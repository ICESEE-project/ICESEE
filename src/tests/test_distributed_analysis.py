"""Tests for bounded-memory execution-mode-3 analysis contracts."""

from types import SimpleNamespace

import numpy as np
import pytest

from src.parallelization.distributed_analysis import (
    assemble_selected_ensemble_rows,
    DistributedObservationLayout,
    StochasticAnalysisProducts,
    contiguous_observation_layout,
    ensemble_transform_from_products,
    iter_local_ensemble_row_blocks,
    validate_observation_partition,
)


class _SpatialComm:
    def __init__(self, descriptors=None):
        self.descriptors = descriptors

    def allgather(self, value):
        return [value] if self.descriptors is None else list(self.descriptors)


class _EnsembleComm:
    def __init__(self, remote_members):
        self.remote_members = {
            int(key): np.asarray(value) for key, value in remote_members.items()
        }
        self.offset = 0

    def allgather(self, value):
        if not isinstance(value, dict):
            remote_size = next(iter(self.remote_members.values())).size
            return [value, remote_size]
        width = len(next(iter(value.values()))) if value else 0
        remote = {
            key: array[self.offset:self.offset + width]
            for key, array in self.remote_members.items()
        }
        self.offset += width
        return [value, remote]


class _ReduceComm:
    def Allreduce(self, local, target):
        target[...] = 2.0 * local

    def allgather(self, value):
        return [value, value]


def _topology(spatial_rank=0, spatial_ranks=3):
    return SimpleNamespace(spatial_rank=spatial_rank, spatial_ranks=spatial_ranks)


def test_observation_layout_is_balanced_and_preserves_global_ids():
    rows = np.array([9, 2, 14, 7, 21, 4, 18])
    layouts = [
        contiguous_observation_layout(rows, _topology(rank, 3))
        for rank in range(3)
    ]
    assert [(item.owned_start, item.owned_stop) for item in layouts] == [
        (0, 3), (3, 5), (5, 7)
    ]
    np.testing.assert_array_equal(
        np.concatenate([item.global_row_ids for item in layouts]), rows
    )


def test_observation_partition_detects_gap_and_overlap():
    layout = DistributedObservationLayout(5, 0, 2, np.array([3, 8]))
    valid = [(0, 2, 5, "observations-v1"), (2, 5, 5, "observations-v1")]
    assert validate_observation_partition(layout, _SpatialComm(valid)) == (
        (0, 2), (2, 5)
    )
    with pytest.raises(ValueError, match="gap"):
        validate_observation_partition(
            layout,
            _SpatialComm([
                (0, 2, 5, "observations-v1"),
                (3, 5, 5, "observations-v1"),
            ]),
        )
    with pytest.raises(ValueError, match="overlaps"):
        validate_observation_partition(
            layout,
            _SpatialComm([
                (0, 3, 5, "observations-v1"),
                (2, 5, 5, "observations-v1"),
            ]),
        )


def test_local_ensemble_blocks_follow_member_ids_and_bound_rows():
    local = {
        0: np.arange(7, dtype=float),
        2: np.arange(7, dtype=float) + 20,
    }
    remote = {
        1: np.arange(7, dtype=float) + 10,
        3: np.arange(7, dtype=float) + 30,
    }
    blocks = list(iter_local_ensemble_row_blocks(
        local, 4, _EnsembleComm(remote), row_chunk_size=3
    ))
    assert [item[1].shape for item in blocks] == [(3, 4), (3, 4), (1, 4)]
    assembled = np.vstack([item[1] for item in blocks])
    expected = np.column_stack([
        np.arange(7, dtype=float) + 10 * member for member in range(4)
    ])
    np.testing.assert_array_equal(assembled, expected)


def test_selected_ensemble_rows_follow_stable_member_ids():
    rows = np.array([5, 1, 6])
    local = {
        0: np.arange(7, dtype=float),
        2: np.arange(7, dtype=float) + 20,
    }
    remote = {
        1: (np.arange(7, dtype=float) + 10)[rows],
        3: (np.arange(7, dtype=float) + 30)[rows],
    }

    assembled = assemble_selected_ensemble_rows(
        local, rows, 4, _EnsembleComm(remote)
    )
    expected = np.column_stack([
        np.arange(7, dtype=float)[rows] + 10 * member for member in range(4)
    ])
    np.testing.assert_array_equal(assembled, expected)


def test_streamed_products_equal_dense_legacy_transform():
    rng = np.random.default_rng(91)
    forecast = rng.normal(size=(13, 5))
    observations = rng.normal(size=13)
    products = StochasticAnalysisProducts.zeros(5)
    for start in range(0, 13, 4):
        products.add_chunk(
            forecast[start:start + 4], observations[start:start + 4],
            error_mode="legacy_prior_anomalies",
        )

    yprime = forecast - forecast.mean(axis=1, keepdims=True)
    factor = 2.0 * yprime
    innovation = observations[:, None] - forecast
    np.testing.assert_allclose(products.cross, yprime.T @ factor)
    np.testing.assert_allclose(products.gram, factor.T @ factor)
    np.testing.assert_allclose(products.rhs, factor.T @ innovation)

    dense_products = StochasticAnalysisProducts(
        yprime.T @ factor,
        factor.T @ factor,
        factor.T @ innovation,
        observation_rows=13,
    )
    np.testing.assert_allclose(
        ensemble_transform_from_products(products),
        ensemble_transform_from_products(dense_products),
        rtol=1.0e-12,
        atol=1.0e-12,
    )


@pytest.mark.parametrize("mode", ["stochastic_R", "generated_R"])
def test_streamed_products_equal_dense_error_factor_modes(mode):
    rng = np.random.default_rng(193)
    forecast = rng.normal(size=(11, 4))
    observations = rng.normal(size=11)
    eta = rng.normal(scale=0.2, size=forecast.shape)
    eta -= eta.mean(axis=1, keepdims=True)
    products = StochasticAnalysisProducts.zeros(4)
    for start in range(0, 11, 3):
        products.add_chunk(
            forecast[start:start + 3], observations[start:start + 3],
            error_mode=mode,
            observation_errors=eta[start:start + 3],
        )

    yprime = forecast - forecast.mean(axis=1, keepdims=True)
    factor = yprime + eta
    innovation = observations[:, None] + eta - forecast
    np.testing.assert_allclose(products.cross, yprime.T @ factor)
    np.testing.assert_allclose(products.gram, factor.T @ factor)
    np.testing.assert_allclose(products.rhs, factor.T @ innovation)


def test_products_reduce_only_fixed_ensemble_matrices():
    products = StochasticAnalysisProducts.zeros(3)
    products.cross[:] = 1.0
    products.gram[:] = 2.0
    products.rhs[:] = 3.0
    products.observation_rows = 7
    reduced = products.reduced(_ReduceComm())
    np.testing.assert_array_equal(reduced.cross, np.full((3, 3), 2.0))
    np.testing.assert_array_equal(reduced.gram, np.full((3, 3), 4.0))
    np.testing.assert_array_equal(reduced.rhs, np.full((3, 3), 6.0))
    assert reduced.observation_rows == 14
