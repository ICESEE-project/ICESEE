"""Tests for mode-3 distributed grouped-local stochastic analysis."""

import numpy as np

from src.parallelization.distributed_local_analysis import (
    apply_distributed_local_patches,
    iter_distributed_local_patches,
)
from src.utils.localization import compute_X5_from_matrices


class _SingleRankComm:
    def Get_size(self):
        return 1

    def allgather(self, value):
        return [value]

    def alltoall(self, values):
        assert len(values) == 1
        return list(values)


def test_distributed_local_patches_match_existing_grouped_analysis():
    rng = np.random.default_rng(20260821)
    target_coords = np.array([[0.0], [0.1], [1.0], [1.1], [4.0]])
    target_rows = np.array([10, 11, 12, 13, 14])
    obs_coords = np.array([[0.0], [0.2], [1.0], [1.2], [8.0]])
    obs_ids = np.array([30, 10, 50, 20, 90])
    yprime = rng.normal(size=(5, 4))
    yprime -= yprime.mean(axis=1, keepdims=True)
    eta = rng.normal(scale=0.1, size=(5, 4))
    eta -= eta.mean(axis=1, keepdims=True)
    innovations = rng.normal(size=(5, 4))

    rounds = list(iter_distributed_local_patches(
        target_global_rows=target_rows,
        target_coordinates=target_coords,
        observation_ids=obs_ids,
        observation_coordinates=obs_coords,
        yprime=yprime,
        eta=eta,
        innovations=innovations,
        radius=0.25,
        spatial_comm=_SingleRankComm(),
        target_chunk_size=3,
    ))
    assert [item[0] for item in rounds] == [slice(0, 3), slice(3, 5)]

    forecast = rng.normal(size=(5, 4))
    analysis = forecast.copy()
    for rows, patches in rounds:
        apply_distributed_local_patches(
            analysis[rows], forecast[rows], patches
        )
        for patch in patches:
            order = np.array([
                int(np.flatnonzero(obs_ids == obs_id)[0])
                for obs_id in patch.observation_ids
            ])
            expected = compute_X5_from_matrices(
                yprime[order], eta[order], innovations[order], 4
            )
            np.testing.assert_allclose(
                patch.transform, expected, rtol=1.0e-12, atol=1.0e-12
            )

    # First pair and second pair are updated by different local neighborhoods;
    # the final point has no observations and retains the global/prior value.
    assert not np.allclose(analysis[:2], forecast[:2])
    assert not np.allclose(analysis[2:4], forecast[2:4])
    np.testing.assert_array_equal(analysis[4], forecast[4])


def test_empty_observation_owner_participates_without_patches():
    rounds = list(iter_distributed_local_patches(
        target_global_rows=np.array([4, 5]),
        target_coordinates=np.array([[0.0, 0.0], [1.0, 1.0]]),
        observation_ids=np.empty(0, dtype=np.int64),
        observation_coordinates=np.empty((0, 2)),
        yprime=np.empty((0, 3)),
        eta=np.empty((0, 3)),
        innovations=np.empty((0, 3)),
        radius=1.0,
        spatial_comm=_SingleRankComm(),
        target_chunk_size=1,
    ))
    assert len(rounds) == 2
    assert all(not patches for _, patches in rounds)
