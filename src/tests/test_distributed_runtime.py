"""Tests for non-selectable execution-mode-3 state-only primitives."""

from types import SimpleNamespace

import numpy as np
import pytest

from src.parallelization.distributed_adapter import (
    DistributedBlockStateLayout,
    DistributedStateBlockLayout,
    DistributedStateLayout,
    contiguous_state_layout,
)
from src.parallelization.distributed_runtime import (
    LocalMemberEnsemble,
    apply_ensemble_transform_local,
    forecast_local_members,
    initialize_local_members,
    members_for_ensemble_slot,
    reconstruct_global_member_for_parity,
    transform_local_members,
)


class _SpatialComm:
    def __init__(self, rank=0, descriptors=None, gathered=None):
        self.rank = rank
        self.descriptors = descriptors
        self.gathered = gathered

    def Get_rank(self):
        return self.rank

    def allgather(self, value):
        return [value] if self.descriptors is None else list(self.descriptors)

    def gather(self, value, root=0):
        if self.rank != root:
            return None
        return [value] if self.gathered is None else list(self.gathered)


class _EnsembleComm:
    def allgather(self, value):
        return [value]


class _Adapter:
    def distributed_state_layout(self, *, topology, icesee_kwargs):
        return contiguous_state_layout(
            icesee_kwargs["nd"], topology, layout_id="test-state"
        )

    def initialize_local_member(
        self, member_id, *, layout, topology, icesee_kwargs
    ):
        values = np.arange(layout.global_size, dtype=float) + 100.0 * member_id
        return values[layout.owned_slice]

    def forecast_local_member(
        self, local_state, member_id, timestep, *, layout, topology, icesee_kwargs
    ):
        return local_state + member_id + timestep

    def observe_local_member(self, *args, **kwargs):
        return np.empty(0)

    def finalize_local_analysis(self, local_forecast, local_analysis, *args, **kwargs):
        return local_analysis


def _topology(slot=0, spatial_rank=0, ensemble_groups=2, spatial_ranks=1, comm=None):
    return SimpleNamespace(
        ensemble_slot=slot,
        spatial_rank=spatial_rank,
        ensemble_groups=ensemble_groups,
        spatial_ranks=spatial_ranks,
        spatial_comm=_SpatialComm() if comm is None else comm,
        ensemble_comm=_EnsembleComm(),
    )


def test_round_robin_member_schedule_is_complete_and_disjoint():
    assignments = [members_for_ensemble_slot(8, 3, slot) for slot in range(3)]
    assert assignments == [(0, 3, 6), (1, 4, 7), (2, 5)]
    assert sorted(member for group in assignments for member in group) == list(range(8))


def test_balanced_contiguous_layout_handles_remainder():
    layouts = [
        contiguous_state_layout(10, _topology(spatial_rank=rank, spatial_ranks=3))
        for rank in range(3)
    ]
    assert [(item.owned_start, item.owned_stop) for item in layouts] == [
        (0, 4),
        (4, 7),
        (7, 10),
    ]


def test_initialize_and_forecast_preserve_member_ids_and_local_shapes():
    topology = _topology(slot=1)
    initial = initialize_local_members(_Adapter(), topology, {"Nens": 5, "nd": 6})
    assert initial.member_ids == (1, 3)
    np.testing.assert_array_equal(initial.members[1], np.arange(6) + 100)

    forecast = forecast_local_members(_Adapter(), initial, 4, topology, {})
    assert forecast.member_ids == initial.member_ids
    np.testing.assert_array_equal(forecast.members[1], np.arange(6) + 105)
    np.testing.assert_array_equal(forecast.members[3], np.arange(6) + 307)


def test_local_member_ensemble_rejects_wrong_slab_shape():
    with pytest.raises(ValueError, match="local state"):
        LocalMemberEnsemble(
            layout=DistributedStateLayout(8, 0, 4),
            members={0: np.zeros(3)},
        )


def test_local_transform_matches_dense_right_multiplication():
    topology = _topology(slot=0, ensemble_groups=1)
    local = LocalMemberEnsemble(
        layout=DistributedStateLayout(5, 0, 5),
        members={
            0: np.arange(5, dtype=float),
            1: np.arange(5, dtype=float) + 10.0,
        },
    )
    transform = np.array([[0.75, 0.20], [0.25, 0.80]])
    updated = apply_ensemble_transform_local(
        _Adapter(), local, transform, 3, topology, {"Nens": 2},
        row_chunk_size=2,
    )
    dense = np.column_stack([local.members[0], local.members[1]]) @ transform
    np.testing.assert_allclose(updated.members[0], dense[:, 0])
    np.testing.assert_allclose(updated.members[1], dense[:, 1])


def test_pure_local_transform_matches_finalized_array_adapter_path():
    topology = _topology(slot=0, ensemble_groups=1)
    local = LocalMemberEnsemble(
        DistributedStateLayout(4, 0, 4),
        {0: np.arange(4.0), 1: np.arange(4.0) + 4.0},
    )
    transform = np.array([[0.7, 0.1], [0.3, 0.9]])
    pure = transform_local_members(
        local, transform, topology, {"Nens": 2}, row_chunk_size=1
    )
    finalized = apply_ensemble_transform_local(
        _Adapter(), local, transform, 1, topology, {"Nens": 2},
        row_chunk_size=1,
    )
    for member_id in local.member_ids:
        np.testing.assert_array_equal(
            pure.members[member_id], finalized.members[member_id]
        )


def test_local_transform_accepts_segmented_multivariable_storage():
    topology = _topology(slot=0, ensemble_groups=1)
    layout = DistributedBlockStateLayout((
        DistributedStateBlockLayout("h", 0, 8, 0, 3),
        DistributedStateBlockLayout("u", 8, 5, 0, 2),
    ))
    local = LocalMemberEnsemble(
        layout,
        {0: np.arange(5, dtype=float), 1: np.arange(5, dtype=float) + 20.0},
    )
    transform = np.array([[0.8, 0.1], [0.2, 0.9]])

    updated = apply_ensemble_transform_local(
        _Adapter(), local, transform, 2, topology, {"Nens": 2}, row_chunk_size=2
    )

    expected = np.column_stack([local.members[0], local.members[1]]) @ transform
    np.testing.assert_allclose(updated.members[0], expected[:, 0])
    np.testing.assert_allclose(updated.members[1], expected[:, 1])


def test_parity_reconstruction_is_root_only_and_checks_partition():
    layout = DistributedStateLayout(6, 0, 2)
    root_comm = _SpatialComm(
        rank=0,
        gathered=[(0, 2, [0.0, 1.0]), (2, 4, [2.0, 3.0]), (4, 6, [4.0, 5.0])],
    )
    reconstructed = reconstruct_global_member_for_parity(
        np.array([0.0, 1.0]), layout, root_comm
    )
    np.testing.assert_array_equal(reconstructed, np.arange(6, dtype=float))

    nonroot = _SpatialComm(rank=1)
    assert reconstruct_global_member_for_parity(
        np.array([0.0, 1.0]), layout, nonroot, root=0
    ) is None

    broken = _SpatialComm(
        rank=0,
        gathered=[(0, 2, [0.0, 1.0]), (3, 6, [3.0, 4.0, 5.0])],
    )
    with pytest.raises(ValueError, match="exact partition"):
        reconstruct_global_member_for_parity(
            np.array([0.0, 1.0]), layout, broken
        )
