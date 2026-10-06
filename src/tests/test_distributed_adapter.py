"""Unit tests for the proposed execution-mode-3 adapter contract."""

from __future__ import annotations

import numpy as np
import pytest

from src.parallelization.distributed_adapter import (
    DistributedAnalysisTargets,
    DistributedBlockStateLayout,
    DistributedObservationCoordinates,
    DistributedStateBlockLayout,
    DistributedStateLayout,
    contiguous_block_state_layout,
    validate_distributed_adapter,
    validate_distributed_analysis_coordinates,
    validate_block_spatial_partition,
    validate_spatial_partition,
)


class _GatherComm:
    def __init__(self, descriptors):
        self.descriptors = descriptors

    def allgather(self, _value):
        return list(self.descriptors)


class _Topology:
    spatial_ranks = 3

    def __init__(self, spatial_rank):
        self.spatial_rank = spatial_rank


def test_layout_normalizes_ghosts_and_exports_checkpoint_metadata():
    layout = DistributedStateLayout(
        global_size=12,
        owned_start=4,
        owned_stop=8,
        ghost_indices=np.array([9, 3, 3]),
        layout_id="mesh-state-v2",
    )

    assert layout.owned_size == 4
    assert layout.owned_slice == slice(4, 8)
    np.testing.assert_array_equal(layout.ghost_indices, [3, 9])
    assert layout.checkpoint_metadata() == {
        "layout_id": "mesh-state-v2",
        "global_size": 12,
        "owned_start": 4,
        "owned_stop": 8,
    }


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"global_size": -1, "owned_start": 0, "owned_stop": 0}, "nonnegative"),
        ({"global_size": 8, "owned_start": 5, "owned_stop": 4}, "owned slab"),
        (
            {
                "global_size": 8,
                "owned_start": 2,
                "owned_stop": 5,
                "ghost_indices": [4],
            },
            "overlap",
        ),
        (
            {
                "global_size": 8,
                "owned_start": 2,
                "owned_stop": 5,
                "ghost_indices": [8],
            },
            "global state",
        ),
    ],
)
def test_layout_rejects_invalid_ownership(kwargs, message):
    with pytest.raises(ValueError, match=message):
        DistributedStateLayout(**kwargs)


def test_block_layout_preserves_variable_major_global_rows_with_local_storage():
    layout = DistributedBlockStateLayout(
        (
            DistributedStateBlockLayout("h", 0, 10, 4, 7),
            DistributedStateBlockLayout("u", 10, 7, 3, 5),
            DistributedStateBlockLayout("bed", 17, 10, 4, 7),
        ),
        layout_id="icepack-pig-v1",
    )

    assert layout.global_size == 27
    assert layout.owned_size == 8
    assert layout.block_names == ("h", "u", "bed")
    assert layout.local_slice("h") == slice(0, 3)
    assert layout.local_slice("u") == slice(3, 5)
    assert layout.local_slice("bed") == slice(5, 8)
    assert list(layout.iter_owned_global_intervals()) == [
        ("h", slice(0, 3), (4, 7)),
        ("u", slice(3, 5), (13, 15)),
        ("bed", slice(5, 8), (21, 24)),
    ]
    np.testing.assert_array_equal(
        layout.local_to_global_rows(np.array([0, 2, 3, 4, 5, 7])),
        np.array([4, 6, 13, 14, 21, 23]),
    )


def test_contiguous_block_layout_balances_each_variable_independently():
    layout = contiguous_block_state_layout(
        {"h": 10, "velocity": 17},
        _Topology(spatial_rank=1),
        layout_id="mixed-functions-v1",
    )

    assert layout.block("h").global_owned_interval == (4, 7)
    assert layout.block("velocity").global_owned_interval == (16, 22)
    assert layout.owned_size == 9


def test_block_partition_accepts_exact_per_variable_coverage():
    layout = DistributedBlockStateLayout(
        (
            DistributedStateBlockLayout("h", 0, 10, 0, 4),
            DistributedStateBlockLayout("u", 10, 7, 0, 3),
        ),
        layout_id="blocks-v1",
    )
    comm = _GatherComm([
        ("blocks-v1", (("h", 0, 10, 0, 4), ("u", 10, 7, 0, 3))),
        ("blocks-v1", (("h", 0, 10, 4, 7), ("u", 10, 7, 3, 5))),
        ("blocks-v1", (("h", 0, 10, 7, 10), ("u", 10, 7, 5, 7))),
    ])

    assert validate_block_spatial_partition(layout, comm) == {
        "h": ((0, 4), (4, 7), (7, 10)),
        "u": ((0, 3), (3, 5), (5, 7)),
    }


def test_block_partition_rejects_different_model_definitions():
    layout = DistributedBlockStateLayout(
        (DistributedStateBlockLayout("h", 0, 10, 0, 5),),
        layout_id="blocks-v1",
    )
    comm = _GatherComm([
        ("blocks-v1", (("h", 0, 10, 0, 5),)),
        ("blocks-v1", (("thickness", 0, 10, 5, 10),)),
    ])

    with pytest.raises(ValueError, match="block definitions"):
        validate_block_spatial_partition(layout, comm)


def test_spatial_partition_accepts_exact_global_coverage():
    layout = DistributedStateLayout(12, 4, 8, layout_id="state-v1")
    comm = _GatherComm(
        [
            (0, 4, 12, "state-v1"),
            (4, 8, 12, "state-v1"),
            (8, 12, 12, "state-v1"),
        ]
    )

    assert validate_spatial_partition(layout, comm) == ((0, 4), (4, 8), (8, 12))


@pytest.mark.parametrize(
    "descriptors, message",
    [
        ([(0, 5, 12, "v1"), (4, 12, 12, "v1")], "overlaps"),
        ([(0, 4, 12, "v1"), (5, 12, 12, "v1")], "gap"),
        ([(0, 4, 12, "v1"), (4, 11, 12, "v1")], "cover"),
        ([(0, 4, 12, "v1"), (4, 13, 13, "v1")], "global state size"),
        ([(0, 4, 12, "v1"), (4, 12, 12, "v2")], "layout_id"),
    ],
)
def test_spatial_partition_rejects_invalid_global_ownership(descriptors, message):
    layout = DistributedStateLayout(12, 0, 4, layout_id="v1")
    with pytest.raises(ValueError, match=message):
        validate_spatial_partition(layout, _GatherComm(descriptors))


class _CompleteAdapter:
    def distributed_state_layout(self, **kwargs):
        return None

    def initialize_local_member(self, *args, **kwargs):
        return None

    def forecast_local_member(self, *args, **kwargs):
        return None

    def observe_local_member(self, *args, **kwargs):
        return None

    def finalize_local_analysis(self, *args, **kwargs):
        return None

    def inverse_local_member(self, *args, **kwargs):
        return None


def test_adapter_validation_separates_state_and_hybrid_contracts():
    validate_distributed_adapter(_CompleteAdapter())
    validate_distributed_adapter(_CompleteAdapter(), require_inversion=True)

    incomplete = _CompleteAdapter()
    incomplete.inverse_local_member = None
    validate_distributed_adapter(incomplete)
    with pytest.raises(TypeError, match="inverse_local_member"):
        validate_distributed_adapter(incomplete, require_inversion=True)


def test_localization_metadata_supports_multidimensional_coordinates():
    layout = DistributedStateLayout(10, 4, 10)
    targets = [DistributedAnalysisTargets(
        name="thickness",
        local_rows=np.array([0, 2, 5]),
        global_rows=np.array([4, 6, 9]),
        coordinates=np.array([[0.0, 1.0], [2.0, 3.0], [5.0, 8.0]]),
        observation_kinds=("thickness", "velocity"),
    )]
    observations = [
        DistributedObservationCoordinates(
            "thickness", np.array([11]), np.array([[0.0, 1.0]])
        ),
        DistributedObservationCoordinates(
            "velocity", np.array([19]), np.array([[2.0, 3.0]])
        ),
    ]

    validate_distributed_analysis_coordinates(targets, observations, layout)


def test_localization_metadata_requires_ids_unique_across_kinds():
    layout = DistributedStateLayout(2, 0, 2)
    targets = [DistributedAnalysisTargets(
        "state", np.array([0]), np.array([0]), np.array([[0.0]]), ("h", "u")
    )]
    observations = [
        DistributedObservationCoordinates("h", np.array([9]), np.array([[0.0]])),
        DistributedObservationCoordinates("u", np.array([9]), np.array([[1.0]])),
    ]
    with pytest.raises(ValueError, match="unique across kinds"):
        validate_distributed_analysis_coordinates(targets, observations, layout)


def test_localization_metadata_supports_segmented_block_rows():
    layout = DistributedBlockStateLayout(
        (
            DistributedStateBlockLayout("h", 0, 10, 4, 7),
            DistributedStateBlockLayout("u", 10, 7, 3, 5),
        )
    )
    targets = [DistributedAnalysisTargets(
        "mixed",
        np.array([0, 2, 3, 4]),
        np.array([4, 6, 13, 14]),
        np.array([[0.0], [1.0], [2.0], [3.0]]),
        ("speed",),
    )]
    observations = [DistributedObservationCoordinates(
        "speed", np.array([2]), np.array([[1.0]])
    )]
    validate_distributed_analysis_coordinates(targets, observations, layout)


def test_localization_metadata_rejects_inconsistent_global_rows():
    layout = DistributedStateLayout(10, 4, 10)
    targets = [DistributedAnalysisTargets(
        "bed", np.array([0]), np.array([5]), np.array([[0.0, 1.0]]), ("bed",)
    )]
    observations = [DistributedObservationCoordinates(
        "bed", np.array([3]), np.array([[0.0, 1.0]])
    )]

    with pytest.raises(ValueError, match="inconsistent global rows"):
        validate_distributed_analysis_coordinates(targets, observations, layout)
