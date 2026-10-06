from pathlib import Path

import numpy as np
import pytest

from src.parallelization.distributed_adapter import (
    DistributedBlockStateLayout,
    DistributedStateBlockLayout,
    DistributedStateLayout,
)
from src.parallelization.distributed_checkpoint import (
    latest_complete_distributed_checkpoint,
    load_distributed_checkpoint,
    save_distributed_checkpoint,
)
from src.parallelization.distributed_runtime import LocalMemberEnsemble


class _SingleComm:
    def allgather(self, value):
        return [value]

    def gather(self, value, root=0):
        return [value]

    def bcast(self, value, root=0):
        return value

    def Barrier(self):
        return None


class _Topology:
    def __init__(self, *, ensemble_groups=1, ensemble_slot=0):
        self.world = _SingleComm()
        self.spatial_comm = _SingleComm()
        self.world_rank = 0
        self.ensemble_groups = ensemble_groups
        self.ensemble_slot = ensemble_slot
        self.spatial_ranks = 1


def test_checkpoint_roundtrip(tmp_path: Path):
    layout = DistributedStateLayout(8, 0, 8, layout_id="model-state-v2")
    ensemble = LocalMemberEnsemble(
        layout,
        {member: np.arange(8, dtype=float) + 10 * member for member in range(4)},
    )
    saved = save_distributed_checkpoint(
        tmp_path,
        7,
        ensemble,
        _Topology(),
        run_id="run-a",
        metadata={"Nens": 4, "analysis_cycle": 3, "random_key": 91},
    )
    loaded, checkpoint = load_distributed_checkpoint(
        saved.path,
        layout,
        _Topology(),
        expected_run_id="run-a",
        number_of_members=4,
    )
    assert checkpoint.timestep == 7
    assert checkpoint.metadata["analysis_cycle"] == 3
    for member in range(4):
        np.testing.assert_array_equal(loaded.members[member], ensemble.members[member])


def test_checkpoint_loads_only_new_local_slab_and_member_assignment(tmp_path: Path):
    source_layout = DistributedStateLayout(8, 0, 8, layout_id="state")
    source = LocalMemberEnsemble(
        source_layout,
        {member: np.arange(8) + 100 * member for member in range(4)},
    )
    saved = save_distributed_checkpoint(
        tmp_path,
        2,
        source,
        _Topology(),
        run_id="layout-change",
        metadata={"Nens": 4},
    )
    target_layout = DistributedStateLayout(8, 3, 7, layout_id="state")
    loaded, _ = load_distributed_checkpoint(
        saved.path,
        target_layout,
        _Topology(ensemble_groups=2, ensemble_slot=1),
        number_of_members=4,
    )
    assert loaded.member_ids == (1, 3)
    np.testing.assert_array_equal(loaded.members[1], np.arange(3, 7) + 100)
    np.testing.assert_array_equal(loaded.members[3], np.arange(3, 7) + 300)


def _block_layout(h_owned, u_owned, *, layout_id="ice-state"):
    return DistributedBlockStateLayout(
        (
            DistributedStateBlockLayout("h", 0, 6, *h_owned),
            DistributedStateBlockLayout("u", 6, 4, *u_owned),
        ),
        layout_id=layout_id,
    )


def test_block_checkpoint_roundtrip(tmp_path: Path):
    layout = _block_layout((0, 6), (0, 4))
    ensemble = LocalMemberEnsemble(
        layout,
        {member: np.arange(10, dtype=float) + 100 * member for member in range(3)},
    )
    saved = save_distributed_checkpoint(
        tmp_path,
        9,
        ensemble,
        _Topology(),
        run_id="block-run",
        metadata={"Nens": 3, "analysis_cycle": 4},
    )
    loaded, checkpoint = load_distributed_checkpoint(
        saved.path,
        layout,
        _Topology(),
        expected_run_id="block-run",
        number_of_members=3,
    )
    assert checkpoint.timestep == 9
    for member in range(3):
        np.testing.assert_array_equal(loaded.members[member], ensemble.members[member])


def test_block_checkpoint_reads_target_local_overlaps(tmp_path: Path):
    source_layout = _block_layout((0, 6), (0, 4))
    source = LocalMemberEnsemble(
        source_layout,
        {member: np.arange(10) + 100 * member for member in range(4)},
    )
    saved = save_distributed_checkpoint(
        tmp_path,
        3,
        source,
        _Topology(),
        run_id="block-layout-change",
        metadata={"Nens": 4},
    )
    target_layout = _block_layout((2, 5), (1, 4))
    loaded, _ = load_distributed_checkpoint(
        saved.path,
        target_layout,
        _Topology(ensemble_groups=2, ensemble_slot=1),
        number_of_members=4,
    )
    assert loaded.member_ids == (1, 3)
    expected_local_rows = np.asarray([2, 3, 4, 7, 8, 9])
    np.testing.assert_array_equal(loaded.members[1], expected_local_rows + 100)
    np.testing.assert_array_equal(loaded.members[3], expected_local_rows + 300)


def test_block_checkpoint_loads_empty_owned_state(tmp_path: Path):
    source_layout = _block_layout((0, 6), (0, 4))
    saved = save_distributed_checkpoint(
        tmp_path,
        4,
        LocalMemberEnsemble(source_layout, {0: np.arange(10, dtype=np.float32)}),
        _Topology(),
        run_id="empty-target",
        metadata={"Nens": 1},
    )
    empty_layout = _block_layout((3, 3), (2, 2))
    loaded, _ = load_distributed_checkpoint(
        saved.path,
        empty_layout,
        _Topology(),
        number_of_members=1,
    )
    assert loaded.members[0].shape == (0,)
    assert loaded.members[0].dtype == np.dtype(np.float32)


def test_block_checkpoint_rejects_changed_block_definition(tmp_path: Path):
    layout = _block_layout((0, 6), (0, 4))
    saved = save_distributed_checkpoint(
        tmp_path,
        1,
        LocalMemberEnsemble(layout, {0: np.arange(10)}),
        _Topology(),
        run_id="blocks",
        metadata={"Nens": 1},
    )
    incompatible = DistributedBlockStateLayout(
        (
            DistributedStateBlockLayout("h", 0, 5, 0, 5),
            DistributedStateBlockLayout("u", 5, 5, 0, 5),
        ),
        layout_id="ice-state",
    )
    with pytest.raises(ValueError, match="block definitions"):
        load_distributed_checkpoint(
            saved.path,
            incompatible,
            _Topology(),
            number_of_members=1,
        )


def test_checkpoint_rejects_wrong_run_and_layout(tmp_path: Path):
    layout = DistributedStateLayout(3, 0, 3, layout_id="state")
    saved = save_distributed_checkpoint(
        tmp_path,
        1,
        LocalMemberEnsemble(layout, {0: np.arange(3)}),
        _Topology(),
        run_id="correct",
        metadata={"Nens": 1},
    )
    with pytest.raises(ValueError, match="run_id"):
        load_distributed_checkpoint(
            saved.path,
            layout,
            _Topology(),
            expected_run_id="wrong",
            number_of_members=1,
        )
    with pytest.raises(ValueError, match="layout_id"):
        load_distributed_checkpoint(
            saved.path,
            DistributedStateLayout(3, 0, 3, layout_id="other"),
            _Topology(),
            number_of_members=1,
        )


def test_failed_shard_write_is_not_published_as_checkpoint(tmp_path: Path):
    layout = DistributedStateLayout(2, 0, 2, layout_id="state")
    unsupported = np.asarray([object(), object()], dtype=object)
    with pytest.raises(RuntimeError, match="cannot commit distributed checkpoint"):
        save_distributed_checkpoint(
            tmp_path,
            4,
            LocalMemberEnsemble(layout, {0: unsupported}),
            _Topology(),
            run_id="failed-write",
            metadata={"Nens": 1},
        )
    assert latest_complete_distributed_checkpoint(tmp_path) is None
    assert not (tmp_path / "checkpoint_00000004").exists()
    assert not any(path.name.endswith(".staging") for path in tmp_path.iterdir())


def test_latest_checkpoint_ignores_incomplete_and_other_runs(tmp_path: Path):
    layout = DistributedStateLayout(2, 0, 2, layout_id="state")
    ensemble = LocalMemberEnsemble(layout, {0: np.arange(2, dtype=float)})
    first = save_distributed_checkpoint(
        tmp_path,
        2,
        ensemble,
        _Topology(),
        run_id="target",
        metadata={"Nens": 1},
    )
    save_distributed_checkpoint(
        tmp_path,
        5,
        ensemble,
        _Topology(),
        run_id="other",
        metadata={"Nens": 1},
    )
    incomplete = tmp_path / "checkpoint_00000009"
    incomplete.mkdir()
    (tmp_path / ".checkpoint_00000011.target.staging").mkdir()

    assert latest_complete_distributed_checkpoint(
        tmp_path, expected_run_id="target"
    ) == first.path
