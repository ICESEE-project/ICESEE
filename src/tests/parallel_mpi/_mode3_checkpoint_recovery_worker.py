# ==============================================================================
# @des: MPI worker for test_mode3_checkpoint_recovery_mpi.py. Two phases,
#       launched as separate jobs (possibly with different rank counts):
#
#       write  -- commit steps 0..5 under the rolling retention policy, then
#                 fail the step-6 transaction with ENOSPC on world rank 1
#                 only (a non-root writer), as on a full laptop disk.
#       resume -- collectively discover the newest valid checkpoint, restore
#                 the owned state under the new topology, verify it, and
#                 continue to the end of the run.
#
#       Rank 0 prints one JSON line with the observations.
# ==============================================================================
from __future__ import annotations

import errno
import json
import sys

import h5py
import numpy as np
from mpi4py import MPI

from ICESEE.src.parallelization import distributed_checkpoint
from ICESEE.src.parallelization.distributed_adapter import (
    DistributedBlockStateLayout,
    DistributedStateBlockLayout,
)
from ICESEE.src.parallelization.distributed_checkpoint import (
    committed_checkpoint_paths,
    load_distributed_checkpoint,
)
from ICESEE.src.parallelization.distributed_runtime import (
    LocalMemberEnsemble,
    members_for_ensemble_slot,
)
from ICESEE.src.parallelization.distributed_topology import create_distributed_topology
from ICESEE.src.parallelization.mode3_checkpointing import (
    CheckpointPolicy,
    Mode3CheckpointManager,
)

NENS = 5
NT = 8
ANALYSIS = [2, 4, 7]
BLOCKS = (("a", 0, 9), ("b", 9, 5))
RUN_ID = "mpi-recovery"


def _layout(topology):
    blocks = []
    for name, offset, size in BLOCKS:
        start = size * topology.spatial_rank // topology.spatial_ranks
        stop = size * (topology.spatial_rank + 1) // topology.spatial_ranks
        blocks.append(DistributedStateBlockLayout(name, offset, size, start, stop))
    return DistributedBlockStateLayout(tuple(blocks), layout_id="mpi-recovery-v1")


def _state(layout, member, step):
    """Deterministic owned values: a function of (member, global row, step)."""

    values = []
    for block in layout.blocks:
        rows = np.arange(block.global_offset + block.owned_start,
                         block.global_offset + block.owned_stop, dtype=float)
        values.append(1000.0 * member + rows + 0.25 * step)
    return np.concatenate(values) if values else np.zeros(0)


def _ensemble(topology, layout, step):
    members = members_for_ensemble_slot(NENS, topology.ensemble_groups, topology.ensemble_slot)
    return LocalMemberEnsemble(layout, {m: _state(layout, m, step) for m in members})


def main():
    phase, root, spatial_ranks = sys.argv[1], sys.argv[2], int(sys.argv[3])
    world = MPI.COMM_WORLD
    topology = create_distributed_topology(world, spatial_ranks=spatial_ranks)
    layout = _layout(topology)
    manager = Mode3CheckpointManager(
        root, run_id=RUN_ID, topology=topology, number_of_members=NENS, nt=NT,
        analysis_steps=ANALYSIS, policy=CheckpointPolicy(),
    )
    report = {"world_size": world.Get_size()}

    if phase == "write":
        manager.write_initial(_ensemble(topology, layout, -1))
        km = 0
        for k in range(6):
            km += k in ANALYSIS
            manager.commit_step(k, _ensemble(topology, layout, k),
                                did_analysis=k in ANALYSIS, analyses_completed=km)
        if world.Get_rank() == 1:
            real_file = h5py.File

            def full_disk(path, mode="r", *args, **kwargs):
                if mode == "w":
                    raise OSError(errno.ENOSPC, "No space left on device")
                return real_file(path, mode, *args, **kwargs)

            distributed_checkpoint.h5py.File = full_disk
        try:
            manager.commit_step(6, _ensemble(topology, layout, 6),
                                did_analysis=False, analyses_completed=km)
            report["step6_error"] = None
        except RuntimeError as error:
            report["step6_error"] = str(error)
        report["committed"] = sorted(s for s, _ in committed_checkpoint_paths(manager.steps_root))

    elif phase == "resume":
        plan = manager.prepare_resume()
        restored, checkpoint = load_distributed_checkpoint(
            plan.checkpoint_path, layout, topology, expected_run_id=RUN_ID,
            number_of_members=NENS,
        )
        expected = _ensemble(topology, layout, plan.completed_step)
        exact = all(
            np.array_equal(restored.members[m], expected.members[m]) for m in expected.member_ids
        )
        report["all_ranks_exact"] = bool(world.allreduce(int(exact), op=MPI.MIN))
        report["plan"] = [plan.completed_step, plan.start_step, plan.analyses_completed]
        report["removed"] = list(plan.removed_artifacts)
        km = plan.analyses_completed
        for k in range(plan.start_step, NT):
            km += k in ANALYSIS
            manager.commit_step(k, _ensemble(topology, layout, k),
                                did_analysis=k in ANALYSIS, analyses_completed=km)
        report["committed"] = sorted(s for s, _ in committed_checkpoint_paths(manager.steps_root))
        report["analyses_completed"] = km

    if world.Get_rank() == 0:
        print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
