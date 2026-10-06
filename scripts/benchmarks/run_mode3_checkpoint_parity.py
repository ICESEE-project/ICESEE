#!/usr/bin/env python3
"""Real-MPI checkpoint/restart parity across two mode-3 process grids."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import sys
import tempfile

from mpi4py import MPI
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from src.parallelization.distributed_adapter import contiguous_state_layout
from src.parallelization.distributed_checkpoint import (
    latest_complete_distributed_checkpoint,
    load_distributed_checkpoint,
    save_distributed_checkpoint,
)
from src.parallelization.distributed_runtime import (
    LocalMemberEnsemble,
    members_for_ensemble_slot,
)
from src.parallelization.distributed_topology import create_distributed_topology


def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nens", type=int, default=7)
    parser.add_argument("--global-size", type=int, default=101)
    parser.add_argument("--source-ensemble-groups", type=int, default=2)
    parser.add_argument("--target-ensemble-groups", type=int, default=1)
    return parser.parse_args()


def _expected(member_id, start, stop):
    rows = np.arange(start, stop, dtype=np.float64)
    return 1000.0 * member_id + rows + 0.001 * rows**2


def main() -> int:
    args = _arguments()
    world = MPI.COMM_WORLD
    if world.Get_size() % args.source_ensemble_groups:
        raise ValueError("world size is incompatible with source process grid")
    if world.Get_size() % args.target_ensemble_groups:
        raise ValueError("world size is incompatible with target process grid")

    source = create_distributed_topology(
        world, ensemble_groups=args.source_ensemble_groups
    )
    source_layout = contiguous_state_layout(
        args.global_size, source, layout_id="checkpoint-parity-v1"
    )
    source_members = members_for_ensemble_slot(
        args.nens, source.ensemble_groups, source.ensemble_slot
    )
    ensemble = LocalMemberEnsemble(
        source_layout,
        {
            member: _expected(
                member, source_layout.owned_start, source_layout.owned_stop
            )
            for member in source_members
        },
    )

    if world.Get_rank() == 0:
        checkpoint_root = tempfile.mkdtemp(prefix="icesee-mode3-checkpoint-")
    else:
        checkpoint_root = None
    checkpoint_root = world.bcast(checkpoint_root, root=0)
    saved = save_distributed_checkpoint(
        checkpoint_root,
        19,
        ensemble,
        source,
        run_id="mode3-checkpoint-parity",
        metadata={"Nens": args.nens, "analysis_cycle": 4, "random_key": 108},
    )

    target = create_distributed_topology(
        world, ensemble_groups=args.target_ensemble_groups
    )
    target_layout = contiguous_state_layout(
        args.global_size, target, layout_id="checkpoint-parity-v1"
    )
    loaded, checkpoint = load_distributed_checkpoint(
        saved.path,
        target_layout,
        target,
        expected_run_id="mode3-checkpoint-parity",
        number_of_members=args.nens,
    )
    local_error = 0.0
    for member, values in loaded.members.items():
        expected = _expected(
            member, target_layout.owned_start, target_layout.owned_stop
        )
        local_error = max(local_error, float(np.max(np.abs(values - expected))))
    maximum_error = world.allreduce(local_error, op=MPI.MAX)
    metadata_ok = (
        checkpoint.timestep == 19
        and checkpoint.metadata.get("analysis_cycle") == 4
        and checkpoint.metadata.get("random_key") == 108
    )
    all_metadata_ok = world.allreduce(bool(metadata_ok), op=MPI.LAND)
    broken_members = dict(ensemble.members)
    if world.Get_rank() == world.Get_size() - 1 and broken_members:
        member = min(broken_members)
        broken = broken_members[member].astype(object)
        if broken.size:
            broken[0] = {"unsupported": "hdf5-object"}
        broken_members[member] = broken
    failure_rejected = False
    try:
        save_distributed_checkpoint(
            checkpoint_root,
            20,
            LocalMemberEnsemble(source_layout, broken_members),
            source,
            run_id="mode3-checkpoint-parity",
            metadata={"Nens": args.nens, "analysis_cycle": 5},
        )
    except RuntimeError:
        failure_rejected = True
    all_failure_rejected = world.allreduce(failure_rejected, op=MPI.LAND)
    latest = latest_complete_distributed_checkpoint(
        checkpoint_root, expected_run_id="mode3-checkpoint-parity"
    )
    latest_ok = latest == saved.path
    all_latest_ok = world.allreduce(latest_ok, op=MPI.LAND)
    passed = bool(
        maximum_error == 0
        and all_metadata_ok
        and all_failure_rejected
        and all_latest_ok
    )
    if world.Get_rank() == 0:
        print("Mode-3 topology-independent checkpoint parity")
        print(
            f"  source grid: {args.source_ensemble_groups} x "
            f"{world.Get_size() // args.source_ensemble_groups}"
        )
        print(
            f"  target grid: {args.target_ensemble_groups} x "
            f"{world.Get_size() // args.target_ensemble_groups}"
        )
        print(f"  maximum state error: {maximum_error:.12g}")
        print(f"  metadata restored: {bool(all_metadata_ok)}")
        print(f"  injected failure rejected: {bool(all_failure_rejected)}")
        print(f"  last complete checkpoint retained: {bool(all_latest_ok)}")
        print(f"  result: {'PASS' if passed else 'FAIL'}")
    world.Barrier()
    if world.Get_rank() == 0:
        shutil.rmtree(checkpoint_root)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
