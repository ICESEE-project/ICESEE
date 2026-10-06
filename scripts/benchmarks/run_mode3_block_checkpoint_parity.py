#!/usr/bin/env python3
"""MPI parity for segmented mode-3 checkpoints across process grids."""

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

from src.parallelization.distributed_adapter import contiguous_block_state_layout
from src.parallelization.distributed_checkpoint import (
    load_distributed_checkpoint,
    save_distributed_checkpoint,
)
from src.parallelization.distributed_runtime import (
    LocalMemberEnsemble,
    members_for_ensemble_slot,
)
from src.parallelization.distributed_topology import create_distributed_topology


BLOCKS = (("h", 103), ("u", 79), ("v", 79), ("s", 97), ("basal_melt", 3))


def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nens", type=int, default=7)
    parser.add_argument("--source-ensemble-groups", type=int, default=2)
    parser.add_argument("--target-ensemble-groups", type=int, default=1)
    return parser.parse_args()


def _expected(member_id, layout):
    pieces = []
    for block_index, block in enumerate(layout.blocks):
        rows = np.arange(block.owned_start, block.owned_stop, dtype=np.float32)
        pieces.append(
            np.float32(10000 * member_id + 1000 * block_index)
            + rows
            + np.float32(0.001) * rows**2
        )
    return np.ascontiguousarray(np.concatenate(pieces), dtype=np.float32)


def main() -> int:
    args = _arguments()
    world = MPI.COMM_WORLD
    for groups in (args.source_ensemble_groups, args.target_ensemble_groups):
        if groups <= 0 or world.Get_size() % groups:
            raise ValueError("world size is incompatible with a process grid")

    source = create_distributed_topology(
        world, ensemble_groups=args.source_ensemble_groups
    )
    source_layout = contiguous_block_state_layout(
        BLOCKS, source, layout_id="checkpoint-block-parity-v2"
    )
    source_members = members_for_ensemble_slot(
        args.nens, source.ensemble_groups, source.ensemble_slot
    )
    ensemble = LocalMemberEnsemble(
        source_layout,
        {member: _expected(member, source_layout) for member in source_members},
    )

    root = tempfile.mkdtemp(prefix="icesee-mode3-block-") if world.rank == 0 else None
    root = world.bcast(root, root=0)
    saved = save_distributed_checkpoint(
        root,
        23,
        ensemble,
        source,
        run_id="mode3-block-checkpoint-parity",
        metadata={"Nens": args.nens, "analysis_cycle": 6},
    )

    target = create_distributed_topology(
        world, ensemble_groups=args.target_ensemble_groups
    )
    target_layout = contiguous_block_state_layout(
        BLOCKS, target, layout_id="checkpoint-block-parity-v2"
    )
    loaded, checkpoint = load_distributed_checkpoint(
        saved.path,
        target_layout,
        target,
        expected_run_id="mode3-block-checkpoint-parity",
        number_of_members=args.nens,
    )
    local_error = 0.0
    dtype_ok = True
    for member, values in loaded.members.items():
        difference = np.abs(values - _expected(member, target_layout))
        local_error = max(local_error, float(np.max(difference, initial=0.0)))
        dtype_ok = dtype_ok and values.dtype == np.dtype(np.float32)
    maximum_error = world.allreduce(local_error, op=MPI.MAX)
    all_dtype_ok = world.allreduce(dtype_ok, op=MPI.LAND)
    metadata_ok = world.allreduce(
        checkpoint.timestep == 23
        and checkpoint.metadata.get("analysis_cycle") == 6,
        op=MPI.LAND,
    )
    passed = maximum_error == 0.0 and all_dtype_ok and metadata_ok

    if world.rank == 0:
        print("Mode-3 segmented checkpoint parity")
        print(
            f"  source grid: {args.source_ensemble_groups} x "
            f"{world.size // args.source_ensemble_groups}"
        )
        print(
            f"  target grid: {args.target_ensemble_groups} x "
            f"{world.size // args.target_ensemble_groups}"
        )
        print(f"  blocks: {dict(BLOCKS)}")
        print(f"  maximum state error: {maximum_error:.12g}")
        print(f"  float32 preserved: {bool(all_dtype_ok)}")
        print(f"  metadata restored: {bool(metadata_ok)}")
        print(f"  result: {'PASS' if passed else 'FAIL'}")
    world.Barrier()
    if world.rank == 0:
        shutil.rmtree(root)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
