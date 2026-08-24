#!/usr/bin/env python3
"""Preflight the bounded-memory invariant for execution mode 3."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from src.parallelization.distributed_memory import estimate_mode3_memory


def _gib(value: int) -> float:
    return float(value) / 1024.0**3


def _arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Estimate ICESEE-owned per-rank mode-3 memory without allocating "
            "the requested global state."
        )
    )
    parser.add_argument("--global-rows", type=int, required=True)
    parser.add_argument("--nens", type=int, default=40)
    parser.add_argument("--ensemble-groups", type=int, default=8)
    parser.add_argument("--spatial-ranks", type=int, default=16)
    parser.add_argument("--state-row-chunk", type=int, default=4096)
    parser.add_argument("--observation-row-chunk", type=int, default=4096)
    parser.add_argument("--dtype", choices=("float32", "float64"), default="float64")
    parser.add_argument(
        "--max-rank-gib",
        type=float,
        default=None,
        help="Fail when estimated ICESEE-owned working memory exceeds this value.",
    )
    return parser.parse_args()


def main() -> int:
    args = _arguments()
    plan = estimate_mode3_memory(
        global_rows=args.global_rows,
        total_members=args.nens,
        ensemble_groups=args.ensemble_groups,
        spatial_ranks=args.spatial_ranks,
        state_row_chunk_size=args.state_row_chunk,
        observation_row_chunk_size=args.observation_row_chunk,
        dtype=np.dtype(args.dtype),
    )
    icesee_peak_gib = _gib(plan.estimated_peak_icesee_bytes)
    minimum_peak_gib = _gib(plan.estimated_minimum_peak_bytes)
    passed = args.max_rank_gib is None or minimum_peak_gib <= args.max_rank_gib

    print("Mode-3 analytical memory preflight")
    print(
        f"  process grid:                  {plan.ensemble_groups} x "
        f"{plan.spatial_ranks}"
    )
    print(f"  global state rows/member:      {plan.global_rows:,}")
    print(f"  total/local scheduled members: {plan.total_members}/{plan.local_members}")
    print(f"  maximum owned rows/rank:       {plan.owned_rows:,}")
    print(f"  dtype bytes/value:             {plan.itemsize}")
    print(f"  owned forecast+analysis:       {_gib(plan.owned_snapshots_bytes):.3f} GiB")
    print(
        "  bounded state transform:       "
        f"{_gib(plan.state_transform_workspace_bytes):.3f} GiB"
    )
    print(
        "  bounded observation workspace: "
        f"{_gib(plan.observation_workspace_bytes):.3f} GiB"
    )
    print(f"  ensemble-space products:       {_gib(plan.ensemble_products_bytes):.3f} GiB")
    print(f"  minimum resident native state: {_gib(plan.minimum_native_state_bytes):.3f} GiB")
    print(f"  estimated ICESEE peak/rank:    {icesee_peak_gib:.3f} GiB")
    print(f"  estimated minimum peak/rank:   {minimum_peak_gib:.3f} GiB")
    print(
        "  whole local members avoided:   "
        f"{_gib(plan.complete_local_members_bytes):.3f} GiB/rank"
    )
    print(f"  complete ensemble avoided:      {_gib(plan.complete_ensemble_bytes):.3f} GiB")
    print("  excludes: extra solver vectors, MPI/HDF5 buffers, Python overhead")
    if args.max_rank_gib is not None:
        print(f"  configured rank budget:         {args.max_rank_gib:.3f} GiB")
        print(f"  result:                         {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
