# ==============================================================================
# @des: Standalone comparator for two Idealized PIG Mode-3 distributed
# checkpoints (2026-09-28, PACE calibration prep, item E). Compares
# h, u, v, s, basal_melt_field for matching member_ids between two
# checkpoint directories written by src/parallelization/distributed_checkpoint.py
# (format "icesee-distributed-block-checkpoint-v2": one manifest.json +
# one shard_*.h5 per spatial rank, each shard holding one or more
# member_<id>/<field> datasets plus that shard's owned_start/owned_stop
# per block).
#
# Reassembly (always valid, no caveats): within ONE checkpoint directory,
# each member's global-order array for a field is built by placing every
# shard's member_<id>/<field> values at [owned_start:owned_stop) -- this
# is exact regardless of how many spatial ranks (P_model) wrote it, since
# the manifest records each shard's owned range explicitly.
#
# Cross-run caveat (READ BEFORE TRUSTING A DIRECT DIFF): Firedrake/PETSc
# mesh partitioning is NOT reproducible across separate process
# invocations (confirmed empirically earlier in this engagement -- same
# mesh file, same rank count, different runs can still get different
# rank-to-DOF assignments). Two checkpoints from the SAME run (e.g. two
# timesteps of one invocation, or comparing this run's own initial vs
# later checkpoint) are safe to diff by raw global DOF index. Two
# checkpoints from SEPARATE invocations (e.g. a P_model=1 run vs a
# P_model=2 run, or a memory-backend run vs an HDF5-backend run) are only
# safe to diff by raw index if you have independently confirmed both runs
# partitioned identically -- otherwise pass --coords-a/--coords-b (two
# .npy files of shape (global_size_per_field, 2), produced by a small
# companion coordinate-dump script) to remap run B onto run A's DOF order
# via a nearest-neighbor coordinate match (reusing
# src/tests/parallel_mpi/_coordinate_aware_compare.py's match_by_coordinate).
# This script refuses a cross-run comparison without --coords-a/--coords-b
# unless --i-confirmed-identical-partitioning is passed explicitly.
# ==============================================================================
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_REPO_ROOT), str(_REPO_ROOT.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

_FIELDS = ("h", "u", "v", "s", "basal_melt_field")


def load_checkpoint(directory: Path) -> dict[int, dict[str, np.ndarray]]:
    """Return {member_id: {field_name: global-order np.ndarray}}."""

    manifest = json.loads((directory / "manifest.json").read_text())
    global_size = {b["name"]: b["global_size"] for b in manifest["blocks"]}
    members: dict[int, dict[str, np.ndarray]] = {}

    for shard in manifest["shards"]:
        shard_path = directory / shard["file"]
        block_by_name = {b["name"]: b for b in shard["blocks"]}
        with h5py.File(shard_path, "r") as f:
            for member_id in shard["member_ids"]:
                group = f[f"member_{member_id}"]
                out = members.setdefault(
                    int(member_id),
                    {name: np.full(size, np.nan) for name, size in global_size.items()},
                )
                for field in _FIELDS:
                    if field not in group:
                        continue
                    block = block_by_name[field]
                    start, stop = block["owned_start"], block["owned_stop"]
                    out[field][start:stop] = group[field][:]
    return members


def compare(
    dir_a: Path, dir_b: Path, *, coords_a: np.ndarray | None, coords_b: np.ndarray | None,
    atol: float, rtol: float,
) -> dict:
    members_a = load_checkpoint(dir_a)
    members_b = load_checkpoint(dir_b)

    common_ids = sorted(set(members_a) & set(members_b))
    if set(members_a) != set(members_b):
        print(f"WARNING: member_id sets differ -- A only: {set(members_a)-set(members_b)}, B only: {set(members_b)-set(members_a)}")

    perm = None
    if coords_a is not None and coords_b is not None:
        from src.tests.parallel_mpi._coordinate_aware_compare import match_by_coordinate
        perm = match_by_coordinate(coords_a, coords_b)

    report = {}
    for member_id in common_ids:
        field_report = {}
        for field in _FIELDS:
            a = members_a[member_id].get(field)
            b = members_b[member_id].get(field)
            if a is None or b is None:
                field_report[field] = {"skipped": "field absent in one checkpoint"}
                continue
            b_reordered = b[perm] if perm is not None else b
            if a.shape != b_reordered.shape:
                field_report[field] = {"error": f"shape mismatch {a.shape} vs {b_reordered.shape}"}
                continue
            diff = np.abs(a - b_reordered)
            within_tol = np.allclose(a, b_reordered, atol=atol, rtol=rtol, equal_nan=True)
            field_report[field] = {
                "max_abs_diff": float(np.nanmax(diff)) if diff.size else 0.0,
                "rms_diff": float(np.sqrt(np.nanmean(diff ** 2))) if diff.size else 0.0,
                "within_tolerance": bool(within_tol),
            }
        report[member_id] = field_report
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint_dir_a", type=Path)
    parser.add_argument("checkpoint_dir_b", type=Path)
    parser.add_argument("--coords-a", type=Path, default=None, help=".npy coordinates for run A (required for a cross-run comparison unless --i-confirmed-identical-partitioning)")
    parser.add_argument("--coords-b", type=Path, default=None)
    parser.add_argument("--i-confirmed-identical-partitioning", action="store_true", help="Skip the coordinate-remap requirement -- only pass this if you have independently verified both runs partitioned the mesh identically (e.g. both are shards from the SAME invocation).")
    parser.add_argument("--atol", type=float, default=1e-10)
    parser.add_argument("--rtol", type=float, default=1e-8)
    args = parser.parse_args()

    if args.coords_a is None or args.coords_b is None:
        if not args.i_confirmed_identical_partitioning:
            print(
                "REFUSING: no --coords-a/--coords-b given, and "
                "--i-confirmed-identical-partitioning not passed. A direct "
                "index comparison across two SEPARATE process invocations "
                "is only valid if partitioning is identical, which is NOT "
                "guaranteed by Firedrake/PETSc (see this file's module "
                "docstring). Pass the flag explicitly if you have confirmed "
                "this (e.g. both checkpoints are from steps of the SAME run).",
                file=sys.stderr,
            )
            sys.exit(2)
        coords_a = coords_b = None
    else:
        coords_a = np.load(args.coords_a)
        coords_b = np.load(args.coords_b)

    report = compare(
        args.checkpoint_dir_a, args.checkpoint_dir_b,
        coords_a=coords_a, coords_b=coords_b, atol=args.atol, rtol=args.rtol,
    )
    all_ok = all(
        entry.get("within_tolerance", False)
        for member in report.values()
        for entry in member.values()
        if "within_tolerance" in entry
    )
    print(json.dumps(report, indent=2))
    print(f"\nALL FIELDS WITHIN TOLERANCE: {all_ok}")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
