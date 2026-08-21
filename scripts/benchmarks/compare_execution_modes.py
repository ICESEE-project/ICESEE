#!/usr/bin/env python3
"""Compare partial-parallel and fully-parallel ICESEE ensemble histories.

The partial-parallel source can be ``icesee_ensemble_data.h5`` or its run
directory.  The fully-parallel source is normally a run directory containing
``icesee_enkf_ens_XXXX.h5`` shards.  Comparison is row-chunked, so the checker
does not reconstruct a large ensemble in memory.

The default gate is deliberately strict numerical parity rather than bitwise
identity.  Modes 1 and 2 reduce floating-point products in different MPI
orders, so bitwise identity is not portable even when their algorithms are
equivalent.  ``--bitwise`` remains available as a diagnostic on a fixed
platform.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np


_SHARD = re.compile(r"icesee_enkf_ens_(\d+)\.h5$")


@dataclass(frozen=True)
class StateLocation:
    path: Path
    dataset: str
    timestep: int
    time_axis: int | None = None


def _resolve_partial(source: Path) -> tuple[Path, h5py.Dataset]:
    path = source / "icesee_ensemble_data.h5" if source.is_dir() else source
    handle = h5py.File(path, "r")
    for name in ("ensemble", "states"):
        if name in handle:
            return path, handle[name]
    handle.close()
    raise KeyError(f"No ensemble/states dataset in {path}")


def _partial_locations(source: Path) -> tuple[list[StateLocation], tuple[int, int]]:
    path, dataset = _resolve_partial(source)
    shape = tuple(dataset.shape)
    dataset_name = dataset.name
    dataset.file.close()
    if len(shape) == 2:
        return [StateLocation(path, dataset_name, 0)], shape
    if len(shape) != 3:
        raise ValueError(f"Expected 2-D or 3-D partial ensemble, got {shape}")
    return (
        [StateLocation(path, dataset_name, k, 2) for k in range(shape[2])],
        shape[:2],
    )


def _full_locations(source: Path) -> tuple[list[StateLocation], tuple[int, int]]:
    directory = source if source.is_dir() else source.parent
    locations = []
    shape = None
    for path in directory.glob("icesee_enkf_ens_*.h5"):
        match = _SHARD.match(path.name)
        if not match:
            continue
        with h5py.File(path, "r") as handle:
            if "states" not in handle:
                continue
            this_shape = tuple(handle["states"].shape)
        if len(this_shape) != 2:
            raise ValueError(f"Expected 2-D state shard, got {this_shape} in {path}")
        if shape is None:
            shape = this_shape
        elif shape != this_shape:
            raise ValueError(f"Inconsistent mode-2 shard shapes: {shape} and {this_shape}")
        locations.append(StateLocation(path, "/states", int(match.group(1))))
    locations.sort(key=lambda item: item.timestep)
    if not locations or shape is None:
        raise FileNotFoundError(f"No mode-2 state shards found in {directory}")
    return locations, shape


def _read_rows(location: StateLocation, start: int, stop: int) -> np.ndarray:
    with h5py.File(location.path, "r") as handle:
        dataset = handle[location.dataset]
        if location.time_axis is None:
            return np.asarray(dataset[start:stop, :], dtype=np.float64)
        return np.asarray(dataset[start:stop, :, location.timestep], dtype=np.float64)


def _load_provenance(source: Path, expected_mode: int):
    directory = source if source.is_dir() else source.parent
    path = directory / "execution_parity_provenance.json"
    if not path.is_file():
        return None
    with path.open("r", encoding="utf-8") as stream:
        value = json.load(stream)
    if int(value.get("execution_mode", -1)) != expected_mode:
        raise ValueError(
            f"Invalid parity provenance in {path}: expected execution mode "
            f"{expected_mode}."
        )
    for name, expected_mtime in value.get("output_mtime_ns", {}).items():
        output = directory / name
        if not output.is_file() or output.stat().st_mtime_ns != expected_mtime:
            raise ValueError(
                f"Parity output {output} was removed or modified after its "
                "benchmark run."
            )
    return value


def _validate_provenance(partial: Path, full: Path, required: bool):
    left = _load_provenance(partial, 1)
    right = _load_provenance(full, 2)
    if required and (left is None or right is None):
        raise ValueError(
            "Fresh-run provenance is missing. Run both modes through "
            "run_execution_mode_parity.py; do not compare an older directory."
        )
    if (left is None) != (right is None):
        raise ValueError(
            "Only one mode has parity provenance; the directories do not "
            "belong to a trustworthy paired run."
        )
    if left is not None and left.get("run_id") != right.get("run_id"):
        raise ValueError(
            "Mode 1 and mode 2 have different parity run IDs; refusing a "
            "cross-run comparison."
        )
    if left is not None:
        if required and (
            int(left.get("schema", 0)) < 2
            or int(right.get("schema", 0)) < 2
        ):
            raise ValueError(
                "Parity provenance predates source-tree fingerprinting. "
                "Run both modes again with the current parity runner."
            )
        for key in (
            "source_config_sha256",
            "application_sha256",
            "git_head",
            "git_diff_sha256",
        ):
            if left.get(key) != right.get(key):
                raise ValueError(
                    f"Mode 1 and mode 2 provenance disagree on {key}; "
                    "refusing to compare outputs produced from different "
                    "inputs or source trees."
                )
        print(f"paired run id: {left['run_id']}")


def compare(partial, full, *, row_chunk, rtol, atol,
            require_provenance=True):
    partial = Path(partial)
    full = Path(full)
    _validate_provenance(partial, full, require_provenance)
    partial_locations, partial_shape = _partial_locations(partial)
    full_locations, full_shape = _full_locations(full)
    if partial_shape != full_shape:
        raise ValueError(f"Ensemble shape mismatch: mode1={partial_shape}, mode2={full_shape}")

    partial_by_t = {item.timestep: item for item in partial_locations}
    full_by_t = {item.timestep: item for item in full_locations}
    if set(partial_by_t) != set(full_by_t):
        missing_full = sorted(set(partial_by_t) - set(full_by_t))
        missing_partial = sorted(set(full_by_t) - set(partial_by_t))
        raise ValueError(
            "Timestep mismatch: "
            f"missing from mode2={missing_full}; missing from mode1={missing_partial}"
        )

    nd, nens = partial_shape
    global_max_abs = 0.0
    global_max_rel = 0.0
    mismatch_count = 0
    value_count = 0
    first_failure = None
    for timestep in sorted(partial_by_t):
        step_max = 0.0
        for start in range(0, nd, row_chunk):
            stop = min(nd, start + row_chunk)
            expected = _read_rows(partial_by_t[timestep], start, stop)
            actual = _read_rows(full_by_t[timestep], start, stop)
            delta = np.abs(actual - expected)
            tolerance = atol + rtol * np.abs(expected)
            failed = ~(delta <= tolerance)
            count = int(np.count_nonzero(failed))
            if count and first_failure is None:
                row, member = np.argwhere(failed)[0]
                first_failure = (
                    timestep, start + int(row), int(member),
                    float(expected[row, member]), float(actual[row, member]),
                )
            mismatch_count += count
            value_count += delta.size
            if delta.size:
                step_max = max(step_max, float(np.nanmax(delta)))
                global_max_abs = max(global_max_abs, float(np.nanmax(delta)))
                denominator = np.maximum(np.abs(expected), np.finfo(float).tiny)
                global_max_rel = max(
                    global_max_rel, float(np.nanmax(delta / denominator))
                )
        print(f"timestep {timestep:6d}: max_abs={step_max:.12g}")

    passed = mismatch_count == 0
    print("\nExecution-mode parity summary")
    print(f"  synchronized timesteps: {len(partial_by_t)}")
    print(f"  ensemble shape:         {nd} x {nens}")
    print(f"  compared values:        {value_count}")
    print(f"  mismatched values:      {mismatch_count}")
    print(f"  max absolute error:     {global_max_abs:.12g}")
    print(f"  max relative error:     {global_max_rel:.12g}")
    print(f"  tolerance:              atol={atol:g}, rtol={rtol:g}")
    if first_failure:
        t, row, member, expected, actual = first_failure
        print(
            "  first mismatch:         "
            f"t={t}, row={row}, member={member}, mode1={expected:.17g}, "
            f"mode2={actual:.17g}"
        )
    print(f"  result:                 {'PASS' if passed else 'FAIL'}")
    return passed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partial", required=True, help="Mode-1 run directory or HDF5 file")
    parser.add_argument("--full", required=True, help="Mode-2 run directory")
    parser.add_argument("--row-chunk", type=int, default=8192)
    parser.add_argument("--rtol", type=float, default=1.0e-8)
    parser.add_argument("--atol", type=float, default=1.0e-10)
    parser.add_argument(
        "--bitwise", action="store_true",
        help="Require exact equality (diagnostic only; not MPI-portable).",
    )
    parser.add_argument(
        "--allow-unprovenanced", action="store_true",
        help=(
            "Unsafe diagnostic override: permit old directories without a "
            "paired-run provenance record. Fresh paired provenance is required "
            "by default."
        ),
    )
    args = parser.parse_args()
    if args.row_chunk <= 0:
        parser.error("--row-chunk must be positive")
    try:
        rtol = 0.0 if args.bitwise else args.rtol
        atol = 0.0 if args.bitwise else args.atol
        passed = compare(
            args.partial, args.full, row_chunk=args.row_chunk,
            rtol=rtol, atol=atol,
            require_provenance=not args.allow_unprovenanced,
        )
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
