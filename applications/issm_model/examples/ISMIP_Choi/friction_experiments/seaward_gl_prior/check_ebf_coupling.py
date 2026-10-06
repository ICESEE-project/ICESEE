#!/usr/bin/env python3
"""Fail unless an EBF result changes both friction and the forecast state."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

EXAMPLE_ROOT = Path(__file__).resolve().parents[2]


def read_mean(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as handle:
        return handle["ensemble_mean"][:]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--wbf",
        type=Path,
        default=EXAMPLE_ROOT
        / (
            "_modelrun_datasets_wbf_2_spinup_dt01_smoke/"
            "icesee_ensemble_data.h5"
        ),
    )
    parser.add_argument(
        "--ebf",
        type=Path,
        default=EXAMPLE_ROOT
        / (
            "_modelrun_datasets_ebf_2_spinup_dt01_smoke/"
            "icesee_ensemble_data.h5"
        ),
    )
    parser.add_argument("--tolerance", type=float, default=1.0e-8)
    parser.add_argument(
        "--initial-tolerance",
        type=float,
        default=1.0e-8,
        help="Maximum permitted WBF/EBF difference at timestep zero.",
    )
    args = parser.parse_args()

    wbf = read_mean(args.wbf)
    ebf = read_mean(args.ebf)
    if wbf.shape[0] != ebf.shape[0] or wbf.shape[0] % 6:
        raise ValueError("WBF and EBF state-vector dimensions are incompatible")

    nt = min(wbf.shape[1], ebf.shape[1])
    nvertices = wbf.shape[0] // 6
    state_stop = 5 * nvertices
    coefficient_start = state_stop

    state_delta = np.abs(ebf[:state_stop, :nt] - wbf[:state_stop, :nt])
    coefficient_delta = np.abs(
        ebf[coefficient_start:, :nt] - wbf[coefficient_start:, :nt]
    )
    max_state_delta = float(np.max(state_delta))
    max_coefficient_delta = float(np.max(coefficient_delta))
    initial_state_delta = float(np.max(state_delta[:, 0]))
    initial_coefficient_delta = float(np.max(coefficient_delta[:, 0]))

    changed_steps = np.flatnonzero(
        np.max(coefficient_delta, axis=0) > args.tolerance
    )
    first_changed = int(changed_steps[0]) if changed_steps.size else None
    coefficient_difference = (
        ebf[coefficient_start:, :nt] - wbf[coefficient_start:, :nt]
    )
    update_steps = np.flatnonzero(
        np.max(
            np.abs(np.diff(coefficient_difference, axis=1)),
            axis=0,
        )
        > args.tolerance
    ) + 1
    post_update_deltas = {
        int(step): float(np.max(state_delta[:, step + 1]))
        for step in update_steps
        if step + 1 < nt
    }
    propagating_steps = [
        step
        for step, delta in post_update_deltas.items()
        if delta > args.tolerance
    ]
    max_post_update_state_delta = max(post_update_deltas.values(), default=0.0)

    print(f"compared timesteps: {nt}")
    print(f"initial dynamic-state difference: {initial_state_delta:.6g}")
    print(f"initial coefficient difference: {initial_coefficient_delta:.6g}")
    print(f"first coefficient difference: {first_changed}")
    print(f"max |EBF-WBF coefficient|: {max_coefficient_delta:.6g}")
    print(f"max |EBF-WBF dynamic state|: {max_state_delta:.6g}")
    print(f"coefficient-update steps: {update_steps.tolist()}")
    print(
        "max dynamic difference one forecast after an update: "
        f"{max_post_update_state_delta:.6g}"
    )
    print(
        "first update with forecast propagation: "
        f"{propagating_steps[0] if propagating_steps else None}"
    )

    if (
        initial_state_delta > args.initial_tolerance
        or initial_coefficient_delta > args.initial_tolerance
    ):
        raise SystemExit(
            "FAIL: WBF and EBF did not start from the same ensemble"
        )
    if max_coefficient_delta <= args.tolerance:
        raise SystemExit("FAIL: EBF did not update friction")
    if not post_update_deltas:
        raise SystemExit(
            "FAIL: result ends before coefficient-to-forecast coupling can be checked"
        )
    if not propagating_steps:
        raise SystemExit(
            "FAIL: friction changed but every following forecast stayed "
            "identical to WBF"
        )

    print("PASS: EBF friction updates propagate into the dynamic state")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
