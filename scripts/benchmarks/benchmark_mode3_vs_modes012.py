#!/usr/bin/env python3
"""Benchmark real-MPI execution-mode-3 primitives against production modes 0-2.

## Why this benchmark compares different things on each side

Mode 3 has no per-application production DA-cycle runner yet. Modes 0, 1, and
2 each implement the complete workflow -- true/nudged state generation,
synthetic-observation generation, ensemble initialization, the forecast-
analysis loop, and output writing -- in `icesee_da_serial.py`,
`icesee_da_partial_parallel.py`, and `icesee_da_full_parallel.py`
respectively. Mode 3's own reference adapter
(`applications/lorenz_model/lorenz_utils/distributed_adapter.py`) is
documented in its own module docstring as "used only for state-only parity
development" and guards against large states with an explicit allgather size
cap -- it is not a production adapter and no production mode-3 runner has been
built on top of it. See `docs/execution-mode-3-design.md`'s phased tracker
and its "Mode 3 remains non-selectable until ... gates pass" acceptance
statement; this script does not change that.

This script therefore runs, side by side, for one matching Lorenz96 case:

- modes 0, 1, 2: the real, full production DA cycle, launched exactly as
  `scripts/ci/run_lorenz96_ci.py` does, through the actual
  `run_da_lorenz96.py` entry point (true state, observations, ensemble init,
  the full forecast/analysis loop, and HDF5 output); and
- "mode 3": the real-MPI, already parity-tested forecast + one stochastic-
  analysis-step primitives used to gate mode-3 development
  (`scripts/benchmarks/run_mode3_lorenz_forecast_parity.py`), run with a
  matching ensemble size, step count, and dt.

This is a wall-clock and correctness benchmark, not a byte-for-byte DA-cycle
parity comparison: the mode-3 side's one synthetic analysis step does not
exercise the same observation schedule, ensemble initialization, or output
format as modes 0-2. Its purpose is to inform whether building a full mode-3
production runner is worth its engineering cost, before that work starts.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
LORENZ_DIR = REPO_ROOT / "applications" / "lorenz_model" / "examples" / "lorenz96"
DATA_DIR = LORENZ_DIR / "_modelrun_datasets"
MODE3_SCRIPT = REPO_ROOT / "scripts" / "benchmarks" / "run_mode3_lorenz_forecast_parity.py"


def _mode_config(mode: int, temporary_directory: Path) -> tuple[Path, Path]:
    source = LORENZ_DIR / "params.yaml"
    with source.open("r", encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    output = DATA_DIR / f"benchmark_mode3_vs_mode_{mode}"
    enkf = config["enkf-parameters"]
    enkf["execution_mode"] = mode
    enkf["data_path"] = str(output)
    enkf["restart_enabled"] = False
    enkf["force_fresh_start"] = True
    destination = temporary_directory / f"lorenz_mode_{mode}.yaml"
    with destination.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(config, stream, sort_keys=False)
    return destination, output, config


def _run_timed(command: list[str], cwd: Path) -> float:
    print(f"+ {' '.join(command)}", flush=True)
    environment = os.environ.copy()
    # ``run_da_lorenz96.py`` is imported as ``ICESEE.applications....``; the
    # package root is this repository's parent directory.
    repo_parent = str(REPO_ROOT.parent)
    existing = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        repo_parent if not existing else f"{repo_parent}{os.pathsep}{existing}"
    )
    started = time.time()
    subprocess.run(command, cwd=cwd, check=True, env=environment)
    return time.time() - started


def _run_mode_0_1_2(mode: int, python: str, mpiexec: str, ranks: int,
                     temporary_directory: Path) -> tuple[float, dict]:
    config, output, source_config = _mode_config(mode, temporary_directory)
    if output.exists():
        shutil.rmtree(output)
    command = [
        python, "-m",
        "ICESEE.applications.lorenz_model.examples.lorenz96.run_da_lorenz96",
        "-F", str(config),
    ]
    if mode in (1, 2):
        command = [mpiexec, "-np", str(ranks)] + command
    elapsed = _run_timed(command, LORENZ_DIR)
    return elapsed, source_config


def _run_mode_3(python: str, mpiexec: str, ensemble_groups: int, spatial_ranks: int,
                 nens: int, steps: int, dt: float, analysis_step: int,
                 atol: float, rtol: float) -> float:
    ranks = ensemble_groups * spatial_ranks
    command = [
        mpiexec, "-np", str(ranks), python, str(MODE3_SCRIPT),
        "--ensemble-groups", str(ensemble_groups),
        "--nens", str(nens),
        "--steps", str(steps),
        "--dt", str(dt),
        "--analysis-step", str(analysis_step),
        "--atol", str(atol),
        "--rtol", str(rtol),
    ]
    return _run_timed(command, REPO_ROOT)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", type=int, default=2,
                        help="MPI ranks for modes 1 and 2, and for mode 3's "
                             "ensemble-groups x spatial-ranks grid")
    parser.add_argument("--spatial-ranks", type=int, default=1,
                        help="Mode-3 spatial ranks per member (Lorenz's "
                             "three-variable state has no useful spatial "
                             "decomposition, so 1 is the default)")
    parser.add_argument("--analysis-step", type=int, default=20)
    parser.add_argument(
        "--mode3-atol", type=float, default=1.0e-8,
        help="Mode-3 primitive gate tolerance. The script's own default "
             "(1e-13) is tuned for its short 20-step demonstration case; a "
             "1000-step case like this benchmark's default accumulates "
             "round-off from per-rank reduction-order differences well past "
             "that, exactly as documented for modes 1 vs 2 in "
             "compare_execution_modes.py.",
    )
    parser.add_argument("--mode3-rtol", type=float, default=1.0e-8)
    parser.add_argument("--mpiexec", default=None)
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()

    if args.ranks < 1 or args.spatial_ranks < 1:
        parser.error("--ranks and --spatial-ranks must be positive")
    if args.ranks % args.spatial_ranks:
        parser.error("--ranks must be divisible by --spatial-ranks")
    ensemble_groups = args.ranks // args.spatial_ranks

    mpiexec = args.mpiexec or shutil.which("mpirun") or shutil.which("mpiexec")
    if mpiexec is None:
        raise RuntimeError("An MPI launcher is required for modes 1, 2, and 3")

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}
    with tempfile.TemporaryDirectory(prefix="icesee-mode3-benchmark-") as temporary:
        temporary_directory = Path(temporary)
        source_config = None
        for mode in (0, 1, 2):
            elapsed, source_config = _run_mode_0_1_2(
                mode, args.python, mpiexec, args.ranks, temporary_directory,
            )
            timings[f"mode {mode}"] = elapsed

        modeling = source_config["modeling-parameters"]
        enkf = source_config["enkf-parameters"]
        nens = int(enkf["Nens"])
        dt = float(modeling["dt"])
        steps = int(round(float(modeling["num_years"]) / dt))
        elapsed = _run_mode_3(
            args.python, mpiexec, ensemble_groups, args.spatial_ranks,
            nens, steps, dt, args.analysis_step,
            args.mode3_atol, args.mode3_rtol,
        )
        timings[
            f"mode 3 (reference primitives, {ensemble_groups}x{args.spatial_ranks} grid)"
        ] = elapsed

    print("\nWall-clock benchmark: Lorenz96, Nens="
          f"{nens}, steps={steps}, dt={dt}")
    for label, elapsed in timings.items():
        print(f"  {label:55s} {elapsed:8.3f} s")
    print(
        "\nNote: mode 3's timing covers real-MPI forecast + one stochastic-"
        "analysis step only (no true/observation generation, ensemble init, "
        "or output writing); see this script's module docstring."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
