# ==============================================================================
# @des: Focused regression test for the true-vs-wrong basal-melt forcing
# experiment (Idealized PIG, 2026-09-28 reconciliation with icesee/features).
#
# Verifies, through the ACTUAL BasalMeltRate() code path (never a duplicated
# formula in this test):
#   - true trajectory (experiment=EXPERIMENT_TRUE): melt_max ramps from 20
#     at the start of the run to ~100 by the end (linear ramp over
#     num_years, per _icepack_model.BasalMeltRate).
#   - wrong/ensemble trajectory (experiment=EXPERIMENT_WRONG): melt_max
#     stays constant at 20 for the entire run.
#   - BasalMeltRate rejects any experiment value other than the two
#     canonical strings (no silent str/bool fallback).
#
# Runs the real computation in a subprocess under `mpirun -n 1` (this
# suite's established convention for anything that touches real Firedrake
# objects -- see src/tests/parallel_mpi/_basal_melt_experiment_worker.py),
# not via direct in-process Firedrake import.
# ==============================================================================
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

_REPO_ROOT = Path(__file__).resolve().parents[2]
_WORKER = Path(__file__).resolve().parent / "parallel_mpi" / "_basal_melt_experiment_worker.py"
# config/_utility_imports.py's module-level parameter loading looks up
# 'params.yaml' relative to the process CWD by default -- run the worker
# from Idealized PIG's own example directory, matching how run_da_icepack.py
# is normally invoked.
_IDEALIZED_PIG_DIR = _REPO_ROOT / "applications" / "icepack_model" / "examples" / "idealized_pig"


def _run_worker() -> dict:
    mpirun = find_compatible_mpi_launcher()
    if mpirun is None:
        pytest.skip("no compatible mpirun available in this environment")
    try:
        import firedrake  # noqa: F401
    except ImportError:
        pytest.skip("firedrake is not importable in this environment")

    env = dict(os.environ)
    env["PYTHONPATH"] = f"{_REPO_ROOT.parent}:{_REPO_ROOT}"
    # A fresh (uncached) PyOP2 JIT compile can pick up an unrelated
    # ISSM-specific mpicc via a stale PATH entry instead of the one
    # matching this Firedrake install's actual Open MPI -- scope the fix
    # to this subprocess only.
    env["PATH"] = f"/opt/homebrew/bin:{env.get('PATH', '')}"

    # The config loader cleans data_path at import time; without an isolated
    # data_path the worker would wipe idealized_pig/_modelrun_datasets.
    data_path = tempfile.mkdtemp(prefix="icesee_basal_melt_test_")
    result = subprocess.run(
        [mpirun, "--oversubscribe", "-n", "1", sys.executable, str(_WORKER),
         f"--data_path={data_path}"],
        capture_output=True, text=True, timeout=180, env=env,
        cwd=str(_IDEALIZED_PIG_DIR),
    )
    assert result.returncode == 0, (
        f"worker non-zero exit ({result.returncode}).\n"
        f"stdout tail:\n{result.stdout[-4000:]}\nstderr tail:\n{result.stderr[-2000:]}"
    )
    # Worker prints exactly one JSON line; be tolerant of any Firedrake/PETSc
    # banner lines that may precede it.
    for line in reversed(result.stdout.strip().splitlines()):
        line = line.strip()
        if line.startswith("{"):
            return json.loads(line)
    raise AssertionError(f"worker produced no JSON line.\nstdout:\n{result.stdout}")


def test_true_trajectory_ramps_20_to_100():
    results = _run_worker()
    melt_max_first, melt_max_last = results["true"]
    assert melt_max_first == pytest.approx(20.0, abs=1e-6)
    assert melt_max_last == pytest.approx(100.0, abs=0.1)


def test_wrong_trajectory_stays_constant_at_20():
    results = _run_worker()
    melt_max_first, melt_max_last = results["wrong"]
    assert melt_max_first == pytest.approx(20.0, abs=1e-6)
    assert melt_max_last == pytest.approx(20.0, abs=1e-6)


def test_true_and_wrong_trajectories_diverge_by_the_end():
    results = _run_worker()
    _true_first, true_last = results["true"]
    _wrong_first, wrong_last = results["wrong"]
    assert true_last - wrong_last > 50.0


def test_basal_melt_rate_rejects_unrecognized_experiment_value():
    """Guards against the exact fragility this reconciliation removed: a
    caller passing a Python bool (or any value other than the two
    canonical strings) must fail loudly, never silently fall into the
    wrong-trajectory branch via a str/bool type mismatch. Exercised inside
    the worker subprocess (real BasalMeltRate/Firedrake objects), not via
    a direct in-process firedrake import."""

    results = _run_worker()
    assert results["invalid_experiment_raises"] is True
