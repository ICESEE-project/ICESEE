# ==============================================================================
# @des: Regression test for the Idealized PIG model_nprocs>1 startup deadlock
# (PACE, 2026-10-03). PyOP2's make_so looks the compiled kernel up in the disk
# cache on each rank independently but compiles collectively, so a model
# group whose ranks see different cache contents hangs forever. ICESEE's
# install_consistent_compile_cache() makes the hit/miss decision collective.
# ==============================================================================
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

pytest.importorskip("firedrake")

_REPO_ROOT = Path(__file__).resolve().parents[2]
_WORKER = Path(__file__).resolve().parent / "parallel_mpi" / "_pyop2_uneven_cache_worker.py"


def _run(tmp_path, fix, timeout):
    mpirun = find_compatible_mpi_launcher()
    if mpirun is None:
        pytest.skip("no compatible mpirun available in this environment")
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{_REPO_ROOT.parent}:{_REPO_ROOT}"
    # A fresh salt per call gives kernels no earlier run has cached.
    salt = f"{time.time_ns() % 100000 / 7.0:.6f}"
    try:
        return subprocess.run(
            [mpirun, "--oversubscribe", "-n", "2", sys.executable, str(_WORKER),
             str(int(fix)), str(tmp_path), salt],
            capture_output=True, text=True, timeout=timeout, env=env,
        )
    except subprocess.TimeoutExpired:
        return None


def test_uneven_private_caches_do_not_deadlock_with_consistent_cache(tmp_path):
    result = _run(tmp_path, fix=True, timeout=300)
    assert result is not None, "2-rank compile deadlocked despite the consistent cache"
    ok = [line for line in result.stdout.splitlines() if line.startswith("OK")]
    assert result.returncode == 0 and len(ok) == 2, result.stdout + result.stderr
    assert len({line.split("value=")[1] for line in ok}) == 1


@pytest.mark.skipif(
    not os.environ.get("ICESEE_RUN_DEADLOCK_REPRO"),
    reason="demonstrates the unpatched PyOP2 deadlock (slow, ends by timeout); "
    "set ICESEE_RUN_DEADLOCK_REPRO=1 to run",
)
def test_unpatched_pyop2_deadlocks_on_uneven_private_caches(tmp_path):
    assert _run(tmp_path, fix=False, timeout=180) is None
