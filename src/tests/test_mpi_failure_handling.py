# ==============================================================================
# @des: Level-2 synthetic MPI test: proves the test harness can actually
#       detect a failed distributed run, and that ICESEE's MPI failure
#       pattern (print diagnostics, then comm.Abort() -- never a further
#       collective) prevents healthy ranks from hanging forever waiting for
#       a rank that raised and will never arrive.
#
#       This exists because a real crashed run was previously observed to
#       exit 0 (see icesee_da_full_parallel.py's exception-handling fix):
#       an unconditional comm_world.Barrier() plus collective cleanup
#       inside the exception handler could itself hang, and even when it
#       didn't, the exception was swallowed with no re-raise and no abort,
#       so the process exited as if nothing had happened. Every later
#       scalability test in this repo depends on "the harness reported
#       PASS" actually meaning the run succeeded; this test is the
#       regression guard for that assumption.
# ==============================================================================
import subprocess
import sys
from pathlib import Path

import pytest

from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

HERE = Path(__file__).resolve().parent
FAILURE_WORKER = HERE / "parallel_mpi" / "_intentional_failure_worker.py"
HEALTHY_WORKER = HERE / "parallel_mpi" / "_healthy_worker.py"

_LAUNCHER = find_compatible_mpi_launcher()

pytestmark = pytest.mark.skipif(
    _LAUNCHER is None, reason="mpirun/mpiexec not available in this environment"
)


def _run(worker_path, timeout=30):
    cmd = [_LAUNCHER, "-n", "2", sys.executable, str(worker_path)]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def test_healthy_two_rank_job_exits_zero():
    """Negative control: proves the harness distinguishes PASS from FAIL --
    a run with no failure must exit 0, not merely 'not hang'."""
    result = _run(HEALTHY_WORKER)
    assert result.returncode == 0, (
        f"healthy job unexpectedly failed:\nstdout={result.stdout}\n"
        f"stderr={result.stderr}"
    )
    assert "completed normally" in result.stdout


def test_one_rank_failure_terminates_with_nonzero_exit_and_visible_traceback():
    """The actual regression guard: one rank raising must not hang the
    job (bounded timeout below turns a hang into a test failure, not an
    indefinite wait), must produce a nonzero process result, and must not
    suppress the original exception's visibility."""
    try:
        result = _run(FAILURE_WORKER, timeout=30)
    except subprocess.TimeoutExpired as exc:
        pytest.fail(
            "one rank raising an exception hung the whole MPI job instead "
            f"of terminating it (timed out after 30s). Partial output: "
            f"{exc.stdout!r} / {exc.stderr!r}"
        )

    combined = result.stdout + result.stderr
    assert result.returncode != 0, (
        "a rank raising an exception must produce a non-zero job result "
        f"(this is exactly the false-success failure mode this test "
        f"guards against); got rc=0.\nstdout={result.stdout}\n"
        f"stderr={result.stderr}"
    )
    assert "intentional test failure on rank 1" in combined, (
        "the original exception's message must remain visible, not be "
        f"suppressed:\n{combined}"
    )
    assert "Fatal error" in combined
    # The defining property under test: rank 0 (healthy) must never reach
    # its Barrier() -- if it did, Abort() failed to prevent the hang this
    # test exists to catch.
    assert "reached the barrier unexpectedly" not in combined
