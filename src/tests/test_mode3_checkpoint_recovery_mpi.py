# ==============================================================================
# @des: Real multi-rank crash consistency and resume for execution-mode-3
#       checkpoints (Level 2, see docs/testing.md). A 4-rank job (2 ensemble
#       groups x 2 spatial ranks) commits rolling checkpoints and then hits
#       ENOSPC on one non-root writer; a second job with a different process
#       grid resumes from the newest valid checkpoint and finishes the run.
# ==============================================================================
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

_REPO_ROOT = Path(__file__).resolve().parents[2]
_WORKER = Path(__file__).resolve().parent / "parallel_mpi" / "_mode3_checkpoint_recovery_worker.py"


def _launch(ranks, *args):
    mpirun = find_compatible_mpi_launcher()
    if mpirun is None:
        pytest.skip("no compatible mpirun available in this environment")
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{_REPO_ROOT.parent}:{_REPO_ROOT}"
    result = subprocess.run(
        [mpirun, "--oversubscribe", "-n", str(ranks), sys.executable, str(_WORKER), *args],
        capture_output=True, text=True, timeout=120, env=env,
    )
    assert result.returncode == 0, (
        f"worker exit {result.returncode}\nstdout:\n{result.stdout[-3000:]}\n"
        f"stderr:\n{result.stderr[-3000:]}"
    )
    assert "Traceback" not in result.stdout + result.stderr
    lines = [line for line in result.stdout.splitlines() if line.startswith("{")]
    assert len(lines) == 1, result.stdout
    return json.loads(lines[0])


def test_enospc_on_one_rank_keeps_previous_checkpoint_and_resume_continues(tmp_path):
    root = str(tmp_path / "_mode3_state_history")

    written = _launch(4, "write", root, "2")
    assert written["world_size"] == 4
    assert "No space left on device" in written["step6_error"]
    # Rolling retention kept the newest two; the failed step 6 left nothing
    # committed and step 5 is intact.
    assert written["committed"] == [4, 5]

    # Resume on a different process grid (1 group x 2 spatial ranks).
    resumed = _launch(2, "resume", root, "2")
    assert resumed["plan"] == [5, 6, 2]  # analyses at 2 and 4 already applied
    assert resumed["all_ranks_exact"] is True
    assert resumed["committed"] == [6, 7]
    assert resumed["analyses_completed"] == 3
