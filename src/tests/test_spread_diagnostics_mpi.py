"""Local numerical MPI checks; no scheduler or numerical-model solve."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher, mpi_test_command


@pytest.mark.parametrize('ranks,spatial', [(3, 1), (4, 2), (4, 4)])
def test_spread_and_diagnostic_ownership_in_real_mpi(ranks, spatial):
    launcher = find_compatible_mpi_launcher()
    if launcher is None:
        pytest.skip('compatible MPI launcher unavailable')
    worker = Path(__file__).parent/'parallel_mpi/_spread_diagnostics_worker.py'
    command = mpi_test_command(launcher, ranks)+[sys.executable, str(worker), str(spatial)]
    result = subprocess.run(command, env=os.environ.copy(), text=True, capture_output=True, timeout=60)
    assert result.returncode == 0, result.stdout+result.stderr
    data = json.loads(next(line for line in result.stdout.splitlines() if line.startswith('{')))
    assert data['ranks'] == ranks and data['spatial_ranks'] == spatial
    assert data['max_member_error'] < 3e-11
