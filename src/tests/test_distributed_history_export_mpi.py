"""Exercise native history finalization across ensemble and spatial ranks."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher, mpi_test_command


@pytest.mark.parametrize('ranks,spatial', [(3, 1), (4, 2)])
def test_distributed_checkpoint_history_finalizes_and_refuses_missing_steps(tmp_path, ranks, spatial):
    launcher = find_compatible_mpi_launcher()
    if launcher is None:
        pytest.skip('compatible MPI launcher not available')
    root = Path(__file__).resolve().parents[2]
    worker = Path(__file__).parent / 'parallel_mpi/_history_export_worker.py'
    env = dict(os.environ, PYTHONPATH=f'{root.parent}:{root}')
    command = mpi_test_command(launcher, ranks) + [sys.executable, str(worker), str(tmp_path), str(spatial)]
    result = subprocess.run(command, env=env, text=True, capture_output=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads(next(line for line in result.stdout.splitlines() if line.startswith('{')))
    assert report == dict(world_size=ranks, spatial_ranks=spatial, exact=True, missing_history_refused=True)
