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


@pytest.mark.parametrize('ranks,spatial', [(3, 1), (4, 2)])
def test_serial_finalizer_recovers_after_actual_mpi_abort(tmp_path, ranks, spatial):
    import h5py
    import numpy as np

    launcher = find_compatible_mpi_launcher()
    if launcher is None:
        pytest.skip('compatible MPI launcher not available')
    root = Path(__file__).resolve().parents[2]
    worker = Path(__file__).parent / 'parallel_mpi/_history_export_worker.py'
    env = dict(os.environ, PYTHONPATH=f'{root.parent}:{root}')
    command = mpi_test_command(launcher, ranks) + [sys.executable, str(worker), str(tmp_path),
                                                  str(spatial), '--abort-after-prefix']
    failed = subprocess.run(command, env=env, text=True, capture_output=True, timeout=90)
    assert failed.returncode != 0, failed.stdout + failed.stderr
    exporter = root / 'src/parallelization/distributed_history_export.py'
    recovered = subprocess.run([sys.executable, str(exporter), str(tmp_path), '--allow-partial',
                                '--simulation-exit-code', str(failed.returncode)],
                               env=env, text=True, capture_output=True, timeout=30)
    assert recovered.returncode == 0, recovered.stdout + recovered.stderr
    rows = np.r_[np.arange(6, -1, -1), 7 + np.arange(2, -1, -1)]
    expected = 100 * rows[:, None, None] + 10 * np.arange(5)[None, :, None] + np.arange(3)[None, None, :]
    with h5py.File(tmp_path / 'icesee_ensemble_data.h5', 'r') as handle:
        assert handle.attrs['export_complete'] and not handle.attrs['complete']
        assert handle.attrs['requested_steps'] == 3
        assert handle.attrs['completed_steps'] == 2
        assert handle.attrs['simulation_exit_code'] == failed.returncode
        np.testing.assert_array_equal(handle['ensemble'][:], expected)
        np.testing.assert_array_equal(handle['t'][:], [0, .05, .1])
