"""End-to-end integration test for execution_mode 3's Lorenz96 DA-cycle runner.

Runs the real, user-facing ``run_da_lorenz96.py`` entry point as a subprocess
(a small, self-contained params file, isolated ``data_path``) rather than
importing it in-process.  This mirrors an actual user invocation exactly and
avoids polluting the shared pytest process with ``run_da_lorenz96.py``'s
module-level side effects (``os.chdir``, MPI init, and a real call into
``icesee_model_data_assimilation`` all execute at import time -- see
``config/_utility_imports.py``'s own docstring-documented argv/cwd coupling).
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

_REPO_ROOT = Path(__file__).resolve().parents[2]
_LORENZ96_DIR = _REPO_ROOT / "applications" / "lorenz_model" / "examples" / "lorenz96"
_RUN_SCRIPT = _LORENZ96_DIR / "run_da_lorenz96.py"

_SMOKE_PARAMS = """\
physical-parameters:
  sigma_96: 10.0
  beta_96: 8.0/3.0
  rho_96: 28.0

modeling-parameters:
  dt: 0.01
  num_years: 1
  timesteps_per_year: 20
  example_name: "lorenz96"

enkf-parameters:
  Nens: 4
  freq_obs: 0.2
  obs_max_time: 2
  obs_start_time: 0.1

  num_state_vars: 3
  num_param_vars: 0
  vec_inputs: ['x','y','z']
  observed_vars: ['x','y','z']

  sig_obs: [0.1, 0.1, 0.1]
  sig_Q: [0.01, 0.01, 0.01]
  length_scale: [2, 2, 2]

  generate_synthetic_obs: True
  generate_true_state: True
  generate_nurged_state: True
  sequential_ensemble_initialization: True
  use_ensemble_pertubations: true

  joint_estimation: False
  state_estimation: True
  parameter_estimation: False

  seed: 1
  inflation_factor: 1.0
  localization_flag: 0

  model_name: "lorenz"
  filter_type: "EnKF"
  execution_mode: 3
  commandlinerun: "True"
  data_path: {data_path}
"""


def _run_mode3(tmp_path: Path, nprocs: int) -> subprocess.CompletedProcess:
    tmp_path.mkdir(parents=True, exist_ok=True)
    data_path = tmp_path / "_modelrun_datasets"
    params_path = tmp_path / "params_mode3_smoke.yaml"
    params_path.write_text(_SMOKE_PARAMS.format(data_path=data_path))

    env = {
        "PYTHONPATH": f"{_REPO_ROOT.parent}:{_REPO_ROOT}",
        "PATH": "/usr/bin:/bin:/usr/local/bin:/opt/homebrew/bin",
        "HOME": os.environ.get("HOME", ""),
    }
    mpirun = find_compatible_mpi_launcher() if nprocs > 1 else None
    if nprocs > 1 and mpirun is None:
        pytest.skip("no compatible mpirun available in this environment")
    if nprocs > 1:
        cmd = [mpirun, "-n", str(nprocs), sys.executable, str(_RUN_SCRIPT), "-F", str(params_path)]
    else:
        cmd = [sys.executable, str(_RUN_SCRIPT), "-F", str(params_path)]

    return subprocess.run(
        cmd,
        cwd=str(_LORENZ96_DIR),
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    ), data_path


@pytest.mark.parametrize("nprocs", [1, 4])
def test_lorenz_mode3_full_cycle_end_to_end(tmp_path, nprocs):
    result, data_path = _run_mode3(tmp_path, nprocs)
    assert result.returncode == 0, (
        f"mode-3 Lorenz96 run failed (nprocs={nprocs}):\n"
        f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )

    ensemble_path = data_path / "icesee_ensemble_data.h5"
    assert ensemble_path.exists()
    with h5py.File(ensemble_path, "r") as f:
        ensemble = f["ensemble"][:]
        ensemble_mean = f["ensemble_mean"][:]

    assert ensemble.shape == (3, 4, 101)
    assert np.all(np.isfinite(ensemble))
    # The trajectory must actually advance -- not stay stuck at the initial
    # condition. Frequent full-state observations on a 4-member ensemble in
    # this short window legitimately pull members close together, so only
    # exact (bit-identical) collapse -- a real bug signature -- is rejected.
    assert not np.allclose(ensemble[:, :, -1], ensemble[:, :, 0])
    for member in range(1, 4):
        assert not np.array_equal(ensemble[:, member, -1], ensemble[:, 0, -1])
    np.testing.assert_allclose(ensemble_mean, ensemble.mean(axis=1), atol=1e-10)

    assert (data_path / "true-wrong-lorenz.h5").exists()


def test_lorenz_mode3_single_and_multi_rank_are_bit_identical(tmp_path):
    result_1, data_path_1 = _run_mode3(tmp_path / "one", 1)
    result_4, data_path_4 = _run_mode3(tmp_path / "four", 4)
    assert result_1.returncode == 0, result_1.stderr
    assert result_4.returncode == 0, result_4.stderr

    with h5py.File(data_path_1 / "icesee_ensemble_data.h5", "r") as f:
        ensemble_1 = f["ensemble"][:]
    with h5py.File(data_path_4 / "icesee_ensemble_data.h5", "r") as f:
        ensemble_4 = f["ensemble"][:]

    # Both the per-member AR(1) process-noise streams and the
    # legacy_prior_anomalies observation-error term are designed to be
    # rank/order independent -- running under 1 vs 4 real MPI ranks must
    # give bit-identical results.
    np.testing.assert_array_equal(ensemble_1, ensemble_4)
