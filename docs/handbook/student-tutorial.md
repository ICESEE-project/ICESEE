# First reproducible ICESEE run

This tutorial uses the small application stored in `applications/lorenz_model/examples/lorenz96`. Despite that directory's historical name, its present three-component equations are the Lorenz-63 equations:

```text
dx/dt = sigma*(y-x)
dy/dt = x*(rho-z)-y
dz/dt = x*y-beta*z
```

The example advances them with classical fourth-order Runge–Kutta. It is useful for learning orchestration, observations, ensembles and output without a glacier dataset. It does not validate distributed Firedrake meshes or glacier physics.

## Reading the execution modes

Before the parallel exercise, read [parallelization modes and the assimilation cycle](parallelization-modes.md). Follow the Mode-0 diagram while running the serial case, then compare Mode 1's root transform and state-row scatter with Mode 2's file-backed row reductions. The worked `P=8`, `Ne=12` examples show why member IDs, MPI ranks, and ranks per model are different quantities. Mode 3 is registered per application and needs its own model/environment procedure; it is not a general replacement for the other commands in this tutorial.

## 1. Establish source and imports

Use the [deployment guide](deployment.md) to activate a compatible Python/MPI environment. In the ICESEE checkout root, record the commit and prepare imports. The checkout directory should be named `ICESEE` for this simple setup:

```bash
export ICESEE_ROOT="$(pwd -P)"
export PYTHONPATH="$(dirname "$ICESEE_ROOT"):$ICESEE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
git rev-parse HEAD
git status --short
python -c 'import ICESEE, sys; print(sys.executable); print(ICESEE.__file__)'
```

The inspected implementation baseline is `3748e416d3a8120cc860119bd7224c6f9ad41ff7` on `features_debug`. The new teaching files and ISSM corrections are included in the immutable snapshot selected by [deployment](deployment.md#2-select-the-scientific-source-explicitly). Do not replace an active research checkout's branch merely to follow the tutorial. Use a clean dedicated checkout or record your actual source and treat differences as part of the experiment.

Create a new absolute working location for this exercise, outside any existing result tree:

```bash
CASE_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/icesee-student.XXXXXX")"
export CASE_ROOT
cp "$ICESEE_ROOT/docs/handbook/examples/lorenz-smoke.yaml" "$CASE_ROOT/params.yaml"
git rev-parse HEAD > "$CASE_ROOT/source.txt"
git diff --binary HEAD > "$CASE_ROOT/source.diff"
```

Temporary storage can be purged; copy successful learning records to your group's approved research storage afterwards. `data_path` is automatically cleaned for a fresh run, so every command below uses a new mode-specific subdirectory.

## 2. Read the small configuration

The supplied fully resolved [smoke YAML](examples/lorenz-smoke.yaml) uses one model-time unit, `dt=0.02`, 50 forecast cycles, eight members and observations at interior times `0.2`, `0.4`, `0.6`, `0.8`. There are 51 saved state columns, including the initial condition. Initial, process and observation uncertainty scales are specified per component. The analysis-factor mode is explicitly `legacy_prior_anomalies` to match the historical example; it is a forecast-dependent factor, not independent Gaussian `R` noise.

The driver uses a prior reference `[2,3,4]` and truth initial state `[1,1,1]`. These are synthetic trajectories used to teach assimilation, not external measurements. State rows are ordered `x,y,z`. `state_distribution: replicated` records that a model group holds copies of this tiny complete state rather than spatial partitions.

## 3. Run the serial baseline

```bash
python "$ICESEE_ROOT/applications/lorenz_model/examples/lorenz96/run_da_lorenz96.py" \
  -F "$CASE_ROOT/params.yaml" --execution_mode=0 \
  --data_path="$CASE_ROOT/mode0" \
  > "$CASE_ROOT/mode0.stdout.log" 2> "$CASE_ROOT/mode0.stderr.log"
printf '%s\n' "$?" > "$CASE_ROOT/mode0.exit-status"
```

Read stdout and stderr before continuing. An exit code is one check, not the whole acceptance test. The log should identify the chosen runner and reach normal completion. The first saved state represents time zero; no extra forecast should be needed just to populate column zero.

## 4. Inspect state and analysis timing

Use ordinary HDF5 reads on the closed file:

```bash
python - <<'PY'
import os
from pathlib import Path
import h5py
import numpy as np
root = Path(os.environ['CASE_ROOT']) / 'mode0'
with h5py.File(root / 'icesee_ensemble_data.h5', 'r') as f:
    print('datasets:', list(f))
    mean = f['ensemble_mean'][:]
    assert mean.shape == (3, 51), mean.shape
    assert np.isfinite(mean).all()
    if 'ensemble' in f:
        members = f['ensemble'][:]
        assert members.shape == (3, 8, 51), members.shape
        assert np.isfinite(members).all()
        np.testing.assert_allclose(mean, members.mean(axis=1), rtol=1e-10, atol=1e-10)
print('Shape and finite-value checks passed')
PY
```

Inspect the saved truth/observation files and their metadata too. Identify the observation state columns corresponding to times `0.2`, `0.4`, `0.6`, `0.8`. In runs emitting `[ICESEE_ANALYSIS]`, each record's `state_column` and `observation_state_column` must match and `model_time` must match the mapped observation time. The exact fields written by the runner, rather than a notebook's assumptions, determine how to interpret the file.

## 5. Run actual parallel modes

Before launching, check the MPI pair and parallel HDF5:

```bash
python -c 'import h5py; print(h5py.__version__, "parallel:", h5py.get_config().mpi)'
python - <<'PY'
from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher
launcher = find_compatible_mpi_launcher()
assert launcher is not None
print('Use this launcher:', launcher)
PY
```

Set `MPIEXEC` to the compatible executable printed above. On an HPC allocation use the scheduler's approved launch form and stay within the allocated rank count. For a local approved two-rank launch:

```bash
"$MPIEXEC" -np 2 python \
  "$ICESEE_ROOT/applications/lorenz_model/examples/lorenz96/run_da_lorenz96.py" \
  -F "$CASE_ROOT/params.yaml" --execution_mode=2 \
  --data_path="$CASE_ROOT/mode2" \
  > "$CASE_ROOT/mode2.stdout.log" 2> "$CASE_ROOT/mode2.stderr.log"
```

The active MPI communicator must have size two; two independent processes each reporting size one are a failed launch. Inspect `time_grid_contract.json`, timestep shards, ensemble mean and final view. The expected mean shape remains `(3,51)`; MPI rank count must not multiply this replicated model's state dimension.

Use the [paired parity harness](testing-and-contributing.md#modes-1-and-2-parity) to compare Modes 1 and 2 from the same YAML into a fresh output root. Its member-wise result is more meaningful than plotting two similar mean curves. Run Mode 3 separately only after reading its application's support limits; a Mode-3 comparison can differ in stochastic initialization/forecast policy and should not be folded into a same-input parity claim without checking those realizations.

## 6. Learning exercises

1. Locate the forecast hook and RK4 routine. Explain which code owns a model step and which code decides when observations are due.
2. Change `Nens` in a copied config, rerun into a new directory, and predict which output dimension changes. Then explain why increasing MPI ranks alone should not change the scientific state dimension.
3. Set `process_noise_schedule: none` in a new case and compare spread/innovations. Record this as a scientific configuration change, even if the job runs faster.
4. Select `stochastic_R` in another copy and explain how `Eta` and innovations differ from the legacy factor. Compare members as well as means.
5. On a disposable run, remove a required shard after preserving the original and verify that an integrity check rejects the incomplete output. Restore the original rather than treating an attractive plot as proof of completeness.
6. With a trained operator, perform the target-model restart exercise described in [operations](operations.md). A serial Lorenz smoke run alone is insufficient preparation for recovering a PIG production job.

For each exercise retain config, source identity, command, exit/logs, output checks and your interpretation. A student has completed onboarding when another student can reproduce that record and explain both the numerical choices and their limitations.
