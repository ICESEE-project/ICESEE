# ISSM Model with ICESEE Data Assimilation

> Coupled Ice Sheet System Model (ISSM) + ICESEE ensemble data assimilation, demonstrated through the ISMIP example. Supports serial and parallel (MPI) execution with multiple EnKF variants, and containerized HPC deployment via Apptainer.

---

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
- [Running the Model](#running-the-model)
  - [Serial Execution](#serial-execution)
  - [Parallel Execution (Local)](#parallel-execution-local)
  - [HPC Deployment with Apptainer](#hpc-deployment-with-apptainer)
- [Bind Mounts Reference](#bind-mounts-reference)
- [Outputs](#outputs)
- [Directory Structure](#directory-structure)
- [Troubleshooting](#troubleshooting)

---

## Overview

The `run_da_issm.py` script couples ISSM with ICESEE to perform data assimilation using ensemble Kalman filter (EnKF) variants. Each ensemble member runs an independent MATLAB-backed ISSM model instance, coordinated via MPI for parallel execution.

Supported filter types (configured in `params.yaml`):

- `EnKF` — standard ensemble Kalman filter
- Additional variants configurable via `filter_type`

---

## Prerequisites

| Requirement | Version / Notes |
|-------------|-----------------|
| Python | 3.8 or higher |
| MPI | OpenMPI or MPICH |
| MATLAB | Required for ISSM execution; must be accessible on `$PATH` |
| Python packages | See `requirements.txt` |

---

## Configuration

Model parameters are defined in `params.yaml`. Edit this file before running.

```yaml
filter_type: EnKF
Nens: 4
model_nprocs: 4
```

| Parameter | Description |
|-----------|-------------|
| `filter_type` | Ensemble filter variant to use |
| `Nens` | Number of ensemble members |
| `model_nprocs` | MPI processes per MATLAB/ISSM instance |

---

## Running the Model

### Serial Execution

For quick testing without MPI:

```bash
python run_da_issm.py --Nens=1 --model_nprocs=1 --verbose
```

### Parallel Execution (Local)

Distribute ensemble members across multiple processors using MPI:

```bash
mpirun -np <N> python run_da_issm.py --Nens=<N> --model_nprocs=<M> --verbose
```

| Flag | Description |
|------|-------------|
| `-np <N>` | Total MPI processes — set equal to `--Nens` |
| `--Nens=<N>` | Number of ensemble members |
| `--model_nprocs=<M>` | Processes per MATLAB instance |
| `--verbose` | Enable detailed debug output |

**Example** — 4 ensemble members, 4 processes each:

```bash
mpirun -np 4 python run_da_issm.py --Nens=4 --model_nprocs=4 --verbose
```

> **Note:** Setting `-np` equal to `--Nens` gives optimal resource utilization. Adjust `--model_nprocs` based on available cores.

---

### HPC Deployment with Apptainer

For reproducible, portable execution on HPC systems.

#### Container Definition

The Apptainer definition file is available at:

```
https://github.com/ICESEE-project/ICESEE-Containers/blob/main/spack-managed/combined-container/combined-env-inbuilt-matlab.def
```

#### Build

```bash
apptainer build combined-env-inbuilt-matlab.sif combined-env-inbuilt-matlab.def
```

#### Launch (SLURM + `srun`)

```bash
srun --mpi=pmix -n 4 apptainer exec \
  -B /path/to/ICESEE:/opt/ICESEE,examples:/opt/ISSM/examples,execution:/opt/ISSM/execution \
  combined-env-inbuilt-matlab.sif \
  with-issm python run_da_issm.py --Nens=4 --model_nprocs=1 --verbose
```

| Component | Description |
|-----------|-------------|
| `srun --mpi=pmix -n 4` | Launch 4 MPI processes via SLURM |
| `combined-env-inbuilt-matlab.sif` | Container image with ICESEE, ISSM, MATLAB runtime, MPI, and dependencies |
| `with-issm` | Activates the hybrid ICESEE + ISSM runtime environment inside the container |
| `--Nens=4` | 4 ensemble members |
| `--model_nprocs=1` | 1 process per ISSM model instance |
| `--verbose` | Enable launcher and server debug output |

#### Quick Environment Check

Verify the container environment before a full run:

```bash
apptainer exec combined-env-inbuilt-matlab.sif \
  with-issm bash -c '
    echo "PY=$(command -v python)"
    echo "MPI=$(command -v mpiexec)"
    echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
  '
```

---

## Bind Mounts Reference

The `-B` flag maps host directories into the container:

```
-B /host/path:/container/path
```

| Host Path | Container Path | Purpose |
|-----------|---------------|---------|
| `/path/to/ICESEE` | `/opt/ICESEE` | Local ICESEE source tree (live edits, no rebuild needed) |
| `examples/` | `/opt/ISSM/examples` | Writable ISSM example directory |
| `execution/` | `/opt/ISSM/execution` | ISSM runtime files and outputs (persisted outside container) |

> **Notes:**
> - Bind-mounting `/opt/ICESEE` is recommended during development — code changes are reflected immediately without rebuilding the `.sif`.
> - Ensure `examples/` and `execution/` exist and are writable before launching.
> - Use real directory paths rather than symlinks where possible.

---

## Outputs

Results are written to `_modelrun_datasets/` and `results/`, named by filter type:

```
results/<filter_type>-model.h5
```

To visualize results, open `read_results.ipynb` in Jupyter and follow the inline instructions to load and plot `.h5` output files. A MATLAB equivalent is available in `read_results.m`.

---

## Directory Structure

```
.
├── run_da_issm.py          # Main data assimilation script
├── params.yaml             # Model and filter configuration
├── requirements.txt        # Python dependencies
├── read_results.ipynb      # Jupyter notebook for result visualization
├── read_results.m          # MATLAB script for result visualization
├── results/                # Output .h5 files
└── _modelrun_datasets/     # Intermediate model run datasets
```

---

## Troubleshooting

**MPI errors**
Confirm the MPI implementation is correctly installed and that `mpirun`/`mpiexec` is on `$PATH`. Make sure `-np` matches `--Nens`.

**MATLAB server issues**
Verify MATLAB is installed and accessible:
```bash
which matlab
matlab -nodisplay -nosplash -r "disp('ok'); exit"
```
Inside a container, check the bind mount and that `$PATH` includes the MATLAB `bin` directory.

**Missing Python packages**
```bash
pip install -r requirements.txt
```

**General runtime issues**
Run with `--verbose` to get detailed launcher and MATLAB server output:
```bash
mpirun -np 4 python run_da_issm.py --Nens=4 --model_nprocs=1 --verbose
```

**Container environment issues**
Use the [quick environment check](#quick-environment-check) above to verify Python, MPI, and library paths inside the container before submitting a full job.