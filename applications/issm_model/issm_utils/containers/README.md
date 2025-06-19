

# ISSM MATLAB Apptainer Image Usage Guide

This guide provides instructions for using the Apptainer (Singularity) image designed for running the **Ice Sheet System Model (ISSM)** with Python and MATLAB integration. The image is optimized for high-performance computing (HPC) clusters and includes MATLAB R2024b, ISSM, Python libraries, and MPI support.

## Prerequisites

- **HPC Cluster**: Access to a cluster with Apptainer/Singularity installed.
- **Apptainer Definition File**: `issm_matlab_container.def` for building the image.
- **Input Data**: Model input files for the `examples` directory.
- **Slurm Scheduler**: Cluster must support Slurm for `srun` commands.
- **MATLAB License**: Valid license accessible via `1711@matlablic.edu`.
- **Build Permissions**: `apptainer build` may require `fakeroot` or root privileges.

## Setup

### 1. Build the Apptainer Image

Build the image using the definition file:

```bash
apptainer build issm_matlab.sif issm_matlab_container.def
```

### 2. Create Directories

Set up directories for input and output:

```bash
mkdir -p examples execution
```

### 3. Populate Directories

- Place ISSM example data and scripts (e.g., model input files) in `examples/`.
- Keep `execution/` empty initially; it will store model outputs.

## Directory Structure

The following host directories are bound to the container:

- `examples/` → `/opt/ISSM/execution`: Contains example data and scripts.
- `execution/` → `/opt/execution`: Stores model outputs.

## Running the ISSM Script

Execute the `run_da_issm.py` script with:

```bash
srun -n 4 apptainer exec \
  -B examples:/opt/ISSM/execution,execution:/opt/execution \
  issm_matlab.sif python run_da_issm.py --Nens=2 --model_nprocs=2
```

### Command Breakdown

- `srun -n 4`: Launches 4 MPI tasks via Slurm.
- `apptainer exec`: Runs a command inside the container.
- `-B examples:/opt/ISSM/execution,execution:/opt/execution`: Binds host directories to container paths.
- `issm_matlab.sif`: The Apptainer image.
- `python run_da_issm.py`: The Python script to run.
- `--Nens=2`: Sets 2 ensemble members for data assimilation.
- `--model_nprocs=2`: Allocates 2 processors for the model simulation.

## Expected Behavior

- The `run_da_issm.py` script runs inside the container.
- Input files are read from `/opt/ISSM/execution` (host: `examples/`).
- Outputs are written to `/opt/execution` (host: `execution/`).
- MATLAB is invoked from Python using the integrated environment.

