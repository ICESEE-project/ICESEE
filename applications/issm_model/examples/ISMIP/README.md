# ISSM Model with ICESEE Data Assimilation

This repository implements the Ice Sheet System Model (ISSM) integrated with ICESEE for data assimilation, demonstrated through the ISMIP example. The code supports both serial and parallel execution, utilizing ensemble Kalman filter (EnKF) variants for data assimilation experiments.

## Overview

The `run_da_issm.py` script combines the ISSM model with ICESEE to perform data assimilation using various EnKF methods. It supports parallel execution via MPI for high-performance computing and serial execution for testing different filter types.

## Prerequisites

- **Python**: Version 3.8 or higher
- **MPI**: An MPI implementation (e.g., OpenMPI or MPICH) for parallel runs
- **MATLAB**: Required for ISSM execution, with a configured MATLAB server
- **Dependencies**: Install Python packages listed in `requirements.txt`:
  ```bash
  pip install -r requirements.txt
  ```

## Running the Model

### Parallel Execution

Parallel runs leverage MPI to distribute tasks across multiple processors, with each MATLAB instance executing the ISSM model.

**Command**:
```bash
mpirun -np <N> python run_da_issm.py --Nens=<N> --model_nprocs=<M> --verbose
```

- `-np <N>`: Number of MPI processes (typically matches `--Nens`).
- `--Nens=<N>`: Number of ensemble members.
- `--model_nprocs=<M>`: Number of processes per MATLAB instance.
- `--verbose`: Optional flag for detailed debugging output.

**Example**:
To run with 4 ensemble members, each MATLAB instance using 4 processes:
```bash
mpirun -np 4 python run_da_issm.py --Nens=4 --model_nprocs=4 --verbose
```

**Notes**:
- Set `-np` equal to `--Nens` for optimal resource utilization.
- Adjust `--model_nprocs` based on available computational resources.
- Each MATLAB instance communicates with Python via the MATLAB server, running `<M>` processes concurrently.

### Serial Execution

Serial execution is ideal for testing EnKF variants, as not all support MPI parallelization.

**Command**:
```bash
python run_da_issm.py --Nens=30
```

**Configuration**:
- Set the `filter_type` in `params.yaml` to select the desired filter:
  - `EnKF`: Ensemble Kalman Filter
  - `DEnKF`: Deterministic Ensemble Kalman Filter
  - `EnTKF`: Ensemble Transform Kalman Filter
  - `EnRSKF`: Ensemble Square Root Kalman Filter

**Example**:
To run with the EnKF filter:
1. Update `params.yaml` to set `filter_type: EnKF`.
2. Run:
   ```bash
   python run_da_issm.py --parallel=serial
   ```

## Configuration

Model parameters are defined in `params.yaml` for easy modification. Key parameters include:
- `filter_type`: Specifies the ensemble filter.
- `Nens`: Number of ensemble members.
- `model_nprocs`: Processes per MATLAB instance (parallel runs only).
- Additional model-specific settings (see `params.yaml` for details).

**Example `params.yaml`**:
```yaml
filter_type: EnKF
Nens: 4
model_nprocs: 4
```

## Outputs

- **Results**: Saved as `.h5` files in the `_modelrun_datasets` and `results/` directories, named as:
  ```
  <filter_type>-model.h5
  ```
- **Visualization**: Use the `read_results.ipynb` Jupyter notebook to load and visualize results.

**Example**:
Open `read_results.ipynb` in Jupyter and follow the instructions to plot data from `.h5` files.

## Directory Structure

- `run_da_issm.py`: Main script for data assimilation.
- `params.yaml`: Configuration file for model parameters.
- `read_results.ipynb`: Jupyter notebook for visualizing results.
- `results/`: Directory for output `.h5` files.
- `requirements.txt`: List of Python dependencies.

## Notes

- Ensure the MATLAB server is running for parallel executions.
- Adjust `--model_nprocs` for large-scale runs to optimize performance.
- Use serial mode for testing and debugging filter variants.

## Troubleshooting

- **MPI Errors**: Confirm that the MPI implementation is correctly installed and configured.
- **MATLAB Server Issues**: Verify MATLAB is installed and accessible via the terminal.
- **Verbose Output**: Use the `--verbose` flag to diagnose runtime issues.


