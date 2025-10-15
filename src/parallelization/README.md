# Parallelization

This directory contains the parallelization infrastructure for distributed data assimilation with ICESEE.

## Overview

ICESEE uses MPI (Message Passing Interface) for parallel execution on HPC systems. This directory provides utilities for parallel ensemble operations, I/O, and communication.

## Structure

- `parallel_mpi/` - MPI-specific implementations
- Core parallel I/O functions
- MPI communication wrappers
- Parallel analysis and forecast functions

## Key Modules

### EnKF_parallel_io.py
High-level parallel I/O operations for ensemble data.

**Features**:
- Parallel HDF5/Zarr reading and writing
- Distributed ensemble storage
- Checkpoint management
- Restart capabilities

### _parallel_i_o.py
Low-level parallel I/O primitives.

### MPI Analysis Functions
- `_mpi_analysis_functions.py` - Parallel EnKF analysis step
- Distributed covariance calculations
- Parallel localization

### MPI Forecast Functions
- `_mpi_forecast_functions.py` - Parallel ensemble forecasting
- Load balancing across processors
- Task distribution

### MPI Initialization
- `_mpi_ensemble_intialization.py` - Parallel ensemble generation
- Distributed initial conditions

### MPI Observations
- `_mpi_generate_synthetic_observations.py` - Parallel observation generation
- `_mpi_generate_true_wrong_state.py` - Truth and background states

## Parallelization Strategy

### Domain Decomposition
Ensemble members are distributed across MPI processes:
- Each process handles a subset of ensemble members
- Communication occurs during analysis step
- Efficient for large ensembles

### Data Parallelism
Model runs are embarrassingly parallel:
- Minimal communication during forecast
- Gather operations for analysis
- Scalable to many processors

## Usage

### MPI Execution
```bash
mpirun -np 48 python run_parallel_da.py --config params.yaml
```

### Configuration
```yaml
parallel_flag: true
n_modeltasks: 48
model_nprocs: 1  # processors per model instance
```

## Performance

### Scaling
- **Strong scaling**: Fixed problem size, increasing processors
- **Weak scaling**: Problem size scales with processors
- Near-ideal scaling for forecast step
- Communication overhead in analysis step

### Optimization Tips
- Use parallel I/O for large datasets
- Balance ensemble size with number of processors
- Consider communication costs
- Use localization to reduce communication

## Requirements

- MPI implementation (OpenMPI, MPICH, Intel MPI)
- mpi4py Python package
- Parallel HDF5 (optional but recommended)
- Sufficient network bandwidth for large ensembles

## Debugging

### Common Issues
- MPI initialization failures
- Deadlocks in communication
- Load imbalance
- I/O bottlenecks

### Tools
```bash
# Check MPI setup
mpirun -np 2 python -c "from mpi4py import MPI; print(MPI.COMM_WORLD.rank)"

# Profile execution
mpirun -np 48 python -m cProfile -o profile.out run_parallel_da.py
```

## References

- [mpi4py Documentation](https://mpi4py.readthedocs.io/)
- [Parallel HDF5](https://docs.h5py.org/en/stable/mpi.html)
