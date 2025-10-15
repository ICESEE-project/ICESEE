# Source (src)

This directory contains the core implementation of the ICESEE data assimilation framework.

## Subdirectories

### EnKF/
Ensemble Kalman Filter implementations:
- `python_enkf/` - Pure Python EnKF implementation
- `cython_enkf/` - Optimized Cython EnKF for performance

### parallelization/
Parallelization infrastructure for distributed computing:
- MPI-based parallel I/O
- Parallel ensemble initialization
- Parallel forecast and analysis functions

### run_model_da/
Main entry points for running data assimilation workflows:
- Serial and parallel execution modes
- Full and partial parallelization strategies
- Localization and error generation functions

### tests/
Test suites and example implementations:
- Model-specific test cases
- Parallel MPI tests
- Zarr storage setup tests

### utils/
Common utility functions and tools used across the framework.

## Architecture

The ICESEE framework follows a modular design:
1. **EnKF Core**: Implements the mathematical operations for ensemble data assimilation
2. **Parallelization Layer**: Handles distributed computing with MPI
3. **Model Interface**: Connects external models to the DA framework
4. **I/O Layer**: Manages efficient data storage and retrieval

For detailed usage, refer to the [ICESEE Wiki](https://github.com/ICESEE-project/ICESEE/wiki).
