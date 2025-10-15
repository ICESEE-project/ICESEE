# ICESEE Copilot Instructions

## Project Overview

ICESEE (ICE ShEet state and parameter Estimator) is a data assimilation software framework designed for coupling with ice sheet models such as ISSM, Icepack, and idealized models like Lorenz-96. It provides a modular, extensible platform for applying ensemble-based data assimilation techniques in glaciological modeling.

## Technology Stack

- **Primary Language**: Python 3.6+
- **Key Dependencies**: NumPy, SciPy, MPI4py, Firedrake, H5py, Zarr
- **Parallelization**: MPI, Dask, Ray, Multiprocessing
- **Testing**: Python unittest/pytest
- **Build System**: setuptools

## Project Structure

```
ICESEE/
├── src/                          # Core source code
│   ├── EnKF/                     # Ensemble Kalman Filter implementations
│   │   ├── python_enkf/          # Python-based EnKF
│   │   └── cython_enkf/          # Cython-optimized EnKF
│   ├── parallelization/          # Parallel computing implementations
│   ├── run_model_da/             # Data assimilation run scripts
│   └── utils/                    # Utility functions
├── applications/                 # Model-specific applications
│   ├── icepack_model/            # Icepack (Firedrake) model integration
│   ├── issm_model/               # ISSM model integration
│   └── lorenz_model/             # Lorenz-96 benchmark model
├── config/                       # Configuration files and flags
└── scripts/                      # Utility scripts
```

## Coding Standards

### Python Style

- Follow standard Python naming conventions:
  - `snake_case` for functions and variables
  - `PascalCase` for class names
  - `UPPER_CASE` for constants
- Include module-level docstrings with author, date, and description
- Use descriptive variable names, especially for scientific variables
- Prefer explicit imports over wildcard imports (avoid `from module import *`)

### Documentation

- Include file headers with:
  ```python
  # ==============================================================================
  # @author: [Author Name]
  # @date: YYYY-MM-DD
  # @description: [Description of module/file purpose]
  # ==============================================================================
  ```
- Document class methods and complex functions with docstrings
- Explain scientific/mathematical operations with comments where appropriate

### Code Organization

- Keep model-specific code in the `applications/` directory
- Core algorithms and utilities go in `src/`
- Use the namespace package structure: `ICESEE.src.`, `ICESEE.applications.`, etc.
- Import from ICESEE using full paths: `from ICESEE.src.EnKF.python_enkf.EnKF import EnsembleKalmanFilter`

## Key Concepts

### Data Assimilation

- **EnKF**: Ensemble Kalman Filter - the core data assimilation method
- **Forecast Step**: Propagate ensemble members forward in time
- **Analysis Step**: Update ensemble based on observations
- **Stochastic vs Deterministic**: Different EnKF implementations (EnKF, EnSRF, DEnKF, EnTKF)

### Ensemble Management

- Ensemble members represent different realizations of model states
- State variables typically include ice thickness (`h`), velocity components (`u`, `v`)
- Parameter estimation can be joint with state estimation

### Configuration

- YAML files are used for configuration parameters
- Command-line flags are defined in `config/README.md`
- Key parameters include:
  - `Nens`: Number of ensemble members
  - `freq_obs`: Observation frequency
  - `sig_obs`: Observation error standard deviation
  - `sig_model`: Model error standard deviation

## Common Tasks

### Adding a New Model

1. Create a new directory under `applications/[model_name]/`
2. Implement the forecast step function following the pattern:
   ```python
   def forecast_step_single(ensemble=None, **kwargs):
       """Advance a single ensemble member forward in time"""
       # Implementation
       return updated_ensemble
   ```
3. Create example configurations in `applications/[model_name]/examples/`

### Working with Parallelization

- Use the `ParallelManager` from `ICESEE.src.parallelization.parallel_mpi`
- Ensure functions to be parallelized are pickleable
- Be mindful of MPI context when using parallel operations
- Set `OMP_NUM_THREADS=1` to avoid thread conflicts with MPI

### Data Storage

- HDF5 (`.h5`) files are used for storing ensemble data
- Zarr is available for alternative storage formats
- Configure chunking and compression parameters in YAML configs
- Follow the pattern: `data_path` for output, `obs_data_path` for observations

## Testing

- Test files are located in `src/tests/`
- Run tests before committing changes
- When modifying core EnKF algorithms, verify with existing example applications
- Test with different parallelization modes if changing parallel code

## Environment Setup

- Use `setup_venv.sh` (Linux/Mac) or `setup_venv.bat` (Windows) to set up the environment
- Install with: `pip install -e .` for development
- Set PYTHONPATH if not using package installation
- Be aware of external dependencies like Firedrake, ISSM that may require separate installation

## Performance Considerations

- Large ensemble sizes require careful memory management
- Use chunked I/O for large datasets
- Consider using Cython implementations for performance-critical code
- Profile MPI communication patterns for scalability

## External Model Integration

- **Icepack**: Uses Firedrake for PDE-based modeling
- **ISSM**: Integrates via MATLAB interface
- **Lorenz-96**: Pure Python implementation for benchmarking
- Each model has specific initialization and simulation routines

## When Making Changes

- Maintain backward compatibility with existing example configurations
- Update documentation if adding new parameters or features
- Consider impact on parallel execution when modifying core algorithms
- Test with at least one example application before considering changes complete
- Be mindful of scientific correctness in data assimilation algorithms
