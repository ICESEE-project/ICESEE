# Ensemble Kalman Filter (EnKF)

This directory contains the core Ensemble Kalman Filter implementations for ICESEE.

## Implementations

### python_enkf/
Pure Python implementation of the EnKF algorithm.

**Files**:
- `EnKF.py` - Main EnKF implementation with various algorithms
- `enkf_class_python.py` - Object-oriented EnKF class interface

**Features**:
- Standard EnKF
- Ensemble Transform Kalman Filter (ETKF)
- Ensemble Square Root Filter (EnSRF)
- Localization support
- Inflation techniques

**Advantages**:
- Easy to understand and modify
- Pure Python (no compilation needed)
- Good for prototyping

### cython_enkf/
Optimized Cython implementation for performance-critical applications.

**Files**:
- `enkf.pyx` - Cython-optimized EnKF implementation
- `setup.py` - Build configuration

**Features**:
- Highly optimized matrix operations
- C-level performance
- Reduced memory footprint
- Same algorithms as Python version

**Advantages**:
- Faster execution
- Better scaling for large ensembles
- Efficient for production runs

## Usage

### Python EnKF
```python
from src.EnKF.python_enkf import EnKF

enkf = EnKF(ensemble_size=50)
analysis = enkf.update(forecast, observations, obs_operator)
```

### Cython EnKF
First compile:
```bash
cd src/EnKF/cython_enkf
python setup.py build_ext --inplace
```

Then use similarly to Python version.

## Algorithm Variants

Both implementations support:
- **Deterministic EnKF**: Uses ensemble mean for analysis
- **Stochastic EnKF**: Perturbed observations
- **ETKF**: Transform-based, deterministic
- **EnSRF**: Square root filter, deterministic

## Localization

Covariance localization is supported through:
- Gaspari-Cohn correlation function
- Distance-based cutoff
- Custom localization matrices

## Performance

For small to medium problems (ensemble size < 100), Python implementation is sufficient. For larger problems or production runs, use the Cython version.

## References

- Evensen, G. (2003). "The Ensemble Kalman Filter: theoretical formulation and practical implementation"
- Hunt et al. (2007). "Efficient data assimilation for spatiotemporal chaos: A local ensemble transform Kalman filter"
