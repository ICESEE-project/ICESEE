# Python EnKF

Pure Python implementation of the Ensemble Kalman Filter for ICESEE.

## Files

### EnKF.py
Main module containing EnKF algorithms and supporting functions.

**Key Functions**:
- Ensemble analysis updates
- Covariance calculations
- Localization functions
- Inflation algorithms
- Observation operator handling

### enkf_class_python.py
Object-oriented interface for the EnKF.

**Features**:
- Clean class-based API
- State management
- Configuration handling
- History tracking

## Usage

### Functional Interface
```python
from src.EnKF.python_enkf.EnKF import enkf_update

X_analysis = enkf_update(
    X_forecast,
    y_obs,
    H,
    R,
    localization=None
)
```

### Class Interface
```python
from src.EnKF.python_enkf.enkf_class_python import EnKFPython

enkf = EnKFPython(config)
X_analysis = enkf.analysis_step(X_forecast, y_obs)
```

## Algorithm Details

### Standard EnKF
The ensemble mean and covariance are used to compute the Kalman gain and update each ensemble member.

### Localization
Covariance localization reduces spurious long-range correlations:
- Gaspari-Cohn function
- Schur product with localization matrix
- Distance-based tapering

### Inflation
Ensemble inflation counters variance underestimation:
- Multiplicative inflation
- Additive inflation
- Adaptive inflation

## Advantages

- **Readable**: Easy to understand algorithm implementation
- **Flexible**: Simple to modify and extend
- **No Compilation**: Works immediately without build step
- **Debugging**: Easier to debug than compiled code

## Performance Considerations

- Suitable for small to medium ensemble sizes (< 100 members)
- For larger problems, consider the Cython implementation
- Numpy vectorization provides reasonable performance

## Testing

Unit tests for this module are located in `src/tests/`.

## References

See the main EnKF README for algorithm references.
