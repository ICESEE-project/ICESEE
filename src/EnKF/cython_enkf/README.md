# Cython EnKF

High-performance Cython implementation of the Ensemble Kalman Filter for ICESEE.

## Files

### enkf.pyx
Cython source code with optimized EnKF implementations.

**Features**:
- C-optimized matrix operations
- Memory-efficient algorithms
- Parallel operation support
- Same functionality as Python version

### setup.py
Build configuration for compiling the Cython extension.

**Dependencies**:
- Cython
- NumPy
- C compiler (gcc, clang, MSVC)

## Building

### Compile the Extension
```bash
cd src/EnKF/cython_enkf
python setup.py build_ext --inplace
```

This creates a compiled module (.so on Linux/Mac, .pyd on Windows) that can be imported like a regular Python module.

### Development Build
For development with debugging symbols:
```bash
python setup.py build_ext --inplace --debug
```

## Usage

After compilation, use like the Python version:

```python
from src.EnKF.cython_enkf import enkf_update

X_analysis = enkf_update(
    X_forecast,
    y_obs,
    H,
    R,
    localization=None
)
```

## Performance

### Speed Improvements
Typical speedups compared to Python version:
- 5-10x for medium ensembles (50-100 members)
- 10-20x for large ensembles (>100 members)
- Better scaling with problem size

### Memory Efficiency
- Reduced memory allocations
- In-place operations where possible
- Efficient array handling

## When to Use

Use the Cython EnKF when:
- Running production data assimilation
- Working with large ensembles (>100 members)
- Performance is critical
- Running many DA cycles

Use the Python EnKF when:
- Prototyping new algorithms
- Debugging
- Working with small problems
- Ease of modification is important

## Compilation Requirements

- Python development headers
- Cython (pip install cython)
- NumPy development headers
- C compiler toolchain

## Troubleshooting

If compilation fails:
1. Ensure Cython is installed: `pip install cython`
2. Check C compiler availability: `gcc --version`
3. Verify NumPy installation: `python -c "import numpy; print(numpy.get_include())"`
4. Check the `out` file for compilation logs

## Testing

Performance benchmarks and correctness tests are in `src/tests/`.

## Maintenance

When updating the algorithm:
1. Modify enkf.pyx
2. Rebuild: `python setup.py build_ext --inplace`
3. Test: Run test suite
4. Update this README if interface changes
