# Flowline EnKF Python-Julia Tests

Test suite for flowline model with Python-Julia interoperability.

## Purpose

This directory contains tests for the flowline model that utilize both Python and Julia implementations, testing the interface between the two languages.

## Overview

These tests validate:
- Python-Julia data exchange
- Julia model execution from Python
- Performance comparisons
- Numerical consistency

## Background

Julia can provide performance benefits for certain numerical operations. These tests ensure that:
- Julia implementations produce correct results
- Data transfers between Python and Julia work correctly
- Performance improvements are realized

## Requirements

- Python with PyJulia installed
- Julia runtime
- Required Julia packages
- ICESEE framework

## Running Tests

### Setup Julia Environment
```bash
# First time setup
julia setup_julia_env.jl
```

### Run Tests
```bash
cd src/tests/flowline_enkf_py_jl
python -m pytest
```

## Performance Testing

Compare Python vs Julia performance:
```bash
python benchmark_python_vs_julia.py
```

## Test Categories

1. **Correctness Tests**: Verify Julia implementation matches Python
2. **Performance Tests**: Measure speedup from Julia
3. **Interface Tests**: Validate data exchange
4. **Integration Tests**: Complete DA workflows

## Troubleshooting

### Julia Not Found
Ensure Julia is in PATH:
```bash
julia --version
```

### PyJulia Issues
Reinstall PyJulia:
```bash
pip install julia
python -c "import julia; julia.install()"
```

### Numerical Differences
Small floating-point differences between Python and Julia are expected. Tests use appropriate tolerances.

## Status

This is an experimental feature exploring Julia integration for performance-critical operations.
