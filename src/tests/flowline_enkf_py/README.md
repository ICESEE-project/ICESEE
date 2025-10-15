# Flowline EnKF Python Tests

Test suite for the Python flowline model with EnKF data assimilation.

## Purpose

This directory contains tests specifically for validating the flowline model integration with the Python EnKF implementation.

## Test Coverage

Tests include:
- Flowline model initialization
- Ensemble generation
- Forward model execution
- EnKF analysis step
- Complete DA cycle
- Configuration handling

## Running Tests

### All Tests
```bash
cd src/tests/flowline_enkf_py
python -m pytest
```

### Specific Test
```bash
python -m pytest test_flowline_enkf.py::test_ensemble_init
```

## Test Structure

Typical test structure:
1. Setup: Initialize model and configuration
2. Execute: Run DA cycle or component
3. Validate: Check results against expected values
4. Cleanup: Reset state for next test

## Expected Results

Tests validate:
- Ensemble dimensions are correct
- Analysis reduces ensemble spread
- RMSE improves with assimilation
- Configuration parameters are respected

## Dependencies

- pytest
- numpy
- ICESEE core modules
- Flowline model implementation

## Adding Tests

When adding new tests:
1. Follow existing naming conventions
2. Use descriptive test names
3. Include docstrings
4. Ensure tests are independent
5. Clean up resources after tests

## Debugging

For verbose output:
```bash
python -m pytest -v -s
```

For specific test debugging:
```bash
python -m pytest --pdb test_file.py::test_name
```
