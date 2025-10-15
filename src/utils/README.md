# Utilities

Common utility functions and tools used throughout the ICESEE framework.

## Modules

### tools.py
General-purpose tools and helper functions.

**Categories**:
- File I/O utilities
- Data conversion functions
- Configuration helpers
- Path management
- Error handling utilities

**Common Functions**:
- File existence checks
- Directory creation
- Data type conversions
- Configuration parsing
- Logging setup

### utils.py
Core utility functions for data assimilation workflows.

**Categories**:
- Ensemble manipulation
- Statistical computations
- Array operations
- Diagnostic functions
- Validation utilities

**Common Functions**:
- Ensemble mean and spread calculations
- RMSE and other metrics
- Data validation
- Array reshaping and indexing
- Time handling

## Usage

Import utilities as needed:

```python
from src.utils import tools, utils

# Use file utilities
tools.ensure_directory_exists(output_path)

# Use ensemble utilities
mean, spread = utils.ensemble_statistics(ensemble)
rmse = utils.compute_rmse(forecast, truth)
```

## Design Philosophy

These utilities follow these principles:
- **Reusable**: Functions are general-purpose
- **Well-tested**: Covered by unit tests
- **Documented**: Clear docstrings
- **Efficient**: Optimized for common operations

## Common Operations

### Ensemble Statistics
```python
ensemble_mean = utils.compute_ensemble_mean(ensemble)
ensemble_spread = utils.compute_ensemble_spread(ensemble)
```

### File Operations
```python
tools.safe_create_directory(path)
tools.archive_old_results(directory)
```

### Data Validation
```python
utils.validate_ensemble_shape(ensemble, expected_shape)
utils.check_configuration(config_dict)
```

## Testing

Unit tests for these utilities are located in `src/tests/`.

## Adding New Utilities

When adding new utility functions:
1. Choose the appropriate module (tools.py or utils.py)
2. Add comprehensive docstrings
3. Include type hints
4. Add unit tests
5. Update this README

## Dependencies

These modules have minimal dependencies:
- numpy (array operations)
- Python standard library

This keeps them lightweight and universally usable across the framework.

## Best Practices

- Check for existing utilities before implementing new ones
- Keep functions focused and single-purpose
- Use descriptive function names
- Handle edge cases gracefully
- Return meaningful error messages
