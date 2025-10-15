# Data Management Scripts

This directory contains scripts for managing ICESEE data files, including retrieval, inspection, post-processing, and aggregation.

## Scripts

### get_data.py
Utilities for retrieving and loading ICESEE data from various storage locations.

**Usage**:
```python
from scripts.data_management.get_data import ...
```

### inspect_files.py
Tools for inspecting ICESEE data files, including metadata examination and data structure validation.

**Features**:
- File format verification
- Metadata display
- Data structure inspection

### post_processing.py
Post-processing workflows for ICESEE output data.

**Capabilities**:
- Data cleaning and validation
- Statistical computations
- Output formatting
- Result extraction

### stack_icesee_data.py
Tools for stacking and aggregating ICESEE data across multiple runs or ensemble members.

**Use Cases**:
- Combining ensemble output
- Multi-run aggregation
- Time series concatenation
- Spatial data stacking

## File Formats

ICESEE primarily uses:
- HDF5 (.h5) for efficient array storage
- Zarr for cloud-optimized storage
- NetCDF for compatibility with geoscience tools

## Usage Examples

```python
# Load data
from scripts.data_management import get_data
data = get_data.load_ensemble_data('path/to/file.h5')

# Inspect file
from scripts.data_management import inspect_files
inspect_files.show_structure('path/to/file.h5')

# Stack data
from scripts.data_management import stack_icesee_data
stacked = stack_icesee_data.stack_runs(['run1.h5', 'run2.h5'])
```

## Dependencies

- h5py (HDF5 support)
- zarr (Zarr support)
- numpy
- ICESEE core utilities
