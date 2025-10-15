# MATLAB Scripts

This directory contains scripts for managing MATLAB processes and interfacing with MATLAB-based models like ISSM.

## Scripts

### kill_matlab_processes.py
Utility for safely terminating MATLAB processes that may be left running after data assimilation runs.

**Purpose**:
- Clean up orphaned MATLAB processes
- Prevent resource leaks
- Manage MATLAB instance lifecycle

**Usage**:
```bash
python kill_matlab_processes.py
```

**Features**:
- Identifies MATLAB processes
- Safe termination procedures
- Process filtering options

## Background

When running MATLAB-based models (like ISSM) with ICESEE, MATLAB processes may occasionally persist after the main workflow completes. These scripts ensure proper cleanup and resource management.

## Safety

These utilities are designed to safely terminate MATLAB processes without corrupting data or affecting other applications. Always ensure important work is saved before running cleanup scripts.

## Related Components

For MATLAB-Python integration, see:
- `applications/issm_model/issm_utils/` - MATLAB interface utilities
- MATLAB Engine API documentation

## Dependencies

- psutil (process management)
- Python 3.x
