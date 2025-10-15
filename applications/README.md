# Applications

This directory contains model-specific implementations and examples for the ICESEE data assimilation framework.

## Supported Models

### Flowline Model
A simple 1D ice flow simulation model used for testing and benchmarking data assimilation workflows.
- **Location**: `flowline_model/`
- **Status**: Integration underway

### Icepack Model
PDE-based ice sheet modeling using the Firedrake finite element framework.
- **Location**: `icepack_model/`
- **Status**: Fully supported
- **Dependencies**: Firedrake

### ISSM (Ice Sheet System Model)
Finite-element ice sheet modeling using MATLAB interface.
- **Location**: `issm_model/`
- **Status**: Development underway
- **Dependencies**: ISSM, MATLAB

### Lorenz96 Model
Idealized nonlinear dynamical system for data assimilation benchmarking and testing.
- **Location**: `lorenz_model/`
- **Status**: Fully supported

## Model Registration

Models are registered in `supported_models.py`, which provides a centralized interface for model discovery and initialization within the ICESEE framework.

## Structure

Each model directory typically contains:
- Example configurations and run scripts
- Model-specific utilities and helper functions
- Integration code for connecting to the ICESEE EnKF framework
