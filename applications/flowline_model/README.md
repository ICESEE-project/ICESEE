# Flowline Model

A simple 1D ice flow simulation model for testing and benchmarking ICESEE data assimilation workflows.

## Overview

The flowline model provides a lightweight ice sheet representation that is computationally efficient for rapid prototyping and testing of data assimilation algorithms. It simulates ice flow along a single flowline using simplified physics.

## Files

- `flowline_enkf.py` - Main EnKF implementation for the flowline model
- `config_loader.py` - Configuration file loader
- `params.yaml` - Model parameters and DA configuration
- `EnKF.ipynb` - Jupyter notebook demonstrating EnKF usage
- `run_flowline_enkf.ipynb` - Complete workflow notebook

## Usage

### Running in Jupyter

Open `run_flowline_enkf.ipynb` to see a complete example of:
1. Model initialization
2. Ensemble generation
3. Data assimilation cycles
4. Results visualization

### Command Line

```bash
python flowline_enkf.py --config params.yaml
```

## Status

Integration with ICESEE framework is currently underway. The model is functional for standalone testing and demonstration purposes.

## Configuration

Model parameters are specified in `params.yaml`, including:
- Ensemble size
- Observation frequency
- Model physics parameters
- Data assimilation settings

Refer to the configuration file for detailed parameter descriptions.
