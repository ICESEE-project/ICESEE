# Lorenz Model

Implementation of the Lorenz96 model for data assimilation testing and benchmarking.

## Overview

The Lorenz96 model is a simplified dynamical system originally designed to study atmospheric predictability. It serves as an excellent testbed for data assimilation algorithms due to its chaotic behavior and computational efficiency.

## Structure

- `examples/` - Example configurations and test cases
- `lorenz_utils/` - Utilities specific to Lorenz96 model integration

## Model Description

The Lorenz96 model is a system of ordinary differential equations with the form:

```
dX_i/dt = (X_{i+1} - X_{i-2}) * X_{i-1} - X_i + F
```

where:
- `X_i` represents the state variable at location i
- `F` is a forcing parameter (typically 8)
- Cyclic boundary conditions are applied

## Use Cases

This model is ideal for:
- Testing new data assimilation algorithms
- Benchmarking EnKF performance
- Understanding nonlinear DA behavior
- Educational demonstrations
- Algorithm parameter tuning

## Status

Fully supported and extensively tested.

## Examples

The `examples/` directory contains various scenarios including:
- Different system sizes (number of variables)
- Various observation configurations
- Perfect and imperfect model experiments
- Parameter estimation tests

## Getting Started

See the examples in `examples/lorenz96/` for complete workflows demonstrating:
- Model initialization
- Ensemble generation
- Observation simulation
- EnKF application
- Results analysis

## Resources

- Original Lorenz96 paper: Lorenz, E. N. (1996). "Predictability: A problem partly solved"
- [ICESEE Wiki](https://github.com/ICESEE-project/ICESEE/wiki)
