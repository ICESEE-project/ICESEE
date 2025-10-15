# Icepack Model

Integration of the Icepack ice sheet model with the ICESEE data assimilation framework.

## Overview

Icepack is a Python library for modeling ice sheets using finite element methods via the Firedrake framework. This directory provides the interface between Icepack and ICESEE's Ensemble Kalman Filter implementation.

## Structure

- `examples/` - Example configurations and test cases for different ice geometries and scenarios
- `icepack_utils/` - Utilities specific to Icepack model integration

## Dependencies

- Firedrake finite element library
- Icepack ice sheet modeling library
- ICESEE core framework

## Examples

The `examples/` directory contains various test cases including:
- Idealized geometries
- Synthetic ice streams
- Realistic glacier scenarios

Each example typically includes:
- Configuration files
- Mesh generation scripts
- Initial condition setup
- Observation generation

## Status

Fully supported and actively used for ice sheet state and parameter estimation.

## Getting Started

Refer to the examples in the `examples/` directory for complete workflows. Each example contains its own README with specific instructions.

## Resources

- [Icepack Documentation](https://icepack.github.io/)
- [Firedrake Documentation](https://www.firedrakeproject.org/)
- [ICESEE Wiki](https://github.com/ICESEE-project/ICESEE/wiki)
