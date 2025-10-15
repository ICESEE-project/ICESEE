# Lorenz Model Examples

This directory contains example configurations and test cases for the Lorenz96 model with ICESEE.

## Overview

These examples demonstrate various data assimilation scenarios using the Lorenz96 model, from basic EnKF applications to advanced parameter estimation.

## Example Structure

Each example typically includes:
- Configuration files (YAML)
- Python scripts or Jupyter notebooks
- Expected results and validation data
- Documentation

## Types of Examples

Examples cover:
- **Basic State Estimation**: Standard EnKF with perfect model
- **Parameter Estimation**: Joint state-parameter estimation
- **Imperfect Model**: Testing robustness to model error
- **Observation Strategies**: Different observation networks
- **Ensemble Sizes**: Performance with varying ensemble sizes
- **Localization**: Covariance localization techniques

## Running Examples

Lorenz96 examples are computationally lightweight and can typically run on a laptop:

```bash
cd lorenz96/
python run_example.py --config params.yaml
```

Or use the provided Jupyter notebooks for interactive exploration.

## Educational Use

These examples are excellent for:
- Learning data assimilation concepts
- Teaching EnKF principles
- Algorithm development and testing
- Benchmarking new methods

## Quick Start

Start with the basic example in `lorenz96/` to understand the fundamental workflow, then explore more advanced scenarios.
