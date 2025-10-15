# ISSM Utilities

This directory contains utilities for ISSM model integration with ICESEE.

## Purpose

These utilities facilitate:
- MATLAB-Python communication
- ISSM data conversion to ICESEE formats
- Process management for MATLAB instances
- Model state and parameter handling
- Results extraction and post-processing

## Key Components

### matlab2python/
Tools for converting MATLAB data structures to Python formats and vice versa.

### containers/
Container configurations for running ISSM with ICESEE in containerized environments.

## Usage

These utilities bridge the gap between ISSM's MATLAB environment and ICESEE's Python framework, enabling seamless data exchange during data assimilation cycles.

## MATLAB Interface

Special consideration is given to:
- Efficient data transfer between MATLAB and Python
- Proper MATLAB process lifecycle management
- Error handling across the MATLAB-Python boundary

## Dependencies

- MATLAB
- ISSM (properly configured in MATLAB)
- Python MATLAB Engine API
- ICESEE core framework

Refer to individual utility modules for specific usage patterns and requirements.
