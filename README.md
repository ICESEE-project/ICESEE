# ICESEE

**ICESEE** (ICE ShEet state and parameter Estimator) is a data assimilation software framework designed for coupling with ice sheet models such as **ISSM**, **Icepack**, and idealized models like **Lorenz-96**. It provides a modular, extensible platform for applying ensemble-based data assimilation techniques in glaciological modeling and beyond.

---

## What is ICESEE?

ICESEE (ICE ShEet state and parameter Estimator) is a scalable ensemble-based data assimilation framework designed for coupling with ice-sheet, climate, and geophysical models such as ISSM, Icepack, and idealized systems like Lorenz-96.

ICESEE provides a modular and extensible scientific computing infrastructure for ensemble forecasting, uncertainty quantification, parameter estimation, and hybrid physics–AI workflows in large-scale environmental modeling.

The framework supports:
- Ensemble-based data assimilation methods
- Distributed-memory HPC execution
- Cloud-enabled scientific workflows
- Cross-language model coupling
- AI-enhanced scientific computing and machine learning integration

---

##  Getting Started

To get started with ICESEE:

- [Installation Guide](https://github.com/ICESEE-project/ICESEE/wiki/1.-Installation)  
- [Using ICESEE](https:https://github.com/ICESEE-project/ICESEE/wiki/2.-Usage)  
- [Build ICESEE as a package](https://github.com/ICESEE-project/ICESEE/wiki/3.-Build-ICESEE-as-a-package)  
- [Developmental notes](https://github.com/ICESEE-project/ICESEE/wiki/4.-Development-Notes)

>For Cluster installation and runs, see [ICESEE-Spack](https://github.com/ICESEE-project/ICESEE-Spack) or [ICESEE-Containers](https://github.com/ICESEE-project/ICESEE-Containers)
 and for cloud runs, see [ICESEE-GHUB](https://github.com/ICESEE-project/ICESEE-GHUB)

---

## AI/ML Integration

ICESEE is being extended with scientific machine learning and AI-enhanced workflows to support next-generation hybrid physics–AI modeling and scalable intelligent simulation systems.

Current and planned AI/ML capabilities include:

- Machine learning–based parameter estimation
- Neural-network observation operators
- AI-enhanced ensemble initialization
- Adaptive covariance and inflation tuning
- Surrogate forward models for accelerated simulations
- Intelligent workflow automation and HPC optimization
- AI-assisted diagnostics and parameter tuning

The AI/ML framework is designed to integrate seamlessly with existing ensemble-based data assimilation workflows while maintaining compatibility with distributed HPC environments and scientific modeling systems.

---

## Supported Models

- `icepack`: PDE-based modeling with Firedrake  
- `issm`: Finite-element ice sheet modeling (via MATLAB interface)  
- `lorenz96`: Idealized nonlinear DA benchmarking  
- `flowline_model`: Simple ice flow simulation  

---

## Documentation

Explore the Wiki to find:

- Configuration and setup tips  
- How to implement new models  
- How to extend or modify filters  
- Debugging common issues  

---

## Key Features

- Modular ensemble-based data assimilation framework
- Cross-language scientific model coupling
- Scalable HPC and distributed-memory workflows
- HDF5/NetCDF-based scientific data pipelines
- Containerized and cloud-enabled execution
- AI/ML-ready scientific workflow infrastructure
- Extensible APIs for integrating external models
- Support for uncertainty quantification and ensemble forecasting

---

## Future Directions

- Integration of scientific machine learning workflows into data assimilation pipelines
- Development of neural surrogate models for accelerated ice-sheet simulations
- AI-enhanced parameter estimation and adaptive ensemble tuning
- Integration with AWS and cloud-native scientific computing environments
- Expansion of distributed HPC workflows for large-scale ensemble forecasting
- Incorporation into the GHUB online ice-sheet platform
- Intelligent workflow orchestration and automated HPC diagnostics

For questions or contributions, please open an issue or pull request on the [GitHub repository](https://github.com/ICESEE-project/ICESEE) or contact me at bkyanjo3@gatech.edu

ICESEE is distributed as free and open-source software under a BSD-style license (see LICENSE). All external dependencies, including ISSM, Icepack, and other coupled models, are governed by their own licenses, which are independent of and do not impose restrictions on the ICESEE license.




