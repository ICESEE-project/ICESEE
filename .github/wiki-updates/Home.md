# Welcome to the ICESEE Wiki

**Ice Sheet State and Parameter Estimator (ICESEE)** is a state-of-the-art data assimilation software package designed for ice sheet models. This advanced software facilitates the creation of an adaptive intelligent wrapper with robust protocols and APIs to seamlessly couple and integrate with various ice sheet models. The primary objective is to simplify the interaction between different models, enabling the adoption of complex data assimilation techniques across multiple frameworks.

This design will be extended to integrate with cloud computing services such as **AWS**, ensuring scalability and efficiency for larger simulations. Eventually, the software will be incorporated into the **GHUB online ice sheet platform**, significantly enhancing its capabilities by including the new features currently under development.

---

## 📚 Documentation Overview

This Wiki provides comprehensive documentation for installing, using, and extending ICESEE. Whether you're a new user or a developer, you'll find guides and references to help you work with ICESEE effectively.

---

## 🚀 Quick Start

New to ICESEE? Start here:

1. **[Installation Guide](1.-Installation)** - Set up ICESEE on your system
2. **[Usage Guide](2.-Usage)** - Learn how to run ICESEE applications and use different EnKF variants
3. **[Configuration Flags](https://github.com/ICESEE-project/ICESEE/blob/main/config/README.md)** - Reference for all command-line and YAML configuration options

---

## 📖 User Guides

### Getting Started
- **[1. Installation](1.-Installation)** - Complete installation instructions for ICESEE and its dependencies
- **[2. Usage](2.-Usage)** - Running ICESEE applications, supported models, and data assimilation methods

### Building and Packaging
- **[4. Build ICESEE as a Package](4.-Build-ICESEE-as-a-package)** - Instructions for building and distributing ICESEE to PyPI

---

## 🔧 Developer Guides

### Integration and Extension
- **[3. Guide to Integrating Models into ICESEE](3.--Guide-to-Integrating-Models-into-the-ICESEE-Framework)** - Comprehensive guide for integrating new models (MPI and non-MPI)
- **[5. Development Notes](5.-Development-Notes)** - Project structure, namespace packages, and development best practices

### Troubleshooting
- **[6. Common Issues and Solutions](6.-Common-Issues-and-solutions)** - Resolving common installation and runtime issues

---

## 🖥️ Platform-Specific Guides

### macOS Installation
- **[7. ISSM MATLAB Installation Guide for macOS](7.-ISSM-Matlab-Installation-Guide-for-macOS)** - Step-by-step guide for installing ISSM with MATLAB on macOS
- **[7.1 ISSM Dual Build Setup](7.1-ISSM-Dual-Build-Setup)** - Building ISSM with both MATLAB and Python interfaces using Makefile

---

## 🎯 Key Features

- **Modular Python Interface**: Easy integration with various ice sheet models
- **Multiple EnKF Variants**: EnKF, DEnKF, EnTKF, EnRSKF with inflation and localization support
- **MPI Parallelization**: Fully parallel and partially parallel modes for scalable computation
- **Supported Models**: Icepack, ISSM, Lorenz-96, and Flowline models
- **Container Support**: Docker and Apptainer recipes for portable deployments
- **HPC Ready**: Designed for high-performance computing environments

---

## 🔗 External Resources

- **[GitHub Repository](https://github.com/ICESEE-project/ICESEE)** - Source code and issue tracking
- **[Main README](https://github.com/ICESEE-project/ICESEE/blob/main/README.md)** - Project overview and quick links
- **[Configuration Documentation](https://github.com/ICESEE-project/ICESEE/blob/main/config/README.md)** - Complete flag reference

---

## 💬 Getting Help

- **Questions?** Open an issue on the [GitHub repository](https://github.com/ICESEE-project/ICESEE/issues)
- **Contributions?** Submit a pull request or contact the maintainer at bkyanjo3@gatech.edu
- **Bug Reports?** Use the GitHub issue tracker with detailed reproduction steps

---

## 📋 Wiki Contents

### Installation and Setup
1. [Installation](1.-Installation)
2. [Usage](2.-Usage)
4. [Build ICESEE as a Package](4.-Build-ICESEE-as-a-package)

### Development and Integration
3. [Guide to Integrating Models into ICESEE](3.--Guide-to-Integrating-Models-into-the-ICESEE-Framework)
5. [Development Notes](5.-Development-Notes)

### Troubleshooting and Platform-Specific
6. [Common Issues and Solutions](6.-Common-Issues-and-solutions)
7. [ISSM MATLAB Installation Guide for macOS](7.-ISSM-Matlab-Installation-Guide-for-macOS)
7.1. [ISSM Dual Build Setup](7.1-ISSM-Dual-Build-Setup)

---

**Last Updated**: January 2025
