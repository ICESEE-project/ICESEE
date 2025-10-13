# ICESEE

**ICESEE** (ICE ShEet state and parameter Estimator) is a data assimilation software framework designed for coupling with ice sheet models such as **ISSM**, **Icepack**, and idealized models like **Lorenz-96**. It provides a modular, extensible platform for applying ensemble-based data assimilation techniques in glaciological modeling and beyond.

---

##  What is ICESEE?

ICESEE simplifies the implementation of advanced data assimilation workflows—such as the Ensemble Kalman Filter (EnKF) and its variants—across a range of geophysical models. It is designed with:

- A modular Python interface  
- Seamless integration with external model codes (MATLAB, Firedrake, ISSM, etc.)  
- Support for high-performance computing and containerized workflows  
- Scalability for future integration with cloud platforms like AWS and portals like GHUB  

---

##  Getting Started

To get started with ICESEE:

- [Installation Guide](https://github.com/ICESEE-project/ICESEE/wiki/1.-Installation)  
- [Using ICESEE](https://github.com/ICESEE-project/ICESEE/wiki/2.-Usage)  
- [Guide to Integrating Models](https://github.com/ICESEE-project/ICESEE/wiki/3.--Guide-to-Integrating-Models-into-the-ICESEE-Framework)  
- [Build ICESEE as a Package](https://github.com/ICESEE-project/ICESEE/wiki/4.-Build-ICESEE-as-a-package)  
- [Development Notes](https://github.com/ICESEE-project/ICESEE/wiki/5.-Development-Notes)

---

## Supported Models

- `icepack`: PDE-based modeling with Firedrake  
- `issm`: Finite-element ice sheet modeling (via MATLAB interface)  
- `lorenz96`: Idealized nonlinear DA benchmarking  
- `flowline_model`: Simple ice flow simulation  

---

## Documentation

Explore the [Wiki](https://github.com/ICESEE-project/ICESEE/wiki) to find:

- **[Installation Guide](https://github.com/ICESEE-project/ICESEE/wiki/1.-Installation)**: Setup instructions for ICESEE and dependencies  
- **[Usage Guide](https://github.com/ICESEE-project/ICESEE/wiki/2.-Usage)**: Running applications and using EnKF variants  
- **[Model Integration Guide](https://github.com/ICESEE-project/ICESEE/wiki/3.--Guide-to-Integrating-Models-into-the-ICESEE-Framework)**: How to implement new models  
- **[Development Notes](https://github.com/ICESEE-project/ICESEE/wiki/5.-Development-Notes)**: Project structure and best practices  
- **[Common Issues](https://github.com/ICESEE-project/ICESEE/wiki/6.-Common-Issues-and-solutions)**: Debugging common issues  
- **[Configuration Flags](config/README.md)**: Complete reference for all command-line and YAML parameters  

---

## Future Plans

- Integration with **AWS** for scalable cloud computing.
- Incorporation into the **GHUB online ice sheet platform** with enhanced features.

For questions or contributions, please open an issue or pull request on the [GitHub repository](https://github.com/ICESEE-project/ICESEE) or contact me at bkyanjo3@gatech.edu






