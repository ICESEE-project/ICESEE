# ICESEE

**ICESEE** (Ice Sheet State and Parameter Estimator) is a data assimilation software framework designed for coupling with ice sheet models such as **ISSM**, **Icepack**, and idealized models like **Lorenz-96**. It provides a modular, extensible platform for applying ensemble-based data assimilation techniques in glaciological modeling and beyond.

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
- [Using ICESEE](https:https://github.com/ICESEE-project/ICESEE/wiki/2.-Usage)  
- [Build ICESEE as a package](https://github.com/ICESEE-project/ICESEE/wiki/3.-Build-ICESEE-as-a-package)  
- [Developmental notes](https://github.com/ICESEE-project/ICESEE/wiki/4.-Development-Notes)

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

## Future Plans

- Integration with **AWS** for scalable cloud computing.
- Incorporation into the **GHUB online ice sheet platform** with enhanced features.

For questions or contributions, please open an issue or pull request on the [GitHub repository](https://github.com/ICESEE-project/ICESEE) or contact me at bkyanjo3@gatech.edu






