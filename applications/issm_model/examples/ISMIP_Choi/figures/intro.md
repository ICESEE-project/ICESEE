# ICESEE on GHUB

ICESEE is an open-source framework for **ice-sheet state estimation** and **parameter inference** using ensemble-based data assimilation. Its core goal is to make transient Ensemble Kalman filtering workflows easier to apply across a range of ice-sheet and idealized models by connecting observations (glaciological and geological) with physics-based simulations.

This GHUB tool packages ICESEE into an interactive **Jupyter Book** format, so users can explore the framework, run lightweight examples, and understand the workflow without needing to set up a full HPC environment.

---

## Start here

::::{grid} 3
:gutter: 2

:::{grid-item}
```{button-link} user_manual
:color: info
:expand:
GHUB User Manual
```
:::

:::{grid-item}
```{button-link} icesee_jupyter_notebooks/run_lorenz96_da
:color: success
:expand:
Run Lorenz-96 Demo
```
:::

:::{grid-item}
```{button-link} https://github.com/ICESEE-project/ICESEE/wiki
:color: primary
:expand:
ICESEE Wiki
```
:::

::::

---

## What ICESEE is designed to do

ICESEE supports ensemble data assimilation workflows (e.g., EnKF-style methods) with an emphasis on:

- **Modular structure** so you can reuse the same DA logic across different models
- **Model coupling** (including external codes/workflows) while keeping the assimilation engine consistent
- **Scalability** toward HPC and cloud-style execution, including GHUB-style environments

For implementation details and broader documentation, see the ICESEE Wiki:
- https://github.com/ICESEE-project/ICESEE/wiki

---

## What you can do in this GHUB book

This book is organized to help you:

1. **Understand the ICESEE workflow** (ensembles, observations, forecasts, and analysis steps)
2. **Run a minimal, self-contained DA example** (Lorenz-96) to validate the environment and illustrate the full cycle
3. **Learn how ICESEE connects to larger ice-sheet models** (e.g., ISSM, Icepack) at a high level, including typical setup requirements

---

## How ICESEE is integrated here

This GHUB-facing repository (**ICESEE-GHUB**) is a *tool wrapper* around the main ICESEE codebase:

- The ICESEE source is pinned from the **ICESEE main branch** for reproducibility while tracking upstream progress
- The Jupyter Book provides documentation and runnable notebooks that rely on that pinned version

---

## Notes on dependencies

Some ICESEE components run with standard scientific Python packages, while others require external modeling stacks (e.g., full ice-sheet model toolchains). This GHUB book prioritizes **lightweight and portable** examples, and clearly separates:

- ✅ **Runs directly in GHUB** (Lorenz-96 DA demo)
- ⚠️ **Requires external installations** (ISSM / Icepack workflows)

---

## Where to go next

```{button-link} icesee_jupyter_notebooks/run_lorenz96_da
:color: success
:expand:
Go to Notebooks → run_lorenz96_da
```

```{button-link} index
:color: secondary
:expand:
Back to Home
```
