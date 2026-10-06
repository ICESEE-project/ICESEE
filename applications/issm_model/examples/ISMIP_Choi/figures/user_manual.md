# User Manual

This section provides practical guidance for using the **ICESEE GHUB Tool**.  
It focuses on what is available in this environment, how to run the included material, and where to look for more advanced workflows.

---

## Using this Jupyter Book

The navigation sidebar organizes the book into:

- **Getting Started** pages for context and orientation  
- A small set of **tutorial notebooks** that run directly in GHUB  
- Reference links to the upstream ICESEE documentation  

Pages and notebooks can be opened directly in the browser without additional setup.

---

## Running notebooks in GHUB

Tutorial notebooks in this tool are designed to be lightweight and runnable in an interactive session.

General recommendations:

- Run cells in order from top to bottom  
- Allow a few seconds for ensemble steps to complete  
- Re-run individual sections to explore parameter sensitivity  

```{admonition} Note
:class: note

Execution time in GHUB is limited compared to HPC environments.  
Examples here prioritize clarity and portability.
```

---

## File and repository layout

The GHUB-facing repository is organized as:

- `index.md` — landing page  
- `intro.md` — overview of ICESEE on GHUB  
- `user_manual.md` — practical usage notes (this page)  
- `icesee_jupyter_notebooks/` — runnable tutorial notebooks  
- `_config.yml`, `_toc.yml` — Jupyter Book configuration  

This structure keeps the book stable while tracking the upstream ICESEE framework.

---

## External workflows

Some ICESEE applications require full ice-sheet modeling stacks (e.g., ISSM or Icepack) and are not executed directly inside GHUB.

For those workflows, this book provides:

- High-level guidance and templates  
- Pointers to upstream documentation  
- A starting point for running on dedicated systems  

Upstream reference:

- https://github.com/ICESEE-project/ICESEE/wiki

---

## Support and feedback

If you encounter issues or have suggestions:

- Use the GitHub issue tracker  
- Consult the ICESEE Wiki for developer notes  

```{button-link} https://github.com/ICESEE-project/ICESEE/issues
:color: primary
:expand:
Report an Issue
```

---

```{button-link} index
:color: secondary
:expand:
Back to Home
```
