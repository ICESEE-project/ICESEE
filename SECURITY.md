# Security Policy

ICESEE is a research data-assimilation framework for ice-sheet and
geophysical models, typically run by a single user or research group on a
workstation or HPC cluster against local or collaborator-provided input
data. It does not run as a network service and has no built-in
authentication, network listener, or multi-tenant execution mode.

## Scope

Security-relevant reports are in scope for:

- code that parses external input (configuration YAML, HDF5/NetCDF/Zarr
  state and observation files, checkpoint/restart files) in a way that could
  execute arbitrary code, corrupt unrelated files, or crash in an unsafe way
  on untrusted input;
- dependency vulnerabilities pulled in via `pyproject.toml` /
  `requirements/` that are reachable from ICESEE's own code paths;
- credential or path handling in CI workflows (`.github/workflows/`) or
  packaging (`config/update_readme.py`, `pyproject.toml`) that could leak
  secrets or write outside the intended directory.

Out of scope: the scientific correctness of the data-assimilation
algorithms (report those as regular issues/PRs, not security reports), and
vulnerabilities in external coupled models (Icepack/Firedrake, ISSM/MATLAB)
themselves — report those upstream.

## Reporting a vulnerability

Please do not open a public GitHub issue for a suspected vulnerability.
Instead, use GitHub's private
["Report a vulnerability"](https://github.com/ICESEE-project/ICESEE/security/advisories/new)
flow on this repository, or contact the maintainer directly:

Brian Kyanjo — bkyanjo3@gatech.edu

Include the affected file(s)/function(s), a minimal reproduction if
possible, and the potential impact. We aim to acknowledge reports within a
reasonable time and will credit reporters in the fix's release notes unless
you ask not to be named.
