# Contributing to ICESEE

ICESEE (Ice Sheet statE and parameter Estimator) is a scientific-computing
framework. Contributions should preserve numerical correctness and existing
scientific behavior as strictly as they preserve code quality.

## Before you start

- Check open issues and the [Wiki](https://github.com/ICESEE-project/ICESEE/wiki)
  for existing design discussion — in particular the
  [model-integration guide](https://github.com/ICESEE-project/ICESEE/wiki/3.--Guide-to-Integrating-Models-into-the-ICESEE-Framework)
  if you're adding or changing a model.
- For architectural changes (new shared abstractions, communicator topology,
  state-management strategy), read `docs/architecture/` first and add an ADR
  there before writing code — see `docs/architecture/README.md`.
- Open an issue describing the change before large refactors or new
  dependencies, so the design can be discussed before code is written.

## Development setup

```bash
python -m pip install -U pip
python -m pip install -e ".[dev,mpi,viz]"
```

`requirements/ci.in` / `requirements/ci.txt` pin the exact environment CI
uses; regenerate with `pip-compile pyproject.toml --extra dev --extra mpi
--extra viz -o requirements/ci.txt` after changing `pyproject.toml`
dependencies (see the comment at the bottom of `pyproject.toml`).

`ICESEE` is imported as a namespace package rooted one directory above the
repository (`import ICESEE.src...`). `scripts/ci/ci_setup.sh` shows the
`PYTHONPATH` layout CI relies on; `Makefile`'s `setup` target does the
equivalent for a local shell.

## The `icesee_kwargs` convention

ICESEE passes exactly one flat runtime dictionary, `icesee_kwargs`, through
configuration loading, model initialization, forecast, analysis, and I/O.
Do not introduce a second, parallel context object (`params`, `model_kwargs`,
plain `kwargs` meant to hold config) or nest one context inside another —
this is the single most important convention in the codebase and is covered
in detail in Wiki page 3.

## Adding or changing a model integration

Follow the existing hook contract (`initialize_ensemble`,
`forecast_step_single`, `generate_true_state`, `generate_nurged_state`,
`inverse_step_single`, `post_analysis_update`, `Obs_fun`/`JObs_fun`/`Cov_Obs_fun`)
documented in Wiki page 3. Shared ICESEE infrastructure (`src/parallelization`,
`src/EnKF`, `src/utils`) must stay model-agnostic — do not add Firedrake,
ISSM/MATLAB, or any other model-specific import there. Model-specific logic
belongs in `applications/<model>/`. See `docs/model-adapter.md` for the
target conceptual interface this hook contract is evolving toward.

## MPI / HPC development guidance

- Prefer buffer-based MPI (`Bcast`/`Gather`/`Allreduce` on numpy arrays with
  an explicit MPI datatype) over the lowercase Python-object API
  (`bcast`/`gather`) for anything sized by state dimension, ensemble size, or
  observation count. Reserve the lowercase API for small, fixed-size metadata.
- Every rank that can reach a collective call must reach it in the same order
  every time — a rank-gated call (`if sub_rank == 0: comm.<collective>(...)`)
  on a communicator sized larger than the gated ranks is a deadlock risk, not
  just a style issue.
- Prefer `Abort()` over a bare `raise` inside code any rank might call
  collectively, so one rank's exception doesn't hang the others at the next
  collective.
- New MPI-parallel code should be exercised at `P < Ne`, `P = Ne`, `P > Ne`,
  and `P` not evenly divisible by `Ne`, not only the convenient case.
- See `docs/architecture/` for the intended `COMM_WORLD` → ensemble
  communicator → model communicator hierarchy new work should target.

## Testing

- Run `pytest -q src/tests` for the serial/non-Firedrake suite before
  opening a PR.
- MPI-dependent tests live under `src/tests/parallel_mpi/`; run them with
  `mpirun -np <N> pytest ...` as appropriate for the test.
- Firedrake/Icepack- and ISSM-dependent tests require those environments and
  are not part of the lightweight CI job in `.github/workflows/ci.yml` —
  run them manually against a real installation before merging a change that
  touches `applications/icepack_model/` or `applications/issm_model/`.
- Any change to shared analysis/state-management code (`src/EnKF/`,
  `src/parallelization/`) that could affect numerical output must include or
  reference a regression comparison (e.g. against the mode-2 parity contract
  in `docs/execution-mode-2-development.md`), not just a passing test suite.
- Performance-motivated changes must not change scientific results within
  the tolerances already established by the relevant parity test; if a
  change's effect on results is uncertain, say so in the PR description
  rather than asserting equivalence.

## Pull requests

- Keep PRs to one logically isolated change (one bug fix, one refactor, one
  feature) — do not mix performance changes with unrelated cleanup.
- Describe what changed and why, which tests were run, and any numerical
  regression evidence, using `.github/pull_request_template.md`.
- Target `develop` unless a maintainer asks for `main` directly.

## Commit messages

Describe the technical change only (what changed and why). Do not add
AI-authorship or AI-assistance trailers to commits, PR descriptions, or
source files — authorship attribution is the repository owner's call, not a
contributor default.
