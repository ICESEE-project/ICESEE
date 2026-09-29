# 0001 — Memory-bounded ensemble analysis (mode 2)

## Status

Accepted. Implemented in `src/parallelization/EnKF_parallel_io.py` and
described in detail in `docs/execution-mode-2-development.md` and
`docs/large-scale-execution.md`.

## Context

Ensemble Kalman Filter analysis naively requires an `Nx x Ne` ensemble
matrix and an `Ny x Ny` (or `Nx x Nx`) covariance/gain structure. Neither
scales to large multidimensional model states (`Nx` in the millions) or
large ensembles. ICESEE's legacy serial analysis path
(`src/EnKF/python_enkf/EnKF.py`) does build these dense structures and does
not scale past small states.

## Decision

Execution mode 2's analysis (`compute_analysis_update`,
`compute_X5_utils_` in `EnKF_parallel_io.py`) partitions state rows across
MPI ranks, streams state and observation rows from HDF5 in
memory-budgeted chunks (`analysis_memory_budget_mb`), and reduces only
`Ne x Ne`-sized matrices (`cross`, `gram`, `rhs`) via `Allreduce`. No rank
ever holds a complete `Nx x Ne` ensemble matrix.

## Consequences

- This design is confirmed correct and scalable by direct code inspection
  (see the scalability audit) and should be preserved, not rewritten, as
  the rest of ICESEE's state handling is brought up to the same standard.
- It sets the bar new shared infrastructure (model-adapter state access,
  forecast-loop I/O) must be measured against: bounded per-rank memory,
  independent of global `Nx`.
- Modes 0 and 1 do not yet meet this bar (see ADR 0004) — this is a known,
  separate gap, not an argument against this decision.
- The opt-in `local_analysis` path is a partial exception: it gathers
  observation-space terms to rank 0, bounded by an explicit memory guard.
  It remains a single-rank ceiling and should be revisited if/when very
  large `Ny` with localization becomes a real workload.
