# 0003 — Hierarchical ensemble/model MPI communicator topology

## Status

Proposed. The scalability audit found the currently implemented topology
(`ParallelManager.icesee_mpi_ens_distribution`,
`src/parallelization/parallel_mpi/icesee_mpi_parallel_manager.py`) conflates
ensemble-count and rank-count into a single `color = rank % X` split, which
means intra-member (spatial) parallelism is an accident of the `P`/`Ne`
ratio rather than a design parameter, and rank assignment is currently
strided rather than block/node-aware.

## Context

Ensemble-level parallelism (how many members advance concurrently) and
spatial/model-level parallelism (how many ranks cooperate on one member's
PDE solve) are logically independent, but today's split makes one
determine the other. With few ensemble members (e.g. `Ne=4`) this can mean
either serialized rounds (`Ne >= P`) or, if `P` is large enough to give
each member its own subcommunicator, a strided rank assignment that
fragments a member's domain-decomposed mesh across compute nodes given
typical `ntasks-per-node` launch configurations.

## Decision (proposed)

Make `ranks_per_model` an explicit, user-suppliable (or auto-derived)
parameter, with `P = Ne x ranks_per_model` as a designed invariant
(including documented, tested behavior for the non-divisible case), and
replace strided (`rank % Ne`) rank assignment with block-contiguous
assignment so each member's ranks stay together for node locality.

## Consequences

- This changes `icesee_mpi_ens_distribution`'s internals but not its
  external contract (`comm_world`, `subcomm`, `color`, `rounds` are still
  provided to callers) — existing callers in `_mpi_forecast_functions.py`
  and mode 1's driver should not need to change, only be re-validated.
- Requires numerical-equivalence regression testing (Lorenz96 parity, then
  a trimmed `idealized_pig` run) before/after, since communicator topology
  changes can change MPI reduction order and therefore floating-point
  roundoff — acceptable within the tolerances already defined in
  `docs/execution-mode-2-development.md`, not silently beyond them.
- Does not by itself fix the per-timestep world `Barrier()`
  (`icesee_da_full_parallel.py:605`) or the `write_forecast` collective-
  batch-rebuild mismatch found in the audit — those are separate, related
  fixes tracked alongside this one.
