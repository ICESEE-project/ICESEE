# Performance instrumentation plan

This document is a design/plan for lightweight, switchable timing and
memory instrumentation, so that later strong/weak/memory scaling studies
(not yet performed — see `scripts/benchmarks/README.md`) don't require
invasive code changes to add. No scaling campaign is run as part of this
document; it defines what will be measured and how.

## Phase decomposition

```
T_total = T_init + T_forecast + T_analysis + T_IO
          + T_sync_comm + T_diagnostics
```

Mapping to current code (mode 2, the most instrumented path today):

| Phase | Current code | Status |
|---|---|---|
| `T_init` | `initialize_model` (mesh, solver, checkpoint load), `icesee_mpi_ens_distribution` | not separately timed today |
| `T_forecast` | `parallel_forecast_step_default_full_parallel_run`, `_mpi_forecast_functions.py` | partially timed (`time_forecast_ensemble_generation`, `time_forecast_noise_generation`, `time_forecast_file_writing` already exist in `icesee_kwargs`) |
| `T_analysis` | `EnKF_parallel_io.compute_analysis_update` / `compute_X5_utils_` | partially timed via existing driver-level timing block (`icesee_da_full_parallel.py:798-871`) |
| `T_IO` | HDF5 batch open/close (`_ensure_batch`, `_create_batch_serial/_parallel`), checkpoint/prune | not separated from forecast/analysis time today — currently folds into whichever phase calls it |
| `T_sync_comm` | `Barrier()`/`Allreduce`/`bcast` wait time | not measured today; requires wrapping collectives, not just timing around them, to separate "waiting for others" from "doing the communication" |
| `T_diagnostics` | flowline-profile-style per-application diagnostics | not measured today; also the phase most exposed to a bug like audit finding D#1 (wrong-frequency execution), so measuring it independently is high-value |

## Design requirements

- **Low, controllable overhead.** A timing wrapper around each phase
  (`MPI.Wtime()` before/after, accumulated into `icesee_kwargs`, exactly as
  the existing `time_forecast_*` counters already do) rather than
  fine-grained per-line profiling. Full profiler-based instrumentation
  (e.g. py-spy, PETSc's own logging) is a separate, opt-in path for
  targeted investigation, not always-on.
- **Switchable.** A single config flag (e.g. `enkf-parameters.profile: true`,
  matching the existing YAML-driven config style) should enable
  per-rank timing/memory collection without code changes at call sites that
  don't need it.
- **Model-agnostic.** Instrumentation lives in the shared driver/DA layer
  (`icesee_da_serial.py`, `icesee_da_partial_parallel.py`,
  `icesee_da_full_parallel.py`) and the shared I/O layer
  (`EnKF_parallel_io.py`), not duplicated per application.
- **Machine-readable output**, suitable for later scaling plots without a
  reformatting step: one row per rank per run, e.g. CSV or JSON Lines
  alongside the existing checkpoint/diagnostic output, with columns for
  each phase's cumulative time, `Ne`, `P`, `ranks_per_model` (once ADR 0003
  lands), `Nx`, `Ny`, and peak memory where captured.
- **Peak memory per rank** via `resource.getrusage(RUSAGE_SELF).ru_maxrss`
  (already a dependency-free stdlib option; `psutil`, already a runtime
  dependency per `pyproject.toml`, gives a cheaper live-sample alternative)
  at phase boundaries, not continuous sampling.
- **Rank idle time** requires wrapping collectives themselves (time spent
  *inside* `Barrier()`/`Allreduce` before all ranks arrive is not
  distinguishable from "the collective's own cost" without a paired
  timestamp scheme) — flagged as the one piece of this plan that needs a
  small amount of new wrapper code, not just accumulation of existing
  `MPI.Wtime()` calls.

## What already exists and should be reused, not replaced

`icesee_da_full_parallel.py` already threads several `time_*` accumulators
through `icesee_kwargs` and has an end-of-run rank-0 timing display/report
step (referenced in the audit's fix history as "centralized timing/
performance display for modes 0 and 3"). The instrumentation plan above
extends that existing mechanism (finer phase separation, `T_IO` and
`T_sync_comm` split out, machine-readable per-rank output) rather than
introducing a parallel timing system.

## Explicitly not part of this document

Formal strong/weak/memory scaling experiments, and any specific numeric
performance claims — see `paper/results/README.md` for where measured
results will be recorded once instrumentation lands and experiments are
actually run.
