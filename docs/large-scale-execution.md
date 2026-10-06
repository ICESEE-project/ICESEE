# Large-scale ensemble execution

ICESEE execution mode 2 (`execution_mode: 2`) distributes model forecasts and
streams the EnKF analysis. It is intended for models whose state and ensemble
history do not fit comfortably in memory.

## Recommended configuration

```yaml
enkf-parameters:
  execution_mode: 2

  # Select `full` only when every ensemble timestep must be retained.
  # `auto` chooses rolling storage when the projected full history is too large.
  ensemble_history_mode: auto       # auto | full | rolling
  batch_size: 2                     # active timestep-file window

  # Zero means derive a row chunk from the memory budget.
  analysis_memory_budget_mb: 256
  analysis_row_chunk_size: 0
  observation_row_chunk_size: 0

  # Avoid dense observation operators and dense synthetic-observation arrays.
  observation_operator_storage: indices
  synthetic_observation_storage: compact

  # Disk storage only; analysis arithmetic remains float64.
  ensemble_storage_dtype: float32   # use float64 when required by the experiment
  storage_safety_factor: 1.25
  fail_on_insufficient_storage: true

  restart_enabled: true
  force_fresh_start: false
  ensemble_finalize_mode: vds       # full history only; skipped in rolling mode
```

## What is bounded

The mode-2 analysis no longer constructs a dense observation operator, a dense
state-by-ensemble analysis matrix, or an observation-by-observation covariance.
It streams observed and state rows and reduces the analysis to `Nens x Nens`
Gram and right-hand-side matrices. Synthetic observations are stored as compact
observed rows rather than as a dense state-sized array.

`ensemble_history_mode: rolling` keeps only the current restartable ensemble
state and a short active write window. Completed ensemble shards are pruned only
after a checkpoint identifying the retained shard has been written. The ensemble
mean time series remains available, but historical member spread cannot be
reconstructed after pruning. `auto` retains full history when projected storage
fits safely and otherwise selects rolling mode.

Mode 2 also passes file-backed true and nudged trajectory targets to applications.
Large application adapters should write those targets in place and return `None`;
returning slices of the full target reconstructs the trajectory in RAM. The
Icepack PIG, Icepack synthetic-stream, and ISSM adapters implement
this contract.

## Capacity model and hard lower bound

With `b` bytes per stored value, a full ensemble history requires approximately

```text
state_dofs * ensemble_members * time_steps * b
```

whereas rolling history is independent of run length and requires approximately
one ensemble snapshot plus the active write window:

```text
state_dofs * ensemble_members * b
```

This snapshot is a physical lower bound for a conventional EnKF. If one model
state is genuinely 37 GB, a 40-member snapshot represents about 1.48 TB before
compression. Rolling history prevents that cost from being multiplied by the
number of timesteps, but it cannot remove the snapshot itself. Such applications
also require one or more of:

- a reduced assimilation state or reduced-order ensemble representation;
- model-native spatial decomposition so no rank holds a complete member;
- validated lower-precision or compressed disk storage;
- fewer simultaneous members or hierarchical/serial member processing.

Choose MPI rank placement from the memory required by active model members, not
only from core count. If one member does not fit on its assigned ranks, the model
adapter must provide domain-decomposed state access; the generic ICESEE I/O layer
cannot compensate for that model-side limit.

## Restart and output behavior

- Rolling mode checkpoints every completed cycle and retains the matching state
  shard. Restart resumes from that shard without repeating its forecast.
- `force_fresh_start: true` intentionally discards reuse of an existing run.
- Full mode can create a zero-copy HDF5 virtual dataset (`vds`) over timestep
  shards. Rolling mode skips this because historical member shards no longer
  exist.
- Use `ensemble_history_mode: full` for experiments requiring member-wise plots
  at every historical timestep, after checking the reported storage estimate.

## Mode-1 parity contract

Execution mode 2 is an alternative storage and distribution strategy, not a
different data-assimilation method. Given the same configuration, observations,
initial ensemble, member identifiers, and `base_seed`, it must reproduce the
mode-1 ensemble at every synchronized timestep. This includes the complete
member ensemble, rather than only the ensemble mean or plotted diagnostics.

The parity contract currently covers:

| Capability | Mode 1 | Mode 2 |
|---|---:|---:|
| Reproducible member initialization independent of MPI placement | yes | yes |
| Reproducible forecast/process noise by member and timestep | yes | yes |
| Global stochastic EnKF analysis | yes | yes, streamed |
| Generated spatial observation-error fields | yes | yes, column-streamed |
| Legacy forecast-anomaly analysis factor | yes | yes, streamed |
| Grouped local analysis and covariance localization | yes | yes, streamed |
| Inflation, frozen variables, relaxation, bounds, and increment caps | yes | yes |
| Sparse/compact observations and bed-observation masks | yes | yes |
| Delayed hybrid inversion and member-wise inversion handoff | yes | yes |
| Checkpoint/restart and explicit restart timestep | yes | yes |
| Full member history | yes | yes, sharded/VDS |
| Run-length-independent rolling history | no | yes |

Use the paired-run harness to create both runs from one configuration and compare
every value without loading either history into memory:

```bash
python scripts/benchmarks/run_execution_mode_parity.py \
  --config applications/lorenz_model/examples/lorenz96/params.yaml \
  --application applications/lorenz_model/examples/lorenz96/run_da_lorenz96.py \
  --output /tmp/icesee-mode-parity \
  --ranks 2 --nens 8
```

Existing runs can be checked independently:

```bash
python scripts/benchmarks/compare_execution_modes.py \
  --partial /path/to/mode1 \
  --full /path/to/mode2
```

The default acceptance gate is deliberately strict numerical parity
(`atol=1e-10`, `rtol=1e-8`). MPI reductions in the two engines need not occur in
the same floating-point order, so bitwise identity is not portable even when the
algorithms are equivalent. Use `--bitwise` only as a fixed-platform diagnostic.
Float64 storage is required for this strict gate. A validated lower-precision
run should instead declare scientifically justified tolerances.
The two runs must use fresh output directories; comparing a restarted run to a
fresh run is a separate restart-equivalence test. New mode-2 optimizations should
not be accepted until this complete-history gate passes for a small deterministic
case and for the target application configuration.

The implementation status, invariants, and acceptance matrix for continued
mode-2 development are maintained in
[`execution-mode-2-development.md`](execution-mode-2-development.md).

## Architecture and future scaling

The current design—member-parallel forecasts, file-backed timestep shards, and a
row-streamed ensemble-space analysis—is an efficient and robust baseline when a
single model member fits on its assigned model ranks. It avoids an in-memory
`state_dofs x Nens x time_steps` history and reduces global EnKF algebra to the
ensemble dimension.

It is not the final optimal architecture for arbitrarily large ice-sheet
problems. The next scaling layer is tracked as proposed `execution_mode: 3` and
should preserve model-native spatial decomposition with a two-dimensional
process topology:
one communicator dimension over ensemble members and another over spatial
subdomains. Analysis should then operate on distributed local patches without
assembling a complete member on any rank. Parallel HDF5 or ADIOS2, asynchronous
checkpointing, topology-aware placement, and optional reduced-rank or multilevel
ensembles are appropriate follow-on improvements. Each changes performance or
the numerical representation and therefore must be introduced behind the parity
gate and validated independently.

See [`execution-mode-3-design.md`](execution-mode-3-design.md) for the proposed
distributed adapter contract, stochastic-analysis semantics, phased rollout,
and acceptance gates. Mode 3 is not yet a selectable runtime mode.
