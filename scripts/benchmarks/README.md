# Benchmark structure and conventions

This describes how ICESEE benchmarks are organized so that model/spatial
scaling, ensemble scaling, DA/analysis scaling, I/O scaling, and memory
scaling can eventually be characterized **separately**, not collapsed into
one aggregate wall-time number. No formal scaling campaign has been run yet
(see `paper/results/README.md`); this is the structure that campaign will
use.

## Existing scripts (already in this directory)

- `run_execution_mode_parity.py`, `compare_execution_modes.py` — correctness
  parity across modes 0/1/2, not performance benchmarks.
- `run_mode3_lorenz_forecast_parity.py`, `run_mode3_checkpoint_parity.py`,
  `run_mode3_block_checkpoint_parity.py`, `run_mode3_local_analysis_parity.py`,
  `run_mode3_local_runtime_parity.py`, `run_mode3_restart_parity.py` —
  mode-3 scaffolding parity gates (mode 3 has no registered production
  application yet).
- `run_mode3_memory_scaling.py` — the one existing memory-scaling script;
  currently mode-3-primitives-only.
- `../ci/run_lorenz96_ci.py` — the CI smoke test (correctness, not timing).

## Benchmark tiers (planned)

### Small — Lorenz-96

Rapid regression and MPI-correctness testing; cheap enough to run on every
PR if needed. Already has working infrastructure
(`applications/lorenz_model/examples/lorenz96/`, `run_lorenz96_ci.py`,
`run_mode3_lorenz_forecast_parity.py`). This tier is the right place to
validate ADR 0003's `ranks_per_model`/block-assignment change first, at
many `Ne`/`P` combinations, before touching a real model.

### Intermediate — a multidimensional model cheap enough for repeated runs

Not yet selected. Candidate: a trimmed/low-resolution `idealized_pig`
configuration (short `num_years`, coarse mesh) or `flowline_model`, small
enough to iterate on in development but exercising real multidimensional
state and real Firedrake/PETSc machinery (unlike Lorenz-96). This tier is
where audit findings D#4/D#6 (intra-member spatial parallelism,
node-locality of rank assignment) should be measured in isolation from
total run length.

### Realistic — Icepack `idealized_pig` (current, ~37 GB input dataset)

The current stress/real-world validation case — not the architectural
upper bound. This is where the audited 8-day runtime was observed and
where before/after comparisons against that baseline belong, once Phase 1
fixes land.

### Future large-scale — continental-scale domains (e.g. Antarctica), ISSM

Not attempted yet. Requires the hierarchical communicator work (ADR 0003)
and, for ISSM, the extensibility questions in `docs/model-adapter.md`
resolved first.

## What each benchmark run should record

See `paper/results/README.md` for the exact ledger format. In short:
machine/environment, `Ne`, `P`, `ranks_per_model`, `Nx`, `Ny`, dataset size,
per-phase timing (`T_init/T_forecast/T_analysis/T_IO/T_sync_comm/T_diagnostics`
per `docs/performance-instrumentation.md`), peak memory per rank, and the
ICESEE commit/version tested. Do not report a single wall-clock number
without this context — it isn't reproducible or comparable across the
tiers above.

## Explicitly not yet done

No strong-scaling, weak-scaling, or memory-scaling sweep has been executed
as part of this planning pass. Architecture stabilization, numerical-
equivalence validation, memory/communication/I-O fixes, and instrumentation
come first, per the agreed implementation strategy.
