# 0002 — Model-adapter boundary between orchestration and model internals

## Status

Proposed. The current implemented mechanism is the function-hook contract
described in Wiki page 3 (`initialize_ensemble`, `forecast_step_single`,
`generate_true_state`, `generate_nurged_state`, `inverse_step_single`,
`post_analysis_update`, `Obs_fun`/`JObs_fun`/`Cov_Obs_fun`, all called with
`**icesee_kwargs`). This ADR proposes the conceptual target the hook
contract should keep evolving toward as state gets larger and more model
backends (ISSM, future MPI and non-MPI models) are added. See
`docs/model-adapter.md` for the full conceptual interface.

## Context

Shared ICESEE code must not import Firedrake, ISSM/MATLAB bindings, or any
other model-specific library. Today this is upheld by convention (model
code lives under `applications/<model>/`, shared code under `src/`) but the
hook contract itself is somewhat Icepack-shaped in places — e.g. it assumes
a flat numpy state vector can be produced/consumed by the model side, and
does not yet have an explicit "give me this rank's locally owned state
without materializing the global vector" hook.

## Decision (proposed)

Formalize the boundary as: ICESEE orchestration (DA drivers,
`EnKF_parallel_io.py`) calls only through adapter-level operations —
initialize, advance forecast, get/set locally owned state, observation
access, checkpoint/restart, diagnostics, finalize — and never assumes how
the model represents state internally (Firedrake `Function`, ISSM/MATLAB
structures, plain numpy, or something else). The current hook contract is
the concrete Python-level implementation of this boundary today; the
proposal is to make the "local, not global, state" part of it explicit and
consistent across all models, rather than something each model integration
reimplements independently.

## Consequences

- Existing model integrations (icepack, issm, lorenz96, flowline_model)
  keep working under the current hook contract; this is an incremental
  evolution, not a breaking API change, per the implementation strategy
  agreed for this work.
- ISSM must be checked against this boundary explicitly as it's designed
  (see ADR 0002's companion question in `docs/model-adapter.md`: "can ISSM
  use this without redesigning ICESEE again?").
- No large abstraction layer should be implemented under this ADR until a
  concrete second consumer (ISSM, or a second Icepack workload) forces
  specific design questions — see `docs/model-adapter.md` for what's
  deliberately left open.
