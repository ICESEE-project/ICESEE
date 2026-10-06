# Model-adapter architecture: conceptual target vs. implemented API

This document distinguishes what ICESEE's model-integration boundary
**conceptually should guarantee** from what is **actually implemented
today**. See ADR [0002](architecture/0002-model-adapter-boundary.md) for the
decision record, and the
[Wiki model-integration guide](https://github.com/ICESEE-project/ICESEE/wiki/3.--Guide-to-Integrating-Models-into-the-ICESEE-Framework)
for the concrete, current API contributors should follow today.

## Implemented today (do not break this)

Every model integration provides an example-specific `_<model>_enkf.py`
with a fixed set of hook functions, called with the single flat
`icesee_kwargs` runtime context:

| Hook | Purpose |
|---|---|
| `initialize_ensemble(ens, **icesee_kwargs)` | Construct/load one prior member |
| `forecast_step_single(ensemble=..., **icesee_kwargs)` | Advance one member one forecast interval |
| `generate_true_state` / `generate_nurged_state` | Reference / unassimilated trajectories |
| `inverse_step_single` | Model-specific deterministic inversion (hybrid workflows) |
| `post_analysis_update` | Enforce model constraints after analysis |
| `Obs_fun`, `JObs_fun`, `Cov_Obs_fun` | Observation operator overrides |

Shared ICESEE code (`src/EnKF/`, `src/parallelization/`, `src/utils/`)
never imports a model-specific library; it only calls these hooks by name,
resolved through `applications/supported_models.py`'s `MODEL_CONFIG`.
This is real, working model-agnosticism today for four models (icepack,
issm, lorenz96, flowline_model) and must be preserved exactly as-is by any
architectural evolution — no breaking change to this contract without a
migration plan for all four.

## Conceptual target (proposed, not yet fully realized)

The audit's core finding is that the *hook contract* is model-agnostic, but
the *state representation crossing it* is not yet uniformly so: several
models (most visibly `idealized_pig`) move whole per-member numpy state
vectors across the boundary, and nothing in the contract currently
distinguishes "give me the whole member's state" from "give me only what
this rank owns." The conceptual target, expressed as capabilities rather
than a mandated class hierarchy:

```
ICESEE orchestration (DA drivers, EnKF_parallel_io.py)
        |
        v
  model adapter  (per-model, e.g. applications/icepack_model/...)
        |  - initialize(comm, config)
        |  - advance(state, dt, ...)              [wraps forecast_step_single]
        |  - local_state(state) -> buffer/view     [NEW: no global array]
        |  - set_local_state(state, array) -> None [NEW]
        |  - state_layout() -> ownership ranges     [NEW]
        |  - observation_operator(state, ...)
        |  - checkpoint(state, path) / restore(path)
        |  - diagnostics(state)
        |  - finalize()
        v
  underlying model (Firedrake/PETSc, ISSM/MATLAB, plain numpy, ...)
```

The three "NEW" capabilities are what today's contract lacks: a way for
ICESEE's analysis code to address a rank's own state slab without the
adapter (or ICESEE) ever building a full `(Nx, Ne)` or even full-`Nx`
per-member array. This is the same principle ADR 0001 already proves works
at the analysis layer; the gap is only at the model/forecast layer.

## Explicitly deferred / out of scope for now

- A formal Python ABC/Protocol class is **not** being introduced in this
  stage. Per the agreed implementation strategy, no large abstraction layer
  should land before a second concrete adapter (ISSM, or a second
  significantly different Icepack workload) exercises it — premature
  interfaces designed against one model tend to be wrong for the second.
- Non-MPI models (plain numpy, e.g. Lorenz-96) already satisfy the target
  trivially (`local_state` is just the whole array, `ranks_per_model=1`);
  no special-casing is anticipated for them.

## ISSM extensibility check

ISSM already has its own multi-rank-per-solve concept (its `params.yaml`
sets `model_nprocs` explicitly, unlike icepack). The proposed
`ranks_per_model` communicator hierarchy (ADR 0003) generalizes exactly
this existing ISSM convention into shared infrastructure rather than
inventing something new — ISSM should be able to adopt the hierarchical
communicator work without a second redesign. The open question to validate
when ISSM integration work resumes: does ISSM's MATLAB-coupled state
representation support a "give me only my local slab" operation
comparable to a Firedrake `Function`'s `.dat.data`, or does the MATLAB
bridge only expose whole-state transfer? This must be answered from ISSM
integration code, not assumed, before `local_state`/`set_local_state` are
implemented as mandatory (vs. optional-with-fallback) adapter capabilities.
