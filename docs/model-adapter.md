# Model integration contracts

This page documents the inspected `features_debug` implementation at commit `3748e416d3a8120cc860119bd7224c6f9ad41ff7` (10 November 2026). Read the [detailed integration guide](handbook/model-integration.md), [API index](handbook/api-reference.md) and [accepted native-runtime decision](architecture/0006-native-runtime-and-time-contract.md) before adding an adapter.

## Legacy hook interface

Modes 0–2 resolve model-specific modules through `applications/supported_models.py` and pass the flat `icesee_kwargs` context. Hooks include `initialize_ensemble`, `forecast_step_single`, truth/unassimilated reference generators, `post_analysis_update`, optional `inverse_step_single`, and observation overrides. The precise spelling, arguments and required results depend on the example; there is no universal class hierarchy to substitute blindly.

State ownership must be explicit. `src/utils/state_ownership.py` distinguishes replicated full state from distributed disjoint pieces. The default distributed interpretation preserves existing model conventions; replicated examples must opt in. Ghost entries are not additional owned state rows. Preserve existing hooks when adding native capabilities.

## Native mode 3

Native mode 3 is implemented and selected through `src/parallelization/distributed_mode3_registry.py`. Its interfaces are in `distributed_adapter.py`, `distributed_fields.py`, `distributed_native_runtime.py`, and `distributed_native_cycle.py`. Native contexts retain model fields between forecasts, expose compact owned snapshots at analysis/checkpoint boundaries, and write analyzed owned slabs back to the native fields. Segmented block layouts describe variable-major global ordering without allocating the gaps between owned pieces.

Spatial communicators join one model group. Ensemble communicators join corresponding spatial positions across model groups. Groups may schedule several stable member IDs. Observation batches require aligned canonical IDs, kinds and ordering; the runtime fingerprints them before collective analysis. A grouped local target combines its accepted observation kinds in one transform.

Initialization, forecast, observation evaluation, projection, inversion, restore and diagnostic hooks are capabilities. Optional capabilities must be advertised only when the application actually implements them. The generic native application runner has narrower restrictions than specialized PIG; registering an adapter does not establish spatial, restart, generated-error, hybrid or scaling support for every combination.

## Implemented support and limits

The specialized Idealized PIG runner wires native Firedrake fields, process noise, bounded analysis, checkpoint and history publication. It is directly selected with `--execution_mode=3`; an extra experimental opt-in is unnecessary at this revision. Read its README for input files, launch scripts and measured validation scope.

Shared numerical analysis should operate on model-independent field/observation contracts. The full codebase still has historical model-specific branches and compatibility paths; a blanket assertion that all shared modules never import or recognize model-specific code would be inaccurate. Keep new model knowledge at the adapter boundary and add a capability rather than extending scattered model-name tests.

Native physical state need not be gathered for PIG forecast/analysis, but its canonical noise bridge retains a guarded global-size stochastic workspace. The tiny Lorenz native example intentionally gathers a small legacy snapshot. Do not advertise a universal zero-global-allocation guarantee.

## Adding or changing an adapter

1. Specify variables, units, owned entries, ghosts and the replicated/distributed contract.
2. Implement and validate exact per-variable ownership with no overlap or gap.
3. Keep member IDs independent of rank placement and test random policy across topology changes.
4. Preserve the selected factor/error/localization/inflation mathematics when changing execution.
5. Respect state-column timing: initial state is column zero; cycle k writes column k+1.
6. Add capability refusal tests, numerical row/product parity tests and real-MPI participation tests.
7. Test model-specific projection and a supported checkpoint replay across an analysis event.
8. Document solver memory, stochastic workspaces, stores and history costs separately.

The [testing guide](handbook/testing-and-contributing.md) explains which tests establish metadata logic, actual MPI behavior or backend integration. The [operations guide](handbook/operations.md) explains run isolation and recovery.
