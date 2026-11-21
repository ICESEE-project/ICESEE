# Integrating and maintaining model adapters

For concrete model walkthroughs, read [Icepack coupling](coupling-icepack.md) and [ISSM coupling](coupling-issm.md). The [application onboarding procedure](application-onboarding.md) and [copyable scaffolds](templates/README.md) connect the contracts below to a new application.
ICESEE has an established generic hook API and a newer native field API. Select the boundary that matches what the model can actually expose. A MATLAB interface returning a complete member is a whole-member boundary even if the solver uses many internal ranks. A Firedrake field can expose owned degrees of freedom and halo operations, enabling a spatial native boundary. The documentation should describe these facts without forcing all models into one abstraction.

## Generic hook resolution

[`SupportedModels.call_model`](../../applications/supported_models.py) recognizes model families `icepack`, `issm`, `lorenz`, and `flowline`. From the application's working-directory basename it first tries `<model>_model.examples.<example>._<model>_enkf`, then the family-level module in `MODEL_CONFIG`. Names must be valid Python identifiers. An import failure inside an existing example module is a real failure, not a reason to run a different model quietly.

This makes the entry script and working directory part of the integration contract. Most example drivers change cwd to their own directory before loading settings. Avoid launching a hook file directly; use the example entrypoint and an absolute configuration path. Choose output paths according to the backend: the generic array examples accept absolute output directories, while the current MATLAB/ISSM bridge requires a relative `data_path` beneath its isolated application workspace. Extending `MODEL_CONFIG` is not enough: you also need a runnable example, physics/input initialization and hooks expected by the selected runner.

| Hook or operation | Required behavior to verify against the runner |
|---|---|
| `initialize_model(**icesee_kwargs)` | Application-specific setup, dimensions, model objects, files and communicator placement; entry scripts may call it before dispatch. |
| `initialize_ensemble(ens, **icesee_kwargs)` | Build/load one prior using explicit member identity. Return block values in the established variable order and distinguish arrays from scalar blocks. |
| `forecast_step_single(ensemble=..., **icesee_kwargs)` | Advance exactly one interval using current member/cycle, return/update the expected variable blocks, and respect the model communicator. |
| `generate_true_state`, `generate_nurged_state` | Produce reference trajectories on the correct initial-state clock. Honor file-backed targets where passed; write in place and return `None` rather than materializing history. |
| `inverse_step_single` or application inverse callback | Only required for the selected hybrid workflow; document inversion inputs, output parameters and internal process/server lifecycle. |
| `post_analysis_update` | Restore model geometry/diagnostic consistency after a mathematical update. Do not silently alter unrelated parameters. |
| `Obs_fun`, `JObs_fun`, `Cov_Obs_fun` | Optional operator overrides where consumed; nonlinear `H` evaluates each member before centering. Check exact signatures at call sites. |

Historical adapters do not use a single formal class for every generic hook and some signatures/return conventions differ. Copy the closest working example only after following its active callers. A hook table is a navigation guide, not permission to replace existing function signatures. Use focused tests to pin the exact integration you add.

## Generic array invariants

The same ordered `vec_inputs` must describe initialization, forecast results, observation indexing, analysis, and output. The number of blocks and dimensions must agree with `num_state_vars`, parameter metadata, `nd`, and any `var_nd`. A shape-compatible state with swapped blocks is scientifically wrong and can escape finite-value checks. Make a tiny state with distinguishable values per block to test this mapping.

Declare `state_distribution: replicated` when model-communicator ranks report identical whole arrays. Declare distributed ownership only for disjoint owned pieces. Ensure a single-element block is represented consistently; `ensure_state_array` normalizes scalars where used. Never gather ghost copies as additional owned state.

In a multi-rank model callback, all ranks in the model communicator must enter solver collectives in the same order. Global observation/truth files need a designated writer. Reusing world rank in place of member ID fails when multiple members are scheduled per slot. Initialization, forecasting, process noise and file naming must all use stable member identities.

## Native fields and persistent members

The active native lifecycle is defined in [`distributed_native_runtime.py`](../../src/parallelization/distributed_native_runtime.py). `NativeDistributedMember` contains a stable `member_id`, a `DistributedFieldRegistry`, and an opaque model context. The model owns persistent vectors for the member's lifetime. ICESEE packs only owned entries for analysis/checkpoints and writes analyzed owned values back into those same vectors.

`DistributedFieldHandle` provides the field's name/layout, owned read/write access and synchronization. The registry validates field compatibility, combines field blocks into an ordered layout, and packs/unpacks local owned values. Native code must not create a whole physical member merely to satisfy this interface. Ghosts belong to model synchronization; the packed analysis state excludes them.

`NativeDistributedModelAdapter` requires these callbacks:

| Callback | Responsibility |
|---|---|
| `initialize_native_member(member_id, *, topology, icesee_kwargs)` | Construct a persistent member on `topology.spatial_comm` and return its field registry/context. |
| `forecast_native_member(member, timestep, *, topology, icesee_kwargs)` | Mutate the same native member in place, including native halo exchange and solver collectives. |
| `observe_native_member(member, observation_rows, *, topology, icesee_kwargs)` | Evaluate only locally owned observation rows and preserve canonical IDs/ordering. |
| `finalize_native_analysis(member, forecast_owned, timestep, *, topology, icesee_kwargs)` | Restore physical/diagnostic constraints after analyzed owned values have been unpacked. |
| Optional `inverse_native_member(...)` | Apply inversion collectively inside the model's spatial communicator without returning a gathered global parameter vector. |

`NativeObservationShard` maps local observation rows back to canonical active positions. This is necessary when mesh partition order differs from the serial observation IDs. PIG maps physical coordinates to distributed rows before routing observations or increments. Numerical row IDs from independently constructed meshes must not be assumed identical just because array lengths agree.

The lower-level `DistributedModelAdapter` protocol in `distributed_adapter.py` operates on packed local arrays (`initialize_local_member`, `forecast_local_member`, etc.). It remains useful for reference and synthetic tests. It is distinct from the persistent native lifecycle; do not conclude that implementing one protocol automatically registers a production policy for the other.

## Localized native analysis metadata

An adapter supporting spatial local analysis supplies owned target rows, stable global rows, physical coordinates and allowed observation kinds. Observation metadata supplies canonical IDs and coordinates per kind. Validators check local/global row mapping, uniqueness, dimensions and referenced kinds without assembling full physical fields.

A target can use multiple observation kinds in one joint collection. Keep the active ordering stable, including sparse/missing observations. Applying separate per-kind transforms sequentially is not generally identical to one joint update. Build a small reference case where the same target sees multiple kinds, then compare full member values with the intended generic policy.

## Registration and explicit capability limits

`register_execution_mode_3(model_name, runner, application=...)` binds a specialized runner. `register_native_application` binds model setup and a native adapter to the shared lifecycle, with an optional analysis policy and declared `supported_features`. The shared default policy uses one spatial rank and rejects requested restart or unsupported scientific features. PIG has a specialized spatial/restart policy; it must not be replaced by that default without a migration and parity study.

Registration occurs during import; dispatch also tries an example-specific bootstrap before reporting an unsupported model. Tests should cover both explicit entrypoint registration and bootstrap when relevant. The registry scopes applications so one Icepack example's runner is not accidentally used for another.

Before exposing a new adapter, publish a matrix covering forecast representation, supported error factors, localization, hybrid inversion, inference plugins, restart, precision, full-history export and allowed topology. Mark untested combinations as unvalidated. If the boundary is whole-member, state the actual memory bound and avoid a strict no-global-member claim.

## Acceptance sequence for a new model

1. Run a deterministic model-only forecast and establish expected units, shape, geometry and boundary behavior. Use no assimilation to isolate model setup.
2. Test initialization with at least two distinguishable members and multiple variables. Verify member IDs, block order and partition coverage.
3. Test observations using known state values, sparse IDs and, if used, a nonlinear operator. Check canonical positions and physical coordinate mapping.
4. Run a tiny analysis against an explicit reference with identical realized observations/noise. Compare members, not just means; compensating errors can leave the mean unchanged.
5. Run actual MPI with at least two groups and, for a distributed model, two spatial ranks per group. Test unequal member/group counts, empty local observation sets and failure on one rank.
6. Interrupt after a committed checkpoint, resume with the same layout and policy, cross an analysis event and compare final values to an uninterrupted reference. A restart that only forecasts further without crossing analysis does not test event bookkeeping.
7. Verify full/rolling storage semantics and postprocessing. Retain a deliberate failing output so the verifier is known to reject wrong dimensions, a missing shard or a changed observation clock.
8. Only then run a representative physical workload and obtain a scientific acceptance decision. Include a resource and validation protocol before large ensembles.

Each stage should preserve logs and a case configuration. Unit tests establish interface behavior; MPI tests establish communicator behavior; model validation establishes physical meaning. None substitutes for the others.
