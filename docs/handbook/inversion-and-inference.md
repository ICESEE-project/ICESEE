# Inversion and parameter inference in ICESEE

We use **parameter inference** for several distinct operations in ICESEE. An ensemble filter can update a parameter through its sampled covariance with observed quantities. A model adapter can instead solve a deterministic inverse problem after an ensemble analysis. Optional post-analysis routines can also regularize bed increments or estimate an apparent mass balance from an analyzed history. These operations have different objectives, inputs and failure boundaries. Their coexistence does not make the resulting procedure an unchanged EnKF posterior.

This chapter follows the ISMIP_Choi ISSM application at scientific-source commit `3748e416d3a8120cc860119bd7224c6f9ad41ff7`, with the narrowly scoped missing-inverse-result guard described below as a later handover correction. The Choi-specific equations and handoff below describe inspected implementation; they are not validated by an executed Choi finite-element inverse. The separate basal-friction teaching exercise has its own [backend evidence](sample-runs.md). Python unit tests use explicit substitutes for the external solver. Read [ISSM coupling](coupling-issm.md) for workspace/server setup, [numerics](numerics.md) for the ensemble update, and [current developments](current-developments.md) for the detailed SMB/bed algorithms and control order.

## 1 Three meanings of an inferred parameter

Let a member contain state variables `x` and parameters `theta`. In joint filtering, we augment the state as `z=(x,theta)` and propagate each member through its configured forecast. The parameter rows of a Kalman-style analysis then have the form

```text
theta_a(i) = theta_f(i) + K_theta [y + eta(i) - H(z_f(i))]
K_theta = sample_cov(theta_f, H(z_f)) [sample_cov(H(z_f)) + R]^{-1}
```

This is a conceptual covariance expression; the implemented stochastic factor, truncation and observation-error policy are specified in the numerics chapter. It explains why an unobserved bed or friction field can change when correlated state observations are assimilated. Its usefulness depends on ensemble spread, identifiability and the observation experiment. Joint estimation is established ice-sheet assimilation practice, including [Bonan et al. (2014)](https://doi.org/10.5194/npg-21-569-2014), [Gillet-Chaulet (2020)](https://doi.org/10.5194/tc-14-811-2020) and [Choi et al. (2025)](https://doi.org/10.5194/tc-19-5423-2025); their experiments do not validate a new ICESEE application automatically.

The **hybrid friction workflow** is different. We first analyze the non-friction variables, then call an ISSM stress-balance control problem for each member and replace its friction coefficient. The deterministic optimizer is inside ISSM, rather than inside ICESEE's sample-covariance calculation. Control methods relate a model's unknown parameters to an observation-misfit objective; [MacAyeal (1993)](https://www.cambridge.org/core/journals/journal-of-glaciology/article/tutorial-on-the-use-of-control-methods-in-icesheet-modeling/CD1F93F5EBF26E7CDAADC376CA2BC6DE) provides foundational context. [Morlighem et al. (2010)](https://doi.org/10.1029/2010GL043853) illustrates why the inferred basal-drag pattern also depends on the chosen ice-flow approximation. We therefore retain the reference model, sliding law and dependency versions as part of an inverse experiment's identity.

Finally, the **inference plugin** applies configured post-analysis rules. Its bed helper smooths an increment proposed elsewhere; it does not solve the friction control problem. Its SMB helper uses thickness tendency and horizontal flux divergence, followed by regularization. Each has separate activation flags. Enabling a parameter in `joint_estimated_params` does not enable these helpers or an ISSM inverse callback.

## 2 Identify the coefficient and the state being exchanged

The Choi application orders six nodal blocks as

```text
vec_inputs = [Thickness, Surface, Vx, Vy, bed, coefficient]
z = (H, s, u, v, b, C)
```

`vel_idx: 2` identifies the first velocity block and `friction_idx: 5` identifies the last coefficient block, using zero-based Python block indices. All six blocks are nodal and equal-sized in this application. The generic Mode 1 reduction and global physical finalizer assume this layout; changing the coefficient position or inserting scalar states requires reviewing that code, rather than changing only a YAML list. Mode 2 has an explicit variable-size resolver, but its compatibility finalization still calls routines with their own equal-block requirements.

The prior-generation HDF5 file `friction_bed_<member>.h5` carries a **coefficient perturbation**. MATLAB initialization adds `mean_friction` to that perturbation. The initialized/forecast/analyzed `ensemble_output_<member>.h5` instead carries the **absolute** `md.friction.coefficient`. In the inverse branch, MATLAB reads this absolute value directly; it does not add `mean_friction` again. Confusing these two file roles shifts every initial guess.

We denote the packed quantity by `C`, matching ISSM's `FrictionCoefficient` control. It is not basal stress, a dimensionless drag multiplier, or automatically the squared quantity `C^2`. For example, `Mismip2.par` constructs its initial coefficient with a square root and the year's seconds conversion, and sets the friction exponents to `p=3`, `q=0`; forecast and inverse branches restore these same exponents. Units and the relationship between coefficient and traction belong to this selected ISSM friction law. Preserve the parameterization and ISSM revision when reporting coefficients; do not relabel the YAML bounds as pascals or transfer them to a different law without deriving the conversion.

At the MATLAB interface, velocity observations follow ISSM's metres-per-year convention; geometry is in metres. The official [`m1qn3inversion` class](https://raw.githubusercontent.com/ISSMteam/ISSM/main/src/m/classes/m1qn3inversion.m) documents the interface and performs velocity conversion during marshalling. This is separate from ICESEE's configured model clock and from internal SI representations. The weight values below are application settings, not an automatic nondimensionalization of all objective terms.

## 3 The actual ISSM objective and optimizer controls

For one member and one active observation time, we write the configured control problem as

```text
minimize over C: J(C) = w_abs J_101(C) + w_log J_103(C) + w_reg J_501(C)
subject to: C_min(j) <= C(j) <= C_max(j) at the control vertices
```

Here `J_101`, `J_103` and `J_501` are ISSM's finite-element response functionals, evaluated through its stress-balance solution. Their official identifiers are absolute surface-velocity misfit, logarithmic surface-speed misfit, and drag-coefficient gradient penalty. The [ISSM Jakobshavn control tutorial](https://issm.jpl.nasa.gov/documentation/tutorials/jks/) and [ISSM formulations guide](https://issm.jpl.nasa.gov/files/issm/userguide1.pdf) identify these responses. The last penalizes coefficient gradients rather than the magnitude of `C` itself. We keep the response identifiers in the equation because quadrature, masks, small-speed treatment, normalization and adjoint implementation are properties of the installed ISSM version; ICESEE does not implement those element integrals.

The YAML key `rel_vel_weight` multiplies **103**, the logarithmic response. It does not select 102, ISSM's relative-velocity response. This naming mismatch is material when interpreting sensitivity experiments. Neither the velocity weights nor the gradient weight are computed from ICESEE's observation covariance `R`. Thus the EnKF likelihood settings and inverse objective weights are independent scientific choices, even when both stages use the same synthetic observations.

| Control | Choi application value and meaning |
|---|---|
| Class and active control | `m1qn3inversion()`, `iscontrol=1`, `control_parameters={'FrictionCoefficient'}` |
| Responses | `[101 103 501]`, with one coefficient column per response at every vertex |
| Response weights | `abs_vel_weight=1`, `rel_vel_weight=1`, `tikhonov_regularization_weight=1e-13` in the base YAML |
| Bounds | `min_friction=1900`, `max_friction=4200`, expanded to nodal arrays |
| Optimizer limits | `maxsteps=40`, `maxiter=40` |
| Convergence controls | `dxmin=0.1`, `gttol=1e-4` |
| Stress-balance tolerances | `restol=0.01`, `reltol=0.1`, `abstol=NaN` |
| Solve | `solve(md,'Stressbalance','runtimename',false)` with ISSM `generic(...,'np',nprocs)` |

These are the active `run_model.m` inverse-branch assignments, not settings inferred from a similarly named `.par` file elsewhere. In the official class documentation, `maxsteps` limits gradient computations, `maxiter` limits function/forward evaluations, `dxmin` compares nearby control vectors in the supremum norm, and `gttol` uses the gradient norm relative to its initial value. Reaching a limit is not evidence that the desired tolerance was achieved. Retain ISSM's termination reason, objective components, gradient information and iteration history when evaluating an actual run.

The application does **not** explicitly assign `incomplete_adjoint`, `control_scaling_factors` or `dfmin_frac` after constructing this object. The official source inspected on 10 November 2026 initializes them to `1`, `1` and `1`, respectively; the first selects the class's linear-viscosity/incomplete-adjoint setting. This documents an external default, not a verified property of an unknown local ISSM build. Record the actual resolved inversion object and solver version before claiming an exact or complete adjoint. [Morlighem et al. (2013)](https://doi.org/10.1002/jgrf.20125) examines exact and incomplete adjoints in an ISSM inverse setting; it is context for that distinction, rather than a validation of this ICESEE wrapper.

## 4 Activation uses the observation clock

Generic Modes 1 and 2 evaluate inversion only when `observation_is_due` selects an analysis. The requested switch is retained as `inversion_enabled`; the effective per-cycle switch is recomputed as

```text
active_inverse = inversion_enabled and observation_time[km] + 1e-12 >= inversion_start_time
```

`resolve_analysis_cycle_time` uses requested `obs_t[km]` when available, and separately returns mapped `t[k+1]` for the log. It falls back to model time when no observation-time entry is available. This matters when an annual observation maps to a nearby state-grid time: a slightly earlier mapped time should not inadvertently defer a year-two inverse. `km` is an observation-column index, not the model timestep `k` and not time in years.

There is no separate inverse-frequency control in these inspected paths. Once its start time is reached, inversion follows every due analysis while enabled. The base Choi YAML requests observations from year 2 through year 40 with `freq_obs: 1` and sets `inversion_start_time: 2`. A different observation schedule changes the actual inverse cadence. A forecast-only step has no inverse call.

Before the start time, ordinary EnKF parameter controls still apply. After activation, generic Modes 1 and 2 reserve the friction update for the inverse stage. `frozen_analysis_vars: [coefficient]` can additionally hold that row during earlier EnKF analyses. Such a freeze is an **analysis-row** rule; it does not disable a later inverse callback, process noise, forecast-side evolution or an independently enabled global helper. Report both the requested flags and the per-cycle activity.

The `friction_experiments` profiles distinguish fixed-background friction, EnKF-updated friction and inversion-based friction. Their `extends` paths must be resolved before a launch. In the inspected checkout, `params_ibf.yaml` refers to `../param.yaml`, whereas the visible base configuration is `../params.yaml`; a similarly named `param.yaml` exists deeper in `method_comparison_40yr`. Do not treat these child profiles as a runnable comparison merely because their scientific switches are informative. Use a reviewed, fully resolved configuration and preserve that file chain.

## 5 Generic Mode 1: analyze five blocks, invert members, replace friction

At an active cycle, the partial-parallel runner saves the complete forecast ensemble to `ensemble_before_analysis_step_<stored-step>.h5`. It removes the last friction block and corresponding observations from the EnKF problem, analyzes the remaining five blocks, and calls the shared finalizer on that reduced state. The writer merges these analyzed/finalized blocks back into the saved full forecast. Its inverse input therefore has the form

```text
z_star(i) = (Finalize[EnKF(H,s,u,v,b)](i), C_forecast(i))
z_next(i) = (z_star_without_C(i), InverseISSM(z_star(i), observed_u_v))
```

This notation describes the handoff, not an assertion that MATLAB uses every geometric field before its solve. That distinction is detailed in Section 7.

The writer gathers the reduced row state on root, constructs the full inverse ensemble, and broadcasts the complete `(nd,Nens)` array to all participating ranks. Rank `r` processes member IDs `r, r+size, ...`; each Python inverse call receives `MPI.COMM_SELF`. The external ISSM calculation still uses the separately configured `nprocs`, derived from the application model-worker context. Python member parallelism and ISSM's internal parallel solve are distinct resource axes.

For each returned mapping, the current writer accepts a recognized friction field and overwrites only those coefficient rows. Its velocity-update branches are commented out. Returned `Vx`, `Vy`, thickness, surface and bed do not replace the already analyzed values here. `None` raises an error, and a non-mapping result is rejected. The base application uses the recognized key `coefficient`; aliases in this older writer do not remove the requirement that the returned key also map to a real configured block.

The output ensemble and mean are written after member updates have been gathered. This is not the same transaction design as Mode 2's temporary shard. A failure can occur after external member files have been written or workers have been killed; root-only output ordering does not roll those effects back. This path also retains full-ensemble memory during inverse preparation, even if forecasts were scheduled member by member.

## 6 Generic Mode 2: preserve the canonical forecast until success

Mode 2 retains a fixed full-state file layout. At an active hybrid cycle it excludes friction observation rows, temporarily freezes the friction analysis rows, and finalizes the same reduced non-friction state used by Mode 1. `finalize_hybrid_analysis_like_mode1` merges the finalized five-block result into the saved forecast, ensuring a pre-inversion EnKF coefficient update cannot leak into the inverse input. It restores the caller's configured freeze policy after the cycle and carries private `_bed_`/`_smb_` inference state back from the reduced mapping.

The complete new analysis remains in an unpublished temporary shard. Workers independently read one full member at a time from it, call `inverse_step_single`, and retain their returned coefficient vectors. They validate the result's mapping, recognized coefficient field and expected length. After all workers report success through an `allgather`, the updated coefficients are written to that temporary file. Root then publishes it with `os.replace(temp_path, canonical)`.

An inverse error therefore leaves the canonical forecast/restart shard unchanged. This is the precise output-publication guarantee exercised by the Python contract tests. It does not restore MATLAB's member HDF5 files, solver scratch, subprocess state, or any external side effect of a completed member solve. Those require the application to stage or recover its own files. Likewise, the optional finalizer may require a separate complete-ensemble memory budget; file-backed analysis is not a blanket memory bound on every enabled inference stage.

## 7 Follow the MATLAB initial guess, observations and returned fields

`inverse_step_single` delegates to `run_model_inverse`. Python writes every configured member block into `ensemble_output_<member>.h5`, selects `fname='inverse_state.mat'`, and replaces the command's `k` with `km` for the MATLAB call. MATLAB uses `km=k+1` to select the one-based observation column. The wrapper's local forecast-time interval is still derived from the original model index. These two uses of an index must not be conflated when extending the interface.

The inverse branch loads `Models/ens_id_0/<reference_data>` to obtain its model, boundary conditions and solver context. It reads target velocity observations from `synthetic_obs.h5`; the analyzed velocities are an initialization, **not** substitutes for those targets. The helper accepts either dense `/hu_obs` or compact `/hu_obs_compact` with zero-based `/obs_indices` and optional `/obs_active`. It translates active global velocity rows to nodal vectors and supports either compact matrix orientation. Missing entries remain `NaN`; having some active rows of each component is not proof of complete spatial coverage or a suitable objective mask. Validate sparse-velocity use with the installed ISSM consistency checks.

Before solving, the branch reads analyzed `Vx`, `Vy` and absolute `coefficient`. It uses analyzed thickness for pressure and analyzed thickness/bed to construct the ocean level set, zeros the coefficient where that level set is negative, sets the selected friction exponents, zero basal-melt forcing and boundary settings, and calls BAMG. It then constructs the inverse object and solves stress balance.

There are two important limits in the inspected sequence:

1. Assignments of analyzed thickness, surface, bed and base to `md.geometry` **before** the inverse solve are commented out. Geometry comes from the loaded/reference model and subsequent BAMG processing, while pressure and grounding-mask calculations use analyzed fields. The branch copies analyzed thickness/surface/bed into geometry **after** solving. We cannot describe this as a stress-balance inverse consistently posed on all analyzed geometry without changing and validating that sequence.
2. After solving, it copies `StressbalanceSolution.FrictionCoefficient` into the coefficient, but explicitly leaves `md.initialization.vx/vy` assigned to themselves. It does not explicitly copy `StressbalanceSolution.Vx/Vy` into the returned velocity fields. Those fields follow initialization/BAMG handling; their correspondence to the optimized stress-balance velocities must be checked rather than assumed.

BAMG's changed mesh also creates an obligation to verify nodal sizes, ordering, observation interpolation and control bounds across the solve. The wrapper expects the same packed field size on return. The Python source/tests do not establish that remeshing preserves this contract in every reference model.

Python reads the six output datasets with Fortran-order flattening and returns bed/coefficient only in its joint-estimation branch. The non-joint branch's local `bed_int`/`coefficient_int` variables are initialized only when `k==0`, so later isolated calls do not establish a safe fixed-parameter inverse path. Use the inspected joint-estimation application configuration unless that alternative is separately repaired and tested. A Python exception in the external call kills the server and returns `None`; successful completion of the wrapper alone does not inspect ISSM convergence diagnostics.

## 8 Native Mode 3: a different application contract

The Choi native route wraps complete nodal member arrays in `WholeMemberNativeAdapter` and requires `spatial_ranks=1`. It does not expose ISSM's internal partition as Python-owned spatial slabs. A native inverse callback packs one member, invokes the same legacy Python/MATLAB inverse and returns all available configured fields. The wrapper applies every returned named field. Thus its overwrite policy differs from generic Modes 1 and 2, which accept only friction.

The Choi runner calls the global native cycle without `analysis_row_controls`, and the registered adapter supplies no global physical-analysis finalizer. It does not automatically reproduce generic five-block reduction, coefficient freezing, localization, SMB/bed regularization or geometry controls. A YAML key common to the generic configuration is insufficient evidence that this native route invokes its policy. In particular, the coefficient can enter the native ensemble transform before inversion unless an actual native policy is bound. Mode numbers identify execution routes, not a proof of scientific equivalence between all application features.

The native runtime snapshots owned fields immediately before inversion, invokes each locally scheduled member, gathers error status and restores those fields on any reported failure before raising collectively. It restores the **pre-inversion analyzed state**, not the pre-analysis forecast. The transaction covers registered fields; external MATLAB files, subprocesses and unregistered model context need adapter-specific staging or rollback.

The handover correction makes the Choi hook raise `RuntimeError` when the legacy inverse returns no mapping, with the stable member ID and timestep in its diagnostic. Previously it forwarded `None`, which the general whole-member protocol validly interprets as successful in-place mutation. That ambiguity could suppress an application solver failure. The generic protocol still permits deliberate in-place callbacks returning `None`; only this application's missing-result conversion changes. Focused tests verify failure propagation through the actual native wrapper and unchanged successful mapping updates, alongside the framework's owned-field rollback tests. They do not execute MATLAB.

The runner writes `did_inversion` from the requested active cycle and writes history checkpoints after the cycle. Such metadata identifies a scheduled operation, not a measured optimizer convergence result. This application has no complete native checkpoint-load/resume branch; shard output alone is not a restartable coupled inverse workflow.

## 9 SMB and bed regularization are separate inference stages

For context, the optional SMB helper estimates, member by member,

```text
r = fitted_time_derivative(H) + divergence(mean_history(H*u), mean_history(H*v))
[(1+lambda_t) I + lambda_s L^T L] q_hat = r + lambda_t q_previous
```

The temporal derivative is an ordinary least-squares slope over at least three analyzed-history times; flux derivatives use a local weighted quadratic reconstruction. The default retained window is five. With nonzero basal ablation, this continuity residual represents apparent/net balance, rather than surface input alone. The optional basis projection, blend and clipping order affect the final result. If SMB is directly observed as a parameter, the helper follows a different posterior-projection route. Neither route calls ISSM's friction optimizer or supplies an automatic uncertainty distribution for its regularized estimate.

For bed, the helper takes the proposed increment `d=b_analysis-b_previous`, optionally masks it, and solves `(I+lambda_b L^T L)d_smooth=d`, followed by projection, optional increment cap, blending, bounds and surface-separation processing. It regularizes an increment instead of flattening the entire prior bed. Default spatial weights/blends differ between the two helpers and remain modeling choices. A freeze in an earlier EnKF row operation is not a guarantee against every independently enabled later helper.

Survey activity is resolved using `km` and `bed_snap_cols`. The global `grounded_only`/`observed_only` support gate acts at survey columns; it is not a universal mask on every continuous nonsurvey update. Deterministic direct-observation anchors further modify the bed. The [current-developments chapter](current-developments.md) gives these equations, defaults, gate order, private-state persistence and limited test evidence in full. That chapter is the source for extending these helpers; selecting a friction inverse does not activate them.

## 10 Reproduce and extend the workflow

A reviewable inverse experiment records the exact ICESEE and ISSM revisions, reference `.mat` model, parameterization, mesh/order, resolved YAML, command-line overrides, observation schema/clock, coefficient convention and all objective/control settings. Keep per-member initial coefficient, targets, optimized coefficient, component costs, solver termination and elapsed solve time. Compare against a fixed-background and an EnKF-only experiment under a matched observation/initial-ensemble design; a successful MPI exit is a separate property from better inferred parameters.

From a configured checkout, the focused Python checks are:

```sh
python -m pytest -q src/tests/test_full_parallel_hybrid_contract.py src/tests/test_full_parallel_large_data.py
python -m pytest -q src/tests/test_distributed_native_cycle.py src/tests/test_distributed_whole_member_adapter.py
python -m pytest -q src/tests/test_issm_ismip_choi_native.py src/tests/test_issm_mode3_runner.py
python -m pytest -q src/tests/test_bed_update_domain.py
```

The first pair checks representative reduced-state handoff and publication contracts. Native tests cover callback translation, rollback on reported errors and stubbed runner scheduling. These commands are test gates for Python behavior, not invitations to infer MATLAB results. Configuration loading in the application has filesystem side effects, so run adapter tests in their documented temporary-output context rather than importing the driver to discover defaults.

For this revision, the first three commands' six modules were run together against the guarded application source: **80 tests passed, none skipped**. Restoring only the original missing-result `return None` branch in memory made both new failure-propagation tests fail with `DID NOT RAISE RuntimeError`; the actual checkout was not changed for that counterfactual. These checks demonstrate regression sensitivity and the healthy mapping contract. They overlap earlier infrastructure checks and are not 80 independent physical inverse experiments.

Before a new coupled scientific result, execute a small isolated MATLAB/ISSM case with known targets, inspect convergence and all six returned blocks, then force an inverse error and verify the promised canonical/owned-state protection. Test mesh correspondence, sparse velocity masks, floating-node bounds and the before/after geometry handoff explicitly. A restart test must include the external solver state and plugin history. Use the [onboarding templates](application-onboarding.md) to declare these capabilities; a placeholder callback is not evidence that the inverse exists.

## Source map

All anchors below refer to the scientific commit identified above; the subsequent missing-result guard is identified by its focused patch/tests rather than silently attributed to that older commit.

| Subject | Reviewed source |
|---|---|
| Base state, settings and inference switches | [`ISMIP_Choi/params.yaml`](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/applications/issm_model/examples/ISMIP_Choi/params.yaml#L43) |
| Prior/absolute coefficient convention | [`run_model.m`, initialization](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/applications/issm_model/examples/ISMIP_Choi/run_model.m#L54) |
| Active MATLAB inverse and observations | [`run_model.m:930`](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/applications/issm_model/examples/ISMIP_Choi/run_model.m#L930), [`observation helper:1074`](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/applications/issm_model/examples/ISMIP_Choi/run_model.m#L1074) |
| Python HDF5 inverse wrapper | [`_issm_model.py:227`](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/applications/issm_model/examples/ISMIP_Choi/_issm_model.py#L227) |
| Mode 1 activation/reduction and inverse writer | [`icesee_da_partial_parallel.py:543`](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/src/run_model_da/icesee_da_partial_parallel.py#L543), [`_parallel_i_o.py:604`](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/src/parallelization/_parallel_i_o.py#L604) |
| Mode 2 reduced finalizer and atomic publication | [`EnKF_parallel_io.py:251`](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/src/parallelization/EnKF_parallel_io.py#L251), [`publication/inversion:898`](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/src/parallelization/EnKF_parallel_io.py#L898) |
| Native application callback and cycle selection | [`_issm_native.py:147`](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/applications/issm_model/examples/ISMIP_Choi/_issm_native.py#L147), [`mode3_runner.py:529`](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/applications/issm_model/examples/ISMIP_Choi/mode3_runner.py#L529) |
| Native rollback and valid in-place protocol | [`distributed_native_runtime.py:711`](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/src/parallelization/distributed_native_runtime.py#L711), [`distributed_whole_member_adapter.py:199`](https://github.com/ICESEE-project/ICESEE/blob/3748e416d3a8120cc860119bd7224c6f9ad41ff7/src/parallelization/distributed_whole_member_adapter.py#L199) |
