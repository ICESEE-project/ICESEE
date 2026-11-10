# Numerical analysis and scientific interpretation

This chapter documents the numerical objects that production operators must recognize in logs, configuration, and code. It separates an implemented algorithm from evidence that a particular physical inference is accurate. A run can conserve shapes, finish all timesteps, and remain finite while using the wrong observations, insufficient ensemble spread, or an unacceptable physical parameter estimate.

## Ensemble convention

Let `X` be an `n by N` forecast ensemble with one member per column. Let `Y` be its model equivalents in active observation space, of shape `m by N`; `Y[:,i] = H(X[:,i])`. `H` can be a model observation callback or an index selection. Form the row-centered anomalies `Yprime = Y - mean(Y, axis=1)`. A nonlinear operator must be applied to each forecast member before subtracting the observation-space mean; `H(X-mean(X))` is generally a different object.

The textbook finite-ensemble cross covariance is `Xprime @ Yprime.T / (N-1)`, and the predicted observation covariance is `Yprime @ Yprime.T / (N-1)`. These definitions explain the statistical roles of spread and observation error. The active stochastic transform code uses an analysis factor and a truncated pseudoinverse, described below; do not replace its definition with a textbook gain expression when reporting exact software behavior.

The common generic state layout stacks variable blocks in `vec_inputs` order. Some row-local controls require equal-sized blocks and explicitly reject `nd % len(vec_inputs) != 0`. Other code supports `var_nd`, scalar inputs, or native block layouts. A variable-size integration must verify every operation it uses; the existence of `var_nd` in the loader does not make equal-block control code applicable.

## Observation-error and analysis-factor modes

The key `enkf_observation_error_mode` selects three distinct analysis constructions in the generic Modes 1/2 paths. Its default is derived from the historically misspelled `use_ensemble_pertubations` option: true selects `legacy_prior_anomalies`, false selects `generated_R`. Set the mode explicitly for a reproducible experiment.

| Mode | Factor `Eta` and innovation matrix `D` | Interpretation |
|---|---|---|
| `stochastic_R` | Draw Gaussian values scaled by active standard deviations; center each observation row across members. `D = y[:,None] + Eta - Y`. | An explicit perturbed-observation construction using configured observation uncertainties. The finite ensemble has a sample factor, rather than an exact prescribed covariance. |
| `generated_R` | Generate spatial fields using the configured field generator, variable scales and deterministic observation/member ordering; select active observed rows and row-center. `D = y[:,None] + Eta - Y`. | A spatially structured generated factor. Match physical coordinates, units, and scale parameters to the intended error model. |
| `legacy_prior_anomalies` | `Eta = Yprime`; `D = y[:,None] - Y`. | A legacy forecast-dependent analysis factor. It is not independent random observation noise and does not become a prescribed `R` just because `sig_obs` exists. |

`active_observation_std` accepts the supported vector/time-array forms of `error_R`, selects active rows, and requires every active standard deviation to be finite and strictly positive. The HDF5 dataset called `R` is used to obtain these standard deviations in the current pipeline. Do not assume every object named `R` is a covariance matrix: inspect its producer and consumer before taking a square root or squaring it.

Generated factors are produced member-major and variable-major by `iter_generated_observation_columns`. The seed starts from `base_seed + 1000003*(observation_slot+1)`. Modes 1 and 2 share this iterator to avoid rank-placement differences and permit column-streamed construction. The native registered runners inspected here advertise `legacy_prior_anomalies` and `stochastic_R`; `generated_R` is rejected by those policies. A generic Mode-2 feature table must not be projected onto every Mode-3 runner.

## The implemented ensemble transform

In [`compute_X5_from_matrices`](../../src/utils/localization.py), define `F = Yprime + Eta`. If `F = U S V.T`, the code retains singular directions accounting for approximately 99.9% of the summed squared singular values. Let `Q` be the truncated pseudoinverse of `F @ F.T`. The transform is

```text
X5 = I + Yprime.T @ Q @ D
Xa = Xf @ X5
```

This is the actual factor-form update. In particular, using `F @ F.T` introduces the cross terms between forecast anomalies and the finite observation factor. These cross terms do not vanish identically in a finite realization. The code's factor formulation should not be described as algebraically identical to every stochastic gain that uses `Yprime @ Yprime.T + (N-1)*R`.

The streamed implementations avoid the `m by m` object by accumulating

```text
C = Yprime.T @ F         # N by N
G = F.T @ F              # N by N
B = F.T @ D              # N by N
X5 = I + C @ Gplus @ Gplus @ B
```

Here `Gplus` is a spectral pseudoinverse. The two pseudoinverse factors follow the identity `(F F.T)+ = F (F.T F)+ (F.T F)+ F.T`; they are intentional. `ensemble_transform_from_products` symmetrizes `G`, orders eigenvalues, clips negative roundoff to zero, retains the configured energy fraction, and discards eigenvalues below a machine-epsilon threshold. Distributed products reduce only ensemble-sized matrices across observation owners. Row chunks then apply `X5` to state rows, avoiding construction of a global covariance or dense `H`.

The loader exposes `analysis_svd_energy` with fallback `0.999` for streamed analysis. Some direct/local transform paths use a hard-coded `0.999`. If changing this value, check both implementations and rerun parity; do not assume one knob controls all numerical paths. In the legacy factor case `F = 2*Yprime`, `legacy_transform_from_gram` implements the corresponding simplified expression.

Reducing covariance storage does not remove finite-ensemble sampling error, long-range spurious correlations, or the need to assess localization. Matrix-free algebra is a computational property, not a statistical argument that localization is unnecessary.

## Localization

The code contains covariance/localization utilities and grouped local analysis. `localization_flag`, `local_analysis`, `global_analysis`, `localized_vars`, radii and taper settings have different roles. Trace the selected runner's analysis branch to see which mechanism is active. Enabling a flag that an application does not implement is not a validated local analysis.

`compute_local_patches_X5` builds neighborhoods using physical node coordinates and a KD tree. It maps observation rows back to block coordinates, finds observations within each target's radius, and groups target nodes with exactly the same sorted observation-neighbor tuple. One transform per unique tuple replaces repeated identical factorizations. For a manual radius, `localization_radius` may be a scalar or a dictionary by variable. Otherwise `estimate_adaptive_radius` estimates a radius from the bounding-box area and observation density with a target count. This density heuristic is not a physically calibrated correlation length.

Nodes without neighbors are omitted from the local patch list. The enclosing analysis branch determines their retained/global update; inspect it before interpreting them as unchanged. An absent coordinate provider can also return no patches. Verify that the expected local groups appear in logs, that observation units match state-coordinate units, and that an intended target variable is actually in `localized_vars` and `vec_inputs`.

The distributed patch machinery validates canonical observation IDs, per-kind coordinates, owned target rows, and cross-kind uniqueness. Joint observation collections preserve ordering across kinds rather than applying arbitrary sequential per-kind filters. This prevents MPI ownership from silently changing a joint stochastic update. A grouped patch is a numerical neighborhood, not a mesh partition or a new scientific observation.

## Controls after analysis

`apply_analysis_controls_local` in [`_mpi_analysis_functions.py`](../../src/parallelization/_mpi_analysis_functions.py) is a shared row-local policy for Modes 1 and 2. It first relaxes analysis increments, optionally caps increments by variable, then inflates anomalies about the controlled row mean, restores frozen variables, and applies the bed-update gate. Defaults for relaxation are a no-op (`1.0`); optional increment limits must be finite and positive. These controls change the posterior calculation and must be recorded as scientific configuration.

`state_inflation_factor`, `param_inflation_factor`, and `bed_inflation_factor` can override the general inflation. Inflation multiplies deviations from the row mean; it does not multiply the absolute state. A frozen variable is restored from forecast values during analysis, but may still evolve during forecast. Bed gates and global inference can subsequently apply additional policy; inspect the full finalization path for the requested feature combination. Model callbacks then restore geometric or diagnostic consistency, for example updating surface or diagnostic velocity in a model-consistent order.

Finite checks reject non-finite analysis values. The distributed collapsed-spread guard rejects an observed ensemble whose maximum anomaly is no larger than `1e-8` times its forecast scale. This prevents dividing innovations by roundoff-sized spread. It does not diagnose all forms of underdispersion. Supply a justified spread source and examine innovations and ensemble spread; do not bypass the guard to complete a job.

## Initialization and process noise

Initialization uncertainty and forecast uncertainty serve different purposes. `sig_model`, `initial_spread_factor`, parameter spreads and adapter initialization determine the prior. `sig_Q`, the random-field method, temporal correlation, and `process_noise_schedule` determine stochastic forecast increments. `sig_obs` and the observation-error mode govern observation factors. They must not be interchanged to obtain a visually smoother trajectory.

`random_streams.py` defines stable namespaces for initialization, canonical process noise, and a separate native helper stream. Member ID, timestep, variable index and `base_seed` determine streams, while MPI rank and round order must not. Active native PIG execution wraps the native adapter with `NativeProcessNoiseAdapter`, which calls the canonical `add_member_process_noise`; the existence of the separate `mode3_process_noise_seed` helper alone does not identify a run's active policy.

The noise schedule defaults to `observations`: apply noise in the forecast cycle immediately preceding that state's analysis. `every_step` applies at all forecast cycles; `none` disables it. The canonical bridge replays deterministic AR history up to a restored step to reconstruct stochastic state without changing restored physical fields. Its checkpoint policy identity includes time-grid and noise configuration. A changed seed, schedule, amplitudes, coordinates, or time grid is a different experiment and can invalidate restart.

Native PIG initialization intentionally differs from the generic initial perturbation policy in the inspected validations. Shared forecast noise and a common mode number therefore do not establish member-wise equivalence across all three modes. Same-input parity requires the same realized prior, observations, stochastic policy, physical setup, and numerical settings.

## Physical inference and hybrid inversion

Joint EnKF estimation places parameters in the state vector so cross covariances can update them. A hybrid inversion instead passes an analyzed state to a model-specific deterministic inverse solve. These are different estimation procedures. `inversion_flag`, delayed activation, inversion callbacks, and frozen variables must be described explicitly in a scientific experiment.

The optional inference plugin provides SMB and bed corrections. With `inference_plugin_enabled=False`, its global hook is a no-op and the bed gate retains legacy freezing. `bed_update_mode` can be `legacy`, `snapshots`, or `continuous`; `snapshots` requires nonempty snapshot columns. Domain gates can limit bed changes to all, forecast-diagnosed grounded, or observation-supported locations. The `observed_only` gate intentionally uses observation support so a misclassified member does not create unsupported floating-ice updates.

SMB correction uses continuity-related history and declared regularization/blending; bed correction applies declared spatial and physical increment constraints. These are prior/model assumptions, not observational facts. Parameter identifiability, physical bounds, geometric consistency, observation coverage and holdout predictions require separate scientific evaluation. The presence of unit tests for plugin behavior does not demonstrate successful real-world bed or SMB recovery.

## Model time and evidence to inspect

`timestep_contract.py` defines `state-columns-v2`. Column zero is the pre-forecast initial state. Cycle `k` advances to state column `k+1`; observation column `j>0` is due after cycle `j-1`. Initial-time observations are refused because a distinct pre-forecast analysis stage is not implemented. Requested observation times can map to a model grid, so archive both requested `obs_t` and mapped `obs_index`/`t`.

Use the emitted `[ICESEE_ANALYSIS]` JSON records and the output verifier to check the actual model time, cycle, observation slot and state column used. Older summaries and `run_schedule` helper comments may describe pre-correction endpoint conventions; they are not stronger evidence than the active contract and recorded events. The time-grid marker and checkpoint schedule protect reuse, but their presence does not prove the observations describe the intended physical date.

For publication and production sign-off retain full-member comparisons, analysis-event timing, prior/posterior spread, innovations, physical bounds and geometry checks, along with errors against synthetic truth or independently withheld observations. A field being finite is a minimal integrity check. An operational wall time is a measurement for that topology, workload and concurrent load; it becomes a scaling result only under a controlled scaling design.
