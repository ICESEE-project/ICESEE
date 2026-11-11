# Numerical analysis and scientific interpretation

In this chapter we connect the ensemble update to the numerical quantities used by ICESEE: forecast members, observation equivalents, error factors, local neighborhoods and post-analysis controls. We describe the implemented calculation and the assumptions needed to interpret it. Correct array shapes, completed timesteps and finite values establish useful computational properties, but they do not establish that observations were applied at the intended time or that a recovered physical parameter is well constrained.

## Ensemble data assimilation in context

Data assimilation combines observations with a dynamical model and a description of uncertainty to estimate an evolving system. The forecast carries information from earlier observations forward; the analysis incorporates information available at the next observation time. This relationship is statistical as well as dynamical: a precise observation should not receive the same weight as a poorly constrained one, and a nearly identical set of members cannot represent every plausible model error. [Carrassi et al. (2018)](https://doi.org/10.1002/wcc.535) provide a broad introduction to these estimation questions and their practical consequences.

We represent forecast uncertainty with an ensemble of model realizations. Each member advances through the nonlinear forecast model, and their sample statistics describe the uncertainty used in the analysis. This construction avoids evolving a full state covariance through a separately linearized forecast model, as developed by [Evensen (1994)](https://doi.org/10.1029/94JC00572). Ensemble representation also gives a natural opportunity to distribute independent forecast work, although the model solver, analysis communication and output still impose their own costs. The practical formulations in [Evensen (2003)](https://doi.org/10.1007/s10236-003-0036-9) provide context for both the ensemble algebra and implementation choices; the exact ICESEE factor is specified below.

The distinction between stochastic and deterministic analysis concerns how we construct the analyzed ensemble. In a perturbed-observation stochastic EnKF, each member uses an observation realization drawn according to an error model. [Burgers et al. (1998)](https://doi.org/10.1175/1520-0493%281998%29126%3C1719%3AASITEK%3E2.0.CO%3B2) explain the role of this construction in the ensemble analysis covariance. Deterministic square-root methods instead update ensemble anomalies to obtain specified analysis statistics without introducing those observation perturbations; [Whitaker and Hamill (2002)](https://doi.org/10.1175/1520-0493%282002%29130%3C1913%3AEDAWPO%3E2.0.CO%3B2) and [Tippett et al. (2003)](https://doi.org/10.1175/1520-0493%282003%29131%3C1485%3AESRF%3E2.0.CO%3B2) describe examples and the nonuniqueness of such transforms. The DEnKF of [Sakov and Oke (2008)](https://doi.org/10.1111/j.1600-0870.2007.00299.x) uses a small-correction approximation with a half-gain anomaly update. These are algorithmic distinctions. Distributing a stochastic transform across MPI ranks does not make it a deterministic square-root filter.

Ice-sheet applications add a distinction between evolving state and uncertain physical parameters. [Bonan et al. (2014)](https://doi.org/10.5194/npg-21-569-2014) studied initial-state and parameter estimation with ETKF in a shallow-ice flowline model. [Gillet-Chaulet (2020)](https://doi.org/10.5194/tc-14-811-2020) examined transient surface observations and basal conditions in a marine flowline experiment. [Choi et al. (2025)](https://doi.org/10.5194/tc-19-5423-2025) considered a two-dimensional setting and observing-system experiments. These studies motivate careful treatment of prior uncertainty and observation design. Their results do not establish the accuracy of a different ICESEE configuration.

## Ensemble convention

We arrange the forecast ensemble `X` as an `n by N` matrix with one member per column. Its model equivalents in active observation space form an `m by N` matrix `Y`, where `Y[:,i] = H(X[:,i])`. `H` can be a model observation callback or an index selection. We compute row-centered anomalies as `Yprime = Y - mean(Y, axis=1, keepdims=True)`. A nonlinear operator must be applied to each forecast member before subtracting the observation-space mean; `H(X-mean(X))` is generally a different object.

The textbook finite-ensemble cross covariance is `Xprime @ Yprime.T / (N-1)`, and the predicted observation covariance is `Yprime @ Yprime.T / (N-1)`. These definitions explain the statistical roles of spread and observation error. The active stochastic transform code uses an analysis factor and a truncated pseudoinverse, described below. A textbook gain expression provides context, but does not by itself specify this exact software calculation.

The common generic state layout stacks variable blocks in `vec_inputs` order. Some row-local controls require equal-sized blocks and explicitly reject `nd % len(vec_inputs) != 0`. Other code supports `var_nd`, scalar inputs, or native block layouts. Variable-size integrations therefore need compatible indexing and controls throughout the selected path; the presence of `var_nd` in the loader does not make equal-block control code applicable.

## Observation-error and analysis-factor modes

The key `enkf_observation_error_mode` selects three distinct analysis constructions in the generic Modes 1/2 paths. Its default is derived from the historically misspelled `use_ensemble_pertubations` option: true selects `legacy_prior_anomalies`, false selects `generated_R`. An explicit setting makes this choice independent of that legacy default and part of the reproducible experiment.

| Mode | Factor `Eta` and innovation matrix `D` | Interpretation |
|---|---|---|
| `stochastic_R` | Draw Gaussian values scaled by active standard deviations; center each observation row across members. `D = y[:,None] + Eta - Y`. | An explicit perturbed-observation construction using configured observation uncertainties. The finite ensemble has a sample factor, rather than an exact prescribed covariance. |
| `generated_R` | Generate spatial fields using the configured field generator, variable scales and deterministic observation/member ordering; select active observed rows and row-center. `D = y[:,None] + Eta - Y`. | A spatially structured generated factor. Physical coordinates, units and scale parameters define the intended spatial error model. |
| `legacy_prior_anomalies` | `Eta = Yprime`; `D = y[:,None] - Y`. | A legacy forecast-dependent analysis factor. It is not independent random observation noise and does not become a prescribed `R` just because `sig_obs` exists. |

`active_observation_std` accepts the supported vector/time-array forms of `error_R`, selects active rows, and requires every active standard deviation to be finite and strictly positive. The HDF5 dataset called `R` is used to obtain these standard deviations in the current pipeline. Thus the name `R` alone does not specify a covariance matrix; its producer and consumer determine whether a value is a variance, a standard deviation or a factor.

Generated factors are produced member-major and variable-major by `iter_generated_observation_columns`. The seed starts from `base_seed + 1000003*(observation_slot+1)`. Modes 1 and 2 share this iterator to avoid rank-placement differences and permit column-streamed construction. The native registered runners inspected here advertise `legacy_prior_anomalies` and `stochastic_R`; `generated_R` is rejected by those policies. The generic Mode-2 capability set consequently differs from the capabilities of these native policies.

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

The loader exposes `analysis_svd_energy` with fallback `0.999` for streamed analysis. Some direct/local transform paths use a hard-coded `0.999`. A change to this value needs matching treatment of both implementations and renewed parity checks, because one configuration key does not control every numerical path. In the legacy factor case `F = 2*Yprime`, `legacy_transform_from_gram` implements the corresponding simplified expression.

Reducing covariance storage does not remove finite-ensemble sampling error, long-range spurious correlations, or the need to assess localization. Matrix-free algebra is a computational property, not a statistical argument that localization is unnecessary.

## Localization

Finite ensembles estimate correlations imperfectly, particularly when the true correlation is weak. Distance-dependent covariance filtering addresses this sampling issue, as illustrated by [Hamill et al. (2001)](https://doi.org/10.1175/1520-0493%282001%29129%3C2776%3ADDFOBE%3E2.0.CO%3B2). [Gaspari and Cohn (1999)](https://doi.org/10.1002/qj.49712555417) provide compactly supported correlation functions used for tapering. We therefore distinguish the cost of evaluating covariance information from the statistical reliability of that information: avoiding a dense matrix does not make distant sample correlations accurate.

Local ensemble analysis is also an algorithmic choice. The LETKF of [Hunt et al. (2007)](https://doi.org/10.1016/j.physd.2006.11.008) is a particular deterministic local transform formulation. ICESEE's grouped local stochastic calculation retains its selected factor and innovation policy; neither grouping neighborhoods nor selecting execution Mode 3 identifies that calculation as LETKF.

The code contains covariance/localization utilities and grouped local analysis. `localization_flag`, `local_analysis`, `global_analysis`, `localized_vars`, radii and taper settings have different roles. The selected runner's analysis branch determines which mechanism is active. A flag alone does not establish that the application implements the requested local analysis.

`compute_local_patches_X5` builds neighborhoods using physical node coordinates and a KD tree. It maps observation rows back to block coordinates, finds observations within each target's radius, and groups target nodes with exactly the same sorted observation-neighbor tuple. One transform per unique tuple replaces repeated identical factorizations. For a manual radius, `localization_radius` may be a scalar or a dictionary by variable. Otherwise `estimate_adaptive_radius` estimates a radius from the bounding-box area and observation density with a target count. This density heuristic is not a physically calibrated correlation length.

Nodes without neighbors are omitted from the local patch list. Their retained or global update is determined by the enclosing analysis branch; omission from the patch list does not by itself imply an unchanged state. An absent coordinate provider can also return no patches. Interpretation of the resulting groups depends on consistent observation and state-coordinate units, the logged neighborhood counts, and the inclusion of the intended targets in `localized_vars` and `vec_inputs`.

The distributed patch machinery validates canonical observation IDs, per-kind coordinates, owned target rows, and cross-kind uniqueness. Joint observation collections preserve ordering across kinds rather than applying arbitrary sequential per-kind filters. This prevents MPI ownership from silently changing a joint stochastic update. A grouped patch is a numerical neighborhood, not a mesh partition or a new scientific observation.

## Controls after analysis

`apply_analysis_controls_local` in [`_mpi_analysis_functions.py`](../../src/parallelization/_mpi_analysis_functions.py) is a shared row-local policy for Modes 1 and 2. It first relaxes analysis increments, optionally caps increments by variable, then inflates anomalies about the controlled row mean, restores frozen variables, and applies the bed-update gate. Defaults for relaxation are a no-op (`1.0`); optional increment limits must be finite and positive. These controls are part of the posterior calculation and therefore part of the scientific configuration.

`state_inflation_factor`, `param_inflation_factor`, and `bed_inflation_factor` can override the general inflation. Inflation multiplies deviations from the row mean; it does not multiply the absolute state. A frozen variable is restored from forecast values during analysis, but may still evolve during forecast. Bed gates and global inference can subsequently apply additional policy, so the full finalization path determines the result for a particular feature combination. An application can install a state-modifying callback to restore its geometric or diagnostic relations. The current native PIG binding has no such callback: its analysis finalizer synchronizes ghosts, and the subsequent forecast recomputes surface/flotation and diagnostic velocity. We therefore do not assume that every saved posterior immediately satisfies those relations.

Finite checks reject non-finite analysis values. The distributed collapsed-spread guard rejects an observed ensemble whose maximum anomaly is no larger than `1e-8` times its forecast scale. This prevents dividing innovations by roundoff-sized spread. It does not diagnose all forms of underdispersion. The scientific question is then whether the prior and forecast uncertainty provide a justified spread source. Suppressing the guard would permit an unstable calculation without resolving that question.

## Initialization and process noise

Initialization uncertainty and forecast uncertainty serve different purposes. `sig_model`, `initial_spread_factor`, parameter spreads and adapter initialization determine the prior. `sig_Q`, the random-field method, temporal correlation, and `process_noise_schedule` determine stochastic forecast increments. `sig_obs` and the observation-error mode govern observation factors. Changing one of these uncertainty sources to compensate for another changes the estimation problem, even when a plotted trajectory becomes smoother.

`random_streams.py` defines stable namespaces for initialization, canonical process noise, and a separate native helper stream. Member ID, timestep, variable index and `base_seed` determine streams, while MPI rank and round order must not. Active native PIG execution wraps the native adapter with `NativeProcessNoiseAdapter`, which calls the canonical `add_member_process_noise`; the existence of the separate `mode3_process_noise_seed` helper alone does not identify a run's active policy.

The noise schedule defaults to `observations`: apply noise in the forecast cycle immediately preceding that state's analysis. `every_step` applies at all forecast cycles; `none` disables it. The canonical bridge replays deterministic AR history up to a restored step to reconstruct stochastic state without changing restored physical fields. Its checkpoint policy identity includes time-grid and noise configuration. A changed seed, schedule, amplitudes, coordinates, or time grid is a different experiment and can invalidate restart.

Native PIG initialization intentionally differs from the generic initial perturbation policy in the inspected validations. Shared forecast noise and a common mode number therefore do not establish member-wise equivalence across all three modes. Same-input parity requires the same realized prior, observations, stochastic policy, physical setup, and numerical settings.

## Physical inference and hybrid inversion

Joint EnKF estimation places parameters in the state vector so cross covariances can update them. A hybrid inversion instead passes an analyzed state to a model-specific deterministic inverse solve. These are different estimation procedures. `inversion_flag`, delayed activation, inversion callbacks, and frozen variables must be described explicitly in a scientific experiment.

The optional inference plugin provides SMB and bed corrections. With `inference_plugin_enabled=False`, its global hook is a no-op and the bed gate retains legacy freezing. `bed_update_mode` can be `legacy`, `snapshots`, or `continuous`; `snapshots` requires nonempty snapshot columns. Domain gates can limit bed changes to all, forecast-diagnosed grounded, or observation-supported locations. The `observed_only` gate intentionally uses observation support so a misclassified member does not create unsupported floating-ice updates.

SMB correction uses continuity-related history and declared regularization/blending; bed correction applies declared spatial and physical increment constraints. These are prior/model assumptions, not observational facts. Parameter identifiability, physical bounds, geometric consistency, observation coverage and holdout predictions require separate scientific evaluation. The presence of unit tests for plugin behavior does not demonstrate successful real-world bed or SMB recovery.

## Model time and evidence to inspect

`timestep_contract.py` defines `state-columns-v2`. Column zero is the pre-forecast initial state. Cycle `k` advances to state column `k+1`; observation column `j>0` is due after cycle `j-1`. Initial-time observations are refused because a distinct pre-forecast analysis stage is not implemented. Requested observation times can map to a model grid, so archive both requested `obs_t` and mapped `obs_index`/`t`.

The emitted `[ICESEE_ANALYSIS]` JSON records and the output verifier identify the actual model time, cycle, observation slot and state column used. Older summaries and `run_schedule` helper comments may describe pre-correction endpoint conventions; they are not stronger evidence than the active contract and recorded events. The time-grid marker and checkpoint schedule protect reuse, but their presence does not prove the observations describe the intended physical date.

A numerical and scientific evaluation therefore includes full-member comparisons, analysis-event timing, prior/posterior spread, innovations, physical bounds and geometry checks, together with errors against synthetic truth or independently withheld observations. A field being finite is a minimal integrity check. An operational wall time is a measurement for that topology, workload and concurrent load; it becomes a scaling result only under a controlled scaling design.

The [annotated reading guide](literature.md) places these references alongside the numerical and implementation questions they help answer.
