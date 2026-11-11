# Current developments and their implementation boundaries

We describe the current inference, configuration, diagnostics and recovery mechanisms as implementation capabilities. We do not assign an introduction date to every feature or claim that all were absent from the earlier ICESEE paper. That paper and its public revisions already discuss filtering alternatives, observation-error modes and grouped/local analysis. The scientific source for this chapter remains commit `3748e416d3a8120cc860119bd7224c6f9ad41ff7`; a feature is usable only through an application/runner that actually calls it.

For the EnKF formulation and execution modes, read [numerics](numerics.md) and [parallelization modes](parallelization-modes.md). The [inversion and inference guide](inversion-and-inference.md) explains ISSM's external friction inverse problem and its generic/native handoffs. Here we distinguish that solver from post-analysis regularization and trace the other supporting developments.

## 1 Locate a capability and its evidence

| Mechanism | Implementation | Current scope and evidence |
|---|---|---|
| Apparent mass-balance inference | `src/utils/stable_smb_inference.py` | Optional global post-analysis helper; source-derived algorithm below; no dedicated manufactured-SMB test found in `src/tests` |
| Bed increment regularization | `src/utils/stable_bed_inference.py` | Optional global helper; constant-increment preservation and related bed tests cover limited contracts |
| Survey cadence, support and anchors | `src/utils/inference_plugin.py` | Explicit observation-column schedule, domain gate and deterministic anchor; dense/compact anchor and mask tests |
| Relaxation, caps, inflation, freezing | `src/parallelization/_mpi_analysis_functions.py` | Shared row controls; native bindings must pass/use the callback; parity tests cover representative layouts |
| Physical finalization | `src/parallelization/_parallel_i_o.py` | Shared generic finalizer, including ISSM geometry; equal-sized variable blocks required |
| YAML inheritance and CLI provenance | `config/config_loader.py`, `config/_utility_imports.py`, `src/utils/run_schedule.py` | Recursive merge and override checks; paths remain application-specific |
| Diagnostic phase reports | `src/utils/performance.py` | One structured/text summary, measured/unmeasured distinctions and registered backend metrics; reporter tests |
| Stochastic geometry cache and JIT coordination | `_error_generation.py`, `pyop2_compile_cache.py` | Cached deterministic graph geometry; communicator-consistent compiled-kernel lookup; separate tests |
| Time-grid and durable histories | `timestep_contract.py`, checkpoint/history modules | Explicit clock identity, owned-field coverage and recovery records; application restore capability remains separate |

These are evidence categories, not a declaration that each target application has completed a PDE or scientific validation. Test names below identify what to exercise when extending the relevant path.

## 2 Apparent mass balance from an analyzed history

`apply_smb_physics_correction` is inactive unless `physics_smb_inference` is true and its caller enables the inference plugin. It expects an SMB block and, for the physics route, equal-sized thickness and horizontal velocity blocks with physical coordinates. Coordinates are multiplied by `mesh_coordinate_scale_to_m`. Thickness is in metres, coordinates in metres and depth-averaged velocities in metres per model-time unit; the result is metres per the same time unit.

For one member, the routine retains a rolling history of analyzed `H`, `u` and `v` and computes

```text
r = least_squares_time_slope(H) + divergence(mean_history(H*u), mean_history(H*v))
```

At least three recorded times are needed. The default history length is five; a repeated final time replaces its previous entry instead of duplicating it. The temporal slope is the ordinary least-squares slope using centered physical times. Spatial flux divergence is obtained from a weighted quadratic moving-least-squares reconstruction on neighboring physical coordinates. The local polynomial uses scaled coordinate offsets, a small ridge and a pseudoinverse fallback. The default divergence-neighbor request is 24, with a minimum stencil chosen to support its five polynomial coefficients.

This residual is **apparent net mass balance**. With a positive basal-loss convention, continuity gives `r = surface_mass_balance - basal_loss`. The routine does not add an independently estimated basal term. Nonzero basal melt, incorrect velocity interpretation, or assimilation jumps therefore prevent interpreting the residual as surface mass balance alone. A windowed analyzed history reduces sensitivity to individual jumps but does not turn every analysis increment into a physical tendency.

The regularized field `q_hat` solves, independently for each member,

```text
[(1 + lambda_time) I + lambda_space L.T L] q_hat
    = r + lambda_time q_previous
```

Here `L` is a symmetric dimensionless graph Laplacian: nearest-neighbor Gaussian adjacency is symmetrized and `D-W` is scaled by the median positive degree. It annihilates a constant field. The equation minimizes a residual penalty, spatial roughness penalty and departure from the previously applied field. The defaults are `smb_spatial_regularization=25`, `smb_temporal_regularization=4` and `smb_graph_neighbors=12`. These are algorithm parameters, not statistically calibrated observation covariances.

The helper blends the current analyzed SMB with the inferred field using `smb_blend_factor` (default 0.35), optionally ramped after `smb_inference_start_time`. Low-order spatial bases are `none`, `x_linear`, `affine` and `quadratic`. It retains what was actually applied as the next temporal reference. Before the configured start, `smb_spinup_hold_factor` can retain the first recorded SMB reference; the default zero leaves the analyzed field unchanged during that hold stage.

### Directly observed SMB follows a different route

When SMB is present in `vec_inputs` and declared in `observed_params`, the helper treats it as already conditioned by the EnKF. It skips the continuity-residual fit, applies only the configured spatial basis and physical clipping, and retains that posterior. Declaring observed SMB without an SMB state block raises an error. This branch avoids silently imposing a second residual-based correction on the directly observed parameter.

The operation order matters. In the physics branch, inferred values are clipped before blending and basis projection; the final projected/blended field is not clipped again. It is therefore inaccurate to describe `smb_physical_bounds` as a guaranteed final bound in every branch. In the directly observed branch, projection precedes the final clip. Record the branch, bounds and basis with the experiment.

## 3 Bed regularization acts on increments

`apply_bed_regularized_correction` stabilizes a bed update produced elsewhere by the filter/inverse. It does not infer bed from the continuity residual or solve ISSM's velocity inverse problem. It recognizes bedrock aliases separately from ice `base`: those surfaces differ beneath floating ice.

Let `b_a` be the proposed bed and `b_p` the previously applied bed. The helper forms `d=b_a-b_p`, multiplies any configured `bed_update_mask` weights clipped to `[0,1]`, and solves

```text
(I + lambda_bed L.T L) d_smooth = d
b_corrected = b_previous + beta_bed * d_smooth
```

Regularization acts on **increments**, so a rough prior bed is not progressively flattened merely because another analysis occurs. A basis projection can retain constant, along-flow linear, affine or quadratic increment structure. An optional `bed_max_update_per_cycle` caps the smoothed increment before blending. Defaults are `bed_spatial_regularization=40`, `bed_graph_neighbors=12` and `bed_update_blend_factor=0.15`. The blend may ramp after `bed_inference_start_time`; `bed_spinup_hold_factor` defaults to one before that start.

The first recorded bed analysis supplies an initial reference. The previous reference can instead come from the forecast reference supplied by the generic finalizer. These are runtime references, not an implicit claim that the helper has loaded a physical time-zero bed. Final operations clip configured bed bounds and then enforce `bed <= surface - separation` when a surface block is present (`bed_min_surface_separation=1` metre by default). This is a geometric safeguard; it does not perform a full flotation/stress-balance solve.

Both regularizers cache sparse factorization data using coordinate summaries and parameter values. Those keys are not complete mesh fingerprints. Retain one mesh/ordering per runtime context and rebuild/reset the private cache when changing it; the separate stochastic graph cache described below uses a content hash.

## 4 Survey cadence and physical support

`resolve_bed_update_schedule` returns the current activity from `km`, the **observation-column** index, and `bed_snap_cols`. `legacy` restores forecast bed outside survey columns when `km` is defined. Explicit `snapshots` requires a nonempty survey list and updates only those columns. `continuous` retains cross-covariance updates at every analysis. Disabling `inference_plugin_enabled` forces legacy scheduling rather than removing the older bed-freezing behavior.

`bed_update_domain` is `all`, `grounded_only` or `observed_only`. In the inspected global gate, restricted domains are applied at bed survey columns. They are not a universal mask on every nonsurvey continuous update. Grounded support is diagnosed per forecast member using `H + bed/di > 0`; observed support uses the stored synthetic-observation mask, with shape and availability checks.

At active finite bed observations, the optional deterministic anchor is

```text
b_member <- b_member + eta * (bed_observation - b_member),  0 <= eta <= 1
```

Pre- and post-stage factors are separate (`bed_observation_nudge_factor` and `bed_observation_post_anchor_factor`, both default zero). This contracts each directly observed innovation by `1-eta` at that operation, while leaving unobserved nodes untouched. It is additional deterministic postprocessing, not an unchanged stochastic EnKF posterior. Dense and compact observation paths are checked for parity, including inactive compact rows.

## 5 Preserve the full analysis-control order

The shared generic row policy first relaxes each variable's proposed increment, optionally caps that increment, inflates anomalies around the controlled row mean, restores `frozen_analysis_vars`, and applies the local bed schedule gate. Relaxation defaults to one; state/parameter/bed inflation controls have separate fallback values. A cap limits an **increment**, not the absolute magnitude of a physically valid velocity or elevation.

Global finalization then applies bed support/anchors, stores a forecast reference, runs bed regularization before geometry, reapplies support and post anchors, restores the application's physical geometry where implemented, and runs SMB correction after geometry. For ISSM the geometry mode preserves either thickness or surface, enforces minimum thickness and diagnoses floating/grounded consistency. Subsequent inverse/model callbacks can change fields further. Freezing in the row policy is not a blanket promise that an independently enabled later inference/geometry callback will never touch that variable.

These shared routines require equal-sized blocks. A native segmented adapter must explicitly bind the relevant row controls or model-specific finalizer. The generic native application rejects unsupported inference flags unless its policy opts in; registration alone does not wire the full global plugin. The specialized PIG native analysis/finalizer does not automatically run these gathered generic inference routines. Verify the active callback instead of transferring a YAML feature claim between modes.

Mode 2 preserves private `_bed_`/`_smb_` state between method calls because `**kwargs` makes a new mapping. Ensemble-coupled inference can require a complete ensemble under a separate finalization memory budget; ordinary member-local geometry can be processed in member slabs. Bounded row analysis is therefore not proof that this optional global stage has the same bound. In-memory temporal/plugin state also needs an explicit restart policy; object-lifetime copying is not itself durable checkpoint restoration.

## 6 Configuration and provenance

The YAML loader resolves `extends` relative to the child file and recursively merges mappings; child non-mapping values replace the corresponding parent value. Cyclic inheritance is rejected. Retain the entire referenced file chain as well as a resolved configuration, because an unchanged child file can otherwise select different parent values later. Verify referenced filenames exist before launching an experiment.

The CLI/YAML import layer flattens physical, modeling and EnKF sections, then reapplies normalized CLI-selected values. Generic `--key=value`/`--key value` overrides are coerced using YAML parsing. The dispatcher checks recorded CLI overrides against the final application context and rejects values an application silently replaced. Record both requested and effective configuration; source/application setup can transform units and time grids. Do not import `_utility_imports.py` merely to inspect defaults: fresh configuration can remove an existing output directory.

Inference timing must also name its clock. `resolve_analysis_cycle_time` returns requested observation time and mapped state-grid time separately for delayed inversion decisions. The generic finalizer currently passes `timestep*dt` to the regularizers. Applications with nonzero initial time or a different convention must verify that mapping rather than assuming every helper uses one absolute-time origin.

## 7 Diagnostics, caches and recoverable output

The generic performance reporter preserves historical metrics while emitting `performance.json` and text from one summary. Wall time is the maximum elapsed rank time; aggregate computational time is a rank sum. Phase timers may overlap and cannot automatically be added into a percentage breakdown. An uninstrumented phase is “not measured”; a registered zero operation count is “no events”; measured zero is a third case. Backend resource plans and I/O counters are registered by their owners rather than inferred from a model name.

The stochastic graph transition cache in `_error_generation.py` stores deterministic geometry work under a 64-MiB per-process limit. Its key includes a content hash of coordinates and graph parameters. Tests compare cached/uncached fields and random state and reject stale reuse after coordinate mutation. It caches geometry, not previously sampled random fields.

The PyOP2 cache wrapper coordinates hit/miss across the compilation communicator before entering a collective compile. This prevents mixed cache hits from leaving only some ranks in the compile broadcast. It retains the compiled-kernel key/write protocol and records JIT timing. It is version-structure sensitive: inspect whether installation succeeded; do not assume every later PyOP2 release has the expected wrappers. Node-local JIT caches and rank-private XDG state remain deployment choices described in [deployment](deployment.md).

Clock markers, ownership fingerprints, atomic checkpoint publication and retained backing histories support recovery. Their existence does not supply every application's restore callback. Preserve old-clock outputs, keep VDS backing files, and test a restart across analysis with the chosen inference/solver state. The [operations guide](operations.md) separates simulation completion, allocation status and output verification.

## 8 Validation and remaining scientific work

Current tests cover bed schedule/support/anchor contracts (`test_bed_update_domain.py`), survey-track masks, bed-background invariance, local/native analysis controls, hybrid finalization and representative YAML inheritance. `test_full_parallel_large_data.py` checks compact metadata and ensemble-coupled finalization behavior. Graph-cache, PyOP2 communicator and performance-report tests cover their corresponding infrastructure. A manufactured-solution SMB suite was not found in the inspected test tree; these source-derived equations should not be presented as an independently validated physical SMB product.

Before scientific use, test a known continuity residual with correct units and basal balance, constant/rough bed increments, projection/bound ordering, survey/non-survey cadence and restart of private temporal references. Then evaluate the target observation experiment using withheld data or a known truth. These checks assess the scientific regularization assumptions; an execution-parity test addresses a separate question.
