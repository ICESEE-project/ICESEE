# ICESEE–ISSM rebuttal experiments: implementation record

This document records the changes made while constructing the reviewer
experiments and translates them into language suitable for the methods,
experimental-design, and limitations sections of the paper. The canonical
configuration is `rebutal_experiments/param.yaml`; WBF, EBF, and IBF inherit
from it and differ only in their basal-friction treatment.

## 1. Reviewer experiment design

Three synchronized experiments were separated explicitly:

- **WBF — wrong basal friction, fixed:** the deliberately incorrect basal
  friction is held fixed. The EnKF updates geometry and velocity but restores
  the coefficient after each analysis.
- **EBF — EnKF basal-friction estimation:** friction is included in the EnKF
  state/parameter analysis and is localized with the same 15 km radius used
  for bed updates.
- **IBF — inversion-based basal-friction estimation:** the EnKF analysis does
  not directly update friction. Instead, member-wise ISSM inversion updates
  the coefficient after the analysis. `frozen_analysis_vars: [coefficient]`
  therefore prevents a second, conflicting EnKF coefficient increment; it
  does not disable the ISSM inversion.

All three cases use the same truth, prior ensemble construction, observations,
assimilation schedule, bed surveys, localization, inflation, and geometry
projection. This makes the friction method the controlled difference.

The present canonical setup runs for 100 years, assimilates annual surface and
velocity observations from year 2 through year 55, and uses sparse bed surveys
at years 2, 8, 14, 20, 24, 30, 36, 40, 46, and 50.

## 2. Configuration and reproducibility

- YAML inheritance was added with recursive deep merging.
- The experiment tree was flattened to one self-contained `param.yaml`, with
  one shallow child for each friction method.
- Each experiment writes to a distinct dataset directory to prevent accidental
  resume from, or overwrite of, another experiment.
- Fresh-start, restart, checkpoint, and partial-run paths were made explicit.
- Initial-state-only generation was added so a prior can be inspected before
  paying for a full ensemble simulation.
- Synthetic observation times and bed-snapshot columns are printed at startup
  to make the data schedule auditable.

## 3. Observation design

- Surface elevation and horizontal velocity are the regularly observed state
  variables.
- Bed elevation is observed only during sparse survey campaigns.
- Bed tracks now use independent strides in both x and y rather than a small
  number of quasi-one-dimensional tracks.
- The observation-support mask is stored and checked against the global mesh.
  Shape mismatches now fail clearly instead of silently applying an incorrect
  mask.
- Bed updates are restricted to grounded ice and localized within 15 km.

This design represents frequent remotely sensed surface/velocity information
and comparatively sparse radar-derived bed information.

## 4. Initial-condition construction

The initial prior is deliberately imperfect but internally physical.

### 4.1 Geometry consistency

After every configured perturbation, geometry is reconstructed so that

`surface = base + thickness`.

Floating base is recomputed from hydrostatic equilibrium and the ocean
level-set mask is recalculated from the flotation function. This eliminated
ISSM failures caused by independently updated surface, base, and thickness.

### 4.2 Heterogeneous thickness and surface

The prior thickness combines a mean multiplicative factor with a smooth,
mixed-sign, deterministic anomaly. The current canonical values are:

- mean thickness scale: 1.05;
- fractional-anomaly standard deviation: 0.10;
- additive-anomaly scale: 120 m;
- permitted additive range: −180 to +180 m;
- floating anomaly factor: 0.25;
- multiplicative factor range: 0.90–1.10.

Surface is not perturbed independently. It follows from the consistent
bed/base and thickness construction.

### 4.3 Bed prior and floating-bed leakage fix

Early implementations accidentally replaced much of the floating bed with the
hidden true bed, producing an artificially perfect shelf and a truth-shaped
edge at the grounding line. That behavior was removed:

- optional masking is diagnosed from the smooth background, not truth;
- masked floating vertices retain the background prior, not the true bed;
- the audit reports the fraction exactly equal to truth.

The accepted bed prior uses the survey-supported upstream kriging-error
structure as a template and stretches it continuously across the full domain.
Its amplitude changes smoothly downstream rather than stitching two fields or
adding Gaussian blobs. Current controls are:

- mean bed offset: −80 m;
- anomaly scale: 120 m;
- smooth bounds: −250 to +150 m;
- source/template limit: x ≤ 300 km;
- downstream anomaly factor: 0.98 in the current long-run profile;
- floating anomaly factor: 0.95;
- floating maximum absolute error: 100 m;
- 40 km grounding-zone transition.

The v6 initial-state audit (which used a stronger downstream taper of 0.60)
gave a grounded-bed RMSE of about 80.5 m, a floating-bed RMSE of about 25.9 m,
and zero vertices copied exactly from truth. The current 0.98 profile must be
identified separately when its statistics are quoted.

## 5. EnKF stabilization and geometry projection

- Variable-specific relaxation and increment limits were introduced.
- The global analysis relaxation is 0.7.
- Bed relaxation is 1.0, but its analysis increment is capped at 20 m per
  cycle.
- Thickness and surface increments are capped at 300 m; velocity components
  at 2500 m/yr.
- Bed updates can be spatially regularized using a 12-neighbor graph.
- Geometry projection uses `preserve_surface`, so accepted bed changes are
  reconciled with thickness without violating `S = B + H`.
- `frozen_analysis_vars` restores specified analysis blocks from the forecast
  ensemble. This is used to distinguish WBF/IBF from EBF cleanly.

## 6. Bed inference

- Bed observations are anchored only at their supported survey locations.
- A configurable observation nudge and post-anchor factor were added.
- Bed-domain gating is applied consistently in global and scattered parallel
  I/O paths.
- Per-cycle bed increments are bounded and the bed is kept below the surface.
- Diagnostics report the number of unique localized observation groups and the
  effective localization radius.

## 7. Friction inversion

- Inversion activation was separated from inversion availability, allowing a
  geometry-only spin-up before inversion in earlier tuned experiments.
- Inversion bounds were widened during tuning to prevent premature saturation.
- Absolute and relative velocity misfit weights and Tikhonov regularization
  are configurable.
- WBF, EBF, and IBF now use explicit method flags rather than relying on
  ambiguous combinations of inversion and joint-estimation settings.

## 8. ISSM coupling and failure handling

- Missing per-member model caches can be reconstructed to permit safe resume.
  This rebuild initializes model structure; it does not substitute truth for
  the evolving ensemble state.
- Forecasts are checked for non-finite or catastrophic velocity before and
  after ISSM.
- A rejected catastrophic member uses one-cycle persistence rather than
  contaminating ensemble covariances with an unphysical forecast.
- The rejection is logged with member ID and maximum speed.
- MATLAB server status handling, command startup, and partial-run resume
  diagnostics were strengthened.

Persistence is a numerical safety mechanism, not a successful model forecast;
the frequency of persistence events should be reported or at least inspected.

## 9. HDF5 output and result reading

- Output-path handling was made consistent with the selected experiment
  directory.
- Resume now checks that the requested ensemble file exists and gives an
  explicit error otherwise.
- `read_results.m` was rewritten to resolve the active dataset, reconstruct
  the correct saved variables, use dynamic truth domains, and distinguish the
  assimilation and forecast periods.
- A guarded `read_results_v1.m` variant was retained for experimental
  difference-map diagnostics.
- Difference maps replaced unstable relative errors near zero-valued fields.

## 10. Diagnostics and figures

- RMSE panels were restored to the solid/no-DA and dotted/assimilated style
  used in the ICESEE EnKF paper.
- Red, blue, and cyan consistently denote grounded, grounded excluding the
  grounding-line zone (or floating for speed), and whole-ice domains.
- Bed and friction no-DA curves are evaluated against fixed parameters rather
  than being interpreted as evolving parameter updates.
- Grounding-line displacement is computed from the actual centerline
  zero-crossing of the flotation function.
- Map panels distinguish true fields from prior-minus-truth differences and
  avoid setting floating truth fields artificially to zero.
- `plot_filter_comparison.m` provides the requested seven-panel WBF/EBF/IBF
  comparison: four friction maps and three time-series metrics.

## 11. New opposite-sign grounding-line robustness experiment

The long result supplied for review starts with the prior grounding line
landward of truth and requires nearly 100 years for convergence. A separate
experiment now tests the opposite initialization.

The new profile is in `rebutal_experiments/seaward_gl_prior/`. It inherits all
accepted settings and changes only two controls:

- `initial_gl_seaward_thickness_m: 600`;
- `initial_gl_seaward_width_m: 50000`.

A compact smoothstep thickness perturbation is applied across both sides of
the initial grounding zone. It first repairs local retreat introduced by the
heterogeneous background prior, then grounds a controlled shelf strip. Its
distance is measured along flow from the locally upstream GL front and is zero
beyond 50 km; the lateral arms of the U-shaped grounding line therefore cannot
spread the perturbation across the shelf. The grounding-line change follows
from the flotation function and internally consistent geometry; no contour is
shifted for plotting. The final initial-state audit must be approved before a
full run.

## 12. Suggested paper language

> We conducted three synchronized data-assimilation experiments that differed
> only in their treatment of basal friction: a fixed wrong-friction control
> (WBF), direct EnKF coefficient estimation (EBF), and member-wise variational
> inversion coupled to the EnKF geometry analysis (IBF). Surface elevation and
> velocity were assimilated annually, whereas bed elevation was assimilated
> only during sparse survey campaigns and only over grounded ice. Bed
> increments were localized within 15 km, spatially regularized, and limited
> to 20 m per analysis cycle.

> Initial geometry errors were spatially heterogeneous and mixed in sign.
> Surface elevation was not perturbed independently; after perturbing bed and
> thickness we reconstructed base, surface, flotation, and the grounded/floating
> mask to satisfy geometric and hydrostatic constraints. Floating-bed priors
> were generated independently of the hidden truth. We additionally performed
> opposite-sign grounding-line initialization tests to assess sensitivity to
> whether the prior grounding line began landward or seaward of the reference.

## 13. Items to report transparently

- State the ensemble size actually used for the final figures.
- Report the assimilation window and total forecast length associated with
  each figure; several earlier tuning runs used shorter windows.
- Report or inspect the number of catastrophic-member persistence fallbacks.
- Describe the seaward-GL experiment as a controlled robustness test, not an
  observationally derived initialization.
- Keep failed/tuning runs out of the quantitative comparison, but retain their
  configuration names for provenance.
