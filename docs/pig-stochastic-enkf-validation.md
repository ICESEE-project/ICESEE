# PIG stochastic EnKF validation

These controls retain the stochastic EnKF and its perturbed observations. They
do not enable LETKF or a serial square-root filter. Historical configurations
remain selectable. No new policy should be promoted on completion/finite values alone.

## Prescribed observation covariance

```yaml
enkf-parameters:
  filter_type: EnKF
  enkf_observation_error_mode: stochastic_R
  enkf_gain_formulation: prescribed_r
  local_analysis: false
  localization_flag: false
```

For forecast anomalies `Y`, diagonal observation covariance `R` and member
innovations `D = d + eta - H(X)`, the streamed products are
`G = Y.T R^-1 Y` and `B = Y.T R^-1 D`. The stochastic transform is
`T = I + [(N-1)I + G]^-1 B`. The implementation removes roundoff in the
ensemble-mean direction. Model rows apply `X @ T`; no observation covariance
matrix of size `m x m` is constructed. Observation perturbations retain their
member/time-keyed streams and ensemble centering.

`sampled_factor` is the historical default. `prescribed_r` requires positive,
finite diagonal variances, at least two members and `stochastic_R`. It is wired
to the global Mode-1 in-memory, Mode-2 file-backed and Mode-3 native resident/store
paths. Mode-1 partitioned I/O and grouped/local analysis reject this policy.
Correlated observation errors require a separate whitening implementation.
`analysis_svd_energy` applies to the historical sampled-factor gain, not this
positive-definite solve. PIG's full observation arrays are still replicated.

Tests compare complete analyzed members against an independent dense Kalman gain,
include many-more-observations-than-members cases and check posterior covariance
in expectation. A single stochastic realization need not have the exact Kalman
posterior covariance. These are algebraic tests, not real-model validation.

## Independent physical controls

`pig_analysis_state_policy: diagnostic_balance` freezes direct `s/u/v` analysis
increments, floors thickness using the existing `hThresh`, derives geometry and
solves diagnostic velocity using the existing solver settings. The same consistency
step runs after geometry-changing native forecast noise; melt-only noise avoids
an unnecessary repeated velocity solve. Returned/checkpointed values are the
finalized fields. Initial physical fields remain unperturbed.

Set `pig_report_analysis: true` to print per-member thickness minima, maximum
floor correction and volume added by projection. Set an explicitly chosen
`pig_max_thickness_projection` (thickness units) to stop before a correction
exceeds that limit. Model-group failures propagate before the next global phase.
Projection is not a mass-conserving constrained optimization; evaluate its volume
effect and do not accept a trajectory merely because projection removed negatives.

## Estimated depth-law forcing

`basal_melt_forecast: estimated_field` retains the previous persistent-field
experiment. The new `estimated_depth_law` is a Mode-3 candidate requiring
`diagnostic_balance`, explicit prior/bounds and `basal_melt_amplitude_trend`.

The predictor projects the analyzed melt field onto the unit-amplitude depth
profile using finite-element integrals: `a = integral(m*b)/integral(b*b)`.
It evaluates the known control/warm chronology on current draft and floating
geometry. The next amplitude is bounded `a + trend*dt`; the saved melt field is
evaluated on the new geometry at the new state time. Checkpoints therefore retain
the actual melt field from which amplitude can be recovered without hidden state.
Zero profile support refuses the experiment because amplitude is unidentifiable.

The drift has units `(m water/year)/year`; it is an explicit sensitivity assumption,
not a learned trend or the truth ramp. No truth values enter the predictor.
Online joint estimation of drift is not implemented. Test an independently
justified drift, including zero, against held-out data before an 82-year prediction.
The scalar amplitude cannot represent arbitrary spatial forcing errors.

For uncertainty in this scalar control, set
`process_noise_spatial_policy: {basal_melt_field: global_scalar}` and
`basal_melt_noise_basis: depth_amplitude`, with `physical_ar1`. Each member/event
gets an independent scalar AR innovation; the adapter multiplies it by the
current unit-amplitude depth profile. Thus `sig_Q` for melt is an amplitude
increment standard deviation, and the covariance is `sigma_a^2 b(x)b(x')`.
This prevents spatial averaging from accidentally suppressing amplitude noise.
Use no independent thickness noise when forcing is the experiment's only model
uncertainty. A Gaussian field length is not required for a scalar amplitude;
geometry supplies its spatial pattern.

## Physical noise

`process_noise_normalization: physical_ar1` requires an explicit nonnegative
`process_noise_correlation_time` in model-time units and
`process_noise_start: stationary` or `zero`. For consecutive noise events separated
by `delta_t`, the AR coefficient is `exp(-delta_t/tau)` (`tau=0` is white).
`stationary` starts with a full random draw at the first forecast-noise event;
it does not perturb the initial state. `zero` starts from zero at the noise origin.
`sig_Q` specifies marginal increment standard deviations under stationary forcing,
not diffusion coefficients. Temporally correlated increments can accumulate;
their amplitudes must be calibrated for the chosen update schedule.

`process_noise_vars` selects independent state blocks; derived fields should be
recomputed rather than independently perturbed. Noise replay preserves member,
variable, coordinate and event identities. Changed controls invalidate restart
identities. `per_step` and `legacy_run_length` retain their previous meanings.

With graph noise, `graph_physical_length: true` makes `length_scale` enter a
Gaussian spatial kernel in **mesh-coordinate units**. Its continuum e-folding
covariance length is approximate on a nonuniform bounded mesh. Measure the actual
correlation. A radius exceeding `graph_max_neighbors` or the 64 MiB sparse-kernel
budget fails rather than silently truncating neighbors. The default graph policy
continues using adjacency/passes and does not interpret `length_scale` this way.

## Reuse references and stop early

`reference_data_path` names a read-only directory containing the completed
`true_nurged_states.h5` and `synthetic_obs.h5`. `reference_params_file` must name
the frozen YAML that generated them. The runner checks the declared forcing
horizon, inputs, physical thresholds, observation configuration and data shapes.
This checks provenance declarations, not a checksum of every large input file.
The output and reference directories must be separate and non-nested; this is
checked before output cleanup. Truth profiles remain in the reference directory.

`run_stop_time` must lie on the configured `dt` grid. It limits execution without
changing `num_years`, `nt`, observations or the truth's ramp. Use checkpoints every
step. A time-zero run stopped early can export a clearly marked partial plotting
history; it can later resume with a later stop time under the same scientific policy.

For a cheap first-analysis test, `diagnostic_fork_checkpoint` names a closed,
matching checkpoint **before the first analysis**, from an observation-only-noise
parent experiment. The new experiment records the parent's run ID and manifest
hash. Normal ownership, member-count and layout-fingerprint checks still apply.
New noise starts at the fork boundary. Parent files are read only. The parent state
is copied into a new checkpoint namespace before advancement.

A fork is not an ordinary compatible restart and does not reproduce the changed
policies' history before the fork. It publishes only a checkpoint window, not a
time-zero aggregate. Set `mode3_export_ensemble_history: false`. An interrupted
fork may resume its own checkpoints with the same fork declaration and policies.
For `dt=.05`, checkpoint 298 is at year 14.95; stopping at 15.05 checks the first
analysis and the following forecast without integrating another 15 years.

## Prepare controlled cases

From `applications/icepack_model/examples/idealized_pig`, use the configuration
builder with the frozen 82-year stochastic-R base. This prepares files only:

```bash
python pace/prepare_diagnostic_config.py /path/to/frozen/params.yaml \
  --variant gain \
  --reference-data /path/to/completed/data \
  --reference-params /path/to/frozen/params.yaml \
  --fork-checkpoint /path/to/history/steps/checkpoint_00000298 \
  --data-path "$PWD/pace_runs/gain-first-analysis/data" \
  --stop-time 15.05 --output "$PWD/gain-first-analysis.yaml"
```

Keep separate outputs for `gain`, then `balance`, then `forcing`, then `noise`.
The builder requires an explicit projection limit for balance, an explicit drift
for forcing, and explicit noise standard deviations/time/length for noise.
Use `--no-analysis` with the same variant to obtain a matched forecast control.
Member count, initial inputs, observation errors, seeds, timestep and forcing
horizon are inherited. `--help` lists the required controls; no calibrated values
are supplied. Submit the selected configuration only after reviewing it and
obtaining allocation approval, using `run_idealized_1.sbatch -F /path/to/case.yaml`.

First require analytical correctness and physical consistency. Then compare
pre/post-analysis and subsequent-forecast errors for thickness, surface, velocity,
active melt and flux/volume against the matched control. Use area-weighted and
floating/grounded-region diagnostics, uncertainty coverage and innovation
covariance. A finite result, narrow ensemble or improved observed surface alone
does not establish scientific acceptance. Test several seeds and longer held-out
windows before promoting one configuration to a full-duration run.


## Spread-relaxation diagnostic

The completed gain/balance case remains strongly underdispersed. At year 15,
after its first saved analysis, thickness RMS member spread is 0.01070 m while
mean RMSE is 8.670 m; at year 30 these are 0.002444 m and 14.074 m. Snapshots at
14.95 and 15 years are different physical times and cannot isolate an analysis
increment by themselves. Same-cycle stage reports address that limitation.

Prepare a paired fork from the same committed pre-first-analysis checkpoint:
one case keeps spread relaxation disabled, and one opts into 0.5 for thickness
and basal melt. Retain the stochastic filter, prior, observations, physical
balance, 1 m projection guard, solver settings, seeds and existing noise. Stop
at 15.05 years to inspect the analysis at 15 and the following forecast. The 0.5
coefficient is an explicit sensitivity choice requiring review, not a tuned
recommendation. The builder options are `--spread-relaxation`, `--spread-fields`
and `--stage-diagnostics`. Forks intentionally produce checkpoint windows rather
than a full history starting at year zero.

Before a longer test, verify baseline replay against its parent, equality of the
raw analyses between cases, preservation of the controlled mean before bounds,
measured spread, finite/admissible finalized fields, consistent diagnostic
velocity and bounded projection. Compare truth errors after physical balancing
and after the next forecast. Larger spread by itself is not an acceptance gate.
Additional seeds and region/area-weighted diagnostics remain necessary for
scientific validation.

RTPS is described in Houtekamer and Zhang (2016), section 4a(3),
[doi:10.1175/MWR-D-15-0440.1](https://doi.org/10.1175/MWR-D-15-0440.1),
following Whitaker and Hamill (2012). This implementation follows the sample-SD
relaxation formula and explicitly documents its prior-direction fallback at
exactly zero analyzed variance. Its default is disabled. Local numerical and MPI
checks do not establish PIG trajectory improvement; that requires an approved
model validation.
