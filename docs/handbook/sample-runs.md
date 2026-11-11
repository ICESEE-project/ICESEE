# Sample runs and the evidence to retain

We begin with a small calculation whose expected dimensions, observation times and active callbacks are explicit. These exercises connect the [parallelization modes](parallelization-modes.md), [Icepack coupling](coupling-icepack.md), [ISSM coupling](coupling-issm.md) and [inversion guide](inversion-and-inference.md) to executable cases. They supply a reproducible teaching route before changing mesh resolution, ensemble size or parameter controls.

An executable command, an executed model and a scientifically accepted result are different evidence. The Lorenz endpoint, coarse Icepack and all three capped basal-friction profiles below have recorded local runs. The Choi case remains unexecuted. Do not interpret a configuration or source check as an executed PDE calculation. The [readiness record](production-readiness.md) retains the broader production decisions.

## 1 Establish the environment and run identity

Activate the numerical stack described in [deployment](deployment.md), verify the actual ICESEE import location, and record the source commit, resolved YAML, parent YAML files, input checksums, Python/MPI/HDF5 versions and exact launch command. Use a new output ID for a fresh attempt. The short ISSM wrappers refuse existing model output, record directories and ISSM workspaces. Run variants sequentially in an isolated writable application checkout because the legacy MATLAB bridge also writes metadata beside the driver.

Choose resources according to the installed adapter. Python MPI ranks, native spatial ranks and an ISSM internal solver's `model_nprocs` are distinct controls. The supplied ISSM profiles request `model_nprocs=1` and set `ICESEE_PERFORMANCE_TEST: true` to retain that request through the existing Mode-1 resource policy; this historical switch name does not make the tutorial a performance benchmark. Without that switch, an earlier one-Python-rank fixed run requested one internal process but used three for reference/initialization and thirteen for forecasts on a twelve-core host. Inspect the actual outgoing solver counts for each run. Python ranks between one and four require a corresponding allocation and a functioning launcher; only the explicitly recorded configurations have backend evidence.

Do not import the full application configuration merely to inspect it: fresh-run setup can remove an existing output directory. The ISSM wrapper's `--check-only` path uses the side-effect-free YAML loader and does not start MATLAB or a model solve.

## 2 Begin with the Lorenz endpoint case

The [endpoint experiment](../../paper/production-report/experiments/README.md) supplies [lorenz-endpoint.yaml](../../paper/production-report/experiments/lorenz-endpoint.yaml), exact commands and [recorded results](../../paper/production-report/experiments/mode-smoke-results.json). Although its historical directory is `lorenz96`, this fixture uses the repository's current three-variable Lorenz-63 equations. It has eight members, fifty forecast intervals and fifty-one stored state columns, including the initial condition. Observations occur at model times 0.2, 0.4, 0.6, 0.8 and 1.0; their state columns are 10, 20, 30, 40 and 50.

The recorded local runs used two Python MPI ranks for Modes 1, 2 and 3. All produced finite `(3,8,51)` histories with means consistent with their members. Modes 1 and 2 have a recorded complete-member comparison: maximum absolute difference `5.3071289585204795e-9` over the history and `9.100062925426755e-10` at the final state, within `rtol=atol=1e-7`. Mode 3 completed five analyses; it is a completion/integrity exercise, without a claim of three-mode member equality. Its initialization and stochastic policy must be accounted for before such a comparison.

Use this small case to learn the forecast-count versus stored-column distinction and confirm that the terminal observation is assimilated. Tiny local timings from this fixture are not scaling results. The paper's experiment README gives the launch form and a new temporary directory for each invocation; follow it before adapting the MPI command to a cluster.

## 3 Advance to a compact Icepack calculation

The new [coarse assimilation case](../../applications/icepack_model/examples/synthetic_ice_stream/tutorials/coarse_assimilation/README.md) delegates to the existing synthetic ice-stream driver and model hooks. It requires no external mesh or spin-up dataset. Its generated rectangle has six by four quadrilateral elements, degree-one fields and four members over one model year with ten forecast steps. The ordered blocks are `h,u,v,smb`; the expected one-spatial-rank state has 140 rows and history `(140,4,11)`. Observations at 0.5 and 1.0 year map to columns five and ten, including the terminal state. The child YAML disables optional inference/localization for this first exercise and explicitly clears the parent's late-start override.

From the repository root, set a fresh absolute output path and inspect the merged configuration without importing a solver:

```bash
CASE=applications/icepack_model/examples/synthetic_ice_stream/tutorials/coarse_assimilation
python "$CASE/run.py" --mode 2 --output /absolute/new/icepack-coarse-mode2 --check
```

After activating the verified Icepack/Firedrake environment and functionally testing its launcher, execute:

```bash
mpiexec -n 2 python "$CASE/run.py" \
  --mode 2 --output /absolute/new/icepack-coarse-mode2
```

The wrapper refuses existing output and preserves the existing application's working-directory discovery. It offers generic Modes 0/1/2 and the registered whole-member native Mode 3, with one model rank per member. Mode 0 uses one Python process; use a fresh output directory for every other mode. This exercise is not a multi-rank spatial validation.

The [recorded Icepack results](../../paper/production-report/experiments/icepack-coarse-results.json) cover actual runs on 10 November 2026: one Python rank for Mode 0 and two MPI ranks for Modes 1–3. All completed ten forecasts and two analyses, with observations at columns five and ten. Their packed state series have shape `(140,4,11)`, finite values and positive thickness; generic stored means are consistent with member means. Modes 1 and 2 agree over the complete history to a maximum absolute difference of `2.007709554163739e-10`, within `rtol=atol=1e-7`. This is a recorded generic-mode comparison and four-mode completion/integrity result, without a claim of four-mode member parity, calibrated uncertainty or scientific skill.

Native Mode 3 writes `native_analysis_0000.h5` through `native_analysis_0010.h5`, each containing a `(140,4)` `states` dataset, plus `_mode3_state_history` manifests and shards. This binding does not export the generic ensemble file or stored `ensemble_mean` by default. Preserve both native output families and use the native verifier below. Its format and restart capability differ from the spatial PIG binding.

The result record retains configuration identity, numerical versions, analysis events and physical block extrema. The first local attempt exposed an MPI compiler mismatch: inherited `mpicc` selected MPICH while PETSc used Open MPI. Successful runs selected the matching compiler in the child environment. Importing Firedrake alone had not established a coherent compile stack; the installed environment was unchanged.

Read the [Icepack coupling guide](coupling-icepack.md) before selecting another generic or native route. Firedrake owns mesh, finite-element fields and communicator-specific solver objects. A shortened time interval or coarser mesh does not make those ownership and physical-finalization contracts optional. The coarse solution and reduced perturbations are illustrative settings rather than calibrated physical recommendations.

Retain physical field extrema by variable, solver errors and observation columns along with complete member arrays. PIG's specialized spatial native adapter and the synthetic ice-stream whole-member native adapter have different history, spatial-memory and restart capabilities. Do not transfer a passing result between them from their shared Icepack dependency alone.

## 4 Compare three friction policies with ISSM

Two new teaching protocols delegate to the existing application drivers:

| Application | Protocol and scientific inputs | Additional requirement |
|---|---|---|
| Basal friction variation | [Short assimilation exercise](../../applications/issm_model/examples/basal_friction_variation/tutorials/short_assimilation/README.md); tracked reference MATLAB model, domain and `.par` file | Working MATLAB/ISSM and the configured Python dependencies; active initialization does not read a bed-kriging file |
| ISMIP Choi | [Short assimilation exercise](../../applications/issm_model/examples/ISMIP_Choi/tutorials/short_assimilation/README.md); existing reference experiment and inference configuration | Approved `bed_kriging_results.h5` in the application's expected location, as well as MATLAB/ISSM |

Each uses four members, two model years and 0.2-year intervals: the intended schedule is ten forecasts, eleven state columns and observations at years one and two (state columns five and ten). They inherit physical parameters, reference mesh, uncertainties, forcing, sliding law, objective weights and friction bounds. The basal case explicitly regenerates truth/background/observations in its fresh output directory. The Choi case retains its existing bed-survey meaning at observation column two. Neither substitutes a fabricated bed prior or surrogate solver.

`params-fixed.yaml` prescribes the prior coefficient by freezing its EnKF update and disabling the inverse step. `params-enkf.yaml` permits joint state–friction EnKF updates with inversion disabled. `params-hybrid.yaml` freezes the EnKF coefficient update and enables the application's external inverse path from year two, the final observation in this short schedule. The fixed treatment is a prescribed prior; the file name does not demonstrate that it has a particular error relative to truth. State-block names differ: Choi uses `coefficient`, while the basal application uses `FrictionCoefficient`.

From the repository root, preflight the more self-contained basal exercise:

```bash
python applications/issm_model/examples/basal_friction_variation/tutorials/short_assimilation/run_case.py \
  --variant all --mode 1 --ranks 1 --run-id basal_student01 --check-only
```

After a working backend and approved compute resources are available, run one treatment:

```bash
python applications/issm_model/examples/basal_friction_variation/tutorials/short_assimilation/run_case.py \
  --variant fixed --mode 1 --ranks 1 --run-id basal_student01
```

Choose a fresh ID for a repeated attempt. `--variant all` runs the three treatments sequentially. Use `--ranks 2` only after a functional two-rank stack probe; `--mode 2` selects the generic file-backed runner. `--mode 0` is limited to one Python process. These tutorials deliberately do not launch Mode 3: the basal driver lacks registration, and Choi's separate whole-member native path needs its own lifecycle/analysis validation.

Both tutorials passed inherited-YAML resolution, Python syntax and negative preflight checks during preparation. The later basal clock correction exports the actual pre-solve reference as column zero, preserves post-step states in columns one through ten, and honors `initial_state_only`. Keep the existing one-year prior spin-up separate from this assimilation clock and inspect the saved time/observation metadata. Read [inversion and inference](inversion-and-inference.md) for the objective, geometry handoff and generic/native return-field differences before interpreting the hybrid treatment.

### Recorded basal backend result and its limits

The [public basal result record](../../paper/production-report/issm-basal-results.json) identifies runtime source `02ac5496d3d6cacf0ded737cd821e3e6f43c9d42`, the exact profile/source hashes and ISSM 4.24/MATLAB R2025b environment. The capped fixed, EnKF-only and hybrid profiles all passed their stated numerical/runtime checks in Mode 1 on one Python rank; the hybrid also passed its inverse-handoff checks. Every recorded outgoing solver command used one internal process. All three `(20082,4,11)` ensembles and `(20082,11)` means were finite, all stored state columns were populated, stored means matched member means, and analysis records occurred at years one/two and columns five/ten. The fixed and EnKF-only runs made no inverse calls, as configured.

In the hybrid run, all four members invoked the external inverse at year two. Each returned 3,347 finite coefficient values, changed from its forecast input, and the returned coefficient matched the final stored ensemble block exactly (maximum absolute handoff error zero). This establishes the executed hybrid lifecycle and coefficient handoff.

Every inverse stopped on `dxmin` during line search. This does not establish gradient-tolerance convergence, objective improvement, calibrated parameter uncertainty or physical skill. Printed objective rows include trial/function evaluations; the last printed trial is not established as the accepted final-control objective. A historical uncapped fixed pass is retained separately with its observed 3/13 internal process counts. None of these basal results validates the unexecuted Choi objective/geometry handoff or a native ISSM route, multi-rank ICESEE scaling or coupled restart.

## 5 Inspect the closed history and retain the result

Wait for the simulation, MPI workers and HDF5 writers to finish. The [sample verifier](../../scripts/docs/verify_sample.py) reads bounded row/time slabs, checks dimensions and finite values, compares stored means with member means, and verifies that referenced VDS files/datasets exist. For an ISSM tutorial history, replace the path below with the actual unique output recorded by the wrapper:

```bash
python scripts/docs/verify_sample.py \
  path/to/icesee_ensemble_data.h5 --members 4 --columns 11
```

For the Lorenz fixture use `--members 8 --columns 51 --rows 3`. The coarse Icepack generic history uses `--members 4 --columns 11 --rows 140`. Its native output requires the directory form:

```bash
python scripts/docs/verify_sample.py /absolute/new/icepack-coarse-mode3 \
  --native-analysis-series --members 4 --columns 11 --rows 140
```

That checks all eleven `states` snapshots for shape and finite values; a stored-mean check is inapplicable because this layout has no `ensemble_mean`. Verify the separate native history manifests/shards as well. Both verifier forms assess integrity. They cannot infer the physical units of packed blocks, establish a correct time grid, count accepted analyses, evaluate inverse convergence or measure scientific skill. A finite VDS fill value is not accepted as evidence of a present backing file.

Retain the wrapper's `requested-run.json`, exact command, stdout/stderr and `execution-status.json`, plus the integrity report, physical time metadata, observation/analysis records and input/source identity. Verify the initial state, terminal observation, geometric admissibility and parameter/inverse targets explicitly. Record whether the result is blocked, completed, integrity-checked or scientifically accepted, and the evidence supporting that statement. Restart and uncertainty calibration are separate experiments; neither follows from a four-member short run.

The [current-developments guide](current-developments.md) explains optional SMB/bed postprocessing, masks, freezing, caches and diagnostic records that can change or qualify these interpretations. Enable one additional mechanism at a time and retain the complete effective configuration.
