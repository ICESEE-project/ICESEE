# A compact basal friction variation exercise

This supplementary case delegates to this application's own `run_da_issm.py`. It uses four members, two model years, 0.2-year steps and observations at years 1 and 2. The intended output has ten forecasts, eleven stored columns and two analyses. It preserves the reference mesh, forcing, parameter uncertainties, sliding law, inverse objective weights and friction bounds inherited from `../../params.yaml`. Fresh truth/background/ensemble generation is enabled so a new data directory does not depend on an earlier output run.

The tracked `data/ISMIP.reference_simulation_0.mat`, `Domain.exp` and `Mismip2.par` supply the scientific model inputs used by the initializer. Unlike the Choi example, this version's active callbacks generate their existing friction perturbations and do not read `bed_kriging_results.h5`; those bed-kriging reads are commented out. No replacement bed prior is fabricated.

Three child profiles use this application's case-sensitive `FrictionCoefficient` block. `params-fixed.yaml` freezes it during EnKF and disables the inverse step; `params-enkf.yaml` permits joint EnKF updates without inversion; `params-hybrid.yaml` freezes the EnKF friction update and enables the existing external inverse path from year 2. In this short case that is the final observation. The named fixed treatment is the prescribed prior, not a measured statement about its error relative to truth.

## Prepare and preflight

Use an isolated writable ICESEE/application checkout and an approved working MATLAB/ISSM environment with `ISSM_DIR`. Python must resolve that source and use a functional MPI/HDF5 stack. The wrapper creates unique relative output and ISSM-example names from `--run-id`; run variants sequentially because metadata files beside the application driver are still shared within a checkout.

`common.yaml` requests `model_nprocs: 1` and sets the existing `ICESEE_PERFORMANCE_TEST: true` resource-policy switch to retain that request. Its historical name does not make this a performance benchmark. Python rank count and internal ISSM process count remain separate; inspect the actual outgoing MATLAB/solver commands. Without the switch, the earlier local fixed-policy run requested one internal process but used three for reference/initialization and thirteen for forecasts on a twelve-core host.

From the repository root:

```bash
python applications/issm_model/examples/basal_friction_variation/tutorials/short_assimilation/run_case.py \
  --variant all --mode 1 --ranks 1 --run-id basal_student01 --check-only
```

This checks inherited YAML, required files, executables and fresh target paths. It does not start MATLAB, test licensing, execute a solve or establish that the legacy callbacks satisfy the intended clock. Missing prerequisites produce a nonzero status.

## Execute the selected legacy path

On permitted compute resources after preflight:

```bash
python applications/issm_model/examples/basal_friction_variation/tutorials/short_assimilation/run_case.py \
  --variant fixed --mode 1 --ranks 1 --run-id basal_student01
```

Choose `--variant enkf`, `--variant hybrid` or `--variant all` deliberately. After a functional two-rank probe, `--ranks 2` exercises member-parallel operation; `--mode 2` selects the file-backed generic route. `--mode 0` is restricted to one Python process. Mode 3 is deliberately absent: the legacy driver's YAML selects it, but that driver does not register a native runner.

The wrapper records requested settings, exact command, stdout/stderr and process exit under `run_records/<run-id>/<variant>`, while model output is `_modelrun_tutorials/<run-id>/<variant>`. Existing records, output and ISSM workspaces are refused. Choose a new ID instead of clearing failed results. The child driver remains responsible for model computation; the tutorial wrapper supplies no surrogate solver or synthetic result. After the child exits, it checks the expected ensemble/observation files and the ensemble/member/column and mean shapes. Missing or malformed output makes `tutorial_exit_code` nonzero even when an older application exception handler returned process exit zero. The separate `output_verification` is a structural completion check.

## Acceptance and current limitations

The [backend record](../../../../../../paper/production-report/issm-basal-results.json) identifies exact source/profile hashes and distinguishes historical runs from the current capped profiles. With runtime source `02ac5496d3d6cacf0ded737cd821e3e6f43c9d42`, the capped fixed, EnKF-only and hybrid profiles all passed their stated numerical/runtime checks on ISSM 4.24 and MATLAB R2025b; the hybrid also passed its inverse-handoff checks. All three used one Python rank in Mode 1 and one internal process for every recorded solver command. Their `(20082,4,11)` ensembles and `(20082,11)` means were finite, all state columns were populated, stored means matched member means, reference endpoints aligned, and two analyses occurred at years one/two and columns five/ten. Fixed and EnKF-only made no inverse calls, as configured.

The hybrid's final analysis made four actual inverse calls, one per member. Each returned 3,347 finite friction coefficients changed from forecast; the maximum absolute difference between returned and stored coefficients was zero. All four inverses stopped on `dxmin` during line search. This demonstrates hybrid execution and the coefficient handoff, without a claim of gradient-tolerance convergence, objective improvement or physical skill. The last printed trial cost is not established as the objective of the accepted final control. The earlier fixed pass used the uncapped 3/13 process counts and is retained separately.

The [runtime correction](https://github.com/ICESEE-project/ICESEE/pull/74), including commit `52c19fad4567148017ba6960cac58c62f7063a8d`, restores the application directory before HDF5 generation returns and exports the actual pre-solve reference as column zero. Post-step reference states map to columns one through ten, including the terminal state, and the basal callbacks now honor `initial_state_only`. The existing one-year prior spin-up remains separate from this assimilation clock. Scalar diagnostic files keep their existing post-step schema. Earlier attempts that exited zero before assimilation or produced shifted/empty reference columns are retained as unaccepted evidence.

The application still sets `HDF5_USE_FILE_LOCKING=FALSE` and retains older error handling in some callbacks. The wrapper's structural check cannot establish finite values, a correct physical clock, geometry, inverse calls or optimizer convergence. Use the closed-file verifier and inspect model/analysis records before accepting a result.

Inspect MATLAB errors, all member/variable arrays, physical time, initial condition, observation columns 5/10, geometry and inverse targets. In particular, establish the actual saved-state clock before accepting parity with another application. Record physical extrema and complete-member comparisons. A four-member exercise does not establish calibrated parameter uncertainty, convergence, restart correctness or scientific skill.

See the [unified sample route](../../../../../../docs/handbook/sample-runs.md), [ISSM coupling](../../../../../../docs/handbook/coupling-issm.md), [inversion](../../../../../../docs/handbook/inversion-and-inference.md) and the more extensively updated [Choi teaching protocol](../../../ISMIP_Choi/tutorials/short_assimilation/README.md) for the distinct paths and required interpretation.
