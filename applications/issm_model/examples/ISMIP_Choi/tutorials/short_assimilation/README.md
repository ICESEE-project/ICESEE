# A short ISSM assimilation and friction-inference exercise

This teaching case delegates to the existing `ISMIP_Choi/run_da_issm.py`. It uses four ensemble members, two model years, the existing 0.2-year step and observations at years 1 and 2. Thus the intended history has ten forecasts, eleven saved state columns and two analyses at state columns 5 and 10. It reuses the prepared Choi mesh and physical inputs; it does not construct a new coarse mesh or change the sliding law, inversion objective, physical forcing or prior amplitudes.

The profiles inherit the real `../../params.yaml` through `common.yaml`. They reduce ensemble size/run duration and choose storage/resources for teaching. These are execution exercises, not a controlled accuracy or scaling experiment. This Choi case has not been executed with the MATLAB/ISSM backend in this documentation revision; the separate basal-friction backend record does not validate its objective, geometry handoff or native route. No Choi output values or scientific improvements are asserted.

## Compare three parameter-update mechanisms

| Profile | EnKF friction block | ISSM inverse step |
|---|---|---|
| `params-fixed.yaml` | Restored from forecast through `frozen_analysis_vars: [coefficient]` | Disabled; friction remains the configured prior/prescribed treatment |
| `params-enkf.yaml` | Joint covariance updates permitted | Disabled |
| `params-hybrid.yaml` | Frozen during EnKF | Enabled at the retained two-year inversion-start threshold, i.e. the final analysis of this short case |

All inherit the same reference, observations and physical setup. Other inherited bed/inference controls remain active according to the base YAML. “Fixed” is the prescribed prior-friction treatment often called WBF in longer experiments; this wrapper does not measure how wrong that prior is. The hybrid case is not a pure EnKF update of friction. Read [the coupling guide](../../../../../../docs/handbook/coupling-issm.md) and [inverse-problem guide](../../../../../../docs/handbook/inversion-and-inference.md) before interpreting differences.

## Prepare the actual model environment and inputs

Activate the approved numerical environment so MATLAB, MPI and ISSM are usable, and set `ISSM_DIR`. Use an isolated writable application checkout named `ICESEE`. Each run also receives a unique `$ISSM_DIR/examples/<example_name>` workspace. Separate output paths alone do not isolate every metadata/IPC file for concurrent ISSM applications; run these variants sequentially as the wrapper does, or use separate application workspaces.

`common.yaml` requests `model_nprocs: 1` and sets `ICESEE_PERFORMANCE_TEST: true`, the existing resource-policy switch that retains the requested count. This historical name does not turn a teaching run into a performance benchmark. Python MPI ranks and internal ISSM solver processes are distinct resources. The generic policy can otherwise redistribute the requested internal count; inspect the actual outgoing solver commands in every accepted run.

Required inputs include `Domain.exp`, `Mismip2.par`, `data/ISMIP.reference_simulation_0.mat` and a prepared `bed_kriging_results.h5` in the Choi directory. The reference `.mat` is tracked in the source snapshot; the bed-kriging file is not present in that inspected checkout and must come from the approved scientific input bundle. It must match the reference mesh and contain the bed ensemble expected by the existing callbacks. The wrapper neither fabricates a replacement prior nor downloads a substitute dataset.

Start with preflight from the repository root:

```bash
python applications/issm_model/examples/ISMIP_Choi/tutorials/short_assimilation/run_case.py \
  --variant all --mode 1 --ranks 1 --run-id student01 --check-only
```

Preflight resolves the inherited YAML through the existing safe loader and checks file/executable paths without importing the application configuration or starting MATLAB. It exits nonzero when prerequisites are missing. A successful preflight is not a MATLAB-license, MPI, PDE or scientific validation.

## Run the prepared exercise

On permitted compute resources, use one Python process initially:

```bash
python applications/issm_model/examples/ISMIP_Choi/tutorials/short_assimilation/run_case.py \
  --variant all --mode 1 --ranks 1 --run-id student01
```

After the appropriate functional MPI/HDF5 probes, use `--ranks 2` for a two-rank Mode-1 exercise, or choose `--mode 2` for the file-backed route. `--mode 0 --ranks 1` selects serial orchestration. This wrapper deliberately offers Modes 0/1/2; Choi's native whole-member route has separate support and does not provide a verified checkpoint resume workflow. It is not silently substituted here.

The wrapper retains requested YAML/settings, exact commands, stdout, stderr and process exit status under `run_records/<run-id>/<variant>` (or `--record-root`). The assimilation output remains in unique relative `_modelrun_tutorials/<run-id>/<variant>` paths under the Choi directory. Existing output, ISSM workspaces or record directories for that ID are refused. Choose a new ID for another attempt and preserve failures.

## Inspect results before accepting them

Check model-time and observation metadata, all four members and six named nodal blocks, finite fields, mean consistency, analysis events and the final model time. Inspect MATLAB logs and returned HDF5, because an outer legacy driver error path can print an error without establishing a trustworthy successful process status. Verify the observed velocity targets and the friction-return contract separately for the hybrid case. Report shape/time checks and measured differences without assuming four members provide a calibrated posterior.

The source-derived expectations above are not stored Choi run results. After execution, `execution-status.json` records the child `simulation_exit_code`, wrapper `tutorial_exit_code` and `output_verification`. The wrapper refuses missing ensemble/observation files or unexpected ensemble/mean dimensions even if an older outer driver returned zero after an error. Passing this structural check does not establish finite fields, a correct physical clock, geometry, inverse calls or convergence. Perform and record the application checks before accepting output. No Choi restart, convergence, performance speedup or physical-skill result is supplied with this unexecuted backend exercise.

Continue through the [unified sample route](../../../../../../docs/handbook/sample-runs.md), [ISSM coupling](../../../../../../docs/handbook/coupling-issm.md) and [inversion guide](../../../../../../docs/handbook/inversion-and-inference.md) for closed-file checks and the scientific interpretation of the three policies.
