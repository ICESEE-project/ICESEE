# Forecast recovery and checkpointed continuation

ICESEE separates numerical step rejection from allocation-level continuation.
Recovery preserves the selected stochastic EnKF, observations, process-noise
schedule and external time grid. Finite values alone do not establish accuracy.

## Numerical transaction

`src/parallelization/forecast_recovery.py` provides `advance_with_recovery`.
Adapters supply capture, restore, advance and validate callbacks on their spatial
communicator. Capture includes owned physical fields, mutable auxiliary state
and adapter-owned random generators. Restore synchronizes halos and discards
failed solver caches. Trials publish no diagnostics or checkpoints. The
controller also restores Python and NumPy legacy RNG states.

After a recoverable rejection, all spatial ranks restore the cycle's starting
state and retry the full interval with twice as many uniform substeps. Maximum
halvings, minimum internal dt and maximum substeps bound the work. Times are
computed from the external interval. An optional model-provided `stable_dt`
bound is reduced across the spatial group and sets the initial subdivision.
The bound's physical derivation belongs to the adapter; no universal ice-sheet
CFL formula is assumed.

The native PIG adapter runs this controller before process noise. It restores
its owned registry and rebuilds FlowSolver with the existing solver settings.
Substeps evaluate forcing at their physical start times, keeping the original
dt and horizon in configuration. Estimated-melt drift uses internal duration.
Profiles and process noise remain outside physical trials. Firedrake
`ConvergenceError` and non-finite forecast fields are recoverable; other
exceptions are fatal. A model stuck inside MPI cannot participate in rejection
agreement and cannot be repaired by this controller.

The policy is off by default. Selecting it changes checkpoint identity. PIG
integrates it only in native Mode 3. The ISSM integration is described below;
other adapters and reference generation need their own verified rollback
integration. An analysis rejection
is not retried as a forecast: recovery does not redraw observations, repeat
assimilation, damp an analysis, lower a physical threshold or change tolerances.
It preserves the last committed state for a reviewed correction.

## Committed boundary stops

`run_control.CycleStopController` observes a declared cycle budget, elapsed-time
budget or signal. All ranks agree after a cycle checkpoint commits. PIG writes
`_continuation/segment.json` atomically with next cycle, analysis count, model time
and reason. Intermediate segments defer the full aggregate; the final successful
segment exports it. Duration, dt and time columns remain unchanged.

Signals are handled while the cycle controller is active. Root-only initialization
and final export can outlast an inadequate reserve. Allow for initialization,
the slowest cycle, checkpoints and export. A hard kill or indefinitely blocked
collective uses the existing manifest-based restart path.

## Submission transactions

`checkpoint_continuation.tick` and PIG `tools/continue_case.py` provide a bounded
supervisor. Preparation freezes resolved YAML inheritance. A source fingerprint
includes tracked edits and untracked source, excluding known run products.
Input identity remains governed by reference/time-grid validation and native
restore; keep input files immutable.

Each submission has a unique job name and an intent record durable before
`sbatch`. After a communication failure, the supervisor reconciles that name
against the user's `squeue` and `sacct` records. Unknown outcomes are held rather
than resubmitted. A directory lock prevents two supervisors. After a supervisor
crash, inspect its recorded lock owner before manually removing a stale lock.

A successor is eligible after an exit-zero completed segment, or a timeout,
node failure or preemption with new committed progress. It retains the original
duration, dt and science, and disables reference regeneration and ensemble
initialization on resume. Scientific failure, cancellation, out-of-memory,
missing progress, changed source/configuration or run identity stop the chain.
Maximum submissions bounds allocations. The supervisor checks closed manifest
metadata and shard existence; native restore validates HDF5 contents. Completion
does not constitute scientific acceptance.

Keep the supervisor running or invoke it again. No daemon or unbounded dependency
chain is installed. Follow site rules for submission hosts, polling and budgets.
Source and dependency roots remain separate; no reinstallation is needed.

## Validation limits

Tests cover rollback after partial progress, RNG restoration, external event
boundaries, limits, fatal errors, PIG forcing time, unchanged solver options,
signals, atomic reports, pins, missing shards, duplicate prevention and uncertain
submissions. Four-process MPI exercises one-rank rejection and collective stops.
The checkpoint failure test injects ENOSPC on a non-root rank, keeps the previous
commit, then resumes with two processes and verifies exact owned-state recovery
and analysis progression.

PACE validation on 24 November 2026 used two members, four ranks and two spatial
ranks per model group, retaining the existing solver settings and dt=0.05 years:

| Check | Evidence | Result |
|---|---|---|
| Physical recovery | Job 16085734, exit 0, 1 min 58 s | No-rejection forecast matched the original path bitwise; an injected rejection matched explicit half-steps bitwise; direct physical checkpoint restore matched bitwise. |
| Common launcher and boundary continuation | Job 16085735, all three steps exit 0, 5 min 46 s | Continuous two-cycle forecast versus a stopped/resumed forecast: maximum absolute owned-field difference 4.547473508864641e-12; stored melt fields matched bitwise. |
| Local fast suite | 983 passed, 9 skipped, 1 deselected | Included scheduler transaction, numerical rollback, noise, analysis and plotting regressions. |

The physical test is reproducible with
`src/tests/parallel_mpi/_pig_forecast_recovery_worker.py` in the normal PIG
environment, using a separate data directory and the lab inputs. It injects a
failure for validation, not into production forecasts. Recovery required typing
the unchanged velocity threshold as a Firedrake Constant for rebuilt solver
contexts; tolerances and iteration settings were retained.

These short checks do not validate an 82-year recovery-enabled trajectory or
analysis/noise equivalence through every event. The launcher comparison ended
at year 0.1, retained the original 82-year reference horizon and encountered no
analysis events. Evaluate all members, physical bounds, observations, diagnostics
and long uninterrupted/resumed histories before accepting a scientific result.

## ISSM whole-member recovery

The `ISMIP_Choi` and `basal_friction_variation` forecast callbacks accept an
explicit `forecast_recovery` policy. A transaction preserves the packed member,
mutable MATLAB model caches, per-member output files and MATLAB RNG. A rejected
trial restores them before retrying the entire external interval with smaller
internal steps. Each trial uses a private ISSM execution directory. Process
noise and analysis remain outside the physical forecast transaction.

The MATLAB bridge reports identified numerical rejections through
`ICESEE:ISSM:ForecastRejected`. The checked forecast rejects non-finite fields,
speeds above the existing 5000 m/year guard, or a returned time that differs from
the requested interval end. Undefined functions, missing files, command timeouts
and unclassified ISSM errors are fatal. They do not trigger repeated solves or
consume stale HDF5 output. Solver tolerances, iteration settings and physical
bounds are unchanged. A MATLAB process that crashes or remains blocked cannot be
repaired by this numerical transaction.

`ISMIP_Choi` native Mode 3 additionally supports committed boundary stops and
complete backend resume. Select `issm_native_continuation: true`, or a recovery
policy, with `checkpoint_every: 1`. Backend MAT files, output snapshots and RNG
are committed before the vector checkpoint that references their checksums.
Resume restores both; a vector-only historical checkpoint is refused. Resident
whole-member execution and the original Python rank count are required. This
does not expose ISSM's internal mesh partition to Python, validate streaming
members or provide rank-topology migration.

On 24 November 2026 a local MATLAB R2025b/ISSM 4.24 check used the packaged Choi
mesh (3347 vertices, 6319 elements), two members, one Python rank and one internal
ISSM rank. The configured reference horizon was two years with dt=0.2 years;
the forecast comparison advanced two cycles to year 0.4. A local-only bed fixture
duplicated the packaged reference bed for the two members; it was not a calibrated
prior. Analysis and inversion were disabled to isolate forecast/restart behavior.

| Check | Result |
|---|---|
| Recovery enabled, no rejected solve | All six fields of both members matched the baseline bitwise. |
| Stop at year 0.2, restart in a fresh MATLAB process, advance to 0.4 | All six fields of both members matched the uninterrupted case bitwise. |
| Inject a numerical rejection after a real solve and MATLAB RNG draw | Retried two 0.1-year substeps matched an explicit half-step reference bitwise; saved MATLAB RNG states matched. |

These are local physical forecast checks, not assimilation/inversion acceptance,
HPC scaling or long-run scientific validation. The common continuation supervisor
can inspect the checkpoint metadata, but an ISSM-specific Slurm preparation and
resubmission workflow has not been validated. The PIG launcher must not be used
to submit an ISSM case.

## Related model practices

- [WRF running guide](https://www2.mmm.ucar.edu/wrf/site/documentation/users_guide/running_wrf.html): adaptive timestep and event/output alignment.
- [CIME running guide](https://esmci.github.io/cime/versions/master/users_guide/running-a-case.html): explicit continuation and resubmission controls.
- [PETSc nonlinear-failure limit](https://petsc.org/release/manualpages/TS/TSSetMaxSNESFailures/): bounded failure handling.

ICESEE still requires model-specific restoration and stability criteria; these
references do not validate a particular PIG timestep or recovery setting.
