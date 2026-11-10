# Production operations and six month handover

Operate ICESEE as a collection of reproducible scientific cases. Keep an approved source/environment baseline available, submit new work into separate run directories, and record who accepts its physical results. A long successful queue run should not require access to the principal developer's shell history, private home directory, or unpublished local files.

## Responsibilities and decisions

| Role to assign locally | Decisions and records |
|---|---|
| Run operator | Select an approved case, preserve provenance, check resource/input readiness, submit, monitor, and run integrity verification. |
| Scientific reviewer | Accept observation timing/units, priors, physical parameters, solver settings, output bounds, posterior behavior and scientific tolerances. |
| Software maintainer | Review code/environment changes, run the relevant regression gates, maintain the wiki, and preserve a known working revision. |
| Backup operator | Independently reproduce the small case and the documented restart drill; be able to retrieve archived inputs/output. |

These are responsibilities, not named appointments. Record the actual people and local escalation route in the group's private operational record. Routine operations need not wait for the absent principal developer when they follow an accepted case and these roles have accepted responsibility. Scientific changes and unresolved integrity failures need the corresponding local reviewer rather than an improvised threshold.

Before the absence, designate one working production baseline and one representative acceptance case for each application actually planned. The initial baseline should include the exact source SHA, dependency environment/digest, input checksums, resource topology and expected verification results. Keep the baseline readable and avoid updating it in place. Documented support for PIG does not automatically cover a different glacier, mesh or observation operator.

## Protecting run directories

The configuration loader has import-time side effects. For a fresh run it calls `shutil.rmtree(data_path, ignore_errors=True)` on the detected root rank, then recreates the directory. It guards the application cwd, the user's home and filesystem root, but that is not protection for a valuable experiment directory. `restart_enabled: true` alone does not prevent this cleanup. A supported restart must explicitly request `resume_from_checkpoint: true` through YAML or CLI.

Use a new absolute `data_path` for every fresh calculation and a separate log/provenance directory outside it. The shipped PIG fresh launcher refuses an existing data directory and takes a launch lock. Keep that guard. Do not rerun a driver with its default `_modelrun_datasets` directory where earlier results are valuable, and do not import a live application merely to list settings.

Do not start a second writer in a directory used by a running job. Read only closed outputs for full validation. Copy or archive the complete needed tree before modifying retention, rebuilding aggregate output, or experimenting with restart. HDF5 file locking should remain enabled; disabling it is not a repair for uncoordinated readers/writers.

## Case record

Keep a run record with these fields:

| Category | Minimum contents |
|---|---|
| Source | Full Git SHA, branch, remote, clean/dirty status, binary diff and hashes of untracked code/config needed by the run |
| Environment | Python executable/version, MPI library and launcher, PETSc/Firedrake/Icepack or ISSM/MATLAB versions, h5py parallel capability, Spack lock/environment commit or container digest |
| Inputs | Mesh, spin-up, observations, forcing, all inherited config files, paths and content checksums |
| Resolved experiment | Actual time vector, observation indices/requested times, variable ordering/dimensions, ensemble size, seed/noise policy, analysis/localization/inference settings |
| Launch | Exact command, scheduler script, job ID, nodes/ranks/threads, ranks per model, working directory and output path |
| Completion | Scheduler parent state, each model/MPI step exit, final committed step and analysis count, verifier exit and machine-readable report |
| Acceptance | Numerical comparison tolerances and rationale, physical diagnostics, reviewer, date and disposition |

Record ordinary version output and approved environment manifests; do not dump all environment variables because they may contain credentials. `git diff` does not include untracked files, so archive any untracked inputs/config/code explicitly.

## Before submission

Activate the approved [deployment environment](deployment.md) and verify the Python/launcher pair functionally. Modes 1/2 need an MPI-enabled h5py for shared `mpio` paths. Verify external input paths and permissions using the scheduled job's environment, not only an interactive login shell. A source-only update normally does not require rebuilding the entire scientific stack, but ABI compatibility and the regression gates still matter.

Calculate full-history and checkpoint capacity, include duplicate export/temporary space, and check filesystem quotas/inodes as well as free bytes. Choose whether retrospective member spread is necessary. For full native PIG aggregate export, every completed cycle must be written and retained (`checkpoint_every=1`, `checkpoint_keep_last=0`); keeping only two restart states cannot reconstruct all history later.

Keep threads explicit, usually one BLAS/OpenMP thread per MPI rank for the provided launchers. Confirm that `world_size` factors into the requested native topology and that memory fits active model members plus analysis/noise workspaces. PIG's shipped scheduler script contains site/account defaults that must match the local allocation; it is a reference, not a portable account assignment.

Run the small case after changing source/environment and before consuming a production allocation. For a scientific configuration change, run a reduced target-model case and obtain an agreed acceptance tolerance before the long run.

## Monitoring and completion

Follow initialization traces, model-step progress, analysis-event JSON, checkpoint commits, memory, CPU and I/O activity. A silent period during code generation or filesystem work is not by itself a deadlock. Compare with the expected phase and inspect rank stack diagnostics where enabled. The PIG launcher can emit timed hang diagnostics; preserve them with stdout/stderr.

On Slurm, inspect both parent and step records with `sacct` after termination. A PDE/MPI step can exit zero while its containing allocation times out during a later verification read. Conversely, a scheduler `COMPLETED` state does not prove output shape, observation timing or scientific validity. Record simulation, finalizer and verifier status separately.

Before accepting output, require:

1. The intended simulation reached its final state and each required step exited successfully, with no traceback/fatal-error markers.
2. The physical time vector has `nt+1` columns including the initial state, the final committed cycle is `nt-1`, and the expected analysis times occurred exactly once.
3. Expected members, variable blocks, ownership ranges, shards, means and required profile files exist with valid shapes and finite values.
4. Means recomputed from member values agree within the declared tolerance. Validate complete history when the experiment depends on it, using bounded slabs rather than an unbounded read.
5. Physical/scientific diagnostics pass their predeclared criteria. Finite values alone do not establish acceptable melt, geometry, velocity, spread or inference accuracy.

The PIG timing verifier in `tools/verify_timing_outputs.py` is tied to particular one/five-year validation schedules and five-field/profile shapes. Use it for those protocols; do not pass a different experiment and assume its hard-coded expectations form a universal verifier.

## Checkpoint inspection and restart

After MPI has stopped, inspect native committed checkpoints without loading model fields:

```bash
python -m ICESEE.src.parallelization.mode3_checkpointing \
  /absolute/run/data/_mode3_state_history
```

Only a committed manifest plus every referenced shard is a restart state. Hidden temporary artifacts, an incomplete manifest, a plotting profile or a VDS are not interchangeable with it. The inspector reports rejected checkpoints and the next cycle. A checkpoint for completed cycle `k` resumes at `k+1`; its analysis is not repeated. The initial checkpoint has `completed_step=-1` and resumes at cycle zero.

For a validated PIG resume, retain the same source or compatible verified implementation, input values, physical setup, mesh/layout, ensemble size, schedule, time grid and stochastic policy. Use the same `data_path` with `--resume_from_checkpoint=True`, and preserve a copy of its previous launch/logs. The common fresh-run scheduler launcher intentionally refuses existing output; make a dedicated resume launcher following the PIG README. Do not remove its fresh-data guard and thereby make ordinary launches destructive or ambiguous.

Native checkpoint metadata rejects incompatible schedule/layout/policy identities. Do not delete or edit identity fields to bypass that check. A different decomposition or environment is a portability experiment requiring a comparison, not a routine continuation. The shared whole-member native policy explicitly rejects resume; use only the application whose runner has a validated binding.

If a run was interrupted before committing its next state, the last complete state is restored and that incomplete cycle is recomputed. Deterministic native noise reconstructs AR history from seeds and schedule. Keep checkpoint retention sufficient to tolerate one corrupt/rejected newest state and verify both a normal restart and a restart crossing an analysis event before production dependence.

## Aggregate export and recovery

PIG's closed native checkpoints can be exported separately using bounded rows:

```bash
python -m ICESEE.src.parallelization.distributed_history_export \
  /absolute/run/data --row-chunk-size 8192
```

A recovery view from a stopped incomplete run must be explicitly marked partial:

```bash
python -m ICESEE.src.parallelization.distributed_history_export \
  /absolute/run/data --row-chunk-size 8192 \
  --allow-partial --simulation-exit-code 1
```

Replace `1` with the actual saved simulation exit code. The exporter requires the initial checkpoint and consecutive committed cycles through the last included cycle; missing middle snapshots are not silently interpolated. The PIG workflow also needs its reference layout. It atomically publishes an aggregate after successful construction and preserves checkpoints and an existing published file on failure. A partial file is a diagnostic/recovery artifact and must not be reported as a completed experiment.

Mode-2 VDS output depends on its source shards and their paths. Archive/move the full dependency tree and verify the moved copy before deleting originals. The native exported aggregate has a different construction and completeness metadata; do not infer its independence from the name `icesee_ensemble_data.h5` alone.

## Failure triage

| Symptom | First independent checks | Safe next action |
|---|---|---|
| Every process reports rank zero | Actual launcher/mpi4py probe and MPI library | Select the matching MPI launcher/environment, rerun the smoke case. |
| HDF5 lacks `mpio` | `h5py.get_config().mpi` in the job Python | Use the approved parallel-HDF5 build for the selected path. |
| Hang near file creation/analysis | Rank stacks, collective participants, open reader/writer phase | Preserve logs; reproduce a small case. Do not disable locking. |
| Wrong state size only when ranks exceed members | Ownership declaration, ghost exclusion, resource plan | Fix the adapter contract and rerun the topology matrix. |
| Collapsed-spread rejection | Realized prior/noise, observed-space anomalies | Review spread configuration scientifically; retain the guard. |
| Resume rejected | Source/layout/run ID, schedule and noise/time-grid identities | Preserve the run and compare provenance; use fresh output for a new experiment. |
| Parent allocation times out after simulation ends | MPI step exit and finalizer/verifier log | Verify closed output separately within an adequate allocation; report both states. |
| Float32 output fails strict parity | Storage dtype and declared tolerance | Repeat strict parity in float64 or justify a separate scientific tolerance. |

For an issue, capture a small config/input case, exact source/environment, command, final progress/checkpoint, error and relevant rank stacks, expected behavior and comparison result. Keep credentials and large private input datasets out of public issue text. Do not fix a solver failure by changing physical parameters or raising iteration limits without recording the scientific/numerical change.

## Maintenance during the six months

At each production batch, verify provenance, capacity, accepted baseline and backup access. After each run, record simulation/verification/scientific disposition and archive needed shards/inputs. Weekly, operators should review failed or stalled cases and storage growth, with no automatic promotion of unreviewed results. Monthly, the backup operator should reproduce the smoke case and inspect a representative checkpoint archive. At each source/dependency update, the maintainer should run relevant regression/parity/restart gates and update the handbook/wiki in the same change.

These are proposed operating practices, not scheduled automations. Avoid routine environment upgrades in the middle of a long experiment series unless needed to fix a verified issue. Keep the prior environment and source recoverable. During handback, transfer the case ledger, unresolved issues, accepted changes, archived manifests and reproducibility notes so the returning developer can understand what ran without rerunning everything.
