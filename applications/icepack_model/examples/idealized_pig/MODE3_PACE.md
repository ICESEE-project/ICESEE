# Idealized PIG: Mode-3 run on PACE

The `features_debug` branch includes the scientific updates, upstream
`main`, and the runtime repairs used in the successful five-year Mode-3 check.
Modes 1 and 2 also completed the five-year confirmation after the history and
finalization repairs. See [RUNTIME_DECISION.md](RUNTIME_DECISION.md) for evidence
and remaining scaling limits.

## Update your scientific checkout

From an existing ICESEE-project/ICESEE checkout, first inspect local changes:

```bash
git status --short
git remote -v
git fetch origin
git switch features_debug
git merge --ff-only origin/features_debug
```

Here `origin` must be `ICESEE-project/ICESEE`, not `ICESEE-Spack`. Keep your
uncommitted work before updating; if Git refuses a fast-forward, inspect the
local commits and merge them deliberately. Do not reset your work. A checkout may name its scientific remote `icesee`; verify the URL.
Updating the scientific checkout does not update the parent installation or
its pinned source automatically. Select the checkout explicitly below.

## Submit the validated configuration

Submit from the example directory in the updated scientific checkout. In the
standard source-only layout, the launcher locates the source from the submission
directory and reuses the sibling dependency installation:

```bash
cd "$HOME/r-arobel3-0/ICESEE-Spack-source/ICESEE/applications/icepack_model/examples/idealized_pig"
sbatch run_idealized_1.sbatch
```

Explicit `ICESEE_SOURCE_ROOT` and `ICESEE_SPACK_ROOT` exports override discovery.
The parameter file comes from `ICESEE_PARAMS_FILE` when set, otherwise the
adjacent `ICESEE-run-configs/idealized_pig.yaml` when present, then the tracked
`params.yaml`. The selected paths and a parameter copy are recorded with logs.

The existing PIG mesh, parameter, SMB and compact spin-up input files must be
available in this checkout's `data/`. They are not downloaded by the launcher.
The compact spin-up input is approximately 29 MB; using the historical full
spin-up history is a separate input choice. The submission uses your configured
PACE account; override `--account` or `--partition` at submission if required by
your allocation.

Defaults are one node, eight MPI ranks, one CPU/thread per rank, 48 GB and a
two-hour limit. Two spatial ranks per model group provide four model groups;
24 members run in six rounds.
For this native Icepack path the total is `model_groups * model_nprocs`;
there is no additional coordinator rank per member. Nodes/ranks are resources,
not the ensemble size. The script checks divisibility and uses block placement
to keep each two-rank model group within a node.

The default allocation matches job 15369255. Both measured allocations use five
years, `dt=0.05`, analyses
at years 2, 3 and 4, and flowline profiles at years 0, 2, 3 and 4. No initial
thickness perturbation, solver-cap increase or physical-parameter adjustment
is introduced. Native process noise and the existing analysis controls are
applied. `--execution_mode=3` selects native Mode 3 directly.

## Mode selection when launching directly

`--execution_mode=3` selects the native PIG runner directly. No additional
opt-in is required. The same topology, adapter, layout, checkpoint and
observation-error validation still applies; invalid configurations fail rather
than falling back to another mode.

From the example directory, the 24-member/eight-rank configuration is:

```bash
DATA_DIR="$PWD/_modelrun_datasets_mode3_$(date +%Y%m%d_%H%M%S)"
if [[ -e "$DATA_DIR" ]]; then
  echo "Preserve existing output or choose another data directory: $DATA_DIR" >&2
else
  mpirun -np 8 python run_da_icepack.py \
    --Nens=24 --execution_mode=3 --model_nprocs=2 \
    --num_years=5 --dt=0.05 \
    --obs_start_time=2 --obs_max_time=4 --freq_obs=1 \
    --data_path="$DATA_DIR"
fi
```

The 24-member/eight-rank case completed on PACE in job 15369255 using only
`--execution_mode=3`: one node, two spatial ranks per model, 100 forecast steps
and three analyses. The job exited zero in 47:09. All eight final checkpoint
shards and all 25 truth/member profile files were checked. Use
`run_idealized_1.sbatch` or `pace/validate_mode3_24.sbatch` for that allocation. This confirms
the PACE run; each laptop's installed Firedrake/MPI stack still needs its own
compatibility check.

## Outputs and repeat runs

By default results go to `_modelrun_datasets/` in the example directory.
The script refuses an existing target, including a symlink, and holds a launch
lock. Keep prior results, or choose a new absolute output path:

```bash
export ICESEE_DATA_PATH="$PWD/_modelrun_datasets_mode3_repeat1"
sbatch run_idealized_1.sbatch
```

Do not reuse an active run's output target. This launcher is for fresh runs;
it deliberately refuses existing restart data. Use the documented checkpoint
resume interface and the same scientific configuration for a restart.

Logs, exact launch arguments, source revision/dirty status and an exit-status
file go to `pace_runs/mode3-JOBID/`. Slurm writes `pig-mode3-JOBID.out` in the
submission directory. Check both exit status and output contents:

```bash
squeue -u "$USER"
sacct -j JOBID --format=JobID,State,ExitCode,Elapsed
```

The dataset directory contains `hs_profiles_true.h5` and
`hs_ensemble_profiles0.h5` through `hs_ensemble_profiles23.h5` for the default
24-member run, plus native checkpoint shards/manifests. These are binary HDF5
files, not folders. Read profiles with h5py. The updated notebook reads the
standalone `icesee_ensemble_data.h5` published after full-history finalization;
older rolling archives need a native checkpoint reader and still lack pruned times.
Native profiles record `written_steps` and `profile_times`; unwritten columns
must not be interpreted as physical zero values. Compilation caches are placed
in node-local scratch, with rank-private XDG state, rather than shared storage.

### Native snapshots and the combined plotting file

The PIG launcher retains the initial ensemble and every completed cycle using
`checkpoint_every=1`, `checkpoint_keep_last=0`. Each rank writes only its owned
rows during the run. On completion, bounded finalization creates standalone
`icesee_ensemble_data.h5` containing `ensemble` (state/member/time),
`ensemble_mean` (state/time), physical `t`, and `analysis_applied`.
Reference-row permutations align exported state with serial truth/observations.
The initial reference mapping and all checkpoint shards remain available for
restart; copying the aggregate alone is sufficient for its plotting datasets.

For five years from time zero with `dt=0.05`, `checkpoint_00000099` is one final
snapshot at year 5, not all timesteps. That example's aggregate has 101 columns:
initial state plus 100 completed cycles numbered 0–99. Column counts follow the
chosen duration and `dt`: `K=int(round(T/dt))` cycles and `K+1` saved times. Use
`ICESEE_DT=0.1 sbatch run_idealized_1.sbatch` for 51 columns over five years.
Actual physical labels come from the run's model times, including when `T/dt`
is not an integer. At an analysis cycle its column
contains the updated ensemble. Flowline diagnostics remain separate and are
not full analyzed-state history.

Earlier rolling-retention runs may retain only initial state and steps 98/99.
Those missing times cannot be reconstructed or exported as complete history.
The new exporter fails explicitly rather than filling the aggregate with zeros
or silently dropping times. See the [README](README.md#8-locate-and-inspect-results)
for manifest inspection, disk-size estimates and finalization recovery.

## Changing allocation for a benchmark

Hold the ensemble and science fixed for strong scaling. For example, four
nodes with eight ranks each keeps two spatial ranks per model and runs two
member rounds:

```bash
export ICESEE_TASKS_PER_NODE=8
export ICESEE_DATA_PATH="$PWD/_modelrun_datasets_mode3_N4"
sbatch --nodes=4 --ntasks=32 --ntasks-per-node=8 run_idealized_1.sbatch
```

For weak scaling, declare which work grows. With fixed mesh, two spatial ranks
per model, and four members per group, doubling the validated allocation to
four nodes would use 64 members:

```bash
export ICESEE_NENS=64 ICESEE_TASKS_PER_NODE=8
export ICESEE_DATA_PATH="$PWD/_modelrun_datasets_mode3_weak_N4"
sbatch --nodes=4 --ntasks=32 --ntasks-per-node=8 run_idealized_1.sbatch
```

These are proposed benchmark commands, not measured results. Changing ensemble
size changes sampling and is not a scientific-equivalence comparison. Other
allocations/topologies need validation. Do not infer speedup from the single
successful run, which shared an allocation with Modes 1 and 2.

## Evidence and limits

Job 15125921.2 exited zero in 34:24 with 32 members, 100 forecast steps and
three analyses. Final checkpoint 99 recorded next step 100 and model time 5.0.
All 16 shards, five fields for every member, and all 33 profile files were
checked; final thickness was 30..4181.29 m. Profiles had 415 points by four
columns, times `[0,2,3,4]`, and written steps `[0,40,60,80]`. Earlier native
restart across an analysis reproduced uninterrupted fields to maximum absolute
difference `4.09e-12` in the two-member check.

These checks establish runtime completion and inspected output integrity, not
independent scientific accuracy or three-mode equivalence. Modes 1/2 use a
different generic initial-perturbation path. The native noise bridge still has
a bounded, leader-local global-state workspace; not every setup is fully
spatially scalable. See the runtime notes for the outstanding work.


## Exporting a failed run's available history

The common batch launcher attempts serial finalization after a failed MPI step,
using only consecutive committed checkpoints. Its `exit-status` remains the
simulation's original failure code; `finalizer-exit-status` and `finalizer.log`
report recovery separately. For a manually launched or forcibly interrupted run,
stop all its processes and run the standalone exporter with `--allow-partial`.
The README gives the exact command and the partial-file schema.

A partial file has `export_complete=True`, `complete=False`, explicit
requested/completed counts and a `t` axis ending at the last committed time.
Initial-only recovery is supported. Missing middle snapshots, missing shards,
invalid numbering and nonfinite values are errors. A complete published file
is protected from partial replacement; `--output` can select another target.
The notebook and generic output reader can use the partial file without the
end-of-run summary. Hard node/allocation termination may prevent automatic
recovery, so manual finalization remains necessary in that case.

The historical full-history run exposed an extra initial forecast and a
one-cycle observation offset. Its initial thickness matched truth at `dt=.05`
instead of zero, and observations at 2/3/4 years updated states at
2.05/3.05/4.05. Those saved values and timestamps remain unchanged.

The corrected clock copies initial physical fields without advancing, analyzes
observation column `j` after cycle `j-1`, and applies configured process noise
before that analysis. Restart metadata and native noise identities distinguish
this contract from older output. Choose a fresh data path when upgrading;
initial-time observations require an unimplemented pre-forecast analysis stage
and are rejected. See the README's time-grid section for details.

PACE job **15432397** completed the one-year, two-member tests in Modes 1–3,
a four-rank noise ownership/replay check, and a native restart. All steps exited
zero. Analyses occurred at 0.5, 0.75 and 1.0 years, including the final state.
Restart from checkpoint 13 through 19 crossed two analysis events; its final
relative L2 difference from the uninterrupted run was 6.57e-13 (maximum absolute
difference 2.99e-9), within the declared `rtol=1e-9`, `atol=1e-7` checks.
This replay used copied completed output, not an interrupted partial-profile test.

The three five-year, 32-member simulations in job **15433074** also exited zero:

| Execution mode | MPI ranks / ranks per model | Simulation step | Elapsed |
| --- | --- | --- | --- |
| 1 | 8 / 1 | 15433074.2 | 42:15 |
| 2 | 8 / 1 | 15433074.0 | 1:38:31 |
| 3 | 16 / 2 | 15433074.1 | 37:11 |

These concurrent runs shared a two-node allocation; these times are operational
measurements, not a controlled scaling comparison. Each output has ensemble
shape `(270615, 32, 101)` at `dt=0.05`, three analyses at years 2/3/4
(columns 40/60/80), and 33 finite truth/member profile files of shape `(415, 4)`
with written times 0/2/3/4. Native output retains the initial condition and all
100 committed forecast steps; checkpoint 99 records next step 100 and time 5.

The original allocation **timed out during post-run verification**, after all
three simulations had finished successfully and the Mode 1–2 output scans had
passed. Verification-only job **15439331** also timed out in repeated virtual
HDF5 reads. These allocation failures are retained in the record. A subsequent
read-only row-slab reader reproduced every short-case check and completed the
remaining full native field/mean/profile/manifest checks. Together those checks
took 67.37 seconds wall time; this is reader validation, not model runtime.
No simulation was rerun to obtain that result.

All final Mode 1–2 members were compared in bounded slabs: the largest per-field
relative L2 difference was 9.12e-13. Native initial thickness and surface exactly
matched truth at zero; initial velocity differences were at most 9.10e-13.
The native five-year thickness range was 29.94–4180.72 m, with final values
30.00–4180.68 m; velocity components stayed within approximately -7013 to 1550
in model units. Assimilation produced transient basal-melt values -0.3664 to
20.0441, returning to 0–20 at the final state; no unreviewed clipping was added.
These checks establish timing, output integrity and the stated restart behavior,
not scientific accuracy. Generic and native initial perturbation paths differ,
so the three modes are not claimed to have identical member trajectories.

To reproduce the prepared cases from this example directory (requires a PACE
allocation and existing lab inputs):

```bash
export ICESEE_SOURCE_ROOT="$(git rev-parse --show-toplevel)"
export ICESEE_SPACK_ROOT="$HOME/r-arobel3-0/ICESEE-Spack"
sbatch pace/test_modes_1_3_year.sbatch
# After checking the short case:
sbatch pace/validate_production_modes.sbatch
```

Both scripts run `tools/verify_timing_outputs.py` on their closed outputs. This
checker targets these fixed observation/profile schedules and geometry; it is
not a general validator for arbitrary PIG configurations.


Job 15384953 (`da6f9dde`, 24 members/eight ranks) completed the full five-year
history path in 1:29:27 with every step exiting zero. The aggregate contains
270615 state rows, 24 members and 101 time columns. Independent reconstruction
compared every one of its 655,970,760 values and every mean exactly with the
checkpoint shards; all 25 profiles passed shape, written-time and finite-value
checks. Export took 2485.64 seconds in that baseline. The revised contiguous
reader/recovery implementation is verified separately on the closed archive;
a new PACE performance result is not implied. See the README for the measured
scope of the historical storage test and the separate corrected-clock validation.
