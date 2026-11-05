# Idealized PIG: student Mode-3 run on PACE

The `features_debug` branch includes the student's scientific updates, upstream
`main`, and the runtime repairs used in the successful five-year Mode-3 check.
Modes 1 and 2 retain unresolved production I/O issues described in
[RUNTIME_DECISION.md](RUNTIME_DECISION.md). Use Mode 3 for this handoff.

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
local commits and merge them deliberately. Do not reset your work. The PACE
integration worktree uses `icesee` for the scientific remote instead.
Updating the scientific checkout does not update the parent installation or
its pinned source automatically. Select the checkout explicitly below.

## Submit the validated configuration

Run these commands from the updated ICESEE checkout containing `src/` and
`applications/`. The dependency installation remains separate:

```bash
export ICESEE_SOURCE_ROOT="$(pwd -P)"
export ICESEE_SPACK_ROOT="$HOME/r-arobel3-0/ICESEE-Spack"
cd "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig"
sbatch run_idealized_1.sbatch
```

The existing PIG mesh, parameter, SMB and compact spin-up input files must be
available in this checkout's `data/`. They are not downloaded by the launcher.
The compact spin-up input is approximately 29 MB; using the historical full
spin-up history is a separate input choice. The submission uses your configured
PACE account; override `--account` or `--partition` at submission if required by
your allocation.

Defaults are two nodes, eight MPI ranks per node, one CPU/thread per rank,
96 GB per node and a two-hour limit. There are 16 MPI ranks total, two spatial
ranks per model group, eight model groups, and 32 members in four rounds.
For this native Icepack path the total is `model_groups * model_nprocs`;
there is no additional coordinator rank per member. Nodes/ranks are resources,
not the ensemble size. The script checks divisibility and uses block placement
to keep each two-rank model group within a node.

The scientific arguments match job 15125921.2: five years, `dt=0.05`, analyses
at years 2, 3 and 4, and flowline profiles at years 0, 2, 3 and 4. No initial
thickness perturbation, solver-cap increase or physical-parameter adjustment
is introduced. Native process noise and the existing analysis controls are
applied. The experimental flag selects native Mode 3 explicitly.

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
`hs_ensemble_profiles0.h5` through `hs_ensemble_profiles31.h5`, plus native
checkpoint shards/manifests. These are binary HDF5 files, not folders. Open
with h5py or the result notebook after selecting this run's data directory.
Native profiles record `written_steps` and `profile_times`; unwritten columns
must not be interpreted as physical zero values. Compilation caches are placed
in node-local scratch, with rank-private XDG state, rather than shared storage.

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
