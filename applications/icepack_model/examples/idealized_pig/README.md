# Idealized Pine Island Glacier: installation, runs and results

This example couples Icepack/Firedrake forecasts to ICESEE ensemble data
assimilation. It includes truth generation, synthetic observations, an optional
unassimilated reference trajectory, ensemble forecasts, assimilation, flowline
diagnostics and restart checkpoints.

The runtime and output changes described here are on `features_debug`. The
completed full-history PACE case used `da6f9dde` on 9 November 2026; the recovery
extension is tested separately below. Record the exact source revision for each
run. Do not assume an older `main`, installation directory or notebook uses the
same code. Other adapters have their own mode and dependency requirements.

Runtime completion and exact output reconstruction do not establish scientific
time alignment. The current initializer advances one extra model step, and
analysis scheduling has a one-cycle observation offset. See
[time alignment under review](#time-alignment-under-review) before interpreting
assimilation plots.

## 1. Choose the source and retain your working environment

You need a working Icepack, Firedrake, PETSc, MPI and HDF5 environment, the ICESEE
source and the PIG inputs. Updating Python source does not require rebuilding
Spack or reinstalling Firedrake. Use the existing environment when it passes the
checks below. For a new environment, follow the current [Firedrake installation
guide](https://www.firedrakeproject.org/install.html) and
[Icepack installation guide](https://icepack.github.io/install/).

### Existing standalone scientific repository

Run inside the ICESEE repository, where `src/` and `applications/` exist:

```bash
git status --short
git remote -v
git fetch origin features_debug
git switch features_debug
git merge --ff-only origin/features_debug
git rev-parse HEAD
```

Here `origin` must point to `ICESEE-project/ICESEE`. If it points to
`ICESEE-Spack`, you are operating in the parent repository. Stop and locate the
scientific repository before fetching or publishing scientific code. If there
are local edits or a fast-forward fails, preserve and review them. Do not use a
hard reset to make an update succeed.

### ICESEE-Spack with a real submodule

The parent and its ICESEE child are separate Git repositories:

```bash
cd /path/to/ICESEE-Spack
git ls-tree HEAD ICESEE
git submodule status ICESEE
git -C ICESEE status --short
git -C ICESEE remote -v
```

`160000 commit` identifies a submodule; `040000 tree` identifies an ordinary
tracked directory, including historical subtree installations. A `.gitmodules`
file alone does not prove that ICESEE is a submodule.

For a clean, current submodule installation, initialize the pinned child:

```bash
git submodule update --init ICESEE
```

To deliberately use the development updates described here, inspect the child
and then select its branch:

```bash
git -C ICESEE fetch origin features_debug
git -C ICESEE switch features_debug
git -C ICESEE merge --ff-only origin/features_debug
```

The child may now differ from the parent pin. This is expected until a reviewed
parent commit records the new child revision. Running `git submodule update`
again can return it to the older pin. See the wiki's
[Source updates and submodules](https://github.com/ICESEE-project/ICESEE/wiki/9.-Source-updates-and-submodules).

### Keep an older subtree installation intact

Do not run a subtree pull into a new submodule or convert a dirty installation
in place. A source-only checkout can reuse the original dependencies:

```bash
cd /path/with/space
git clone --no-recurse-submodules https://github.com/ICESEE-project/ICESEE-Spack.git ICESEE-Spack-source
cd ICESEE-Spack-source
git submodule update --init ICESEE
git -C ICESEE switch features_debug
git -C ICESEE pull --ff-only origin features_debug
```

Inspect remotes and clean status first. Initialize only the ICESEE child when
reusing installed dependencies; there is no reason to download and build a
second Spack stack for a source update. Keep the original installation and
inputs until the new run is verified.

## 2. Activate dependencies and verify imports

On PACE, set both roots. The source root contains this example; the dependency
root contains `scripts/activate.sh` and `.spack-env/icesee`:

```bash
export ICESEE_SOURCE_ROOT=/absolute/path/to/ICESEE-Spack-source/ICESEE
export ICESEE_SPACK_ROOT=/absolute/path/to/existing/ICESEE-Spack
source "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig/pace/resolve_environment.sh"
```

The helper changes to the PIG example directory. Verify the active interpreter,
package location and parallel HDF5 support before a scheduler submission:

```bash
python - <<'PYTHON'
import sys, ICESEE, firedrake, icepack, h5py
from mpi4py import MPI
print("Python:", sys.executable)
print("ICESEE paths:", list(ICESEE.__path__))
print("Icepack:", icepack.__file__)
print("Firedrake:", firedrake.__file__)
print("Parallel h5py:", h5py.get_config().mpi)
print("MPI:", MPI.Get_library_version())
PYTHON
git -C "$ICESEE_SOURCE_ROOT" rev-parse HEAD
git -C "$ICESEE_SOURCE_ROOT" status --short
```

Import checks do not validate a multi-rank run. Keep the launcher, mpi4py,
PETSc and h5py in the same MPI environment; do not mix a system `mpirun` with
Python packages linked to another MPI library. Avoid stacking unrelated Conda,
virtualenv and Spack environments.

For a laptop, activate its working Firedrake/Icepack environment and make the
ICESEE checkout importable using that environment's installation procedure.
The PACE resolver is specific to ICESEE-Spack; it is not a laptop activation
script. PACE validation does not certify a different laptop software stack.

## 3. Provide and identify inputs

Use the lab's existing approved PIG inputs from the retained ICESEE-Spack
installation. The new scientific checkout can link to those files while using
the same installed dependencies. Users joining the project need access to that
lab input directory before launching this example. The current configuration uses:

| Input | Purpose |
|---|---|
| `data/extended_beta1000yrs_compact_idx20000.h5` | Compact final spin-up checkpoint, approximately 29 MB |
| `data/PigFull2017GeomFull.exp` | Argus geometry used by mesh setup |
| `data/OLS_Trend_plus_Resid_9b9.tif` | Surface mass-balance raster |
| `data/PigFull2017GeomFull.yaml`, if present | Optional mesh settings read beside the `.exp` file |

The historical `extended_beta1000yrs.h5` is about 39.5 GB and contains the full
spin-up history. The compact checkpoint avoids reading that history for each
initialization. `tools/build_compact_initialization.py` extracts the required
state. Use its `--help` and verify the extracted fields against the source;
a different compact file is not automatically the same initial condition.
The legacy `paramsFile` YAML key is not consumed by the active PIG initialization
path at this revision. Do not create an empty input file to silence a path error.

When the approved inputs already exist in another installation, links avoid a
second copy. From this example directory, for example:

```bash
PIG_INPUTS="$ICESEE_SPACK_ROOT/ICESEE/applications/icepack_model/examples/idealized_pig/data"
if [[ ! -d "$PIG_INPUTS" ]]; then
  PIG_INPUTS="$ICESEE_SPACK_ROOT/applications/icepack_model/examples/idealized_pig/data"
fi
printf 'Existing lab input directory: %s\n' "$PIG_INPUTS"
mkdir -p data
for name in extended_beta1000yrs_compact_idx20000.h5 PigFull2017GeomFull.exp OLS_Trend_plus_Resid_9b9.tif; do
  if [[ ! -f "$PIG_INPUTS/$name" ]]; then
    echo "Missing approved input: $PIG_INPUTS/$name" >&2
  elif [[ -e "data/$name" || -L "data/$name" ]]; then
    echo "Keep and inspect existing input: data/$name"
  else
    ln -s "$PIG_INPUTS/$name" "data/$name"
  fi
done
```

The two path forms support a nested ICESEE installation and the historical
standalone layout. If your lab stores inputs elsewhere, set `PIG_INPUTS` to that
approved directory. Inspect any reported missing files or existing targets before
launching; the commands preserve existing destinations. If the original mesh
has an optional adjacent `.yaml`, preserve that input too. For an intentionally
different approved dataset, use a separate directory/configuration and record
its identity.

Inputs can remain in their original storage location. Use explicit configuration
paths or links to approved existing inputs instead of duplicating the full
history. Record hashes of files actually used. On Linux:

```bash
sha256sum data/extended_beta1000yrs_compact_idx20000.h5 \
  data/PigFull2017GeomFull.exp data/OLS_Trend_plus_Resid_9b9.tif
```

Keep input filenames, hashes, source revision and environment information with
run results. Identical names do not establish identical contents. Different HDF5 file hashes can also reflect container metadata rather than changed numerical fields; inspect dataset contents and writer metadata before concluding that the physical initial state differs.

## 4. Keep personal scientific settings outside tracked source

Copy the reviewed YAML once, then edit your own copy:

```bash
mkdir -p "$HOME/ICESEE-run-configs"
cp params.yaml "$HOME/ICESEE-run-configs/idealized_pig.yaml"
```

Use `-F` or `--force-params` to select it. Relative input paths are resolved from
the example's working directory; absolute input paths avoid ambiguity across
checkouts. Explicit command-line overrides select the run's member count,
mode, duration, timestep and output directory. Keep those arguments in the
recorded command as well as retaining the YAML.

Important experiment settings include:

- `generate_true_state` and `generate_synthetic_obs` generate the reference and observations.
- `generate_nurged_state=1` also creates the full wrong/no-assimilation trajectory. It is optional and adds work; the validated reference run used `0`.
- `initialize_ensemble` controls fresh ensemble initialization. Resume restores saved members instead.
- `save_steps` specifies requested profile years; the initial profile is saved separately.
- `seed`, `sig_Q`, `sig_obs`, inflation and localization are scientific choices.
- Initial native PIG thickness is unperturbed. Basal-melt forcing is the intended experiment discrepancy. Do not add thickness noise as a performance fix.

Do not change solver caps, timesteps, forcing, initial uncertainty or analysis
options to make a scaling comparison complete. Such changes define a different
scientific experiment and require separate validation.

## 5. Select an execution mode and resource layout

| Mode | PIG execution | Important distinction |
|---|---|---|
| 0 | Serial runner | Useful as a small reference where supported |
| 1 | Partially parallel forecast/analysis workflow | Generic state/history path |
| 2 | Parallel forecast with distributed I/O/analysis workflow | Generic state/history path |
| 3 | Native distributed Icepack state and ensemble-by-space execution | Requires the registered native adapter and valid ownership/topology |

Select the runner with **`--execution_mode=3`** (or 1/2). There is no additional
Mode-3 opt-in flag. Invalid native layouts or unsupported options fail explicitly;
the code does not silently fall back to Mode 2. The legacy positional scheduling
argument is a separate setting; omit it in the commands below.

For native Mode 3:

```text
MPI ranks = concurrent model groups × spatial ranks per model group
model groups = MPI ranks / model_nprocs
```

Ensemble members are scheduled over those groups. Eight ranks with
`model_nprocs=2` provide four groups; 24 members require six rounds. There is
no additional MPI coordinator rank per member in this native Icepack path.
External MATLAB/ISSM processes have separate resource requirements.

On multiple nodes, choose uniform ranks per node divisible by `model_nprocs`
and use block placement so a model group stays within a node. The launcher
checks allocation sizes, and rank logs record host placement for inspection. `Nens`
can exceed the number of groups; it is not a request for the same number of MPI
ranks. Strong scaling holds `Nens`, mesh and scientific settings fixed.

## 6. Submit a fresh PACE run

Submit directly from `applications/icepack_model/examples/idealized_pig` in the
checkout you want to execute. The launcher discovers that checkout from
`SLURM_SUBMIT_DIR`. In the source-only layout it reuses the sibling
`ICESEE-Spack` installation; explicit `ICESEE_SOURCE_ROOT` and
`ICESEE_SPACK_ROOT` exports take precedence.

The default allocation is 24 members, one node, eight MPI ranks, two spatial
ranks per group, 48 GB and a two-hour cap:

```bash
cd "$HOME/r-arobel3-0/ICESEE-Spack-source/ICESEE/applications/icepack_model/examples/idealized_pig"
sbatch run_idealized_1.sbatch
```

The default result directory is this example's `_modelrun_datasets`. It must
not already exist, including as a symlink. Keep completed results and select a
new path for each additional run:

```bash
export ICESEE_DATA_PATH="$PWD/_modelrun_datasets_trial02"
sbatch run_idealized_1.sbatch
```

The parameter file is selected in this order: explicit `ICESEE_PARAMS_FILE`,
`$HOME/r-arobel3-0/ICESEE-run-configs/idealized_pig.yaml` in the standard lab
layout (beside the two installation directories), then this example's tracked
`params.yaml`. The launcher prints the selected paths and saves a copy of the
parameters in `pace_runs/mode3-JOBID/params.yaml`. To select another configuration:

```bash
export ICESEE_PARAMS_FILE="$HOME/r-arobel3-0/ICESEE-run-configs/idealized_pig.yaml"
sbatch run_idealized_1.sbatch
```

The launcher passes that file through `--force-params`; editing the tracked
script or copying it elsewhere is unnecessary. Arguments after the script name
are not forwarded to Python. Defaults are five years, `dt=0.05`, observations at
years 2/3/4 and requested profiles at years 0/2/3/4. Select another timestep at
submission with `ICESEE_DT`; `ICESEE_NUM_YEARS` selects duration. A changed duration
still needs a valid observation window.

```bash
export ICESEE_DT=0.1
export ICESEE_DATA_PATH="$PWD/_modelrun_datasets_dt01_trial01"
sbatch run_idealized_1.sbatch
```

The launch record prints resolved forecast steps and history columns. Explicit
launcher arguments override the YAML timestep/duration, so inspect `launch.txt`.

For the two-node, 32-member allocation:

```bash
export ICESEE_NENS=32
export ICESEE_DATA_PATH="$PWD/_modelrun_datasets_32members_trial01"
sbatch --nodes=2 --ntasks=16 --ntasks-per-node=8 --mem=96G run_idealized_1.sbatch
```

Check account/partition access before submitting. Email notifications can be
requested with `sbatch --mail-user=you@gatech.edu --mail-type=BEGIN,END,FAIL
run_idealized_1.sbatch`; keep personal addresses outside tracked scripts.

For a standalone MPI run on resources where you are allowed to run it:

```bash
DATA_DIR="$PWD/_modelrun_datasets_local_trial01"
if [[ -e "$DATA_DIR" || -L "$DATA_DIR" ]]; then
  echo "Choose a fresh output directory: $DATA_DIR" >&2
else
  mpirun -np 8 python run_da_icepack.py \
    -F "$HOME/ICESEE-run-configs/idealized_pig.yaml" \
    --execution_mode=3 --Nens=24 --model_nprocs=2 \
    --num_years=5 --dt=0.05 \
    --obs_start_time=2 --obs_max_time=4 --freq_obs=1 \
    '--save_steps=[2,3,4]' --data_path="$DATA_DIR"
fi
```

Use Slurm for computational runs on PACE; do not execute this MPI workload on
a login node. A direct fresh CLI run can recreate the configured output path;
keep the guard above and never point it at previous results.

For Mode 1 or 2, select that mode explicitly and use its supported rank layout.
`pace/validate_modes_1_2.sbatch` is a two-mode confirmation script, not a
controlled scaling benchmark. The generic and native initialization paths
differ; changing the mode does not by itself establish member-wise equivalence.

## 7. Cache placement and monitoring

`pace/launch_rank.sh` uses node-local `$TMPDIR`, falling back to `/tmp` for
recognized shared filesystems. It assigns job/step-local PyOP2 and TSFC cache
paths and rank-private `XDG_CACHE_HOME` for SQLite state. The launch log records
host/cache paths. A larger shared project directory can hold inputs, results
and checkpoints, but available storage alone does not make it suitable for
concurrent JIT/SQLite cache traffic. Keep thread counts at one per MPI rank
unless a separately validated hybrid layout is intended.

```bash
squeue -u "$USER"
sacct -j JOBID --format=JobID,State,ExitCode,Elapsed
cat pace_runs/mode3-JOBID/launch.txt
cat pace_runs/mode3-JOBID/source.txt
tail -n 40 pace_runs/mode3-JOBID/stdout.log
tail -n 40 pace_runs/mode3-JOBID/stderr.log
```

Truth and optional wrong-state generation currently run on world rank zero;
other ranks wait during that phase. Idle peers then do not prove Mode 3 was
ignored. Watchdog timeout stack dumps are diagnostics, not automatically a
failed solve. Check changing stacks, CPU/I/O progress and checkpoint commits
before declaring a hang. Preserve failed logs and data.

A nonlinear convergence error is different from an MPI or storage deadlock.
Record the failing truth/member/step, pre/post-analysis fields and solver
residuals. Do not silently increase the solver iteration cap or disable MPI
file locking. The history fixes retain normal write safety.

## 8. Locate and inspect results

The end-of-run `performance.json` records phase timings and available numerical-package versions. Retain it with the scheduler accounting and input provenance. Phase timers can overlap: `forecast_step/with_analysis` is a subset of forecast-cycle time, not an isolated additive analysis timer. Rank-summed time is not wall time, and logical checkpoint throughput is not a raw filesystem benchmark.

The common launcher writes logs, source state, exact arguments and
`exit-status` to `pace_runs/mode3-JOBID/`. Data go to `ICESEE_DATA_PATH`, or
`_modelrun_datasets/` by default. Custom launchers can choose another log root;
read their `launch.txt` rather than assuming a path.

Typical native output includes:

```text
DATA_DIR/
  hs_profiles_true.h5
  hs_ensemble_profiles0.h5
  ...
  hs_ensemble_profiles23.h5       # for a 24-member run
  synthetic_obs.h5
  true_nurged_states.h5
  performance.json
  icesee_ensemble_data.h5       # published after successful full-history finalization
  _mode3_state_history/
    reference_layout.h5        # serial-reference to native-row mapping
    initial_condition/checkpoint_00000000/
    steps/checkpoint_00000099/
      manifest.json
      shard_00000000.h5
      ...
```

`hs_profiles_true.h5` is a binary HDF5 file, **not a folder** or a text document.
The member filenames are `hs_ensemble_profiles{member}.h5`, not
`hs_profiles_true{member}.h5`. Older files may lack an extension. Inspect the
selected data directory and the code revision before concluding generation
failed. A text editor's binary-file warning is expected.

```bash
python - "$ICESEE_DATA_PATH" <<'PYTHON'
from pathlib import Path
import sys, h5py, numpy as np
root=Path(sys.argv[1])
for path in sorted(root.glob('hs*.h5')):
    with h5py.File(path, 'r') as f:
        print(path.name)
        for key in ('h_profiles', 's_profiles'):
            a=f[key][:]
            print(key, a.shape, 'finite=', np.isfinite(a).all())
        for key in ('written_steps', 'profile_times'):
            print(key, f[key][:] if key in f else f.attrs.get(key))
PYTHON
```

Check all members, requested/written columns, finite values and physical ranges.
An unwritten profile column is not evidence of physically zero thickness; use
metadata and the simulated interval. Older generic output can contain zero
placeholders. Native output records written columns explicitly. A finite array
can still be physically unreasonable.

Native member profiles are sampled after forecast advancement, before any analysis in that cycle. They are flowline forecast diagnostics, not a complete analyzed-state history. Use the committed checkpoint to inspect analyzed fields.

### Distributed checkpoints during the run, one plotting file at completion

The PIG launcher and example YAML retain every native timestep. During the
run, each rank writes its owned state rows to checkpoint shards. This avoids
first gathering the entire ensemble onto one rank. Each committed manifest
describes member/field ownership; keep it with all referenced shards for restart.

After all cycles complete, the runner materializes **`icesee_ensemble_data.h5`**
from those closed checkpoints. It reads bounded row chunks and one member at a
time, computes the ensemble mean, and atomically publishes the finished file.
A failed finalization preserves checkpoints and any previously published file;
the run reports failure instead of publishing an incomplete aggregate. Its
`analysis_file_io/ensemble_history_export` timer is a subset of total
`analysis_file_io`, allowing benchmarks to separate finalization from forecasts.

The combined file is standalone HDF5, with no external links to the shards:

| Dataset | Shape | Content |
|---|---|---|
| `ensemble` | `(state, member, time)` | Initial ensemble and every completed cycle |
| `ensemble_mean` | `(state, time)` | Member mean at each saved time |
| `t` | `(time,)` | Physical model years |
| `checkpoint_steps` | `(time,)` | `-1` for the initial state, then zero-based cycle indices |
| `analysis_applied` | `(time,)` | Whether that cycle included an analysis update |

Rows are reordered into the **serial reference numbering used by truth and
observations**. Native MPI row numbers alone are not guaranteed to identify the
same physical node. The coordinate-derived permutations are saved once in
`_mode3_state_history/reference_layout.h5`; finalization refuses incompatible
numbering, incomplete member/field ownership or missing snapshots. Initial
thickness, forcing, random streams and filtering/solver settings are not altered
by this storage path.

The aggregate includes the updated ensemble on analysis cycles. It does not
add a separately saved pre-analysis ensemble at the same time. Flowline files
remain forecast diagnostics sampled at their requested years.

Older completed runs made with rolling retention may have only the initial
checkpoint and the newest two steps, and no aggregate file. Changing readers
cannot recreate their pruned snapshots. The full-history behavior applies to
new runs made with the updated launcher/configuration.

### Time columns follow the chosen timestep

The exporter reads the resolved `nt` and physical `model_time` values from
checkpoint manifests. It does not hard-code 100 cycles or 101 columns. For this
PIG driver, requested duration `T` and user-selected `dt` resolve to
`K = int(round(T / dt))` cycles, `K + 1` columns including the initial state,
and `t[j] = j * dt`. Final checkpoint index is `K - 1`.

| Five-year example | Forecast cycles | Time columns |
|---|---|---|
| `dt=0.1` | 50 | 51 |
| `dt=0.05` | 100 | 101 |
| `dt=0.025` | 200 | 201 |

If `T / dt` is not an integer, the driver's rounded cycle count means the
actual final time is `K * dt`. The file preserves that actual time grid; it does
not stretch the labels to end at the requested duration. Check the launch
schedule when choosing such values. `ICESEE_DT=0.1 sbatch run_idealized_1.sbatch`
selects the first example without editing the script.

### Does `checkpoint_00000099` contain every timestep?

**No. Each checkpoint contains one ensemble snapshot.** Its shards partition
that snapshot across ranks; they are not different times or necessarily one
file per member. The final aggregate is what combines snapshots over time.

For a run starting at year zero, five years with `dt=0.05` means 100 cycles:

| Location under `_mode3_state_history/` | Snapshot time | Aggregate time column |
|---|---|---|
| `initial_condition/checkpoint_00000000/` | 0 years, before forecasting | 0 |
| `steps/checkpoint_00000000/` | 0.05 years, after cycle 0 | 1 |
| `steps/checkpoint_00000098/` | 4.95 years | 99 |
| `steps/checkpoint_00000099/` | 5.0 years, final state | 100 |

The initial and step checkpoint zero have different meanings and live in
separate directories. With all cycles retained, the combined file has **101
columns: the initial state plus 100 completed steps**. The final manifest has
`timestep=99`, `metadata.next_step=100`, `metadata.nt=100` and
`metadata.model_time=5.0`. `analysis_events_completed=3` is a cumulative count;
it does not mean that three earlier states are inside checkpoint 99.

Consult each manifest for other run lengths. To list retained snapshots without
opening their state shards:

```bash
python - "$ICESEE_DATA_PATH" <<'PYTHON'
from pathlib import Path
import json, sys
root = Path(sys.argv[1]) / '_mode3_state_history' / 'steps'
for path in sorted(root.glob('checkpoint_*/manifest.json')):
    manifest = json.loads(path.read_text())
    metadata = manifest.get('metadata', {})
    print(path.parent.name, 'complete=', manifest.get('complete'),
          'time=', metadata.get('model_time'),
          'next_step=', metadata.get('next_step'), 'nt=', metadata.get('nt'))
PYTHON
```

### Full-history settings and recovery of finalization

The PIG launcher passes `--checkpoint_every=1 --checkpoint_keep_last=0
--mode3_export_ensemble_history=True`; the example YAML has the same defaults.
No extra Mode-3 opt-in is needed. Other native applications keep their own
checkpoint policies. Enabling this exporter requires both writing every cycle
and keeping every checkpoint; incompatible retention is rejected before model
setup. `mode3_export_row_chunk_size` bounds the number of reference rows in a
finalization buffer (default 8192).

Keep the checkpoint tree for restart. The standalone aggregate can be copied or
downloaded by itself for plotting. Retaining both duplicates state payload on
disk: for 270615 state rows, 24 members, float64 and 101 times, raw checkpoint
state is about 5.25 GB, the aggregate ensemble about 5.25 GB and its mean about
0.22 GB, before HDF5/metadata overhead. This is a size calculation, not measured
throughput or a new scaling benchmark. The earlier measured PACE cases below
predate full-history finalization; they do not establish its production-scale
I/O time. Export tests cover reordered/split fields, complete time coverage,
atomic failure, standalone reads and collective error handling.

If all checkpoints and `reference_layout.h5` exist but finalization was
interrupted, rerun only the export in an environment permitted for that I/O:

```bash
python -m ICESEE.src.parallelization.distributed_history_export   /absolute/path/to/completed-data --row-chunk-size=8192
```

The default command requires the complete requested history. For a failed or
interrupted simulation, stop its processes first, then export the committed
prefix without rerunning the model:

```bash
python "$ICESEE_SOURCE_ROOT/src/parallelization/distributed_history_export.py" \
  /absolute/path/to/failed-data --allow-partial
```

Run this in the existing Python environment with NumPy and h5py; it does not
initialize MPI or import Firedrake. The output contains the initial ensemble
and every consecutive committed cycle through the last available checkpoint.
An uncommitted trailing checkpoint is ignored. Missing middle snapshots,
incompatible row numbering, missing shards and nonfinite fields are errors;
recovery never fabricates missing time columns. Old pruned archives cannot be
reconstructed.

The common batch script attempts this recovery automatically **after `srun`
returns a nonzero status**. It records `finalizer.log` and
`finalizer-exit-status` beside the original logs. The job and `exit-status`
retain the original simulation failure even when recovery succeeds. A hard
node failure, forced termination or exhausted allocation may prevent that
recovery process from running; use the standalone command afterward in an
environment permitted for its I/O and memory use. For manual `mpirun` launches,
run the standalone command after the failed MPI processes have exited.

Partial output uses the same dataset names and serial-reference row order as
a complete history. Read `t` for its actual available time columns. Attributes
`history_status="partial"`, `complete=False`, `export_complete=True`,
`requested_steps`, `completed_steps`, `last_completed_step` and
`last_model_time` describe its extent. `simulation_exit_code` is recorded by
batch recovery when available. Here `export_complete` describes a finished
file export; it does not say that the simulation succeeded. If all cycles
were saved before a later error, `complete=True` can coexist with a nonzero
simulation exit code.

Partial recovery refuses to replace a complete, unrelated or longer published
history. Use `--output /absolute/path/to/recovered-history.h5` for a separate
file when necessary. A later successful full export can replace a partial
file atomically. The checkpoint tree and previous published file are preserved
if export fails. The common `sbatch` script does not forward arbitrary trailing
Python flags.

The notebook and `load_da_outputs` accept finalized partial histories and
restrict their time axes and observation columns to the available prefix.
End-of-run summary metadata is optional for native exports; the initial
checkpoint records the observation plan. A missing optional background
trajectory is reported and shown as unavailable, not replaced with zeros.

## 9. Resume without destroying output

Use the same scientific YAML, ensemble identities, inputs and supported topology.
A restart must restore native state and RNG progression instead of regenerating
the ensemble. Only a committed checkpoint with its manifest and all referenced
shards is a restart candidate. Keep the original truth/observations/profiles.

From an allocated compute environment, the direct native interface is:

```bash
mpirun -np 8 python run_da_icepack.py \
  -F /absolute/path/to/the-original-run.yaml \
  --execution_mode=3 --Nens=24 --model_nprocs=2 \
  --num_years=5 --dt=0.05 \
  --obs_start_time=2 --obs_max_time=4 --freq_obs=1 \
  '--save_steps=[2,3,4]' --data_path=/absolute/path/to/existing-data \
  --resume_from_checkpoint=True
```

For Slurm, create a dedicated resume launcher retaining the same environment,
resource binding and provenance, with an existing-data check appropriate for
resume. **Do not reuse the common fresh-run launcher unchanged:** it correctly
refuses an existing data target. Do not remove its guards for routine fresh runs.

The generic native policy can retain only the newest two step checkpoints;
the PIG plotting workflow instead sets `checkpoint_every=1` and
`checkpoint_keep_last=0` to retain the full history. Do not enable aggregate
export for an older rolling archive with missing times. These settings affect
storage and restart availability. Arbitrary mesh/topology migration is not established by
the short same-topology restart test.

## 10. Plot and compare runs

Use the data path for the run you actually launched. New native aggregate output
has the familiar `ensemble`/`ensemble_mean` datasets. Load only the requested
rows/times rather than copying the entire ensemble history into memory. For portable
flowline plots, h5py/NumPy/Matplotlib can read `h_profiles`, `s_profiles`,
`distances`, `bed_values`, `profile_times` and `written_steps` directly without
loading Firedrake or a full spin-up history.

For a complete ensemble time series at one serial-reference state row:

```python
from pathlib import Path
import h5py
import matplotlib.pyplot as plt

data = Path("/absolute/path/to/completed-data")
row = 0  # Choose the desired row within the documented state blocks.
with h5py.File(data / "icesee_ensemble_data.h5", "r") as f:
    if not f.attrs.get("export_complete", f.attrs.get("complete", False)):
        raise ValueError("Ensemble history export did not finish")
    if not f.attrs.get("complete", False):
        print("Partial history through", f.attrs["last_model_time"], "years")
    t = f["t"][:]
    members = f["ensemble"][row, :, :].T
    mean = f["ensemble_mean"][row, :]
    analysis = f["analysis_applied"][:]
plt.plot(t, members, color="gray", alpha=0.3)
plt.plot(t, mean, color="blue", label="Ensemble mean")
plt.scatter(t[analysis], mean[analysis], label="Analysis cycles")
plt.xlabel("Time (years)")
plt.legend()
```

The existing notebook's ensemble-loading cells can use those same dataset names.
Mesh-map cells must use the actual PIG mesh and compatible row ordering; a
rectangular template mesh is not evidence of a correct physical map.

The standalone comparison utility needs NumPy, h5py and Matplotlib, not
Firedrake. Run it on **completed** native outputs:

```bash
python tools/compare_mode3_runs.py \
  /absolute/path/to/reference-data /absolute/path/to/candidate-data \
  --output /absolute/path/to/comparison-figures
```

It reads bounded final-field slabs, checks ownership fingerprints and matching
member/time metadata, and writes `comparison.json`, `field_metrics.csv`,
`h_profiles_comparison.png` and `s_profiles_comparison.png`. Profile figures show
truth, ensemble mean ± one ensemble standard deviation, and difference panels.
A numbering or coordinate mismatch stops comparison rather than subtracting
unrelated rows. Passing the same data directory twice is a reader/plot smoke
check only; it is not an independent validation run.

For a comparison, first match source revision, YAML/CLI science, seed, inputs
by content, member IDs, MPI topology and profile coordinates/times. Compare
complete final members in bounded field slabs; report maximum absolute error,
RMSE and relative L2 per field. Compare physical flowline profiles and plot
ensemble means with ensemble spread and truth. Confirm ownership/coordinate
fingerprints before subtracting raw checkpoint rows from different environments.
If numbering differs, resolve the physical mapping; equal vector lengths are
not sufficient. Set acceptance tolerances for the scientific quantities under
review rather than declaring success from a similar-looking mean curve.

## 11. Scaling and validation limits

Strong scaling fixes mesh, members, duration, timestep, observations, filtering,
noise and storage policy while increasing allocated resources. Weak scaling
must state what grows; increasing members changes sampling and the scientific
experiment. Report cold/warm cache conditions, repeated timings, all-rank memory,
phase timings, logical I/O and failures. A laptop/cluster wall-time comparison
is not controlled scaling when inputs, CPU, environment or settings differ.

Proposed native allocation commands, to submit only within your approved
allocation and project budget:

```bash
# Strong scaling: retain 32 members and the common launcher's scientific settings.
export ICESEE_NENS=32 ICESEE_MODEL_NPROCS=2 ICESEE_TASKS_PER_NODE=8
export ICESEE_DATA_PATH="$PWD/_modelrun_datasets_strong_N4_trial01"
sbatch --nodes=4 --ntasks=32 --ntasks-per-node=8 run_idealized_1.sbatch

# Weak scaling: four members per group, increasing to 64 members on four nodes.
export ICESEE_NENS=64
export ICESEE_DATA_PATH="$PWD/_modelrun_datasets_weak_N4_trial01"
sbatch --nodes=4 --ntasks=32 --ntasks-per-node=8 run_idealized_1.sbatch
```

These are unmeasured benchmark proposals. Keep output targets distinct and record
resource/settings changes between trials. The weak case changes ensemble size;
it is not a member-wise equivalence test against 32 members.

Verified PACE cases at the published baseline:

| Case | Result | Scope |
|---|---|---|
| Mode 3, 32 members, 16 ranks, 5 years | Job 15125921.2, exit 0, 34:24 | All final fields and 33 profiles inspected |
| Mode 3, 24 members, 8 ranks, 5 years | Job 15369255, exit 0, 47:09 parent / 46:48 MPI step | All eight final shards and 25 profiles inspected |
| Mode 3, second PACE installation, 24 members, 8 ranks, 5 years | Job 15371824, exit 0, 46:14 parent / 45:54 MPI step | All eight final shards and 25 profiles inspected |
| Mode 1, 32 members, 8 ranks, 5 years | Job 15142488.1, exit 0, 42:01 | Final members, means and profiles inspected |
| Mode 2, 32 members, 8 ranks, 5 years | Job 15142488.0, exit 0, 1:43:05 | Final members, means and profiles inspected |

The second-installation 24-member case used identical runtime/configuration
content, input datasets, generated observations and observation-error arrays.
Truth profiles matched exactly. All final member fields and 25 profiles were
finite with expected metadata; final thickness ranged from 30 to 4181.58 m.
Across all final thickness members/rows, maximum absolute difference was
0.211592 m, RMSE 0.0151454 m and relative L2 `8.57429e-6`; maximum member
flowline thickness difference was 0.0135394 m. Thus these runs are close but
not bitwise equivalent. Profiles before the first analysis differed by at most
`1.59e-12` m; larger differences appeared afterward. Python patch versions and
Open MPI versions differed (3.11.14/3.11.15 and 5.0.9/5.0.10 respectively),
but those differences have not been isolated as the cause. No scientific
acceptance threshold is implied. Version strings alone do not identify every
compiler, BLAS or numerical build option.

Mode 1/2 final per-field member relative L2 differences were at most `1.54e-12`.
That is a measured comparison, not a universal acceptance threshold. Modes 1/2
and native Mode 3 use different initial perturbation paths; these checks do not
establish three-mode member-wise scientific equivalence. Native observation
setup still contains replicated reads, and native noise retains a capped
leader-local global-state workspace. Alternative filters and optional member
streaming need their own correctness and performance validation; do not switch
filters to hide an I/O bottleneck. See [RUNTIME_DECISION.md](RUNTIME_DECISION.md)
and [MODE3_PACE.md](MODE3_PACE.md).

## 12. Report a reproducible problem

Include source hash/dirty status, import paths and package versions, exact
command and selected YAML, input sizes/hashes, allocation/rank placement,
cache locations, scheduler and process exits, last committed checkpoint,
relevant log/traceback, profile metadata and field extrema. Remove credentials
and private input contents. Preserve the failing run; do not clear HDF5
consistency flags or overwrite outputs to make a reader proceed.


### Time alignment under review

A full-history validation exposed a pre-existing one-cycle observation offset
in the PIG native schedule: the forecast advances before the cycle index is
compared with the observation time-grid column. With `dt=0.05`, observation
columns at years 2/3/4 were applied to saved states at 2.05/3.05/4.05. The
exporter preserves those actual model times and records `analysis_applied`;
it does not relabel or change the trajectory. This output validation is not
evidence that observation timing is scientifically correct. Boundary cases
(initial, interior and final observations) require a separate reviewed
scheduling correction and validation across execution modes.


There is a separate initialization offset. Both the existing generic
`initialize_ensemble` and native `initialize_member` call `Icepack` for one
`dt` before the ensemble is saved at metadata time zero. In job 15384953,
initial ensemble thickness differed from truth at zero by up to 40.123212 m,
but matched truth at 0.05 years within 1.14e-12 m; velocity and surface matched
that first truth step within 3e-12 in their stored units. Members had identical
initial thickness, so this is an extra initial advance, not added thickness
noise. The aggregate preserves the checkpoint's declared times and values;
it does not correct either scheduling convention. Changing initialization and
observation alignment requires a separate scientific review and validation.

### Full-history output verification

PACE job **15384953**, source `da6f9dde`, completed with 24 members, eight MPI
ranks, two spatial ranks/model, five years and `dt=0.05`. Parent and all steps
exited zero (parent 1:29:27; MPI step 1:28:54). It retained initial state and all
100 cycles, ending at checkpoint 99/next step 100/model time 5.0 with three
analyses and eight final shards. The standalone ensemble has shape
`(270615,24,101)` and mean `(270615,101)`; its file size is 5,634,162,632 bytes.

An independent reader compared **all 655,970,760 ensemble values** with the
source shards in serial-reference order: every value matched exactly, as did
all sequential-member means. All 25 truth/member profiles had finite `h`/`s`
arrays of shape `(415,4)`, written steps `[0,40,60,80]` and profile times
`[0,2,3,4]`. Final thickness was 30..4181.5861 m, with basal melt 0..20 in model
units. These are output-integrity results; the time-alignment issues above
remain unresolved.

The original point-selection exporter took 2485.64 seconds in that PACE run;
truth generation took 371.96 seconds and the forecast-loop timer 2347.68
seconds. Phase timers are not an additive partition of the entire job. The
revised reader uses bounded contiguous HDF5 reads instead of long point lists.
A local Mac replay and recovery test use the same closed checkpoint archive;
their timings are not a PACE speedup measurement. No new PDE trajectory or
solver setting is introduced by replaying finalization.


Recovery verification includes 132 focused checks covering HDF5 export, batch
status propagation, output readers, runner wiring and checkpoint behavior. The
revised full reader replayed the closed production archive in 54.72 seconds on
a Mac, and an independent shard reader again matched every value exactly.
A read-only 12-cycle prefix fixture exported 13 time columns through recorded
time 0.6, with all 84,431,880 values matching the corresponding full-history
prefix. That fixture is not a newly failed model run; separate MPI-abort tests
exercise recovery after an actual process failure. Local replay timings do not
measure cold-cache performance or PACE speedup.
