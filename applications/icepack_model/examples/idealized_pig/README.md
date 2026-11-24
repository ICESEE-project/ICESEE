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

Initialization and observation scheduling now use a shared state-column clock.
The initial state is saved without an extra forecast, and observations and
configured process noise are aligned with their actual model times. See
[time grid, observations and restart](#time-grid-observations-and-restart) for
validation, supported boundaries and updating older runs.

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
- The supplied YAML selects `process_noise_normalization: per_step`, which removes total-run-length normalization from forecast noise. Older YAML files that omit this key retain `legacy_run_length` for historical reproducibility. Use a fresh data directory when changing it; see [noise scaling and restart rules](../../../../docs/handbook/numerics.md#initialization-and-process-noise).
- Initial native PIG thickness is unperturbed. Basal-melt forcing is the intended experiment discrepancy. Do not add thickness noise as a performance fix.

Do not change solver caps, timesteps, forcing, initial uncertainty or analysis
options to make a scaling comparison complete. Such changes define a different
scientific experiment and require separate validation.

### Estimated basal-melt forcing

The default `basal_melt_forecast: prescribed_wrong` rebuilds the imposed wrong
forcing every forecast. Carrying `basal_melt_field` in a checkpoint alone does
not make that forcing respond to its analyzed value.

An explicit `basal_melt_forecast: estimated_field` instead propagates the
analyzed melt field on the current floating-ice mask. Initialization draws a
member-keyed multiplier from the explicitly supplied
`basal_melt_prior_scale_bounds`, applied to the nominal wrong-forcing field.
`basal_melt_bounds` supplies the allowed field interval in m/year. These
settings must be provided; their calibration is part of the experiment.
The truth generator continues to use its configured reference forcing.

This policy defaults `initial_perturbation_vars` to an empty list, so generic
initialization introduces no thickness, velocity or surface perturbations.
The melt prior is generated by the model initializer. Forecast noise remains
a separate choice. Identical member IDs and seeds give the same initial melt
multiplier regardless of MPI rank count or member scheduling.

The application derives `analysis_bounds` for `basal_melt_field` from those
explicit bounds. The shared analysis controls apply them after inflation;
the predictor also enforces them before applying the binary floating mask.
The generic and native PIG predictors use the same estimated-field provider.
The forcing policy, prior, observation-error mode and analysis controls enter
the experiment identity. A changed experiment requires a fresh output path.

Set `enkf_observation_error_mode` explicitly. `stochastic_R` uses the configured
observation standard deviations; `legacy_prior_anomalies` uses forecast
anomalies as its error factor and has different statistical semantics.
Neither changing the filter nor increasing rank count creates physically
missing ensemble uncertainty. Check error against truth, ensemble spread,
forcing persistence and restart equivalence before interpreting a new
estimated-forcing experiment scientifically.

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

The default scheduler allocation is one node, eight MPI ranks, 48 GB and a
two-hour cap. Member count, execution mode and ranks/model come from the YAML:

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

The launcher always reads this example's `params.yaml`, unless you explicitly
select another file with `-F` (or `--force-params`):

```bash
sbatch run_idealized_1.sbatch -F "$HOME/r-arobel3-0/ICESEE-run-configs/idealized_pig.yaml"
```

There is no automatic external-file fallback. Legacy `ICESEE_PARAMS_FILE`,
`ICESEE_NENS`, `ICESEE_MODEL_NPROCS`, `ICESEE_DT` and `ICESEE_NUM_YEARS` exports
are ignored with a message; edit the selected YAML instead. Other script
arguments are rejected rather than silently ignored.

Set `execution_mode`, `Nens`, `model_nprocs` and the observation schedule in
`enkf-parameters`. Set duration, profile years and `timesteps_per_year` (this
example's **dt in years**) in `modeling-parameters`. The launcher no longer
forces Mode 3 or the old 24-member/five-year observation settings. For Mode 3,
set `execution_mode: 3`; use `model_nprocs: 2` for two spatial ranks/model.
No scientific parameter defaults in the repository YAML are changed.

Each job records the selected path in `launch.txt`, the original selected file
as `params.source.yaml`, and the resolved configuration as `params.yaml` in
`pace_runs/modeMODE-JOBID/`. YAML inheritance is resolved once. The model reads
that saved configuration, so later edits to the source YAML cannot change an
already started job. The resolved copy records the absolute output path,
resolved ranks/model and timestep count. Other scientific values retain the
selected YAML settings. Choose `ICESEE_DATA_PATH` explicitly for a different
output directory; it is also recorded in the saved configuration.

For a two-node, 32-member Mode-3 run, first set `Nens: 32`, `execution_mode: 3`
and `model_nprocs: 2` in the YAML, then submit:

```bash
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
schedule when choosing such values. setting `modeling-parameters.timesteps_per_year: 0.1` in the selected YAML
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

For Slurm, use the common launcher with a separate copy of the original YAML.
Set `resume_from_checkpoint: true` and the original absolute `data_path` in
`enkf-parameters`; disable ensemble initialization and truth/observation/background
generation in that copy. Keep the original duration, dt, ensemble, observations,
noise, forcing, input paths and physical policies. Then submit
`sbatch run_idealized_1.sbatch -F /absolute/path/to/resume.yaml` from this example.
The launcher requires an existing native history for resume and retains the
fresh-run existing-target guard.

The generic native policy can retain only the newest two step checkpoints;
the PIG plotting workflow instead sets `checkpoint_every=1` and
`checkpoint_keep_last=0` to retain the full history. Do not enable aggregate
export for an older rolling archive with missing times. These settings affect
storage and restart availability. Arbitrary mesh/topology migration is not established by
the short same-topology restart test.

### Bounded recovery and automatic continuation

The native PIG adapter supports optional forecast subcycling. A shared controller
is reusable by adapters that can restore their full model state. It retries
explicitly identified nonlinear forecast failures with smaller internal steps,
without changing external `dt`, observation times or history columns. Failed
attempts are rolled back; process noise is applied once after a successful
external forecast. PIG profiles are written after the accepted forecast. Solver
settings, thickness floor, analysis acceptance limit and filter are unchanged.

For external `dt=0.05`, this sensitivity configuration permits three halvings:

```yaml
enkf-parameters:
  forecast_recovery:
    max_halvings: 3
    min_dt: 0.00625
    max_substeps: 8
```

These values illustrate the controls, rather than prescribing a validated PIG
configuration. Recovery is disabled when the key is absent. Selecting or changing
it changes checkpoint identity; an incompatible existing run is refused. Validate
it in a fresh case. Analysis rejections, compiler/input errors and process loss
are not numerical forecast retries. PIG integration currently requires Mode 3;
other adapters must supply their own complete rollback and substep callbacks.

For long runs, `run_segment_cycles` or `run_walltime_seconds` requests a stop after
an accepted cycle and its checkpoint. Initialization counts against the wall-time
budget. Reserve scheduler time for the slowest cycle, checkpoint, cleanup and final
aggregate. `checkpoint_every: 1` is required. Signals `SIGUSR1` and `SIGTERM` request
the same stop while the controller is active; a blocked collective or hard kill
still needs a later restart. Intermediate segments retain shards and defer the
full aggregate until completion.

The continuation tool freezes expanded YAML and source identity. It preserves the
forcing horizon and scientific settings, records each submission before contacting
Slurm, reconciles uncertain outcomes and enforces an explicit submission budget.
Scientific failure stops the chain. Transient scheduler failures are eligible for
continuation only if a new committed cycle exists. `prepare` submits nothing;
`tick` may submit one segment and `watch` submits permitted successors:

```bash
# Run in this example, using the normal activated dependency environment.
# The YAML must select Mode 3 and checkpoint_every=1.
python tools/continue_case.py prepare pace_runs/continued-case \
  -F /absolute/path/to/experiment.yaml \
  --data "$PWD/pace_runs/continued-case-data" \
  --dependencies "$HOME/r-arobel3-0/ICESEE-Spack" \
  --nodes 1 --ranks 8 --memory 48G --walltime 02:00:00 \
  --run-seconds 5400 --segment-cycles 100 --max-submissions 20

python tools/continue_case.py tick pace_runs/continued-case
python tools/continue_case.py watch pace_runs/continued-case
```

Choose budgets from measured cycle/output costs and site allocation rules. The
example permits twenty two-hour allocations; this is not a default budget.
Keep one supervisor per case. If its terminal closes, invoking `watch` again
reconciles the recorded job. Inspect a stale `.supervisor-lock/owner.json` and
confirm no supervisor is alive before removing that lock manually. Unknown
submission outcomes are held for reconciliation. Do not edit active-case source.

See the [recovery contracts and validation limits](../../../../docs/forecast-recovery-and-continuation.md).

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

Open `read_results.ipynb` from this example directory. In its configuration cell,
set `results_dir` to the completed run's `data` directory, then restart the kernel
and run all cells. `ICESEE_DATA_PATH`, if set in the notebook environment, overrides
that default. No simulation configuration import, MPI initialization or Firedrake
mesh reconstruction is needed. The notebook provides:

- Member time series, means and 10–90% member ranges for thickness, surface,
  velocity components, speed and basal melt.
- Local absolute error versus sample ensemble standard deviation. These are
  point diagnostics, not domain-wide or area-weighted RMSE.
- Bed, surface, base and thickness profiles at actual saved flowline years.
- Spatial snapshots on verified PIG coordinates, with an explicitly labelled
  DOF-index fallback if coordinates are unavailable.

Use `point_index` for a scalar-space DOF or `point_xy_m=(x, y)` for the nearest
physical location in metres. `time_stride=1` uses every saved timestep; increase
it for quicker exploration of long histories. Snapshot and profile years are
selected from the available saved times. Missing backgrounds are omitted, and
observation markers appear only for variables identified as observed in the
run's saved `params.yaml`. Surface elevation is measured in metres, not SMB units.
The default point is a reproducible active-melt example; it is not selected by
agreement with truth and does not represent the whole domain.

The reader in `ICESEE.src.postprocessing.pig` keeps state/ensemble arrays as
HDF5 datasets. Its `extract_datasets_from_h5` convenience function is for small
metadata files and rejects large histories before allocating them. Existing
custom plots can use `results.point(...)`, `results.snapshot(...)`, or explicit
slices of `ensemble_vec_full` and `ensemble_vec_mean`. Do not apply `[:]` to a
complete ensemble or truth history. `reference_dir` can locate shared truth and
observations if a copied run's original `reference_data_path` is unavailable.

Physical maps require `dof_coordinates.h5` in `serial_reference` order. For
native PIG checkpoints whose stored coordinate fingerprint matches the compact
initial checkpoint, export these coordinates without starting a model:

```python
from ICESEE.src.postprocessing.pig import export_reference_coordinates
export_reference_coordinates(data, "data/extended_beta1000yrs_compact_idx20000.h5")
```

The exporter verifies the native coordinate fingerprint and all field row maps,
refuses incompatible layouts and never replaces an existing coordinate file.
Its verification currently requires the full native coordinate vector to match
the saved shard fingerprints; an arbitrary repartitioned layout can be refused.
The maps show node samples rather than inventing rectangular cells or filling
domain holes. A rectangular template mesh is not a correct PIG map.

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


## Time grid, observations and restart

The initial physical model fields are copied without advancing the forecast.
Column zero therefore represents the initial state at `t[0]`. The existing
configured generic initial perturbations are separate from this correction;
native PIG initial physical fields remain unperturbed.

Forecast cycle `k` advances from `t[k]` to `t[k+1]`. An observation indexed by
state column `j` is used after cycle `j-1`. Modes 1–3 apply configured observation-scheduled
process noise at that same boundary, before analysis. The analyzed
state is saved in column `j`. With `dt=0.05`, observations at years 2, 3 and 4
therefore update state columns 40, 60 and 80, rather than 41, 61 and 81. Base seeds, member identities and coordinate
keying are unchanged. Moving a noise event to its correct cycle changes its
timestep-keyed draw, so corrected runs are not bitwise replicas of old-clock
trajectories.
`analysis_applied` marks these columns in native full-history output. Completed
analysis logs record both the state and observation columns and model time.

Observations at the final simulated state are supported. Initial-time
observations (`obs_index=0`) require a pre-forecast analysis stage, which is not
implemented; choose observation times strictly after `t[0]`. A time-zero
physical initial state and a time-zero observation are distinct concepts.

New runs record `time_grid_contract.json` (`state-columns-v2`). Native run and
noise identities, and generic checkpoint metadata, also distinguish the new
clock. Preserve existing output and choose a fresh `data_path` when upgrading:
older state cannot be silently resumed under the corrected schedule. Different
time grids are rejected on reuse. Previously saved files remain readable and
are not relabelled or repaired by export. A restart consumes observations
already included in the checkpoint and replays the configured stochastic
history before the next due event.

### Historical timing findings

Full-history job 15384953 (`da6f9dde`) preserved a pre-existing initialization
and observation offset. Its initial thickness differed from truth at zero by
up to 40.123212 m but matched truth at 0.05 years within 1.14e-12 m. Its
2/3/4-year observations updated states at 2.05/3.05/4.05. Those archived outputs
remain unchanged. Exact export of their saved values did not establish correct
scientific timing; use a fresh run with the corrected clock for new analyses.

### Corrected-clock validation

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
units. These are historical output-integrity results from before the clock
correction, not validation of that trajectory's time alignment.

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

### Ensemble spread plots

A single arbitrary row is a useful debugging trace, but it does not summarize a
spatial field. Choose reproducible physical locations (for example inland ice,
the grounding zone and floating ice), then plot each variable's ensemble mean
with member quantile bands through time. A packed state row is a degree of
freedom, not necessarily a mesh element; use the field's block metadata and its
reference mesh to select the row. Keep the same physical location across runs.

Also plot per-variable spatial standard-deviation maps at selected times and a
spatial RMS spread curve, `sqrt(sum(w * var_members) / sum(w))`, using sample
variance (`ddof=1`). Use area/lumped-mass weights on a nonuniform finite-element
mesh; uniform weights give a nodal summary, not an area-weighted one. Keep
thickness, velocity components, surface and basal melt separate because their
units differ. In twin experiments compare spread with the ensemble-mean error
against truth on the same locations, times and weights; a narrow band alone
is not evidence of accuracy or calibration.

Use `icesee_ensemble_data.h5` for every variable and timestep. The flowline files
contain only thickness/surface at the configured profile years. Inspect native
`blocks`, `t` and `analysis_applied` metadata rather than assuming equal block
lengths or inferring time from column numbers. A point trace is inexpensive:

```python
import json
import h5py
import numpy as np
import matplotlib.pyplot as plt

# row is the reference-mesh DOF selected for this variable/location.
with h5py.File(data_dir / "icesee_ensemble_data.h5", "r") as f:
    blocks = json.loads(f.attrs["blocks"])
    block = next(b for b in blocks if b["name"] == "h")
    assert 0 <= row < block["global_size"]
    samples = f["ensemble"][block["global_offset"] + row, :, :]
    times = f["t"][:]
mean = samples.mean(axis=0)
lower, upper = np.quantile(samples, [0.1, 0.9], axis=0)
plt.plot(times, mean)
plt.fill_between(times, lower, upper, alpha=0.25, label="10–90% member range")
```

The quantile band describes ensemble members, not a confidence interval for the
mean. Native PIG deliberately starts from identical physical fields, so zero
initial physical spread is expected before the configured uncertainty evolves.

### Physically constrained native analyses

The opt-in `pig_analysis_state_policy: diagnostic_balance` (Mode 3) analyzes
thickness and the basal-melt field, freezes direct analysis increments to
`u`, `v`, and `s`, enforces the existing `hThresh` after every observation update,
then derives surface, flotation/friction and diagnostic velocity from that
geometry. Estimated melt is bounded and masked on the resulting geometry.
Solver settings and initial physical fields are unchanged. Checkpoints/history
record the finalized state, rather than the unconstrained analysis array.

This changes the scientific analysis and requires a fresh experiment identity.
The default `independent_fields` preserves existing runs. A 32-member, 30-year
PACE diagnostic completed all eight analyses and finalized its history. However,
its unweighted nodal thickness RMSE at year 30 was 14.07 m with RMS member spread
only 0.00244 m; speed RMSE was 49.50 m/year. This is runtime and output evidence,
not scientific acceptance. Physical consistency, uncertainty calibration and
forcing-model error require separate checks.

### Short stochastic EnKF comparisons

The [stochastic EnKF validation guide](../../../../docs/pig-stochastic-enkf-validation.md)
describes the prescribed-covariance gain, physical noise controls, depth-law melt
predictor, immutable reference reuse and short diagnostic forks. The filter remains
the stochastic EnKF. These policies are explicit experiment choices; the local
checks do not establish a successful 82-year scientific trajectory.

Use `pace/prepare_diagnostic_config.py --help` to prepare separate configurations.
The `gain`, `balance`, `forcing` and `noise` variants apply changes cumulatively.
`--no-analysis` prepares a matched forecast control. Noise amplitudes, correlation
time, coordinate length, drift and projection limits have no calibrated defaults.
The tool does not submit jobs or overwrite the example's `params.yaml`.


### Local basal melt versus the prescribed maximum

The truth simulation ramps the **maximum amplitude** from 20 to 100 m/year over
the configured forcing horizon. This does not prescribe a straight-line melt
history at every point. `BasalMeltRate` also uses the evolving ice-base draft
`z = surface - thickness`, a depth-dependent factor and a floating-ice mask:

```text
local_melt = maximum_amplitude * depth_factor(draft, scenario) * floating_mask
```

For the `control` scenario, melt is zero for drafts at or above -450 m and reaches
the full amplitude at or below -500 m, with a linear transition between them.
The `warm` thresholds are -400 m and -450 m. The truth generator switches between
these scenarios at its prescribed times. Consequently, local melt can jump when
scenario thresholds change and can fall to zero as thinning raises the ice base
above the shallow cutoff, even if the point remains floating. Zero melt alone
cannot establish that a point has grounded. Plotting a fixed physical point does
not prevent its grounding state or draft from evolving.

Melt enters the thickness equation through net accumulation
`(surface_mass_balance - basal_melt) * water_to_ice`. Removing local melt reduces
that direct thinning contribution, while ice-flux divergence and surface mass
balance can still change thickness. Thickness/surface changes also affect the
diagnostic velocity solution. These feedbacks are distinct from numerical blow-up;
check extrema, conservation and timestep sensitivity before drawing that conclusion.

For early-stop diagnostics, `run_stop_time` does not shorten the original forcing
horizon or accelerate its amplitude ramp. A saved post-step truth column contains
the melt used over the preceding integration interval, evaluated from that
interval's starting geometry. Compare forcing, geometry and times accordingly.
If the intended experiment requires strictly linear local melt, removing the depth
cutoff or changing scenario switches is a scientific change requiring a separate
truth/observation dataset. A curve matching the implemented law does not by itself
establish that the law matches the intended experiment.

### Fixed initial melt pattern

`basal_melt_spatial_policy: evolving_depth` is the default behavior described
above. The optional `fixed_initial_depth` policy freezes the initial **control**
depth factor and applies the current floating-ice mask:

```text
local_melt(t) = maximum_amplitude(t) * initial_control_depth_factor * floating_mask(t)
```

For an 82-year forcing horizon, the truth amplitude increases from 20 to 100
m/year; the wrong-forcing reference retains amplitude 20 m/year. Both use the
same initial spatial pattern. Scenario switches and subsequent draft changes
do not change that pattern. Local melt increases linearly while the point stays
floating, but grounding still sets it to zero. Locations with zero initial depth
factor retain zero melt.

Changing this policy requires new truth and synthetic observations in a separate
output directory. Restarting an old trajectory at year 60 retains the first 60
years of response to the old forcing; that is a forcing-intervention experiment,
not the fixed-pattern baseline from year zero. Reuse the existing spin-up and
installed dependencies, but regenerate the reference trajectories from their
initial state. Independent truth and unassimilated wrong-forcing/background
trajectories can run concurrently in separate processes and output directories;
neither may overwrite the other's open HDF5 files. The new policy is under
validation; its availability does not establish long-run scientific acceptance.

To prepare independent jobs without allocating an ensemble, resolve a base YAML
into separate fresh configurations:

```bash
python pace/prepare_reference_jobs.py params.yaml /absolute/path/to/new-reference-run
sbatch pace/run_reference.sbatch -F /absolute/path/to/new-reference-run/truth/params.yaml
sbatch pace/run_reference.sbatch -F /absolute/path/to/new-reference-run/background/params.yaml
```

Set `ICESEE_SOURCE_ROOT` to the source checkout containing this implementation
and `ICESEE_SPACK_ROOT` to the installed dependency environment. The prepared
configs retain the base duration, timestep, physical and observation settings;
they select the fixed initial melt pattern. Each job uses the direct
`tools/generate_reference.py` driver on one MPI rank, without assimilation
dispatch or ensemble allocation. The truth job also
generates observations; the background job generates only `nurged_state`.
Their output directories are independent. A `reference_manifest.json` is
published only after generation returns and its trajectory shape is checked.
It records generation completion, not scientific acceptance. Check exit status,
all time columns, finite values and physical ranges before using new references.

### Optional relaxation to prior spread

`analysis_spread_relaxation` and `analysis_spread_relaxation_factors` control RTPS
within the existing stochastic EnKF. They are separate from
`analysis_relaxation_factor`, which damps analysis increments. RTPS preserves the
row's analyzed mean before subsequent bounds/model finalization and targets
`SD_new = (1-alpha)*SD_analysis + alpha*SD_forecast`. Coefficients must be in
`[0, 1]`; the default is zero, preserving existing behavior. Per-field settings
are available through the shared bounded-row analysis controls in Modes 1–3.

For a **sensitivity test**, not a calibrated production configuration:

```yaml
enkf-parameters:
  analysis_spread_relaxation: 0.0
  analysis_spread_relaxation_factors:
    h: 0.5
    basal_melt_field: 0.5
  analysis_stage_diagnostics: true
```

The spread step follows existing multiplicative inflation, then existing bounds,
frozen-field rules and physical finalization remain authoritative. Thus the final
physical ensemble can have different moments from the unconstrained RTPS result.
If analyzed spread is exactly zero, the centered forecast direction supplies the
requested spread; if the forecast also has none, no artificial noise is created.
No RNG draws or rank/member reassignments occur. Enabling/changing RTPS changes
the PIG scientific experiment identity and cannot silently resume an old run.

`analysis_stage_diagnostics` currently requires persistent native Mode 3. It
records `forecast`, `raw_analysis`, `controlled_analysis` and `balanced_analysis`
for each analysis at the **same model time**, under
`<data_path>/_analysis_diagnostics/analysis_<cycle>.json`. Reports contain per-field
nodal RMSE against cached truth, RMS sample member spread, extrema and nonfinite
counts. They read one truth snapshot on ensemble roots and use bounded member
rows; no full ensemble history is loaded. These are uniform owned-row summaries,
not area-weighted ice-sheet metrics. The balanced stage precedes any optional
inversion. Pair a report with its committed checkpoint; a report alone does not
prove the cycle was committed.

Use `pace/prepare_diagnostic_config.py --spread-relaxation 0.5
--spread-fields h basal_melt_field --stage-diagnostics` with the other required
arguments to create a separate test configuration. For paired short runs sharing
one Slurm allocation, `ICESEE_LOG_PATH` selects a fresh per-case log directory;
its default remains `pace_runs/mode<mode>-<jobid>`. Logs must sit outside the
recreated dataset directory. No existing YAML or running source is changed by
preparing a diagnostic configuration.
