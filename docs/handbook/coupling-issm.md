# Coupling ISSM to ICESEE

We couple ICESEE to ISSM through Python application adapters, persistent MATLAB subprocesses, and member-specific HDF5 state exchange. ISSM performs its own finite-element forecast and inversion; ICESEE schedules ensemble work, constructs observations and analysis updates, and writes the assimilation history. The Python-visible exchange is a complete member. It does not expose the partition of ISSM's internal parallel solve as owned spatial slabs.

This chapter describes inspected source at scientific commit `3748e416d3a8120cc860119bd7224c6f9ad41ff7`. The later [compact basal-friction runs](sample-runs.md) execute the actual MATLAB/ISSM backend with the separately identified directory/clock corrections. Their case-specific evidence does not validate every row below or the Choi inverse. Existing unit tests cover field translation, resource planning and stubbed native orchestration; these remain distinct from a coupled backend run. For the ISSM model itself see [Larour et al. (2012)](https://doi.org/10.1029/2011JF002140). The [parallelization chapter](parallelization-modes.md) explains the shared runners; [deployment](deployment.md) describes numerical environments and containers.

## 1 Select the actual ISSM application

`applications/supported_models.py` resolves `model_name: issm`. It first tries an example-local `_issm_enkf` module using the current application-directory name, then falls back to `issm_model/issm_utils/_issm_enkf.py`. The example drivers change into their own directory before this lookup. A model name alone therefore does not identify the active callbacks; record the driver and resolved adapter path.

| Application | Driver and distinguishing features | Native registration |
|---|---|---|
| `ISMIP_Choi` | Explicit example-local imports; six nodal state/parameter blocks; graph noise/coordinates; delayed friction inversion and geometry controls | Driver imports its `mode3_runner`, which registers the ISSM whole-member route |
| `basal_friction_variation` | Legacy example-local wrapper; different parameter names/case; explicit `var_nd`; additional scalar diagnostics | Current driver does not import a native registration |
| `ISMIP` | Older example retained in the repository | Inspect its own files; it is not the maintained Choi native path |

The basal-friction example's checked-in YAML currently says `execution_mode: 3`, while its driver does not register that route. Select a reviewed generic mode explicitly when using this application. Copying the Choi mode number or importing its runner elsewhere does not translate the basal example's fields or policies. Existing README snippets using `--parallel=serial`, Python 3.8, or `Nens=1` for an EnKF are historical; the package requires Python 3.11+, the runner selector is `--execution_mode`, and covariance-based ensemble calculations require an appropriate nontrivial ensemble.

## 2 Trace initialization and the MATLAB boundary

Both main drivers follow this sequence:

1. Load the selected YAML and CLI overrides into `icesee_kwargs`, initialize the Python MPI resource context, and obtain the initial member identifier.
2. Resolve `ISSM_DIR`, add ISSM's interface paths, and prepare `$ISSM_DIR/examples/<example_name>` through `setup_example_directory`.
3. Convert scientific/model settings, construct the physical time grid, and write a MATLAB-compatible `icesee_kwargs_<ens_id>.mat` through `matlab_icesee_kwargs` and `scipy.io.savemat`.
4. Copy `issm_env.m`, `matlab_server.m`, `Domain.exp`, the selected `ParamFile` and metadata into the ISSM example workspace.
5. Construct a `MatlabServer`, invoke the application's `initialize_model`, and read the initialized variable size from `ensemble_init_<ens_id>.h5`.
6. Return to the Python application directory and dispatch the selected ICESEE runner.

`initialize_model` calls `setup_reference_data` and `setup_ensemble_data`. These helpers contain `MPI.COMM_WORLD.Barrier()` operations, even when a caller supplies a smaller communicator. Initialization must therefore occur once in a matched world-rank sequence, not independently inside a per-member loop reached by an arbitrary rank subset. Root prepares the base reference/`Models/ens_id_0` data and metadata for all ensemble identifiers.

The initialization MATLAB file in the Choi application loads configuration but also sets the active branch to `steps = [55]`, the reference-model initialization branch. Changing only the YAML `steps` value does not select an arbitrary mesh-generation branch in that file. The referenced `.mat` model must exist; the example is not a self-contained download-and-mesh command.

### The server uses files rather than network ports

The active implementation is in `issm_utils/matlab2python/mat2py_utils.py`. `MatlabServer` lazily creates `_MatlabServer` when a delegated operation is first used. Each server launches headless MATLAB as a subprocess and communicates through absolute files named `cmdfile_<color>.txt` and `statusfile_<color>.txt` in its workspace. There is no TCP port or `nservers` setting in this interface.

Python atomically publishes a command and a `running` status. MATLAB reads and evaluates the command, atomically publishes `done` or `error`, and deletes the command file. Python verifies `done` before accepting a returned forecast, so deletion alone cannot make a failed MATLAB command appear successful. Startup also checks readiness and whether the subprocess exited. The command timeout default in the implementation is very long; do not infer a short deadlock detector from its older docstring.

`shutdown()` sends `exit`, waits for the child and can terminate its process group if necessary. The current process cleanup restricts stale-server termination to command/status paths owned by that server. Preserve these paths and the MATLAB process identity in incident records. Do not run a blanket user-wide MATLAB kill command to stop one failed member.

## 3 Separate Python ranks from ISSM solve resources

We distinguish three controls:

| Control | Meaning in the inspected implementation |
|---|---|
| Python launcher rank count, such as `mpiexec -n 4` | Number of ICESEE Python processes and their associated application/server setup |
| `ranks_per_model` | Generic Modes 1/2 model-communicator grouping; ISSM registers no verified multi-rank Python state decomposition |
| `model_nprocs` | Passed to MATLAB and used as ISSM's `md.cluster = generic(..., 'np', nprocs)` for its internal solve |

Increasing `model_nprocs` does not make ISSM return a local state slab to each Python rank. It can launch additional solver processes from MATLAB. Budget Python workers, MATLAB workers, internal solver processes and threads against the actual allocation. A launch with four Python ranks and four internal solver processes per active solve needs substantially more resources than four cores; the exact overlap depends on the selected solve and site launcher.

`MatlabServer._check_conditions` requires `Nens >= MPI.COMM_WORLD.size` when a communicator is provided. Generic resource planning can schedule additional members in rounds; equality is not an established optimum. For a first backend test request one Python rank per active server and `model_nprocs=1`, then inspect the effective outgoing solver commands. The generic Mode-1 resource policy can increase this value from available host/SLURM resources and increase it again after ensemble initialization: a real one-Python-rank basal test requested one but used three for reference generation and thirteen for forecasts on a twelve-core host. The teaching profiles use the existing `ICESEE_PERFORMANCE_TEST: true` resource-policy switch to retain the requested solver count; the historical name does not turn their timings into a benchmark. Verify the observed count before increasing nested parallelism.

The generic planner accepts `issm_server_count` as a cap on active model groups. This is the actual key corresponding to a requested server-count limit; `nservers` is not an alias. The planner's spare-rank tests establish topology behavior. The ISSM example drivers still initialize servers/meshes unconditionally using the returned communicator and member identifier, so a cap that creates spare Python ranks requires application-level validation before production use. Native Choi constructs its own `spatial_ranks=1` topology and does not use this generic cap to limit its groups. Do not treat a parsed key as evidence that every application policy implements it.

## 4 State packing and physical coordinates

In `ISMIP_Choi`, the configured order is:

```text
Thickness | Surface | Vx | Vy | bed | coefficient
```

Each block has `Nv = md.mesh.numberofvertices` nodal values. The driver sets `nd = Nv * total_state_param_vars`; `icesee_get_index` maps the variable-major blocks. MATLAB writes complete arrays, and the Python return path flattens them with Fortran ordering. `coefficient` is the model's friction-coefficient field; its interpretation follows the configured ISSM sliding law, not a generic basal-stress unit assumption.

The basal-friction application instead configures:

```text
Thickness | Surface | Vx | Vy | Bed | FrictionCoefficient
```

These names are case-sensitive dataset and index keys. Its `scalar_inputs` include `IceVolume` and `IceVolumeAboveFloatation`, but those names are absent from `vec_inputs`: they are separately read diagnostics, not one-element state blocks. Its `var_nd` construction therefore gives all six packed blocks the same nodal size. `test_issm_basal_friction_variation_layout.py` covers this distinction. If a new application genuinely adds a scalar state block, set both `var_nd` and `nd = sum(block_sizes)` consistently; changing `scalar_inputs` alone does not update the driver's uniform dimension formula.

The shared coordinate provider reads `fric_x` and `fric_y` from `mesh_idxy_0.h5`, validates their sizes and returns coordinates in nodal state order. Choi uses this physical graph for random fields and localization. Preserve coordinate order, units and mesh identity when replacing a reference model. An observation at one vertex index is not transferable to another mesh merely because both arrays have the same length.

## 5 Forecast and inverse callbacks

`forecast_step_single` sets `tinitial=t[k]` and `tfinal=t[k+1]`, then calls `run_model` with one packed member. Python writes that member's named blocks to `ensemble_output_<ens_id>.h5`, asks MATLAB to run `run_model.m`, and reads the complete returned member. MATLAB selects its branch using the requested state filename, updates the ISSM model, solves the transient and writes the expected datasets. `output_frequency=1` is part of that exchange: the Python cycle expects the corresponding model result.

Choi's forecast wrapper raises on a failed MATLAB command and refuses to reuse a stale output file. The basal-friction wrapper has older error-return behavior and needs separate failure-path validation. Example-level outer exception handlers also print errors and call `exit()`; collect simulation-step logs and inspected output rather than accepting an exit code alone as backend success.

`inverse_step_single`/`run_model_inverse` provide the hybrid friction-inversion path. The analyzed geometry and friction fields are written back to MATLAB; the configured velocity observations, objective weights, friction bounds and regularization define the inverse solve. The driver preserves `inversion_enabled` separately from the cycle-dependent `inversion_flag`, with `inversion_start_time` selecting when inversion starts. Choi native uses the observation/model cycle-time helper before setting this flag. Retain `km` and the observation identity when invoking an inverse callback.

Geometry projection, increment relaxation/caps, frozen variables and bed-update domains affect the analyzed state. They are part of the scientific configuration. In particular, changing bed, surface and thickness independently can violate geometry or flotation relationships unless the model's configured consistency step restores them. The Choi MATLAB initialization includes a geometry-consistency construction; retain and test it when designing a different prior. Callback availability does not establish that every generic analysis or native path applies identical inference/localization controls.

## 6 Coupling across all execution modes

| Mode | Inspected ISSM route | Validation and limits |
|---|---|---|
| 0 serial | Generic serial runner calls the selected ISSM callbacks through one Python process/server | Suitable starting architecture for a small coupled reference; serial filter branches need their own numerical tests |
| 1 member parallel | Generic in-memory orchestration; default EnKF factor is rooted and row updates are distributed | Complete-member MATLAB exchange remains; nested ISSM processes are separate from Python MPI grouping |
| 2 file backed | Generic `default_run` forecast/I/O and streamed analysis with ISSM callbacks | Reduces ensemble/history allocations in ICESEE; does not remove the complete-member exchange or ISSM solve memory |
| 3 native Choi | Registered Choi runner, `WholeMemberNativeAdapter`, `spatial_ranks=1`, scheduled members and native global-analysis/checkpoint calls | Complete fields remain on each assigned Python worker; no Python-visible spatial ISSM partition; not the PIG spatial policy |

All four rows describe source-inspected routes. The later compact basal-friction backend evidence concerns its selected generic configuration; it does not validate every row or the Choi native path. The Choi native tests replace heavy dependencies with stubs and verify registration, orchestration, field translation and delayed inversion; they are not a finite-element accuracy or production-scale test.

### Native Choi boundaries

The driver imports its native runner before dispatch. That runner initializes an additional server/setup context, constructs native member fields, evaluates observations and calls `run_native_global_analysis_cycle`. It accepts `legacy_prior_anomalies` and `stochastic_R`, and rejects other error modes including `generated_R`. It calls the global-analysis cycle even when a generic YAML also contains local-analysis settings; do not infer the generic localized/inference policy from those keys alone.

The whole-member adapter has no `reactivate_native_member`. Therefore `initialize_native_member_pool` uses the persistent pool containing all locally assigned member states; `use_member_streaming` does not automatically give it a streaming implementation. Persistent MATLAB-worker count and resident NumPy-member count are different memory measures. Historical headers claiming a strict group-count bound on every resident member are not a verified property of this implementation.

The runner writes an initial checkpoint and a checkpoint after every cycle under `_mode3_state_history`, but the inspected loop has no checkpoint-load/resume branch. Writing the shared checkpoint format does not provide the PIG runner's restore workflow. Preserve these shards for inspection; do not set resume flags expecting native Choi to reconstruct its MATLAB model and continue. A coupled restart needs an application restore policy and a real interrupted-run test.

## 7 Prepare an isolated application workspace

The coupling writes in two places: the Python example directory (metadata and member exchange) and `$ISSM_DIR/examples/<example_name>` (server command/status files and `Models/ens_id_N`). A unique `data_path` does not isolate the latter, and even a unique `example_name` does not isolate metadata files written in the Python application directory. Use a separate source/application workspace and ISSM example workspace for concurrent runs.

Many wrapper paths are built as `icesee_path + '/' + data_path`. An absolute `data_path` is not consistently honored across this bridge. Use a unique **relative** data path under the isolated application, or arrange a reviewed parent link/mount to shared project storage. Record the resolved host/container path. Preserve valuable output before loading fresh-run configuration, which can delete the selected data directory.

Required items include a compatible `ISSM_DIR`, MATLAB on the executable path and its license configuration, the reference model from `reference_data_dir/reference_data`, `Domain.exp`, `ParamFile`, required prior files such as Choi's `bed_kriging_results.h5`, and writable ISSM examples/execution storage. Initializer and model scripts may reference additional prepared model files; retain their hashes with the selected experiment. Source updates do not generate missing spin-up results.

A container must bind the isolated Python source, writable ISSM examples/execution paths and all referenced input/output paths. The combined container's ISSM and Icepack wrappers select different MPI stacks; verify the actual Python MPI library as well as ISSM's internal executable before using a nested launch. The [deployment chapter](deployment.md) gives these environment boundaries.

## 8 Launch procedures using existing drivers

These are source-checked launch forms, not a newly executed ISSM validation. Begin with a reviewed personal copy of the full example YAML and its required data; the small Lorenz YAML is not an ISSM template.

```bash
# From an isolated checkout directory named ICESEE, after stack activation:
export ICESEE_ROOT="$(pwd -P)"
export PYTHONPATH="$(dirname "$ICESEE_ROOT"):$ICESEE_ROOT:${PYTHONPATH:-}"
: "${ISSM_DIR:?Activate the approved ISSM installation}"
command -v matlab
command -v mpiexec
CASE_CONFIG_DIR="$(mktemp -d "${TMPDIR:-/tmp}/icesee-issm-config.XXXXXX")"
cp applications/issm_model/examples/ISMIP_Choi/params.yaml "$CASE_CONFIG_DIR/params.yaml"
export ISSM_PARAMS="$CASE_CONFIG_DIR/params.yaml"
```

Edit this personal YAML before launch: retain the intended state/order and reference inputs; choose the experiment length, observation schedule, ensemble spread and inversion policy; set a unique ISSM `example_name` and relative `data_path`. The Choi driver takes its time increment from `timesteps_per_year`, despite that key's historical name. A value `0.2` is used as a 0.2-year step, not five such steps implied by an integer steps-per-year interpretation.

Serial and in-memory examples use the actual Choi driver:

```bash
cd "$ICESEE_ROOT/applications/issm_model/examples/ISMIP_Choi"
python run_da_issm.py -F "$ISSM_PARAMS" --execution_mode=0 \
  --Nens=4 --model_nprocs=1 --data_path=_issm_runs/serial_case01
mpiexec -n 2 python run_da_issm.py -F "$ISSM_PARAMS" --execution_mode=1 \
  --Nens=4 --model_nprocs=1 --data_path=_issm_runs/member_case01
```

The file-backed and native launch forms are:

```bash
mpiexec -n 2 python run_da_issm.py -F "$ISSM_PARAMS" --execution_mode=2 \
  --Nens=4 --model_nprocs=1 --data_path=_issm_runs/file_case01
mpiexec -n 2 python run_da_issm.py -F "$ISSM_PARAMS" --execution_mode=3 \
  --Nens=4 --model_nprocs=1 --data_path=_issm_runs/native_case01
```

Run each in its isolated workspace and inspect the recorded source, server status, MATLAB logs, time grid, member datasets and model outputs. The native command establishes the registered whole-member route, not spatial ISSM ownership or restart support. To use `basal_friction_variation`, copy its own full YAML, invoke its own driver directory, choose Mode 0/1/2 explicitly and retain its case-sensitive field names. Do not add a Choi registration import as a shortcut around application capability checks.

## 9 Validate before extending the coupling

The minimum coupled backend check initializes the intended mesh/reference, exchanges one complete member, advances one forecast interval, performs one supported analysis/inversion and verifies physical state/time and field ordering. Include a MATLAB command that fails deliberately and confirm that the launcher reports failure without consuming stale HDF5. Compare Modes 1/2 using aligned initialization, observations and complete members before interpreting a performance difference. Native whole-member consistency needs its own comparison to the matched generic experiment.

Run lightweight tests separately from backend tests. `test_issm_ismip_choi_native.py`, `test_issm_mode3_runner.py`, `test_issm_basal_friction_variation_layout.py` and the resource-plan tests answer structural questions with stubs or small arrays. A new application's [onboarding sequence](application-onboarding.md) should add shape/index contracts first, then real MPI/backend/failure tests. Record which tests were skipped because MATLAB, ISSM or input files were unavailable.

From the checkout root, the focused structural tests can be invoked in the configured test environment:

```bash
python -m pytest -q \
  src/tests/test_issm_ismip_choi_native.py \
  src/tests/test_issm_mode3_runner.py \
  src/tests/test_issm_basal_friction_variation_layout.py
```

For a genuinely spatial ISSM adapter, the Python boundary would need stable owned-row/global-ID metadata and model-supported local exchange. Increasing internal `model_nprocs`, changing an ownership label or wrapping a full array in a native registry does not create that interface. Such work requires an explicit design and model-level validation rather than a documentation claim.

## Source reading map

- `applications/issm_model/examples/ISMIP_Choi/run_da_issm.py`, `_issm_model.py`, `_issm_enkf.py`, `_issm_native.py`, `mode3_runner.py`, `initialize_model.m`, `run_model.m` and `params.yaml`
- `applications/issm_model/examples/basal_friction_variation/` for the separate driver, field names and diagnostics
- `applications/issm_model/issm_utils/matlab2python/mat2py_utils.py`, `matlab_server.m`, `issm_env.m` and `server_utils.py`
- `applications/issm_model/issm_utils/_coordinates.py` for physical nodal coordinates
- `src/parallelization/parallel_mpi/{resource_plan,model_capabilities,icesee_mpi_parallel_manager}.py`
- `src/parallelization/distributed_whole_member_adapter.py`, `distributed_native_runtime.py` and `distributed_native_cycle.py`

## Corrections exposed by the teaching runs

The later [ISSM runtime correction](https://github.com/ICESEE-project/ICESEE/pull/74) restores the basal application directory before HDF5 generation returns and exports the actual pre-solve reference state as column zero. Subsequent MATLAB states map to identical Python columns, including the terminal state. It preserves the transient solves and the existing prior spin-up. The [sample-run record](sample-runs.md) reports actual coupled acceptance separately from these source and unit-test checks. Choi's native missing-inverse-result guard is covered by the same correction; it is distinct from MATLAB optimizer validation.
