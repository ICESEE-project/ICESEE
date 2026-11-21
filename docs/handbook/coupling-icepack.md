# Coupling Icepack to ICESEE

We couple Icepack at two boundaries. The generic runners exchange packed arrays with an example's forecast hooks; the spatial native PIG runner retains Firedrake fields and exchanges owned entries with the analysis. Both boundaries require a precise definition of the physical fields, observation numbering and model time. Sharing a model name does not make their initialization, memory use or available scientific policies identical.

This guide describes source inspected at runtime commit `3748e416d3a8120cc860119bd7224c6f9ad41ff7`. **Source inspected** means the route and its checks exist in that source. **Measured** refers to the retained PACE experiments and the separately recorded local teaching fixture. A source route or a backend test file is not evidence that every configuration has been run successfully. The [Icepack model paper](https://doi.org/10.5194/gmd-14-4593-2021) and [Firedrake paper](https://doi.org/10.1145/2998441) describe the forecast software; the ICESEE adapter and its validation are additional layers.

## Examples and execution boundaries

| Application and mode | Boundary inspected in source | Evidence and limitation |
| --- | --- | --- |
| PIG, Mode 0 | Serial runner and example-specific packed hooks | Source inspected. The corrected-clock five-year PACE comparison did not include this serial route. |
| PIG, Mode 1 | In-memory ensemble orchestration with packed member forecasts | Measured: eight ranks, one rank/model, 32 members and five years in job 15433074. |
| PIG, Mode 2 | File-backed ensemble orchestration and streamed analysis; forecast hook still exchanges member arrays | Measured: eight ranks, one rank/model for the same case. Bounded analysis does not make every model-side temporary spatially distributed. |
| PIG, Mode 3 | Persistent native Firedrake fields on spatial model groups, owned-state analysis and sharded checkpoints | Measured: sixteen ranks, two ranks/model for the same case, plus short restart evidence. This runner uses the global native analysis or its store-streaming counterpart. A localized native PIG policy is not wired by that selection. |
| Synthetic ice stream, Modes 0–2 | Example-specific generic hooks and rectangular-mesh initialization | Source inspected and the new coarse four-member, ten-cycle teaching case executed; see the [sample runs](sample-runs.md) for its exact scope and all-mode evidence. |
| Synthetic ice stream, Mode 3 | `WholeMemberNativeAdapter` with `FileBackedNativeAnalysis` and one spatial rank per member | The new coarse case executed this binding with two MPI ranks, one spatial rank per member and two analyses. It retains a complete member on its rank and refuses multi-rank spatial execution and requested resume through the shared policy. It is distinct from PIG's spatial native runtime. |

All modes select execution and storage policy, rather than implicitly selecting a new filter. PIG's native runner accepts `legacy_prior_anomalies` and `stochastic_R`; it rejects `generated_R`. Its call sites do not select the lower-level local-analysis machinery. Synthetic Mode 3 registers the shared file-backed analysis policy and declares localization, inference and inversion features; those declarations are not a five-year PIG validation or proof of every feature combination. The [parallelization guide](parallelization-modes.md) and [numerical chapter](numerics.md) explain these separate choices.

## Entrypoint, setup and communicator lifetime

The PIG entrypoint is [`run_da_icepack.py`](../../applications/icepack_model/examples/idealized_pig/run_da_icepack.py). Its configuration is loaded before Firedrake, so ICESEE's CLI arguments are consumed before PETSc initializes. Thread limits are also set before NumPy imports. The entrypoint imports `mode3_runner` to register the specialized policy. PIG does not perform the same early `chdir` as the synthetic driver, so launching from its example directory is important to legacy hook discovery and relative inputs.

For Modes 0–2, `ParallelManager.icesee_mpi_init` supplies the model communicator and the entrypoint builds the legacy context. [`initialize_model`](../../applications/icepack_model/examples/idealized_pig/_icepack_model.py) calls `initializeMesh`, constructs the Icepack model and solver, and calls `initializeRun`. The context includes mesh/function spaces, physical fields, forcing, reference geometry and the solver. The generic callbacks subsequently receive this context through `icesee_kwargs`.

Mode 3 follows a different setup sequence. World root builds a serial reference using Firedrake `COMM_SELF`, generates or reuses truth and observation artifacts, and establishes serial reference coordinates. The runner then creates the spatial/ensemble topology. [`_icepack_native.py`](../../applications/icepack_model/examples/idealized_pig/_icepack_native.py) constructs a shared model context once per spatial group using `topology.spatial_comm`, and keeps the members' physical functions alive on that context. Meshes, checkpoint reads and solvers must use the communicator that actually owns those fields. A world communicator substituted inside a per-group callback can produce unmatched collectives.

`initializeMesh` reads the checkpoint mesh and passes the selected communicator to mesh construction. `initializeRun` opens Firedrake `CheckpointFile` collectively on that communicator. With `compact_initialization: true`, it loads named final fields without a history index; the historical input uses `idx=20000` for velocity, thickness and surface. A single-rank context uses Firedrake's managed `COMM_SELF`. PIG also installs a consistent PyOP2 compile-cache policy before model initialization to address unequal compiled-cache state across model ranks.

Native PIG validates node-local model groups before group checkpoint loading. Source, input paths and local scratch/cache placement therefore remain part of a correct launch. The context cache is associated with the topology and has an explicit release function; a new application must not retain model objects whose communicator has already been freed. The entrypoint's `PETSC_CONFIGURE_OPTIONS` assignment is not evidence that a running environment was rebuilt; the loaded MPI/PETSc/Firedrake libraries must be recorded independently.

## Packed state and persistent fields

PIG's default packed order is `h`, `u`, `v`, `s`, `basal_melt_field`. The two velocity components occupy separate scalar-sized blocks even though Firedrake stores velocity in a vector function. Thickness and surface have length units; velocity, basal melt, surface mass balance and density conversion use the application's configured physical units. A new forcing file needs a documented unit conversion, rather than a matching array shape alone.

The default PIG YAML counts all five entries in `num_state_vars`, sets `num_param_vars: 0`, and places surface in `observed_params: ['s']`. These are existing configuration labels, not a physical classification of surface as a static parameter. The combined observation list determines selected rows. Basal melt remains a packed slot when joint estimation is disabled; the state definition and the joint-estimation nudge are separate decisions. The generic initializer returns a dictionary of blocks; forecasting unpacks the array into Firedrake functions and returns updated blocks in the same order.

The spatial bridge is [`_distributed_fields.py`](../../applications/icepack_model/icepack_utils/_distributed_fields.py). `scalar_ownership(function)` reads the scalar PETSc layout vector's global size and ownership range and checks it against the owned `Dat` view. `build_icepack_field_registry` creates named `h`, `u`, `v`, and optional `s`/basal-melt handles. Each block has its own global offset and owned interval. A local packed state concatenates only owned pieces; it does not contain the gaps between globally separated variable blocks.

Firedrake ghost entries support local model computation and are not additional global state. `read_owned` and `write_owned` use the owned views; synchronization refreshes halos after updates. The two velocity handles share a backing `Dat`, and the registry avoids redundant synchronization using the handle's synchronization key. A new finite-element space must satisfy the component and ownership assumptions of this bridge; the current use of scalar thickness ownership for velocity components is not a blanket mixed-space guarantee.

## Observations and mesh numbering

Truth and synthetic observations use a serial reference numbering. A distributed mesh can represent the same physical nodes in a different numerical order. PIG matches owned node coordinates to the reference, builds a reference-to-layout permutation, translates observation row IDs, and routes them to the rank that owns each row. It retains `obs_positions` so the observation values and error factors are subset in the same canonical order.

Equal vector lengths do not establish equal numbering. The required checks include unique/gap-free ownership, matching physical coordinates, variable offsets and preserved observation ordering. The same canonical mapping is used for configured process-noise increments. Reference coordinate setup and the canonical noise bridge still involve global-sized metadata/workspace; native physical member fields remain distributed, but total application memory is not described by the owned-field size alone.

## Forecast, analysis and recovery

The generic PIG hook has signature `forecast_step_single(ensemble=None, **icesee_kwargs)` and delegates to `run_model`. The native bridge instead mutates persistent member functions through `forecast_native_member(member, timestep, *, topology, icesee_kwargs)`. Both should advance exactly one resolved interval. PIG interprets `timesteps_per_year` as a step length in years, whereas the synthetic driver uses a count per year and computes its reciprocal. A copied YAML field with the same name can therefore change the time grid unless the driver resolves it explicitly.

PIG forecasting recomputes forcing for the configured true/wrong experiment, advances thickness prognostically, applies the forecast thickness policy, recomputes surface/flotation and grounding-line-dependent friction, then solves diagnostic velocity. `forecast_member` recomputes the floating mask from the current member surface and bed rather than retaining the bootstrap mask indefinitely. The reference inflow thickness belongs to the shared setup; evolving member fields belong to the member. The current native forecast regenerates basal melt through `_basal_melt_for_step` before the advance. A posterior melt slot can therefore be checkpointed without being carried unchanged as the next cycle's forcing. A parameter-estimation experiment needs to account for that application-specific evolution policy.

An important current boundary is that PIG does **not** supply a state-modifying `finalize_analysis` callback to `IcepackNativeAdapter`. Its finalization method synchronizes ghosts. An analysis can therefore change packed thickness, surface and velocity without immediately rerunning the physical geometry/velocity construction. The subsequent forecast performs those operations. This is a scientific property of the present binding, not evidence that every posterior is already geometrically projected. The observed transient posterior ranges and physical consistency need their own interpretation.

Column zero is the physical initial state. Cycle `k` advances to column `k+1`; scheduled process noise precedes the analysis associated with that completed state. Native PIG starts with unperturbed physical fields, while generic initialization retains its perturbation policy. Sharing a noise generator does not establish equal realized priors across modes. Native checkpoints store owned fields, progress and compatibility metadata. Restore allocates native functions without rerunning initialization perturbations or an extra forecast. Full native history export requires `checkpoint_every: 1` and `checkpoint_keep_last: 0`; the runner rejects incompatible retention/export settings.

## A reproducible coupling exercise

A small case uses a separate reviewed YAML copied from the intended example, valid local input paths, a short duration, at least two members, and interior/final observation and profile times. PIG's retained short validation used one year, `dt=0.05`, two members and observations at 0.5/0.75/1.0 years. The following commands describe source-inspected launch patterns; they are not additional validation results and submit no jobs until executed by an operator in an appropriate allocation.

```bash
cd "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig"
# CONFIG is an absolute path to a reviewed short PIG configuration.
# OUTPUT_ROOT is a new writable run directory, never an accepted old output.
python run_da_icepack.py -F "$CONFIG" --execution_mode=0 --data_path="$OUTPUT_ROOT/mode0"
mpiexec -n 2 python run_da_icepack.py -F "$CONFIG" --execution_mode=1 --model_nprocs=1 --data_path="$OUTPUT_ROOT/mode1"
mpiexec -n 2 python run_da_icepack.py -F "$CONFIG" --execution_mode=2 --model_nprocs=1 --data_path="$OUTPUT_ROOT/mode2"
mpiexec -n 4 python run_da_icepack.py -F "$CONFIG" --execution_mode=3 --model_nprocs=2 --data_path="$OUTPUT_ROOT/mode3"
```

The MPI launcher must first pass a functional compatibility probe in the activated environment. Input files must be readable by every model group. Native groups must fit the permitted node placement. Four-rank PIG commands require the matching Firedrake stack; they are not an invitation to run PDE workloads on an HPC login node. The existing `pace/test_modes_1_3_year.sbatch` and `pace/validate_production_modes.sbatch` provide the recorded PACE case structure, including a fixed-case output verifier.

Before a new coupled experiment, the relevant existing test gates include:

```bash
cd "$ICESEE_SOURCE_ROOT"
python -m pytest -q src/tests/test_icepack_distributed_fields.py src/tests/test_icepack_reference_numbering.py src/tests/test_timestep_contract.py
python -m pytest -q src/tests/test_icepack_checkpoint_communicator.py src/tests/test_icepack_mesh_communicator.py src/tests/test_icepack_idealized_pig_native.py src/tests/test_icepack_mode3_runner.py
```

Backend-dependent tests can skip when Firedrake/Icepack is unavailable; a skip does not establish that backend's behavior. Successful interface checks are followed by model-only forecasts, a small actual coupled run, member-wise comparisons for genuinely matching priors, profile/event checks and a supported checkpoint replay crossing analysis. Common failures arise from mismatched imports (`src...` versus `ICESEE.src...` class identities), a wrong working directory, a world-communicator mesh inside a subgroup, ghost duplication, untranslated observation IDs, stale output reuse, and treating a whole-member synthetic adapter as PIG's spatial adapter.

The next step for a new case is the [application onboarding guide](application-onboarding.md). The [ISSM coupling guide](coupling-issm.md) explains a different boundary, where external solver parallelism does not by itself imply Python-visible spatial state ownership.

The [coarse synthetic exercise](../../applications/icepack_model/examples/synthetic_ice_stream/tutorials/coarse_assimilation/README.md) also ran the registered whole-member native Mode 3. Its finite snapshots and two-analysis evidence do not extend the spatial or restart capabilities of that binding. Modes 1/2 agree across this fixture to a maximum absolute difference of `2.007709554163739e-10`.
