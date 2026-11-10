# Architecture and execution paths

ICESEE combines model forecasts with ensemble analysis through a flat runtime context and application hooks. Its scalability depends on three different objects: the state needed by one model member, the rows needed by one analysis block, and the ensemble history retained over time. Keeping one of these bounded does not automatically bound the others. This distinction is essential when sizing production jobs and when describing the updated architecture.

## Code map

| Responsibility | Active implementation and reading notes |
|---|---|
| Configuration and CLI | [`config/_utility_imports.py`](../../config/_utility_imports.py) builds `icesee_kwargs` at import time; [`config_loader.py`](../../config/config_loader.py) resolves YAML inheritance |
| Application selection | [`applications/supported_models.py`](../../applications/supported_models.py) resolves a current-example hook module, then a model-family fallback |
| Runner dispatch | [`run_models_da.py`](../../src/run_model_da/run_models_da.py) validates CLI requests and maps canonical integer modes to lazily imported runners |
| Serial lifecycle | [`icesee_da_serial.py`](../../src/run_model_da/icesee_da_serial.py); useful for small cases and serial filter methods |
| Partially parallel lifecycle | [`icesee_da_partial_parallel.py`](../../src/run_model_da/icesee_da_partial_parallel.py); member/model forecast parallelism and a shared analysis transform |
| File-backed lifecycle | [`icesee_da_full_parallel.py`](../../src/run_model_da/icesee_da_full_parallel.py) and [`EnKF_parallel_io.py`](../../src/parallelization/EnKF_parallel_io.py) |
| Mode-3 dispatch and registration | [`icesee_da_distributed.py`](../../src/run_model_da/icesee_da_distributed.py), [`distributed_mode3_registry.py`](../../src/parallelization/distributed_mode3_registry.py) |
| Generic analysis | [`_mpi_analysis_functions.py`](../../src/parallelization/_mpi_analysis_functions.py), [`utils/localization.py`](../../src/utils/localization.py), [`python_enkf/`](../../src/EnKF/python_enkf/) |
| Distributed algebra and patches | [`distributed_analysis.py`](../../src/parallelization/distributed_analysis.py), [`distributed_local_analysis.py`](../../src/parallelization/distributed_local_analysis.py) |
| Native fields and members | [`distributed_fields.py`](../../src/parallelization/distributed_fields.py), [`distributed_native_runtime.py`](../../src/parallelization/distributed_native_runtime.py), [`distributed_native_cycle.py`](../../src/parallelization/distributed_native_cycle.py) |
| File-backed native analysis | [`native_file_analysis.py`](../../src/parallelization/native_file_analysis.py), member-store implementations |
| MPI planning | [`parallel_mpi/resource_plan.py`](../../src/parallelization/parallel_mpi/resource_plan.py), [`model_capabilities.py`](../../src/parallelization/parallel_mpi/model_capabilities.py), [`distributed_topology.py`](../../src/parallelization/distributed_topology.py) |
| Random streams and noise | [`random_streams.py`](../../src/utils/random_streams.py), [`_mpi_forecast_functions.py`](../../src/parallelization/_mpi_forecast_functions.py), [`native_process_noise.py`](../../src/parallelization/native_process_noise.py) |
| Physical time and restart | [`timestep_contract.py`](../../src/utils/timestep_contract.py), [`mode3_checkpointing.py`](../../src/parallelization/mode3_checkpointing.py), [`distributed_checkpoint.py`](../../src/parallelization/distributed_checkpoint.py) |
| Output interchange and views | [`postprocessing/contract.py`](../../src/postprocessing/contract.py), [`outputs.py`](../../src/postprocessing/outputs.py), [`distributed_history_export.py`](../../src/parallelization/distributed_history_export.py) |
| Diagnostics | [`utils/performance.py`](../../src/utils/performance.py), [`init_trace.py`](../../src/utils/init_trace.py), [`pyop2_compile_cache.py`](../../src/utils/pyop2_compile_cache.py) |

Files named `*_v0.py`, older duplicate I/O modules, notebooks, and comment headers are not authoritative runner selectors. Follow imports from the dispatcher and current application's entry script. The generated [API index](api-reference.md) provides the wider module inventory.

## A run from configuration to output

1. The application selects its working directory and imports the configuration loader. The loader parses named flags, reads/inherits YAML, normalizes selected values, flattens sections, applies generic overrides, and prepares the output directory. This phase has side effects; use a unique directory for a fresh run.
2. The application initializes model-specific metadata and resolves its actual `nt`, `dt`, state layout, input paths, and physical time vector. Values derived by an application can differ from loader defaults. Dispatch checks whether explicitly requested generic CLI values were replaced.
3. `icesee_model_data_assimilation` chooses the runner solely from `execution_mode`. A mode is a storage/orchestration strategy, while scientific settings remain in the same context.
4. The selected path constructs or loads truth/reference data, observations, initial members, communicators, and storage. Stable member IDs determine columns and random-stream identities.
5. Each forecast cycle advances the member from column `k` to column `k+1`. If an observation belongs to the new state column, analysis uses that state, applies declared controls, and writes/checkpoints the completed state.
6. The normal completion path closes and publishes output and emits performance information. A recoverable checkpoint and a final history view are distinct artifacts. A VDS exposes source datasets; it is not another complete copy of their values.

The context is mutable because setup, forecasts, and analysis add runtime objects and indices. It is not safe to serialize the entire dictionary blindly: it can contain MPI communicators, model modules, live server objects, arrays, and path-dependent fields. Archive YAML and a purposeful resolved metadata record.

## Execution modes and actual support

| Mode | Main representation | Appropriate use and limits |
|---|---|---|
| `0` serial | Small in-memory/member data, serial lifecycle | A useful debugging baseline; serial filter variants live here. It does not validate multi-rank behavior. |
| `1` partial parallel | Member/model forecast groups; analysis transform distributed to update work | The reference for same-input Mode-2 parity. Model memory and full-history costs can still matter. |
| `2` fully parallel bounded-memory | File-backed timestep shards; streamed observation/state rows; full or rolling retention | Avoids dense global history in RAM and reduces analysis products to ensemble size. Individual forecasts and application finalization may still need whole-member workspaces. |
| `3` application-specific native/distributed | Registered runner, native field pool or a declared whole-member boundary | PIG exposes model-native spatial ownership. Other registered applications have different representations and features; inspect their runner before production. |

`default_run`, `sequential_run`, `even_distribution`, the positional distribution selector, and `execution_flag` are scheduling controls within a runner. They do not replace `execution_mode`. `ranks_per_model` is a generic resource-planning setting; a specialized native runner may construct topology from `model_nprocs`. Inspect the application-specific startup summary to see what actually ran.

Mode-3 registration is present for `lorenz`, Icepack `idealized_pig`, Icepack `synthetic_ice_stream`, `flowline_1d`, and ISSM `ISMIP_Choi`. Registration is scoped by model and example, and import-time registration or the dispatcher bootstrap loads the implementation. It is not a promise of identical support for restart, localization, inference, error modes, or multi-rank physical state. An unsupported application raises a named `NotImplementedError`; there is no license to silently run another mode instead.

The PIG native runner uses persistent Firedrake fields and can split each model solve over spatial ranks. The Lorenz runner intentionally gathers its tiny output snapshot for compatibility. ISSM's `WholeMemberNativeAdapter` crosses the Python/MATLAB boundary with complete members, so its round-robin scheduling cannot claim the strict no-complete-member invariant. The shared `native_application.py` policy requires one spatial rank and rejects restart and scientific options not declared in `supported_features`. These boundaries matter more than the common mode number.

## State ownership

The generic `StateOwnership` helper distinguishes `replicated` and `distributed`. In a replicated model, every model rank reports the same complete state and its global size equals one local size. In a distributed model, disjoint local pieces combine and the global size is the sum of owned local sizes. Adding all replicated copies produces a wrong state dimension; selecting only one piece of a distributed field loses data. This was a real source of multi-rank failures, so an adapter must declare its semantics deliberately.

The `lorenz96` example sets `state_distribution: replicated`. The generic default remains `distributed` for backward compatibility. This flag describes the existing generic array boundary; it is not the full native-field interface. Native layouts additionally record per-variable offsets, owned intervals, and ghost indices. Native ownership validators require exact coverage with no gaps or overlaps, consistent layout identity, and stable block definitions. Ghost entries are synchronized through the model adapter and excluded from owned analysis slabs.

## Two-dimensional MPI topology

`create_distributed_topology` builds a rectangular grid with `P = E*S`, where `E` is the number of ensemble groups and `S` the spatial ranks per group. Rank order is ensemble-major: `ensemble_slot, spatial_rank = divmod(world_rank, S)`. `spatial_comm` connects ranks within one member solve; `ensemble_comm` connects the same spatial coordinate across groups. Stable member IDs are scheduled round-robin across ensemble slots, so a slot can process multiple members.

For example, `P=8`, `S=2` gives four groups of two ranks. With 12 members, each group processes three members. An analysis block exchanges corresponding owned rows across ensemble groups, while the native forecast uses spatial halo exchange within its group. The grid must be rectangular and ownership metadata must agree before large state data move. Do not substitute world collectives for one of these dimensions: a seemingly harmless world barrier inside a member callback can deadlock when groups process different members.

## Storage, memory, and the limits of bounded execution

For state dimension `n`, ensemble size `N`, saved columns `T+1`, and stored item size `b`, a full ensemble history is approximately `n*N*(T+1)*b` bytes before metadata, compression, temporary copies, means, truth, checkpoints, and diagnostics. Rolling retention removes the growth with run length; it cannot remove the conventional EnKF snapshot cost `n*N*b`. A row block with `r` rows needs multiple `r*N*8` arithmetic buffers even if on-disk storage uses float32.

Mode 2 streams observations into ensemble-sized products, then reads state rows to apply the transform. HDF5 windows must be opened/extended by all required communicator ranks in matching order. Subgroup-root writes cannot secretly create a collective window. Analysis reads use immutable forecasts; private serial shards and bounded assembly avoid mixing reader and writer handles unsafely. Full-history VDS finalization is metadata over shards; deleting or moving those shards breaks it. Rolling pruning happens only after the retained restart state is committed.

Mode 3 avoids gathering physical PIG member fields during native forecast/analysis, but not every auxiliary object is fully spatially distributed. The current PIG process-noise bridge generates a canonical global increment on a model-group leader and scatters owned entries. It has an explicit 256 MiB stochastic-workspace guard and an `O(n)` workspace; it does not gather physical member state. Reducing analysis row chunks will not remove this lower bound. Treat an exceeded guard as a design limit requiring another noise policy, not as permission to increase memory blindly.

## Failure boundaries

Healthy collectives require all communicator participants to enter in a consistent order. If one rank raises unexpectedly, a collective cleanup barrier is unsafe because its peers may be blocked elsewhere. The full-parallel error convention prints the original traceback and rank, performs only best-effort rank-local checkpoint work, and calls `Abort(1)`. Do not add collective close/finalization in that handler. Reconstruct a recovery view explicitly after the job has stopped and the files have been inspected.

In CI, success requires more than exit zero: completion markers, expected columns/shapes, finite values, and member-wise parity are independently checked. In production, scheduler accounting, MPI-step exit, verifier exit, and scientific acceptance remain separate records.
