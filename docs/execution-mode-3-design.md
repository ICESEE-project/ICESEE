# Execution mode 3: spatially distributed design

## Status and purpose

`execution_mode: 3` is a proposed large-scale execution architecture. It is not
yet accepted by the runtime configuration. Mode 3 will distribute both ensemble
members and each member's spatial state so that no rank must hold or reconstruct
a complete high-dimensional member.

Execution mode and filter choice are orthogonal. The default mode-3 analysis
must preserve ICESEE's existing stochastic EnKF equations, observation-factor
mode, grouped local analysis, localization, inflation, constraints, and hybrid
inversion schedule. LETKF, EnSRF, reduced-order, and multilevel methods may be
added later only as explicitly named filter algorithms with their own scientific
validation; selecting mode 3 must not silently select one of them.

## Process topology

Mode 3 uses a two-dimensional MPI topology:

- an **ensemble communicator** joins ranks owning the same spatial slab for
  different ensemble members; and
- a **spatial communicator** joins the subdomains belonging to one ensemble
  member.

For `P_e` ensemble groups and `P_s` spatial ranks per member, the useful process
grid is `P_e x P_s`. Multiple members may be scheduled sequentially within an
ensemble group when `Nens > P_e`; member identity, not rank identity, continues
to key all random streams.

## Distributed model-adapter contract

A mode-3-capable model adapter must provide:

1. a stable global state layout and ownership map;
2. local state slabs plus global offsets or indices;
3. halo/ghost exchange needed by the forecast model;
4. a local observation operator and ownership of observation rows;
5. distributed geometry projection and analysis-finalization hooks;
6. serialization metadata sufficient to restart with a different compatible
   rank placement; and
7. for hybrid workflows, an inversion handoff that does not gather a complete
   member on rank 0.

Adapters that cannot expose distributed state remain supported by modes 0--2,
but cannot claim mode-3 scalability.

The initial API is implemented in
`src/parallelization/distributed_adapter.py`. A single-vector adapter may
describe one contiguous owned interval. A multidimensional model instead uses
a segmented block layout: every named variable has its own global size and
contiguous PETSc-style owned interval, while the global state retains ICESEE's
existing variable-major ordering. Local arrays concatenate only the owned
pieces and never allocate the gaps between them. ICESEE validates collectively
that every variable block is covered exactly once, using metadata proportional
to the number of variables and spatial ranks rather than to the global state
size.

`src/parallelization/distributed_fields.py` defines the model-agnostic field
contract and registry that maps those owned pieces to native model vectors.
`src/parallelization/distributed_native_runtime.py` retains those vectors and
their model context for the lifetime of each locally scheduled member. It
creates compact owned snapshots only at analysis or checkpoint boundaries and
writes analyzed slabs back into the persistent native vectors before the next
forecast. `src/parallelization/distributed_native_cycle.py` composes that
lifecycle with bounded observation batches, the existing stochastic analysis,
transactional member-wise inversion, topology-portable restore, and
rank-sharded checkpoints. Observation batches are fingerprinted across the
ensemble communicator before analysis collectives so a disagreement in row
identity or ordering fails clearly instead of hanging MPI.
The first native bridge is
`applications/icepack_model/icepack_utils/_distributed_fields.py`: it reads
PETSc ownership from Firedrake scalar spaces, treats velocity components as
separate state blocks, writes only owned degrees of freedom, and performs one
halo refresh per backing Firedrake `Dat`. This bridge establishes distributed
storage semantics and now provides a persistent member context plus callback
wrappers for native initialization, forecast, observation, and post-analysis
finalization. Applications may additionally supply inversion and checkpoint-
restore callbacks. Those optional capabilities are exposed only when the
application configures them, so state-only adapters cannot accidentally claim
hybrid or restart support. Wiring the production PIG physics callbacks and
passing its forecast/restart/scale gates remain necessary before that
application can select mode 3.

Hybrid adapters expose a distributed inversion callback; the inversion
requirement is not imposed on state-only models. ICESEE snapshots locally owned
fields before inversion and restores them collectively when any rank reports a
failure. Model adapters must similarly stage auxiliary solver side effects
until collective success, or provide their own rollback. An optional
localization adapter contract describes locally owned analysis targets and
observation coordinates in one or more physical dimensions. It also declares
which observation kinds affect each target variable. This keeps mesh geometry
and variable coupling in the model adapter while the stochastic analysis
remains model agnostic.

Canonical observation IDs are unique across observation kinds. For each target
block, all active accepted kinds are concatenated into one bounded local
observation collection and enter one stochastic transform. They are not
applied as sequential per-kind filters. This is required for mode-1 scientific
parity when, for example, thickness and velocity observations both update the
same geometry block.

## Stochastic patch-local analysis

The first mode-3 analysis implementation preserves the existing stochastic
EnKF. The implemented global-analysis primitives stream local observation rows
into fixed-size `Nens x Nens` products and apply the resulting ICESEE ensemble
transform to bounded local-state row blocks. At most
`row_chunk_size x Nens` state values are assembled across ranks owning the same
spatial slab; no complete member is reconstructed. The implementation preserves
the current `stochastic_R`, `generated_R`, and `legacy_prior_anomalies` factors.

Grouped-local analysis now extends the same mechanism to observation groups.
Each observation owner retains a local coordinate tree. Target owners send
bounded coordinate chunks only to ranks whose observation bounding boxes
intersect the localization radius, and receive only matching canonical
observation IDs and stochastic-analysis rows. Stable observation IDs make
patch grouping and SVD ordering independent of message arrival or rank
placement. Only the ensemble-slot-zero spatial communicator performs these
sparse queries; the resulting `Nens x Nens` transforms can then be shared
across the orthogonal ensemble communicator. No global coordinate tree,
observation table, patch registry, or complete member is formed. A real
four-rank parity gate matches the existing grouped-local stochastic analysis
to machine precision (`6.22e-15` maximum absolute error).

This is a distribution of the existing stochastic analysis, not LETKF or EnSRF.
Any alternative square-root or deterministic update must use a separate
`filter_type` and a separate scientific benchmark.

## Checkpointing and storage

The first checkpoint backend uses rank-sharded serial HDF5, which works even
when `h5py` was not built with MPI support. Each shard contains only locally
owned slabs and is indexed by stable member IDs and global block intervals. It
preserves the state dtype, supports spatial ranks that own zero rows of a small
block, and does not repack segmented state into a complete vector. The
manifest is committed last, so an interrupted write cannot be mistaken for a
complete restart. A compatible restart may change ensemble or spatial rank
placement and reads only overlaps needed by its new local ownership. Parallel
HDF5 and ADIOS2 remain planned optional backends for platforms where collective
or asynchronous I/O is advantageous. Checkpoint sidecar metadata carries
random-stream keys, the completed cycle, observation metadata, inversion state,
and other model-agnostic restart state. No backend may require a rank-0
full-member gather.

Hybrid MPI plus controlled threaded kernels or BLAS is permitted within each
spatial rank. Thread counts must be explicit to prevent oversubscription. Dask
may support preprocessing and workflow scheduling, but is not part of the MPI
collective analysis path.

## Phased implementation

1. **Adapter API and topology:** add non-selectable interfaces, ownership tests,
   and communicator construction. Contiguous and segmented block contracts,
   native-field packing, and topology construction are implemented; no runtime
   mode is enabled by them.
2. **State-only forecast parity:** deterministic round-robin member scheduling,
   local-slab initialization and forecast validation, and a parity-only root
   reconstruction helper are implemented in
   `src/parallelization/distributed_runtime.py`. A deliberately small Lorenz
   reference adapter now exercises spatial ownership and distributed RK4
   forecast parity. Its complete-state allgather is guarded and is permitted
   only for this three-variable reference; large adapters must use native halo
   exchange. A real-MPI state-only parity gate is provided by
   `scripts/benchmarks/run_mode3_lorenz_forecast_parity.py`; it reconstructs
   complete states only through the parity-only helper and compares them with
   the whole-member RK4 recurrence. The model-agnostic persistent native cycle
   and Icepack/Firedrake callback bridge are implemented. The next step is to
   connect the production Icepack physics callbacks, remove its whole-array
   forecast path, and compare reconstructed histories against mode 2.
3. **Distributed global stochastic analysis:** bounded-memory observation
   products and local-state ensemble transforms are implemented. A real four-
   rank Lorenz gate (two ensemble groups by two spatial ranks) passes forecast
   plus stochastic-analysis parity exactly.
4. **Grouped-local stochastic analysis:** multidimensional adapter metadata,
   sparse distributed coordinate queries, canonical observation grouping, and
   bounded target-row updates are implemented. The four-rank numerical parity
   gate is `scripts/benchmarks/run_mode3_local_analysis_parity.py`. Integration
   across the full two-dimensional process grid is exercised by
   `scripts/benchmarks/run_mode3_local_runtime_parity.py`; the initial
   two-ensemble-group by two-spatial-rank gate matches the existing grouped-
   local stochastic update to machine precision (`3.44e-15` maximum absolute
   error). Integration into a selectable runner and a distributed ice-model
   adapter remain.
5. **Parallel checkpoint/restart:** the first rank-sharded, manifest-last HDF5
   backend and rank-placement-independent loader are implemented. A four-rank
   gate now writes on a `2 x 2` process grid, restores exactly on a `1 x 4`
   grid, rejects an injected shard-write failure collectively, and retains the
   preceding complete checkpoint. An interrupted/restarted stochastic
   trajectory also matches the uninterrupted trajectory exactly. A segmented
   multi-variable topology-portability gate exercises the block checkpoint
   format separately.
6. **Native restart and hybrid lifecycle:** topology-portable restore directly
   into persistent fields and transactional member-wise inversion are
   implemented and unit tested. The ISSM hybrid gate must still validate
   geometry finalization, delayed inversion, member-wise handoff, and the
   forecast following inversion using model-native distributed state.
7. **Scale gate:** `scripts/benchmarks/run_mode3_memory_scaling.py` analytically
   preflights a state too large to allocate on the launch machine. It reports
   compact snapshots, bounded state/observation chunks, ensemble-space
   products, and one minimum resident native-state copy separately. The budget
   intentionally excludes additional nonlinear-solver vectors, MPI/HDF5
   buffers, Python overhead, and filesystem cache; every production adapter
   must add a measured application-specific peak before mode 3 is enabled.

The immediate production path is therefore: retain persistent model-native
distributed fields for each locally scheduled member; pack only owned pieces
for analysis; call the model with native halo exchange; evaluate observations
from local fields; and checkpoint owned member/block shards. The required
memory invariant is that no rank allocates a complete member, the complete
ensemble, or the complete time history. Analysis work arrays must remain
bounded by configured state-row and observation-row chunk sizes.

Mode 3 remains non-selectable until the state-only ice-model, restart, hybrid
inversion, and scale gates pass. Passing the small algebraic gates alone is not
sufficient to expose it as a production execution mode.

For example, a state with approximately 37 GiB per member can be preflighted
without allocating it:

```bash
python scripts/benchmarks/run_mode3_memory_scaling.py \
  --global-rows 4966055936 --nens 40 \
  --ensemble-groups 8 --spatial-ranks 16 \
  --state-row-chunk 4096 --observation-row-chunk 4096 \
  --max-rank-gib 40
```

Increasing `spatial-ranks` reduces the owned state and native-state terms;
reducing either chunk size bounds temporary analysis storage. Increasing
`ensemble-groups` reduces the number of persistent members scheduled on each
process-grid column. The process count is
`ensemble-groups * spatial-ranks`, and the model itself must support the chosen
spatial communicator.

## Acceptance gates

Mode 3 must pass all mode-2 scientific controls plus the following:

- complete-member, complete-history parity against mode 2 for Lorenz;
- complete-state parity for a spatially decomposed application;
- grouped-local analysis and localization parity by global state index;
- fresh versus interrupted/resumed equivalence;
- ISSM state-only parity followed by a shortened hybrid inversion parity run;
- rank-layout invariance for at least two compatible `P_e x P_s` topologies;
- deterministic behavior for member-keyed and observation-keyed random streams;
- atomic checkpoint recovery after an injected writer or rank failure; and
- measured per-rank peak memory that scales with the local slab, not the global
  member or full ensemble history.

Scientific parity uses the same default numerical gate as mode 2
(`atol=1e-10`, `rtol=1e-8`) unless a storage backend declares and justifies a
different precision tolerance.

## Bounded-memory member lifecycle (implemented, 2026-09-26/27)

The sections above describe mode 3's original design target. This section
documents what is actually implemented and real-Firedrake-verified as of
2026-09-27, for the Icepack (idealized_pig) application specifically.

**Lifecycle**: SHARED MODEL CONTEXT (mesh, function spaces, solver, static
forcing fields -- `_SharedPigContext`, cached once per ensemble group,
never duplicated by rounds) → ACTIVE NATIVE MEMBER (at most one member's
live Firedrake state resident at a time) → INACTIVE MEMBER STORE (a
packed owned-state array, backend-neutral) → DURABLE CHECKPOINT (a
completely separate, unrelated mechanism -- `distributed_checkpoint.py`).

**Two `InactiveMemberStore` backends** (`src/parallelization/
distributed_member_store.py`, `distributed_member_store_hdf5.py`):
`MemoryInactiveMemberStore` (default, unchanged behavior) and
`HDF5MemberMajorStore` (opt-in, one HDF5 file per rank regardless of
ensemble size, dataset-handle-cached, real-Firedrake-verified bit-for-bit
equivalent to the memory backend across P=1, genuine spatial
decomposition, and real rounds scheduling). Selected via
`build_inactive_member_store(icesee_kwargs, world_rank=...)` -- the only
place the generic runtime imports a concrete backend; applications never
import a backend directly. Config: `member_store_backend` (`"memory"` |
`"hdf5_member_major"`, default `"memory"`), `member_store_root` (default
`<data_path>/_mode3_member_store/<run_id>/`).

**Store-streaming analysis** (`run_native_store_streaming_analysis_cycle`,
a sibling to the original `run_native_global_analysis_cycle` -- the
original is kept, unmodified, as the reference implementation): the
ensemble-space transform (`X_I^a = X_I^f @ T`) is applied one state-row
block at a time, pulling/pushing through the store's `get_ensemble_rows`/
`put_ensemble_rows` bulk primitives instead of requiring every
round-assigned member's full array simultaneously resident. **Important
optimization**: this pass is skipped entirely on any timestep with no
scheduled analysis event, since the ensemble transform is then provably
(verified bit-for-bit) exactly the identity matrix -- a real Icepack run
showed this cuts HDF5 operations by 12.5x for a typical sparse-observation
schedule (1675 -> 134 bulk calls, formula validated in
`src/parallelization/mode3_projections.py`).

**Temporary-store lifecycle**: created lazily on pool construction; each
rank cleans up only its own file, after the run's existing final barrier
(no new synchronization). On abnormal termination (exception before that
point), the store file is deliberately left on disk for debugging --
**restart must use the last durable checkpoint, never the temporary
inactive-state store**, which is not a recovery artifact.

**What this does NOT yet claim**: real-scale (Ne~1000, ~37GB/member)
performance or memory savings -- only local, small-Ne correctness and a
measured-and-validated set of operation-count/traffic formulas used to
project (not claim) larger-scale behavior. See
`src/parallelization/mode3_projections.py` for the formulas and
`scripts/benchmarks/mode3_large_state_benchmark.py` for the synthetic
benchmark harness. PACE validation is required before any scalability
claim beyond what is stated here.
