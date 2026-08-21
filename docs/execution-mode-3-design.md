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

## Stochastic patch-local analysis

The first mode-3 analysis implementation will be the existing grouped-local
stochastic EnKF expressed over distributed state and observation patches. Each
patch reads only owned state rows, its localized observation group, and required
halo metadata. Observation perturbations remain deterministically keyed by
cycle, observation, and ensemble-member identifiers. Collective operations are
restricted to the relevant ensemble communicator and operate on patch matrices
or `Nens x Nens` factors rather than complete members.

This is a distribution of the existing stochastic analysis, not LETKF or EnSRF.
Any alternative square-root or deterministic update must use a separate
`filter_type` and a separate scientific benchmark.

## Checkpointing and storage

Parallel HDF5 is the first checkpoint backend because ICESEE already uses HDF5
metadata and datasets. ADIOS2 may be added as an optional backend for platforms
where it provides better collective or asynchronous I/O. Checkpoints contain
distributed state slabs, ownership metadata, member IDs, random-stream state,
the completed cycle, observation metadata, inversion state, and sufficient
information for atomic recovery. No backend may require a rank-0 full-member
gather.

Hybrid MPI plus controlled threaded kernels or BLAS is permitted within each
spatial rank. Thread counts must be explicit to prevent oversubscription. Dask
may support preprocessing and workflow scheduling, but is not part of the MPI
collective analysis path.

## Phased implementation

1. **Adapter API and topology:** add non-selectable interfaces, ownership tests,
   and communicator construction.
2. **State-only forecast parity:** implement Lorenz and one distributed model
   adapter; compare reconstructed results against mode 2.
3. **Distributed global stochastic analysis:** preserve exact analysis controls
   before optimizing patch locality.
4. **Grouped-local stochastic analysis:** distribute observation groups and
   patch updates without full-member reconstruction.
5. **Parallel checkpoint/restart:** add rank-placement-independent restart and
   failure cleanup.
6. **ISSM hybrid gate:** validate geometry finalization, delayed inversion,
   member-wise handoff, and the forecast following inversion.
7. **Scale gate:** demonstrate bounded per-rank memory and I/O for a state that
   cannot fit as one complete member on any rank.

Only after phases 1--3 pass should `normalize_execution_mode` accept mode 3.

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
