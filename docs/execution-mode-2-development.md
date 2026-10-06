# Execution mode 2: parity and development contract

This document is the engineering contract for ICESEE's fully parallel engine
(`execution_mode: 2`). Mode 2 is not a separate filter. It must preserve the
scientific semantics of partial parallel mode (`execution_mode: 1`) while
changing data placement, I/O, and reduction strategy so that large ensembles
do not require a dense in-memory history.

## Non-negotiable parity invariants

For one configuration, seed, observation file, initial ensemble, and set of
member identifiers, modes 1 and 2 must use the same:

1. observation times, active observation rows, values, and uncertainties;
2. member-keyed initialization and process-noise streams;
3. member forecast inputs and model adapter;
4. observation-error/analysis-factor mode;
5. EnKF controls, including localization, grouped local analysis, inflation,
   relaxation, frozen variables, bounds, and increment caps;
6. geometry projection and application-specific analysis finalization;
7. delayed inversion schedule and member-wise inversion handoff; and
8. restart state, timestep, and output-time interpretation.

The only accepted differences are floating-point roundoff caused by a different
MPI reduction order, or a separately declared storage-precision tolerance.
Comparing only ensemble means or plotted RMSE is insufficient: acceptance uses
every member value at every retained synchronized timestep.

## Architecture

Mode 2 currently follows five bounded-memory phases:

1. **Forecast by member.** Member IDs, rather than rank order, key all random
   streams. The shared member-advance helper is also used by mode 1.
2. **Timestep shard publication.** A complete restartable ensemble snapshot is
   stored as a state-by-member HDF5 shard. Publication is atomic.
3. **Streamed observation-space analysis.** State and observation rows are read
   in bounded chunks. MPI reductions operate on ensemble-space Gram,
   cross-product, and right-hand-side matrices.
4. **Streamed state update and finalization.** The transform is applied in row
   or member slabs. Application constraints use the same finalizer as mode 1.
   Hybrid inversion is completed member by member before the analysis shard is
   atomically published.
5. **Checkpoint/history policy.** Full mode retains shards and can expose them
   through an HDF5 virtual dataset. Rolling mode checkpoints the matching shard
   before pruning older member history.

## Capability and validation matrix

"Implemented" does not by itself mean that a target application has passed its
paired integration benchmark. The last column records the strongest current
gate.

| Capability | Mode-2 implementation | Current validation |
|---|---|---|
| Deterministic member initialization and forecast | member-keyed RNG and shared forecast helper | paired Lorenz complete-history PASS |
| Standard stochastic observation errors | compact seeded observation-sized factor | paired Lorenz complete-history PASS |
| Generated spatial observation errors | deterministic member columns, streamed to temporary HDF5 | exact unit parity and determinism |
| Legacy forecast-anomaly factor | streamed ensemble-space products | algebra/unit tests |
| Global analysis and analysis controls | row-streamed transform and common controls | unit tests; paired Lorenz PASS |
| Sparse observations and bed masks | compact row-index representation | unit/integration tests |
| Grouped local analysis/localization | observation terms distributed; patch factors currently assembled on root within a hard budget | unit tested; large distributed patch redesign pending |
| ISSM geometry finalization | bounded member-slab common finalizer | exact mode-1/mode-2 unit parity |
| Delayed hybrid inversion | frozen EnKF friction block, member-wise inversion, atomic publish | paired shortened ISSM hybrid complete-history PASS |
| Restart and explicit restart timestep | shard/checkpoint synchronization and cleanup of newer state | implementation/unit tests; paired interrupted-vs-fresh run required |
| Full history | sharded files plus optional VDS/materialization | unit/integration tests |
| Rolling history | checkpoint-before-prune, retained mean history | unit tests |

The paired ISSM inversion gate has passed. Interrupted/resumed ISSM equivalence
and longer survey-campaign runs remain separate release gates; they are not
reasons to weaken the numerical contract.

### Verified ISSM hybrid parity gate

On 2026-08-21, a fresh paired ISMIP-Choi run (paired run ID
`28340b56-9e18-42cd-beb7-87eee9436b5a`) passed the complete-member comparison:

- 12 synchronized timesteps;
- ensemble shape `20082 x 4`;
- 963,936 compared values and zero mismatches;
- maximum absolute difference `1.17133822641e-08`; and
- acceptance tolerances `atol=1e-10`, `rtol=1e-8`.

The shortened run includes the first joint state/bed analysis and hybrid
inversion at year 2, followed by the next forecast. It therefore validates the
mode-2 reduced-state finalization, frozen-friction handoff, member-wise ISSM
inversion, publication, and post-inversion forecast against mode 1. It does not
by itself validate interrupted restart equivalence or later observation
campaigns.

## Acceptance suite

Run focused bounded-memory and semantics tests:

```bash
pytest -q \
  src/tests/test_full_parallel_large_data.py \
  src/tests/test_random_field_coordinates.py
```

Run both engines from a fresh common configuration and compare all members:

```bash
python scripts/benchmarks/run_execution_mode_parity.py \
  --config applications/lorenz_model/examples/lorenz96/params.yaml \
  --application applications/lorenz_model/examples/lorenz96/run_da_lorenz96.py \
  --output /tmp/icesee-mode-parity \
  --ranks 2 --nens 8
```

The default gate is `atol=1e-10`, `rtol=1e-8`. Use `--bitwise` only to diagnose
one fixed MPI/platform combination. For every new model adapter or hybrid
feature, add a shortened paired run that includes at least one observation
cycle, one forecast after analysis, and—when enabled—one inversion cycle.

Application-specific command-line arguments can be forwarded by repeating
`--application-arg`. For example, a shortened ISSM hybrid benchmark can be run
with:

```bash
/Users/bkyanjo3/venv-firedrake/bin/python \
  scripts/benchmarks/run_execution_mode_parity.py \
  --config applications/issm_model/examples/ISMIP_Choi/params.yaml \
  --application applications/issm_model/examples/ISMIP_Choi/run_da_issm.py \
  --output /private/tmp/icesee-issm-mode-parity \
  --mpiexec /opt/homebrew/bin/mpirun \
  --python /Users/bkyanjo3/venv-firedrake/bin/python \
  --ranks 2 --nens 4 --num-years 2.4 \
  --application-arg=--model_nprocs=1
```

With the checked-in ISMIP-Choi schedule, year 2 is both the first state/bed
analysis and the first enabled inversion, and 2.4 years retains the following
forecast. This is the inexpensive regression gate for the snapshot-metadata
and inversion-handoff paths. A longer application benchmark is still required
before claiming parity for later survey campaigns, forecasts, and restart
behavior. In general, an ISSM parity configuration must be short enough for
routine testing but long enough to include at least one observation analysis,
one successful inversion, and the forecast immediately after that inversion.
A state-only comparison is not sufficient evidence for hybrid parity.

Restart equivalence is a separate three-run test: create a fresh mode-1
reference, a fresh uninterrupted mode-2 run, and an interrupted/resumed mode-2
run. Both mode-2 histories after the restart point must satisfy the same parity
gate against the mode-1 reference.

## Scaling limits and next development stage

The current member-parallel, streamed ensemble-space design is a strong generic
baseline when one member fits on its assigned model ranks and `Nens x Nens`
analysis matrices are small. It removes the run-length multiplier from rolling
storage, but cannot remove the physical size of one ensemble snapshot. A 37-GB
member and 40 members imply about 1.48 TB per uncompressed snapshot.

With the first ISSM hybrid parity gate passed, the preferred scaling upgrade is
tracked as proposed `execution_mode: 3`: a two-dimensional MPI topology with
one communicator dimension over ensemble members and one over model spatial
subdomains. Mode 3 is an execution architecture, not a new filter, and is not
selectable until its adapter and parity contracts are implemented. Model
adapters should expose distributed state slabs and local observation operators
so no rank reconstructs a whole member. The analysis should preserve the same
grouped-local stochastic EnKF update on distributed state and observation
slabs, with member-keyed observation perturbations and the same analysis
equations as mode 1. Parallel HDF5 or ADIOS2 checkpointing can remove rank-0 I/O
bottlenecks, while a two-dimensional MPI topology can distribute both ensemble
members and model spatial subdomains.
Hybrid MPI plus threaded model kernels or BLAS is appropriate when thread counts
are controlled to avoid oversubscription. Dask can help with preprocessing,
postprocessing, and coarse workflow scheduling, but it should not be inserted
into the MPI collective analysis hot path because doing so can duplicate data
and introduce a second scheduler.

LETKF, EnSRF, reduced-order, and multilevel methods remain useful optional
research directions, but they are different filters or approximations. They
must be exposed under distinct algorithm names and validated scientifically;
they must not replace the stochastic EnKF merely because execution mode 2 is
selected.

The proposed architecture, adapter interface, rollout phases, and acceptance
gates are specified in
[`execution-mode-3-design.md`](execution-mode-3-design.md).
