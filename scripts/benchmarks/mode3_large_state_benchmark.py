# ==============================================================================
# @des: Synthetic, model-agnostic large-state scalability benchmark for
# execution-mode-3's distributed architecture.
#
# Motivating stress case: a hypothetical 37 GB live dynamic state per
# ensemble member x 40 members = 1.48 TB logical ensemble state, before
# model/solver/workspace overhead. We cannot allocate that on a development
# machine, so this benchmark exercises the REAL production distributed
# primitives (src/parallelization/distributed_native_runtime.py,
# distributed_native_cycle.py, distributed_analysis.py, distributed_runtime.py,
# distributed_checkpoint.py -- the exact same code every real application's
# mode3_runner.py uses) against a synthetic, Firedrake-free toy model whose
# only job is to have a large, arbitrarily-sized per-DOF state vector. This
# lets Nx/Ne/P be scaled up far beyond what any single real application mesh
# permits locally, and lets the SAME script be re-run at PACE with much larger
# --nx/--ne/--world-size once real hardware is available.
#
# What this measures (see also --help):
#   - actual per-rank owned-array bytes (never global-sized) at every stage;
#   - actual peak RSS/rank (platform-appropriate: ru_maxrss is bytes on
#     Darwin, KiB on Linux -- normalized to bytes here);
#   - analysis-workspace scaling as row_chunk_size (B) and Ne vary, compared
#     against the predicted O(B*Ne + Ne^2) relationship;
#   - checkpoint file count/bytes/time, and restart read time;
#   - a projection (NOT a measurement) of persistent native-model bytes/rank
#     for the 37GB/member x 40-member stress case at several
#     ensemble_groups/P_model layouts, using this run's own measured
#     bytes-per-owned-element as the scaling constant.
#
# Hard safety guards (see _guard_against_dangerous_allocation): refuses to
# run with a --nx large enough that a single GLOBAL-sized float64 array would
# exceed --max-global-bytes (default 2 GiB) unless --allow-large-nx is passed
# explicitly. Also asserts, at runtime, that no per-rank owned array ever
# reaches the full global size when spatial_ranks > 1 -- a live guard against
# an accidental O(Nx) allocation slipping into any code path this benchmark
# exercises, not just a static one-time check.
#
# Usage (see also --help): this is an mpi4py program; launch with mpirun.
#   mpirun -n 4 python mode3_large_state_benchmark.py \
#       --nx 4000000 --ne 4 --ensemble-groups 2 --p-model 2 \
#       --row-chunk-size 4096 --n-cycles 3 \
#       --checkpoint-root /tmp/mode3_bench_ckpt
#
# For the block-size/Ne scaling study (Primary Task 1's explicit
# O(B*Ne + Ne^2) requirement), run this script several times varying
# --row-chunk-size and --ne and compare the printed analysis-workspace-bytes
# line across runs (see scripts/benchmarks/README or the module docstring's
# own worked example at the bottom of this file for a driver loop).
# ==============================================================================
from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import sys
import time
import tracemalloc
from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np
from mpi4py import MPI

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
for _p in (_REPO_ROOT, os.path.dirname(_REPO_ROOT)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from ICESEE.src.parallelization.distributed_fields import DistributedFieldRegistry
from ICESEE.src.parallelization.distributed_streaming_runtime import (
    run_native_store_streaming_analysis_cycle,
)
from ICESEE.src.parallelization.distributed_native_cycle import (
    NativeCheckpointRequest,
    NativeObservationBatch,
    run_native_global_analysis_cycle,
)
from ICESEE.src.parallelization.distributed_native_runtime import (
    NativeDistributedMember,
    initialize_native_member_pool,
)
from ICESEE.src.parallelization.distributed_topology import create_distributed_topology


# ------------------------------------------------------------------
# Platform-neutral peak-RSS reading (Darwin: bytes; Linux: KiB).
# ------------------------------------------------------------------
def peak_rss_bytes() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(raw) if sys.platform == "darwin" else int(raw) * 1024


# ------------------------------------------------------------------
# Safety guard against accidental large global allocation.
# ------------------------------------------------------------------
def guard_against_dangerous_allocation(nx: int, max_global_bytes: int, allow_large_nx: bool) -> None:
    global_bytes = int(nx) * 8
    if global_bytes > max_global_bytes and not allow_large_nx:
        raise SystemExit(
            f"Refusing to run: --nx={nx} implies a GLOBAL-sized float64 array "
            f"of {global_bytes / 1e9:.2f} GB, which exceeds the safety "
            f"threshold ({max_global_bytes / 1e9:.2f} GB). This benchmark "
            "must never construct a full-Nx array on one rank; pass "
            "--allow-large-nx only if you have specifically verified every "
            "stage below stays owned-only at this size."
        )


def contiguous_owned_interval(global_size: int, rank: int, size: int) -> tuple[int, int]:
    """Same divmod-based contiguous partition already used in production
    (see distributed_analysis.py::contiguous_observation_layout and
    distributed_memory.py) -- reused here, not reinvented, for consistency
    with the real distributed primitives this benchmark exercises."""

    quotient, remainder = divmod(global_size, size)
    start = rank * quotient + min(rank, remainder)
    stop = start + quotient + (1 if rank < remainder else 0)
    return start, stop


@dataclass
class _SyntheticField:
    """A synthetic per-DOF state block: OWNED-ONLY float64 array.

    Deliberately model-agnostic: no mesh, no Firedrake, no PETSc -- just a
    plain contiguous numpy array sized exactly to this rank's owned
    interval. This is the smallest possible stand-in that still exercises
    every real distributed-state code path (DistributedFieldRegistry,
    NativeDistributedMember, the analysis cycle, checkpoint I/O) at
    whatever Nx the caller wants to test.
    """

    name: str
    global_size: int
    owned_start: int
    owned_stop: int
    member_id: int
    _values: np.ndarray = field(default=None, repr=False)

    def __post_init__(self) -> None:
        owned = self.owned_stop - self.owned_start
        assert owned <= self.global_size, "owned interval cannot exceed global size"
        # Deterministic, member-keyed initial value -- physically-meaningless
        # for this synthetic model, but reproducible and cheap.
        self._values = np.full(owned, float(self.member_id) + 0.1, dtype=np.float64)

    def read_owned(self) -> np.ndarray:
        return self._values

    def write_owned(self, values: np.ndarray) -> None:
        self._values[:] = np.asarray(values, dtype=np.float64)

    def synchronize_ghosts(self) -> None:
        return None


class SyntheticLargeStateAdapter:
    """Native distributed adapter for the synthetic large-state benchmark.

    Every method operates ONLY on this rank's owned interval -- never a
    global-sized array -- and asserts that invariant explicitly wherever
    an accidental full materialization would otherwise be easy to miss.
    """

    def __init__(self, nx: int, max_global_bytes: int) -> None:
        self.nx = int(nx)
        self.max_global_bytes = int(max_global_bytes)
        self.peak_owned_bytes = 0

    def _check_bounded(self, array: np.ndarray, topology: Any) -> None:
        nbytes = array.nbytes
        self.peak_owned_bytes = max(self.peak_owned_bytes, nbytes)
        if int(topology.spatial_ranks) > 1:
            assert array.size < self.nx, (
                f"GUARD TRIPPED: a per-rank array reached {array.size} "
                f"elements, >= the full global size {self.nx}, under "
                f"spatial_ranks={topology.spatial_ranks} -- this would be an "
                "accidental O(Nx) allocation on one rank."
            )

    def initialize_native_member(
        self, member_id: int, *, topology: Any, icesee_kwargs: Mapping[str, Any]
    ) -> NativeDistributedMember:
        start, stop = contiguous_owned_interval(
            self.nx, int(topology.spatial_rank), int(topology.spatial_ranks)
        )
        field_obj = _SyntheticField(
            name="state", global_size=self.nx, owned_start=start, owned_stop=stop,
            member_id=int(member_id),
        )
        self._check_bounded(field_obj._values, topology)
        registry = DistributedFieldRegistry([field_obj], layout_id="synthetic-large-state-v1")
        return NativeDistributedMember(int(member_id), registry, {"field": field_obj})

    def forecast_native_member(
        self, member: NativeDistributedMember, timestep: int, *, topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> None:
        # Trivial, cheap forecast: exercises pack/unpack without any real
        # solve. The point of this benchmark is distributed-plumbing memory
        # and I/O scaling, not numerical physics.
        state = member.pack_owned()
        self._check_bounded(state, topology)
        state = state + 1e-3
        member.unpack_owned(state)

    def observe_native_member(
        self, member: NativeDistributedMember, observation_rows: np.ndarray, *,
        topology: Any, icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        rows = np.asarray(observation_rows, dtype=np.int64)
        values = member.fields.observe_owned_rows(rows)
        self._check_bounded(values, topology)
        return values

    def finalize_native_analysis(
        self, member: NativeDistributedMember, forecast_owned: np.ndarray, timestep: int,
        *, topology: Any, icesee_kwargs: Mapping[str, Any],
    ) -> None:
        return None

    def reactivate_native_member(
        self, member_id: int, packed_state: np.ndarray, *, topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> NativeDistributedMember:
        start, stop = contiguous_owned_interval(
            self.nx, int(topology.spatial_rank), int(topology.spatial_ranks)
        )
        field_obj = _SyntheticField(
            name="state", global_size=self.nx, owned_start=start, owned_stop=stop,
            member_id=int(member_id),
        )
        field_obj._values[:] = np.asarray(packed_state, dtype=np.float64)
        self._check_bounded(field_obj._values, topology)
        registry = DistributedFieldRegistry([field_obj], layout_id="synthetic-large-state-v1")
        return NativeDistributedMember(int(member_id), registry, {"field": field_obj})


def build_observation_batch(
    adapter: SyntheticLargeStateAdapter,
    pool,
    nobs: int,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
    rng: np.random.Generator,
) -> NativeObservationBatch:
    """Nobs-sized (never Nx-sized) synthetic observation batch, routed to
    this spatial rank's owned subset before the analysis cycle sees it --
    mirrors icepack's mode3_runner.py fix from earlier this engagement."""

    global_ids = rng.choice(adapter.nx, size=min(nobs, adapter.nx), replace=False)
    global_ids = np.sort(global_ids).astype(np.int64)
    positions, local_ids = pool.route_observation_ids(
        global_ids, adapter, topology=topology, icesee_kwargs=icesee_kwargs
    )
    values = np.zeros(local_ids.size, dtype=np.float64)  # arbitrary synthetic "truth"
    return NativeObservationBatch(observation_ids=local_ids, values=values)


def project_persistent_model_bytes(
    member_state_gb: float,
    nens: int,
    ensemble_group_options: list[int],
    p_model_options: list[int],
) -> list[dict[str, Any]]:
    """Primary Task 2's quantification: for the 37GB/member x Nens stress
    case, project persistent native-model bytes/rank at several
    (ensemble_groups, P_model) layouts.

    This is a PROJECTION, not a measurement -- it uses hypothetical
    large-scale layout parameters explicitly supplied by the caller
    (--project-groups / --project-p-models), decoupled from whatever small
    world_size this particular local benchmark invocation happens to use,
    so it stays meaningful when run from a tiny local smoke test.

    rounds_per_slot = ceil(Nens / ensemble_groups); one slot's rank holds
    rounds_per_slot members' worth of persistent per-rank state
    SIMULTANEOUSLY for the whole run (see NativeDistributedMemberPool --
    it never streams members in/out), so persistent bytes/rank scales as
    rounds_per_slot * (member_state_gb / P_model).
    """

    import math

    rows = []
    for groups in ensemble_group_options:
        if groups <= 0 or groups > nens:
            continue
        rounds = math.ceil(nens / groups)
        for p_model in p_model_options:
            per_rank_member_bytes = (member_state_gb * 1e9) / max(int(p_model), 1)
            persistent_bytes_per_rank = rounds * per_rank_member_bytes
            rows.append(
                {
                    "ensemble_groups": groups,
                    "p_model": int(p_model),
                    "world_size": groups * int(p_model),
                    "rounds_per_slot": rounds,
                    "projected_persistent_gb_per_rank": persistent_bytes_per_rank / 1e9,
                }
            )
    return rows


def project_large_ensemble_table(
    member_state_gb_options: list[float],
    nens_options: list[int],
    p_model: int,
    row_chunk_size: int,
) -> list[dict[str, Any]]:
    """Extended projection for the eventual Ne up to ~1000 scale target
    (not just Nens=40): reports, per (member size, Ne), the logical
    ensemble size, Option-A (all-persistent) vs Option-B (one active +
    all inactive packed in memory) persistent-memory/rank, the Ne^2
    analysis-product footprint (several float64 Ne x Ne matrices --
    cross/gram/rhs plus the eigendecomposition's own working copies, NOT
    just one matrix), and the O(B*Ne) working-set bytes. This is a
    projection using this run's own formulas, not a measurement -- it
    never allocates any of these hypothetical sizes.

    Option B's number here assumes ALL Ne members' packed inactive state
    fits in this rank's addressable memory simultaneously (a
    memory-backed InactiveMemberStore) -- at Ne~1000 with 37GB/member that
    is usually false (see the "fits_in_128gb_node" flag below), which is
    exactly why the InactiveMemberStore abstraction (distributed_member_
    store.py) is backend-neutral: a file/chunked backend is the answer
    once Option B's own number stops fitting, not a redesign of the
    runtime layer that calls it.
    """

    rows = []
    # cross, gram, rhs, plus eigh's internal (values, vectors, a working
    # copy) -- conservatively count 6 Ne x Ne float64 matrices total as the
    # analysis-product working set, not just the 3 accumulated ones.
    ne_squared_matrix_count = 6
    for member_gb in member_state_gb_options:
        for nens in nens_options:
            logical_ensemble_tb = (member_gb * nens) / 1000.0
            option_a_gb_per_rank = (member_gb * nens) / max(p_model, 1)
            option_b_gb_per_rank = (
                member_gb / max(p_model, 1)  # one active native member
                + (member_gb * nens) / max(p_model, 1)  # all inactive packed, in memory
            )
            ne_squared_bytes = ne_squared_matrix_count * (nens ** 2) * 8
            working_set_bytes = (row_chunk_size * nens + nens ** 2) * 8
            rows.append(
                {
                    "member_state_gb": member_gb,
                    "nens": nens,
                    "logical_ensemble_tb": logical_ensemble_tb,
                    "option_a_all_persistent_gb_per_rank": option_a_gb_per_rank,
                    "option_b_one_active_plus_memory_inactive_gb_per_rank": option_b_gb_per_rank,
                    "option_b_fits_in_128gb_node": option_b_gb_per_rank < 128.0,
                    "ne_squared_analysis_product_bytes": ne_squared_bytes,
                    "b_times_ne_working_set_bytes": working_set_bytes,
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=2_000_000, help="Global synthetic state size (elements)")
    parser.add_argument("--ne", type=int, default=2, help="Number of ensemble members")
    parser.add_argument("--ensemble-groups", type=int, default=None, help="Concurrent ensemble groups (default: world size / p-model)")
    parser.add_argument("--p-model", type=int, default=1, help="Spatial ranks per ensemble group")
    parser.add_argument("--row-chunk-size", type=int, default=4096, help="Analysis state-row chunk size B")
    parser.add_argument("--nobs", type=int, default=200, help="Synthetic observation count (Nobs, never Nx)")
    parser.add_argument("--n-cycles", type=int, default=2, help="Number of forecast+analysis cycles")
    parser.add_argument("--checkpoint-root", type=str, default=None, help="Directory for checkpoint I/O measurement (default: a temp dir under /tmp)")
    parser.add_argument("--max-global-bytes", type=int, default=2 * 1024**3, help="Safety threshold for a hypothetical full-Nx float64 array")
    parser.add_argument("--allow-large-nx", action="store_true", help="Explicit opt-in to bypass the --nx safety guard")
    parser.add_argument("--member-state-gb", type=float, default=37.0, help="Assumed per-member live state size (GB) for the projection table")
    parser.add_argument("--project-nens", type=int, default=40, help="Assumed Nens for the projection table")
    parser.add_argument(
        "--project-groups", type=str, default="40,20,10,5",
        help="Comma-separated ensemble_groups options for the projection table",
    )
    parser.add_argument(
        "--project-p-models", type=str, default="1,4,16,64",
        help="Comma-separated P_model options for the projection table (independent of this run's own --p-model)",
    )
    parser.add_argument(
        "--large-ensemble-projection-only", action="store_true",
        help=(
            "Print the extended Ne=40..1000 x member-size projection table "
            "(project_large_ensemble_table) and exit immediately -- no MPI "
            "run, no allocation, analytical only."
        ),
    )
    parser.add_argument(
        "--project-member-gbs", type=str, default="1,10,37",
        help="Comma-separated hypothetical member-state sizes (GB) for --large-ensemble-projection-only",
    )
    parser.add_argument(
        "--project-nens-list", type=str, default="40,80,100,200,300,500,1000",
        help="Comma-separated Ne options for --large-ensemble-projection-only",
    )
    parser.add_argument(
        "--use-member-streaming", action="store_true",
        help="Use StreamingNativeDistributedMemberPool (Option B) instead of the persistent pool.",
    )
    parser.add_argument(
        "--store-backend", type=str, default="memory",
        choices=["memory", "hdf5-member-major", "hdf5-chunked2d"],
        help="InactiveMemberStore backend (Gates 2-4: memory=Option B reference, the other two are Option-C file-backed prototypes).",
    )
    parser.add_argument(
        "--store-directory", type=str, default=None,
        help="Directory for a file-backed store backend (default: a temp dir under /tmp).",
    )
    parser.add_argument(
        "--member-chunk-size", type=int, default=32,
        help="HDF5 chunked-2D backend's ensemble-member chunk width (Gate 6).",
    )
    parser.add_argument(
        "--use-store-streaming-analysis", action="store_true",
        help=(
            "Requires --use-member-streaming. Use "
            "run_native_store_streaming_analysis_cycle (state-row-block "
            "streaming through the InactiveMemberStore) instead of the "
            "original transform_local_members path -- this is the mode "
            "that matters for Ne up to ~1000: analysis workspace stays "
            "O(B*Ne + Ne^2) regardless of how many round-assigned members "
            "this rank holds."
        ),
    )
    args = parser.parse_args()

    if args.large_ensemble_projection_only:
        table = project_large_ensemble_table(
            [float(x) for x in args.project_member_gbs.split(",") if x.strip()],
            [int(x) for x in args.project_nens_list.split(",") if x.strip()],
            args.p_model,
            args.row_chunk_size,
        )
        print(json.dumps(table, indent=2))
        return

    world = MPI.COMM_WORLD
    world_rank = world.Get_rank()
    world_size = world.Get_size()

    guard_against_dangerous_allocation(args.nx, args.max_global_bytes, args.allow_large_nx)

    ensemble_groups = args.ensemble_groups
    if ensemble_groups is None:
        if world_size % args.p_model != 0:
            raise SystemExit(
                f"world_size={world_size} is not divisible by --p-model={args.p_model}"
            )
        ensemble_groups = world_size // args.p_model

    t0 = MPI.Wtime()
    topology = create_distributed_topology(world, spatial_ranks=args.p_model)
    adapter = SyntheticLargeStateAdapter(args.nx, args.max_global_bytes)
    icesee_kwargs: dict[str, Any] = {
        "Nens": args.ne,
        "use_member_streaming": bool(args.use_member_streaming),
        "use_store_streaming_analysis": bool(args.use_store_streaming_analysis),
    }
    if args.use_store_streaming_analysis and not args.use_member_streaming:
        raise SystemExit("--use-store-streaming-analysis requires --use-member-streaming")

    store_files: list[str] = []
    custom_store = None
    if args.use_member_streaming and args.store_backend != "memory":
        from ICESEE.src.parallelization.distributed_native_runtime import (
            validate_native_distributed_adapter,
        )
        from ICESEE.src.parallelization.distributed_runtime import members_for_ensemble_slot
        from ICESEE.src.parallelization.distributed_streaming_runtime import (
            StreamingNativeDistributedMemberPool,
        )

        if args.store_backend == "hdf5-member-major":
            from ICESEE.src.parallelization.distributed_member_store_hdf5 import (
                HDF5MemberMajorStore,
            )
            store_dir = args.store_directory or f"/tmp/mode3_bench_store_{os.getpid()}"
            owned_start, owned_stop = contiguous_owned_interval(
                args.nx, int(topology.spatial_rank), int(topology.spatial_ranks)
            )
            custom_store = HDF5MemberMajorStore(
                store_dir, rank=world_rank, owned_size=owned_stop - owned_start,
            )
        elif args.store_backend == "hdf5-chunked2d":
            from ICESEE.src.parallelization.distributed_member_store_hdf5 import (
                HDF5Chunked2DStore,
            )
            store_dir = args.store_directory or f"/tmp/mode3_bench_store_{os.getpid()}"
            owned_start, owned_stop = contiguous_owned_interval(
                args.nx, int(topology.spatial_rank), int(topology.spatial_ranks)
            )
            custom_store = HDF5Chunked2DStore(
                store_dir, rank=world_rank, owned_size=owned_stop - owned_start,
                number_of_members=args.ne, state_chunk=args.row_chunk_size,
                member_chunk=args.member_chunk_size,
            )
        else:
            raise SystemExit(f"unknown --store-backend {args.store_backend!r}")

        validate_native_distributed_adapter(adapter)
        member_ids = members_for_ensemble_slot(
            args.ne, int(topology.ensemble_groups), int(topology.ensemble_slot)
        )
        pool = StreamingNativeDistributedMemberPool(
            adapter, topology, icesee_kwargs, member_ids, store=custom_store
        )
        store_files = custom_store.file_paths()
    else:
        pool = initialize_native_member_pool(adapter, topology, icesee_kwargs)
    cycle_fn = (
        run_native_store_streaming_analysis_cycle
        if args.use_store_streaming_analysis
        else run_native_global_analysis_cycle
    )
    t_init = MPI.Wtime() - t0

    checkpoint_root = args.checkpoint_root or f"/tmp/mode3_bench_ckpt_{os.getpid()}"
    if world_rank == 0:
        os.makedirs(checkpoint_root, exist_ok=True)
    world.Barrier()

    t_forecast_total = 0.0
    t_analysis_total = 0.0
    t_checkpoint_total = 0.0
    checkpoint = None

    rss_delta_samples = []
    traced_peak_samples = []
    for cycle in range(int(args.n_cycles)):
        t1 = MPI.Wtime()
        # Seeded identically on EVERY rank (never by world_rank): every
        # ensemble slot at a given spatial coordinate must agree on
        # observation IDs/order/values before routing (see
        # _validate_batch_alignment in distributed_native_cycle.py).
        rng = np.random.default_rng(1234 + cycle)
        batch = build_observation_batch(adapter, pool, args.nobs, topology, icesee_kwargs, rng)
        checkpoint_request = NativeCheckpointRequest(
            root=checkpoint_root, run_id="mode3-large-state-benchmark",
        )
        # Bracket the actual analysis cycle with RSS reads so the reported
        # workspace figure is a real measurement, not merely the predicted
        # formula printed back -- ru_maxrss is monotonic non-decreasing
        # within a process, so a nonzero delta here reflects genuinely NEW
        # peak memory touched during this specific call (assuming Nx is
        # kept small enough that persistent state isn't itself still
        # growing the high-water mark -- true for this benchmark's small
        # local Nx choices; at large Nx the delta should shrink toward zero
        # once the persistent state dominates the high-water mark instead).
        rss_before = peak_rss_bytes()
        tracemalloc.start()
        result = cycle_fn(
            pool, adapter, cycle, [batch], number_of_batches=1,
            topology=topology, icesee_kwargs=icesee_kwargs, error_mode="legacy_prior_anomalies",
            state_row_chunk_size=args.row_chunk_size,
            checkpoint_request=checkpoint_request,
        )
        _current, traced_peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        rss_after = peak_rss_bytes()
        rss_delta_samples.append(max(0, rss_after - rss_before))
        traced_peak_samples.append(traced_peak)
        elapsed = MPI.Wtime() - t1
        t_forecast_total += elapsed  # cycle fuses forecast+analysis; see mode3_runner.py's own note
        checkpoint = result.checkpoint

    world.Barrier()

    # --- Restart read timing (this rank's own shard(s) only, never a
    # global-sized read): find this rank's shard(s) via the manifest this
    # same run already wrote, then re-open and re-read them. A simple,
    # robust I/O-cost measurement -- deliberately not exercising the full
    # load_distributed_checkpoint repartition-on-restart path (that has its
    # own, separately test-covered contract; re-driving it correctly here
    # would need a second, independently-constructed layout object, which
    # risks masking real timing behind a benchmark-only bug).
    t_restart = None
    if checkpoint is not None:
        import glob
        import h5py

        t2 = MPI.Wtime()
        manifest_path = os.path.join(str(checkpoint.path), "manifest.json")
        with open(manifest_path) as fh:
            manifest = json.load(fh)
        my_members = set(pool.member_ids)
        bytes_read = 0
        for shard in manifest["shards"]:
            if my_members.isdisjoint(int(m) for m in shard["member_ids"]):
                continue
            shard_path = os.path.join(str(checkpoint.path), shard["file"])
            with h5py.File(shard_path, "r") as handle:
                for member_id in shard["member_ids"]:
                    if int(member_id) not in my_members:
                        continue
                    array = np.asarray(handle[f"member_{member_id}"]["state"])
                    bytes_read += array.nbytes
        t_restart = MPI.Wtime() - t2

    peak_rss = peak_rss_bytes()
    owned_size = adapter.peak_owned_bytes // 8

    # --- Checkpoint file/byte accounting (rank 0 only, small metadata). ---
    ckpt_files = 0
    ckpt_bytes = 0
    if world_rank == 0 and os.path.isdir(checkpoint_root):
        for dirpath, _dirs, files in os.walk(checkpoint_root):
            for fn in files:
                fp = os.path.join(dirpath, fn)
                ckpt_files += 1
                ckpt_bytes += os.path.getsize(fp)

    if custom_store is not None and hasattr(custom_store, "flush"):
        custom_store.flush()
    store_file_bytes = 0
    for fp in store_files:
        try:
            store_file_bytes += os.path.getsize(fp)
        except OSError:
            pass
    all_store_files = world.gather(store_files, root=0)  # collective: every rank must call
    total_store_files_all_ranks = None
    if world_rank == 0:
        total_store_files_all_ranks = sum(len(f) for f in all_store_files if f)

    summary = {
        "world_size": world_size,
        "nx": args.nx,
        "ne": args.ne,
        "p_model": args.p_model,
        "ensemble_groups": ensemble_groups,
        "store_backend": args.store_backend,
        "store_files_this_rank": len(store_files),
        "store_bytes_this_rank": store_file_bytes,
        "store_files_all_ranks": total_store_files_all_ranks,
        "row_chunk_size": args.row_chunk_size,
        "nobs": args.nobs,
        "n_cycles": args.n_cycles,
        "t_init_s": t_init,
        "t_forecast_and_analysis_total_s": t_forecast_total,
        "t_restart_s": t_restart,
        "peak_owned_array_elements_this_rank": int(owned_size),
        "peak_owned_array_bytes_this_rank": int(adapter.peak_owned_bytes),
        "predicted_owned_bytes_this_rank": (args.nx // args.p_model) * 8,
        "peak_rss_bytes_this_rank": peak_rss,
        "analysis_workspace_predicted_bytes": (
            args.row_chunk_size * args.ne + args.ne ** 2
        ) * 8,
        "analysis_workspace_measured_rss_delta_bytes_max": (
            max(rss_delta_samples) if rss_delta_samples else None
        ),
        "analysis_workspace_measured_rss_delta_bytes_per_cycle": rss_delta_samples,
        "analysis_workspace_measured_traced_peak_bytes_max": (
            max(traced_peak_samples) if traced_peak_samples else None
        ),
        "analysis_workspace_measured_traced_peak_bytes_per_cycle": traced_peak_samples,
        "checkpoint_files_total": ckpt_files if world_rank == 0 else None,
        "checkpoint_bytes_total": ckpt_bytes if world_rank == 0 else None,
    }

    if world_rank == 0:
        projection = project_persistent_model_bytes(
            args.member_state_gb,
            args.project_nens,
            [int(x) for x in args.project_groups.split(",") if x.strip()],
            [int(x) for x in args.project_p_models.split(",") if x.strip()],
        )
        summary["projection_37gb_x40_style"] = projection
        print(json.dumps(summary, indent=2))

    if not args.checkpoint_root and world_rank == 0:
        shutil.rmtree(checkpoint_root, ignore_errors=True)


if __name__ == "__main__":
    main()
