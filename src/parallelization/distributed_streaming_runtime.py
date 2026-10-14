# ==============================================================================
# @des: Bounded-memory round scheduling for execution-mode-3 native members
# (Option B: one active native member per ensemble slot; inactive
# round-assigned members retained as packed distributed-memory state via
# InactiveMemberStore, never as live native/Firedrake/PETSc objects).
# @date: 2026-09-26
# ==============================================================================
"""Drop-in, opt-in alternative to ``NativeDistributedMemberPool``.

Motivation (see docs/execution-mode-3-design.md's fix history and this
engagement's large-state audit): ``NativeDistributedMemberPool`` builds and
holds every member assigned to one ensemble slot for the whole run, so
persistent native-model memory/rank scales as
``rounds_per_slot * native_member_bytes`` whenever ``Nens > ensemble_groups``.
For a hypothetical 37 GB/member x 40-member stress case, that multiplier is
the dominant remaining large-state risk (the observation/analysis side is
already bounded -- see ``distributed_analysis.py``/``distributed_runtime.py``
and the large-state benchmark).

``StreamingNativeDistributedMemberPool`` keeps AT MOST ONE member's native
model objects resident at a time. Every other round-assigned member's owned
state lives only as a plain packed array in an ``InactiveMemberStore``
(``distributed_member_store.py`` -- memory-backed by default; a
checkpoint-backed store could implement the same protocol without any
change here, by design). This targets
persistent native-model memory/rank == O(Nx/P_model) for one active member,
instead of O(rounds_per_slot * Nx/P_model) -- but does NOT eliminate the
underlying O(rounds_per_slot * Nx/P_model) STORAGE requirement for the
packed arrays themselves; it only avoids paying native/PETSc overhead on
top of that for every round simultaneously. Quantify the actual difference
with the large-state benchmark's Option-B mode before assuming it is large.

Strictly backward compatible: nothing here is used unless
``initialize_native_member_pool`` (distributed_native_runtime.py) is asked
to select it, which itself only happens when the adapter implements
``reactivate_native_member`` (duck-typed, exactly like the existing
``partition_native_observations`` optional-method convention) AND the
caller explicitly opts in. Models that implement no such method are
completely unaffected and continue using ``NativeDistributedMemberPool``.

Public method surface intentionally mirrors ``NativeDistributedMemberPool``
exactly (``layout``, ``member_ids``, ``snapshot_owned``, ``restore_owned``,
``restore_and_finalize``, ``forecast_with_adapter``, ``observe_with_adapter``,
``route_observation_ids``, ``observe_partitioned_with_adapter``) so
``distributed_native_cycle.py``'s orchestration functions and every
application's ``mode3_runner.py`` work UNCHANGED against either pool type --
this module makes no changes to that shared cycle-orchestration code at all.

Known, measured tradeoff (see Phase 14 in the overnight report this module
was built for): ``run_native_global_analysis_cycle`` calls
``forecast_with_adapter`` (all members) and, later, ``observe_with_adapter``
(all members) as two SEPARATE top-level calls. Since forecasting a member
here deactivates it immediately afterward, and observation of an inactive
member is answered directly from its packed array (see
``extract_rows_from_packed`` below -- no reactivation needed for
observation), this pool does NOT pay a double activation cost for the
common case. It DOES pay one activation+deactivation per member per
timestep for forecasting, which is the fundamental, unavoidable cost of
bounded-memory rounds; the large-state benchmark measures exactly how much.
"""

from __future__ import annotations

from typing import Any, Iterable, Mapping

import numpy as np

from .distributed_adapter import DistributedBlockStateLayout, DistributedStateLayout
from .distributed_analysis import (
    StochasticAnalysisProducts,
    assemble_selected_ensemble_rows,
    ensemble_transform_from_products,
)
from .distributed_checkpoint import save_distributed_checkpoint
from .distributed_fields import DistributedFieldRegistry
from .distributed_member_store import InactiveMemberStore, MemoryInactiveMemberStore
from .distributed_native_cycle import (
    NativeCheckpointRequest,
    NativeCycleResult,
    NativeObservationBatch,
    _assemble_observation_errors,
    _validate_batch_alignment,
)
from .distributed_native_runtime import (
    NativeDistributedMember,
    NativeObservationShard,
    _layout_signature,
)
from .distributed_runtime import LocalMemberEnsemble, members_for_ensemble_slot


def extract_rows_from_packed(
    layout: DistributedBlockStateLayout,
    packed_array: np.ndarray,
    global_rows: np.ndarray,
) -> np.ndarray:
    """Read identity observations directly from a PACKED owned-state array.

    Mirrors ``DistributedFieldRegistry.observe_owned_rows`` exactly, but
    operates on a plain numpy array (an inactive member's stored state)
    instead of live Firedrake-backed fields -- so observing an inactive
    member never requires reactivating it.
    """

    rows = np.asarray(global_rows, dtype=np.int64).ravel()
    packed = np.asarray(packed_array)
    if rows.size and (rows.min() < 0 or rows.max() >= layout.global_size):
        raise ValueError("global observation rows leave the distributed state")
    result = np.empty(rows.size, dtype=packed.dtype)
    matched = np.zeros(rows.size, dtype=bool)
    for _name, local, (global_start, global_stop) in layout.iter_owned_global_intervals():
        selected = (rows >= global_start) & (rows < global_stop)
        if not np.any(selected):
            continue
        result[selected] = packed[local.start + (rows[selected] - global_start)]
        matched[selected] = True
    if not np.all(matched):
        missing = rows[~matched]
        preview = ", ".join(str(int(row)) for row in missing[:5])
        suffix = "..." if missing.size > 5 else ""
        raise ValueError(
            "observation rows are not owned by this spatial rank: "
            f"{preview}{suffix}"
        )
    return np.ascontiguousarray(result)


def partition_owned_rows_from_layout(
    layout: DistributedBlockStateLayout,
    global_rows: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Layout-only equivalent of ``DistributedFieldRegistry.partition_owned_rows``.

    Ownership is a pure function of the shared layout (identical for every
    member of one ensemble group -- they all share one Firedrake mesh
    partition), so this never needs a live member either.
    """

    raw_rows = np.asarray(global_rows)
    if raw_rows.ndim != 1:
        raise ValueError("global observation rows must be one-dimensional")
    if not np.issubdtype(raw_rows.dtype, np.integer):
        raise TypeError("global observation rows must use an integer dtype")
    rows = raw_rows.astype(np.int64, copy=False)
    if rows.size and (rows.min() < 0 or rows.max() >= layout.global_size):
        raise ValueError("global observation rows leave the distributed state")
    owned = np.zeros(rows.size, dtype=bool)
    for _name, _local, (start, stop) in layout.iter_owned_global_intervals():
        owned |= (rows >= start) & (rows < stop)
    positions = np.flatnonzero(owned).astype(np.int64, copy=False)
    return positions, rows[positions]


class StreamingNativeDistributedMemberPool:
    """Bounded-memory (at most one live native member) member pool.

    See module docstring for the full design rationale and the exact
    tradeoff this makes.
    """

    def __init__(
        self,
        adapter: Any,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
        member_ids: tuple[int, ...],
        *,
        store: InactiveMemberStore | None = None,
        member_factory: Any = None,
    ) -> None:
        reactivate = getattr(adapter, "reactivate_native_member", None)
        if not callable(reactivate):
            raise TypeError(
                f"{type(adapter).__name__} does not implement "
                "reactivate_native_member; it cannot be used with "
                "StreamingNativeDistributedMemberPool"
            )
        self._adapter = adapter
        self._topology = topology
        self._icesee_kwargs = icesee_kwargs
        self._member_ids = tuple(sorted(int(m) for m in member_ids))
        if not self._member_ids:
            raise ValueError("streaming member pool requires at least one member")
        self._store = store if store is not None else MemoryInactiveMemberStore()
        self._active_id: int | None = None
        self._active: NativeDistributedMember | None = None
        self._layout: DistributedStateLayout | DistributedBlockStateLayout | None = None
        self._layout_signature: tuple[Any, ...] | None = None

        # Build each member ONCE via the adapter's existing initial-state
        # constructor (unchanged contract), immediately pack + store, then
        # release it -- construction itself never holds more than one live
        # native member at a time either. A restart passes a
        # ``member_factory`` that allocates members without initializing
        # them; restore_owned then overwrites the stored state.
        factory = (
            member_factory
            if member_factory is not None
            else adapter.initialize_native_member
        )
        for member_id in self._member_ids:
            member = factory(
                member_id, topology=topology, icesee_kwargs=icesee_kwargs
            )
            if not isinstance(member.fields, DistributedFieldRegistry):
                raise TypeError("native members must expose a DistributedFieldRegistry")
            signature = _layout_signature(member.fields.layout)
            if self._layout is None:
                self._layout = member.fields.layout
                self._layout_signature = signature
            elif signature != self._layout_signature:
                raise ValueError("native members must share one distributed layout")
            self._store.put(member_id, member.pack_owned())
            del member

    @property
    def layout(self):
        return self._layout

    @property
    def member_ids(self) -> tuple[int, ...]:
        return self._member_ids

    def _deactivate(self) -> None:
        if self._active is not None:
            self._store.put(self._active_id, self._active.pack_owned())
            self._active = None
            self._active_id = None

    def _activate(self, member_id: int, *, seed_state: np.ndarray | None = None) -> NativeDistributedMember:
        member_id = int(member_id)
        if self._active_id == member_id and seed_state is None:
            return self._active
        self._deactivate()
        packed = seed_state if seed_state is not None else self._store.get(member_id)
        member = self._adapter.reactivate_native_member(
            member_id,
            np.asarray(packed, dtype=np.float64),
            topology=self._topology,
            icesee_kwargs=self._icesee_kwargs,
        )
        self._active = member
        self._active_id = member_id
        return member

    def _packed_state_for(self, member_id: int) -> np.ndarray:
        if member_id == self._active_id:
            return self._active.pack_owned()
        return self._store.get(member_id)

    def snapshot_owned(self, *, dtype=None) -> LocalMemberEnsemble:
        members = {}
        for member_id in self._member_ids:
            packed = self._packed_state_for(member_id)
            members[member_id] = np.asarray(packed, dtype=dtype) if dtype else packed
        return LocalMemberEnsemble(self._layout, members)

    def restore_owned(self, local_ensemble: LocalMemberEnsemble) -> None:
        if _layout_signature(local_ensemble.layout) != self._layout_signature:
            raise ValueError("local ensemble layout is incompatible with native members")
        if set(local_ensemble.member_ids) != set(self._member_ids):
            raise ValueError("local ensemble member IDs are incompatible")
        for member_id in self._member_ids:
            self._store.put(member_id, np.asarray(local_ensemble.members[member_id]))
        # Any currently active member's live state is now stale relative to
        # what was just written into the store -- drop it so the next
        # activation reads the freshly restored array.
        self._active = None
        self._active_id = None

    def restore_and_finalize(
        self,
        local_analysis: LocalMemberEnsemble,
        local_forecast: LocalMemberEnsemble,
        timestep: int,
        adapter: Any,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> None:
        if set(local_forecast.member_ids) != set(self._member_ids):
            raise ValueError("forecast member IDs are incompatible")
        for member_id in self._member_ids:
            analyzed = np.asarray(local_analysis.members[member_id])
            member = self._activate(member_id, seed_state=analyzed)
            result = adapter.finalize_native_analysis(
                member,
                np.asarray(local_forecast.members[member_id]),
                int(timestep),
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            if result is not None:
                raise TypeError(
                    "native analysis finalizers mutate distributed fields and "
                    "must return None"
                )
            self._deactivate()

    def forecast_with_adapter(
        self,
        timestep: int,
        adapter: Any,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> LocalMemberEnsemble:
        forecast_members: dict[int, np.ndarray] = {}
        for member_id in self._member_ids:
            member = self._activate(member_id)
            result = adapter.forecast_native_member(
                member, int(timestep), topology=topology, icesee_kwargs=icesee_kwargs
            )
            if result is not None:
                raise TypeError(
                    "native forecast callbacks mutate distributed fields and "
                    "must return None"
                )
            forecast_members[member_id] = member.pack_owned()
            self._deactivate()
        return LocalMemberEnsemble(self._layout, forecast_members)

    def forecast_streaming(
        self,
        timestep: int,
        adapter: Any,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> None:
        """Bounded-memory counterpart of ``forecast_with_adapter``: advances
        every round-assigned member exactly the same way (activate, forecast,
        pack, deactivate -- ``_deactivate`` already persists into the store),
        but never accumulates a ``{member_id: full_owned_array}`` return
        dict. Use this instead of ``forecast_with_adapter`` when the caller
        will read forecast values back via the store (``get``/``get_rows``)
        rather than needing them all resident at once -- see
        ``run_native_store_streaming_analysis_cycle``.
        """

        for member_id in self._member_ids:
            member = self._activate(member_id)
            result = adapter.forecast_native_member(
                member, int(timestep), topology=topology, icesee_kwargs=icesee_kwargs
            )
            if result is not None:
                raise TypeError(
                    "native forecast callbacks mutate distributed fields and "
                    "must return None"
                )
            self._deactivate()

    def finalize_all_from_store(
        self,
        adapter: Any,
        timestep: int,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> None:
        """Bounded-memory counterpart of ``restore_and_finalize``: reactivates
        each round-assigned member directly from its (already analyzed, via
        ``transform_members_via_store``) stored state, calls
        ``finalize_native_analysis``, and deactivates -- no external
        ``LocalMemberEnsemble`` needed for either the analyzed or forecast
        arrays.

        KNOWN LIMITATION (documented, not silently worked around): the
        generic ``finalize_native_analysis`` contract also accepts a
        ``forecast_owned`` array for finalizers that need the PRE-analysis
        forecast value (e.g. a positivity clamp referencing the forecast).
        Since the streaming path never keeps that forecast snapshot
        resident once written to the store, this passes an empty array as a
        placeholder. Icepack's own ``finalize_analysis`` callback is
        ``None`` (unused, confirmed dead code for this argument), so this is
        safe for icepack today; a future model whose finalizer genuinely
        needs the forecast value would need a second, forecast-snapshot
        store namespace -- flagged, not implemented, since no current
        model needs it.
        """

        placeholder_forecast = np.empty(0, dtype=np.float64)
        for member_id in self._member_ids:
            analyzed = self._store.get(member_id)
            member = self._activate(member_id, seed_state=analyzed)
            result = adapter.finalize_native_analysis(
                member,
                placeholder_forecast,
                int(timestep),
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            if result is not None:
                raise TypeError(
                    "native analysis finalizers mutate distributed fields and "
                    "must return None"
                )
            self._deactivate()

    def observe_with_adapter(
        self,
        observation_rows: np.ndarray,
        adapter: Any,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> dict[int, np.ndarray]:
        rows = np.asarray(observation_rows, dtype=np.int64).ravel()
        observed: dict[int, np.ndarray] = {}
        for member_id in self._member_ids:
            packed = self._packed_state_for(member_id)
            values = extract_rows_from_packed(self._layout, packed, rows)
            if values.size != rows.size:
                raise ValueError(
                    "native observation callback must return one value per "
                    "requested local row"
                )
            observed[member_id] = values
        return observed

    def route_observation_ids(
        self,
        global_observation_ids: np.ndarray,
        adapter: Any,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> tuple[np.ndarray, np.ndarray]:
        partition = getattr(adapter, "partition_native_observations", None)
        if callable(partition):
            # Needs a live representative member for a nonlinear/non-identity
            # observation operator; activate one transiently (deactivated
            # again immediately, so this never leaves more than one member
            # resident).
            representative = self._activate(self._member_ids[0])
            routed = partition(
                representative,
                np.asarray(global_observation_ids, dtype=np.int64),
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            self._deactivate()
            if not isinstance(routed, tuple) or len(routed) != 2:
                raise TypeError(
                    "partition_native_observations must return "
                    "(canonical_positions, local_observation_ids)"
                )
            return routed
        return partition_owned_rows_from_layout(self._layout, global_observation_ids)

    def observe_partitioned_with_adapter(
        self,
        global_observation_ids: np.ndarray,
        adapter: Any,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> NativeObservationShard:
        positions, local_ids = self.route_observation_ids(
            global_observation_ids, adapter, topology=topology, icesee_kwargs=icesee_kwargs
        )
        values = self.observe_with_adapter(
            local_ids, adapter, topology=topology, icesee_kwargs=icesee_kwargs
        )
        return NativeObservationShard(positions, local_ids, values)


def transform_members_via_store(
    store: InactiveMemberStore,
    member_ids: tuple[int, ...],
    owned_size: int,
    transform: np.ndarray,
    ensemble_comm: Any,
    number_of_members: int,
    *,
    row_chunk_size: int,
) -> None:
    """Store-streaming counterpart of ``distributed_runtime.transform_local_members``.

    Mathematically identical: for each state-row block ``I``, this computes
    ``X_I^a = X_I @ transform`` (a pure per-block right-multiplication, no
    dependence on rows outside ``I`` -- see this session's derivation of
    ``iter_local_ensemble_row_blocks``, which established the SAME equation
    and confirmed row-block separability holds for the existing ensemble
    transform). The only difference is where each block's data comes from
    (``store.get_rows`` instead of a pre-resident ``{member_id: full_array}``
    dict) and where the analyzed block goes (``store.put_rows`` instead of
    an in-memory output dict) -- so this never requires this rank's
    round-assigned members' full owned arrays to be simultaneously resident.

    The cross-rank ``ensemble_comm.allgather`` below exists for the SAME
    reason as in ``iter_local_ensemble_row_blocks``: different ensemble
    GROUPS (not different spatial ranks -- ``ensemble_comm`` connects ranks
    at the SAME spatial coordinate across ensemble groups) hold different
    members' data, so building one ``block_rows x Ne`` block requires
    combining contributions from every ensemble group at this spatial
    coordinate. This is classification (A) from the derivation: gathering
    different ensemble MEMBERS for the SAME owned state rows -- never a
    spatial gather, never a global-state reconstruction.
    """

    owned_size = int(owned_size)
    row_chunk_size = int(row_chunk_size)
    number_of_members = int(number_of_members)
    bulk_get = getattr(store, "get_ensemble_rows", None)
    bulk_put = getattr(store, "put_ensemble_rows", None)
    use_bulk = callable(bulk_get) and callable(bulk_put)
    for start in range(0, owned_size, row_chunk_size):
        stop = min(owned_size, start + row_chunk_size)
        row_slice = slice(start, stop)
        if use_bulk:
            # ONE physical store operation for all of THIS rank's own
            # round-assigned members, instead of len(member_ids) separate
            # ones -- see Gate 7: a filesystem backend doing one small
            # read per member per row block scales terribly at Ne~1000
            # even though it is free for an in-memory backend.
            local_block = np.asarray(bulk_get(member_ids, row_slice))
            payload = {
                int(member_id): np.ascontiguousarray(local_block[:, i])
                for i, member_id in enumerate(member_ids)
            }
        else:
            payload = {
                int(member_id): np.ascontiguousarray(
                    store.get_rows(int(member_id), row_slice)
                )
                for member_id in member_ids
            }
        gathered = ensemble_comm.allgather(payload)
        columns: dict[int, np.ndarray] = {}
        for rank_payload in gathered:
            for member_id, values in rank_payload.items():
                member_id = int(member_id)
                if member_id in columns:
                    raise ValueError("ensemble member is owned by multiple slots")
                columns[member_id] = np.asarray(values)
        expected = set(range(number_of_members))
        if set(columns) != expected:
            missing = sorted(expected - set(columns))
            extra = sorted(set(columns) - expected)
            raise ValueError(
                f"ensemble member ownership is incomplete; missing={missing}, "
                f"extra={extra}"
            )
        block_dtype = np.result_type(np.float64, *[value.dtype for value in columns.values()])
        block = np.empty((stop - start, number_of_members), dtype=block_dtype)
        for member_id in range(number_of_members):
            values = columns[member_id]
            if values.shape != (stop - start,):
                raise ValueError("gathered member chunk has an invalid shape")
            block[:, member_id] = values
        analysis_block = block @ transform
        if use_bulk:
            local_analyzed = np.stack(
                [analysis_block[:, int(member_id)] for member_id in member_ids], axis=1
            )
            bulk_put(member_ids, row_slice, local_analyzed)
        else:
            for member_id in member_ids:
                store.put_rows(
                    int(member_id), row_slice, analysis_block[:, int(member_id)]
                )


def run_native_store_streaming_analysis_cycle(
    pool: StreamingNativeDistributedMemberPool,
    adapter: Any,
    timestep: int,
    observation_batches: Iterable[NativeObservationBatch],
    *,
    number_of_batches: int,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
    error_mode: str,
    state_row_chunk_size: int = 4096,
    energy_fraction: float = 0.999,
    checkpoint_request: NativeCheckpointRequest | None = None,
) -> NativeCycleResult:
    """Store-streaming sibling of ``distributed_native_cycle.run_native_global_analysis_cycle``.

    Kept as a SEPARATE function, not a modification of the original --
    per this session's directive, the original persistent/Option-B analysis
    path (``transform_local_members``) is preserved unmodified as a
    correctness reference. This function requires a
    ``StreamingNativeDistributedMemberPool`` specifically (it uses
    ``pool.forecast_streaming``/``pool.finalize_all_from_store``, and reaches
    into ``pool._store`` for the streaming transform) -- it is not meant to
    be a generic drop-in for the persistent pool.

    Observation handling is IDENTICAL to the original (byte-for-byte the
    same code): ``pool.observe_with_adapter`` already reads directly from
    the store for a streaming pool (see ``StreamingNativeDistributedMemberPool.
    observe_with_adapter``/``extract_rows_from_packed``), so no change was
    needed there. Only the state-row transform step is replaced with
    ``transform_members_via_store``, and only the forecast/finalize calls
    use the streaming (non-dict-accumulating) pool methods.

    KNOWN REMAINING COST (documented, not hidden): ``pool.snapshot_owned()``
    for the return value and for ``checkpoint_request`` still assembles one
    ``{member_id: full_owned_array}`` dict for this rank's round-assigned
    members, same as the original path -- checkpoint writing itself is not
    yet made to stream row-blocks directly from the store. This is a
    separate, Option-C-adjacent piece of work, not attempted here.
    """

    timestep = int(timestep)
    if timestep < 0:
        raise ValueError("timestep must be nonnegative")
    number_of_batches = int(number_of_batches)
    if number_of_batches < 0:
        raise ValueError("number_of_batches must be nonnegative")
    counts = topology.ensemble_comm.allgather(number_of_batches)
    if any(int(value) != number_of_batches for value in counts):
        raise ValueError("ensemble slots disagree on local observation batch count")

    number_of_members = int(icesee_kwargs.get("Nens", 1))
    pool.forecast_streaming(
        timestep, adapter, topology=topology, icesee_kwargs=icesee_kwargs
    )

    products = StochasticAnalysisProducts.zeros(number_of_members)
    batches = iter(observation_batches)
    for batch_index in range(number_of_batches):
        try:
            batch = next(batches)
        except StopIteration as error:
            raise ValueError(
                "observation source ended before number_of_batches"
            ) from error
        if not isinstance(batch, NativeObservationBatch):
            raise TypeError("observation source must yield NativeObservationBatch")
        _validate_batch_alignment(batch, topology.ensemble_comm)
        forecast_values = pool.observe_with_adapter(
            batch.observation_ids,
            adapter,
            topology=topology,
            icesee_kwargs=icesee_kwargs,
        )
        rows = np.arange(batch.observation_ids.size, dtype=np.int64)
        forecast_matrix = assemble_selected_ensemble_rows(
            forecast_values,
            rows,
            number_of_members,
            topology.ensemble_comm,
        ).astype(np.float64, copy=False)
        error_matrix = _assemble_observation_errors(
            batch,
            number_of_members=number_of_members,
            ensemble_comm=topology.ensemble_comm,
        )
        if topology.is_ensemble_root:
            products.add_chunk(
                forecast_matrix,
                batch.values,
                error_mode=error_mode,
                observation_errors=error_matrix,
            )

    try:
        next(batches)
    except StopIteration:
        pass
    else:
        raise ValueError("observation source exceeded number_of_batches")

    if topology.is_ensemble_root:
        global_products = products.reduced(topology.spatial_comm)
        transform = ensemble_transform_from_products(
            global_products, energy_fraction=energy_fraction
        )
        observation_rows = int(global_products.observation_rows)
    else:
        transform = None
        observation_rows = None
    transform = np.asarray(
        topology.ensemble_comm.bcast(transform, root=0), dtype=np.float64
    )
    observation_rows = int(topology.ensemble_comm.bcast(observation_rows, root=0))

    # Safe, provable optimization (this session's Phase-1 operation-count
    # derivation, validated against a real Icepack run): when
    # number_of_batches == 0, no chunk ever reaches products.add_chunk, so
    # cross/gram/rhs stay exactly zero and ensemble_transform_from_products
    # returns exactly np.eye(Ne) (verified bit-for-bit, not assumed) --
    # X_I @ eye(Ne) == X_I for every row block, so the full store-streaming
    # read/transform/write pass is a mathematically exact no-op on every
    # non-analysis timestep. Skipping it here leaves the store holding
    # exactly the forecast state transform_members_via_store would have
    # written back unchanged, so finalize_all_from_store below still runs
    # unconditionally (preserving every adapter's own finalize/ghost-sync
    # semantics exactly, including for a future model whose finalizer is
    # NOT a no-op) -- only the row-block MATH pass is skipped, never the
    # per-model finalize. `number_of_batches` is already validated
    # identical across every rank in `topology.ensemble_comm` above, and is
    # computed identically on every world rank by every caller (never
    # rank-dependent), so this branch is collective-safe with no new
    # communication. This is exactly the wasted-I/O pattern this session's
    # operation-count derivation found: 25 timesteps x 67 row blocks =
    # 1675 bulk HDF5 operations measured for a run with only 2 real
    # analysis events -- 23 of those 25 passes were provably no-ops.
    if number_of_batches > 0:
        transform_members_via_store(
            pool._store,
            pool.member_ids,
            pool.layout.owned_size,
            transform,
            topology.ensemble_comm,
            number_of_members,
            row_chunk_size=state_row_chunk_size,
        )
    pool.finalize_all_from_store(
        adapter, timestep, topology=topology, icesee_kwargs=icesee_kwargs
    )

    local_analysis = pool.snapshot_owned()

    checkpoint = None
    if checkpoint_request is not None:
        metadata = dict(checkpoint_request.metadata or {})
        metadata.update({
            "cycle_timestep": timestep,
            "observation_rows": observation_rows,
            "error_mode": str(error_mode),
        })
        checkpoint = save_distributed_checkpoint(
            checkpoint_request.root,
            timestep,
            pool.snapshot_owned(),
            topology,
            run_id=checkpoint_request.run_id,
            metadata=metadata,
        )

    return NativeCycleResult(
        local_forecast=local_analysis,
        local_analysis=local_analysis,
        transform=transform,
        observation_rows=observation_rows,
        checkpoint=checkpoint,
    )
