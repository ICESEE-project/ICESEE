"""Persistent model-native member lifecycle for execution mode 3.

The numerical model owns its distributed vectors for the lifetime of a local
ensemble member.  ICESEE packs only owned degrees of freedom when analysis or
checkpoint code needs a compact array, then writes the owned analysis back to
the same native vectors.  This avoids rebuilding model objects and prevents a
whole global member from entering the Python analysis layer.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping, Protocol, runtime_checkable

import numpy as np

from .distributed_fields import DistributedFieldRegistry
from .distributed_adapter import (
    DistributedBlockStateLayout,
    DistributedStateLayout,
    validate_block_spatial_partition,
    validate_spatial_partition,
)
from .distributed_runtime import LocalMemberEnsemble, members_for_ensemble_slot


def _layout_signature(layout: Any) -> tuple[Any, ...]:
    """Return an array-safe identity for a distributed ownership layout."""

    blocks = getattr(layout, "blocks", None)
    if blocks is not None:
        return (
            "blocks",
            str(layout.layout_id),
            tuple(
                (
                    str(block.name),
                    int(block.global_offset),
                    int(block.global_size),
                    int(block.owned_start),
                    int(block.owned_stop),
                    tuple(np.asarray(block.ghost_indices, dtype=np.int64).tolist()),
                )
                for block in blocks
            ),
        )
    return (
        "contiguous",
        str(layout.layout_id),
        int(layout.global_size),
        int(layout.owned_start),
        int(layout.owned_stop),
        tuple(np.asarray(layout.ghost_indices, dtype=np.int64).tolist()),
    )


@dataclass
class NativeDistributedMember:
    """Persistent native model state and its ICESEE field registry."""

    member_id: int
    fields: DistributedFieldRegistry
    model_context: Any = None

    def __post_init__(self) -> None:
        self.member_id = int(self.member_id)
        if self.member_id < 0:
            raise ValueError("ensemble member IDs must be nonnegative")
        if not isinstance(self.fields, DistributedFieldRegistry):
            raise TypeError("fields must be a DistributedFieldRegistry")

    def pack_owned(self, *, dtype=None) -> np.ndarray:
        return self.fields.pack_owned(dtype=dtype)

    def unpack_owned(self, local_state: np.ndarray) -> None:
        self.fields.unpack_owned(local_state, synchronize=True)


@dataclass(frozen=True)
class NativeObservationShard:
    """Observation rows and model equivalents owned by one spatial rank.

    ``canonical_positions`` maps the local rows back to the original active
    observation ordering.  This permits arbitrary mesh ownership and sparse
    observation patterns without gathering a global model state.
    """

    canonical_positions: np.ndarray
    observation_ids: np.ndarray
    member_values: Mapping[int, np.ndarray]

    def __post_init__(self) -> None:
        positions = np.asarray(self.canonical_positions)
        rows = np.asarray(self.observation_ids)
        if positions.ndim != 1 or rows.ndim != 1:
            raise ValueError("native observation shard rows must be one-dimensional")
        if not np.issubdtype(positions.dtype, np.integer):
            raise TypeError("canonical observation positions must be integers")
        if not np.issubdtype(rows.dtype, np.integer):
            raise TypeError("native observation IDs must be integers")
        if positions.size != rows.size:
            raise ValueError("native observation positions and IDs must align")
        normalized: dict[int, np.ndarray] = {}
        for member_id, values in self.member_values.items():
            array = np.asarray(values)
            if array.ndim != 1 or array.size != rows.size:
                raise ValueError(
                    "native member observations must match the local row count"
                )
            normalized[int(member_id)] = np.ascontiguousarray(array)
        object.__setattr__(
            self,
            "canonical_positions",
            np.ascontiguousarray(positions, dtype=np.int64),
        )
        object.__setattr__(
            self, "observation_ids", np.ascontiguousarray(rows, dtype=np.int64)
        )
        object.__setattr__(self, "member_values", normalized)


@runtime_checkable
class NativeDistributedModelAdapter(Protocol):
    """Model-owned lifecycle required by the scalable mode-3 path.

    The adapter creates persistent native model objects and advances them in
    place.  Core ICESEE code receives only compact owned snapshots and local
    observation values; a complete member is never part of this interface.
    """

    def initialize_native_member(
        self,
        member_id: int,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> NativeDistributedMember:
        """Create one persistent member on its model-native communicator."""

    def forecast_native_member(
        self,
        member: NativeDistributedMember,
        timestep: int,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> None:
        """Advance one persistent member in place, including halo exchange."""

    def observe_native_member(
        self,
        member: NativeDistributedMember,
        observation_rows: np.ndarray,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        """Evaluate only the observation rows owned by this spatial rank."""

    def finalize_native_analysis(
        self,
        member: NativeDistributedMember,
        forecast_owned: np.ndarray,
        timestep: int,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> None:
        """Apply constraints after analyzed owned values have been restored."""


@runtime_checkable
class NativeDistributedInversionAdapter(Protocol):
    """Optional extension for a distributed hybrid inversion workflow.

    The callback operates collectively on the member's spatial communicator,
    mutates model-native distributed fields in place, and must not return a
    gathered coefficient or state vector.
    """

    def inverse_native_member(
        self,
        member: NativeDistributedMember,
        timestep: int,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> None:
        """Apply one member-wise inversion to native distributed fields."""


_REQUIRED_NATIVE_CALLBACKS = (
    "initialize_native_member",
    "forecast_native_member",
    "observe_native_member",
    "finalize_native_analysis",
)


def validate_native_distributed_adapter(
    adapter: Any,
    *,
    require_inversion: bool = False,
) -> None:
    missing = [
        name for name in _REQUIRED_NATIVE_CALLBACKS
        if not callable(getattr(adapter, name, None))
    ]
    if require_inversion and not callable(
        getattr(adapter, "inverse_native_member", None)
    ):
        missing.append("inverse_native_member")
    if missing:
        raise TypeError(
            "native distributed adapter is missing required callbacks: "
            + ", ".join(missing)
        )


class NativeDistributedMemberPool:
    """Members scheduled on one ensemble group, retained in native storage."""

    def __init__(self, members: Mapping[int, NativeDistributedMember]) -> None:
        normalized = {int(member_id): member for member_id, member in members.items()}
        if not normalized:
            raise ValueError("native member pool requires at least one member")
        for member_id, member in normalized.items():
            if member_id != member.member_id:
                raise ValueError("native member mapping key disagrees with member_id")
        first_member = next(iter(normalized.values()))
        layout = first_member.fields.layout
        signature = _layout_signature(layout)
        if any(
            _layout_signature(member.fields.layout) != signature
            for member in normalized.values()
        ):
            raise ValueError("native members must share one distributed layout")
        self._members = normalized
        self._layout = layout
        self._layout_signature = signature

    @property
    def layout(self):
        return self._layout

    @property
    def member_ids(self) -> tuple[int, ...]:
        return tuple(sorted(self._members))

    def member(self, member_id: int) -> NativeDistributedMember:
        try:
            return self._members[int(member_id)]
        except KeyError as error:
            raise KeyError(f"native member {member_id} is not scheduled here") from error

    def storage_plan(self, *, dtype=np.float64) -> dict[str, int | float]:
        """Describe core owned-state storage without touching global vectors."""

        itemsize = int(np.dtype(dtype).itemsize)
        local_entries = int(self.layout.owned_size) * len(self._members)
        global_entries = int(self.layout.global_size) * len(self._members)
        return {
            "local_members": len(self._members),
            "global_entries_per_member": int(self.layout.global_size),
            "owned_entries_per_member": int(self.layout.owned_size),
            "owned_snapshot_bytes": local_entries * itemsize,
            "global_snapshot_bytes_avoided": global_entries * itemsize,
            "owned_fraction": (
                float(self.layout.owned_size) / float(self.layout.global_size)
                if self.layout.global_size else 0.0
            ),
        }

    def snapshot_owned(self, *, dtype=None) -> LocalMemberEnsemble:
        """Pack a bounded local snapshot for analysis or checkpointing."""

        return LocalMemberEnsemble(
            self.layout,
            {
                member_id: self._members[member_id].pack_owned(dtype=dtype)
                for member_id in self.member_ids
            },
        )

    def restore_owned(self, local_ensemble: LocalMemberEnsemble) -> None:
        """Write local analysis/restart slabs back into persistent vectors."""

        if _layout_signature(local_ensemble.layout) != self._layout_signature:
            raise ValueError("local ensemble layout is incompatible with native members")
        if set(local_ensemble.member_ids) != set(self.member_ids):
            raise ValueError("local ensemble member IDs are incompatible")
        for member_id in self.member_ids:
            self._members[member_id].unpack_owned(local_ensemble.members[member_id])

    def restore_and_finalize(
        self,
        local_analysis: LocalMemberEnsemble,
        local_forecast: LocalMemberEnsemble,
        timestep: int,
        adapter: NativeDistributedModelAdapter,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> None:
        """Restore analyzed slabs, then enforce model-native consistency."""

        if set(local_forecast.member_ids) != set(self.member_ids):
            raise ValueError("forecast member IDs are incompatible")
        self.restore_owned(local_analysis)
        for member_id in self.member_ids:
            result = adapter.finalize_native_analysis(
                self._members[member_id],
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

    def forecast(
        self,
        timestep: int,
        forecast_member: Callable[..., None],
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> LocalMemberEnsemble:
        """Advance native members in place and return only owned snapshots.

        The callback receives the persistent ``NativeDistributedMember`` and
        must operate on its model-native distributed fields.  It must not
        return a global array.
        """

        if not callable(forecast_member):
            raise TypeError("forecast_member must be callable")
        for member_id in self.member_ids:
            result = forecast_member(
                self._members[member_id],
                int(timestep),
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            if result is not None:
                raise TypeError(
                    "native forecast callbacks mutate distributed fields and "
                    "must return None"
                )
        return self.snapshot_owned()

    def forecast_with_adapter(
        self,
        timestep: int,
        adapter: NativeDistributedModelAdapter,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> LocalMemberEnsemble:
        return self.forecast(
            timestep,
            adapter.forecast_native_member,
            topology=topology,
            icesee_kwargs=icesee_kwargs,
        )

    def observe_with_adapter(
        self,
        observation_rows: np.ndarray,
        adapter: NativeDistributedModelAdapter,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> dict[int, np.ndarray]:
        """Evaluate one bounded local observation slab for every member."""

        rows = np.asarray(observation_rows, dtype=np.int64).ravel()
        observed: dict[int, np.ndarray] = {}
        for member_id in self.member_ids:
            values = np.asarray(
                adapter.observe_native_member(
                    self._members[member_id],
                    rows,
                    topology=topology,
                    icesee_kwargs=icesee_kwargs,
                )
            )
            if values.ndim != 1 or values.size != rows.size:
                raise ValueError(
                    "native observation callback must return one value per "
                    "requested local row"
                )
            observed[member_id] = np.ascontiguousarray(values)
        return observed

    def route_observation_ids(
        self,
        global_observation_ids: np.ndarray,
        adapter: NativeDistributedModelAdapter,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> tuple[np.ndarray, np.ndarray]:
        """Partition a canonical global observation-id array to this rank.

        Returns ``(canonical_positions, local_observation_ids)`` without
        evaluating the observation operator. By default, observation IDs are
        variable-major global state rows and the field registry performs
        ownership routing. A nonlinear or non-identity model operator may
        instead implement ``partition_native_observations`` and return
        ``(positions, local_ids)`` itself. In either case storage is
        proportional to the local observation stencil, and the canonical
        positions allow downstream code to select matching data/error
        metadata without a complete-member gather.

        Callers that build a :class:`NativeObservationBatch` for
        ``run_native_global_analysis_cycle`` (whose ``observation_ids`` must
        already be spatially local -- see that function's own docstring) call
        this once per canonical observation set, then subset their own
        values/error arrays by the returned ``canonical_positions`` before
        constructing the batch. ``observe_partitioned_with_adapter`` uses this
        same routing internally when a caller wants routing and evaluation in
        one step instead.
        """

        raw_ids = np.asarray(global_observation_ids)
        if raw_ids.ndim != 1:
            raise ValueError("global observation IDs must be one-dimensional")
        if not np.issubdtype(raw_ids.dtype, np.integer):
            raise TypeError("global observation IDs must be integers")
        observation_ids = raw_ids.astype(np.int64, copy=False)

        representative = self._members[self.member_ids[0]]
        partition = getattr(adapter, "partition_native_observations", None)
        if callable(partition):
            routed = partition(
                representative,
                observation_ids,
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            if not isinstance(routed, tuple) or len(routed) != 2:
                raise TypeError(
                    "partition_native_observations must return "
                    "(canonical_positions, local_observation_ids)"
                )
            positions, local_ids = routed
            positions = np.asarray(positions)
            local_ids = np.asarray(local_ids)
            if positions.ndim != 1 or local_ids.ndim != 1:
                raise ValueError("routed observation arrays must be one-dimensional")
            if not np.issubdtype(positions.dtype, np.integer):
                raise TypeError("canonical observation positions must be integers")
            if not np.issubdtype(local_ids.dtype, np.integer):
                raise TypeError("local observation IDs must be integers")
            positions = positions.astype(np.int64, copy=False)
            local_ids = local_ids.astype(np.int64, copy=False)
            if positions.size != local_ids.size:
                raise ValueError("routed observation positions and IDs must align")
            if positions.size and (
                positions.min() < 0 or positions.max() >= observation_ids.size
            ):
                raise ValueError("routed positions leave the observation array")
            return positions, local_ids
        return representative.fields.partition_owned_rows(observation_ids)

    def observe_partitioned_with_adapter(
        self,
        global_observation_ids: np.ndarray,
        adapter: NativeDistributedModelAdapter,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> NativeObservationShard:
        """Route and evaluate only observations owned by this spatial rank.

        See :meth:`route_observation_ids` for the routing contract. This
        method both routes and evaluates in one step; callers that need to
        route once per canonical observation set (e.g. because ownership is
        static across an entire run) but evaluate every analysis step should
        call :meth:`route_observation_ids` directly instead.
        """

        positions, local_ids = self.route_observation_ids(
            global_observation_ids,
            adapter,
            topology=topology,
            icesee_kwargs=icesee_kwargs,
        )
        values = self.observe_with_adapter(
            local_ids,
            adapter,
            topology=topology,
            icesee_kwargs=icesee_kwargs,
        )
        return NativeObservationShard(positions, local_ids, values)


def native_layout_fingerprint(
    adapter: Any,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
) -> str | None:
    """Optional rank-local identity of the adapter's DOF numbering.

    Adapters whose global state numbering can change between process
    invocations (e.g. a mesh partitioned at runtime) may implement
    ``native_layout_fingerprint(topology=..., icesee_kwargs=...)`` so that
    checkpoints record it and a restart can refuse a mismatched numbering.
    Returns ``None`` when the adapter does not provide one.
    """

    fingerprint = getattr(adapter, "native_layout_fingerprint", None)
    if not callable(fingerprint):
        return None
    value = fingerprint(topology=topology, icesee_kwargs=icesee_kwargs)
    return None if value is None else str(value)


def _restart_member_factory(adapter: Any):
    """Member constructor for a restart: allocate, never initialize.

    ``allocate_native_member`` (optional adapter capability) creates a
    member's native storage without perturbing or forecasting it.  Without
    it the adapter's initializer is the only constructor available; its
    state is then completely overwritten by the checkpoint, so the restart
    is still exact, only slower.
    """

    allocate = getattr(adapter, "allocate_native_member", None)
    if callable(allocate):
        return allocate, True
    return adapter.initialize_native_member, False


def restore_native_member_pool(
    adapter: NativeDistributedModelAdapter,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
    checkpoint_path: Any,
    *,
    expected_run_id: str | None = None,
):
    """Rebuild this rank's member pool from a committed checkpoint.

    Unlike :func:`initialize_native_member_pool` this never generates an
    initial ensemble: members are allocated (see
    ``_restart_member_factory``), then every owned state row is overwritten
    with the checkpoint's values for the members scheduled on this ensemble
    slot under the *current* topology.  Returns ``(pool, checkpoint)``.
    """

    from .distributed_checkpoint import load_distributed_checkpoint

    validate_native_distributed_adapter(adapter)
    number_of_members = int(icesee_kwargs.get("Nens", 1))
    member_ids = members_for_ensemble_slot(
        number_of_members,
        int(topology.ensemble_groups),
        int(topology.ensemble_slot),
    )
    if not member_ids:
        raise ValueError(
            "mode-3 process grid has an ensemble slot with no scheduled member"
        )
    factory, allocates = _restart_member_factory(adapter)
    if not allocates and int(topology.world_rank) == 0:
        print(
            "[ICESEE] restart: adapter has no allocate_native_member; members are "
            "constructed with the initializer and then overwritten from the "
            "checkpoint",
            flush=True,
        )

    if icesee_kwargs.get("use_member_streaming") and callable(
        getattr(adapter, "reactivate_native_member", None)
    ):
        from .distributed_member_store import build_inactive_member_store
        from .distributed_streaming_runtime import StreamingNativeDistributedMemberPool

        store = build_inactive_member_store(
            icesee_kwargs, world_rank=int(topology.world_rank)
        )
        pool = StreamingNativeDistributedMemberPool(
            adapter,
            topology,
            icesee_kwargs,
            member_ids,
            store=store,
            member_factory=factory,
        )
    else:
        pool = NativeDistributedMemberPool(
            {
                member_id: factory(
                    member_id, topology=topology, icesee_kwargs=icesee_kwargs
                )
                for member_id in member_ids
            }
        )
    if isinstance(pool.layout, DistributedStateLayout):
        validate_spatial_partition(pool.layout, topology.spatial_comm)
    elif isinstance(pool.layout, DistributedBlockStateLayout):
        validate_block_spatial_partition(pool.layout, topology.spatial_comm)
    else:  # pragma: no cover - registry construction already guarantees this
        raise TypeError("native member returned an unsupported distributed layout")

    local_ensemble, checkpoint = load_distributed_checkpoint(
        checkpoint_path,
        pool.layout,
        topology,
        expected_run_id=expected_run_id,
        number_of_members=number_of_members,
        expected_layout_fingerprint=native_layout_fingerprint(
            adapter, topology, icesee_kwargs
        ),
    )
    pool.restore_owned(local_ensemble)
    restore = getattr(adapter, "restore_native_checkpoint", None)
    if callable(restore) and isinstance(pool, NativeDistributedMemberPool):
        for member_id in pool.member_ids:
            result = restore(
                pool.member(member_id),
                checkpoint,
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            if result is not None:
                raise TypeError(
                    "native checkpoint restore callbacks mutate model context "
                    "and must return None"
                )
    return pool, checkpoint


def initialize_native_member_pool(
    adapter: NativeDistributedModelAdapter,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
):
    """Create scheduled members and validate spatial ownership.

    Returns a bounded-memory ``StreamingNativeDistributedMemberPool``
    (src/parallelization/distributed_streaming_runtime.py) instead of the
    default persistent ``NativeDistributedMemberPool`` when BOTH:
    (a) ``icesee_kwargs.get("use_member_streaming")`` is truthy (explicit
    opt-in -- never automatic, so every existing caller is unaffected by
    default), and (b) the adapter implements ``reactivate_native_member``
    (duck-typed, exactly like the existing ``partition_native_observations``
    optional-method convention). Any adapter that does not implement
    streaming continues to use the persistent pool unconditionally, so this
    is fully backward compatible.
    """

    validate_native_distributed_adapter(adapter)
    member_ids = members_for_ensemble_slot(
        int(icesee_kwargs.get("Nens", 1)),
        int(topology.ensemble_groups),
        int(topology.ensemble_slot),
    )
    if not member_ids:
        raise ValueError(
            "mode-3 process grid has an ensemble slot with no scheduled member"
        )

    if icesee_kwargs.get("use_member_streaming") and callable(
        getattr(adapter, "reactivate_native_member", None)
    ):
        # Local import: distributed_streaming_runtime.py imports several
        # names FROM this module, so importing it at module load time here
        # would be circular. Deferred to call time instead.
        from .distributed_member_store import build_inactive_member_store
        from .distributed_streaming_runtime import StreamingNativeDistributedMemberPool

        store = build_inactive_member_store(
            icesee_kwargs, world_rank=int(topology.world_rank)
        )
        pool = StreamingNativeDistributedMemberPool(
            adapter, topology, icesee_kwargs, member_ids, store=store
        )
        if isinstance(pool.layout, DistributedStateLayout):
            validate_spatial_partition(pool.layout, topology.spatial_comm)
        elif isinstance(pool.layout, DistributedBlockStateLayout):
            validate_block_spatial_partition(pool.layout, topology.spatial_comm)
        else:  # pragma: no cover - registry construction already guarantees this
            raise TypeError("native member returned an unsupported distributed layout")
        return pool

    pool = NativeDistributedMemberPool(
        {
            member_id: adapter.initialize_native_member(
                member_id,
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            for member_id in member_ids
        }
    )
    if isinstance(pool.layout, DistributedStateLayout):
        validate_spatial_partition(pool.layout, topology.spatial_comm)
    elif isinstance(pool.layout, DistributedBlockStateLayout):
        validate_block_spatial_partition(pool.layout, topology.spatial_comm)
    else:  # pragma: no cover - registry construction already guarantees this
        raise TypeError("native member returned an unsupported distributed layout")
    return pool


def apply_native_memberwise_inversion(
    pool: NativeDistributedMemberPool,
    adapter: NativeDistributedInversionAdapter,
    timestep: int,
    *,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
) -> LocalMemberEnsemble:
    """Apply hybrid inversion transactionally without gathering a member.

    Every spatial rank invokes the adapter for the same locally scheduled
    member IDs.  If any rank reports an error, all owned model fields are
    restored to their pre-inversion values before a collective exception is
    raised.  Model adapters remain responsible for making their native solver
    call collective on ``topology.spatial_comm``.
    """

    inverse = getattr(adapter, "inverse_native_member", None)
    if not callable(inverse):
        raise TypeError("native distributed adapter does not provide inversion")
    before = pool.snapshot_owned()
    local_error = None
    try:
        for member_id in pool.member_ids:
            result = inverse(
                pool.member(member_id),
                int(timestep),
                topology=topology,
                icesee_kwargs=icesee_kwargs,
            )
            if result is not None:
                raise TypeError(
                    "native inversion callbacks mutate distributed fields and "
                    "must return None"
                )
    except Exception as error:  # synchronize failures before publication
        local_error = f"{type(error).__name__}: {error}"

    failures = [value for value in topology.world.allgather(local_error) if value]
    if failures:
        pool.restore_owned(before)
        raise RuntimeError(
            "native member-wise inversion failed; pre-inversion state restored: "
            + "; ".join(failures)
        )
    return pool.snapshot_owned()
