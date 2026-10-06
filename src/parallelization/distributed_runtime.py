"""Non-selectable state-only runtime primitives for execution mode 3.

These functions exercise distributed model adapters without changing ICESEE's
runtime dispatch.  They intentionally contain no analysis implementation.  A
mode-3 runner may use them only after forecast parity has been established.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from .distributed_adapter import (
    DistributedBlockStateLayout,
    DistributedStateLayout,
    validate_distributed_adapter,
    validate_block_spatial_partition,
    validate_spatial_partition,
)
from .distributed_analysis import iter_local_ensemble_row_blocks


def members_for_ensemble_slot(
    number_of_members: int,
    ensemble_groups: int,
    ensemble_slot: int,
) -> tuple[int, ...]:
    """Return stable round-robin member IDs assigned to an ensemble group."""

    number_of_members = int(number_of_members)
    ensemble_groups = int(ensemble_groups)
    ensemble_slot = int(ensemble_slot)
    if number_of_members < 0:
        raise ValueError("number_of_members must be nonnegative")
    if ensemble_groups <= 0:
        raise ValueError("ensemble_groups must be positive")
    if not 0 <= ensemble_slot < ensemble_groups:
        raise ValueError("ensemble_slot must identify an ensemble group")
    return tuple(range(ensemble_slot, number_of_members, ensemble_groups))


@dataclass(frozen=True)
class LocalMemberEnsemble:
    """Owned local slabs for the members scheduled on one ensemble group."""

    layout: DistributedStateLayout | DistributedBlockStateLayout
    members: Mapping[int, np.ndarray]

    def __post_init__(self) -> None:
        normalized: dict[int, np.ndarray] = {}
        for member_id, local_state in self.members.items():
            member_id = int(member_id)
            if member_id < 0:
                raise ValueError("ensemble member IDs must be nonnegative")
            state = np.asarray(local_state)
            if state.ndim != 1 or state.size != self.layout.owned_size:
                raise ValueError(
                    f"member {member_id} local state must have shape "
                    f"({self.layout.owned_size},)"
                )
            normalized[member_id] = np.ascontiguousarray(state)
        object.__setattr__(self, "members", normalized)

    @property
    def member_ids(self) -> tuple[int, ...]:
        """Scheduled member IDs in deterministic order."""

        return tuple(sorted(self.members))


def initialize_local_members(
    adapter: Any,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
) -> LocalMemberEnsemble:
    """Initialize locally owned slabs for this rank's scheduled members."""

    validate_distributed_adapter(adapter)
    layout = adapter.distributed_state_layout(
        topology=topology,
        icesee_kwargs=icesee_kwargs,
    )
    if isinstance(layout, DistributedStateLayout):
        validate_spatial_partition(layout, topology.spatial_comm)
    elif isinstance(layout, DistributedBlockStateLayout):
        validate_block_spatial_partition(layout, topology.spatial_comm)
    else:
        raise TypeError(
            "distributed_state_layout must return DistributedStateLayout or "
            "DistributedBlockStateLayout"
        )

    number_of_members = int(icesee_kwargs.get("Nens", 1))
    member_ids = members_for_ensemble_slot(
        number_of_members,
        topology.ensemble_groups,
        topology.ensemble_slot,
    )
    local_members = {
        member_id: adapter.initialize_local_member(
            member_id,
            layout=layout,
            topology=topology,
            icesee_kwargs=icesee_kwargs,
        )
        for member_id in member_ids
    }
    return LocalMemberEnsemble(layout=layout, members=local_members)


def forecast_local_members(
    adapter: Any,
    local_ensemble: LocalMemberEnsemble,
    timestep: int,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
) -> LocalMemberEnsemble:
    """Advance every scheduled member while preserving local ownership."""

    forecast = {
        member_id: adapter.forecast_local_member(
            local_state.copy(),
            member_id,
            int(timestep),
            layout=local_ensemble.layout,
            topology=topology,
            icesee_kwargs=icesee_kwargs,
        )
        for member_id, local_state in local_ensemble.members.items()
    }
    return LocalMemberEnsemble(layout=local_ensemble.layout, members=forecast)


def transform_local_members(
    local_forecast: LocalMemberEnsemble,
    transform: np.ndarray,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
    *,
    row_chunk_size: int,
) -> LocalMemberEnsemble:
    """Right-multiply locally owned member slabs by an ensemble transform.

    The only dense state object is one ``row_chunk_size x Nens`` block shared
    by ranks with the same spatial slab.  This preserves the current
    right-multiplication ``X^a = X^f X5`` without assembling a complete member
    or complete distributed ensemble.  Model-consistency finalization is kept
    outside this pure storage/communication primitive so array adapters and
    persistent model-native adapters use the identical numerical kernel.
    """

    number_of_members = int(icesee_kwargs.get("Nens", 1))
    transform = np.asarray(transform, dtype=np.float64)
    expected = (number_of_members, number_of_members)
    if transform.shape != expected:
        raise ValueError(f"analysis transform must have shape {expected}")

    output = {
        member_id: np.empty(local_forecast.layout.owned_size, dtype=np.float64)
        for member_id in local_forecast.member_ids
    }
    for rows, forecast_block in iter_local_ensemble_row_blocks(
        local_forecast.members,
        number_of_members,
        topology.ensemble_comm,
        row_chunk_size=row_chunk_size,
    ):
        analysis_block = forecast_block @ transform
        for member_id in local_forecast.member_ids:
            output[member_id][rows] = analysis_block[:, member_id]

    return LocalMemberEnsemble(layout=local_forecast.layout, members=output)


def apply_ensemble_transform_local(
    adapter: Any,
    local_forecast: LocalMemberEnsemble,
    transform: np.ndarray,
    timestep: int,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
    *,
    row_chunk_size: int,
) -> LocalMemberEnsemble:
    """Transform local slabs and invoke the existing array-adapter finalizer."""

    transformed = transform_local_members(
        local_forecast,
        transform,
        topology,
        icesee_kwargs,
        row_chunk_size=row_chunk_size,
    )

    finalized = {
        member_id: adapter.finalize_local_analysis(
            local_forecast.members[member_id],
            local_analysis,
            member_id,
            int(timestep),
            layout=local_forecast.layout,
            topology=topology,
            icesee_kwargs=icesee_kwargs,
        )
        for member_id, local_analysis in transformed.members.items()
    }
    return LocalMemberEnsemble(layout=local_forecast.layout, members=finalized)


def reconstruct_global_member_for_parity(
    local_state: np.ndarray,
    layout: DistributedStateLayout,
    spatial_comm: Any,
    *,
    root: int = 0,
) -> np.ndarray | None:
    """Reconstruct one member on ``root`` solely for scientific parity tests.

    Production mode-3 forecast, analysis, and checkpoint paths must not call
    this helper.  It exists so early distributed adapters can be compared with
    mode 2 before the global analysis path is implemented.
    """

    if not isinstance(layout, DistributedStateLayout):
        raise TypeError(
            "contiguous parity reconstruction does not accept segmented layouts"
        )

    state = np.asarray(local_state)
    if state.ndim != 1 or state.size != layout.owned_size:
        raise ValueError("local_state shape does not match the owned state slab")
    payload = (layout.owned_start, layout.owned_stop, np.ascontiguousarray(state))
    gathered = spatial_comm.gather(payload, root=int(root))
    if int(spatial_comm.Get_rank()) != int(root):
        return None

    global_state = np.empty(layout.global_size, dtype=state.dtype)
    cursor = 0
    for start, stop, slab in sorted(gathered, key=lambda item: int(item[0])):
        start, stop = int(start), int(stop)
        slab = np.asarray(slab)
        if start != cursor or stop < start or slab.shape != (stop - start,):
            raise ValueError("gathered state slabs do not form an exact partition")
        global_state[start:stop] = slab
        cursor = stop
    if cursor != layout.global_size:
        raise ValueError("gathered state slabs do not cover the global member")
    return global_state
