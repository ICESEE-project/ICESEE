"""Bounded-memory planning utilities for execution mode 3.

The planner is deliberately analytical: it can evaluate a state that is too
large to allocate on the current machine.  It estimates arrays owned by the
ICESEE analysis/runtime layer and keeps model-native solver vectors and third-
party MPI/HDF5 buffers separate because those costs are adapter and platform
dependent.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


def _positive_integer(name: str, value: int) -> int:
    normalized = int(value)
    if normalized <= 0:
        raise ValueError(f"{name} must be positive")
    return normalized


@dataclass(frozen=True)
class Mode3MemoryPlan:
    """Conservative per-rank memory estimate for the native mode-3 cycle."""

    global_rows: int
    owned_rows: int
    total_members: int
    ensemble_groups: int
    spatial_ranks: int
    local_members: int
    state_row_chunk_size: int
    observation_row_chunk_size: int
    itemsize: int
    owned_snapshots_bytes: int
    state_transform_workspace_bytes: int
    observation_workspace_bytes: int
    ensemble_products_bytes: int
    minimum_native_state_bytes: int

    @property
    def estimated_peak_icesee_bytes(self) -> int:
        return (
            self.owned_snapshots_bytes
            + self.state_transform_workspace_bytes
            + self.observation_workspace_bytes
            + self.ensemble_products_bytes
        )

    @property
    def estimated_minimum_peak_bytes(self) -> int:
        """ICESEE work arrays plus one resident native state copy.

        Real model peaks can be higher because nonlinear solvers commonly own
        additional vectors.  Adapters should measure that extra contribution
        during their application-specific scale gate.
        """

        return self.estimated_peak_icesee_bytes + self.minimum_native_state_bytes

    @property
    def complete_local_members_bytes(self) -> int:
        """Memory that whole-member scheduling would require on this rank."""

        return self.global_rows * self.local_members * self.itemsize

    @property
    def complete_ensemble_bytes(self) -> int:
        return self.global_rows * self.total_members * self.itemsize


def estimate_mode3_memory(
    *,
    global_rows: int,
    total_members: int,
    ensemble_groups: int,
    spatial_ranks: int,
    state_row_chunk_size: int = 4096,
    observation_row_chunk_size: int = 4096,
    dtype=np.float64,
) -> Mode3MemoryPlan:
    """Estimate the maximum ICESEE-owned working memory on one rank.

    The estimate assumes balanced spatial ownership and round-robin member
    scheduling.  It is conservative for the current native cycle:

    * two compact owned snapshots coexist (forecast and analysis);
    * a state transform holds input and output chunks across all members;
    * an observation batch holds model equivalents and perturbations; and
    * four ``Nens x Nens`` analysis products/transforms coexist.

    Persistent model-native fields, solver work vectors, Python object
    overhead, filesystem caches, and MPI implementation buffers are excluded
    and must be added by the application-specific deployment plan.
    """

    global_rows = _positive_integer("global_rows", global_rows)
    total_members = _positive_integer("total_members", total_members)
    ensemble_groups = _positive_integer("ensemble_groups", ensemble_groups)
    spatial_ranks = _positive_integer("spatial_ranks", spatial_ranks)
    state_row_chunk_size = _positive_integer(
        "state_row_chunk_size", state_row_chunk_size
    )
    observation_row_chunk_size = _positive_integer(
        "observation_row_chunk_size", observation_row_chunk_size
    )
    if ensemble_groups > total_members:
        raise ValueError("ensemble_groups cannot exceed total_members")

    itemsize = int(np.dtype(dtype).itemsize)
    owned_rows = int(math.ceil(global_rows / spatial_ranks))
    local_members = int(math.ceil(total_members / ensemble_groups))
    state_chunk = min(owned_rows, state_row_chunk_size)

    owned_snapshots = 2 * owned_rows * local_members * itemsize
    state_workspace = 2 * state_chunk * total_members * itemsize
    observation_workspace = (
        2 * observation_row_chunk_size * total_members * itemsize
    )
    ensemble_products = 4 * total_members * total_members * itemsize
    minimum_native_state = owned_rows * local_members * itemsize

    return Mode3MemoryPlan(
        global_rows=global_rows,
        owned_rows=owned_rows,
        total_members=total_members,
        ensemble_groups=ensemble_groups,
        spatial_ranks=spatial_ranks,
        local_members=local_members,
        state_row_chunk_size=state_row_chunk_size,
        observation_row_chunk_size=observation_row_chunk_size,
        itemsize=itemsize,
        owned_snapshots_bytes=owned_snapshots,
        state_transform_workspace_bytes=state_workspace,
        observation_workspace_bytes=observation_workspace,
        ensemble_products_bytes=ensemble_products,
        minimum_native_state_bytes=minimum_native_state,
    )
