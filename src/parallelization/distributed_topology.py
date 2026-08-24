"""Two-dimensional MPI topology for the proposed execution mode 3.

This module intentionally does not register execution mode 3.  It provides the
first independently testable building block for distributing work across both
ensemble groups and model spatial subdomains.  Runtime selection remains
disabled until the distributed adapter and analysis parity contracts exist.
"""

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class DistributedTopology:
    """Communicators and coordinates in an ensemble-by-space process grid."""

    world: Any
    spatial_comm: Any
    ensemble_comm: Any
    world_rank: int
    world_size: int
    ensemble_slot: int
    spatial_rank: int
    ensemble_groups: int
    spatial_ranks: int

    @property
    def is_spatial_root(self) -> bool:
        """Whether this rank is the root of its member's spatial communicator."""

        return self.spatial_rank == 0

    @property
    def is_ensemble_root(self) -> bool:
        """Whether this rank is first in its same-spatial-slab communicator."""

        return self.ensemble_slot == 0


def create_distributed_topology(
    world: Any,
    *,
    ensemble_groups: int | None = None,
    spatial_ranks: int | None = None,
) -> DistributedTopology:
    """Split ``world`` into orthogonal ensemble and spatial communicators.

    Exactly one dimension may be inferred from the world size.  The initial
    contract requires a rectangular process grid with no idle ranks; later
    schedulers may map multiple ensemble members onto each ensemble group.

    Rank ordering is ensemble-major.  For a process-grid coordinate ``(e, s)``,
    ``spatial_comm`` contains all ``s`` for fixed ``e`` and ``ensemble_comm``
    contains all ``e`` for fixed ``s``.
    """

    world_size = int(world.Get_size())
    world_rank = int(world.Get_rank())
    if world_size <= 0:
        raise ValueError("MPI world communicator must contain at least one rank")

    if ensemble_groups is None and spatial_ranks is None:
        raise ValueError("set ensemble_groups, spatial_ranks, or both")

    if ensemble_groups is not None:
        ensemble_groups = int(ensemble_groups)
        if ensemble_groups <= 0:
            raise ValueError("ensemble_groups must be positive")
    if spatial_ranks is not None:
        spatial_ranks = int(spatial_ranks)
        if spatial_ranks <= 0:
            raise ValueError("spatial_ranks must be positive")

    if ensemble_groups is None:
        if world_size % spatial_ranks:
            raise ValueError("world size must be divisible by spatial_ranks")
        ensemble_groups = world_size // spatial_ranks
    elif spatial_ranks is None:
        if world_size % ensemble_groups:
            raise ValueError("world size must be divisible by ensemble_groups")
        spatial_ranks = world_size // ensemble_groups

    if ensemble_groups * spatial_ranks != world_size:
        raise ValueError(
            "ensemble_groups * spatial_ranks must equal the MPI world size"
        )

    ensemble_slot, spatial_rank = divmod(world_rank, spatial_ranks)
    spatial_comm = world.Split(color=ensemble_slot, key=spatial_rank)
    ensemble_comm = world.Split(color=spatial_rank, key=ensemble_slot)

    return DistributedTopology(
        world=world,
        spatial_comm=spatial_comm,
        ensemble_comm=ensemble_comm,
        world_rank=world_rank,
        world_size=world_size,
        ensemble_slot=ensemble_slot,
        spatial_rank=spatial_rank,
        ensemble_groups=ensemble_groups,
        spatial_ranks=spatial_ranks,
    )
