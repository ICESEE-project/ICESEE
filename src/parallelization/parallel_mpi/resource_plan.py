# ==============================================================================
# @des: Centralized, MPI-free resource planning for ICESEE's hierarchical
#       ensemble/model MPI topology.
#
#       Separates two independent dimensions of parallelism that prior
#       code conflated into a single `Nens >= size_world` branch condition:
#       ensemble parallelism (how many members run concurrently) and
#       model/spatial parallelism (how many ranks cooperate on one
#       member). This module computes *only* the plan -- which world rank
#       belongs to which model group, which member each group handles in
#       each round -- as a pure function of (world_size, Nens,
#       ranks_per_model). It does not import mpi4py and constructs no
#       communicators, so it is fully unit-testable without launching MPI.
#       Communicator construction from a plan lives in
#       icesee_mpi_parallel_manager.py.
# ==============================================================================
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass(frozen=True)
class ResourcePlan:
    """A deterministic hierarchical ensemble/model execution plan.

    ``COMM_WORLD`` (size ``world_size``) is partitioned into
    ``num_model_groups`` block-contiguous groups of ``ranks_per_model``
    ranks each (group ``g`` = world ranks
    ``[g*ranks_per_model, (g+1)*ranks_per_model)``), plus
    ``spare_ranks`` leftover world ranks belonging to no group. Over
    ``num_rounds`` rounds, group ``g`` handles ensemble member
    ``schedule[r][g]`` in round ``r`` (``None`` if that group has no member
    in that round -- possible only in the final round when ``Nens`` is not
    a multiple of ``num_model_groups``). Every member in
    ``range(nens)`` appears in ``schedule`` exactly once.
    """

    world_size: int
    nens: int
    ranks_per_model: int
    num_model_groups: int
    num_rounds: int
    active_ranks: int
    spare_ranks: int
    schedule: Tuple[Tuple[Optional[int], ...], ...]

    def group_id(self, rank: int) -> Optional[int]:
        """This world rank's model-group id, or ``None`` if it is spare."""
        if rank < 0 or rank >= self.world_size:
            raise ValueError(f"rank {rank} out of range [0, {self.world_size})")
        if rank >= self.active_ranks:
            return None
        return rank // self.ranks_per_model

    def rank_in_group(self, rank: int) -> Optional[int]:
        """This world rank's rank within its model group, or ``None`` if spare."""
        if rank < 0 or rank >= self.world_size:
            raise ValueError(f"rank {rank} out of range [0, {self.world_size})")
        if rank >= self.active_ranks:
            return None
        return rank % self.ranks_per_model

    def is_spare(self, rank: int) -> bool:
        return self.group_id(rank) is None

    def member_for(self, round_id: int, group: int) -> Optional[int]:
        if round_id < 0 or round_id >= self.num_rounds:
            raise ValueError(f"round {round_id} out of range [0, {self.num_rounds})")
        if group < 0 or group >= self.num_model_groups:
            raise ValueError(f"group {group} out of range [0, {self.num_model_groups})")
        return self.schedule[round_id][group]

    def members_covered(self):
        """Every member id that appears anywhere in the schedule, for
        coverage tests (no missing / no duplicate members)."""
        seen = []
        for row in self.schedule:
            for member in row:
                if member is not None:
                    seen.append(member)
        return seen

    def summary(self) -> str:
        """One-line diagnostic, e.g. for a rank-0-only print at startup."""
        return (
            f"ICESEE topology: world={self.world_size}, ensembles={self.nens}, "
            f"ranks/model={self.ranks_per_model}, groups={self.num_model_groups}, "
            f"rounds={self.num_rounds}, spare={self.spare_ranks}"
        )

    def describe_all_ranks(self) -> str:
        """Multi-line, one-row-per-world-rank topology dump: group id, rank
        within group, spare status, and the member each rank's group
        handles in every round. Intended for rank-0-only, verbose/debug-
        gated output (``icesee_kwargs.get("verbose")``) -- not printed by
        default, to avoid per-rank spam at any real world size.
        """
        lines = [self.summary(), "  rank | group | rank_in_group | members_by_round"]
        for rank in range(self.world_size):
            group = self.group_id(rank)
            if group is None:
                lines.append(f"  {rank:4d} | spare | spare          | -")
                continue
            members = [self.member_for(r, group) for r in range(self.num_rounds)]
            lines.append(
                f"  {rank:4d} | {group:5d} | {self.rank_in_group(rank):14d} | {members}"
            )
        return "\n".join(lines)


def plan_resources(
    world_size: int,
    nens: int,
    requested_ranks_per_model=None,
    max_model_groups=None,
) -> ResourcePlan:
    """Compute a deterministic hierarchical resource plan.

    ``max_model_groups`` (Stage 4D.2): an optional, independent cap on
    ``num_model_groups``, applied after the usual
    ``min(raw_groups, nens)`` computation. ``None`` (the default) changes
    nothing -- every existing caller/configuration is unaffected. This
    exists for models like ISSM where a "model group" is a persistent,
    expensive external worker (one MATLAB server, in turn launching its
    own separately-scheduled nested MPI solve via ``model_nprocs``) and a
    caller may want to deliberately bound how many run concurrently
    (``issm_server_count`` in ISSM's own config) even when more world
    ranks are available -- the extra ranks become ordinary spare ranks,
    reusing this module's existing spare-rank machinery unchanged rather
    than inventing a second one. A value below 1 is rejected the same way
    an invalid ``ranks_per_model`` is.

    ``requested_ranks_per_model`` semantics:

    - ``None`` (key absent from configuration) -- **legacy/default**.
      Reproduces the pre-hierarchical-topology behavior exactly:
      ``ranks_per_model = 1`` whenever ``nens >= world_size`` (the regime
      every currently shipped configuration, including ISSM, actually
      runs in -- ISSM's own integration hard-requires ``nens >= world_size``
      and must see ``ranks_per_model == 1`` unchanged). When
      ``world_size > nens``, applies the same conservative floor-division
      policy as ``"auto"`` below (this differs from the old code's
      *strided*, sometimes-uneven-group-size assignment for that regime --
      see ADR 0003/0005 and the Stage 4 CHANGELOG entry for why that
      change is safe: no shipped configuration set or depended on the old
      formula's specific group sizes or rank layout, only on every member
      being processed correctly, which this plan still guarantees).
    - ``"auto"`` -- ``ranks_per_model = max(1, world_size // nens)``, the
      same conservative floor-division policy, explicitly requested rather
      than implied by omission.
    - a positive integer -- the user requests exactly that many ranks per
      model group.

    Raises ``ValueError`` for any invalid/insufficient configuration
    (fail fast, at the planning boundary, rather than deep inside MPI
    setup) -- ``world_size`` or ``nens`` below 1, a requested
    ``ranks_per_model`` below 1, or ``world_size < ranks_per_model`` (not
    enough ranks to form even one model group).
    """
    world_size = int(world_size)
    nens = int(nens)
    if world_size < 1:
        raise ValueError(f"world_size must be >= 1, got {world_size}")
    if nens < 1:
        raise ValueError(f"nens must be >= 1, got {nens}")

    if requested_ranks_per_model is None:
        ranks_per_model = 1 if nens >= world_size else max(1, world_size // nens)
    elif isinstance(requested_ranks_per_model, str) and requested_ranks_per_model.strip().lower() == "auto":
        ranks_per_model = max(1, world_size // nens)
    else:
        ranks_per_model = int(requested_ranks_per_model)
        if ranks_per_model < 1:
            raise ValueError(
                f"ranks_per_model must be >= 1, got {ranks_per_model}"
            )

    if world_size < ranks_per_model:
        raise ValueError(
            f"world_size ({world_size}) is smaller than the requested "
            f"ranks_per_model ({ranks_per_model}); not enough ranks to "
            "form even one model group. Launch with more ranks or reduce "
            "ranks_per_model."
        )

    raw_groups = world_size // ranks_per_model
    # Never form more groups than there are members to ever assign -- a
    # group with no member in any round would just be permanently idle
    # ranks wearing a group identity instead of being cleanly spare.
    num_model_groups = min(raw_groups, nens)

    if max_model_groups is not None:
        max_model_groups = int(max_model_groups)
        if max_model_groups < 1:
            raise ValueError(
                f"max_model_groups must be >= 1, got {max_model_groups}"
            )
        num_model_groups = min(num_model_groups, max_model_groups)

    if num_model_groups < 1:
        raise ValueError(
            f"Computed 0 model groups for world_size={world_size}, "
            f"ranks_per_model={ranks_per_model}, nens={nens}; this should "
            "be unreachable given the checks above."
        )

    active_ranks = num_model_groups * ranks_per_model
    spare_ranks = world_size - active_ranks
    num_rounds = math.ceil(nens / num_model_groups)

    schedule = []
    member_id = 0
    for _round in range(num_rounds):
        row = []
        for _group in range(num_model_groups):
            if member_id < nens:
                row.append(member_id)
                member_id += 1
            else:
                row.append(None)
        schedule.append(tuple(row))

    return ResourcePlan(
        world_size=world_size,
        nens=nens,
        ranks_per_model=ranks_per_model,
        num_model_groups=num_model_groups,
        num_rounds=num_rounds,
        active_ranks=active_ranks,
        spare_ranks=spare_ranks,
        schedule=tuple(schedule),
    )
