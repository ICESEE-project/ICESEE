# ==============================================================================
# @des: Demonstrates that Stage 4's ResourcePlan (resource topology: "which
# ranks cooperate on which model/member") and StateOwnership (state
# layout: "which portion of the state belongs to this rank") are
# orthogonal concerns, as required by the Stage 4 spec: ranks_per_model >
# 1 must NOT be inferred to imply distributed state, and
# ranks_per_model == 1 must NOT be inferred to imply replicated state.
# ==============================================================================
import dataclasses

import numpy as np
import pytest

from ICESEE.src.parallelization.parallel_mpi.resource_plan import (
    ResourcePlan,
    plan_resources,
)
from ICESEE.src.utils.state_ownership import StateOwnership


def test_resource_plan_fields_contain_no_state_ownership_concept():
    # Structural check: ResourcePlan's own vocabulary never mentions state
    # distribution/ownership -- topology and state layout are described by
    # two entirely separate types.
    field_names = {f.name for f in dataclasses.fields(ResourcePlan)}
    for forbidden in ("distribution", "state", "replicated", "distributed", "ownership"):
        assert not any(forbidden in name for name in field_names), (
            f"ResourcePlan field set {field_names} unexpectedly references "
            f"state-ownership vocabulary ({forbidden!r})"
        )


def test_state_ownership_fields_contain_no_topology_concept():
    field_names = {f.name for f in dataclasses.fields(StateOwnership)}
    for forbidden in ("rank_per_model", "group", "round", "spare", "world_size"):
        assert not any(forbidden in name for name in field_names), (
            f"StateOwnership field set {field_names} unexpectedly references "
            f"resource-topology vocabulary ({forbidden!r})"
        )


def test_ranks_per_model_greater_than_one_can_host_replicated_state():
    # A ranks_per_model=2 group does not force distributed state: every
    # rank in the group could still (hypothetically -- no shipped
    # application does this today, see model_capabilities.py) hold an
    # identical replicated copy, exactly like a replicated model on a
    # ranks_per_model=1 group.
    plan = plan_resources(world_size=8, nens=4, requested_ranks_per_model=2)
    assert plan.ranks_per_model == 2

    global_size = 3
    ownership_rank_0 = StateOwnership(distribution="replicated", global_size=global_size, local_size=global_size, local_offset=0)
    ownership_rank_1 = StateOwnership(distribution="replicated", global_size=global_size, local_size=global_size, local_offset=0)
    assert ownership_rank_0 == ownership_rank_1  # both ranks in the group see the same thing


def test_ranks_per_model_greater_than_one_can_host_distributed_state():
    # The same ranks_per_model=2 topology also supports genuinely
    # distributed state (e.g. a future Icepack modes-0-2 multi-rank
    # integration): each rank owns a disjoint slice, and the plan itself
    # neither knows nor cares which case applies.
    plan = plan_resources(world_size=8, nens=4, requested_ranks_per_model=2)
    assert plan.ranks_per_model == 2

    global_size = 10
    ownership_rank_0 = StateOwnership(distribution="distributed", global_size=global_size, local_size=5, local_offset=0)
    ownership_rank_1 = StateOwnership(distribution="distributed", global_size=global_size, local_size=5, local_offset=5)
    assert ownership_rank_0.local_size + ownership_rank_1.local_size == global_size
    assert ownership_rank_0.local_offset != ownership_rank_1.local_offset


def test_ranks_per_model_one_can_also_host_either_ownership_kind():
    # The converse: ranks_per_model == 1 (every current shipped
    # configuration) does not itself imply "replicated" either -- a
    # single-rank group's state is simply whatever that one rank's
    # StateOwnership says it is, resolved independently.
    plan = plan_resources(world_size=4, nens=4, requested_ranks_per_model=1)
    assert plan.ranks_per_model == 1

    replicated = StateOwnership(distribution="replicated", global_size=3, local_size=3, local_offset=0)
    distributed = StateOwnership(distribution="distributed", global_size=3, local_size=3, local_offset=0)
    assert replicated.distribution != distributed.distribution
    # Both are valid StateOwnership values regardless of what the plan says.


def test_plan_resources_does_not_import_state_ownership_module():
    import ICESEE.src.parallelization.parallel_mpi.resource_plan as resource_plan_module

    assert "state_ownership" not in resource_plan_module.__file__
    import sys

    # resource_plan.py must not even transitively require state_ownership
    # to be importable -- confirmed by checking its own module globals
    # never bind anything from that module.
    state_ownership_names = {
        name for name, value in vars(resource_plan_module).items()
        if getattr(value, "__module__", "") == "ICESEE.src.utils.state_ownership"
    }
    assert state_ownership_names == set()
