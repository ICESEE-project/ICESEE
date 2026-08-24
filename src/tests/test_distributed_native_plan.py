from __future__ import annotations

from dataclasses import dataclass

import pytest

from src.parallelization.distributed_native_plan import (
    build_native_execution_plan,
    resolve_native_distributed_adapter,
)


class _World:
    def __init__(self, gathered=None):
        self.gathered = gathered

    def allgather(self, value):
        return list(self.gathered) if self.gathered is not None else [value]


@dataclass
class _Topology:
    ensemble_groups: int = 2
    spatial_ranks: int = 4
    world: object = None

    def __post_init__(self):
        if self.world is None:
            self.world = _World()


class _Adapter:
    def initialize_native_member(self, member_id, **kwargs):
        raise NotImplementedError

    def forecast_native_member(self, member, timestep, **kwargs):
        return None

    def observe_native_member(self, member, observation_rows, **kwargs):
        raise NotImplementedError

    def finalize_native_analysis(self, member, forecast_owned, timestep, **kwargs):
        return None

    def distributed_global_state_size(self, **kwargs):
        return 8000


class _LocalAdapter(_Adapter):
    def distributed_analysis_targets(self, **kwargs):
        return ()

    def distributed_observation_coordinates(self, **kwargs):
        return ()


class _InversionAdapter(_LocalAdapter):
    def inverse_native_member(self, member, timestep, **kwargs):
        return None


def test_plan_resolves_direct_adapter_and_applies_bounded_memory_estimate():
    adapter = _Adapter()
    plan = build_native_execution_plan(
        {
            "distributed_model_adapter": adapter,
            "Nens": 8,
            "mode3_max_rank_memory_gib": 1.0,
            "observation_error_mode": "stochastic_R",
        },
        _Topology(),
    )
    assert plan.adapter is adapter
    assert plan.analysis_kind == "global_stochastic"
    assert plan.error_mode == "stochastic_r"
    assert plan.memory.global_rows == 8000
    assert plan.memory.local_members == 4


def test_plan_supports_context_factory_without_model_specific_imports():
    adapter = _LocalAdapter()
    seen = {}

    def factory(*, topology, icesee_kwargs):
        seen["topology"] = topology
        seen["context"] = icesee_kwargs
        return adapter

    context = {
        "distributed_model_adapter_factory": factory,
        "Nens": 4,
        "local_analysis": True,
        "mode3_global_state_rows": 400,
    }
    topology = _Topology()
    plan = build_native_execution_plan(context, topology)
    assert plan.analysis_kind == "grouped_local_stochastic"
    assert seen["topology"] is topology
    assert seen["context"] is context
    assert "distributed_analysis_targets" in plan.capabilities


def test_plan_rejects_local_analysis_without_model_metadata_callbacks():
    with pytest.raises(TypeError, match="localization callbacks"):
        build_native_execution_plan(
            {
                "distributed_model_adapter": _Adapter(),
                "Nens": 4,
                "local_analysis": True,
                "mode3_global_state_rows": 400,
            },
            _Topology(),
        )


def test_plan_requires_distributed_inversion_callback_when_enabled():
    with pytest.raises(TypeError, match="inverse_native_member"):
        build_native_execution_plan(
            {
                "distributed_model_adapter": _LocalAdapter(),
                "Nens": 4,
                "inversion_flag": True,
                "mode3_global_state_rows": 400,
            },
            _Topology(),
        )
    plan = build_native_execution_plan(
        {
            "distributed_model_adapter": _InversionAdapter(),
            "Nens": 4,
            "inversion_flag": True,
            "mode3_global_state_rows": 400,
        },
        _Topology(),
    )
    assert plan.inversion_enabled


def test_plan_fails_before_allocation_when_memory_budget_is_unsafe():
    with pytest.raises(MemoryError, match="increase spatial_ranks"):
        build_native_execution_plan(
            {
                "distributed_model_adapter": _Adapter(),
                "Nens": 40,
                "mode3_global_state_rows": 5_000_000_000,
                "mode3_max_rank_memory_gib": 1.0,
            },
            _Topology(ensemble_groups=2, spatial_ranks=2),
        )


def test_adapter_capability_mismatch_across_ranks_is_rejected():
    adapter = _Adapter()
    identity = (type(adapter).__module__, type(adapter).__qualname__, ())
    topology = _Topology(world=_World([identity, ("other", "Adapter", ())]))
    with pytest.raises(RuntimeError, match="inconsistent adapter"):
        resolve_native_distributed_adapter(
            {"distributed_model_adapter": adapter}, topology
        )
