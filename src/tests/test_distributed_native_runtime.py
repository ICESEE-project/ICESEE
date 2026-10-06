from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytest

from src.parallelization.distributed_fields import DistributedFieldRegistry
from src.parallelization.distributed_native_runtime import (
    NativeDistributedMember,
    NativeDistributedMemberPool,
    NativeObservationShard,
    initialize_native_member_pool,
    validate_native_distributed_adapter,
)
from src.parallelization.distributed_runtime import LocalMemberEnsemble


@dataclass
class _Field:
    name: str
    global_size: int
    owned_start: int
    owned_stop: int
    values: np.ndarray
    syncs: int = 0

    def read_owned(self):
        return self.values

    def write_owned(self, values):
        self.values[:] = values

    def synchronize_ghosts(self):
        self.syncs += 1


def _member(member_id):
    h = _Field("h", 12, 3, 6, np.array([1, 2, 3.0]) + member_id)
    u = _Field("u", 8, 2, 4, np.array([10, 11.0]) + member_id)
    registry = DistributedFieldRegistry([h, u], layout_id="native-pool-v1")
    return NativeDistributedMember(member_id, registry, {"h": h, "u": u})


def test_native_pool_round_trip_and_forecast_keep_local_shapes():
    pool = NativeDistributedMemberPool({0: _member(0), 2: _member(2)})
    initial = pool.snapshot_owned(dtype=np.float32)
    assert initial.layout.global_size == 20
    assert initial.layout.owned_size == 5
    assert initial.members[0].dtype == np.float32

    updated = LocalMemberEnsemble(
        pool.layout,
        {member_id: values + 5 for member_id, values in initial.members.items()},
    )
    pool.restore_owned(updated)
    np.testing.assert_array_equal(pool.member(0).pack_owned(), [6, 7, 8, 15, 16])

    def forecast(member, timestep, **kwargs):
        state = member.pack_owned()
        member.unpack_owned(state + timestep)

    forecasted = pool.forecast(3, forecast, topology=object(), icesee_kwargs={})
    np.testing.assert_array_equal(forecasted.members[0], [9, 10, 11, 18, 19])
    assert forecasted.members[0].size == pool.layout.owned_size


def test_native_pool_rejects_global_return_and_incompatible_restore():
    pool = NativeDistributedMemberPool({0: _member(0)})
    with pytest.raises(TypeError, match="must return None"):
        pool.forecast(
            1,
            lambda *args, **kwargs: np.zeros(20),
            topology=object(),
            icesee_kwargs={},
        )

    other = _member(1)
    with pytest.raises(ValueError, match="member IDs"):
        pool.restore_owned(LocalMemberEnsemble(pool.layout, {1: other.pack_owned()}))


class _NativeAdapter:
    def initialize_native_member(self, member_id, **kwargs):
        return _member(member_id)

    def forecast_native_member(self, member, timestep, **kwargs):
        member.unpack_owned(member.pack_owned() + timestep)

    def observe_native_member(self, member, observation_rows, **kwargs):
        return member.pack_owned()[observation_rows]

    def finalize_native_analysis(self, member, forecast_owned, timestep, **kwargs):
        member.model_context["finalized_from"] = forecast_owned.copy()


def test_native_adapter_forecast_observe_finalize_and_storage_plan():
    adapter = _NativeAdapter()
    validate_native_distributed_adapter(adapter)
    pool = NativeDistributedMemberPool({0: _member(0), 2: _member(2)})

    before = pool.snapshot_owned()
    after = pool.forecast_with_adapter(
        2, adapter, topology=object(), icesee_kwargs={}
    )
    observations = pool.observe_with_adapter(
        np.array([0, 4]), adapter, topology=object(), icesee_kwargs={}
    )
    np.testing.assert_array_equal(observations[0], after.members[0][[0, 4]])

    analysis = LocalMemberEnsemble(
        pool.layout,
        {member_id: values - 1 for member_id, values in after.members.items()},
    )
    pool.restore_and_finalize(
        analysis,
        before,
        2,
        adapter,
        topology=object(),
        icesee_kwargs={},
    )
    np.testing.assert_array_equal(pool.member(0).pack_owned(), analysis.members[0])
    np.testing.assert_array_equal(
        pool.member(0).model_context["finalized_from"], before.members[0]
    )

    plan = pool.storage_plan(dtype=np.float32)
    assert plan["owned_entries_per_member"] == 5
    assert plan["global_entries_per_member"] == 20
    assert plan["owned_snapshot_bytes"] == 2 * 5 * 4
    assert plan["global_snapshot_bytes_avoided"] == 2 * 20 * 4


class _GlobalRowNativeAdapter(_NativeAdapter):
    def observe_native_member(self, member, observation_rows, **kwargs):
        return member.fields.observe_owned_rows(observation_rows)


def test_native_pool_routes_irregular_global_state_observations_locally():
    adapter = _GlobalRowNativeAdapter()
    pool = NativeDistributedMemberPool({0: _member(0), 2: _member(2)})

    shard = pool.observe_partitioned_with_adapter(
        np.asarray([0, 4, 14, 3, 19], dtype=np.int64),
        adapter,
        topology=object(),
        icesee_kwargs={},
    )
    assert isinstance(shard, NativeObservationShard)
    np.testing.assert_array_equal(shard.canonical_positions, [1, 2, 3])
    np.testing.assert_array_equal(shard.observation_ids, [4, 14, 3])
    np.testing.assert_array_equal(shard.member_values[0], [2, 10, 1])
    np.testing.assert_array_equal(shard.member_values[2], [4, 12, 3])


class _NonlinearObservationAdapter(_NativeAdapter):
    def partition_native_observations(self, member, observation_ids, **kwargs):
        # Model-defined observation IDs need not be packed state rows.
        return np.asarray([1, 3]), observation_ids[[1, 3]]

    def observe_native_member(self, member, observation_rows, **kwargs):
        return observation_rows.astype(float) + member.member_id


def test_native_pool_accepts_model_defined_observation_partition():
    adapter = _NonlinearObservationAdapter()
    pool = NativeDistributedMemberPool({0: _member(0), 2: _member(2)})
    shard = pool.observe_partitioned_with_adapter(
        np.asarray([100, 101, 102, 103], dtype=np.int64),
        adapter,
        topology=object(),
        icesee_kwargs={},
    )
    np.testing.assert_array_equal(shard.canonical_positions, [1, 3])
    np.testing.assert_array_equal(shard.observation_ids, [101, 103])
    np.testing.assert_array_equal(shard.member_values[2], [103, 105])


class _Collective:
    def allgather(self, value):
        return [value]


@dataclass
class _NativeTopology:
    ensemble_groups: int = 2
    ensemble_slot: int = 0
    spatial_comm: object = field(default_factory=_Collective)


def _complete_member(member_id):
    h = _Field("h", 3, 0, 3, np.array([1, 2, 3.0]) + member_id)
    registry = DistributedFieldRegistry([h], layout_id="complete-native-v1")
    return NativeDistributedMember(member_id, registry, {"h": h})


class _InitializingAdapter(_NativeAdapter):
    def __init__(self):
        self.initialized = []

    def initialize_native_member(self, member_id, **kwargs):
        self.initialized.append(member_id)
        return _complete_member(member_id)


def test_initialize_native_pool_schedules_members_without_complete_ensemble():
    adapter = _InitializingAdapter()
    pool = initialize_native_member_pool(
        adapter,
        _NativeTopology(),
        {"Nens": 5},
    )
    assert pool.member_ids == (0, 2, 4)
    assert adapter.initialized == [0, 2, 4]
