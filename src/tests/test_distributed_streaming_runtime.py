# ==============================================================================
# @des: Fast, Firedrake-free regression tests for the bounded-memory
# StreamingNativeDistributedMemberPool (Option B) and its
# InactiveMemberStore, using a synthetic fake adapter -- the generic
# orchestration logic (activate/deactivate/forecast/observe/restore) is
# entirely model-agnostic and should not need real Firedrake to validate.
# Real-Firedrake correctness (reactivation round-trip exactness, scientific
# equivalence against the persistent pool) is covered separately in
# test_icepack_member_streaming.py.
# ==============================================================================
from __future__ import annotations

import types

import numpy as np
import pytest

from ICESEE.src.parallelization.distributed_fields import DistributedFieldRegistry
from ICESEE.src.parallelization.distributed_member_store import MemoryInactiveMemberStore
from ICESEE.src.parallelization.distributed_native_cycle import (
    NativeObservationBatch,
    run_native_global_analysis_cycle,
)
from ICESEE.src.parallelization.distributed_native_runtime import NativeDistributedMember
from ICESEE.src.parallelization.distributed_streaming_runtime import (
    StreamingNativeDistributedMemberPool,
    extract_rows_from_packed,
    partition_owned_rows_from_layout,
    run_native_store_streaming_analysis_cycle,
)


class _FakeField:
    def __init__(self, values):
        self.name = "state"
        values = np.asarray(values, dtype=float)
        self.global_size = values.size
        self.owned_start = 0
        self.owned_stop = values.size
        self._values = values.copy()

    def read_owned(self):
        return self._values

    def write_owned(self, values):
        self._values[:] = np.asarray(values, dtype=float)

    def synchronize_ghosts(self):
        return None


def _topology():
    return types.SimpleNamespace(spatial_comm=None, spatial_ranks=1, spatial_rank=0)


class _FakeAdapter:
    def __init__(self, size=3):
        self.size = size
        self.init_calls = []
        self.reactivate_calls = []
        self.forecast_calls = []

    def initialize_native_member(self, member_id, *, topology, icesee_kwargs):
        self.init_calls.append(int(member_id))
        field = _FakeField([float(member_id)] * self.size)
        registry = DistributedFieldRegistry([field], layout_id="fake-streaming-v1")
        return NativeDistributedMember(int(member_id), registry, {"member_id": member_id})

    def reactivate_native_member(self, member_id, packed_state, *, topology, icesee_kwargs):
        self.reactivate_calls.append(int(member_id))
        field = _FakeField(packed_state)
        registry = DistributedFieldRegistry([field], layout_id="fake-streaming-v1")
        return NativeDistributedMember(int(member_id), registry, {"member_id": member_id})

    def forecast_native_member(self, member, timestep, *, topology, icesee_kwargs):
        self.forecast_calls.append((member.member_id, int(timestep)))
        state = member.pack_owned()
        member.unpack_owned(state + 1.0)

    def observe_native_member(self, member, rows, *, topology, icesee_kwargs):
        return member.fields.observe_owned_rows(np.asarray(rows, dtype=np.int64))

    def finalize_native_analysis(self, member, forecast_owned, timestep, *, topology, icesee_kwargs):
        return None


class _FakeAdapterWithoutStreaming:
    def initialize_native_member(self, member_id, *, topology, icesee_kwargs):
        field = _FakeField([0.0])
        registry = DistributedFieldRegistry([field], layout_id="no-stream-v1")
        return NativeDistributedMember(int(member_id), registry, {})

    def forecast_native_member(self, member, timestep, *, topology, icesee_kwargs):
        return None


def test_memory_inactive_member_store_put_get_release():
    store = MemoryInactiveMemberStore()
    assert not store.has(0)
    store.put(0, np.array([1.0, 2.0, 3.0]))
    assert store.has(0)
    np.testing.assert_array_equal(store.get(0), [1.0, 2.0, 3.0])
    assert store.total_bytes() == 24
    store.release(0)
    assert not store.has(0)
    with pytest.raises(KeyError):
        store.get(0)


def test_pool_requires_reactivate_native_member():
    with pytest.raises(TypeError, match="reactivate_native_member"):
        StreamingNativeDistributedMemberPool(
            _FakeAdapterWithoutStreaming(), _topology(), {}, (0,)
        )


def test_pool_builds_each_member_once_and_releases_immediately():
    adapter = _FakeAdapter()
    pool = StreamingNativeDistributedMemberPool(adapter, _topology(), {}, (0, 1, 2))
    assert adapter.init_calls == [0, 1, 2]
    assert pool._active is None  # nothing resident after construction
    assert pool.member_ids == (0, 1, 2)


def test_pool_forecast_activates_and_deactivates_every_member():
    adapter = _FakeAdapter()
    pool = StreamingNativeDistributedMemberPool(adapter, _topology(), {}, (0, 1))
    result = pool.forecast_with_adapter(5, adapter, topology=_topology(), icesee_kwargs={})
    assert set(result.member_ids) == {0, 1}
    np.testing.assert_array_equal(result.members[0], [1.0, 1.0, 1.0])
    np.testing.assert_array_equal(result.members[1], [2.0, 2.0, 2.0])
    assert adapter.reactivate_calls == [0, 1]
    assert adapter.forecast_calls == [(0, 5), (1, 5)]
    assert pool._active is None  # deactivated after the last member too


def test_pool_observe_does_not_reactivate():
    """Observing an inactive member must read directly from its packed
    state -- no reactivation needed, unlike forecasting."""
    adapter = _FakeAdapter()
    pool = StreamingNativeDistributedMemberPool(adapter, _topology(), {}, (0, 1))
    pool.forecast_with_adapter(0, adapter, topology=_topology(), icesee_kwargs={})
    calls_before = len(adapter.reactivate_calls)
    observed = pool.observe_with_adapter(
        np.array([0, 2], dtype=np.int64), adapter, topology=_topology(), icesee_kwargs={}
    )
    assert len(adapter.reactivate_calls) == calls_before  # unchanged
    np.testing.assert_array_equal(observed[0], [1.0, 1.0])
    np.testing.assert_array_equal(observed[1], [2.0, 2.0])


def test_pool_restore_and_finalize_updates_stored_state():
    adapter = _FakeAdapter()
    pool = StreamingNativeDistributedMemberPool(adapter, _topology(), {}, (0, 1))
    local_forecast = pool.forecast_with_adapter(0, adapter, topology=_topology(), icesee_kwargs={})
    from ICESEE.src.parallelization.distributed_runtime import LocalMemberEnsemble
    analyzed = LocalMemberEnsemble(
        pool.layout, {0: np.array([9.0, 9.0, 9.0]), 1: np.array([8.0, 8.0, 8.0])}
    )
    pool.restore_and_finalize(analyzed, local_forecast, 0, adapter, topology=_topology(), icesee_kwargs={})
    snap = pool.snapshot_owned()
    np.testing.assert_array_equal(snap.members[0], [9.0, 9.0, 9.0])
    np.testing.assert_array_equal(snap.members[1], [8.0, 8.0, 8.0])
    assert pool._active is None


def test_pool_route_observation_ids_partitions_by_layout():
    adapter = _FakeAdapter()
    pool = StreamingNativeDistributedMemberPool(adapter, _topology(), {}, (0,))
    positions, local_ids = pool.route_observation_ids(
        np.array([0, 1, 2], dtype=np.int64), adapter, topology=_topology(), icesee_kwargs={}
    )
    np.testing.assert_array_equal(positions, [0, 1, 2])
    np.testing.assert_array_equal(local_ids, [0, 1, 2])


def test_extract_rows_from_packed_matches_registry_observe_owned_rows():
    field = _FakeField([10.0, 20.0, 30.0])
    registry = DistributedFieldRegistry([field], layout_id="check-v1")
    packed = registry.pack_owned()
    from_registry = registry.observe_owned_rows(np.array([0, 2], dtype=np.int64))
    from_packed = extract_rows_from_packed(registry.layout, packed, np.array([0, 2], dtype=np.int64))
    np.testing.assert_array_equal(from_registry, from_packed)


def test_partition_owned_rows_from_layout_matches_registry():
    field = _FakeField([1.0, 2.0, 3.0])
    registry = DistributedFieldRegistry([field], layout_id="check-v2")
    rows = np.array([0, 1, 2], dtype=np.int64)
    expected = registry.partition_owned_rows(rows)
    actual = partition_owned_rows_from_layout(registry.layout, rows)
    np.testing.assert_array_equal(expected[0], actual[0])
    np.testing.assert_array_equal(expected[1], actual[1])


class _Comm:
    """Minimal single-rank fake MPI communicator (matches the existing
    pattern in test_distributed_native_cycle.py)."""

    def Get_rank(self):
        return 0

    def Get_size(self):
        return 1

    def allgather(self, value):
        return [value]

    def bcast(self, value, root=0):
        return value

    def Allreduce(self, source, target):
        target[...] = source


def _single_rank_topology():
    comm = _Comm()
    return types.SimpleNamespace(
        world=comm, spatial_comm=comm, ensemble_comm=comm,
        world_rank=0, world_size=1, ensemble_slot=0, spatial_rank=0,
        ensemble_groups=1, spatial_ranks=1, is_ensemble_root=True,
    )


class _CycleFakeAdapter:
    """Fake adapter exercising a full analysis cycle: deterministic
    forecast (add timestep), identity observation, no-op finalize."""

    def initialize_native_member(self, member_id, *, topology, icesee_kwargs):
        field = _FakeField(np.arange(5.0) + 10.0 * member_id)
        registry = DistributedFieldRegistry([field], layout_id="cycle-fake-v1")
        return NativeDistributedMember(int(member_id), registry, {})

    def reactivate_native_member(self, member_id, packed_state, *, topology, icesee_kwargs):
        field = _FakeField(packed_state)
        registry = DistributedFieldRegistry([field], layout_id="cycle-fake-v1")
        return NativeDistributedMember(int(member_id), registry, {})

    def forecast_native_member(self, member, timestep, *, topology, icesee_kwargs):
        member.unpack_owned(member.pack_owned() + float(timestep))

    def observe_native_member(self, member, observation_rows, *, topology, icesee_kwargs):
        return member.pack_owned()[observation_rows]

    def finalize_native_analysis(self, member, forecast_owned, timestep, *, topology, icesee_kwargs):
        return None


def test_store_streaming_cycle_matches_original_cycle_bit_for_bit():
    """The store-streaming analysis path (transform_members_via_store,
    run_native_store_streaming_analysis_cycle) must produce EXACTLY the
    same analyzed state as the original transform_local_members path, given
    the same pool type (StreamingNativeDistributedMemberPool held constant
    -- this isolates the transform-step change from the already-separately-
    verified pool-type change)."""

    adapter = _CycleFakeAdapter()
    topology = _single_rank_topology()
    icesee_kwargs = {"Nens": 2}

    pool_old = StreamingNativeDistributedMemberPool(adapter, topology, icesee_kwargs, (0, 1))
    observed = np.array([2.0, 7.0])
    batch_old = NativeObservationBatch(np.array([0, 4]), observed)
    result_old = run_native_global_analysis_cycle(
        pool_old, adapter, 1, [batch_old], number_of_batches=1,
        topology=topology, icesee_kwargs=icesee_kwargs, error_mode="legacy_prior_anomalies",
        state_row_chunk_size=2,
    )

    pool_new = StreamingNativeDistributedMemberPool(adapter, topology, icesee_kwargs, (0, 1))
    batch_new = NativeObservationBatch(np.array([0, 4]), observed)
    result_new = run_native_store_streaming_analysis_cycle(
        pool_new, adapter, 1, [batch_new], number_of_batches=1,
        topology=topology, icesee_kwargs=icesee_kwargs, error_mode="legacy_prior_anomalies",
        state_row_chunk_size=2,
    )

    np.testing.assert_array_equal(result_old.transform, result_new.transform)
    for member_id in (0, 1):
        np.testing.assert_array_equal(
            result_old.local_analysis.members[member_id],
            result_new.local_analysis.members[member_id],
        )
    # And directly from each pool's own store, post-cycle.
    np.testing.assert_array_equal(pool_old.snapshot_owned().members[0], pool_new.snapshot_owned().members[0])
    np.testing.assert_array_equal(pool_old.snapshot_owned().members[1], pool_new.snapshot_owned().members[1])


def test_store_streaming_cycle_skips_transform_pass_when_no_batches_scheduled():
    """The Phase-1 operation-count optimization: number_of_batches=0 must
    skip transform_members_via_store entirely (a provable no-op, since the
    ensemble transform is exactly eye(Ne) with no observation contributions)
    while still calling finalize_all_from_store -- and the result must be
    bit-for-bit identical to running the SAME cycle WITHOUT the
    optimization (i.e. to the original run_native_global_analysis_cycle
    reference path on a persistent pool, given no batches either)."""

    adapter = _CycleFakeAdapter()
    topology = _single_rank_topology()
    icesee_kwargs = {"Nens": 2}

    from ICESEE.src.parallelization.distributed_native_runtime import NativeDistributedMemberPool

    reference_pool = NativeDistributedMemberPool(
        {
            0: adapter.initialize_native_member(0, topology=topology, icesee_kwargs=icesee_kwargs),
            1: adapter.initialize_native_member(1, topology=topology, icesee_kwargs=icesee_kwargs),
        }
    )
    result_reference = run_native_global_analysis_cycle(
        reference_pool, adapter, 1, [], number_of_batches=0,
        topology=topology, icesee_kwargs=icesee_kwargs, error_mode="legacy_prior_anomalies",
        state_row_chunk_size=2,
    )

    streaming_pool = StreamingNativeDistributedMemberPool(adapter, topology, icesee_kwargs, (0, 1))
    result_streaming = run_native_store_streaming_analysis_cycle(
        streaming_pool, adapter, 1, [], number_of_batches=0,
        topology=topology, icesee_kwargs=icesee_kwargs, error_mode="legacy_prior_anomalies",
        state_row_chunk_size=2,
    )

    np.testing.assert_array_equal(result_reference.transform, np.eye(2))
    np.testing.assert_array_equal(result_streaming.transform, np.eye(2))
    for member_id in (0, 1):
        np.testing.assert_array_equal(
            result_reference.local_analysis.members[member_id],
            result_streaming.local_analysis.members[member_id],
        )
