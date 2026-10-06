from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from src.parallelization.distributed_analysis import (
    StochasticAnalysisProducts,
    ensemble_transform_from_products,
)
from src.parallelization.distributed_fields import DistributedFieldRegistry
from src.parallelization.distributed_native_cycle import (
    NativeCheckpointRequest,
    NativeLocalAnalysisRequest,
    NativeObservationBatch,
    _validate_batch_alignment,
    restore_native_pool_from_checkpoint,
    run_native_grouped_local_analysis_cycle,
    run_native_global_analysis_cycle,
)
from src.parallelization.distributed_adapter import (
    DistributedAnalysisTargets,
    DistributedObservationCollection,
    DistributedObservationCoordinates,
)
from src.parallelization.distributed_native_runtime import (
    NativeDistributedMember,
    NativeDistributedMemberPool,
    apply_native_memberwise_inversion,
)


class _Comm:
    def Get_rank(self):
        return 0

    def Get_size(self):
        return 1

    def allgather(self, value):
        return [value]

    def gather(self, value, root=0):
        return [value]

    def bcast(self, value, root=0):
        return value

    def Allreduce(self, source, target):
        target[...] = source

    def alltoall(self, value):
        return value

    def Barrier(self):
        return None


class _Field:
    def __init__(self, values):
        self.name = "state"
        self.global_size = len(values)
        self.owned_start = 0
        self.owned_stop = len(values)
        self.values = np.asarray(values, dtype=float)

    def read_owned(self):
        return self.values

    def write_owned(self, values):
        self.values[:] = values

    def synchronize_ghosts(self):
        return None


def _member(member_id):
    field = _Field(np.arange(5.0) + 10.0 * member_id)
    registry = DistributedFieldRegistry([field], layout_id="native-cycle-v1")
    return NativeDistributedMember(member_id, registry, {"field": field})


class _Adapter:
    def forecast_native_member(self, member, timestep, **kwargs):
        member.unpack_owned(member.pack_owned() + float(timestep))

    def observe_native_member(self, member, observation_rows, **kwargs):
        return member.pack_owned()[observation_rows]

    def finalize_native_analysis(self, member, forecast_owned, timestep, **kwargs):
        member.model_context["finalized"] = int(timestep)

    def inverse_native_member(self, member, timestep, **kwargs):
        member.unpack_owned(member.pack_owned() + 3.0)
        member.model_context["inverted"] = int(timestep)


def _topology():
    comm = _Comm()
    return SimpleNamespace(
        world=comm,
        spatial_comm=comm,
        ensemble_comm=comm,
        world_rank=0,
        world_size=1,
        ensemble_slot=0,
        spatial_rank=0,
        ensemble_groups=1,
        spatial_ranks=1,
        is_ensemble_root=True,
    )


def test_native_cycle_matches_dense_stochastic_transform_and_checkpoints(tmp_path):
    pool = NativeDistributedMemberPool({0: _member(0), 1: _member(1)})
    adapter = _Adapter()
    topology = _topology()
    observed = np.array([2.0, 7.0])
    batch = NativeObservationBatch(np.array([0, 4]), observed)

    result = run_native_global_analysis_cycle(
        pool,
        adapter,
        1,
        [batch],
        number_of_batches=1,
        topology=topology,
        icesee_kwargs={"Nens": 2},
        error_mode="legacy_prior_anomalies",
        state_row_chunk_size=2,
        checkpoint_request=NativeCheckpointRequest(
            tmp_path, "native-cycle-test", {"random_stream": "member-keyed"}
        ),
    )

    dense_forecast = np.column_stack([
        np.arange(5.0) + 1.0,
        np.arange(5.0) + 11.0,
    ])
    products = StochasticAnalysisProducts.zeros(2)
    products.add_chunk(
        dense_forecast[[0, 4], :], observed,
        error_mode="legacy_prior_anomalies",
    )
    expected_transform = ensemble_transform_from_products(products)
    expected_analysis = dense_forecast @ expected_transform

    np.testing.assert_allclose(result.transform, expected_transform)
    np.testing.assert_allclose(result.local_analysis.members[0], expected_analysis[:, 0])
    np.testing.assert_allclose(result.local_analysis.members[1], expected_analysis[:, 1])
    np.testing.assert_allclose(pool.member(0).pack_owned(), expected_analysis[:, 0])
    assert result.observation_rows == 2
    assert pool.member(0).model_context["finalized"] == 1
    assert result.checkpoint is not None
    assert (result.checkpoint.path / "manifest.json").is_file()


def test_native_observation_batch_rejects_misaligned_values():
    try:
        NativeObservationBatch(np.array([1, 2]), np.array([1.0]))
    except ValueError as error:
        assert "align" in str(error)
    else:  # pragma: no cover
        raise AssertionError("misaligned observation batch was accepted")


def test_native_cycle_applies_native_inversion_without_gathering():
    pool = NativeDistributedMemberPool({0: _member(0), 1: _member(1)})
    result = run_native_global_analysis_cycle(
        pool,
        _Adapter(),
        1,
        [NativeObservationBatch(np.array([0]), np.array([2.0]))],
        number_of_batches=1,
        topology=_topology(),
        icesee_kwargs={"Nens": 2},
        error_mode="legacy_prior_anomalies",
        apply_inversion=True,
    )
    np.testing.assert_allclose(
        pool.member(0).pack_owned(), result.local_analysis.members[0]
    )
    assert pool.member(0).model_context["inverted"] == 1


def test_native_grouped_local_cycle_updates_only_requested_rows(tmp_path):
    pool = NativeDistributedMemberPool({0: _member(0), 1: _member(1)})
    target = DistributedAnalysisTargets(
        "state",
        np.array([1, 3]),
        np.array([1, 3]),
        np.array([[0.25], [0.75]]),
        ("state",),
    )
    observations = DistributedObservationCoordinates(
        "state",
        np.array([0, 4]),
        np.array([[0.2], [0.8]]),
    )
    result = run_native_grouped_local_analysis_cycle(
        pool,
        _Adapter(),
        1,
        [NativeLocalAnalysisRequest(target, observations, np.array([2.0, 7.0]), 0.1)],
        topology=_topology(),
        icesee_kwargs={"Nens": 2},
        error_mode="legacy_prior_anomalies",
        target_chunk_size=1,
        checkpoint_request=NativeCheckpointRequest(tmp_path, "local-native-test"),
    )

    forecast = np.column_stack([np.arange(5.0) + 1.0, np.arange(5.0) + 11.0])
    analysis = np.column_stack([
        result.local_analysis.members[0],
        result.local_analysis.members[1],
    ])
    np.testing.assert_allclose(analysis[[0, 2, 4]], forecast[[0, 2, 4]])
    assert np.max(np.abs(analysis[[1, 3]] - forecast[[1, 3]])) > 0.0
    assert result.target_blocks == 1
    assert result.observation_row_uses == 2
    assert result.checkpoint is not None


def test_native_local_request_combines_multiple_observation_kinds_once():
    target = DistributedAnalysisTargets(
        "geometry",
        np.array([1, 3]),
        np.array([1, 3]),
        np.array([[0.25], [0.75]]),
        ("thickness", "speed"),
    )
    observations = DistributedObservationCollection((
        DistributedObservationCoordinates(
            "thickness", np.array([10]), np.array([[0.2]])
        ),
        DistributedObservationCoordinates(
            "speed", np.array([20]), np.array([[0.8]])
        ),
    ))
    request = NativeLocalAnalysisRequest(
        target, observations, np.array([2.0, 7.0]), 0.1
    )
    assert request.observations.kinds == ("thickness", "speed")
    np.testing.assert_array_equal(request.observations.observation_ids, [10, 20])


def test_native_checkpoint_restores_persistent_fields(tmp_path):
    pool = NativeDistributedMemberPool({0: _member(0), 1: _member(1)})
    result = run_native_global_analysis_cycle(
        pool,
        _Adapter(),
        1,
        [NativeObservationBatch(np.array([0]), np.array([2.0]))],
        number_of_batches=1,
        topology=_topology(),
        icesee_kwargs={"Nens": 2},
        error_mode="legacy_prior_anomalies",
        checkpoint_request=NativeCheckpointRequest(tmp_path, "restore-test"),
    )
    expected = pool.member(0).pack_owned().copy()
    pool.member(0).unpack_owned(np.full(5, -99.0))
    checkpoint = restore_native_pool_from_checkpoint(
        pool,
        result.checkpoint.path,
        _topology(),
        number_of_members=2,
        expected_run_id="restore-test",
    )
    np.testing.assert_allclose(pool.member(0).pack_owned(), expected)
    assert checkpoint.timestep == 1


class _MismatchedComm:
    def allgather(self, value):
        return [value, b"different"]


def test_native_batch_alignment_rejects_different_ensemble_slot_order():
    batch = NativeObservationBatch(np.array([4, 2]), np.array([1.0, 3.0]))
    try:
        _validate_batch_alignment(batch, _MismatchedComm())
    except ValueError as error:
        assert "disagree" in str(error)
    else:  # pragma: no cover
        raise AssertionError("inconsistent observation ordering was accepted")


class _FailingInversionAdapter(_Adapter):
    def inverse_native_member(self, member, timestep, **kwargs):
        member.unpack_owned(member.pack_owned() + 1000.0)
        raise ValueError("synthetic inversion failure")


def test_native_inversion_failure_rolls_back_owned_fields_collectively():
    pool = NativeDistributedMemberPool({0: _member(0), 1: _member(1)})
    before = pool.snapshot_owned()
    try:
        apply_native_memberwise_inversion(
            pool,
            _FailingInversionAdapter(),
            8,
            topology=_topology(),
            icesee_kwargs={"Nens": 2},
        )
    except RuntimeError as error:
        assert "pre-inversion state restored" in str(error)
    else:  # pragma: no cover
        raise AssertionError("failed native inversion was accepted")
    for member_id in pool.member_ids:
        np.testing.assert_array_equal(
            pool.member(member_id).pack_owned(), before.members[member_id]
        )
