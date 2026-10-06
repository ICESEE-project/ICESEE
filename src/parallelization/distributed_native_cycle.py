"""Model-agnostic bounded-memory cycle for execution-mode-3 development.

This module composes already parity-tested mode-3 primitives while runtime
selection remains disabled.  Forecasts and model consistency operate on
persistent model-native distributed fields.  Observations are consumed in
spatially local batches, analysis products remain ``Nens x Nens``, and state
transforms operate on bounded owned-row chunks.  No step reconstructs a full
member or the full observation table.

The stochastic analysis is ICESEE's existing analysis.  This module changes
ownership and communication, not the filter or its scientific semantics.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from os import PathLike
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from .distributed_analysis import (
    StochasticAnalysisProducts,
    assemble_selected_ensemble_rows,
    ensemble_transform_from_products,
)
from .distributed_checkpoint import (
    DistributedCheckpoint,
    load_distributed_checkpoint,
    save_distributed_checkpoint,
)
from .distributed_adapter import (
    DistributedAnalysisTargets,
    DistributedObservationCollection,
    DistributedObservationCoordinates,
)
from .distributed_local_runtime import apply_grouped_local_analysis_local
from .distributed_native_runtime import (
    NativeDistributedMemberPool,
    NativeDistributedInversionAdapter,
    NativeDistributedModelAdapter,
    apply_native_memberwise_inversion,
)
from .distributed_runtime import LocalMemberEnsemble, transform_local_members


@dataclass(frozen=True)
class NativeObservationBatch:
    """One spatial rank's bounded portion of an observation cycle.

    ``observation_ids`` are interpreted by the model adapter and need not be
    packed state rows.  ``values`` are the matching measured values.  For the
    stochastic-R/generated-R modes, ``member_errors`` contains only errors for
    members scheduled on this ensemble slot.  The ensemble communicator joins
    those member-keyed columns without gathering model state.
    """

    observation_ids: np.ndarray
    values: np.ndarray
    member_errors: Mapping[int, np.ndarray] | None = None

    def __post_init__(self) -> None:
        raw_ids = np.asarray(self.observation_ids)
        values = np.asarray(self.values, dtype=np.float64)
        if raw_ids.ndim != 1 or not np.issubdtype(raw_ids.dtype, np.integer):
            raise TypeError("observation_ids must be a one-dimensional integer array")
        if values.ndim != 1 or values.size != raw_ids.size:
            raise ValueError("observation values must align with observation_ids")
        errors = None
        if self.member_errors is not None:
            errors = {}
            for member_id, member_values in self.member_errors.items():
                array = np.asarray(member_values, dtype=np.float64)
                if array.ndim != 1 or array.size != raw_ids.size:
                    raise ValueError(
                        "member observation errors must align with observation_ids"
                    )
                errors[int(member_id)] = np.ascontiguousarray(array)
        object.__setattr__(
            self, "observation_ids", np.ascontiguousarray(raw_ids, dtype=np.int64)
        )
        object.__setattr__(self, "values", np.ascontiguousarray(values))
        object.__setattr__(self, "member_errors", errors)


@dataclass(frozen=True)
class NativeCheckpointRequest:
    """Optional rank-sharded checkpoint emitted after a completed cycle."""

    root: str | PathLike[str]
    run_id: str
    metadata: Mapping[str, Any] | None = None


@dataclass(frozen=True)
class NativeCycleResult:
    """Bounded local products returned by one native forecast-analysis cycle."""

    local_forecast: LocalMemberEnsemble
    local_analysis: LocalMemberEnsemble
    transform: np.ndarray
    observation_rows: int
    checkpoint: DistributedCheckpoint | None = None


@dataclass(frozen=True)
class NativeLocalAnalysisRequest:
    """One target block and its joint locally owned observation set.

    ``observations`` may describe one kind or a collection of kinds. All rows
    participate in one stochastic transform, preserving mode-1 semantics.
    """

    target: DistributedAnalysisTargets
    observations: (
        DistributedObservationCoordinates | DistributedObservationCollection
    )
    values: np.ndarray
    radius: float
    member_errors: Mapping[int, np.ndarray] | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.target, DistributedAnalysisTargets):
            raise TypeError("target must be DistributedAnalysisTargets")
        if not isinstance(
            self.observations,
            (DistributedObservationCoordinates, DistributedObservationCollection),
        ):
            raise TypeError(
                "observations must be DistributedObservationCoordinates or "
                "DistributedObservationCollection"
            )
        values = np.asarray(self.values, dtype=np.float64).ravel()
        if values.size != self.observations.observation_ids.size:
            raise ValueError("local observation values and metadata must align")
        kinds = (
            (self.observations.kind,)
            if isinstance(self.observations, DistributedObservationCoordinates)
            else self.observations.kinds
        )
        rejected = sorted(set(kinds) - set(self.target.observation_kinds))
        if rejected:
            raise ValueError(
                "analysis target does not accept observation kinds: "
                + ", ".join(rejected)
            )
        radius = float(self.radius)
        if not np.isfinite(radius) or radius <= 0.0:
            raise ValueError("localization radius must be finite and positive")
        errors = None
        if self.member_errors is not None:
            errors = {}
            for member_id, member_values in self.member_errors.items():
                array = np.asarray(member_values, dtype=np.float64).ravel()
                if array.size != values.size:
                    raise ValueError(
                        "member observation errors must align with observations"
                    )
                errors[int(member_id)] = np.ascontiguousarray(array)
        object.__setattr__(self, "values", np.ascontiguousarray(values))
        object.__setattr__(self, "radius", radius)
        object.__setattr__(self, "member_errors", errors)


@dataclass(frozen=True)
class NativeLocalCycleResult:
    """Local products from one persistent grouped-local stochastic cycle."""

    local_forecast: LocalMemberEnsemble
    local_analysis: LocalMemberEnsemble
    target_blocks: int
    observation_row_uses: int
    checkpoint: DistributedCheckpoint | None = None


def _batch_fingerprint(batch: NativeObservationBatch) -> bytes:
    """Return a stable identity for batch ordering and measured values."""

    digest = hashlib.blake2b(digest_size=20)
    digest.update(np.asarray(batch.observation_ids, dtype="<i8").tobytes())
    digest.update(np.asarray(batch.values, dtype="<f8").tobytes())
    return digest.digest()


def _validate_batch_alignment(batch: NativeObservationBatch, ensemble_comm: Any) -> None:
    """Reject inconsistent observation order before analysis collectives."""

    fingerprint = _batch_fingerprint(batch)
    fingerprints = ensemble_comm.allgather(fingerprint)
    if any(value != fingerprint for value in fingerprints):
        raise ValueError(
            "ensemble slots disagree on observation IDs, order, or values"
        )


def _local_request_fingerprint(request: NativeLocalAnalysisRequest) -> bytes:
    digest = hashlib.blake2b(digest_size=20)
    digest.update(request.target.name.encode("utf-8"))
    for kind in request.target.observation_kinds:
        digest.update(kind.encode("utf-8"))
        digest.update(b"\0")
    digest.update(np.asarray(request.target.local_rows, dtype="<i8").tobytes())
    digest.update(np.asarray(request.target.global_rows, dtype="<i8").tobytes())
    digest.update(np.asarray(request.target.coordinates, dtype="<f8").tobytes())
    kinds = (
        (request.observations.kind,)
        if isinstance(request.observations, DistributedObservationCoordinates)
        else request.observations.kinds
    )
    for kind in kinds:
        digest.update(kind.encode("utf-8"))
        digest.update(b"\0")
    digest.update(np.asarray(request.observations.observation_ids, dtype="<i8").tobytes())
    digest.update(np.asarray(request.observations.coordinates, dtype="<f8").tobytes())
    digest.update(np.asarray(request.values, dtype="<f8").tobytes())
    digest.update(np.asarray([request.radius], dtype="<f8").tobytes())
    return digest.digest()


def _validate_local_request_alignment(
    request: NativeLocalAnalysisRequest,
    ensemble_comm: Any,
) -> None:
    fingerprint = _local_request_fingerprint(request)
    fingerprints = ensemble_comm.allgather(fingerprint)
    if any(value != fingerprint for value in fingerprints):
        raise ValueError(
            "ensemble slots disagree on local-analysis metadata, order, or values"
        )


def _assemble_observation_errors(
    batch: NativeObservationBatch,
    *,
    number_of_members: int,
    ensemble_comm: Any,
) -> np.ndarray | None:
    if batch.member_errors is None:
        return None
    rows = np.arange(batch.observation_ids.size, dtype=np.int64)
    return assemble_selected_ensemble_rows(
        batch.member_errors,
        rows,
        number_of_members,
        ensemble_comm,
    ).astype(np.float64, copy=False)


def run_native_global_analysis_cycle(
    pool: NativeDistributedMemberPool,
    adapter: NativeDistributedModelAdapter,
    timestep: int,
    observation_batches: Iterable[NativeObservationBatch],
    *,
    number_of_batches: int,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
    error_mode: str,
    state_row_chunk_size: int = 4096,
    energy_fraction: float = 0.999,
    checkpoint_request: NativeCheckpointRequest | None = None,
    apply_inversion: bool = False,
) -> NativeCycleResult:
    """Advance, analyze, finalize, and optionally checkpoint native members.

    Every ensemble slot at one spatial coordinate must expose the same batch
    count and observation IDs.  Batch *contents* are spatially local and may
    differ across spatial ranks.  Only ensemble slot zero accumulates and
    spatially reduces observation products, avoiding duplicate observations;
    the resulting small transform is broadcast across ensemble slots.
    """

    timestep = int(timestep)
    if timestep < 0:
        raise ValueError("timestep must be nonnegative")
    number_of_batches = int(number_of_batches)
    if number_of_batches < 0:
        raise ValueError("number_of_batches must be nonnegative")
    counts = topology.ensemble_comm.allgather(number_of_batches)
    if any(int(value) != number_of_batches for value in counts):
        raise ValueError("ensemble slots disagree on local observation batch count")

    number_of_members = int(icesee_kwargs.get("Nens", 1))
    local_forecast = pool.forecast_with_adapter(
        timestep,
        adapter,
        topology=topology,
        icesee_kwargs=icesee_kwargs,
    )
    products = StochasticAnalysisProducts.zeros(number_of_members)
    batches = iter(observation_batches)
    local_observation_rows = 0
    for batch_index in range(number_of_batches):
        try:
            batch = next(batches)
        except StopIteration as error:
            raise ValueError(
                "observation source ended before number_of_batches"
            ) from error
        if not isinstance(batch, NativeObservationBatch):
            raise TypeError("observation source must yield NativeObservationBatch")
        _validate_batch_alignment(batch, topology.ensemble_comm)
        forecast_values = pool.observe_with_adapter(
            batch.observation_ids,
            adapter,
            topology=topology,
            icesee_kwargs=icesee_kwargs,
        )
        rows = np.arange(batch.observation_ids.size, dtype=np.int64)
        forecast_matrix = assemble_selected_ensemble_rows(
            forecast_values,
            rows,
            number_of_members,
            topology.ensemble_comm,
        ).astype(np.float64, copy=False)
        error_matrix = _assemble_observation_errors(
            batch,
            number_of_members=number_of_members,
            ensemble_comm=topology.ensemble_comm,
        )
        if topology.is_ensemble_root:
            products.add_chunk(
                forecast_matrix,
                batch.values,
                error_mode=error_mode,
                observation_errors=error_matrix,
            )
            local_observation_rows += int(batch.values.size)

    try:
        next(batches)
    except StopIteration:
        pass
    else:
        raise ValueError("observation source exceeded number_of_batches")

    if topology.is_ensemble_root:
        global_products = products.reduced(topology.spatial_comm)
        transform = ensemble_transform_from_products(
            global_products,
            energy_fraction=energy_fraction,
        )
        observation_rows = int(global_products.observation_rows)
    else:
        transform = None
        observation_rows = None
    transform = np.asarray(
        topology.ensemble_comm.bcast(transform, root=0), dtype=np.float64
    )
    observation_rows = int(
        topology.ensemble_comm.bcast(observation_rows, root=0)
    )

    local_analysis = transform_local_members(
        local_forecast,
        transform,
        topology,
        icesee_kwargs,
        row_chunk_size=state_row_chunk_size,
    )
    pool.restore_and_finalize(
        local_analysis,
        local_forecast,
        timestep,
        adapter,
        topology=topology,
        icesee_kwargs=icesee_kwargs,
    )
    if apply_inversion:
        local_analysis = apply_native_memberwise_inversion(
            pool,
            adapter,
            timestep,
            topology=topology,
            icesee_kwargs=icesee_kwargs,
        )

    checkpoint = None
    if checkpoint_request is not None:
        metadata = dict(checkpoint_request.metadata or {})
        metadata.update({
            "cycle_timestep": timestep,
            "observation_rows": observation_rows,
            "error_mode": str(error_mode),
        })
        checkpoint = save_distributed_checkpoint(
            checkpoint_request.root,
            timestep,
            pool.snapshot_owned(),
            topology,
            run_id=checkpoint_request.run_id,
            metadata=metadata,
        )

    return NativeCycleResult(
        local_forecast=local_forecast,
        local_analysis=local_analysis,
        transform=transform,
        observation_rows=observation_rows,
        checkpoint=checkpoint,
    )


def run_native_grouped_local_analysis_cycle(
    pool: NativeDistributedMemberPool,
    adapter: NativeDistributedModelAdapter,
    timestep: int,
    requests: Sequence[NativeLocalAnalysisRequest],
    *,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
    error_mode: str,
    target_chunk_size: int = 1024,
    energy_fraction: float = 0.999,
    checkpoint_request: NativeCheckpointRequest | None = None,
    apply_inversion: bool = False,
) -> NativeLocalCycleResult:
    """Run grouped-local stochastic analysis on persistent native fields.

    Requests contain only state targets and observations owned by the current
    spatial rank.  Sparse localization communication is delegated to the
    parity-tested grouped-local runtime.  Analysis blocks are applied
    sequentially to one compact owned snapshot and native geometry is
    finalized once after every requested block has been restored.
    """

    timestep = int(timestep)
    if timestep < 0:
        raise ValueError("timestep must be nonnegative")
    requests = tuple(requests)
    counts = topology.ensemble_comm.allgather(len(requests))
    if any(int(value) != len(requests) for value in counts):
        raise ValueError("ensemble slots disagree on local-analysis request count")

    local_forecast = pool.forecast_with_adapter(
        timestep,
        adapter,
        topology=topology,
        icesee_kwargs=icesee_kwargs,
    )
    local_analysis = local_forecast
    local_row_uses = 0
    for request in requests:
        if not isinstance(request, NativeLocalAnalysisRequest):
            raise TypeError("requests must contain NativeLocalAnalysisRequest")
        _validate_local_request_alignment(request, topology.ensemble_comm)
        forecast_values = pool.observe_with_adapter(
            request.observations.observation_ids,
            adapter,
            topology=topology,
            icesee_kwargs=icesee_kwargs,
        )
        local_analysis = apply_grouped_local_analysis_local(
            local_analysis,
            request.target,
            request.observations,
            forecast_values,
            request.values,
            radius=request.radius,
            topology=topology,
            icesee_kwargs=icesee_kwargs,
            error_mode=error_mode,
            local_observation_errors=request.member_errors,
            target_chunk_size=target_chunk_size,
            energy_fraction=energy_fraction,
        )
        local_row_uses += int(request.values.size)

    pool.restore_and_finalize(
        local_analysis,
        local_forecast,
        timestep,
        adapter,
        topology=topology,
        icesee_kwargs=icesee_kwargs,
    )
    if apply_inversion:
        local_analysis = apply_native_memberwise_inversion(
            pool,
            adapter,
            timestep,
            topology=topology,
            icesee_kwargs=icesee_kwargs,
        )

    if topology.is_ensemble_root:
        spatial_uses = sum(
            int(value) for value in topology.spatial_comm.allgather(local_row_uses)
        )
    else:
        spatial_uses = None
    observation_row_uses = int(
        topology.ensemble_comm.bcast(spatial_uses, root=0)
    )

    checkpoint = None
    if checkpoint_request is not None:
        metadata = dict(checkpoint_request.metadata or {})
        metadata.update({
            "cycle_timestep": timestep,
            "analysis_kind": "grouped_local_stochastic",
            "target_blocks": len(requests),
            "observation_row_uses": observation_row_uses,
            "error_mode": str(error_mode),
        })
        checkpoint = save_distributed_checkpoint(
            checkpoint_request.root,
            timestep,
            pool.snapshot_owned(),
            topology,
            run_id=checkpoint_request.run_id,
            metadata=metadata,
        )

    return NativeLocalCycleResult(
        local_forecast=local_forecast,
        local_analysis=local_analysis,
        target_blocks=len(requests),
        observation_row_uses=observation_row_uses,
        checkpoint=checkpoint,
    )


def restore_native_pool_from_checkpoint(
    pool: NativeDistributedMemberPool,
    checkpoint_path: str | PathLike[str],
    topology: Any,
    *,
    number_of_members: int,
    expected_run_id: str | None = None,
    adapter: NativeDistributedModelAdapter | None = None,
    icesee_kwargs: Mapping[str, Any] | None = None,
) -> DistributedCheckpoint:
    """Restore target-local checkpoint overlaps into persistent model fields.

    The loader supports a compatible process-grid change and reads only the
    member/block intervals owned under the current topology.  An adapter may
    optionally implement ``restore_native_checkpoint`` to recover auxiliary
    native solver state represented by identifiers in checkpoint metadata.
    """

    local_ensemble, checkpoint = load_distributed_checkpoint(
        checkpoint_path,
        pool.layout,
        topology,
        expected_run_id=expected_run_id,
        number_of_members=int(number_of_members),
    )
    pool.restore_owned(local_ensemble)
    restore = getattr(adapter, "restore_native_checkpoint", None)
    if callable(restore):
        kwargs = {} if icesee_kwargs is None else icesee_kwargs
        for member_id in pool.member_ids:
            result = restore(
                pool.member(member_id),
                checkpoint,
                topology=topology,
                icesee_kwargs=kwargs,
            )
            if result is not None:
                raise TypeError(
                    "native checkpoint restore callbacks mutate model context "
                    "and must return None"
                )
    return checkpoint
