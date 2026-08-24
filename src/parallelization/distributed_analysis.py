"""Bounded-memory stochastic-analysis primitives for execution mode 3.

The routines in this module preserve ICESEE's existing stochastic EnKF
analysis.  They change data ownership and communication only; they do not
introduce an LETKF, EnSRF, or another filter.  A spatial rank may assemble a
small row block across all ensemble members, but no operation reconstructs a
complete model member.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterator, Mapping

import numpy as np


@dataclass(frozen=True)
class DistributedObservationLayout:
    """Contiguous ownership of active observation rows.

    Observation rows use their canonical active-row ordering, independently
    of MPI rank placement.  The stable global IDs are retained so model
    adapters can map a row back to a state variable, mesh node, or nonlinear
    local observation operator.
    """

    global_count: int
    owned_start: int
    owned_stop: int
    global_row_ids: np.ndarray
    layout_id: str = "observations-v1"

    def __post_init__(self) -> None:
        count = int(self.global_count)
        start = int(self.owned_start)
        stop = int(self.owned_stop)
        rows = np.asarray(self.global_row_ids, dtype=np.int64)
        if count < 0 or not 0 <= start <= stop <= count:
            raise ValueError("invalid distributed observation interval")
        if rows.ndim != 1 or rows.size != stop - start:
            raise ValueError("global_row_ids must match the owned interval")
        if rows.size and np.unique(rows).size != rows.size:
            raise ValueError("global observation row IDs must be unique locally")
        if not str(self.layout_id):
            raise ValueError("layout_id must be nonempty")
        object.__setattr__(self, "global_count", count)
        object.__setattr__(self, "owned_start", start)
        object.__setattr__(self, "owned_stop", stop)
        object.__setattr__(self, "global_row_ids", rows)
        object.__setattr__(self, "layout_id", str(self.layout_id))

    @property
    def owned_count(self) -> int:
        return self.owned_stop - self.owned_start

    def checkpoint_metadata(self) -> dict[str, int | str]:
        return {
            "layout_id": self.layout_id,
            "global_count": self.global_count,
            "owned_start": self.owned_start,
            "owned_stop": self.owned_stop,
        }


def contiguous_observation_layout(
    global_row_ids: np.ndarray,
    topology: Any,
    *,
    layout_id: str = "observations-v1",
) -> DistributedObservationLayout:
    """Partition canonical active observations over spatial ranks.

    The same layout is present in every ensemble group.  Analysis products
    must therefore be accumulated on one ensemble slot (normally slot zero)
    and reduced over that slot's ``spatial_comm`` to avoid double counting.
    """

    rows = np.asarray(global_row_ids, dtype=np.int64).ravel()
    rank = int(topology.spatial_rank)
    size = int(topology.spatial_ranks)
    if size <= 0 or not 0 <= rank < size:
        raise ValueError("topology contains invalid spatial coordinates")
    quotient, remainder = divmod(rows.size, size)
    start = rank * quotient + min(rank, remainder)
    stop = start + quotient + (rank < remainder)
    return DistributedObservationLayout(
        global_count=rows.size,
        owned_start=start,
        owned_stop=stop,
        global_row_ids=rows[start:stop],
        layout_id=layout_id,
    )


def validate_observation_partition(
    layout: DistributedObservationLayout,
    spatial_comm: Any,
) -> tuple[tuple[int, int], ...]:
    """Validate gap-free active-row ownership without gathering row arrays."""

    descriptor = (
        layout.owned_start,
        layout.owned_stop,
        layout.global_count,
        layout.layout_id,
    )
    descriptors = spatial_comm.allgather(descriptor)
    if not descriptors:
        raise ValueError("spatial communicator returned no observation layouts")
    if {int(item[2]) for item in descriptors} != {layout.global_count}:
        raise ValueError("spatial ranks disagree on global observation count")
    if {str(item[3]) for item in descriptors} != {layout.layout_id}:
        raise ValueError("spatial ranks disagree on observation layout_id")
    intervals = sorted((int(item[0]), int(item[1])) for item in descriptors)
    cursor = 0
    for start, stop in intervals:
        if start < cursor:
            raise ValueError("distributed observation ownership overlaps")
        if start > cursor:
            raise ValueError("distributed observation ownership contains a gap")
        if stop < start or stop > layout.global_count:
            raise ValueError("a rank reported an invalid observation interval")
        cursor = stop
    if cursor != layout.global_count:
        raise ValueError("observation ownership does not cover active rows")
    return tuple(intervals)


def iter_local_ensemble_row_blocks(
    local_members: Mapping[int, np.ndarray],
    number_of_members: int,
    ensemble_comm: Any,
    *,
    row_chunk_size: int,
) -> Iterator[tuple[slice, np.ndarray]]:
    """Assemble bounded local-row blocks across all ensemble members.

    Each yielded matrix has shape ``(local_rows_in_chunk, Nens)``.  The
    collective runs only across ranks holding the same spatial slab, so this
    is not a complete-member gather.  Member IDs, rather than rank order,
    determine columns and make the result invariant to ensemble scheduling.
    """

    number_of_members = int(number_of_members)
    row_chunk_size = int(row_chunk_size)
    if number_of_members <= 0:
        raise ValueError("number_of_members must be positive")
    if row_chunk_size <= 0:
        raise ValueError("row_chunk_size must be positive")
    normalized: dict[int, np.ndarray] = {}
    local_size = None
    dtype = None
    for member_id, values in local_members.items():
        member_id = int(member_id)
        array = np.asarray(values)
        if array.ndim != 1:
            raise ValueError("local member slabs must be one-dimensional")
        if local_size is None:
            local_size = array.size
            dtype = array.dtype
        elif array.size != local_size:
            raise ValueError("local member slabs must have equal lengths")
        if member_id in normalized:
            raise ValueError("duplicate local ensemble member ID")
        normalized[member_id] = array
    if local_size is None:
        # Ranks with no scheduled members still need the slab size.  Exchange
        # scalar sizes only, then participate in every block collective.
        sizes = ensemble_comm.allgather(None)
        known = [int(value) for value in sizes if value is not None]
        if not known or len(set(known)) != 1:
            raise ValueError("cannot determine a common local slab size")
        local_size = known[0]
        dtype = np.dtype(np.float64)
    else:
        sizes = ensemble_comm.allgather(int(local_size))
        known = [int(value) for value in sizes if value is not None]
        if len(set(known)) != 1:
            raise ValueError("ensemble ranks disagree on local slab size")

    for start in range(0, int(local_size), row_chunk_size):
        stop = min(int(local_size), start + row_chunk_size)
        payload = {
            member_id: np.ascontiguousarray(values[start:stop])
            for member_id, values in normalized.items()
        }
        gathered = ensemble_comm.allgather(payload)
        columns: dict[int, np.ndarray] = {}
        for rank_payload in gathered:
            for member_id, values in rank_payload.items():
                member_id = int(member_id)
                if member_id in columns:
                    raise ValueError("ensemble member is owned by multiple slots")
                columns[member_id] = np.asarray(values)
        expected = set(range(number_of_members))
        if set(columns) != expected:
            missing = sorted(expected - set(columns))
            extra = sorted(set(columns) - expected)
            raise ValueError(
                f"ensemble member ownership is incomplete; missing={missing}, "
                f"extra={extra}"
            )
        block_dtype = np.result_type(dtype, *[value.dtype for value in columns.values()])
        block = np.empty((stop - start, number_of_members), dtype=block_dtype)
        for member_id in range(number_of_members):
            values = columns[member_id]
            if values.shape != (stop - start,):
                raise ValueError("gathered member chunk has an invalid shape")
            block[:, member_id] = values
        yield slice(start, stop), block


def assemble_selected_ensemble_rows(
    local_members: Mapping[int, np.ndarray],
    selected_rows: np.ndarray,
    number_of_members: int,
    ensemble_comm: Any,
) -> np.ndarray:
    """Assemble selected local rows across members, never a complete member.

    This is the grouped-local counterpart of
    :func:`iter_local_ensemble_row_blocks`.  Callers bound memory by passing a
    single target or observation chunk.  Stable member IDs determine columns,
    so scheduling and MPI rank placement do not affect the result.
    """

    rows = np.asarray(selected_rows, dtype=np.int64).ravel()
    number_of_members = int(number_of_members)
    if number_of_members <= 0:
        raise ValueError("number_of_members must be positive")
    payload: dict[int, np.ndarray] = {}
    for member_id, values in local_members.items():
        member_id = int(member_id)
        array = np.asarray(values)
        if array.ndim != 1:
            raise ValueError("local member slabs must be one-dimensional")
        if rows.size and (rows.min() < 0 or rows.max() >= array.size):
            raise ValueError("selected rows leave the local member slab")
        payload[member_id] = np.ascontiguousarray(array[rows])
    gathered = ensemble_comm.allgather(payload)
    columns: dict[int, np.ndarray] = {}
    for rank_payload in gathered:
        for member_id, values in rank_payload.items():
            member_id = int(member_id)
            if member_id in columns:
                raise ValueError("ensemble member is owned by multiple slots")
            columns[member_id] = np.asarray(values)
    expected = set(range(number_of_members))
    if set(columns) != expected:
        missing = sorted(expected - set(columns))
        extra = sorted(set(columns) - expected)
        raise ValueError(
            f"ensemble member ownership is incomplete; missing={missing}, extra={extra}"
        )
    dtype = np.result_type(*[value.dtype for value in columns.values()])
    block = np.empty((rows.size, number_of_members), dtype=dtype)
    for member_id in range(number_of_members):
        if columns[member_id].shape != (rows.size,):
            raise ValueError("gathered selected member rows have an invalid shape")
        block[:, member_id] = columns[member_id]
    return block


@dataclass
class StochasticAnalysisProducts:
    """Fixed-size ensemble-space products accumulated from row chunks."""

    cross: np.ndarray
    gram: np.ndarray
    rhs: np.ndarray
    observation_rows: int = 0

    @classmethod
    def zeros(cls, number_of_members: int) -> "StochasticAnalysisProducts":
        shape = (int(number_of_members), int(number_of_members))
        if shape[0] <= 0:
            raise ValueError("number_of_members must be positive")
        return cls(*(np.zeros(shape, dtype=np.float64) for _ in range(3)))

    def add_chunk(
        self,
        forecast_observations: np.ndarray,
        observations: np.ndarray,
        *,
        error_mode: str,
        observation_errors: np.ndarray | None = None,
    ) -> None:
        """Accumulate one observation-row chunk using current ICESEE modes."""

        forecast = np.asarray(forecast_observations, dtype=np.float64)
        observed = np.asarray(observations, dtype=np.float64).ravel()
        if forecast.ndim != 2 or forecast.shape[0] != observed.size:
            raise ValueError("forecast and observation chunk shapes disagree")
        if forecast.shape[1] != self.cross.shape[0]:
            raise ValueError("forecast ensemble size disagrees with products")
        yprime = forecast - forecast.mean(axis=1, keepdims=True)
        mode = str(error_mode).lower()
        if mode == "legacy_prior_anomalies":
            eta = yprime
            innovations = observed[:, None] - forecast
        elif mode in {"stochastic_r", "generated_r"}:
            if observation_errors is None:
                raise ValueError(f"{error_mode} requires observation_errors")
            eta = np.asarray(observation_errors, dtype=np.float64)
            if eta.shape != forecast.shape:
                raise ValueError("observation_errors must match forecast shape")
            innovations = observed[:, None] + eta - forecast
        else:
            raise ValueError(
                "error_mode must be stochastic_R, generated_R, or "
                "legacy_prior_anomalies"
            )
        factor = yprime + eta
        self.cross += yprime.T @ factor
        self.gram += factor.T @ factor
        self.rhs += factor.T @ innovations
        self.observation_rows += observed.size

    def reduced(self, communicator: Any) -> "StochasticAnalysisProducts":
        """Sum products over observation owners with fixed memory use."""

        reduced_arrays = []
        for local in (self.cross, self.gram, self.rhs):
            target = np.empty_like(local)
            communicator.Allreduce(local, target)
            reduced_arrays.append(target)
        row_counts = communicator.allgather(int(self.observation_rows))
        return StochasticAnalysisProducts(
            *reduced_arrays,
            observation_rows=sum(int(value) for value in row_counts),
        )


def ensemble_transform_from_products(
    products: StochasticAnalysisProducts,
    *,
    energy_fraction: float = 0.999,
) -> np.ndarray:
    """Compute ICESEE's existing stochastic transform from streamed products."""

    cross = np.asarray(products.cross, dtype=np.float64)
    gram = np.asarray(products.gram, dtype=np.float64)
    rhs = np.asarray(products.rhs, dtype=np.float64)
    if cross.shape != gram.shape or gram.shape != rhs.shape:
        raise ValueError("analysis products must have identical square shapes")
    gram = 0.5 * (gram + gram.T)
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = np.maximum(values[order], 0.0)
    vectors = vectors[:, order]
    inverse = np.zeros_like(values)
    total = float(values.sum())
    if total > 0.0:
        cumulative = np.cumsum(values)
        nkeep = int(np.searchsorted(
            cumulative, float(energy_fraction) * total, side="left"
        ) + 1)
        tolerance = np.finfo(values.dtype).eps * max(gram.shape) * values[0]
        keep = (np.arange(values.size) < nkeep) & (values > tolerance)
        inverse[keep] = 1.0 / values[keep]
    gram_pinv = (vectors * inverse) @ vectors.T
    return np.eye(gram.shape[0]) + cross @ gram_pinv @ gram_pinv @ rhs
