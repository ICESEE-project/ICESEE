"""Distributed grouped-local stochastic EnKF analysis for execution mode 3.

This module distributes ICESEE's existing patch-based stochastic analysis.  It
does not implement an LETKF or a square-root filter.  Observation owners keep
local coordinate trees, target owners issue bounded spatial queries, and only
the observation terms needed by those targets cross the spatial communicator.
No rank gathers a global observation table or broadcasts a global patch map.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterator

import numpy as np
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class DistributedLocalPatch:
    """One transform shared by locally owned state rows."""

    local_rows: np.ndarray
    global_rows: np.ndarray
    observation_ids: np.ndarray
    transform: np.ndarray

    def __post_init__(self) -> None:
        local_rows = np.asarray(self.local_rows, dtype=np.int64)
        global_rows = np.asarray(self.global_rows, dtype=np.int64)
        observation_ids = np.asarray(self.observation_ids, dtype=np.int64)
        transform = np.asarray(self.transform, dtype=np.float64)
        if local_rows.ndim != 1 or global_rows.shape != local_rows.shape:
            raise ValueError("local_rows and global_rows must be matching vectors")
        if observation_ids.ndim != 1:
            raise ValueError("observation_ids must be one-dimensional")
        if transform.ndim != 2 or transform.shape[0] != transform.shape[1]:
            raise ValueError("transform must be square")
        object.__setattr__(self, "local_rows", local_rows)
        object.__setattr__(self, "global_rows", global_rows)
        object.__setattr__(self, "observation_ids", observation_ids)
        object.__setattr__(self, "transform", transform)


def _analysis_transform(yprime: np.ndarray, eta: np.ndarray,
                        innovations: np.ndarray, energy_fraction: float) -> np.ndarray:
    """Match the existing observation-space SVD used by grouped local analysis."""

    yprime = np.asarray(yprime, dtype=np.float64)
    eta = np.asarray(eta, dtype=np.float64)
    innovations = np.asarray(innovations, dtype=np.float64)
    if yprime.ndim != 2 or eta.shape != yprime.shape or innovations.shape != yprime.shape:
        raise ValueError("local observation terms must have matching 2-D shapes")
    number_of_members = yprime.shape[1]
    factor = yprime + eta
    left, singular_values, _ = np.linalg.svd(factor, full_matrices=False)
    eigenvalues = singular_values * singular_values
    inverse = np.zeros_like(eigenvalues)
    total = float(eigenvalues.sum())
    if total > 0.0:
        cumulative = np.cumsum(eigenvalues)
        nkeep = int(np.searchsorted(
            cumulative, float(energy_fraction) * total, side="left"
        ) + 1)
        tolerance = (
            np.finfo(eigenvalues.dtype).eps
            * max(factor.shape)
            * eigenvalues[0]
        )
        keep = (np.arange(eigenvalues.size) < nkeep) & (eigenvalues > tolerance)
        inverse[keep] = 1.0 / eigenvalues[keep]
    retained = left[:, :inverse.size]
    correction = yprime.T @ (
        retained @ (inverse[:, None] * (retained.T @ innovations))
    )
    return np.eye(number_of_members, dtype=np.float64) + correction


def _normalize_coordinates(values: np.ndarray, name: str) -> np.ndarray:
    coords = np.asarray(values, dtype=np.float64)
    if coords.ndim == 1:
        coords = coords[:, None]
    if coords.ndim != 2 or coords.shape[1] < 1:
        raise ValueError(f"{name} must have shape (n, spatial_dimensions)")
    if not np.all(np.isfinite(coords)):
        raise ValueError(f"{name} contains non-finite coordinates")
    return coords


def _bbox_intersects(points: np.ndarray, lower: np.ndarray, upper: np.ndarray,
                     radius: float) -> np.ndarray:
    below = np.maximum(lower[None, :] - points, 0.0)
    above = np.maximum(points - upper[None, :], 0.0)
    return np.sum((below + above) ** 2, axis=1) <= radius * radius


def iter_distributed_local_patches(
    *,
    target_global_rows: np.ndarray,
    target_coordinates: np.ndarray,
    observation_ids: np.ndarray,
    observation_coordinates: np.ndarray,
    yprime: np.ndarray,
    eta: np.ndarray,
    innovations: np.ndarray,
    radius: float,
    spatial_comm: Any,
    target_chunk_size: int = 1024,
    energy_fraction: float = 0.999,
) -> Iterator[tuple[slice, tuple[DistributedLocalPatch, ...]]]:
    """Yield exact grouped transforms for bounded chunks of owned target rows.

    Observation IDs must be stable canonical active-row IDs.  Sorting those IDs
    before each SVD makes the result independent of rank placement and message
    arrival order.  Every rank participates in the same number of query rounds,
    including ranks with no target rows or no observations.
    """

    targets = np.asarray(target_global_rows, dtype=np.int64).ravel()
    target_coords = _normalize_coordinates(target_coordinates, "target_coordinates")
    obs_ids = np.asarray(observation_ids, dtype=np.int64).ravel()
    obs_coords = _normalize_coordinates(
        observation_coordinates, "observation_coordinates"
    )
    yprime = np.asarray(yprime, dtype=np.float64)
    eta = np.asarray(eta, dtype=np.float64)
    innovations = np.asarray(innovations, dtype=np.float64)
    radius = float(radius)
    target_chunk_size = int(target_chunk_size)
    if targets.size != target_coords.shape[0]:
        raise ValueError("target rows and coordinates have different lengths")
    if obs_ids.size != obs_coords.shape[0]:
        raise ValueError("observation IDs and coordinates have different lengths")
    if yprime.ndim != 2 or yprime.shape[0] != obs_ids.size:
        raise ValueError("yprime rows must match locally owned observations")
    if eta.shape != yprime.shape or innovations.shape != yprime.shape:
        raise ValueError("eta and innovations must match yprime")
    if obs_ids.size and np.unique(obs_ids).size != obs_ids.size:
        raise ValueError("observation IDs must be unique on each owner")
    if radius <= 0.0 or not np.isfinite(radius):
        raise ValueError("radius must be finite and positive")
    if target_chunk_size <= 0:
        raise ValueError("target_chunk_size must be positive")
    if target_coords.shape[1] != obs_coords.shape[1]:
        raise ValueError("target and observation coordinate dimensions differ")

    dimension = target_coords.shape[1]
    if obs_ids.size:
        bounds = (obs_coords.min(axis=0), obs_coords.max(axis=0), obs_ids.size)
        tree = cKDTree(obs_coords)
    else:
        bounds = (np.full(dimension, np.inf), np.full(dimension, -np.inf), 0)
        tree = None
    all_bounds = spatial_comm.allgather(bounds)
    local_rounds = (targets.size + target_chunk_size - 1) // target_chunk_size
    rounds = max(int(value) for value in spatial_comm.allgather(local_rounds))
    communicator_size = int(spatial_comm.Get_size())

    for round_index in range(rounds):
        start = round_index * target_chunk_size
        stop = min(targets.size, start + target_chunk_size)
        chunk_coords = target_coords[start:stop]
        query_ids = np.arange(start, stop, dtype=np.int64)
        requests = []
        for lower, upper, count in all_bounds:
            if not count or not query_ids.size:
                requests.append((np.empty(0, dtype=np.int64),
                                 np.empty((0, dimension), dtype=np.float64)))
                continue
            selected = _bbox_intersects(
                chunk_coords, np.asarray(lower), np.asarray(upper), radius
            )
            requests.append((query_ids[selected], chunk_coords[selected]))
        if len(requests) != communicator_size:
            raise ValueError("spatial communicator size changed during analysis")

        received_requests = spatial_comm.alltoall(requests)
        replies = []
        for source_query_ids, source_coords in received_requests:
            source_query_ids = np.asarray(source_query_ids, dtype=np.int64)
            source_coords = np.asarray(source_coords, dtype=np.float64).reshape(-1, dimension)
            pair_queries = []
            pair_local_obs = []
            if tree is not None and source_query_ids.size:
                for query_id, neighbors in zip(
                    source_query_ids, tree.query_ball_point(source_coords, r=radius)
                ):
                    if neighbors:
                        pair_queries.extend([int(query_id)] * len(neighbors))
                        pair_local_obs.extend(int(value) for value in neighbors)
            local_indices = np.asarray(pair_local_obs, dtype=np.int64)
            replies.append({
                "query_ids": np.asarray(pair_queries, dtype=np.int64),
                "observation_ids": obs_ids[local_indices],
                "yprime": yprime[local_indices],
                "eta": eta[local_indices],
                "innovations": innovations[local_indices],
            })

        received_replies = spatial_comm.alltoall(replies)
        by_query: dict[int, list[tuple[int, np.ndarray, np.ndarray, np.ndarray]]] = {
            int(query_id): [] for query_id in query_ids
        }
        for reply in received_replies:
            reply_query_ids = np.asarray(reply["query_ids"], dtype=np.int64)
            reply_obs_ids = np.asarray(reply["observation_ids"], dtype=np.int64)
            reply_yprime = np.asarray(reply["yprime"], dtype=np.float64)
            reply_eta = np.asarray(reply["eta"], dtype=np.float64)
            reply_innovations = np.asarray(reply["innovations"], dtype=np.float64)
            for index, query_id in enumerate(reply_query_ids):
                by_query[int(query_id)].append((
                    int(reply_obs_ids[index]), reply_yprime[index],
                    reply_eta[index], reply_innovations[index],
                ))

        groups: dict[tuple[int, ...], list[int]] = {}
        terms: dict[tuple[int, ...], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for query_id in query_ids:
            records = sorted(by_query[int(query_id)], key=lambda item: item[0])
            if not records:
                continue
            signature = tuple(item[0] for item in records)
            if len(set(signature)) != len(signature):
                raise ValueError("an observation row is owned by multiple spatial ranks")
            groups.setdefault(signature, []).append(int(query_id))
            if signature not in terms:
                terms[signature] = (
                    np.vstack([item[1] for item in records]),
                    np.vstack([item[2] for item in records]),
                    np.vstack([item[3] for item in records]),
                )

        patches = []
        for signature, local_rows in groups.items():
            local_rows_array = np.asarray(local_rows, dtype=np.int64) - start
            local_terms = terms[signature]
            patches.append(DistributedLocalPatch(
                local_rows=local_rows_array,
                global_rows=targets[start + local_rows_array],
                observation_ids=np.asarray(signature, dtype=np.int64),
                transform=_analysis_transform(
                    *local_terms, energy_fraction=float(energy_fraction)
                ),
            ))
        yield slice(start, stop), tuple(patches)


def apply_distributed_local_patches(
    analysis_block: np.ndarray,
    forecast_block: np.ndarray,
    patches: tuple[DistributedLocalPatch, ...],
) -> np.ndarray:
    """Apply locally owned transforms without a global row-membership scan."""

    analysis = np.asarray(analysis_block)
    forecast = np.asarray(forecast_block)
    if analysis.shape != forecast.shape or analysis.ndim != 2:
        raise ValueError("analysis and forecast blocks must have matching 2-D shapes")
    for patch in patches:
        analysis[patch.local_rows, :] = (
            forecast[patch.local_rows, :] @ patch.transform
        )
    return analysis
