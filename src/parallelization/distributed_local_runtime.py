"""Bounded grouped-local analysis runtime for execution mode 3.

This module composes the two orthogonal mode-3 communicators without changing
ICESEE's stochastic filter.  Ensemble ranks assemble only locally owned rows;
the first ensemble slot performs sparse spatial localization; and the resulting
small transforms are broadcast to ranks holding the same spatial slab.
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from .distributed_adapter import (
    DistributedAnalysisTargets,
    DistributedObservationCollection,
    DistributedObservationCoordinates,
)
from .distributed_analysis import assemble_selected_ensemble_rows
from .distributed_local_analysis import (
    apply_distributed_local_patches,
    iter_distributed_local_patches,
)
from .distributed_runtime import LocalMemberEnsemble


def _local_observation_terms(
    local_forecast_observations: Mapping[int, np.ndarray],
    observed_values: np.ndarray,
    local_observation_errors: Mapping[int, np.ndarray] | None,
    *,
    number_of_members: int,
    ensemble_comm: Any,
    error_mode: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = np.arange(np.asarray(observed_values).size, dtype=np.int64)
    forecast = assemble_selected_ensemble_rows(
        local_forecast_observations, rows, number_of_members, ensemble_comm
    ).astype(np.float64, copy=False)
    observed = np.asarray(observed_values, dtype=np.float64).ravel()
    if forecast.shape[0] != observed.size:
        raise ValueError("forecast-observation rows and observed values disagree")
    yprime = forecast - forecast.mean(axis=1, keepdims=True)
    mode = str(error_mode).lower()
    if mode == "legacy_prior_anomalies":
        eta = yprime
        innovations = observed[:, None] - forecast
    elif mode in {"stochastic_r", "generated_r"}:
        if local_observation_errors is None:
            raise ValueError(f"{error_mode} requires member-keyed observation errors")
        eta = assemble_selected_ensemble_rows(
            local_observation_errors, rows, number_of_members, ensemble_comm
        ).astype(np.float64, copy=False)
        innovations = observed[:, None] + eta - forecast
    else:
        raise ValueError(
            "error_mode must be stochastic_R, generated_R, or "
            "legacy_prior_anomalies"
        )
    return yprime, eta, innovations


def apply_grouped_local_analysis_local(
    local_forecast: LocalMemberEnsemble,
    target: DistributedAnalysisTargets,
    observation_metadata: (
        DistributedObservationCoordinates | DistributedObservationCollection
    ),
    local_forecast_observations: Mapping[int, np.ndarray],
    observed_values: np.ndarray,
    *,
    radius: float,
    topology: Any,
    icesee_kwargs: Mapping[str, Any],
    error_mode: str,
    local_observation_errors: Mapping[int, np.ndarray] | None = None,
    target_chunk_size: int = 1024,
    energy_fraction: float = 0.999,
) -> LocalMemberEnsemble:
    """Apply one model-declared localized state update with bounded memory.

    The metadata contains all active observation kinds accepted by this target
    for the current cycle. They are analyzed jointly, matching mode 1; applying
    one sequential update per observation kind would define a different filter.
    A higher-level runner calls this function once per variable/state block,
    then invokes the model's geometry-consistency finalizer exactly once.
    """

    number_of_members = int(icesee_kwargs.get("Nens", 1))
    if target.local_rows.size and (
        target.local_rows.min() < 0
        or target.local_rows.max() >= local_forecast.layout.owned_size
    ):
        raise ValueError("analysis targets leave the locally owned state slab")
    if observation_metadata.observation_ids.size != np.asarray(observed_values).size:
        raise ValueError("observation metadata and values have different lengths")
    yprime, eta, innovations = _local_observation_terms(
        local_forecast_observations,
        observed_values,
        local_observation_errors,
        number_of_members=number_of_members,
        ensemble_comm=topology.ensemble_comm,
        error_mode=error_mode,
    )

    output = {
        member_id: np.asarray(values, dtype=np.float64).copy()
        for member_id, values in local_forecast.members.items()
    }
    local_rounds = (
        target.local_rows.size + int(target_chunk_size) - 1
    ) // int(target_chunk_size)
    if topology.is_ensemble_root:
        spatial_rounds = max(
            int(value) for value in topology.spatial_comm.allgather(local_rounds)
        )
    else:
        spatial_rounds = None
    rounds = int(topology.ensemble_comm.bcast(spatial_rounds, root=0))

    patch_iterator = None
    if topology.is_ensemble_root:
        patch_iterator = iter_distributed_local_patches(
            target_global_rows=target.global_rows,
            target_coordinates=target.coordinates,
            observation_ids=observation_metadata.observation_ids,
            observation_coordinates=observation_metadata.coordinates,
            yprime=yprime,
            eta=eta,
            innovations=innovations,
            radius=radius,
            spatial_comm=topology.spatial_comm,
            target_chunk_size=target_chunk_size,
            energy_fraction=energy_fraction,
        )

    for _ in range(rounds):
        root_payload = next(patch_iterator) if patch_iterator is not None else None
        target_slice, patches = topology.ensemble_comm.bcast(root_payload, root=0)
        selected = target.local_rows[target_slice]
        forecast_block = assemble_selected_ensemble_rows(
            local_forecast.members,
            selected,
            number_of_members,
            topology.ensemble_comm,
        )
        analysis_block = forecast_block.copy()
        apply_distributed_local_patches(analysis_block, forecast_block, patches)
        for member_id in local_forecast.member_ids:
            output[member_id][selected] = analysis_block[:, member_id]

    return LocalMemberEnsemble(layout=local_forecast.layout, members=output)
