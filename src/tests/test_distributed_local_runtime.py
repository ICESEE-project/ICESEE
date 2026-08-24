"""Tests for the bounded grouped-local execution-mode-3 runtime."""

from types import SimpleNamespace

import numpy as np

from src.parallelization.distributed_adapter import (
    DistributedAnalysisTargets,
    DistributedObservationCoordinates,
    DistributedStateLayout,
)
from src.parallelization.distributed_local_runtime import (
    apply_grouped_local_analysis_local,
)
from src.parallelization.distributed_runtime import LocalMemberEnsemble
from src.utils.localization import compute_X5_from_matrices


class _SingleComm:
    def allgather(self, value):
        return [value]

    def alltoall(self, value):
        return value

    def bcast(self, value, root=0):
        return value

    def Get_size(self):
        return 1

    def Get_rank(self):
        return 0


def test_grouped_local_runtime_matches_existing_stochastic_transform():
    comm = _SingleComm()
    topology = SimpleNamespace(
        ensemble_comm=comm,
        spatial_comm=comm,
        is_ensemble_root=True,
    )
    rng = np.random.default_rng(811)
    nens = 5
    state = rng.normal(size=(9, nens))
    forecast_observations = rng.normal(size=(4, nens))
    observed = rng.normal(size=4)
    members = {member: state[:, member] for member in range(nens)}
    obs_members = {
        member: forecast_observations[:, member] for member in range(nens)
    }
    local = LocalMemberEnsemble(DistributedStateLayout(9, 0, 9), members)
    target = DistributedAnalysisTargets(
        "thickness",
        np.array([1, 4, 7]),
        np.array([1, 4, 7]),
        np.array([[0.1], [0.5], [0.9]]),
        ("thickness",),
    )
    observations = DistributedObservationCoordinates(
        "thickness",
        np.array([10, 20, 30, 40]),
        np.array([[0.0], [0.45], [0.55], [1.0]]),
    )

    updated = apply_grouped_local_analysis_local(
        local,
        target,
        observations,
        obs_members,
        observed,
        radius=0.16,
        topology=topology,
        icesee_kwargs={"Nens": nens},
        error_mode="legacy_prior_anomalies",
        target_chunk_size=2,
    )
    result = np.column_stack([updated.members[index] for index in range(nens)])
    expected = state.copy()
    yprime = forecast_observations - forecast_observations.mean(axis=1, keepdims=True)
    eta = yprime
    innovations = observed[:, None] - forecast_observations
    obs_coords = observations.coordinates[:, 0]
    for state_row, coordinate in zip(target.local_rows, target.coordinates[:, 0]):
        selected = np.flatnonzero(np.abs(obs_coords - coordinate) <= 0.16)
        if selected.size:
            transform = compute_X5_from_matrices(
                yprime[selected], eta[selected], innovations[selected], nens
            )
            expected[state_row] = state[state_row] @ transform
    np.testing.assert_allclose(result, expected, rtol=1.0e-12, atol=1.0e-12)
