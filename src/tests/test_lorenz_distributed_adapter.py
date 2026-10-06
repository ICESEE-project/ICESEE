"""State-only parity tests for the tiny Lorenz distributed adapter."""

from types import SimpleNamespace

import numpy as np
import pytest

from applications.lorenz_model.lorenz_utils.distributed_adapter import (
    LorenzDistributedAdapter,
    lorenz_rk4_step,
)
from src.parallelization.distributed_adapter import (
    DistributedStateLayout,
    validate_distributed_adapter,
)


class _Comm:
    def __init__(self, rank, descriptors, state_pieces):
        self.rank = rank
        self.descriptors = descriptors
        self.state_pieces = state_pieces

    def Get_rank(self):
        return self.rank

    def allgather(self, value):
        if len(value) == 4:
            return list(self.descriptors)
        return list(self.state_pieces)


def _topology(rank, comm):
    return SimpleNamespace(
        ensemble_slot=0,
        spatial_rank=rank,
        ensemble_groups=1,
        spatial_ranks=3,
        spatial_comm=comm,
    )


def _kwargs():
    return {
        "Nens": 1,
        "nd": 3,
        "u0b": np.array([1.0, 2.0, 3.0]),
        "sigma_96": 10.0,
        "beta_96": 8.0 / 3.0,
        "rho_96": 28.0,
        "dt": 0.01,
    }


def test_lorenz_adapter_satisfies_state_only_contract():
    validate_distributed_adapter(LorenzDistributedAdapter())


def test_distributed_lorenz_forecast_matches_serial_rk4():
    adapter = LorenzDistributedAdapter()
    kwargs = _kwargs()
    descriptors = [
        (0, 1, 3, adapter.layout_id),
        (1, 2, 3, adapter.layout_id),
        (2, 3, 3, adapter.layout_id),
    ]
    pieces = [(0, 1, [1.0]), (1, 2, [2.0]), (2, 3, [3.0])]
    expected = lorenz_rk4_step(kwargs["u0b"], kwargs)

    forecast = []
    for rank in range(3):
        topology = _topology(rank, _Comm(rank, descriptors, pieces))
        layout = adapter.distributed_state_layout(
            topology=topology, icesee_kwargs=kwargs
        )
        local = adapter.initialize_local_member(
            0, layout=layout, topology=topology, icesee_kwargs=kwargs
        )
        forecast.append(
            adapter.forecast_local_member(
                local,
                0,
                0,
                layout=layout,
                topology=topology,
                icesee_kwargs=kwargs,
            )[0]
        )

    np.testing.assert_allclose(forecast, expected, rtol=0.0, atol=1.0e-14)


def test_lorenz_reference_adapter_refuses_large_state_gather():
    adapter = LorenzDistributedAdapter()
    kwargs = _kwargs()
    topology = _topology(
        0,
        _Comm(
            0,
            [(0, 1, 3, adapter.layout_id)],
            [(0, 1, [1.0])],
        ),
    )
    layout = DistributedStateLayout(
        global_size=2048,
        owned_start=0,
        owned_stop=1,
        layout_id=adapter.layout_id,
    )
    with pytest.raises(RuntimeError, match="model-native halo exchange"):
        adapter.forecast_local_member(
            np.array([1.0]),
            0,
            0,
            layout=layout,
            topology=topology,
            icesee_kwargs=kwargs,
        )
