"""Small reference adapter for execution-mode-3 forecast parity.

The three-variable Lorenz system is intentionally assembled on every spatial
rank during a forecast because all three variables are mutually coupled.  That
is acceptable only for this tiny parity reference and is guarded by a maximum
state size.  Large model adapters must instead perform model-native halo
exchange and must never reconstruct a complete member.
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from src.parallelization.distributed_adapter import (
    DistributedStateLayout,
    contiguous_state_layout,
)


def lorenz_rk4_step(state: np.ndarray, icesee_kwargs: Mapping[str, Any]) -> np.ndarray:
    """Advance the configured three-variable Lorenz system by one RK4 step."""

    state = np.asarray(state, dtype=float)
    if state.shape != (3,):
        raise ValueError("the Lorenz reference state must have shape (3,)")
    sigma = float(icesee_kwargs["sigma_96"])
    beta = float(icesee_kwargs["beta_96"])
    rho = float(icesee_kwargs["rho_96"])
    dt = float(icesee_kwargs["dt"])

    def rhs(values: np.ndarray) -> np.ndarray:
        x, y, z = values
        return np.array(
            [sigma * (y - x), x * (rho - z) - y, x * y - beta * z],
            dtype=float,
        )

    k1 = rhs(state)
    k2 = rhs(state + 0.5 * dt * k1)
    k3 = rhs(state + 0.5 * dt * k2)
    k4 = rhs(state + dt * k3)
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


class LorenzDistributedAdapter:
    """Mode-3 reference adapter used only for state-only parity development."""

    layout_id = "lorenz-three-variable-v1"

    def distributed_state_layout(
        self,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> DistributedStateLayout:
        global_size = int(icesee_kwargs.get("nd", 3))
        if global_size != 3:
            raise ValueError("LorenzDistributedAdapter requires nd=3")
        return contiguous_state_layout(
            global_size,
            topology,
            layout_id=self.layout_id,
        )

    def initialize_local_member(
        self,
        member_id: int,
        *,
        layout: DistributedStateLayout,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        initial_ensemble = icesee_kwargs.get("distributed_initial_ensemble")
        if initial_ensemble is None:
            global_state = np.asarray(icesee_kwargs["u0b"], dtype=float)
        else:
            initial_ensemble = np.asarray(initial_ensemble, dtype=float)
            expected = (layout.global_size, int(icesee_kwargs.get("Nens", 1)))
            if initial_ensemble.shape != expected:
                raise ValueError(
                    "distributed_initial_ensemble must have shape " f"{expected}"
                )
            global_state = initial_ensemble[:, int(member_id)]
        if global_state.shape != (layout.global_size,):
            raise ValueError("u0b must contain the complete three-variable state")
        return np.ascontiguousarray(global_state[layout.owned_slice])

    def forecast_local_member(
        self,
        local_state: np.ndarray,
        member_id: int,
        timestep: int,
        *,
        layout: DistributedStateLayout,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        maximum = int(icesee_kwargs.get("mode3_reference_gather_max_state", 1024))
        if layout.global_size > maximum:
            raise RuntimeError(
                "reference allgather is disabled for large states; implement "
                "model-native halo exchange"
            )
        payload = (
            layout.owned_start,
            layout.owned_stop,
            np.ascontiguousarray(local_state),
        )
        pieces = topology.spatial_comm.allgather(payload)
        global_state = np.empty(layout.global_size, dtype=float)
        cursor = 0
        for start, stop, slab in sorted(pieces, key=lambda item: int(item[0])):
            start, stop = int(start), int(stop)
            slab = np.asarray(slab, dtype=float)
            if start != cursor or slab.shape != (stop - start,):
                raise ValueError("Lorenz spatial slabs do not form an exact partition")
            global_state[start:stop] = slab
            cursor = stop
        if cursor != layout.global_size:
            raise ValueError("Lorenz spatial slabs do not cover the complete state")
        forecast = lorenz_rk4_step(global_state, icesee_kwargs)
        return np.ascontiguousarray(forecast[layout.owned_slice])

    def observe_local_member(
        self,
        local_state: np.ndarray,
        observation_rows: np.ndarray,
        *,
        layout: DistributedStateLayout,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        rows = np.asarray(observation_rows, dtype=np.int64)
        if rows.ndim != 1:
            raise ValueError("observation_rows must be one-dimensional")
        if np.any((rows < layout.owned_start) | (rows >= layout.owned_stop)):
            raise ValueError("observation_rows must be owned by this spatial rank")
        return np.asarray(local_state)[rows - layout.owned_start].copy()

    def finalize_local_analysis(
        self,
        local_forecast: np.ndarray,
        local_analysis: np.ndarray,
        member_id: int,
        timestep: int,
        *,
        layout: DistributedStateLayout,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        analysis = np.asarray(local_analysis, dtype=float)
        if analysis.shape != (layout.owned_size,):
            raise ValueError("local analysis shape does not match the owned slab")
        return np.ascontiguousarray(analysis)
