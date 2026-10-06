# ==============================================================================
# @des: Production execution-mode-3 native adapter for the Lorenz96 model.
# @date: 2026-08-24
# @author: Brian Kyanjo
# ==============================================================================
"""Production execution-mode-3 native adapter for the Lorenz96 application.

Lorenz96's state (x, y, z) has no useful spatial decomposition -- a
3-variable ODE system cannot be split across spatial ranks in any way that
means something physically.  This adapter is registered as mode 3's
production Lorenz96 adapter anyway because mode 3's actual,
application-independent contract is bounded per-rank *ensemble* memory (no
rank ever holds every scheduled member's full time history at once), not
spatial decomposition -- see ``docs/execution-mode-3-design.md``.  Every
native member therefore owns its entire 3-element state locally (no ghosts,
no spatial communication), and this adapter requires
``topology.spatial_ranks == 1``.

The forecast step calls the exact same ``run_model`` function used by
execution modes 0-2
(``applications/lorenz_model/examples/lorenz96/_lorenz96_model.py``), so the
noise-free trajectory contribution is bit-identical across every execution
mode for the same input state.

Modes 0-2's actual EnKF forecast step is not noise-free, though:
``EnsembleKalmanFilter.forecast_step`` (``src/EnKF/python_enkf/EnKF.py``)
adds AR(1)-correlated process noise to every member on every timestep via
``generate_enkf_field``, using the ``alpha``/``rho`` decorrelation
parameters computed once in ``icesee_da_serial.py`` from ``sig_Q``/``dt``.
Skipping that noise would silently under-disperse the mode-3 ensemble
relative to modes 0-2. That serial backend draws noise inside one
sequential Python loop over members -- as of the modes-0/2 process-noise
correction, each member's own AR(1) state is correctly isolated (see
``EnsembleKalmanFilter._process_noise_state`` and
``src/utils/random_streams.py::process_noise_seed``), but it is still one
member's forecast computed *after* another's within a single timestep,
not genuinely independent calls the way this adapter's members are. This
adapter gives every member its own independent, reproducible AR(1) noise
stream keyed by ``(member_id, timestep, variable_index)`` via
``src/utils/random_streams.py::mode3_process_noise_rng`` -- a *separate*
stream from ``process_noise_seed``/``process_noise_rng`` (different
argument order, different namespace constant; see that module for why the
two must not be merged or aliased), statistically equivalent (matching
per-step variance and decorrelation) to modes 0-2's now-corrected stream,
but not bit-identical to it. The AR(1) state is persisted per member
across timesteps in ``NativeDistributedMember.model_context["noise_state"]``,
since native members are long-lived by design.
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from ICESEE.applications.lorenz_model.examples.lorenz96._lorenz96_model import (
    run_model,
)
from ICESEE.src.parallelization.distributed_fields import DistributedFieldRegistry
from ICESEE.src.parallelization.distributed_native_runtime import (
    NativeDistributedMember,
)
from ICESEE.src.run_model_da._error_generation import generate_enkf_field
from ICESEE.src.utils.random_streams import mode3_process_noise_rng


def ar1_parameters(icesee_kwargs: Mapping[str, Any]) -> tuple[float, float]:
    """Return ``(alpha, rho)`` exactly as ``icesee_da_serial.py`` computes them.

    Reused here (rather than requiring the mode-3 runner to inject
    already-computed values) so this adapter stays independently correct
    and testable. Mirrors ``icesee_da_serial.py`` lines ~325-339 exactly:
    ``min_tau=200, max_tau=500``, ``tau=max(max_tau, max(min_tau, dt))``,
    ``alpha=1-dt/tau`` (clamped to ``(0, 1]``, default ``0.5`` if invalid).
    """

    dt = float(icesee_kwargs["dt"])
    min_tau = 200
    max_tau = 500
    tau = max(max_tau, max(min_tau, dt))
    alpha = 1 - dt / tau
    if alpha <= 0 or alpha > 1:
        alpha = 0.5
    n = int(icesee_kwargs["nt"])
    rho = np.sqrt(
        (1 / dt)
        * ((1 - alpha) ** 2)
        * (1 / (n - (2 * alpha) - (n * alpha**2) + (2 * alpha ** (n + 1))))
    )
    return float(alpha), float(rho)


class _LorenzStateField:
    """The full ``[x, y, z]`` state, fully owned locally (no spatial split)."""

    def __init__(self, values: np.ndarray) -> None:
        self.name = "lorenz_xyz"
        self.global_size = 3
        self.owned_start = 0
        self.owned_stop = 3
        self._values = np.asarray(values, dtype=np.float64).reshape(3).copy()

    def read_owned(self) -> np.ndarray:
        return self._values

    def write_owned(self, values: np.ndarray) -> None:
        self._values[:] = np.asarray(values, dtype=np.float64).reshape(3)

    def synchronize_ghosts(self) -> None:
        return None


class LorenzNativeAdapter:
    """Native distributed adapter registered for ``model_name: "lorenz"``."""

    def initialize_native_member(
        self,
        member_id: int,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> NativeDistributedMember:
        if int(topology.spatial_ranks) != 1:
            raise ValueError(
                "LorenzNativeAdapter requires spatial_ranks == 1: Lorenz96's "
                "3-variable state has no useful spatial decomposition"
            )
        initial_state = icesee_kwargs["native_initial_ensemble"][int(member_id)]
        field = _LorenzStateField(initial_state)
        registry = DistributedFieldRegistry([field], layout_id="lorenz96-xyz-v1")
        # Persistent per-member AR(1) process-noise state (one value per
        # state variable), seeded independently of the initial-condition
        # stream -- see module docstring and
        # `random_streams.mode3_process_noise_rng`.
        noise_state = np.zeros(3, dtype=np.float64)
        return NativeDistributedMember(
            int(member_id),
            registry,
            {"field": field, "noise_state": noise_state},
        )

    def forecast_native_member(
        self,
        member: NativeDistributedMember,
        timestep: int,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> None:
        state = member.pack_owned()
        updated = run_model(state, **icesee_kwargs)
        new_state = np.array(
            [updated["x"], updated["y"], updated["z"]], dtype=np.float64
        ).reshape(3)
        new_state += self._process_noise_increment(member, timestep, icesee_kwargs)
        member.unpack_owned(new_state)

    def _process_noise_increment(
        self,
        member: NativeDistributedMember,
        timestep: int,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        """AR(1)-correlated process noise, one independent stream per member.

        See the module docstring for why this is statistically equivalent
        to, but not bit-identical to, modes 0-2's sequential-shared-state
        stream.
        """

        sig_Q = icesee_kwargs.get("sig_Q")
        if sig_Q is None or len(sig_Q) == 0:
            return np.zeros(3, dtype=np.float64)
        alpha, rho = ar1_parameters(icesee_kwargs)
        dt = float(icesee_kwargs["dt"])
        base_seed = icesee_kwargs.get("base_seed", icesee_kwargs.get("seed", 42))
        vec_inputs = icesee_kwargs.get("vec_inputs", ["x", "y", "z"])
        Lx = icesee_kwargs.get("Lx", 1.0)
        Ly = icesee_kwargs.get("Ly", 1.0)
        noise_state = member.model_context["noise_state"]
        increment = np.zeros(3, dtype=np.float64)
        for ii, _key in enumerate(vec_inputs):
            if ii >= len(sig_Q):
                continue
            rng = mode3_process_noise_rng(base_seed, member.member_id, timestep, ii)
            field_kwargs = dict(icesee_kwargs)
            field_kwargs.update(
                {
                    "ii_sig": ii,
                    "Lx_dim": np.sqrt(Lx * Ly),
                    "noise_dim": 1,
                    "hdim": 1,
                    "num_vars": len(vec_inputs),
                    "rng": rng,
                }
            )
            draw = np.asarray(generate_enkf_field(**field_kwargs)).reshape(-1)[0]
            new_noise = alpha * noise_state[ii] + np.sqrt(1 - alpha**2) * draw
            noise_state[ii] = new_noise
            increment[ii] = np.sqrt(dt) * sig_Q[ii] * rho * new_noise
        return increment

    def observe_native_member(
        self,
        member: NativeDistributedMember,
        observation_rows: np.ndarray,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> np.ndarray:
        return member.pack_owned()[observation_rows]

    def finalize_native_analysis(
        self,
        member: NativeDistributedMember,
        forecast_owned: np.ndarray,
        timestep: int,
        *,
        topology: Any,
        icesee_kwargs: Mapping[str, Any],
    ) -> None:
        # Lorenz96 has no post-analysis constraint (no bed/positivity
        # clamp); modes 0-2 also have no `post_analysis_update` for this
        # model (see applications/lorenz_model/lorenz_utils/_lorenz_enkf.py).
        return None
