"""Deterministic random streams shared by ICESEE execution modes.

Member streams are keyed by scientific identity rather than MPI rank or
execution order.  Native distributed adapters can therefore initialize a
member with the same stream used by modes 1 and 2 without depending on their
private runner implementation.
"""

from __future__ import annotations

import numpy as np


def initialization_seed(
    base_seed: int,
    member_id: int,
    variable_index: int = 0,
) -> int:
    """Return ICESEE's stable ensemble-initialization seed.

    The word sequence and namespace constant intentionally preserve the
    established mode-1/mode-2 realization exactly.  Do not include a rank,
    communicator size, or scheduling index in this key.
    """

    return int(
        np.random.SeedSequence(
            [
                int(base_seed) & 0xFFFFFFFF,
                int(member_id) & 0xFFFFFFFF,
                int(variable_index) & 0xFFFFFFFF,
                0x1CE511,
            ]
        ).generate_state(1, dtype=np.uint32)[0]
    )


def initialization_rng(
    base_seed: int,
    member_id: int,
    variable_index: int = 0,
) -> np.random.Generator:
    """Create a rank/order-independent initialization generator."""

    return np.random.default_rng(
        initialization_seed(base_seed, member_id, variable_index)
    )


def process_noise_seed(
    base_seed: int,
    timestep: int,
    member_id: int,
    variable_index: int = 0,
) -> int:
    """Return ICESEE's stable process/background-noise (AR(1)) seed.

    Keyed only by (base_seed, timestep, member_id, variable_index) --
    never by MPI rank, communicator size, or scheduling/round order -- so
    every execution mode derives the identical innovation stream for the
    same scientific identity. The word sequence and namespace constant
    intentionally preserve the realization already established by mode 2's
    forecast path (see ``_process_noise_seed`` in
    ``src/parallelization/_mpi_forecast_functions.py``, which now delegates
    here instead of keeping its own copy), and modes 0 and 2 both now use
    this exact stream (see ``EnsembleKalmanFilter.forecast_step``'s
    "serial" backend in ``src/EnKF/python_enkf/EnKF.py``).

    NOTE on naming: this is deliberately NOT the same stream as
    ``mode3_process_noise_seed`` below, despite the similar name and
    purpose. Native mode-3 adapters (which give every member a fully
    independent forecast, unlike modes 0/2's sequential per-member loop)
    use a different argument order and namespace constant, documented on
    that function, and are not bit-compatible with this one. Do not merge
    the two or have one call the other.
    """

    return int(
        np.random.SeedSequence(
            [
                int(base_seed) & 0xFFFFFFFF,
                int(timestep) & 0xFFFFFFFF,
                int(member_id) & 0xFFFFFFFF,
                int(variable_index) & 0xFFFFFFFF,
                0x1CE5EE,
            ]
        ).generate_state(1, dtype=np.uint32)[0]
    )


def mode3_process_noise_seed(
    base_seed: int,
    member_id: int,
    timestep: int,
    variable_index: int = 0,
) -> int:
    """Return a rank/order-independent forecast-step process-noise seed,
    for native mode-3 adapters only -- NOT the same stream as
    ``process_noise_seed`` above, despite the similar name and purpose.

    Modes 0/2's "serial" forecast backend (``EnsembleKalmanFilter.
    forecast_step`` in ``src/EnKF/python_enkf/EnKF.py``) draws AR(1)
    process noise sequentially, member by member, with each member's own
    state now correctly isolated (fixed in the modes-0/2 process-noise
    correction) but still computed inside one sequential Python loop. A
    native mode-3 adapter instead gives every member its own independent
    forecast call, with no shared loop at all, so it needs its own
    deterministic seed keyed by ``(member_id, timestep, variable_index)``
    -- statistically equivalent (same per-step variance and decorrelation)
    to ``process_noise_seed``'s stream, but never bit-identical to it, and
    must never be made to alias it. Uses a namespace constant distinct
    from both ``initialization_seed``'s ``0x1CE511`` and
    ``process_noise_seed``'s ``0x1CE5EE`` so none of the three stream
    families can ever collide.

    Note the argument order differs from ``process_noise_seed``
    (``member_id`` before ``timestep`` here, the reverse there) -- this
    mirrors each stream's own historical realization exactly; do not
    "fix" one to match the other's order, as that would change the
    seed each already produces.
    """

    return int(
        np.random.SeedSequence(
            [
                int(base_seed) & 0xFFFFFFFF,
                int(member_id) & 0xFFFFFFFF,
                int(timestep) & 0xFFFFFFFF,
                int(variable_index) & 0xFFFFFFFF,
                0x9011CE,
            ]
        ).generate_state(1, dtype=np.uint32)[0]
    )


def mode3_process_noise_rng(
    base_seed: int,
    member_id: int,
    timestep: int,
    variable_index: int = 0,
) -> np.random.Generator:
    """Create a rank/order-independent per-member forecast process-noise
    generator for native mode-3 adapters. See ``mode3_process_noise_seed``
    for why this is a distinct stream from ``process_noise_seed``/no RNG
    helper exists for that one directly (modes 0/2 pass ``rng=`` explicitly
    at the call site instead; see ``src/EnKF/python_enkf/EnKF.py``).
    """

    return np.random.default_rng(
        mode3_process_noise_seed(base_seed, member_id, timestep, variable_index)
    )
