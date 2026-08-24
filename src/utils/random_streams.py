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
