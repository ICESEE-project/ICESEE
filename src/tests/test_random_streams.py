from __future__ import annotations

import numpy as np

from src.utils.random_streams import initialization_rng, initialization_seed


def _legacy_initialization_seed(base_seed, member_id, variable_index=0):
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


def test_initialization_seed_preserves_existing_mode_1_mode_2_stream():
    for member_id in (0, 1, 17, 2**32 + 3):
        for variable_index in (0, 1, 9):
            assert initialization_seed(42, member_id, variable_index) == (
                _legacy_initialization_seed(42, member_id, variable_index)
            )


def test_initialization_rng_is_member_keyed_and_order_independent():
    expected = {
        member_id: initialization_rng(81, member_id, 4).normal(size=8)
        for member_id in (0, 3, 9)
    }
    actual = {
        member_id: initialization_rng(81, member_id, 4).normal(size=8)
        for member_id in (9, 0, 3)
    }
    for member_id in expected:
        np.testing.assert_array_equal(actual[member_id], expected[member_id])
    assert not np.array_equal(expected[0], expected[3])
