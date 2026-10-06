from __future__ import annotations

import numpy as np

from src.utils.random_streams import (
    initialization_rng,
    initialization_seed,
    mode3_process_noise_rng,
    mode3_process_noise_seed,
    process_noise_seed,
)


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


def _legacy_process_noise_seed(base_seed, timestep, ens_id, variable_index):
    # The exact arithmetic mode 2's _mpi_forecast_functions.py used to keep
    # privately before it was generalized into process_noise_seed here.
    words = np.random.SeedSequence(
        [
            int(base_seed) & 0xFFFFFFFF,
            int(timestep) & 0xFFFFFFFF,
            int(ens_id) & 0xFFFFFFFF,
            int(variable_index) & 0xFFFFFFFF,
            0x1CE5EE,
        ]
    ).generate_state(1, dtype=np.uint32)
    return int(words[0])


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


def test_process_noise_seed_preserves_existing_mode_2_stream():
    for timestep in (0, 1, 17):
        for ens_id in (0, 2, 5):
            for variable_index in (0, 1, 2):
                assert process_noise_seed(42, timestep, ens_id, variable_index) == (
                    _legacy_process_noise_seed(42, timestep, ens_id, variable_index)
                )


def test_process_noise_seed_is_reproducible():
    assert process_noise_seed(42, 3, 1, 0) == process_noise_seed(42, 3, 1, 0)


def test_process_noise_seed_varies_by_member():
    seeds = {ens_id: process_noise_seed(42, 3, ens_id, 0) for ens_id in range(4)}
    assert len(set(seeds.values())) == len(seeds)


def test_process_noise_seed_varies_by_timestep():
    seeds = {k: process_noise_seed(42, k, 1, 0) for k in range(4)}
    assert len(set(seeds.values())) == len(seeds)


def test_process_noise_seed_varies_by_variable():
    seeds = {ii: process_noise_seed(42, 3, 1, ii) for ii in range(4)}
    assert len(set(seeds.values())) == len(seeds)


def test_process_noise_seed_does_not_depend_on_rank_or_call_order():
    # No rank/scheduling argument exists in the signature at all; calling
    # it in any order for the same (base_seed, k, member, variable) must
    # give the same result, mirroring initialization_seed's invariant.
    forward = [process_noise_seed(42, 5, ens_id, 0) for ens_id in range(4)]
    backward = [process_noise_seed(42, 5, ens_id, 0) for ens_id in reversed(range(4))]
    assert forward == list(reversed(backward))


# --- mode3_process_noise_seed: a deliberately separate stream, used only
# by native mode-3 adapters (e.g. Lorenz96's distributed_native_adapter.py).
# Not bit-compatible with process_noise_seed above, and must stay that way.

def _legacy_mode3_process_noise_seed(base_seed, member_id, timestep, variable_index):
    # The exact historical arithmetic from test/ICESEE's random_streams.py,
    # preserved here so a future edit can't silently change this stream's
    # realization without a test catching it.
    words = np.random.SeedSequence(
        [
            int(base_seed) & 0xFFFFFFFF,
            int(member_id) & 0xFFFFFFFF,
            int(timestep) & 0xFFFFFFFF,
            int(variable_index) & 0xFFFFFFFF,
            0x9011CE,
        ]
    ).generate_state(1, dtype=np.uint32)
    return int(words[0])


def test_mode3_process_noise_seed_preserves_historical_stream():
    for member_id in (0, 1, 5):
        for timestep in (0, 3, 17):
            for variable_index in (0, 1, 2):
                assert mode3_process_noise_seed(42, member_id, timestep, variable_index) == (
                    _legacy_mode3_process_noise_seed(42, member_id, timestep, variable_index)
                )


def test_mode3_process_noise_seed_is_reproducible():
    assert mode3_process_noise_seed(42, 1, 3, 0) == mode3_process_noise_seed(42, 1, 3, 0)


def test_mode3_process_noise_seed_varies_by_member_timestep_and_variable():
    by_member = {m: mode3_process_noise_seed(42, m, 3, 0) for m in range(4)}
    by_timestep = {k: mode3_process_noise_seed(42, 1, k, 0) for k in range(4)}
    by_variable = {v: mode3_process_noise_seed(42, 1, 3, v) for v in range(4)}
    assert len(set(by_member.values())) == 4
    assert len(set(by_timestep.values())) == 4
    assert len(set(by_variable.values())) == 4


def test_mode3_process_noise_rng_matches_its_own_seed():
    seed = mode3_process_noise_seed(42, 2, 7, 1)
    expected = np.random.default_rng(seed).normal(size=5)
    actual = mode3_process_noise_rng(42, 2, 7, 1).normal(size=5)
    np.testing.assert_array_equal(actual, expected)


def test_mode3_stream_never_aliases_the_modes_0_2_stream():
    # The whole point of the two separate names: same inputs, different
    # functions, must not collide -- neither in value nor in argument
    # order (deliberately swapped: member before timestep here, timestep
    # before member in process_noise_seed).
    base_seed, a, b, variable_index = 42, 3, 5, 0
    assert mode3_process_noise_seed(base_seed, a, b, variable_index) != (
        process_noise_seed(base_seed, a, b, variable_index)
    )
    assert mode3_process_noise_seed is not process_noise_seed
