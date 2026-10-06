# ==============================================================================
# @des: Level-1 unit tests for generate_initial_member_increment
# (src/run_model_da/_error_generation.py), the shared initial-ensemble
# perturbation generator used by every execution mode (0/1/2). Covers the
# exact properties the modes-0/1 reconciliation fix was required to
# establish: reproducibility, per-member distinctness, nonzero ensemble
# spread, and correct sig_Q scaling response.
# ==============================================================================
import numpy as np
import pytest

from ICESEE.src.run_model_da._error_generation import generate_initial_member_increment


def _base_kwargs(sig_Q):
    return {
        "total_state_param_vars": len(sig_Q),
        "sig_Q": sig_Q,
        "Lx": 10.0,
        "Ly": 10.0,
        "base_seed": 42,
        "random_field_method": "fft",
    }


def test_reproducible_for_same_base_seed_and_member_id():
    hdim = 8
    kwargs = _base_kwargs([0.01, 0.02, 0.03])

    increment_1, raw_1 = generate_initial_member_increment(hdim, kwargs, ensemble_id=3)
    increment_2, raw_2 = generate_initial_member_increment(hdim, kwargs, ensemble_id=3)

    np.testing.assert_array_equal(increment_1, increment_2)
    np.testing.assert_array_equal(raw_1, raw_2)


def test_distinct_perturbation_for_distinct_member_id():
    hdim = 8
    kwargs = _base_kwargs([0.01, 0.02, 0.03])

    increment_member_0, _ = generate_initial_member_increment(hdim, kwargs, ensemble_id=0)
    increment_member_1, _ = generate_initial_member_increment(hdim, kwargs, ensemble_id=1)

    assert not np.allclose(increment_member_0, increment_member_1)


def test_ensemble_spread_is_nonzero_across_members():
    hdim = 8
    kwargs = _base_kwargs([0.05, 0.05])
    nens = 6

    members = np.stack(
        [
            generate_initial_member_increment(hdim, kwargs, ensemble_id=ens)[0]
            for ens in range(nens)
        ],
        axis=1,
    )

    spread = members.std(axis=1)
    assert np.all(spread > 0.0)


def test_zero_sig_q_gives_zero_increment_but_nonzero_raw_field():
    hdim = 8
    kwargs = _base_kwargs([0.0, 0.0])

    increment, raw = generate_initial_member_increment(hdim, kwargs, ensemble_id=0)

    np.testing.assert_array_equal(increment, np.zeros_like(increment))
    assert np.any(raw != 0.0)


def test_increment_scales_linearly_with_sig_q():
    hdim = 8
    small = _base_kwargs([0.01, 0.01])
    large = _base_kwargs([0.02, 0.02])

    increment_small, raw_small = generate_initial_member_increment(hdim, small, ensemble_id=2)
    increment_large, raw_large = generate_initial_member_increment(hdim, large, ensemble_id=2)

    # Same underlying field (seed depends only on base_seed/member/variable
    # index, not on sig_Q), so doubling sig_Q must exactly double the
    # scaled increment while leaving the raw field unchanged.
    np.testing.assert_array_equal(raw_small, raw_large)
    np.testing.assert_allclose(increment_large, 2.0 * increment_small)


def test_per_variable_blocks_use_independent_seeds():
    # A single-variable field repeated across two "variables" with
    # identical sig_Q would be indistinguishable only if the seed collapsed
    # across variable_index; assert the two hdim-sized blocks differ.
    hdim = 8
    kwargs = _base_kwargs([1.0, 1.0])

    _, raw = generate_initial_member_increment(hdim, kwargs, ensemble_id=0)
    block_0, block_1 = raw[:hdim], raw[hdim:]

    assert not np.allclose(block_0, block_1)


def test_rejects_vector_size_not_a_multiple_of_hdim():
    hdim = 8
    kwargs = _base_kwargs([0.01])

    with pytest.raises(ValueError):
        generate_initial_member_increment(hdim, kwargs, ensemble_id=0, vector_size=hdim + 1)
