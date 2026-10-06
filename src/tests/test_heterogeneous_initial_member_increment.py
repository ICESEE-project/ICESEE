# ==============================================================================
# @des: Tests for heterogeneous-state-layout support in
# generate_initial_member_increment / resolve_variable_block_sizes
# (src/run_model_da/_error_generation.py).
#
# Fixes the confirmed-live bug: flowline_1d's default (Mode 0) config sets
# scalar_inputs=['xg'] over vec_inputs=['h','u','xg'] (state layout
# h:NX + u:NX + xg:1, total 2*NX+1) and reached a special branch in
# src/EnKF/_ensemble_initialization.py that computed hdim as a naive
# average (state_size // total_state_param_vars) and produced a noise
# vector far too short to broadcast into the ensemble -- an outright
# crash on the very first real run, not just wrong values.
#
# resolve_variable_block_sizes generalizes the var_nd convention already
# used independently by two other call sites (Lorenz96's and ISSM's
# mode3_runner.py/run_da_issm.py: var_nd = {var: (1 if var in
# scalar_inputs else variable_size) for var in vec_inputs}) into the
# single source of truth for per-variable block length, and
# generate_initial_member_increment now uses it whenever scalar_inputs/
# var_nd are configured -- while staying bit-identical to before for the
# uniform case (see test_generate_initial_member_increment.py, unchanged
# and still passing).
# ==============================================================================
import numpy as np
import pytest

from ICESEE.src.run_model_da._error_generation import (
    generate_initial_member_increment,
    resolve_variable_block_sizes,
)


# --- resolve_variable_block_sizes -----------------------------------------

def test_uniform_layout_no_scalar_or_var_nd():
    sizes = resolve_variable_block_sizes(["a", "b", "c"], 12)
    assert sizes == [4, 4, 4]


def test_scalar_plus_field_flowline_like_layout():
    # h:NX + u:NX + xg:1, NX=4 -> vector_size = 9
    sizes = resolve_variable_block_sizes(
        ["h", "u", "xg"], 9, scalar_inputs=["xg"]
    )
    assert sizes == [4, 4, 1]


def test_multiple_scalar_variables():
    sizes = resolve_variable_block_sizes(
        ["field", "s1", "s2"], 6, scalar_inputs=["s1", "s2"]
    )
    assert sizes == [4, 1, 1]


def test_heterogeneous_field_sizes_via_var_nd():
    # field A: N=5, scalar B: 1, field C: M=3 (N != M)
    sizes = resolve_variable_block_sizes(
        ["A", "B", "C"], 9, scalar_inputs=["B"], var_nd={"C": 3}
    )
    assert sizes == [5, 1, 3]


def test_var_nd_alone_without_scalar_inputs():
    sizes = resolve_variable_block_sizes(
        ["A", "B"], 7, var_nd={"B": 3}
    )
    assert sizes == [4, 3]


def test_all_variables_explicitly_sized_must_match_exactly():
    sizes = resolve_variable_block_sizes(
        ["A", "B"], 5, var_nd={"A": 2, "B": 3}
    )
    assert sizes == [2, 3]

    with pytest.raises(ValueError):
        resolve_variable_block_sizes(["A", "B"], 6, var_nd={"A": 2, "B": 3})


def test_mismatched_layout_raises_clear_error_not_silent_wrong_state():
    # 2 regular vars, explicit total from 1 scalar = 1, remaining = 9-1=8,
    # 8 % 2 == 0 so this actually resolves; use an odd remainder instead.
    with pytest.raises(ValueError):
        resolve_variable_block_sizes(
            ["h", "u", "xg"], 8, scalar_inputs=["xg"]
        )  # remaining=7, 2 regular vars, 7 % 2 != 0


def test_offsets_do_not_overlap():
    sizes = resolve_variable_block_sizes(
        ["A", "B", "C"], 9, scalar_inputs=["B"], var_nd={"C": 3}
    )
    offsets = np.cumsum([0] + sizes)
    assert list(offsets) == [0, 5, 6, 9]
    # Every block's index range is disjoint and contiguous.
    covered = np.zeros(9, dtype=bool)
    for start, stop in zip(offsets[:-1], offsets[1:]):
        assert not covered[start:stop].any()
        covered[start:stop] = True
    assert covered.all()


# --- generate_initial_member_increment: heterogeneous layouts --------------

def _flowline_like_kwargs(nx=4, sig_q=(0.02, 0.02, 0.0), base_seed=42):
    return {
        "total_state_param_vars": 3,
        "vec_inputs": ["h", "u", "xg"],
        "scalar_inputs": ["xg"],
        "sig_Q": list(sig_q),
        "Lx": 10.0,
        "Ly": 10.0,
        "base_seed": base_seed,
        "random_field_method": "fft",
    }


def test_flowline_like_layout_produces_correctly_sized_increment():
    nx = 4
    kwargs = _flowline_like_kwargs(nx=nx)
    increment, raw = generate_initial_member_increment(
        nx, kwargs, ensemble_id=0, vector_size=2 * nx + 1
    )
    assert increment.shape == (2 * nx + 1,)
    assert raw.shape == (2 * nx + 1,)
    # h and u blocks (sig_Q nonzero) perturbed; xg block (sig_Q=0.0) is zero.
    assert np.any(increment[:nx] != 0.0)
    assert np.any(increment[nx:2 * nx] != 0.0)
    np.testing.assert_array_equal(increment[2 * nx:], [0.0])


def test_scalar_block_gets_nonzero_raw_field_when_sig_q_nonzero():
    nx = 4
    kwargs = _flowline_like_kwargs(nx=nx, sig_q=(0.0, 0.0, 0.05))
    increment, raw = generate_initial_member_increment(
        nx, kwargs, ensemble_id=0, vector_size=2 * nx + 1
    )
    # h/u blocks zeroed (sig_Q=0), xg block (the scalar) nonzero and scaled.
    np.testing.assert_array_equal(increment[: 2 * nx], np.zeros(2 * nx))
    assert increment[2 * nx] != 0.0
    assert raw[2 * nx] != 0.0
    np.testing.assert_allclose(increment[2 * nx], 0.05 * raw[2 * nx])


def test_zero_sig_q_scalar_does_not_affect_neighboring_blocks():
    nx = 4
    kwargs = _flowline_like_kwargs(nx=nx, sig_q=(0.02, 0.02, 0.0))
    increment, _ = generate_initial_member_increment(
        nx, kwargs, ensemble_id=1, vector_size=2 * nx + 1
    )
    assert increment[2 * nx] == 0.0
    assert not np.allclose(increment[:2 * nx], 0.0)


def test_reproducible_for_same_seed_member_and_layout():
    nx = 4
    kwargs = _flowline_like_kwargs(nx=nx)
    inc_1, _ = generate_initial_member_increment(nx, kwargs, ensemble_id=3, vector_size=2 * nx + 1)
    inc_2, _ = generate_initial_member_increment(nx, kwargs, ensemble_id=3, vector_size=2 * nx + 1)
    np.testing.assert_array_equal(inc_1, inc_2)


def test_distinct_perturbation_for_distinct_members():
    nx = 4
    kwargs = _flowline_like_kwargs(nx=nx)
    inc_a, _ = generate_initial_member_increment(nx, kwargs, ensemble_id=0, vector_size=2 * nx + 1)
    inc_b, _ = generate_initial_member_increment(nx, kwargs, ensemble_id=1, vector_size=2 * nx + 1)
    assert not np.allclose(inc_a, inc_b)


def test_variable_blocks_use_independent_streams():
    nx = 4
    kwargs = _flowline_like_kwargs(nx=nx, sig_q=(1.0, 1.0, 1.0))
    _, raw = generate_initial_member_increment(nx, kwargs, ensemble_id=0, vector_size=2 * nx + 1)
    h_block, u_block = raw[:nx], raw[nx:2 * nx]
    assert not np.allclose(h_block, u_block)


def test_sig_q_scaling_is_per_variable_and_linear():
    nx = 4
    small = _flowline_like_kwargs(nx=nx, sig_q=(0.01, 0.01, 0.01))
    large = _flowline_like_kwargs(nx=nx, sig_q=(0.02, 0.02, 0.02))
    inc_small, _ = generate_initial_member_increment(nx, small, ensemble_id=2, vector_size=2 * nx + 1)
    inc_large, _ = generate_initial_member_increment(nx, large, ensemble_id=2, vector_size=2 * nx + 1)
    np.testing.assert_allclose(inc_large, 2.0 * inc_small)


def test_size_mismatch_between_layout_and_vector_size_raises_clear_error():
    nx = 4
    kwargs = _flowline_like_kwargs(nx=nx)
    with pytest.raises(ValueError):
        # Wrong vector_size: declares one fewer element than h+u+xg needs.
        generate_initial_member_increment(nx, kwargs, ensemble_id=0, vector_size=2 * nx)


def test_uniform_case_is_unaffected_when_scalar_inputs_and_var_nd_absent():
    # Same vec_inputs present, but no scalar_inputs/var_nd configured:
    # must fall back to the exact pre-existing uniform-hdim behavior.
    nx = 4
    kwargs = {
        "total_state_param_vars": 3,
        "vec_inputs": ["h", "u", "xg"],
        "sig_Q": [0.02, 0.02, 0.02],
        "Lx": 10.0,
        "Ly": 10.0,
        "base_seed": 42,
        "random_field_method": "fft",
    }
    increment, _ = generate_initial_member_increment(
        nx, kwargs, ensemble_id=0, vector_size=3 * nx
    )
    assert increment.shape == (3 * nx,)
