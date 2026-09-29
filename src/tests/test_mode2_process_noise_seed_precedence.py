# ==============================================================================
# @des: Regression test for the Mode-2 seed-precedence fix in
# add_member_process_noise (src/parallelization/_mpi_forecast_functions.py).
#
# generate_enkf_field's own fallback (_error_generation.py) prefers a
# pre-existing "base_seed" entry in icesee_kwargs over a freshly-computed
# "seed" when both are present. config/_utility_imports.py unconditionally
# sets icesee_kwargs["base_seed"] (default 42 if not configured) for every
# application -- so this ambiguity was live for every Mode-2 run, not just
# configs that set base_seed explicitly in their own params.yaml. Fixed by
# also passing "rng" explicitly, mirroring the identical fix already applied
# to Mode 0's forecast step (EnKF.py) and generate_initial_member_increment.
# ==============================================================================
import numpy as np
import pytest

from ICESEE.src.parallelization._mpi_forecast_functions import add_member_process_noise


def _base_kwargs(tmp_path, base_seed=42, sig_q=0.05, num_state_vars=2, hdim=4):
    vec_inputs = [f"v{i}" for i in range(num_state_vars)]
    return {
        "k": 0,
        "process_noise_schedule": "every_step",
        "num_state_vars": num_state_vars,
        "total_state_param_vars": num_state_vars,
        "joint_estimation": False,
        "localization_flag": False,
        "vec_inputs": vec_inputs,
        "sig_Q": [sig_q] * num_state_vars,
        "alpha": 0.5,
        "rho": 1.0,
        "dt": 1.0,
        "Lx": 10.0,
        "Ly": 10.0,
        "base_seed": base_seed,  # always present in real icesee_kwargs, default 42
        "random_field_method": "fft",
        "data_path": str(tmp_path),
        "sub_rank": 0,
    }


def test_members_receive_distinct_noise_with_default_base_seed_present(tmp_path):
    # base_seed is present (as it always is via config loading) -- this is
    # exactly the configuration that previously collapsed every member's
    # process noise to the same field.
    ensemble_a = np.zeros(8)
    ensemble_b = np.zeros(8)
    add_member_process_noise(ensemble_a, 0, _base_kwargs(tmp_path / "a"))
    add_member_process_noise(ensemble_b, 1, _base_kwargs(tmp_path / "b"))
    assert not np.allclose(ensemble_a, ensemble_b)


def test_reproducible_for_same_seed_timestep_member(tmp_path):
    result_1 = np.zeros(8)
    result_2 = np.zeros(8)
    add_member_process_noise(result_1, 0, _base_kwargs(tmp_path / "a"))
    add_member_process_noise(result_2, 0, _base_kwargs(tmp_path / "b"))
    np.testing.assert_array_equal(result_1, result_2)


def test_varies_by_timestep(tmp_path):
    result_k0 = np.zeros(8)
    result_k1 = np.zeros(8)
    add_member_process_noise(result_k0, 0, _base_kwargs(tmp_path / "a"))
    kwargs_k1 = _base_kwargs(tmp_path / "b")
    kwargs_k1["k"] = 5
    add_member_process_noise(result_k1, 0, kwargs_k1)
    assert not np.allclose(result_k0, result_k1)


def test_zero_sig_q_leaves_ensemble_unchanged(tmp_path):
    result = np.zeros(8)
    add_member_process_noise(result, 0, _base_kwargs(tmp_path, sig_q=0.0))
    np.testing.assert_array_equal(result, np.zeros(8))


def test_amplitude_scales_linearly_with_sig_q(tmp_path):
    small = np.zeros(8)
    large = np.zeros(8)
    add_member_process_noise(small, 0, _base_kwargs(tmp_path / "small", sig_q=0.02))
    add_member_process_noise(large, 0, _base_kwargs(tmp_path / "large", sig_q=0.04))
    np.testing.assert_allclose(large, 2.0 * small)


def test_reseeding_works_even_when_base_seed_is_explicitly_configured(tmp_path):
    # Same property as the default case above, but with a non-default
    # base_seed explicitly set (e.g. ISSM/icepack synthetic_ice_stream-style
    # configs) -- confirms the fix does not depend on base_seed happening
    # to equal the module's own internal default.
    result_a = np.zeros(8)
    result_b = np.zeros(8)
    add_member_process_noise(result_a, 0, _base_kwargs(tmp_path / "a", base_seed=123))
    add_member_process_noise(result_b, 1, _base_kwargs(tmp_path / "b", base_seed=123))
    assert not np.allclose(result_a, result_b)
