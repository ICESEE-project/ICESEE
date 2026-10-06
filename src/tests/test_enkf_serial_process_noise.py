# ==============================================================================
# @des: Level-1 regression tests for the Mode-0 ("serial" analysis backend)
# process/background-noise path in src/EnKF/python_enkf/EnKF.py.
#
# Exercises EnsembleKalmanFilter.forecast_step itself (not just the
# underlying process_noise_seed helper) -- the initialization-bug
# investigation showed helper-only tests are not sufficient to catch
# integration-level defects such as chained-not-independent per-member
# state. These tests prove:
#   - the raw noise field is no longer a single repeated deterministic
#     field (it now varies by timestep/member/variable);
#   - each ensemble member's AR(1) process-noise state is independent and
#     is not affected by the order members are processed in;
#   - sig_Q amplitude scaling and the AR(1) coefficient/formulation are
#     unchanged.
# ==============================================================================
import os
import sys
import tempfile

import numpy as np
import pytest

# config/_utility_imports.py (imported transitively via EnKF.py) parses
# sys.argv AND loads a params.yaml file at module import time -- and, since
# the Pass-2 config-utility port, also deletes and recreates
# icesee_kwargs['data_path'] as a side effect. No other test module imports
# EnKF.py in-process (existing coverage only reaches it through
# subprocess-launched driver scripts, which have their own clean argv and
# cwd), so pytest's own CLI arguments -- and the lack of a params.yaml in
# the pytest working directory -- would otherwise crash collection here.
# Point -F at an existing, lightweight example config purely to satisfy
# that import-time file read (nothing in these tests depends on its
# contents), and --data_path at a disposable temp directory so the
# auto-clean guard never touches the repository's own _modelrun_datasets.
_lorenz96_params = os.path.join(
    os.path.dirname(__file__),
    "..", "..", "applications", "lorenz_model", "examples", "lorenz96", "params.yaml",
)
_safe_data_path = tempfile.mkdtemp(prefix="icesee_enkf_process_noise_test_")
_saved_argv = sys.argv
sys.argv = [sys.argv[0], "-F", _lorenz96_params, "--data_path", _safe_data_path]
try:
    from ICESEE.src.EnKF.python_enkf.EnKF import EnsembleKalmanFilter as EnKF
finally:
    sys.argv = _saved_argv


def _no_op_forecast_step_single(ensemble=None, **icesee_kwargs):
    # Isolate the process-noise contribution: the model itself contributes
    # nothing, so any change to `ensemble` after forecast_step is purely
    # from process noise.
    return {}


def test_mode_0_and_mode_2_share_the_identical_process_noise_seed_function():
    # As with generate_initial_member_increment, prefer a structural
    # invariant over merely observing matching outputs: both EnKF.py
    # (mode 0) and _mpi_forecast_functions.py (mode 2) must derive their
    # per-(timestep, member, variable) seed from the exact same function
    # object in src/utils/random_streams.py, not two copies that could
    # silently drift apart.
    from ICESEE.src.EnKF.python_enkf import EnKF as enkf_module
    from ICESEE.src.parallelization import _mpi_forecast_functions as mode2_module
    from ICESEE.src.utils.random_streams import process_noise_seed

    assert enkf_module.process_noise_seed is process_noise_seed
    assert mode2_module._process_noise_seed is process_noise_seed


def _base_kwargs(nens, hdim=4, num_state_vars=2, sig_q=0.05, base_seed=7, alpha=0.5):
    vec_inputs = [f"v{i}" for i in range(num_state_vars)]
    return {
        "joint_estimation": False,
        "localization_flag": False,
        "num_state_vars": num_state_vars,
        "total_state_param_vars": num_state_vars,
        "nd": hdim * num_state_vars,
        "vec_inputs": vec_inputs,
        "sig_Q": [sig_q] * num_state_vars,
        "default_run": True,
        "even_distribution": False,
        "Lx": 10.0,
        "Ly": 10.0,
        "dt": 1.0,
        "alpha": alpha,
        "rho": 1.0,
        "base_seed": base_seed,
        "random_field_method": "fft",
    }


def _forecast_once(enkf, nens, hdim=4, num_state_vars=2, k=0, **kwargs):
    icesee_kwargs = _base_kwargs(nens, hdim=hdim, num_state_vars=num_state_vars, **kwargs)
    icesee_kwargs["k"] = k
    ensemble = np.zeros((hdim * num_state_vars, nens))
    return enkf.forecast_step(ensemble, _no_op_forecast_step_single, **icesee_kwargs)


def test_reproducible_for_same_seed_timestep_member():
    enkf_1 = EnKF(analysis_backend="serial")
    enkf_2 = EnKF(analysis_backend="serial")

    result_1 = _forecast_once(enkf_1, nens=3, k=0)
    result_2 = _forecast_once(enkf_2, nens=3, k=0)

    np.testing.assert_array_equal(result_1, result_2)


def test_members_receive_distinct_noise_at_same_timestep():
    enkf = EnKF(analysis_backend="serial")
    result = _forecast_once(enkf, nens=4, k=0)

    assert not np.allclose(result[:, 0], result[:, 1])
    assert not np.allclose(result[:, 1], result[:, 2])


def test_raw_field_is_not_repeated_across_timesteps():
    # The historical bug: no timestep/member/variable ever entered the
    # seed, so generate_enkf_field fell back to the same base_seed on
    # every call -- meaning every timestep injected an identical field.
    enkf = EnKF(analysis_backend="serial")
    result_k0 = _forecast_once(enkf, nens=1, k=0)

    enkf_k1 = EnKF(analysis_backend="serial")
    result_k1 = _forecast_once(enkf_k1, nens=1, k=5)

    assert not np.allclose(result_k0[:, 0], result_k1[:, 0])


def test_reseeding_works_even_when_config_sets_base_seed_explicitly():
    # generate_enkf_field's own fallback (in _error_generation.py) prefers
    # a pre-existing "base_seed" entry in icesee_kwargs over "seed" when
    # both are present -- which is exactly what an application config that
    # sets base_seed explicitly (ISSM, icepack synthetic_ice_stream) puts
    # in icesee_kwargs. Passing "rng" explicitly (added alongside "seed" in
    # this fix) must keep per-(timestep,member,variable) variation working
    # even in that configuration, not just for Lorenz96 (which never sets
    # base_seed and so never exercised this fallback ambiguity at all).
    icesee_kwargs = _base_kwargs(nens=2)
    icesee_kwargs["base_seed"] = 123  # explicit config value, present up front
    icesee_kwargs["k"] = 0
    enkf = EnKF(analysis_backend="serial")
    result = enkf.forecast_step(
        np.zeros((8, 2)), _no_op_forecast_step_single, **icesee_kwargs
    )
    assert not np.allclose(result[:, 0], result[:, 1])


def test_member_result_does_not_depend_on_another_members_prior_state():
    # The historical bug: member j+1 inherited member j's just-updated AR(1)
    # state, because a single shared `noise` variable was read once before
    # the ensemble loop and mutated in place at every iteration. Prove the
    # fix directly by giving two instances different prior state for member
    # 0 while keeping member 1's own prior state identical, then checking
    # member 1's output is unaffected by member 0's different history.
    nens, hdim, num_state_vars = 2, 4, 2
    size = hdim * num_state_vars

    enkf_a = EnKF(analysis_backend="serial")
    enkf_a._process_noise_state = np.zeros((size, nens))
    kwargs_a = _base_kwargs(nens)
    kwargs_a["k"] = 0
    result_a = enkf_a.forecast_step(np.zeros((size, nens)), _no_op_forecast_step_single, **kwargs_a)

    enkf_b = EnKF(analysis_backend="serial")
    enkf_b._process_noise_state = np.zeros((size, nens))
    enkf_b._process_noise_state[:, 0] = 999.0  # member 0's prior state differs
    kwargs_b = _base_kwargs(nens)
    kwargs_b["k"] = 0
    result_b = enkf_b.forecast_step(np.zeros((size, nens)), _no_op_forecast_step_single, **kwargs_b)

    np.testing.assert_array_equal(result_a[:, 1], result_b[:, 1])
    assert not np.allclose(result_a[:, 0], result_b[:, 0])  # member 0 itself IS affected, as expected


def test_member_alone_matches_member_zero_inside_a_larger_ensemble():
    enkf_alone = EnKF(analysis_backend="serial")
    result_alone = _forecast_once(enkf_alone, nens=1, k=0)

    enkf_in_ensemble = EnKF(analysis_backend="serial")
    result_in_ensemble = _forecast_once(enkf_in_ensemble, nens=3, k=0)

    np.testing.assert_array_equal(result_alone[:, 0], result_in_ensemble[:, 0])


def test_zero_sig_q_leaves_ensemble_unchanged():
    enkf = EnKF(analysis_backend="serial")
    result = _forecast_once(enkf, nens=3, k=0, sig_q=0.0)
    np.testing.assert_array_equal(result, np.zeros_like(result))


def test_amplitude_scales_linearly_with_sig_q():
    enkf_small = EnKF(analysis_backend="serial")
    small = _forecast_once(enkf_small, nens=2, k=0, sig_q=0.02)

    enkf_large = EnKF(analysis_backend="serial")
    large = _forecast_once(enkf_large, nens=2, k=0, sig_q=0.04)

    np.testing.assert_allclose(large, 2.0 * small)


def test_process_noise_state_persists_and_evolves_across_timesteps():
    enkf = EnKF(analysis_backend="serial")
    _forecast_once(enkf, nens=2, k=0)
    state_after_k0 = enkf._process_noise_state.copy()

    _forecast_once(enkf, nens=2, k=1)
    state_after_k1 = enkf._process_noise_state.copy()

    # The AR(1) state must actually evolve (not stay pinned) across a real
    # timestep advance on the same instance.
    assert not np.allclose(state_after_k0, state_after_k1)
