# ==============================================================================
# @des: Level-1 integration test for src/EnKF/_ensemble_initialization.py
# (the modes-0/1 ensemble initializer). Exercises the actual call site --
# not just the underlying generate_initial_member_increment helper -- to
# confirm the historical zero-spread / unscaled-noise bug (a single shared
# RNG created once outside the per-member loop, raw noise added with no
# sig_Q scaling) stays fixed at the integration boundary: ensemble spread
# across members must be nonzero, reproducible for a fixed base_seed, and
# must equal exactly what generate_initial_member_increment produces
# directly for the same (hdim, base_seed, member_id).
# ==============================================================================
import numpy as np

from ICESEE.src.EnKF._ensemble_initialization import ensemble_initialization
from ICESEE.src.run_model_da._error_generation import generate_initial_member_increment


class _ZeroInitModel:
    """Stub model_module: every member starts from an all-zero state, so the
    only contribution to ensemble_vec is the initial-perturbation noise."""

    @staticmethod
    def initialize_ensemble(ens, **icesee_kwargs):
        hdim = icesee_kwargs["nd"] // len(icesee_kwargs["vec_inputs"])
        return {var: np.zeros(hdim) for var in icesee_kwargs["vec_inputs"]}


def _run(tmp_path, base_seed, nens=5, hdim=4, num_state_vars=2, sig_Q=(0.02, 0.02)):
    vec_inputs = [f"v{i}" for i in range(num_state_vars)]
    nd = hdim * num_state_vars
    tmp_path.mkdir(parents=True, exist_ok=True)
    icesee_kwargs = {
        "model_module": _ZeroInitModel(),
        "vec_inputs": vec_inputs,
        "nd": nd,
        "Nens": nens,
        "num_state_vars": num_state_vars,
        "total_state_param_vars": num_state_vars,
        "sig_Q": list(sig_Q),
        "Lx": 10.0,
        "Ly": 10.0,
        "base_seed": base_seed,
        "random_field_method": "fft",
        "data_path": str(tmp_path),
        "nt": 1,
        "default_run": True,
    }
    _, ensemble_vec, *_ = ensemble_initialization(**icesee_kwargs)
    return ensemble_vec, hdim


def test_reproducible_for_fixed_base_seed(tmp_path):
    vec_1, _ = _run(tmp_path / "a", base_seed=11)
    vec_2, _ = _run(tmp_path / "b", base_seed=11)
    np.testing.assert_array_equal(vec_1, vec_2)


def test_ensemble_spread_is_nonzero(tmp_path):
    ensemble_vec, _ = _run(tmp_path, base_seed=11, nens=6)
    spread = ensemble_vec.std(axis=1)
    assert np.all(spread > 0.0)


def test_matches_shared_generator_directly(tmp_path):
    hdim = 4
    num_state_vars = 2
    ensemble_vec, _ = _run(
        tmp_path, base_seed=99, nens=3, hdim=hdim, num_state_vars=num_state_vars
    )

    icesee_kwargs = {
        "total_state_param_vars": num_state_vars,
        "sig_Q": [0.02, 0.02],
        "Lx": 10.0,
        "Ly": 10.0,
        "base_seed": 99,
        "random_field_method": "fft",
    }
    for ens in range(3):
        expected_increment, _ = generate_initial_member_increment(
            hdim, icesee_kwargs, ensemble_id=ens
        )
        np.testing.assert_array_equal(ensemble_vec[:, ens], expected_increment)


def test_increment_scales_with_sig_q(tmp_path):
    small, _ = _run(tmp_path / "small", base_seed=5, sig_Q=(0.01, 0.01))
    large, _ = _run(tmp_path / "large", base_seed=5, sig_Q=(0.02, 0.02))
    np.testing.assert_allclose(large, 2.0 * small)


def test_modes_0_1_and_mode_2_share_the_identical_increment_function():
    # Modes 0/1 (src/EnKF/_ensemble_initialization.py) and mode 2
    # (src/parallelization/_mpi_ensemble_intialization.py) must use
    # mathematically consistent initial-ensemble semantics. Rather than
    # relying on two independent implementations happening to agree, both
    # modules import the exact same function object -- so a member's
    # initial perturbation cannot silently diverge between execution
    # modes as long as this identity holds.
    from ICESEE.src.EnKF import _ensemble_initialization as modes01
    from ICESEE.src.parallelization import _mpi_ensemble_intialization as mode2

    assert modes01.generate_initial_member_increment is mode2.generate_initial_member_increment
