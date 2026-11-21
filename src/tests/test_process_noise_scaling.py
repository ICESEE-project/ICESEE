"""Noise amplitude contracts independent of a numerical model."""
import numpy as np
import pytest

from ICESEE.src.utils.process_noise import process_noise_coefficients


@pytest.mark.parametrize("nt", [20, 100, 1640])
def test_legacy_coefficients_reproduce_previous_formula(nt):
    dt = .05
    alpha = 1-dt/500
    expected = np.sqrt((1/dt)*((1-alpha)**2)*(1/(nt-2*alpha-nt*alpha**2+2*alpha**(nt+1))))
    actual_alpha, rho, tau = process_noise_coefficients(dict(dt=dt,nt=nt))
    assert actual_alpha == alpha
    assert rho == expected
    assert tau == 500


@pytest.mark.parametrize("dt", [.01, .05, .1])
def test_per_step_coefficients_do_not_depend_on_final_run_length(dt):
    a = process_noise_coefficients(dict(dt=dt,nt=100,process_noise_normalization="per_step"))
    b = process_noise_coefficients(dict(dt=dt,nt=1640,process_noise_normalization="per_step"))
    assert a == b
    assert a[1] == 1
    # Fresh AR state is still zero; no stationary-state or initial-field draw.
    increment_std = np.sqrt(dt)*.08*np.sqrt(1-a[0]**2)
    assert increment_std > 0


@pytest.mark.parametrize("dt", [0, -1, np.nan, np.inf])
def test_invalid_dt_is_rejected(dt):
    with pytest.raises(ValueError,match="dt"):
        process_noise_coefficients(dict(dt=dt,nt=100))


def test_unknown_policy_is_not_silently_replaced():
    with pytest.raises(ValueError,match="process_noise_normalization"):
        process_noise_coefficients(dict(dt=.05,nt=100,process_noise_normalization="perstep"))


def test_serial_and_parallel_per_step_noise_use_identical_increments(monkeypatch, tmp_path):
    from ICESEE.src.EnKF.python_enkf import EnKF as serial_module
    from ICESEE.src.parallelization import _mpi_forecast_functions as parallel
    def white(**kw):
        return kw["rng"].standard_normal(kw["noise_dim"])
    monkeypatch.setattr(serial_module, "generate_enkf_field", white)
    monkeypatch.setattr(parallel, "generate_enkf_field", white)
    serial = serial_module.EnsembleKalmanFilter(analysis_backend="serial")
    options = dict(dt=.05, nt=1640, process_noise_normalization="per_step",
                   process_noise_schedule="every_step", joint_estimation=False,
                   localization_flag=False, num_state_vars=2, total_state_param_vars=2,
                   nd=8, vec_inputs=["v0","v1"], sig_Q=[.8,.08], base_seed=42,
                   default_run=True, even_distribution=False, Lx=10., Ly=10.,
                   data_path=str(tmp_path), sub_rank=0)
    alpha, rho, _ = process_noise_coefficients(options)
    options.update(alpha=alpha, rho=rho)
    a, b = np.zeros((8,2)), np.zeros((8,2))
    for k in range(3):
        options["k"] = k
        a = serial.forecast_step(a, lambda **kw:{}, **options)
        for mid in range(2):
            b[:,mid] = parallel.add_member_process_noise(b[:,mid], mid, options)
        np.testing.assert_array_equal(a,b)


def test_pig_template_explicitly_selects_per_step_noise():
    from pathlib import Path
    import yaml
    root = Path(__file__).resolve().parents[2]
    config = yaml.safe_load((root/"applications/icepack_model/examples/idealized_pig/params.yaml").read_text())
    options = config["enkf-parameters"]
    assert options["process_noise_normalization"] == "per_step"
    assert process_noise_coefficients(dict(options,dt=.05,nt=1640))[1] == 1.
