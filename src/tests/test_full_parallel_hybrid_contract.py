import numpy as np

from ICESEE.src.parallelization.EnKF_parallel_io import (
    finalize_hybrid_analysis_like_mode1,
)


def test_hybrid_finalization_excludes_friction_then_merges_into_forecast():
    forecast = np.arange(24.0).reshape(6, 4)
    analysis = forecast + 1000.0
    calls = []

    def finalizer(*, analysis_vec, forecast_vec, timestep, icesee_kwargs):
        calls.append((analysis_vec.copy(), forecast_vec, timestep, icesee_kwargs))
        assert analysis_vec.shape == (4, 4)
        assert forecast_vec is None
        assert icesee_kwargs["vec_inputs"] == ["h", "bed"]
        assert icesee_kwargs["nd"] == 4
        assert icesee_kwargs["nd_old"] == 6
        assert icesee_kwargs["excluded_indices"] == [2]
        return analysis_vec + 7.0

    kwargs = {
        "inversion_flag": True,
        "vec_inputs": ["h", "bed", "friction"],
        "vec_inputs_old": ["h", "bed", "friction"],
        "var_nd": {"h": 2, "bed": 2, "friction": 2},
        "var_nd_old": {"h": 2, "bed": 2, "friction": 2},
        "bed_update_domain": "all",
    }
    result = finalize_hybrid_analysis_like_mode1(
        analysis, forecast, 10, kwargs, finalizer
    )

    assert len(calls) == 1
    np.testing.assert_array_equal(result[:4], analysis[:4] + 7.0)
    # Mode 1 does not admit the analyzed friction block before inversion.
    np.testing.assert_array_equal(result[4:], forecast[4:])


def test_hybrid_finalization_passes_reduced_forecast_for_bed_gate():
    forecast = np.arange(24.0).reshape(6, 4)
    analysis = forecast + 100.0

    def finalizer(*, analysis_vec, forecast_vec, timestep, icesee_kwargs):
        np.testing.assert_array_equal(forecast_vec, forecast[:4])
        return analysis_vec

    result = finalize_hybrid_analysis_like_mode1(
        analysis,
        forecast,
        10,
        {
            "inversion_flag": True,
            "vec_inputs": ["h", "bed", "friction"],
            "var_nd": 2,
            "bed_update_domain": "grounded_observed",
        },
        finalizer,
    )
    np.testing.assert_array_equal(result[:4], analysis[:4])
    np.testing.assert_array_equal(result[4:], forecast[4:])
