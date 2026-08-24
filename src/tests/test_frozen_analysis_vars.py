import numpy as np

from ICESEE.src.utils.localization import restore_frozen_analysis_vars


def test_restore_frozen_analysis_vars_restores_only_requested_block():
    forecast = np.arange(24, dtype=float).reshape(6, 4)
    analysis = forecast + 100.0
    rows = np.arange(6)

    result = restore_frozen_analysis_vars(
        analysis.copy(),
        forecast,
        rows,
        ["Thickness", "bed", "coefficient"],
        hdim=2,
        frozen_vars=["coefficient"],
    )

    np.testing.assert_allclose(result[:4], analysis[:4])
    np.testing.assert_allclose(result[4:], forecast[4:])


def test_restore_frozen_analysis_vars_handles_partitioned_rows():
    forecast = np.arange(8, dtype=float).reshape(2, 4)
    analysis = forecast + 50.0

    result = restore_frozen_analysis_vars(
        analysis.copy(),
        forecast,
        global_rows=np.array([4, 5]),
        vec_inputs=["Thickness", "bed", "coefficient"],
        hdim=2,
        frozen_vars=["Coefficient"],
    )

    np.testing.assert_allclose(result, forecast)
