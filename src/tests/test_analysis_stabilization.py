import numpy as np
import pytest

from ICESEE.src.parallelization._mpi_analysis_functions import (
    stabilize_analysis_increments,
)


def test_analysis_increment_relaxation_and_block_limits():
    forecast = np.zeros((6, 2))
    analysis = np.full((6, 2), 100.0)
    rows = np.arange(6)
    kwargs = {
        "analysis_relaxation_factor": 0.5,
        "analysis_relaxation_factors": {"bed": 1.0},
        "analysis_increment_limits": {"Thickness": 20.0, "bed": 60.0},
    }

    result = stabilize_analysis_increments(
        analysis, forecast, rows, ["Thickness", "bed"], 3, kwargs
    )

    np.testing.assert_allclose(result[:3], 20.0)
    np.testing.assert_allclose(result[3:], 60.0)


def test_analysis_stabilization_rejects_invalid_relaxation():
    with pytest.raises(ValueError, match="must be in"):
        stabilize_analysis_increments(
            np.ones((2, 1)),
            np.zeros((2, 1)),
            np.arange(2),
            ["Thickness"],
            2,
            {"analysis_relaxation_factor": 0.0},
        )
