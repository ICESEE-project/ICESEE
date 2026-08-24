import numpy as np
import pytest

from ICESEE.src.utils.utils import cross_flow_track_mask


def test_cross_flow_tracks_leave_unobserved_gaps():
    x_coord_m = np.arange(0.0, 61_000.0, 1_000.0)
    mask = cross_flow_track_mask(
        x_coord_m, stride_km=30.0, half_width_m=1_000.0
    )

    expected = np.isin(
        x_coord_m,
        [0, 1_000, 29_000, 30_000, 31_000, 59_000, 60_000],
    )
    np.testing.assert_array_equal(mask, expected)
    assert mask.sum() < mask.size


def test_cross_flow_tracks_reject_invalid_stride():
    with pytest.raises(ValueError, match="stride"):
        cross_flow_track_mask([0.0, 1.0], stride_km=0.0)
