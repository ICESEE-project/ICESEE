"""Small, explicit diagnostics used by the result views.

Every function here is a direct transcription of a computation in one of the
examples' ``read_results.ipynb`` notebooks (cited per function), so the
numbers the launcher shows are the numbers the notebook prints.
"""
from __future__ import annotations

import numpy as np


def block(index: int, size: int) -> slice:
    """State rows of variable block ``index`` when every block has ``size``
    rows -- the notebooks' ``ensemble_true_state[k*hdim:(k+1)*hdim]``."""
    return slice(index * size, (index + 1) * size)


def magnitude(u, v):
    """``np.sqrt(u**2 + v**2)`` -- velocity magnitude from its components."""
    return np.sqrt(np.asarray(u) ** 2 + np.asarray(v) ** 2)


def nearest_steps(t, times) -> list[int]:
    """``[int(np.argmin(np.abs(t - year))) for year in desired_years]``
    (synthetic_ice_stream/read_results.ipynb, difference maps)."""
    t = np.asarray(t, dtype=float)
    return [int(np.argmin(np.abs(t - float(x)))) for x in times]


def error_stats(truth, estimate) -> dict:
    """RMSE / max |error| / mean error of ``estimate - truth``
    (synthetic_ice_stream/read_results.ipynb, error tables)."""
    error = np.asarray(estimate) - np.asarray(truth)
    return {
        "rmse": float(np.sqrt(np.mean(error ** 2))),
        "max_error": float(np.max(np.abs(error))),
        "mean_error": float(np.mean(error)),
    }


def rmse_over_time(truth, estimate) -> np.ndarray:
    """Domain RMSE at every time step: ``sqrt(mean((estimate - truth)**2,
    axis=0))`` (synthetic_ice_stream/read_results.ipynb, SMB RMSE)."""
    n = min(np.asarray(truth).shape[1], np.asarray(estimate).shape[1])
    diff = np.asarray(estimate)[:, :n] - np.asarray(truth)[:, :n]
    return np.sqrt(np.mean(diff ** 2, axis=0))


def grid_point_index(nx: int, ny: int, hdim: int, profile: str = "middle") -> int:
    """The ``h_indx`` the Icepack/ISSM notebooks compute for a named profile
    location, using ICESEE's own ``get_grid_dimensions`` / ``midindices``."""
    from ICESEE.src.utils.tools import get_grid_dimensions, midindices

    mx, my = get_grid_dimensions(nx, ny, hdim)
    ix, iy = midindices(mx, my)
    if profile == "begining" or profile == "beginning":
        return 0
    if profile == "end":
        return (mx - 1) * (my - 1) + 0
    return ((ix // 2) * (my - 1)) + (iy) // 2


def grid_point_xy(h_indx: int, nx: int, ny: int, hdim: int, Lx: float, Ly: float) -> tuple[float, float]:
    """The notebooks' approximate physical location of ``h_indx`` on the
    nodal grid (``row, col = h_indx // (my-1), h_indx % (mx-1)``), used only
    to mark the profile point on maps."""
    from ICESEE.src.utils.tools import get_grid_dimensions

    mx, my = get_grid_dimensions(nx, ny, hdim)
    row, col = h_indx // (my - 1), h_indx % (mx - 1)
    dx, dy = Lx / (my - 1), Ly / (mx - 1)
    return col * dx, row * dy
