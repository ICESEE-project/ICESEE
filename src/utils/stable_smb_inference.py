"""Stable, drop-in SMB inference from assimilated H, u, and v.

Expected units after applying ``mesh_coordinate_scale_to_m``:
    coordinates : m
    thickness   : m
    u, v        : m / model-time-unit (normally m / yr)
    SMB         : m / model-time-unit (normally m / yr)

The routine estimates *apparent* mass balance

    dH/dt + div(H * u_bar)

so u and v should be depth-averaged velocities. If basal mass balance is
nonzero, the result cannot be interpreted as surface mass balance alone.
"""

import numpy as np
from scipy.spatial import cKDTree
from scipy import sparse
from scipy.sparse.linalg import factorized

try:
    # Normal package layout: src/utils/stable_smb_inference.py
    from .localization import get_mesh_coordinates as _registry_coordinates
except ImportError:
    try:
        # Also supports running the module directly from src/utils.
        from localization import get_mesh_coordinates as _registry_coordinates
    except ImportError:
        _registry_coordinates = None


def _state_slice(vec_inputs, name, hdim):
    """Return the slice for an equal-sized, block-ordered state variable."""
    aliases = {
        "h": {"h", "thickness", "ice_thickness"},
        "u": {"u", "vx", "velocity_x", "vel_x", "v_x"},
        "v": {"v", "vy", "velocity_y", "vel_y", "v_y"},
        "smb": {"smb", "surface_mass_balance", "mass_balance"},
    }
    lowered = [str(v).lower() for v in vec_inputs]
    for i, key in enumerate(lowered):
        if key in aliases[name]:
            return slice(i * hdim, (i + 1) * hdim)
    return None


def _configured_smb_source(vec_inputs, observed_params):
    """Select the SMB source from the declared state/observation configuration."""
    aliases = {"smb", "surface_mass_balance", "mass_balance"}
    state_names = {str(value).lower() for value in vec_inputs}
    observed_names = {str(value).lower() for value in observed_params}
    has_smb_state = bool(state_names & aliases)
    smb_is_observed = bool(observed_names & aliases)

    if smb_is_observed and not has_smb_state:
        raise ValueError(
            "SMB is listed in observed_params but is absent from vec_inputs"
        )
    if not has_smb_state:
        return "disabled"
    return "observations" if smb_is_observed else "physics"


def _weighted_meshless_divergence(coords_m, fx, fy, k=24):
    """Weighted quadratic MLS estimate of d(fx)/dx + d(fy)/dy.

    ``fx`` and ``fy`` may be (n,) or (n, nens). The returned array has the
    corresponding shape. A quadratic local reconstruction reduces the
    one-sided truncation error that otherwise appears along mesh boundaries.
    """
    coords_m = np.asarray(coords_m, dtype=float)
    fx = np.asarray(fx, dtype=float)
    fy = np.asarray(fy, dtype=float)
    squeeze = fx.ndim == 1
    if squeeze:
        fx = fx[:, None]
        fy = fy[:, None]

    n = coords_m.shape[0]
    if n < 3:
        out = np.zeros_like(fx)
        return out[:, 0] if squeeze else out

    # Five polynomial coefficients are fitted, so use at least eight points.
    kk = min(max(int(k), 8) + 1, n)
    tree = cKDTree(coords_m)
    dist, nbr = tree.query(coords_m, k=kk)
    dist = dist[:, 1:]
    nbr = nbr[:, 1:]

    dx = coords_m[nbr, 0] - coords_m[:, None, 0]
    dy = coords_m[nbr, 1] - coords_m[:, None, 1]

    # Scale every stencil before forming its quadratic terms. Coordinates in
    # metres would otherwise make the quadratic system poorly conditioned.
    bandwidth = np.maximum(dist[:, -1], np.finfo(float).eps)
    xn = dx / bandwidth[:, None]
    yn = dy / bandwidth[:, None]
    design = np.stack((xn, yn, xn * xn, xn * yn, yn * yn), axis=2)
    weight = np.exp(-4.0 * (dist / bandwidth[:, None]) ** 2)

    normal = np.einsum("nki,nk,nkj->nij", design, weight, design)
    trace = np.trace(normal, axis1=1, axis2=2)
    ridge = 1.0e-9 * np.maximum(trace / 5.0, np.finfo(float).eps)
    normal += ridge[:, None, None] * np.eye(5)[None, :, :]

    dfx = fx[nbr, :] - fx[:, None, :]
    dfy = fy[nbr, :] - fy[:, None, :]
    rhs_fx = np.einsum("nki,nk,nke->nie", design, weight, dfx)
    rhs_fy = np.einsum("nki,nk,nke->nie", design, weight, dfy)

    try:
        coef_fx = np.linalg.solve(normal, rhs_fx)
        coef_fy = np.linalg.solve(normal, rhs_fy)
    except np.linalg.LinAlgError:
        inverse = np.linalg.pinv(normal, rcond=1.0e-10)
        coef_fx = inverse @ rhs_fx
        coef_fy = inverse @ rhs_fy

    dfx_dx = coef_fx[:, 0, :] / bandwidth[:, None]
    dfy_dy = coef_fy[:, 1, :] / bandwidth[:, None]
    out = dfx_dx + dfy_dy
    return out[:, 0] if squeeze else out


def _normalized_graph_laplacian(coords_m, k=12):
    """Build a symmetric, dimensionless, mean-preserving Laplacian.

    Unlike a symmetric normalized Laplacian, this form annihilates a constant
    field even where node degree changes along a boundary.
    """
    n = coords_m.shape[0]
    kk = min(max(int(k), 2) + 1, n)
    tree = cKDTree(coords_m)
    dist, nbr = tree.query(coords_m, k=kk)
    dist = dist[:, 1:]
    nbr = nbr[:, 1:]

    positive = dist[dist > 0.0]
    ell = np.median(positive) if positive.size else 1.0
    weights = np.exp(-(dist / max(ell, np.finfo(float).eps)) ** 2)

    rows = np.repeat(np.arange(n), nbr.shape[1])
    cols = nbr.reshape(-1)
    vals = weights.reshape(-1)
    wmat = sparse.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    wmat = wmat.maximum(wmat.T)
    degree = np.asarray(wmat.sum(axis=1)).ravel()
    scale = np.median(degree[degree > 0.0]) if np.any(degree > 0.0) else 1.0
    return (sparse.diags(degree) - wmat) / scale


def _regularization_solver(coords_m, lambda_space, lambda_time, graph_k):
    """Factorize the regularized least-squares matrix once."""
    lap = _normalized_graph_laplacian(coords_m, graph_k)
    n = coords_m.shape[0]
    matrix = (
        (1.0 + lambda_time) * sparse.eye(n, format="csc")
        + lambda_space * (lap.T @ lap).tocsc()
    )
    return factorized(matrix)


def _project_smb_to_basis(smb, coords_m, mode="none"):
    """Project member-wise SMB onto a small, physically chosen spatial basis."""
    mode = str(mode).lower()
    if mode in {"none", "off", "false", ""}:
        return smb

    x = coords_m[:, 0]
    y = coords_m[:, 1]
    # Center and scale coordinates to keep the least-squares solve conditioned.
    xs = (x - np.mean(x)) / max(np.std(x), np.finfo(float).eps)
    ys = (y - np.mean(y)) / max(np.std(y), np.finfo(float).eps)

    if mode in {"x_linear", "linear_x"}:
        basis = np.column_stack((np.ones_like(xs), xs))
    elif mode in {"affine", "linear_xy"}:
        basis = np.column_stack((np.ones_like(xs), xs, ys))
    elif mode in {"quadratic", "quadratic_xy"}:
        basis = np.column_stack(
            (np.ones_like(xs), xs, ys, xs * xs, xs * ys, ys * ys)
        )
    else:
        raise ValueError(
            "smb_projection_basis must be one of: none, x_linear, affine, quadratic"
        )

    coefficients, *_ = np.linalg.lstsq(basis, smb, rcond=None)
    return basis @ coefficients


def _append_history(icesee_kwargs, time_now, h, u, v, max_length):
    """Store a bounded rolling history in icesee_kwargs."""
    history = icesee_kwargs.setdefault("_smb_inference_history", [])

    # Restart/repeated-write protection: replace, rather than duplicate, time.
    if history and np.isclose(history[-1][0], time_now):
        history[-1] = (float(time_now), h.copy(), u.copy(), v.copy())
    else:
        history.append((float(time_now), h.copy(), u.copy(), v.copy()))

    if len(history) > max_length:
        del history[:-max_length]
    return history


def apply_smb_physics_correction(
    analysis_vec,
    vec_inputs,
    hdim,
    icesee_kwargs,
    timestep=None,
    model_time=None,
):
    """Infer a stable SMB field and insert it into each ensemble member.

    This is a drop-in replacement for the original function. Pass the physical
    model time explicitly when it is available.

    Important configuration keys in ``icesee_kwargs``:

    physics_smb_inference        bool, default False
    smb_history_length           int,  default 5 (minimum 3)
    smb_divergence_neighbors     int,  default 24
    smb_graph_neighbors          int,  default 12
    smb_spatial_regularization   float, default 25
    smb_temporal_regularization  float, default 4
    smb_blend_factor             float, default 0.35
    smb_inference_start_time     float, default 0
    smb_blend_ramp_time          float, default 0
    smb_spinup_hold_factor       float, default 0; use 1 to hold initial SMB
    smb_projection_basis         str, default 'none'
    mesh_coordinate_scale_to_m   float, default 1.0; use 1000 for km input
    smb_physical_bounds          optional (lower, upper), in m/yr
    dt                           time per model step; falls back to icesee_kwargs['dt']
    """
    if not icesee_kwargs.get("physics_smb_inference", False):
        return analysis_vec

    source = _configured_smb_source(
        vec_inputs, icesee_kwargs.get("observed_params", [])
    )
    icesee_kwargs["_smb_last_source"] = source
    if source == "disabled":
        return analysis_vec

    smb_slice = _state_slice(vec_inputs, "smb", hdim)

    # This supports either pasting the function beside the project's existing
    # get_mesh_coordinates helper, or importing this file as a module and
    # supplying coordinates/getter through icesee_kwargs.
    coords = icesee_kwargs.get("mesh_coordinates")
    if coords is None:
        coordinate_getter = icesee_kwargs.get("mesh_coordinate_getter")
        if coordinate_getter is None:
            coordinate_getter = _registry_coordinates
        if callable(coordinate_getter):
            coords = coordinate_getter(icesee_kwargs)
    if coords is None:
        return analysis_vec
    coords = np.asarray(coords, dtype=float)
    if coords.shape != (hdim, 2):
        raise ValueError(
            f"Expected mesh coordinates with shape {(hdim, 2)}, got {coords.shape}"
        )

    coordinate_scale = float(icesee_kwargs.get("mesh_coordinate_scale_to_m", 1.0))
    coords_m = coords * coordinate_scale

    dt = float(icesee_kwargs.get("dt", 1.0))
    if dt <= 0.0:
        raise ValueError("dt must be positive for SMB inference")

    if model_time is not None:
        time_now = float(model_time)
    elif timestep is None:
        # Backward-compatible fallback, though passing timestep is preferred.
        time_now = float(icesee_kwargs.get("_smb_inference_call_count", 0)) * dt
        icesee_kwargs["_smb_inference_call_count"] = (
            icesee_kwargs.get("_smb_inference_call_count", 0) + 1
        )
    else:
        time_now = float(timestep) * dt

    smb_now = np.asarray(analysis_vec[smb_slice, :], dtype=float)

    # When SMB is declared observed, the EnKF has already conditioned this
    # ensemble on those observations. Do not reuse the same likelihood in the
    # physics solve. Enforce only the predeclared spatial parameterization and
    # independent physical bounds, then retain that posterior.
    if source == "observations":
        corrected = _project_smb_to_basis(
            smb_now,
            coords_m,
            icesee_kwargs.get("smb_projection_basis", "none"),
        )
        bounds = icesee_kwargs.get("smb_physical_bounds")
        if bounds is not None:
            lower, upper = map(float, bounds)
            if lower >= upper:
                raise ValueError("smb_physical_bounds must satisfy lower < upper")
            corrected = np.clip(corrected, lower, upper)
        analysis_vec[smb_slice, :] = corrected
        icesee_kwargs["_smb_regularized_previous"] = corrected.copy()
        return analysis_vec

    h_slice = _state_slice(vec_inputs, "h", hdim)
    u_slice = _state_slice(vec_inputs, "u", hdim)
    v_slice = _state_slice(vec_inputs, "v", hdim)
    if any(state_slice is None for state_slice in (h_slice, u_slice, v_slice)):
        return analysis_vec

    h_now = np.asarray(analysis_vec[h_slice, :], dtype=float)
    u_now = np.asarray(analysis_vec[u_slice, :], dtype=float)
    v_now = np.asarray(analysis_vec[v_slice, :], dtype=float)

    spinup_reference = icesee_kwargs.get("_smb_spinup_reference")
    if spinup_reference is None or np.shape(spinup_reference) != np.shape(smb_now):
        spinup_reference = smb_now.copy()
        icesee_kwargs["_smb_spinup_reference"] = spinup_reference

    history_length = max(int(icesee_kwargs.get("smb_history_length", 5)), 3)
    history = _append_history(
        icesee_kwargs, time_now, h_now, u_now, v_now, history_length
    )
    if len(history) < 3:
        return analysis_vec

    # Early in an assimilation the large corrections to H and velocity are not
    # physical tendencies. Do not interpret that spin-up transient as SMB.
    start_time = float(icesee_kwargs.get("smb_inference_start_time", 0.0))
    if time_now <= start_time:
        hold = np.clip(
            float(icesee_kwargs.get("smb_spinup_hold_factor", 0.0)), 0.0, 1.0
        )
        if hold > 0.0:
            held_smb = (
                (1.0 - hold) * smb_now + hold * spinup_reference
            )
            analysis_vec[smb_slice, :] = _project_smb_to_basis(
                held_smb,
                coords_m,
                icesee_kwargs.get("smb_projection_basis", "none"),
            )
        return analysis_vec

    times = np.asarray([item[0] for item in history], dtype=float)
    centered_times = times - times.mean()
    time_denominator = np.dot(centered_times, centered_times)
    if time_denominator <= np.finfo(float).eps:
        return analysis_vec

    h_stack = np.stack([item[1] for item in history], axis=0)
    u_stack = np.stack([item[2] for item in history], axis=0)
    v_stack = np.stack([item[3] for item in history], axis=0)

    # Least-squares temporal slope is much less sensitive to individual EnKF
    # analysis jumps than a two-snapshot finite difference.
    dhdt = np.tensordot(centered_times, h_stack, axes=(0, 0)) / time_denominator

    # A time-mean flux is consistent with the windowed thickness tendency and
    # reduces spatial derivative noise.
    mean_fx = np.mean(h_stack * u_stack, axis=0)
    mean_fy = np.mean(h_stack * v_stack, axis=0)
    div_flux = _weighted_meshless_divergence(
        coords_m,
        mean_fx,
        mean_fy,
        k=int(icesee_kwargs.get("smb_divergence_neighbors", 24)),
    )
    raw_smb = dhdt + div_flux

    lambda_space = max(
        float(icesee_kwargs.get("smb_spatial_regularization", 25.0)), 0.0
    )
    lambda_time = max(
        float(icesee_kwargs.get("smb_temporal_regularization", 4.0)), 0.0
    )
    graph_k = int(icesee_kwargs.get("smb_graph_neighbors", 12))

    cache_key = (
        coords_m.shape[0],
        float(np.sum(coords_m)),
        float(np.sum(coords_m * coords_m)),
        lambda_space,
        lambda_time,
        graph_k,
    )
    cache = icesee_kwargs.get("_smb_regularization_cache")
    if cache is None or cache[0] != cache_key:
        solver = _regularization_solver(
            coords_m, lambda_space, lambda_time, graph_k
        )
        icesee_kwargs["_smb_regularization_cache"] = (cache_key, solver)
    else:
        solver = cache[1]

    previous = icesee_kwargs.get("_smb_regularized_previous")
    if previous is None or np.shape(previous) != np.shape(smb_now):
        previous = smb_now.copy()

    inferred = np.empty_like(raw_smb)
    for member in range(raw_smb.shape[1]):
        rhs = raw_smb[:, member] + lambda_time * previous[:, member]
        inferred[:, member] = solver(rhs)

    finite = np.isfinite(inferred)
    inferred = np.where(finite, inferred, previous)

    bounds = icesee_kwargs.get("smb_physical_bounds")
    if bounds is not None:
        lower, upper = map(float, bounds)
        if lower >= upper:
            raise ValueError("smb_physical_bounds must satisfy lower < upper")
        inferred = np.clip(inferred, lower, upper)

    blend = np.clip(float(icesee_kwargs.get("smb_blend_factor", 0.35)), 0.0, 1.0)
    ramp_time = max(float(icesee_kwargs.get("smb_blend_ramp_time", 0.0)), 0.0)
    if ramp_time > 0.0:
        blend *= np.clip((time_now - start_time) / ramp_time, 0.0, 1.0)
    corrected = (1.0 - blend) * smb_now + blend * inferred
    corrected = _project_smb_to_basis(
        corrected,
        coords_m,
        icesee_kwargs.get("smb_projection_basis", "none"),
    )
    analysis_vec[smb_slice, :] = corrected

    # Anchor the next temporal solve to what was actually applied. This avoids
    # accumulating an unapplied raw correction during the blend ramp.
    icesee_kwargs["_smb_regularized_previous"] = corrected.copy()
    return analysis_vec
