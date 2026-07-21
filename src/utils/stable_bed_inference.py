"""Stable post-analysis regularization for an inferred glacier/ice-sheet bed.

Unlike SMB, bed is not obtained from the mass-continuity residual. This routine
stabilizes the bed increments produced by an augmented-state filter or inverse
step. It regularizes *increments*, not absolute bed elevation, so real prior
topography is not progressively flattened.
"""

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.sparse.linalg import factorized

try:
    from .localization import get_mesh_coordinates as _registry_coordinates
except ImportError:
    try:
        from localization import get_mesh_coordinates as _registry_coordinates
    except ImportError:
        _registry_coordinates = None


def _bed_state_slice(vec_inputs, aliases, hdim):
    lowered = [str(value).lower() for value in vec_inputs]
    for index, value in enumerate(lowered):
        if value in aliases:
            return slice(index * hdim, (index + 1) * hdim)
    return None


def _bed_graph_laplacian(coords_m, k=12):
    """Symmetric, dimensionless graph Laplacian satisfying L @ 1 == 0."""
    n = coords_m.shape[0]
    if n < 2:
        return sparse.csr_matrix((n, n))

    kk = min(max(int(k), 2) + 1, n)
    distance, neighbors = cKDTree(coords_m).query(coords_m, k=kk)
    distance = distance[:, 1:]
    neighbors = neighbors[:, 1:]

    positive = distance[distance > 0.0]
    length_scale = np.median(positive) if positive.size else 1.0
    weights = np.exp(
        -(distance / max(length_scale, np.finfo(float).eps)) ** 2
    )

    rows = np.repeat(np.arange(n), neighbors.shape[1])
    columns = neighbors.reshape(-1)
    values = weights.reshape(-1)
    adjacency = sparse.coo_matrix(
        (values, (rows, columns)), shape=(n, n)
    ).tocsr()
    adjacency = adjacency.maximum(adjacency.T)

    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    scale = np.median(degree[degree > 0.0]) if np.any(degree > 0.0) else 1.0
    return (sparse.diags(degree) - adjacency) / scale


def _bed_increment_solver(coords_m, spatial_strength, graph_neighbors):
    laplacian = _bed_graph_laplacian(coords_m, graph_neighbors)
    n = coords_m.shape[0]
    matrix = (
        sparse.eye(n, format="csc")
        + spatial_strength * (laplacian.T @ laplacian).tocsc()
    )
    return factorized(matrix)


def _project_bed_increment(increment, coords_m, mode="none"):
    """Optionally retain only identifiable large-scale increment modes."""
    mode = str(mode).lower()
    if mode in {"none", "off", "false", ""}:
        return increment

    x = coords_m[:, 0]
    y = coords_m[:, 1]
    xs = (x - np.mean(x)) / max(np.std(x), np.finfo(float).eps)
    ys = (y - np.mean(y)) / max(np.std(y), np.finfo(float).eps)

    if mode in {"constant", "offset"}:
        basis = np.ones((coords_m.shape[0], 1))
    elif mode in {"x_linear", "linear_x"}:
        basis = np.column_stack((np.ones_like(xs), xs))
    elif mode in {"affine", "linear_xy"}:
        basis = np.column_stack((np.ones_like(xs), xs, ys))
    elif mode in {"quadratic", "quadratic_xy"}:
        basis = np.column_stack(
            (np.ones_like(xs), xs, ys, xs * xs, xs * ys, ys * ys)
        )
    else:
        raise ValueError(
            "bed_projection_basis must be none, constant, x_linear, affine, "
            "or quadratic"
        )

    coefficients, *_ = np.linalg.lstsq(basis, increment, rcond=None)
    return basis @ coefficients


def _bed_coordinates(model_kwargs, hdim):
    coordinates = model_kwargs.get("mesh_coordinates")
    if coordinates is None:
        getter = model_kwargs.get("mesh_coordinate_getter")
        if getter is None:
            getter = _registry_coordinates
        if callable(getter):
            coordinates = getter(model_kwargs)
    if coordinates is None:
        return None

    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.shape != (hdim, 2):
        raise ValueError(
            f"Expected mesh coordinates with shape {(hdim, 2)}, "
            f"got {coordinates.shape}"
        )
    scale = float(model_kwargs.get("mesh_coordinate_scale_to_m", 1.0))
    return coordinates * scale


def apply_bed_regularized_correction(
    analysis_vec,
    vec_inputs,
    hdim,
    model_kwargs,
    timestep=None,
    model_time=None,
):
    """Stabilize member-wise bed updates produced by the analysis.

    Configuration keys
    ------------------
    physics_bed_inference          bool, default False
    bed_inference_start_time       float, default 0
    bed_spinup_hold_factor         float in [0, 1], default 1
    bed_blend_ramp_time            float, default 0
    bed_update_blend_factor        float in [0, 1], default 0.15
    bed_spatial_regularization     float, default 40
    bed_graph_neighbors            int, default 12
    bed_max_update_per_cycle       optional positive metres
    bed_projection_basis           none/constant/x_linear/affine/quadratic
    bed_physical_bounds            optional (minimum, maximum) metres
    bed_enforce_below_surface      bool, default True when surface is present
    bed_min_surface_separation     float, default 1 metre
    bed_update_mask                optional (hdim,) weights in [0, 1]
    """
    if not model_kwargs.get("physics_bed_inference", False):
        return analysis_vec

    bed_slice = _bed_state_slice(
        vec_inputs,
        # Do not alias "base" to bedrock: on floating ice, the ice base and
        # the subglacial bed are distinct physical surfaces.
        {"bed", "bedrock", "bedtopography", "bed_elevation"},
        hdim,
    )
    if bed_slice is None:
        return analysis_vec

    surface_slice = _bed_state_slice(
        vec_inputs,
        {"surface", "ice_surface", "s", "surface_elevation"},
        hdim,
    )
    coords_m = _bed_coordinates(model_kwargs, hdim)
    if coords_m is None:
        return analysis_vec

    bed_analysis = np.asarray(analysis_vec[bed_slice, :], dtype=float)

    initial_reference = model_kwargs.get("_bed_initial_reference")
    if initial_reference is None or np.shape(initial_reference) != np.shape(bed_analysis):
        initial_reference = bed_analysis.copy()
        model_kwargs["_bed_initial_reference"] = initial_reference

    previous = model_kwargs.get("_bed_previous_applied")
    if previous is None or np.shape(previous) != np.shape(bed_analysis):
        forecast_reference = model_kwargs.get("_bed_forecast_reference")
        if np.shape(forecast_reference) == np.shape(bed_analysis):
            previous = np.asarray(forecast_reference, dtype=float).copy()
        else:
            previous = initial_reference.copy()

    params = model_kwargs.get("params", {})
    dt = float(model_kwargs.get("dt", params.get("dt", 1.0)))
    if model_time is not None:
        time_now = float(model_time)
    elif timestep is not None:
        time_now = float(timestep) * dt
    else:
        call = int(model_kwargs.get("_bed_inference_call_count", 0))
        time_now = call * dt
        model_kwargs["_bed_inference_call_count"] = call + 1

    start_time = float(model_kwargs.get("bed_inference_start_time", 0.0))
    if time_now <= start_time:
        hold = np.clip(
            float(model_kwargs.get("bed_spinup_hold_factor", 1.0)), 0.0, 1.0
        )
        corrected = (1.0 - hold) * bed_analysis + hold * initial_reference
    else:
        raw_increment = bed_analysis - previous

        # An optional sensitivity/localization mask can completely suppress bed
        # updates where H/u/v provide little information.
        update_mask = model_kwargs.get("bed_update_mask")
        if update_mask is not None:
            update_mask = np.asarray(update_mask, dtype=float).reshape(hdim, 1)
            raw_increment *= np.clip(update_mask, 0.0, 1.0)

        spatial_strength = max(
            float(model_kwargs.get("bed_spatial_regularization", 40.0)), 0.0
        )
        graph_neighbors = int(model_kwargs.get("bed_graph_neighbors", 12))
        cache_key = (
            hdim,
            float(np.sum(coords_m)),
            float(np.sum(coords_m * coords_m)),
            spatial_strength,
            graph_neighbors,
        )
        cache = model_kwargs.get("_bed_regularization_cache")
        if cache is None or cache[0] != cache_key:
            solver = _bed_increment_solver(
                coords_m, spatial_strength, graph_neighbors
            )
            model_kwargs["_bed_regularization_cache"] = (cache_key, solver)
        else:
            solver = cache[1]

        smooth_increment = np.column_stack(
            [solver(raw_increment[:, member]) for member in range(raw_increment.shape[1])]
        )
        smooth_increment = _project_bed_increment(
            smooth_increment,
            coords_m,
            model_kwargs.get("bed_projection_basis", "none"),
        )

        maximum_update = model_kwargs.get("bed_max_update_per_cycle")
        if maximum_update is not None:
            maximum_update = float(maximum_update)
            if maximum_update <= 0.0:
                raise ValueError("bed_max_update_per_cycle must be positive")
            smooth_increment = np.clip(
                smooth_increment, -maximum_update, maximum_update
            )

        blend = np.clip(
            float(model_kwargs.get("bed_update_blend_factor", 0.15)), 0.0, 1.0
        )
        ramp_time = max(float(model_kwargs.get("bed_blend_ramp_time", 0.0)), 0.0)
        if ramp_time > 0.0:
            blend *= np.clip((time_now - start_time) / ramp_time, 0.0, 1.0)
        corrected = previous + blend * smooth_increment

    bounds = model_kwargs.get("bed_physical_bounds")
    if bounds is not None:
        lower, upper = map(float, bounds)
        if lower >= upper:
            raise ValueError("bed_physical_bounds must satisfy lower < upper")
        corrected = np.clip(corrected, lower, upper)

    if (
        surface_slice is not None
        and model_kwargs.get("bed_enforce_below_surface", True)
    ):
        surface = np.asarray(analysis_vec[surface_slice, :], dtype=float)
        separation = max(
            float(model_kwargs.get("bed_min_surface_separation", 1.0)), 0.0
        )
        corrected = np.minimum(corrected, surface - separation)

    analysis_vec[bed_slice, :] = corrected
    model_kwargs["_bed_previous_applied"] = corrected.copy()
    return analysis_vec
