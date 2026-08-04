"""Model-agnostic coordinate registry and patch-based local EnKF analysis."""

import numpy as np
from scipy.spatial import cKDTree


_COORD_PROVIDERS = {}


def register_coord_provider(model_name, fn):
    """Register a model-specific callable returning node coordinates."""
    _COORD_PROVIDERS[str(model_name).lower()] = fn


def get_mesh_coordinates(model_kwargs):
    """Return and cache physical coordinates for the current model nodes.

    Providers may return one- or multi-dimensional coordinates, but their
    leading dimension must follow the model state-vector node ordering.  A
    failed/early lookup is deliberately not cached so that file-backed model
    providers (for example ISSM) can be retried after model initialization.
    """
    if model_kwargs.get("mesh_coords") is not None:
        return model_kwargs["mesh_coords"]

    model_name = str(model_kwargs.get("model_name", "")).lower()
    provider = _COORD_PROVIDERS.get(model_name)
    coords = None

    if provider is None:
        print(
            f"[ICESEE][coordinates] No coordinate provider for model "
            f"'{model_name}'; coordinate-dependent features are unavailable."
        )
    else:
        try:
            coords = provider(model_kwargs)
        except Exception as exc:
            print(
                f"[ICESEE][coordinates] Coordinate provider for "
                f"'{model_name}' failed: {exc}"
            )

    if coords is None:
        return None

    coords = np.asarray(coords, dtype=float)
    if coords.ndim == 1:
        coords = coords[:, None]
    if coords.ndim != 2 or coords.shape[1] < 1:
        raise ValueError(
            "Mesh coordinates must have shape (n_nodes,) or "
            f"(n_nodes, n_spatial_dims); got {coords.shape}"
        )
    if not np.all(np.isfinite(coords)):
        raise ValueError("Mesh coordinates contain non-finite values")

    model_kwargs["mesh_coords"] = coords
    return coords


def prepare_random_field_coordinates(model_kwargs, expected_nodes=None):
    """Pack registered mesh coordinates for graph random-field generation.

    This is called after model initialization and before ensemble
    initialization.  FFT fields require no coordinates and therefore retain
    their historical behavior.  In graph mode, silently falling back to an
    index-space chain would defeat the purpose of selecting the graph method,
    so a missing provider or incompatible node count is reported immediately.
    """
    method = str(
        model_kwargs.get(
            "random_field_method",
            model_kwargs.get("enkf_field_method", "fft"),
        )
    ).strip().lower()
    if method not in {"fft", "graph"}:
        raise ValueError(
            f"Unsupported random_field_method={method!r}; expected 'fft' or 'graph'"
        )

    # Store the normalized spelling so every downstream noise-generation path
    # (initial ensemble, forecast noise, and stochastic observations) agrees.
    model_kwargs["random_field_method"] = method
    if method == "fft":
        return None

    coords = get_mesh_coordinates(model_kwargs)
    if coords is None:
        model_name = model_kwargs.get("model_name", "model")
        raise ValueError(
            f"random_field_method='graph' requires registered node coordinates "
            f"for model {model_name!r}"
        )
    if expected_nodes is not None and coords.shape[0] != int(expected_nodes):
        raise ValueError(
            "Registered mesh-coordinate count does not match the spatial state "
            f"dimension: {coords.shape[0]} != {int(expected_nodes)}"
        )

    model_kwargs["mesh_coords"] = coords
    if model_kwargs.get("verbose", False):
        print(
            f"[ICESEE] graph random fields use {coords.shape[0]} registered "
            f"mesh coordinates ({coords.shape[1]}D)"
        )
    return coords


def restore_frozen_analysis_vars(
    analysis_vec,
    forecast_vec,
    global_rows,
    vec_inputs,
    hdim,
    frozen_vars,
):
    """Restore selected state-vector blocks to their forecast values.

    This provides a clean fixed-parameter control while retaining assimilation
    updates for the remaining state and parameter blocks.
    """
    if not frozen_vars:
        return analysis_vec

    frozen_names = {str(value).lower() for value in frozen_vars}
    global_rows = np.asarray(global_rows, dtype=int)

    for block_index, key in enumerate(vec_inputs):
        if str(key).lower() not in frozen_names:
            continue
        start = block_index * hdim
        stop = start + hdim
        mask = (global_rows >= start) & (global_rows < stop)
        if np.any(mask):
            analysis_vec[mask, :] = forecast_vec[mask, :]

    return analysis_vec


def active_observation_std(model_kwargs, k_obs, obs_indices):
    """Return configured standard deviations for active observation rows."""
    error_r = model_kwargs.get("error_R")
    if error_r is None:
        raise ValueError(
            "enkf_observation_error_mode='stochastic_R' requires error_R"
        )

    error_r = np.asarray(error_r, dtype=float)
    obs_indices = np.asarray(obs_indices, dtype=int)
    k_obs = int(k_obs)
    max_index = int(obs_indices.max()) if obs_indices.size else -1

    if error_r.ndim == 1:
        sigma = error_r[obs_indices]
    elif (
        error_r.ndim == 2
        and error_r.shape[0] > k_obs
        and error_r.shape[1] > max_index
    ):
        sigma = error_r[k_obs, obs_indices]
    elif (
        error_r.ndim == 2
        and error_r.shape[1] > k_obs
        and error_r.shape[0] > max_index
    ):
        sigma = error_r[obs_indices, k_obs]
    else:
        raise ValueError(
            f"error_R shape {error_r.shape} is incompatible with observation "
            f"column {k_obs} and {obs_indices.size} active rows"
        )

    sigma = np.asarray(sigma, dtype=float).ravel()
    invalid = ~np.isfinite(sigma) | (sigma <= 0.0)
    if np.any(invalid):
        bad_rows = obs_indices[invalid][:8].tolist()
        raise ValueError(
            "Active observations require finite positive standard deviations; "
            f"invalid rows include {bad_rows}"
        )
    return sigma


def stochastic_observation_terms(HA, d, sigma, seed):
    """Build consistent perturbed-observation terms from configured ``R``."""
    HA = np.asarray(HA, dtype=float)
    d = np.asarray(d, dtype=float).ravel()
    sigma = np.asarray(sigma, dtype=float).ravel()
    if HA.shape[0] != d.size or d.size != sigma.size:
        raise ValueError(
            f"Observation shapes disagree: HA={HA.shape}, d={d.shape}, "
            f"sigma={sigma.shape}"
        )

    rng = np.random.default_rng(int(seed))
    eta = rng.standard_normal(HA.shape) * sigma[:, None]
    eta -= np.mean(eta, axis=1, keepdims=True)
    ha_prime = HA - np.mean(HA, axis=1, keepdims=True)
    d_prime = d[:, None] + eta - HA
    return ha_prime, eta, d_prime


def build_obs_coords(obs_indices, node_coords, vec_inputs, hdim):
    """Map global observation rows to per-block physical coordinates."""
    del vec_inputs  # retained in the public signature for compatibility
    obs_indices = np.asarray(obs_indices, dtype=int)
    if obs_indices.size == 0:
        return np.zeros((0, 2))
    return node_coords[obs_indices % hdim]


def compute_X5_from_matrices(HAprime, Eta, Dprime, Nens):
    """Compute the Evensen (2003) ensemble-space analysis transform."""
    m = Dprime.shape[0]
    nrmin = min(m, Nens)

    h_aprime_eta = HAprime + Eta
    U, singular_values, _ = np.linalg.svd(h_aprime_eta, full_matrices=False)
    eigenvalues = singular_values**2

    retained_sum = np.sum(eigenvalues[:nrmin])
    cumulative = 0.0
    for index in range(nrmin):
        if retained_sum > 0.0 and cumulative / retained_sum < 0.999:
            cumulative += eigenvalues[index]
            eigenvalues[index] = 1.0 / eigenvalues[index]
        else:
            eigenvalues[index:nrmin] = 0.0
            break

    x1 = eigenvalues[:nrmin, None] * U[:, :nrmin].T
    x2 = x1 @ Dprime
    x3 = U[:, :nrmin] @ x2
    x4 = HAprime.T @ x3
    return x4 + np.eye(Nens)


def estimate_adaptive_radius(obs_coords, target_count=None):
    """Estimate a localization radius from observation point density."""
    obs_coords = np.asarray(obs_coords, dtype=float)
    n_obs = obs_coords.shape[0]
    if n_obs < 2:
        return 1.0

    if target_count is None:
        target_count = min(max(20, int(np.sqrt(n_obs))), n_obs)

    minimum = obs_coords.min(axis=0)
    maximum = obs_coords.max(axis=0)
    area = max(np.prod(maximum - minimum), 1.0e-12)
    density = n_obs / area
    radius = np.sqrt(target_count / (np.pi * density))

    if not np.isfinite(radius) or radius <= 0.0:
        tree = cKDTree(obs_coords)
        k = min(target_count, n_obs)
        distances, _ = tree.query(obs_coords, k=k)
        radius = float(np.mean(distances[:, -1]))

    return radius


def compute_local_patches_X5(
    vec_inputs,
    hdim,
    HAprime,
    Eta,
    Dprime,
    Nens,
    obs_indices,
    model_kwargs,
):
    """Compute exact-grouped local transforms for patch-based analysis."""
    if not model_kwargs.get("local_analysis", False):
        return {}

    requested = model_kwargs.get("localized_vars", [])
    target_vars = (
        [value for value in requested if value in vec_inputs]
        if requested
        else list(vec_inputs)
    )

    node_coords = get_mesh_coordinates(model_kwargs)
    if node_coords is None:
        return {}

    obs_coords = build_obs_coords(obs_indices, node_coords, vec_inputs, hdim)
    if obs_coords.shape[0] == 0:
        _log_once(
            model_kwargs,
            "no_obs",
            "[ICESEE][local_analysis] No active observations this cycle; "
            "skipping.",
        )
        return {}

    manual_radius = model_kwargs.get("localization_radius")
    obs_tree = cKDTree(obs_coords)
    results = {}

    for key in target_vars:
        block_index = vec_inputs.index(key)
        start = block_index * hdim
        global_indices = start + np.arange(hdim)

        if manual_radius is not None:
            radius = (
                manual_radius[key]
                if isinstance(manual_radius, dict)
                else manual_radius
            )
            mode_string = f"radius={radius} (manual)"
        else:
            radius = estimate_adaptive_radius(
                obs_coords,
                target_count=model_kwargs.get("target_local_obs_count"),
            )
            mode_string = f"radius={radius:.2f} (auto, density-based)"

        neighbor_lists = obs_tree.query_ball_point(node_coords, r=radius)
        groups = {}
        for node_index, observation_neighbors in enumerate(neighbor_lists):
            if not observation_neighbors:
                continue
            signature = tuple(sorted(observation_neighbors))
            groups.setdefault(signature, []).append(node_index)

        results[key] = []
        for observation_tuple, node_group in groups.items():
            observation_rows = np.asarray(observation_tuple, dtype=int)
            local_transform = compute_X5_from_matrices(
                HAprime[observation_rows, :],
                Eta[observation_rows, :],
                Dprime[observation_rows, :],
                Nens,
            )
            node_group = np.asarray(node_group, dtype=int)
            results[key].append(
                (global_indices[node_group], local_transform)
            )

        group_count = len(results[key])
        _log_once(
            model_kwargs,
            f"groups_{key}",
            f"[ICESEE][local_analysis] '{key}': {group_count} unique "
            f"local-observation groups ({mode_string})",
            value=group_count,
        )

    return results


def _log_once(model_kwargs, tag, message, value=None):
    """Write a diagnostic only on first use or a meaningful value change."""
    from tqdm import tqdm

    cache = model_kwargs.setdefault("_log_once_cache", {})
    previous = cache.get(tag)
    should_print = previous is None
    if not should_print and value is not None:
        denominator = max(abs(previous), 1)
        should_print = abs(value - previous) / denominator > 0.10

    if should_print:
        communicator = model_kwargs.get("comm_world")
        rank = communicator.Get_rank() if communicator is not None else 0
        if rank == 0:
            tqdm.write(message)
        cache[tag] = value if value is not None else True


def apply_local_patches(analysis_vec, prior_vec, global_rows, local_patches):
    """Apply local transforms to distributed rows, in place."""
    if not local_patches:
        return analysis_vec

    for patches in local_patches.values():
        for patch_rows_global, local_transform in patches:
            mask = np.isin(global_rows, patch_rows_global)
            if np.any(mask):
                analysis_vec[mask, :] = prior_vec[mask, :] @ local_transform

    return analysis_vec
