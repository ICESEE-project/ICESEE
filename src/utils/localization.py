"""Model-agnostic coordinate registry and patch-based local EnKF analysis."""

import numpy as np
from scipy.spatial import cKDTree


_COORD_PROVIDERS = {}


def register_coord_provider(model_name, fn):
    """Register a model-specific callable returning node coordinates."""
    _COORD_PROVIDERS[str(model_name).lower()] = fn


def get_mesh_coordinates(model_kwargs):
    """Return and cache ``(n_dof, 2)`` coordinates for the current model."""
    if "mesh_coords" in model_kwargs:
        return model_kwargs["mesh_coords"]

    model_name = str(model_kwargs.get("model_name", "")).lower()
    provider = _COORD_PROVIDERS.get(model_name)
    coords = None

    if provider is None:
        print(
            f"[ICESEE][localization] No coordinate provider for model "
            f"'{model_name}'; localization disabled for this run."
        )
    else:
        try:
            coords = provider(model_kwargs)
        except Exception as exc:
            print(
                f"[ICESEE][localization] Coordinate provider for "
                f"'{model_name}' failed: {exc}"
            )

    model_kwargs["mesh_coords"] = coords
    return coords


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
