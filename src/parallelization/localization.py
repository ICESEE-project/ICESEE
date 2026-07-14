"""
localization.py

Static, variable-agnostic parameter-space localization for ICESEE's
partial-parallel stochastic EnKF.

Design:
    - Increment tapering (cheap): analysis increment for any block named in
      params["localization"]["localized_vars"] is scaled by a static,
      distance-based taper rho in [0, 1] before being applied. State-variable
      blocks and the global X5 computation are completely untouched.
    - Taper is computed ONCE per run (static), from each node's distance to
      the nearest node that is EVER observed for that variable, across the
      whole simulation.
    - Fully automated: everything needed is derived from existing
      model_kwargs / params.yaml entries and cached; no manual wiring beyond
      exposing a mesh coordinate source per model (see get_mesh_coordinates).

params.yaml:
    localization:
      enabled: false
      localized_vars: []          # e.g. ["bed", "c"]
      radius: 50.0                # km; scalar or {var: radius}
      taper_type: "gaspari_cohn"  # or "gaussian"
"""

import numpy as np
from scipy.spatial import distance


# ============================================================
# Taper kernels
# ============================================================

def gaspari_cohn(dist, c):
    """Gaspari-Cohn 5th-order taper. dist, c in same units. -> 0 beyond 2c."""
    dist = np.asarray(dist, dtype=float)
    r = np.abs(dist) / max(c, 1e-12)
    taper = np.zeros_like(r)

    m1 = r <= 1.0
    r1 = r[m1]
    taper[m1] = 1 - (5/3)*r1**2 + (5/8)*r1**3 + 0.5*r1**4 - 0.25*r1**5

    m2 = (r > 1.0) & (r <= 2.0)
    r2 = r[m2]
    with np.errstate(divide="ignore", invalid="ignore"):
        taper[m2] = (4 - 5*r2 + (5/3)*r2**2 + (5/8)*r2**3
                     - 0.5*r2**4 + (1/12)*r2**5 - 2/(3*r2))

    return np.clip(np.nan_to_num(taper), 0.0, 1.0)


def gaussian_taper(dist, c):
    dist = np.asarray(dist, dtype=float)
    return np.exp(-0.5 * (dist / max(c, 1e-12)) ** 2)


_TAPER_FUNCS = {
    "gaspari_cohn": gaspari_cohn,
    "gaussian": gaussian_taper,
}


# ============================================================
# Mesh coordinate source (per-model, cached, soft-fails)
# ============================================================

def get_mesh_coordinates(model_kwargs):
    """
    Return (n,2) array of x,y node coordinates [km] for the current model's
    mesh, or None if unavailable. Cached in model_kwargs['mesh_coords'].
    """
    if "mesh_coords" in model_kwargs:
        return model_kwargs["mesh_coords"]

    model_name = str(model_kwargs.get("model_name", "")).lower()
    coords = None

    if model_name == "issm":
        icesee_path = model_kwargs.get("icesee_path")
        data_path = model_kwargs.get("data_path")
        mesh_file = f"{icesee_path}/{data_path}/mesh_idxy_0.h5"
        try:
            import h5py
            with h5py.File(mesh_file, "r") as f:
                x = np.asarray(f["/fric_x"][:], dtype=float).ravel() / 1000.0
                y = np.asarray(f["/fric_y"][:], dtype=float).ravel() / 1000.0
            coords = np.column_stack([x, y])
        except (FileNotFoundError, KeyError, OSError) as e:
            print(f"[ICESEE][localization] ISSM mesh coords unavailable: {e}")

    elif model_name == "icepack":
        mesh = model_kwargs.get("mesh")
        if mesh is not None:
            try:
                coords = np.asarray(mesh.coordinates.dat.data_ro, dtype=float)
            except AttributeError as e:
                print(f"[ICESEE][localization] Could not read Icepack mesh coords: {e}")
        else:
            print("[ICESEE][localization] Icepack mesh not found in model_kwargs "
                  "(expected model_kwargs['mesh']); localization disabled.")

    else:
        print(f"[ICESEE][localization] No coordinate source for model "
              f"'{model_name}'; localization disabled.")

    model_kwargs["mesh_coords"] = coords
    return coords


# ============================================================
# Observation-anchor node derivation (automated from existing masks)
# ============================================================

def _observed_node_mask(key, n_nodes, model_kwargs):
    """
    Boolean mask (n_nodes,) of nodes ever used as an observation anchor
    for `key`, across the whole run. Derived entirely from structures
    already built in _create_synthetic_observations.
    """
    bed_mask_cols = model_kwargs.get("bed_mask_map_cols", {})
    bed_mask_static = model_kwargs.get("bed_mask_map_static", {})
    obs_set = set(model_kwargs.get("observed_vars", []) +
                  model_kwargs.get("observed_params", []))

    if key in bed_mask_cols:
        return np.any(bed_mask_cols[key], axis=1)
    if key in bed_mask_static:
        return np.asarray(bed_mask_static[key], dtype=bool)
    if key in obs_set:
        return np.ones(n_nodes, dtype=bool)
    return np.zeros(n_nodes, dtype=bool)


# ============================================================
# Taper builder (static, cached, fully automated)
# ============================================================

def build_localization_tapers(vec_inputs, indx_map, model_kwargs):
    """
    Compute and cache a static taper array per localized variable.
    No-op (returns {}) if localization is disabled or unconfigured.
    Safe to call every cycle — only computes once.
    """
    if "localization_taper" in model_kwargs:
        return model_kwargs["localization_taper"]

    loc_cfg = model_kwargs.get("params", {}).get("localization", {}) \
        if "params" in model_kwargs else model_kwargs.get("localization", {})

    if not loc_cfg.get("enabled", False):
        model_kwargs["localization_taper"] = {}
        return {}

    localized_vars = [v for v in loc_cfg.get("localized_vars", []) if v in vec_inputs]
    if not localized_vars:
        model_kwargs["localization_taper"] = {}
        return {}

    coords = get_mesh_coordinates(model_kwargs)
    if coords is None:
        print("[ICESEE][localization] No mesh coordinates; localization disabled for this run.")
        model_kwargs["localization_taper"] = {}
        return {}

    radius_cfg = loc_cfg.get("radius", 50.0)
    taper_type = loc_cfg.get("taper_type", "gaspari_cohn")
    taper_fn = _TAPER_FUNCS.get(taper_type, gaspari_cohn)

    tapers = {}
    n_nodes = coords.shape[0]

    for key in localized_vars:
        if key not in indx_map:
            print(f"[ICESEE][localization] '{key}' not in indx_map, skipping.")
            continue

        obs_mask = _observed_node_mask(key, n_nodes, model_kwargs)
        radius = radius_cfg[key] if isinstance(radius_cfg, dict) else radius_cfg

        if not np.any(obs_mask):
            print(f"[ICESEE][localization] '{key}': no obs anchors found; taper=0 everywhere.")
            tapers[key] = np.zeros(n_nodes)
            continue

        dmin = distance.cdist(coords, coords[obs_mask]).min(axis=1)
        tapers[key] = taper_fn(dmin, radius)

        print(f"[ICESEE][localization] '{key}': radius={radius} km, "
              f"{obs_mask.sum()}/{n_nodes} anchor nodes, "
              f"mean_rho={tapers[key].mean():.3f}")

    model_kwargs["localization_taper"] = tapers
    return tapers


# ============================================================
# Application (called inside analysis_enkf_update)
# ============================================================

def apply_localization(analysis_vec, prior_vec, global_rows, vec_inputs, hdim, model_kwargs):
    """
    Taper the analysis increment in-place for any block with a cached taper.
    Fully automated no-op if localization is disabled/unconfigured.
    """
    tapers = model_kwargs.get("localization_taper", {})
    if not tapers:
        return analysis_vec

    for ii, key in enumerate(vec_inputs):
        if key not in tapers:
            continue

        start, end = ii * hdim, (ii + 1) * hdim
        local_mask = (global_rows >= start) & (global_rows < end)
        if not np.any(local_mask):
            continue

        local_node_idx = global_rows[local_mask] - start
        rho_local = tapers[key][local_node_idx].reshape(-1, 1)

        increment = analysis_vec[local_mask, :] - prior_vec[local_mask, :]
        analysis_vec[local_mask, :] = prior_vec[local_mask, :] + rho_local * increment

    return analysis_vec