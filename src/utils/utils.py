# =============================================================================
# @Author: Brian Kyanjo
# @Date: 2024-09-24
# @Description: This script includes the some of the utility functions used in the
#               EnKF data assimilation scheme. 
# =============================================================================

# import libraries
import h5py
import numpy as np
import re
import sys
import traceback
from collections.abc import Iterable
from scipy.stats import norm
from scipy.interpolate import interp1d
from scipy.spatial.distance import cdist


# import utility functions
from ICESEE.src.utils.tools import icesee_get_index


# --- helper functions ---
def isiterable(obj):
    return isinstance(obj, Iterable)

class UtilsFunctions:
    def __init__(self, params=None, model_kwargs=None, ensemble=None):
        self.params = params
        self.ensemble = ensemble
        self.model_kwargs = model_kwargs 

    # ------------------------------ helpers ------------------------------
    @staticmethod
    def _as_list(x):
        import numpy as np
        if x is None:
            return []
        if isinstance(x, (list, tuple)):
            return list(x)
        if isinstance(x, np.ndarray):
            return list(x.ravel())
        return [x]

    def _get_obs_col_from_model_step(self, k_model: int):
        """
        Map a model time index k_model (0..nt-1) -> obs column km (0..m_obs-1),
        using mapping produced by _create_synthetic_observations.
        """
        m = self.model_kwargs.get("obs_model_to_col", None)
        if m is None:
            return None
        return m.get(int(k_model), None)

    # ------------------------------ H -----------------------------------
    def H_matrix(self, n_model, km=None,obs_mask_full=None):
        """
        Dense observation operator H.

        If km is not None (obs column index), bed rows are included ONLY when:
          - km is a bed snapshot column
          - and the bed mask for that column is True (tracks ∩ grounded)

        Non-bed observed variables are included normally (they are not time-gated here).
        """

        import numpy as np

        params = self.params
        observed = params["all_observed"]                 # e.g., ['h','u','v','smb','bed']
        vec_inputs = self.model_kwargs["vec_inputs"]      # e.g., ['h','s','u','v','bed','fric','smb']

        vecs, indx_map, _ = icesee_get_index(**self.model_kwargs)

        bed_aliases = {'bed', 'bedrock', 'bed_topography', 'bedtopo', 'bedtopography'}
        is_bed_key = lambda k: (str(k).lower() in bed_aliases)

        bed_mask_static = self.model_kwargs.get("bed_mask_map_static", {})  # key -> (n_bed,)
        bed_mask_cols   = self.model_kwargs.get("bed_mask_map_cols", {})    # key -> (n_bed, m_obs)
        bed_snap_cols   = set(self.model_kwargs.get("bed_snap_cols", []))

        all_obs_indices = []

        for key in observed:
            idx = np.asarray(indx_map[key], dtype=int)

            if is_bed_key(key):
                # If we are building H for a specific obs column, include bed only if that column is a snapshot.
                if km is not None:
                    if int(km) not in bed_snap_cols:
                        # no bed observations at this obs time
                        continue

                    if key in bed_mask_cols:
                        mask = np.asarray(bed_mask_cols[key][:, int(km)], dtype=bool)  # (n_bed,)
                    elif key in bed_mask_static:
                        mask = np.asarray(bed_mask_static[key], dtype=bool)
                    else:
                        continue

                    if mask.size != idx.size:
                        raise ValueError(
                            f"[H_matrix] bed mask size mismatch for '{key}': "
                            f"mask={mask.size}, idx={idx.size}"
                        )
                    idx = idx[mask]

                else:
                    # legacy time-independent H (not recommended for grounded-only snapshots)
                    if key in bed_mask_static:
                        mask = np.asarray(bed_mask_static[key], dtype=bool)
                        if mask.size != idx.size:
                            raise ValueError(
                                f"[H_matrix] bed static mask size mismatch for '{key}': "
                                f"mask={mask.size}, idx={idx.size}"
                            )
                        idx = idx[mask]
                    # else: include all bed points

            # NEW: apply missing-data mask (NaN-driven)
            if obs_mask_full is not None:
                obs_mask_full = np.asarray(obs_mask_full, dtype=bool).ravel()
                idx = idx[obs_mask_full[idx]]

            if idx.size > 0:
                all_obs_indices.append(idx)   

            # all_obs_indices.append(idx)

        if len(all_obs_indices) == 0:
            return np.zeros((0, n_model), dtype=float)

        obs_indices = np.concatenate(all_obs_indices).astype(int)
        if obs_indices.size == 0:
            return np.zeros((0, n_model), dtype=float)

        if obs_indices.max() >= n_model:
            raise ValueError(
                f"[H_matrix] obs index {obs_indices.max()} >= state size {n_model}"
            )

        m_obs = obs_indices.size
        H = np.zeros((m_obs, n_model), dtype=float)
        H[np.arange(m_obs), obs_indices] = 1.0
        return H

    def Obs_fun(self, virtual_obs, H=None, km=None):
        n = 1 if np.isscalar(virtual_obs) else virtual_obs.shape[0]
        if H is None:
            hu_obs = self.model_kwargs.get("hu_obs_loaded", None)
            #  read hu_obs from file
            # _synthetic_obs = self.model_kwargs.get("synthetic_obs_file", None)
            # with h5py.File(_synthetic_obs, 'r') as f:
            #     hu_obs  = f['hu_obs'][:]
            obs_mask_full = None
            if (hu_obs is not None) and (km is not None):
                obs_mask_full = ~np.isnan(hu_obs[:, int(km)])
            H = self.H_matrix(n, km=km, obs_mask_full=obs_mask_full)
        return H @ virtual_obs

    def JObs_fun(self, n_model):
        k_model = self.model_kwargs.get("k", None)
        km = None
        if k_model is not None:
            km = self._get_obs_col_from_model_step(k_model)

        hu_obs = self.model_kwargs.get("hu_obs_loaded", None)
        #  read hu_obs from file
        # _synthetic_obs = self.model_kwargs.get("synthetic_obs_file", None)
        # with h5py.File(_synthetic_obs, 'r') as f:
        #     hu_obs  = f['hu_obs'][:]

        obs_mask_full = None
        if (hu_obs is not None) and (km is not None):
            obs_mask_full = ~np.isnan(hu_obs[:, int(km)])

        return self.H_matrix(n_model, km=km, obs_mask_full=obs_mask_full)

    # -------------------------- observation schedule --------------------------
    def generate_observation_schedule(self, **kwargs):
        import numpy as np
        import traceback, sys

        try:
            t = np.asarray(kwargs["t"], dtype=float)
            if t.ndim != 1 or t.size == 0:
                raise ValueError("`t` must be a 1D non-empty array of times.")
            t_min, t_max = float(t[0]), float(t[-1])

            freq_obs = float(self.params["freq_obs"])
            obs_start = float(self.params["obs_start_time"])
            obs_max_cfg = float(self.params["obs_max_time"])

            obs_start = max(obs_start, t_min)
            obs_max = min(obs_max_cfg, t_max)

            if freq_obs <= 0.0 or obs_start > obs_max:
                return np.array([]), np.array([], dtype=int), 0

            n_obs = int(np.floor((obs_max - obs_start) / freq_obs)) + 1
            obs_t_req = obs_start + np.arange(n_obs, dtype=float) * freq_obs

            dt_grid = np.min(np.diff(t)) if len(t) > 1 else 1.0

            obs_idx = []
            for tobs in obs_t_req:
                i = int(np.argmin(np.abs(t - tobs)))
                if abs(t[i] - tobs) <= 0.5 * dt_grid:
                    obs_idx.append(i)

            obs_idx = np.array(sorted(set(obs_idx)), dtype=int)
            return obs_t_req, obs_idx, len(obs_idx)

        except Exception as e:
            print(f"Error occurred in generate_observation_schedule: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
            return np.array([]), np.array([], dtype=int), 0

    # --------------------- consistent synthetic observation maker ---------------------
    def _create_synthetic_observations(self, **kwargs):
        """
        Create synthetic observations hu_obs consistent with H_matrix(..., km).

        Returns:
            hu_obs:      (nd, m_obs)  full-state-sized with zeros at unobserved entries
            error_R:     (m_obs*2+1, nd)  (kept same orientation as your return; adjust if needed)
            bed_masks:   dict with:
                          - 'static' : bed_mask_map_static[key] -> (n_bed,)
                          - 'cols'   : bed_mask_map_cols[key]   -> (n_bed, m_obs)
        Side-effects:
            Updates model_kwargs with:
              bed_mask_map_static, bed_mask_map_cols, bed_snap_cols,
              ind_m, obs_t, obs_model_to_col mapping.
        """
        import numpy as np

        statevec_true = kwargs.get("statevec_true", None)
        if statevec_true is None:
            raise ValueError("statevec_true is required")

        params = kwargs.get("params", self.params)
        vec_inputs = list(kwargs["vec_inputs"])

        # Observation schedule
        obs_t_req, ind_m, m_obs = self.generate_observation_schedule(**kwargs)
        ind_m = np.asarray(ind_m, dtype=int)  # model time indices (0-based)
        obs_t_req = np.asarray(obs_t_req, dtype=float)

        print(f"[ICESEE] observation times requested: {obs_t_req}")
        print(f"[ICESEE] observation model indices ind_m: {ind_m}, total m_obs={m_obs}")

        # map model time index -> obs column
        obs_model_to_col = {int(k_model): int(col) for col, k_model in enumerate(ind_m)}

        # Bed snapshot times (in same units as obs_t_req)
        bed_snaps = np.asarray(kwargs.get("bed_obs_snapshot", []), dtype=float).ravel()
        bed_snap_cols = []
        if bed_snaps.size > 0 and obs_t_req.size > 0:
            for bed_time in bed_snaps:
                j = int(np.argmin(np.abs(obs_t_req - bed_time)))
                bed_snap_cols.append(j)
        bed_snap_cols = sorted(set(bed_snap_cols))

        print("[ICESEE] bed_snaps:", bed_snaps)
        print("[ICESEE] bed_snap_cols:", bed_snap_cols)

        # Indices in state vector
        vecs, indx_map, _ = icesee_get_index(statevec_true, **kwargs)

        # preallocate obs in *state* size (so you can keep your old hu_obs layout)
        nd_full = statevec_true.shape[0]
        # hu_obs = np.zeros((nd_full, m_obs), dtype=float)
        hu_obs = np.full((nd_full, m_obs), np.nan, dtype=float)  # use NaN for unobserved

        # error_R: keep your block-structure logic
        total_state_param_vars = params["total_state_param_vars"]
        hdim = nd_full // total_state_param_vars

        # error_R = np.zeros((nd_full, m_obs), dtype=float)
        error_R = np.full((nd_full, m_obs), np.nan, dtype=float)  # use NaN for unobserved
        sig_obs = params["sig_obs"]
        # for i, sig in enumerate(sig_obs):
        #     a = i * hdim
        #     b = a + hdim
        #     error_R[a:b, :] = float(sig)

        bed_aliases = {'bed', 'bedrock', 'bed_topography', 'bedtopo', 'bedtopography'}
        key_is_bed = {k: (str(k).lower() in bed_aliases) for k in vec_inputs}

        # ---------------- build STATIC bed mask (tracks/spacing/user mask) ----------------
        Lx = kwargs.get("Lx", self.params.get("Lx", None))
        Ly = kwargs.get("Ly", self.params.get("Ly", None))
        model_name = kwargs.get("model_name", None)

        bed_stride_km = kwargs.get("bed_obs_stride", None)
        bed_spacing_pts = kwargs.get("bed_obs_spacing", None)
        bed_indices_user = kwargs.get("bed_obs_indices", None)
        bed_mask_user = kwargs.get("bed_obs_mask", None)

        bed_mask_map_static = {}
        bed_mask_map_cols = {}

        # Precompute key indices
        key_idx_map = {k: np.asarray(indx_map[k], dtype=int) for k in vec_inputs if k in indx_map}

        for k in vec_inputs:
            if not key_is_bed.get(k, False):
                continue

            bed_idx = key_idx_map[k]       # global indices for bed
            n_bed = bed_idx.size

            # default: observe everywhere (then time/grounded gate will reduce further)
            mask = np.ones(n_bed, dtype=bool)

            # Priority 1: explicit boolean mask
            if isinstance(bed_mask_user, (list, np.ndarray)):
                um = np.asarray(bed_mask_user, dtype=bool).ravel()
                if um.size == n_bed:
                    mask = um
                else:
                    # safe fallback
                    mask = np.ones(n_bed, dtype=bool)

            # Priority 2: explicit indices
            elif isinstance(bed_indices_user, (list, np.ndarray)):
                mask = np.zeros(n_bed, dtype=bool)
                idxs = np.asarray(bed_indices_user, dtype=int).ravel()
                idxs = idxs[(idxs >= 0) & (idxs < n_bed)]
                mask[idxs] = True

            # Priority 3: spacing in points
            elif isinstance(bed_spacing_pts, (int, np.integer)) and int(bed_spacing_pts) > 1:
                n = int(bed_spacing_pts)
                mask = np.zeros(n_bed, dtype=bool)
                mask[::n] = True

            # Priority 4: stride in km (ISSM mesh uses x/y files)
            elif (bed_stride_km is not None) and (Lx is not None) and (Ly is not None):
                import re, h5py

                if re.match(r"(?i)^issm$", str(model_name)):
                    icesee_path = kwargs.get("icesee_path")
                    data_path = kwargs.get("data_path")
                    mesh_file = f"{icesee_path}/{data_path}/mesh_idxy_0.h5"
                    with h5py.File(mesh_file, "r") as f:
                        x_param = np.asarray(f["/fric_x"][:], dtype=float).ravel() / 1000.0
                        y_param = np.asarray(f["/fric_y"][:], dtype=float).ravel() / 1000.0

                    # tracks perpendicular to flow along x (as in your kriging script)
                    x_min, x_max = x_param.min(), x_param.max()
                    stride = float(bed_stride_km)  # km
                    x_lines = np.arange(x_min, x_max + 1e-6, stride)

                    if x_lines.size > 1:
                        dx_nom = (x_max - x_min) / max(x_lines.size - 1, 1)
                    else:
                        dx_nom = stride
                    band = 0.5 * dx_nom

                    mask = np.zeros(n_bed, dtype=bool)
                    for xl in x_lines:
                        mask |= (np.abs(x_param - xl) <= band)

                    print(f"[ICESEE<-ISSM] bed static track mask '{k}': {mask.sum()} of {n_bed}")

            bed_mask_map_static[k] = mask
            bed_mask_map_cols[k] = np.zeros((n_bed, m_obs), dtype=bool)

        # ---------------- build observations per obs column ----------------
        obs_set = set(kwargs.get("observed_vars", []) + kwargs.get("observed_params", []))
        for ii, key in enumerate(vec_inputs):
            if key not in obs_set:
                continue
            if key_is_bed.get(key, False):
                continue  # bed handled later due to snapshots/masks

            error_R[indx_map[key], :] = float(sig_obs[ii])

        # locate thickness key once (needed for grounded-only bed obs)
        thickness_candidates = ["h", "thickness", "ice_thickness"]
        thickness_key = None
        for key in indx_map.keys():
            if key.lower() in thickness_candidates:
                thickness_key = key  
                break

        di = float(kwargs.get("di", 0.8930))

        for km, k_model in enumerate(ind_m):
            # Non-bed observations
            for ii, key in enumerate(vec_inputs):
                if key not in indx_map:
                    continue
                idx = np.asarray(indx_map[key], dtype=int)

                # make error_R
                # if (key in obs_set) and (not key_is_bed.get(key, False)):
                #     sigma = float(sig_obs[vec_inputs.index(key)])
                #     error_R[idx, km] = sigma

                if (key in obs_set) and (not key_is_bed.get(key, False)):
                    sigma = error_R[idx, km]
                    hu_obs[idx, km] = statevec_true[idx, k_model] + np.random.normal(0.0, sigma, size=idx.size)
                else:
                    #  leave NaN for unobserved entries
                    pass

            # Bed observations only at snapshot columns
            if km in bed_snap_cols:
                for ii, key in enumerate(vec_inputs):
                    if not key_is_bed.get(key, False):
                        continue
                    if key not in indx_map:
                        continue
                    if key not in obs_set:
                        continue

                    bed_idx = np.asarray(indx_map[key], dtype=int)  # global indices for bed
                    static_mask = bed_mask_map_static.get(key, None)
                    if static_mask is None:
                        static_mask = np.ones(bed_idx.size, dtype=bool)

                    # grounded-only: need thickness at same nodes
                    if thickness_key is None:
                        raise ValueError("Cannot find thickness variable (h/thickness/ice_thickness) required for grounded-only bed observations.")

                    h_idx = np.asarray(indx_map[thickness_key], dtype=int)
                    if h_idx.size != bed_idx.size:
                        raise ValueError(
                            f"Thickness vector length {h_idx.size} != bed vector length {bed_idx.size}. "
                            "Grounded-only bed obs assumes pointwise alignment."
                        )

                    bed_local = statevec_true[bed_idx, k_model]
                    h_local = statevec_true[h_idx, k_model]

                    # grounded mask (common choice): ocean_levelset = h + bed/di > 0
                    ocean_levelset = h_local + bed_local / di
                    grounded = (ocean_levelset > 0.0)

                    final_mask = static_mask & grounded
                    bed_mask_map_cols[key][:, km] = final_mask

                    idx_obs = bed_idx[final_mask]
                    if idx_obs.size > 0:
                        # sigma_obs = error_R[idx_obs, km]
                        sigma_bed = float(sig_obs[ii])
                        error_R[idx_obs, km] = sigma_bed
                        hu_obs[idx_obs, km] = statevec_true[idx_obs, k_model] + np.random.normal(
                            0.0, sigma_bed, size=idx_obs.size
                        )

        # ---------------- publish masks + mapping so H_matrix can match ----------------
        kwargs["obs_t"] = obs_t_req
        kwargs["ind_m"] = ind_m
        kwargs["m_obs"] = m_obs
        kwargs["obs_model_to_col"] = obs_model_to_col

        kwargs["bed_snap_cols"] = bed_snap_cols
        kwargs["bed_mask_map_static"] = bed_mask_map_static
        kwargs["bed_mask_map_cols"] = bed_mask_map_cols

        bed_masks = {"static": bed_mask_map_static, "cols": bed_mask_map_cols}
        return hu_obs, error_R.T, bed_masks, kwargs