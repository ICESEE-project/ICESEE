# =============================================================================
# @Author: Brian Kyanjo
# @Date: 2024-09-24
# @Description: This script includes the some of the utility functions used in the
#               EnKF data assimilation scheme. 
# =============================================================================

# import libraries
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
    def __init__(self, params=None,model_kwargs=None, ensemble=None):
        """
        Initialize the utility functions with model parameters.
        
        Parameters:
        params (dict): Model parameters, including those used for bed topography.
        """
        self.params   = params
        self.ensemble = ensemble
        self.model_kwargs = model_kwargs

    def H_matrix(self, n_model):
        """
        Build dense observation operator H.
        H has shape (m_obs, n_model).
        H(i,j) = 1 if observation i corresponds to state index j, else 0.
        """

        params = self.params
        observed = params["all_observed"]         # ['h','u','v','smb']
        vec_inputs = params["vec_inputs"]         # ['h','s','u','v','bed','fric','smb']

        # Recompute index map for the *current full state vector*
        vecs, indx_map, _ = icesee_get_index(**self.model_kwargs)

        # COLLECT OBSERVATION INDICES
        obs_indices = np.concatenate([indx_map[key] for key in observed]).astype(int)

        # SAFETY CHECK
        if obs_indices.max() >= n_model:
            raise ValueError(
                f"H_matrix error: obs index {obs_indices.max()} >= state size {n_model}. "
                "likely vec_inputs or nd inconsistent."
            )

        m_obs = obs_indices.size

        # H: zero everywhere except H[i, obs_indices[i]] = 1
        H = np.zeros((m_obs, n_model))
        H[np.arange(m_obs), obs_indices] = 1.0

        return H

    def H_matrix__(self, n_model):
        """
        Build observation operator H for multi-variable ice-sheet DA
        WITHOUT changing the outer pipeline.
        H returns a dense matrix of shape (m_obs, n_model).
        """

        # unpack parameters
        params = self.params
        observed = params['all_observed']       # e.g., ['h','u','v','smb']
        vec_inputs = params['vec_inputs']       # full list like ['h','s','u','v','bed','fric','smb']

        # Full state indexing
        vecs, indx_map, dim_per_proc = icesee_get_index(**self.model_kwargs)

        # Build consistent obs index list
        obs_indices = []
        for key in observed:
            obs_indices.append(indx_map[key])
        obs_indices = np.concatenate(obs_indices)

        # Number of observations
        m_obs = obs_indices.size

        # Allocate H
        H = np.zeros((m_obs, n_model))

        # Fill H with identity rows at observation positions
        H[np.arange(m_obs), obs_indices] = 1.0
        print(f"observed indices: {obs_indices}")

        return H

    def H_matrix_(self, n_model):
        """ observation operator matrix
        """
        n = n_model

        # Initialize the H matrix
        H = np.zeros((self.params["number_obs_instants"] * 2 + 1, n))

        # Calculate distance between measurements
        di = int((n - 2) / (2 * self.params["number_obs_instants"]))

        # Fill the H matrix
        for i in range(1,self.params["number_obs_instants"]+1):
            H[i-1, i * di-1] = 1
            H[self.params["number_obs_instants"] + i - 1, int((n - 2) / 2) + i * di -1] = 1

        H[self.params["number_obs_instants"] * 2, n - 2] = 1  # Final element

        # check if we have parameter estimation
        if self.params.get('joint_estimation', False):
            ndim = n // self.params["total_state_param_vars"]
            state_variables_size = ndim*self.params["num_state_vars"]
            # parameters are not required to observe the state variables
            num_params_size = n - state_variables_size
            # H_param = np.zeros(num_params_size)
            # H[:,state_variables_size:] = H_param

        # lets have zeros for unobserved variables and parameters
        all_observed = self.params['all_observed']
        vec_inputs = self.params['vec_inputs']
        ndim = n // self.params["total_state_param_vars"]
        vecs, indx_map, dim_per_proc = icesee_get_index(**self.model_kwargs)
        for ii, key in enumerate(vec_inputs):
            if key in all_observed:
                H[:, indx_map[key]] = H[:, indx_map[key]]
            else:
                H[:, indx_map[key]] = 0

        return H

    def Obs_fun(self, virtual_obs):
        """
        Observation operator that reduces the full observation vector into a smaller subset.
        
        Parameters:
        huxg_virtual_obs (numpy array): The virtual observation vector (n-dimensional).
        number_obs_instants (int): The number of observations.
        
        Returns:
        numpy array: The reduced observation vector.
        """
        # check if virtual_obs is a scalar
        if np.isscalar(virtual_obs):
            n = 1
        else:
            n = virtual_obs.shape[0]

        return np.dot(self.H_matrix(n), virtual_obs)

    def JObs_fun(self, n_model):
        """
        Jacobian of the observation operator.
        
        Parameters:
        n_model (int): The size of the model state vector.
        number_obs_instants (int): The number of observations.
        
        Returns:
        numpy array: The Jacobian matrix of the observation operator.
        """

        return self.H_matrix(n_model)

    
    # def generate_observation_schedule(self, **kwargs):
    #     try:
    #         t = np.array(kwargs["t"])
    #         freq_obs = self.params["freq_obs"]
    #         obs_start_time = self.params["obs_start_time"]
    #         obs_max_time = self.params["obs_max_time"]

    #         max_t = np.max(t)
    #         obs_max_time = min(obs_max_time, max_t)

    #         obs_t = np.arange(obs_start_time, obs_max_time + freq_obs, freq_obs)
    #         obs_t = obs_t[obs_t <= obs_max_time]

    #         obs_idx = []
    #         for time in obs_t:
    #             idx = np.argmin(np.abs(t - time))
    #             obs_idx.append(idx)
    #         obs_idx = np.array(obs_idx, dtype=int)

    #         num_observations = len(obs_idx)
    #         return obs_t, obs_idx, num_observations
    #     except Exception as e:
    #         print(f"Error occurred in generate_observation_schedule: {e}")
    #         tb_str = "".join(traceback.format_exception(*sys.exc_info()))
    #         print(f"Traceback details:\n{tb_str}")
    #         # self.mpi_comm.Abort(1)

    def generate_observation_schedule(self, **kwargs):
        try:
            import numpy as np

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
                return np.array([]), np.array([], dtype=int), 0, np.array([])

            # --- Build ideal observation times ---
            n_obs = int(np.floor((obs_max - obs_start) / freq_obs)) + 1
            obs_t_req = obs_start + np.arange(n_obs, dtype=float) * freq_obs

            # --- Match model time points to observation times ---
            dt_grid = np.min(np.diff(t)) if len(t) > 1 else 1.0
            tol = 1e-6 * dt_grid

            # For each t in the model time grid, check if it’s close to any obs time
            obs_idx = []
            for i, ti in enumerate(t):
                if np.any(np.abs(ti - obs_t_req) < tol):
                    obs_idx.append(i)

            obs_idx = np.array(obs_idx, dtype=int)
            obs_t_aligned = t[obs_idx]
            num_observations = len(obs_idx)

            return obs_t_req, obs_idx, num_observations
        except Exception as e:
            print(f"Error occurred in generate_observation_schedule: {e}")
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")
    
    # --- Create synthetic observations ---
    # def _create_synthetic_observations(self,**kwargs):
    #     """create synthetic observations"""
    #     statevec_true = kwargs.get('statevec_true', None)
    #     nd, nt = statevec_true.shape

    #     obs_t, ind_m, m_obs = self.generate_observation_schedule(**kwargs)

    #     vecs, indx_map, _ = icesee_get_index(statevec_true, **kwargs)

    #     # create synthetic observations
    #     hu_obs = np.zeros((nd,self.params["number_obs_instants"]))

    #     # check if params["sig_obs"] is a scalar
    #     # if isinstance(self.params["sig_obs"], (int, float)):
    #     #     self.params["sig_obs"] = np.ones(self.params["nt"]+1) * self.params["sig_obs"]
    #     if kwargs.get('joint_estimation', False) or self.params.get('localization_flag', False):
    #         hdim = statevec_true.shape[0] // self.params["total_state_param_vars"]
    #     else:
    #         hdim = statevec_true.shape[0] // self.params["total_state_param_vars"]
    #         nd = hdim * self.params["num_state_vars"]


    #     error_R = np.zeros((nd, m_obs * 2 + 1))
    #     for i, sig in enumerate(self.params["sig_obs"]):
    #         start_idx = i*hdim
    #         end_idx = start_idx + hdim
    #         error_R[start_idx:end_idx,:] = np.ones((hdim,1)) * sig

    #     # print(f"[ICESEE] vec_inputs: {kwargs['vec_inputs']}")
    #     bed_flag = (key == 'bed' or key == 'bedrock' or key == 'bed_topography' or key == 'bedtopo' or key == 'bedtopography')
    #     km = 0
    #     km_temp = 0
    #     for step in range(nt):
    #         if (km<m_obs) and (step+1 == ind_m[km]):
    #             # for key in kwargs['vec_inputs']:
    #             for ii, key in enumerate(kwargs['vec_inputs']):
    #                 if (ii < kwargs['num_state_vars'] or key in kwargs.get('observed_params', [])) and not bed_flag:
    #                     hu_obs[indx_map[key],km] = statevec_true[indx_map[key],step+1] + np.random.normal(0,error_R[indx_map[key],km],len(indx_map[key]))
    #                 else:
    #                     # fill with zeros for parameters
    #                     hu_obs[indx_map[key],km] = np.zeros(len(indx_map[key]))

    #                 if key == 'bed' or key == 'bedrock' or key == 'bed_topography' or key == 'bedtopo' or key == 'bedtopography':
    #                     if step+1 == kwargs.get('bed_obs_snapshot', 0)[km_temp]:
    #                         hu_obs[indx_map[key],km] = statevec_true[indx_map[key],step+1] + np.random.normal(0,error_R[indx_map[key],km],len(indx_map[key]))
    #                         km_temp += 1

    #             km += 1

    #     return hu_obs, error_R.T

    def _create_synthetic_observations(self, **kwargs):
        """create synthetic observations (same logic; bed snapshots are sparse by spacing/mask)"""
        import numpy as np

        statevec_true = kwargs.get('statevec_true', None)
        assert statevec_true is not None, "statevec_true is required"
        # nd, nt = statevec_true.shape
        params = kwargs.get('params')
        nd = kwargs.get('nd', params.get('nd', statevec_true.shape[0]))
        nt = kwargs.get('nt', params.get('nt', statevec_true.shape[1]))

        # Observation schedule
        obs_t, ind_m, m_obs = self.generate_observation_schedule(**kwargs)
        ind_m = np.asarray(ind_m, dtype=int)  # 1-based
        print(f"[ICESEE] observation times: {obs_t}, indices: {ind_m}, total: {m_obs}")

        # Index maps
        vecs, indx_map, _ = icesee_get_index(statevec_true, **kwargs)
        vec_inputs = list(kwargs['vec_inputs'])

        # Preallocate observations
        # hu_obs = np.zeros((nd, self.params["number_obs_instants"]))
        hu_obs = np.zeros((nd, m_obs))

        # hdim / nd handling (unchanged)
        total_state_param_vars = self.params["total_state_param_vars"]
        hdim = statevec_true.shape[0] // total_state_param_vars
        if not (kwargs.get('joint_estimation', False) or self.params.get('localization_flag', False)):
            nd = hdim * self.params["num_state_vars"]

        # Build error_R (nd, m_obs*2+1), unchanged
        error_R = np.zeros((nd, m_obs * 2 + 1))
        sig_obs = self.params["sig_obs"]
        for i, sig in enumerate(sig_obs):
            start_idx = i * hdim
            end_idx = start_idx + hdim
            error_R[start_idx:end_idx, :] = sig  # broadcast

        # Fast memberships & cached indices
        observed_params = set(kwargs.get('observed_params', []))
        bed_aliases = {'bed', 'bedrock', 'bed_topography', 'bedtopo', 'bedtopography'}
        key_is_bed = {k: (k in bed_aliases) for k in vec_inputs}
        key_idx_map = {k: np.asarray(indx_map[k], dtype=int) for k in vec_inputs}

        # ---- Bed snapshot control (same timing logic) ----
        bed_snaps = kwargs.get('bed_obs_snapshot', [])
        bed_snaps = list(bed_snaps) if isinstance(bed_snaps, (list, np.ndarray)) else []
        km = 0
        km_temp = 0

        # ---- Build a sparse-observation mask for each bed-like key ----
        # Options (use whichever you pass in):
        #   - bed_obs_indices: exact indices into bed subvector (0..len(bed)-1)
        #   - bed_obs_mask: boolean mask over bed subvector
        #   - bed_obs_stride_km: spacing (km) to convert to every-nth point along x
        #   - bed_obs_spacing: spacing (points) = every n-th point
        # Needs Lx for stride_km->points if you use it. Pull from kwargs or params.
        Lx = kwargs.get('Lx', self.params.get('Lx', None))
        bed_stride_km = kwargs.get('bed_obs_stride_km', None)
        bed_spacing_pts = kwargs.get('bed_obs_spacing', None)
        bed_indices_user = kwargs.get('bed_obs_indices', None)
        bed_mask_user = kwargs.get('bed_obs_mask', None)

        bed_mask_map = {}
        for k in vec_inputs:
            if not key_is_bed[k]:
                continue
            bed_idx = key_idx_map[k]               # global indices into statevec
            local_len = bed_idx.size               # size of the bed subvector

            # Default: observe all (original behavior)
            mask = np.ones(local_len, dtype=bool)

            # Priority 1: explicit mask
            if isinstance(bed_mask_user, (list, np.ndarray)):
                mask = np.asarray(bed_mask_user, dtype=bool)
                if mask.size != local_len:
                    # Try to safely broadcast smaller masks by stepping
                    if mask.ndim == 1 and mask.size > 0:
                        step = int(np.ceil(local_len / mask.size))
                        rep = int(np.ceil(local_len / mask.size))
                        mask = np.tile(mask, rep)[:local_len]
                    else:
                        mask = np.ones(local_len, dtype=bool)

            # Priority 2: explicit indices
            elif isinstance(bed_indices_user, (list, np.ndarray)):
                mask = np.zeros(local_len, dtype=bool)
                idxs = np.asarray(bed_indices_user, dtype=int)
                idxs = idxs[(idxs >= 0) & (idxs < local_len)]
                mask[idxs] = True

            # Priority 3: spacing in points
            elif isinstance(bed_spacing_pts, (int, np.integer)) and bed_spacing_pts > 1:
                n = int(bed_spacing_pts)
                mask = np.zeros(local_len, dtype=bool)
                mask[::n] = True

            # Priority 4: spacing in km (needs Lx)
            elif (bed_stride_km is not None) and (Lx is not None):
                # Convert stride in km → every-nth point, assuming uniform x-grid
                # Use (hdim-1) intervals to estimate dx; fall back safely.
                intervals = max(hdim - 1, 1)
                dx_m = float(Lx) / intervals
                # n = max(int(round((bed_stride_km * 1000.0) / max(dx_m, 1e-12))), 1)
                n = max(int(round((bed_stride_km) / max(dx_m, 1e-12))), 1)
                mask = np.zeros(local_len, dtype=bool)
                mask[::n] = True

            # Save per-key mask (over the bed subvector)
            bed_mask_map[k] = mask

        obs_set = set(kwargs["observed_vars"] + kwargs["observed_params"])
        # Map each bed snapshot time to the obs column
        bed_time_to_col = {t: col for col, t in enumerate(ind_m)}

        for step in range(nt):
            if (km < m_obs) and (step + 1 == ind_m[km]):
                for key in vec_inputs:
                    idx = key_idx_map[key]
                    bed_flag = key_is_bed[key]

                    # ---------- STANDARD VARIABLES ----------
                    if key in obs_set and not bed_flag:
                        sigma = error_R[idx, km]
                        hu_obs[idx, km] = (
                            statevec_true[idx, step+1] +
                            np.random.normal(0.0, sigma, size=idx.size)
                        )
                    else:
                        hu_obs[idx, km] = 0.0

                    # ---------- BED SPECIAL CASE --------------
                    if bed_flag and (step+1 in bed_time_to_col):
                        col = bed_time_to_col[step+1]
                        mask = bed_mask_map[key]
                        idx_obs = idx[mask]
                        if idx_obs.size > 0:
                            sigma_obs = error_R[idx_obs, col]
                            hu_obs[idx_obs, col] = (
                                statevec_true[idx_obs, step+1] +
                                np.random.normal(0.0, sigma_obs, size=idx_obs.size)
                            )
                        # print("Nonzero bed obs:", np.nonzero(hu_obs[bed_ind,:]))

                km += 1

        return hu_obs, error_R.T

    
    def bed(self, x):
        """
        Bed topography function, which computes the bed shape based on input x and model parameters.
        
        Parameters:
        x (jax.numpy array): Input spatial grid points.
        
        Returns:
        jax.numpy array: The bed topography values at each x location.
        """
        import jax.numpy as jnp
        
        # Ensure parameters are floats
        params     = self.params
        sillamp    = float(params['sillamp'])
        sillsmooth = float(params['sillsmooth'])
        xsill      = float(params['xsill'])

        # Compute the bed topography
        b = sillamp * (-2 * jnp.arccos((1 - sillsmooth) * jnp.sin(jnp.pi * x / (2 * xsill))) / jnp.pi - 1)
        return b
    
    def compute_smb_mask(self,  k, km,  state_block_size, hu_obs=None, smb_init=None, smb_clim=None, model_kwargs=None):
        """
        Compute a robust SMB mask based on observations (if available) or ensemble statistics.
        
        Parameters:
        - statevec_ens: np.array, current state ensemble
        - state_block_size: int, starting index for SMB-related states
        - hu_obs: np.array or None, observed SMB values (optional)
        - smb_init: np.array, initial SMB state
        - smb_clim: np.array or None, climatological SMB values (optional)
        - model_kwargs: dict, model parameters containing time `t`
        - params: dict, additional model parameters (e.g., `nt`)

        Returns:
        - Updated `statevec_ens` with a robust SMB mask applied.
        """
        statevec_ens = self.ensemble
        params = self.params

        # Define the SMB state
        smb = statevec_ens[state_block_size:, :]

        # Time information
        t = model_kwargs["t"]
        nt = params["nt"]
        time_factor = np.clip(t[k] / (t[nt - 1] - t[0]), 0, 1)  # Prevent division by zero

        # If SMB observations exist, use them to set a dynamic threshold
        # if hu_obs is not None:
        if np.max(hu_obs[state_block_size:, km]) > 0:
            smb_crit = np.percentile(hu_obs[state_block_size:, km], 95, axis=0) * 1.25
        # elif smb_clim is not None:
        elif np.max(smb_clim) > 0:
            # If climatological data exists, use it as a reference
            smb_crit = smb_clim + np.std(smb, axis=1)
        else:
            # Default to ensemble statistics if no observations or climatology
            smb_mean = np.mean(smb, axis=1)
            smb_std = np.std(smb, axis=1)
            smb_crit = smb_mean + 2.0 * smb_std

        # Compute lower and upper bounds
        smb_crit_upper = smb_crit
        smb_crit_lower = smb_crit * 0.8  # Allow some variability

        # Logical mask for significant deviations
        smb_flag = np.any((smb > smb_crit_upper) | (smb < smb_crit_lower), axis=1)

        if np.any(smb_flag):  # If any deviations exist
            alpha = 0.05  # Smoothing factor
            correction_factor = (1 - np.exp(-alpha * time_factor))

            # Apply smooth correction
            statevec_ens[state_block_size:, :] = smb_init + (smb - smb_init) * correction_factor

        return statevec_ens



    