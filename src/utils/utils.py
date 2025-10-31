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
    def __init__(self, params,ensemble=None):
        """
        Initialize the utility functions with model parameters.
        
        Parameters:
        params (dict): Model parameters, including those used for bed topography.
        """
        self.params   = params
        self.ensemble = ensemble
    
    def H_matrix(self, n_model):
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
            H_param = np.zeros(num_params_size)
            H[:,state_variables_size:] = H_param
            
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
    def _create_synthetic_observations(self,**kwargs):
        """create synthetic observations"""
        statevec_true = kwargs.get('statevec_true', None)
        nd, nt = statevec_true.shape

        obs_t, ind_m, m_obs = self.generate_observation_schedule(**kwargs)

        vecs, indx_map, _ = icesee_get_index(statevec_true, **kwargs)

        # create synthetic observations
        hu_obs = np.zeros((nd,self.params["number_obs_instants"]))

        # check if params["sig_obs"] is a scalar
        # if isinstance(self.params["sig_obs"], (int, float)):
        #     self.params["sig_obs"] = np.ones(self.params["nt"]+1) * self.params["sig_obs"]
        if kwargs.get('joint_estimation', False) or self.params.get('localization_flag', False):
            hdim = statevec_true.shape[0] // self.params["total_state_param_vars"]
        else:
            hdim = statevec_true.shape[0] // self.params["total_state_param_vars"]
            nd = hdim * self.params["num_state_vars"]


        error_R = np.zeros((nd, m_obs * 2 + 1))
        for i, sig in enumerate(self.params["sig_obs"]):
            start_idx = i*hdim
            end_idx = start_idx + hdim
            error_R[start_idx:end_idx,:] = np.ones((hdim,1)) * sig

        # print(f"[ICESEE] vec_inputs: {kwargs['vec_inputs']}")
        km = 0
        for step in range(nt):
            if (km<m_obs) and (step+1 == ind_m[km]):
                # hu_obs[:,km] = statevec_true[:,step+1] + norm(loc=0,scale=self.params["sig_obs"][step+1]).rvs(size=nd)
                # hu_obs[:,km] = statevec_true[:,step+1] + np.random.normal(0,self.params["sig_obs"][step+1],nd)
                # TODO: start from here tomorrow.
                for key in kwargs['vec_inputs']:
                    hu_obs[indx_map[key],km] = statevec_true[indx_map[key],step+1] + np.random.normal(0,error_R[indx_map[key],km],len(indx_map[key]))

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



    