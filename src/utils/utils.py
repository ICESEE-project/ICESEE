# =============================================================================
# @Author: Brian Kyanjo
# @Date: 2024-09-24
# @Description: This script includes the observation operator and its Jacobian 
#               for the EnKF data assimilation scheme. It also includes the bed
#               topography function used in the model.
# =============================================================================

# import libraries
import numpy as np
import re
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

    
    def generate_observation_schedule(self,**kwargs):
        """
        Generate observation times and indices from a given array of time points.

        Parameters:
            t (list or np.ndarray): Array of time points.
            freq_obs (int): Frequency of observations in the same unit as `t`.
            obs_max_time (int): Maximum observation time in the same unit as `t`.

        Returns:
            obs_t (list): Observation times.
            obs_idx (list): Indices corresponding to observation times in `t`.
        """
        # unpack kwargs
        t = kwargs["t"]

        # Convert input to a numpy array for easier manipulation
        t = np.array(t)
        
        # Generate observation times
        obs_t = np.arange(self.params["obs_start_time"], self.params["obs_max_time"] + self.params["freq_obs"], self.params["freq_obs"])
        # obs_t = np.linspace(obs_start_time, obs_max_time, int(obs_max_time/freq_obs)+1)
        
        # Find indices of observation times in the original array
        obs_idx = np.array([np.where(t == time)[0][0] for time in obs_t if time in t]).astype(int)

        # print(f"Number of observation instants: {len(obs_idx)} at times: {t[obs_idx]}")
        
        # number of observation instants
        num_observations = len(obs_idx)

        return obs_t, obs_idx, num_observations
    
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

        error_R = np.zeros((nd, m_obs * 2 + 1))
        for i, sig in enumerate(self.params["sig_obs"]):
            start_idx = i*hdim
            end_idx = start_idx + hdim
            error_R[start_idx:end_idx,:] = np.ones((hdim,1)) * sig

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
    
    def inflate_ensemble(self,in_place=True):
        """
        Inflate ensemble members by a given factor.
        
        Args: 
            ensemble: ndarray (n x N) - The ensemble matrix of model states (n is state size, N is ensemble size).
            inflation_factor: float - scalar or iterable length equal to model states
            in_place: bool - whether to update the ensemble in place
        Returns:
            inflated_ensemble: ndarray (n x N) - The inflated ensemble.
        """
        # check if the inflation factor is scalar
        if np.isscalar(self.params['inflation_factor']):
            _scalar = True
            _inflation_factor = float(self.params['inflation_factor'])
        elif isiterable(self.params['inflation_factor']):
            if len(self.params['inflation_factor']) == self.ensemble.shape[0]:
                _inflation_factor[:] = self.params['inflation_factor'][:]
                _scalar = False
            else:
                raise ValueError("Inflation factor length must be equal to the state size")
        
        # check if we need inflation
        if _scalar:
            if _inflation_factor == 1.0:
                return self.ensemble
            elif _inflation_factor < 0.0:
                raise ValueError("Inflation factor must be positive scalar")
        else:
            _inf = False
            for i in _inflation_factor:
                if i>1.0:
                    _inf = True
                    break
            if not _inf:
                return self.ensemble
        
        ens_size = self.ensemble.shape[1]
        mean_vec = np.mean(self.ensemble, axis=1)
        if in_place:
            inflated_ensemble = self.ensemble
            for ens_idx in range(ens_size):
                state = inflated_ensemble[:, ens_idx]
                if _scalar:
                    state = (state - mean_vec) * _inflation_factor
                else:
                    state = (state - mean_vec) * _inflation_factor

                inflated_ensemble[:, ens_idx] = state + mean_vec
        else:
            inflated_ensemble = np.zeros(self.ensemble.shape)
            for ens_idx in range(ens_size):
                state = self.ensemble[:, ens_idx].copy()
                if _scalar:
                    state = (state - mean_vec) * _inflation_factor
                else:
                    state = (state - mean_vec) * _inflation_factor

                inflated_ensemble[:, ens_idx] = state + mean_vec
        
        return inflated_ensemble
    
    def _inflate_ensemble(self,rescale=False):
        """inflate ensemble members by a given factor"""

        _inflation_factor = float(self.params['inflation_factor'])
        x = np.mean(self.ensemble, axis=0, keepdims=True)
        X = self.ensemble - x

        # rescale the ensemble to correct the variance
        if rescale:
            N, M = self.ensemble.shape
            X *= np.sqrt(N/(N-1))

        x = x.squeeze(axis=0)

        if _inflation_factor == 1.0:
            return self.ensemble
        else:
            return x + _inflation_factor * X


    def _localization_matrix(self,euclidean_distance, localization_radius, loc_type='Gaspari-Cohn'):     
        """
        Calculate the localization matrix based on the localization type, euclean_distance and radius
        of influence.
        
        Parameters:
        euclidean_distance (numpy array): The Euclidean distance between the observation and state
        localization_radius (float or numpy array): Distance beyond which the localization matrix is tapered to zero.
        method (str): The localization method.
        
        Returns:
        numpy array: The localization matrix (same size as the Euclidean distance).
        """

        # Get original shape
        dist_size = euclidean_distance.shape

        # Gaspari-Cohn localization
        if re.match(r'\Agaspari(_|-)*cohn\Z', loc_type, re.IGNORECASE):
            # Normalize distances relative to localization radius
            radius = euclidean_distance.flatten() / (0.5 * localization_radius)

            # Initialize localization matrix with zeros
            localization_matrix = np.zeros_like(radius)

            # Gaspari-Cohn function
            mask0 = radius < 1
            mask1 = (radius >= 1) & (radius < 2)

            # Compute values where radius < 1
            loc_func0 = (((-0.25 * radius + 0.5) * radius + 0.625) * radius - 5.0 / 3.0) * radius**2 + 1
            localization_matrix[mask0] = loc_func0[mask0]

            # Compute values where 1 <= radius < 2
            radius_safe = np.where(radius == 0, 1e-10, radius)  # Avoid division by zero
            loc_func1 = ((((1.0 / 12.0 * radius_safe - 0.5) * radius_safe + 0.625) * radius_safe + 5.0 / 3.0) * radius_safe - 5.0) * radius_safe + 4.0 - 2.0 / 3.0 / radius_safe
            localization_matrix[mask1] = loc_func1[mask1]
            return localization_matrix.reshape(dist_size)
        # Gaussian localization
        elif re.match(r'\Agaussian\Z', loc_type, re.IGNORECASE):
            return np.exp(-0.5 * (euclidean_distance / localization_radius)**2)

        else:
            raise ValueError(f"Unknown localization type: {loc_type}")

    import numpy as np

    def compute_sample_correlations_vectorized(self, shuffled_ens, forward_ens):
        """
        Compute sample correlations between shuffled_ens and forward_ens in a vectorized manner.
        
        Parameters:
            shuffled_ens (np.ndarray): Array of shape (n_members, n_variables) representing the shuffled ensemble.
            forward_ens (np.ndarray): Array of shape (n_members, n_variables) representing the forward ensemble.
        
        Returns:
            np.ndarray: An array of correlation coefficients (one per variable).
        """
        # Number of ensemble members
        Nens = self.ensemble.shape[1]

        # Compute means for each variable (column-wise)
        mean_shuffled = np.mean(shuffled_ens, axis=0)
        mean_forward = np.mean(forward_ens, axis=0)

        # Center the ensembles by subtracting the means
        centered_shuffled = shuffled_ens - mean_shuffled
        centered_forward = forward_ens - mean_forward

        # Compute the covariance for each variable (element-wise multiplication, then sum over rows)
        cov = np.sum(centered_shuffled * centered_forward, axis=0) / (Nens - 1)

        # Compute the standard deviations for each variable with Bessel's correction (ddof=1)
        std_shuffled = np.std(shuffled_ens, axis=0, ddof=1)
        std_forward = np.std(forward_ens, axis=0, ddof=1)

        # Calculate the correlation coefficient for each variable
        correlations = cov / (std_shuffled * std_forward)

        return correlations

    
    def _adaptive_localization(self, euclidean_distance=None, 
                              localization_radius=None, ensemble_init=None, loc_type='Gaspari-Cohn'):
        """Adaptively calculates the radius of influence for each observation density
           which is then used to dynamically compute the localization matrix.
           returns: adaptive localization matrix
        @reference: See https://doi.org/10.1016/j.petrol.2019.106559 for more details
        """

        # get the shape of the ensemble size
        nd, Nens = self.ensemble.shape

        # if localization radius is not provided, use the adaptive method
        if localization_radius is None:
            # correlation based localization
            if Nens >= 30:
                # random shuffle the initial ensemble
                np.random.shuffle(ensemble_init)
                # ensemble members after forward simulation
                forward_ens = self.ensemble

                # get initial sample correlation btn the shuffled and forward ens
                # sample_ind
                # sample_correlations = self.compute_sample_correlations_vectorized(shuffled_ens, forward_ens)
                sample_correlations = np.corrcoef(ensemble_init, forward_ens, rowvar=False)
                
                # # subsitute noise field of sample_correlations 
                # sample_correlations[np.isnan(sample_correlations)] = 0

                # use the MAD rule to estimate noise levels; sig_gs = median(abs(eta_gs))/0.6745
                sig_gs = np.median(np.abs(sample_correlations), axis=0) / 0.6745

                # use the universal rule to subsitute noise fields; theta_gs = sqrt(2*ln(number of rho_gs))*sig_gs
                theta_gs = np.sqrt(2 * np.log(Nens)) * sig_gs

                # construct the tapering matrix by applying the the estimated noise levels 
                # to the sample correlations
                tapering_matrix = np.exp(-0.5 * (sample_correlations / theta_gs)**2)

            # distance based localization
            else:
                # if the dist between the model variable and the observation is zero, then the weight is 1
                if np.any(euclidean_distance == 0): 
                    localization_matrix = np.ones(self.ensemble.shape[0])
                    return localization_matrix

                # use a type based on variance  
                var = np.var(self.ensemble,axis=0)
                avg_var = np.mean(var)
                localization_radius = self.params['base_radius'] * np.sqrt(1 + self.params['scaling_factor'] * np.sqrt(avg_var))

                # call the localization matrix function
                localization_matrix = self._localization_matrix(euclidean_distance, localization_radius)
                return localization_matrix
        else:
            # call the localization matrix function
            localization_matrix = self._localization_matrix(euclidean_distance, localization_radius)
            return localization_matrix

    import numpy as np

    def _adaptive_localization_v2(self, cutoff_distance):
        """
        Compute an adaptive localization matrix based on ensemble correlations.

        Parameters:
        cutoff_distance (numpy array): Predefined cutoff distances for localization.

        Returns:
        numpy array: The computed localization matrix.
        """
        # Get ensemble size
        nd, Nens = self.ensemble.shape

        # Compute correlation matrix
        R = np.corrcoef(self.ensemble, rowvar=False)

        # Compute threshold for localization radius
        rad_flag = 1 / np.sqrt(Nens - 1)

        # Find the first occurrence where correlation drops below threshold
        mask = R < rad_flag  # Boolean mask

        # Get the first (i, j) index where R[i, j] < rad_flag
        indices = np.argwhere(mask)  # Get all (i, j) pairs that satisfy the condition
        
        if indices.size > 0:
            first_i = indices[0, 0]  # First valid row index
            radius = cutoff_distance[first_i]  # Assign corresponding cutoff distance

            # Call the localization matrix function with the scalar radius
            localization_matrix = self._localization_matrix(cutoff_distance, radius)
            return localization_matrix

        # Return None if no valid index found (handle this case as needed)
        return None

    def rmse(self,truth, estimate):
        """
        Calculate the Root Mean Squared Error (RMSE) between the true and estimated values.
        
        Parameters:
        truth (numpy array): The true values.
        estimate (numpy array): The estimated values.
        
        Returns:
        float: The RMSE value.
        """
        return np.sqrt(np.mean((truth - estimate) ** 2))

    def compute_euclidean_distance(self, grid_x, grid_y):
        """
        Compute the Euclidean distance matrix between all grid points.

        Parameters:
        grid_x (numpy array): X-coordinates of the grid points (1D array).
        grid_y (numpy array): Y-coordinates of the grid points (1D array).

        Returns:
        numpy array: Euclidean distance matrix (NxN, where N = number of grid points).
        """
        # Stack X, Y coordinates into (N, 2) array where N is the number of points
        grid_points = np.column_stack((grid_x.ravel(), grid_y.ravel()))

        # Compute pairwise Euclidean distances
        distance_matrix = cdist(grid_points, grid_points, metric='euclidean')

        return distance_matrix
    
    def gaspari_cohn(self,r):
        """
        Compute the Gaspari-Cohn localization function.
        
        Parameters:
        r (numpy array): Normalized distance (d / r0), where d is the Euclidean distance 
                        and r0 is the localization radius.
        
        Returns:
        numpy array: Localization weights corresponding to r.
        """
        gc = np.zeros_like(r)  # Initialize localization weights

        # Case 0 <= r < 1
        mask1 = (r >= 0) & (r < 1)
        gc[mask1] = (((-0.25 * r[mask1] + 0.5) * r[mask1] + 0.625) * r[mask1] - 5.0 / 3.0) * r[mask1]**2 + 1

        # Case 1 <= r < 2
        mask2 = (r >= 1) & (r < 2)
        gc[mask2] = ((((1.0 / 12.0 * r[mask2] - 0.5) * r[mask2] + 0.625) * r[mask2] + 5.0 / 3.0) * r[mask2] - 5.0) * r[mask2] + 4.0 - 2.0 / (3.0 * np.where(r[mask2] == 0, 1e-10, r[mask2]))

        # Case r >= 2 (default to 0)
        return gc
    
    def create_tapering_matrix(self,grid_x, grid_y, localization_radius):
        """
        Create a tapering matrix using the Gaspari-Cohn localization function.

        Parameters:
        grid_x (numpy array): X-coordinates of grid points (1D array).
        grid_y (numpy array): Y-coordinates of grid points (1D array).
        localization_radius (float): Cutoff radius beyond which correlations are zero.

        Returns:
        numpy array: Tapering matrix (NxN), where N = number of grid points.
        """
        # Compute Euclidean distance matrix
        distance_matrix = self.compute_euclidean_distance(grid_x, grid_y)
        # print(distance_matrix)

        # Normalize distances by the localization radius
        # if is radius is a scalar
        if np.isscalar(localization_radius):
            r = distance_matrix / (0.5*localization_radius)
        else:
            if localization_radius.shape[0] == distance_matrix.shape[0]:
                r = distance_matrix / 0.5*localization_radius[:, None]
                # r = np.ones_like(distance_matrix)*localization_radius[:, None]
            elif localization_radius.shape[0] > distance_matrix.shape[0]:  
                obs_indices = np.arange(distance_matrix.shape[0])  # Select only the required points
                r = distance_matrix / 0.5*localization_radius[obs_indices, None]
                # r = np.ones_like(distance_matrix)*localization_radius[obs_indices, None]

        # Normalize distances by the localization radius
        # r = distance_matrix / localization_radius
        if False:
            # create a localization matrix without distance
            r = np.ones_like(distance_matrix)*localization_radius

        # Compute tapering matrix using Gaspari-Cohn function
        tapering_matrix = self.gaspari_cohn(r)
        # print(f"tapering matrix: {tapering_matrix}")

        return tapering_matrix

    def compute_adaptive_localization_radius(self, grid_x, grid_y, base_radius=2.0, method='variance'):
        """
        Compute an adaptive localization radius for each grid point.

        Parameters:
        ensemble (numpy array): Ensemble state matrix (N_grid x N_ens).
        grid_x (numpy array): X-coordinates of grid points (1D array).
        grid_y (numpy array): Y-coordinates of grid points (1D array).
        base_radius (float): Default radius before adaptation.
        method (str): 'variance', 'observation_density', or 'correlation'.

        Returns:
        numpy array: Adaptive localization radius for each grid point.
        """
        num_points, Nens = self.ensemble.shape  # Get grid size and ensemble size
        adaptive_radius = np.full(num_points, base_radius)  # Default radius

        if method == 'variance':
            # Compute ensemble variance at each grid point
            ensemble_variance = np.var(self.ensemble, axis=1)

            # Normalize variance (relative to max spread)
            normalized_variance = ensemble_variance / np.max(ensemble_variance)

            # Scale localization radius based on variance
            adaptive_radius *= (1 + normalized_variance)

        elif method == 'observation_density':
            # Compute observation density (using a Gaussian kernel approach)
            grid_points = np.column_stack((grid_x.ravel(), grid_y.ravel()))
            obs_density = np.sum(np.exp(-cdist(grid_points, grid_points, 'euclidean')**2 / base_radius**2), axis=1)

            # Normalize observation density
            normalized_density = obs_density / np.max(obs_density)

            # Decrease localization radius in high-density regions
            adaptive_radius *= (1 - normalized_density)

        elif method == 'correlation':
            # Compute correlation matrix from the ensemble
            correlation_matrix = np.corrcoef(self.ensemble, rowvar=True)

            # Set radius where correlation drops below 1/sqrt(Nens-1)
            threshold = 1 / np.sqrt(Nens - 1)
            for i in range(num_points):
                below_threshold = np.where(correlation_matrix[i, :] < threshold)[0]
                if below_threshold.size > 0:
                    adaptive_radius[i] = base_radius * np.min(below_threshold) / num_points  # Scale adaptively

        else:
            raise ValueError("Invalid method. Choose 'variance', 'observation_density', or 'correlation'.")

        return adaptive_radius

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



    