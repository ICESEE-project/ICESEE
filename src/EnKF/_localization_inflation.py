# =============================================================================
# @Author: Brian Kyanjo
# @Date: 2024-09-24
# @Description: This script includes localization and inflation functions used 
#               in the EnKF data assimilation scheme.
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
from scipy.stats import pearsonr
from scipy.spatial import distance_matrix
from scipy.optimize import curve_fit
from scipy.sparse import csr_matrix
from scipy.ndimage import uniform_filter


# import utility functions
from ICESEE.src.utils.tools import icesee_get_index


# --- helper functions ---
def isiterable(obj):
    return isinstance(obj, Iterable)

class LocalizationInflationUtils:
    """
    A class containing utility functions for localization and inflation in the
    Ensemble Kalman Filter (EnKF) data assimilation scheme.
    """
    def __init__(self, params,ensemble=None):
        """
        Initialize the utility functions with model parameters.
        
        Parameters:
        params (dict): Model parameters, including those used for bed topography.
        """
        self.params   = params
        self.ensemble = ensemble

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
        
    def gaspari_cohn(self, r):
        """
        Gaspari-Cohn taper function for localization in EnKF.
        Defined for 0 <= r <= 2.
        """
        r = np.abs(r)
        taper = np.zeros_like(r)
    
        mask1 = (r >= 0) & (r <= 1)
        mask2 = (r > 1) & (r <= 2)

        taper[mask1] = (((-0.25 * r[mask1] + 0.5) * r[mask1] + 0.625) * r[mask1] - 5/3) * r[mask1]**2 + 1
        taper[mask2] = ((((1/12 * r[mask2] - 0.5) * r[mask2] + 0.625) * r[mask2] + 5/3) * r[mask2] - 5) * r[mask2]**2 + 4 - 2/(3 * r[mask2])

    def localization(self, Lx,Ly,nx, ny, n_points, n_vars, n_ens, state_size):
        # Generate grid points
        x = np.linspace(0, Lx, int(np.sqrt(n_points * Lx / Ly)))
        y = np.linspace(0, Ly, int(np.sqrt(n_points * Ly / Lx)))
        X, Y = np.meshgrid(x, y)
        grid_points = np.vstack([X.ravel(), Y.ravel()]).T[:n_points]
        # grid_points = np.pad(grid_points, ((0, n_points - grid_points.shape[0]), (0, 0)), 'constant', constant_values=0)
        # Extrapolate missing values using the last row values (close to Lx and Ly)
        missing_rows = n_points - grid_points.shape[0]
        if missing_rows > 0:
            last_row = grid_points[-1]  # Get the last available row
            extrapolated_rows = np.tile(last_row, (missing_rows, 1))  # Repeat last row
            grid_points = np.vstack([grid_points, extrapolated_rows])  # Append extrapolated rows


        # Distance matrix (only valid points)
        dist_matrix = distance_matrix(grid_points[:n_points], grid_points[:n_points])  # 425 × 425
        max_distance = np.max(dist_matrix)

        # Generate ensemble
        np.random.seed(42)
        ensemble_vec = np.random.randn(state_size, n_ens)  # 1700 × 24

        # Function to compute spatially varying L
        def compute_spatial_L(ensemble, grid_points, threshold=0.1, min_L=0, max_L=max_distance):
            """
            Compute spatially varying localization length scale L for each point or region.
            
            Parameters:
            - ensemble: Ensemble matrix (state_size × n_ens)
            - grid_points: Grid point coordinates (n_points × 2)
            - threshold: Correlation threshold below which to define L
            - min_L, max_L: Bounds for L
            
            Returns:
            - L_array: Array of L values for each point (n_points,)
            """
            n_points = grid_points.shape[0]
            n_ens = ensemble.shape[1]
            n_vars = state_size // n_points
            
            # Compute ensemble mean and anomalies for the first variable block
            ens_mean = np.mean(ensemble, axis=1)
            ens_anom = ensemble - ens_mean[:, np.newaxis]  # Anomalies (state_size × n_ens)
            ens_block = ens_anom[:n_points, :]  # 425 × 24 (first variable block)
            
            # Compute correlations for each point with all others
            L_array = np.zeros(n_points)
            for i in range(n_points):
                correlations = np.zeros(n_points)
                for j in range(n_points):
                    if i != j:
                        corr, _ = pearsonr(ens_block[i, :], ens_block[j, :])
                        correlations[j] = corr if not np.isnan(corr) else 0
                    else:
                        correlations[j] = 1.0
                
                # Take absolute correlations
                correlations = np.abs(correlations)
                
                # Sort distances and correlations for point i
                distances = dist_matrix[i, :n_points]
                mask = distances > 0
                dists = distances[mask]
                corrs = correlations[mask]
                
                sorted_pairs = sorted(zip(dists, corrs))
                dists, corrs = zip(*sorted_pairs)
                dists, corrs = np.array(dists), np.array(corrs)
                
                # Find distance where correlation drops below threshold
                if np.max(corrs) <= threshold:
                    L = min_L
                else:
                    L = dists[np.where(corrs <= threshold)[0][0]] if threshold in corrs else dists[-1]
                    L = max(min_L, min(max_L, L))  # Clip to reasonable range
                
                L_array[i] = L
            
            # Optionally smooth L_array spatially (e.g., using a moving average)
            L_array = np.clip(L_array, min_L, max_L)  # Ensure bounds
            return L_array

        # Compute spatially varying L
        L_array = compute_spatial_L(ensemble_vec, grid_points)
        print(f"Spatially Varying Localization Length Scales L (min, max): {np.min(L_array):.2f}, {np.max(L_array):.2f} meters")

        # Generate localization matrix with spatially varying L
        def gaspari_cohn_spatial(r, L):
            """Gaspari-Cohn localization function with spatially varying L."""
            r = np.abs(r)
            out = np.zeros_like(r)
            idx1 = r <= 1.0
            idx2 = (r > 1.0) & (r <= 2.0)
            out[idx1] = 1 - 5/3 * r[idx1]**2 + 5/8 * r[idx1]**3 + 1/2 * r[idx1]**4 - 1/4 * r[idx1]**5
            out[idx2] = -5/3 * r[idx2] + 5/8 * r[idx2]**2 + 1/2 * r[idx2]**3 - 1/4 * r[idx2]**4 + 1/12 * (2/r[idx2] - 1/r[idx2]**4)
            return np.where(r > 2, 0, out)  # Zero beyond 2L

        # Create spatially varying localization matrix
        loc_matrix_spatial = np.zeros((n_points, n_points))
        for i in range(n_points):
            for j in range(n_points):
                dist = dist_matrix[i, j]
                L_i = L_array[i]  # Use L for point i (could average or use j's L)
                loc_matrix_spatial[i, j] = gaspari_cohn_spatial(dist / L_i, L_i)

        # Expand to full state space
        loc_matrix = np.zeros((state_size, state_size))
        for var_i in range(n_vars):
            for var_j in range(n_vars):
                start_i, start_j = var_i * n_points, var_j * n_points
                loc_matrix[start_i:start_i + n_points, start_j:start_j + n_points] = loc_matrix_spatial

        from scipy.sparse import csr_matrix
        # loc_matrix = csr_matrix(loc_matrix)
        return loc_matrix

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
    
    def gaspari_cohn_v0(self,r):
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