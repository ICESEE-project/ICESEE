# ==============================================================================
# @des: This file contains run functions for any localization functions used in the EnKF.
# @date: 2025-07-30
# @author: Brian Kyanjo
# ==============================================================================

import numpy as np
from scipy.stats import pearsonr
from scipy.spatial import distance_matrix
from scipy.optimize import curve_fit
from scipy.sparse import csr_matrix
from scipy.ndimage import uniform_filter

def gaspari_cohn(r):
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

    return np.maximum(taper, 0)  # Ensure non-negative values

def localization(Lx,Ly,nx, ny, n_points, n_vars, n_ens, state_size):
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
