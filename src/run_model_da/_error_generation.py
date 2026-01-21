# ==============================================================================
# @des: This file contains run functions for any error generation
# @date: 2025-07-30
# @author: Brian Kyanjo
# ==============================================================================
    
# --- Imports ---
import os
import sys
import gc # garbage collector to free up memory
import copy
import re
import time
import h5py
import numpy as np
from tqdm import tqdm 
import bigmpi4py as BM # BigMPI for large data transfer and communication
from scipy.sparse import csr_matrix
from scipy.sparse import block_diag
from scipy.stats import multivariate_normal
from scipy.spatial import distance_matrix
from scipy.linalg import cholesky, solve_triangular


def compute_Q_err_random_fields(hdim, num_blocks, sig_Q, rho, len_scale):
    """
    """
    import numpy as np
    import gstools as gs

    gs.config.USE_GSTOOLS_CORE = True

    pos = np.arange(hdim).reshape(-1, 1)
    model = gs.Gaussian(dim=1, var=1, len_scale=len_scale)

    sig_Q_sq = [s**2 for s in sig_Q]
    C = np.zeros((num_blocks, num_blocks))
    outer = np.outer(sig_Q, sig_Q)       # shape: (num_blocks, num_blocks)
    C = rho * outer                      # initialize with off-diagonal terms
    np.fill_diagonal(C, sig_Q_sq)       # set the diagonal elements

    try:
        L_C = np.linalg.cholesky(C)
    except np.linalg.LinAlgError:
        eps = 1e-6
        C  += np.eye(C.shape[0]) * eps
        L_C = np.linalg.cholesky(C)

    return pos, model, L_C

def compute_noise_random_fields(k, hdim, pos, model, num_blocks, L_C):
    import numpy as np
    import gstools as gs

    Y = np.zeros((hdim, num_blocks))
    for i in range(num_blocks):
        srf_i = gs.SRF(model, seed=k * num_blocks + i)
        Y[:, i] = srf_i(pos).flatten()
    X = Y @ L_C.T
    total_noise_k = X.flatten()
    # all_noise.append(total_noise_k)
    return total_noise_k

# -----> debuging
def generate_pseudo_random_field_1d(N, Lx, rh, grid_extension=2, verbose=False):
    """
    Generate a 1D pseudo-random field with zero mean, unit variance, and specified covariance.
    
    Parameters:
    - N: Number of grid points
    - Lx: Physical domain size
    - rh: Decorrelation length for covariance
    - grid_extension: Factor to extend grid to avoid periodicity (default=2)
    - verbose: If True, print diagnostic information (default=False)
    
    Returns:
    - q: 1D array of shape (N,) containing the random field
    """

    import numpy as np
    from scipy.optimize import brentq
    import warnings
    
    # Grid spacing
    dx = Lx / N
    # dx = Lx/nx
    
    # rh = min(min(Lx)/10,rh)

    # Validate parameters
    if rh < dx:
        warnings.warn(f"Decorrelation length rh={rh} is smaller than grid spacing dx={dx}. "
                      "Consider increasing rh, decreasing Lx, or increasing N.")
    
    # Extended grid to avoid periodicity
    N_ext = int(N * grid_extension)
    
    # Wave numbers
    kx = np.fft.fftfreq(N_ext, d=dx) * 2 * np.pi
    
    # Delta k for Fourier summation
    dk = 2 * np.pi / (N_ext * dx)
    
    # Compute sigma by solving the covariance equation
    def covariance_eq(sigma):
        k2 = kx**2
        exp_term = np.exp(-2 * k2 / sigma**2)
        numerator = np.sum(exp_term * np.cos(kx * rh))
        denominator = np.sum(exp_term)
        return numerator / denominator - np.exp(-1)
    
    # Dynamically find a bracketing interval
    a, b = 1e-6, 100
    fa = covariance_eq(a)
    fb = covariance_eq(b)
    
    if verbose:
        print(f"[ICESEE] covariance_eq at sigma={a}: {fa}")
        print(f"[ICESEE] covariance_eq at sigma={b}: {fb}")
    
    # Try expanding the interval if signs are the same
    if fa * fb > 0:
        warnings.warn("Initial interval [1e-6, 100] does not bracket a root. Trying to find a new interval.")
        sigma_values = np.logspace(-6, 6, 25)  # Test a wider, finer range
        f_values = [covariance_eq(s) for s in sigma_values]
        
        if verbose:
            print("[ICESEE] Testing sigma values:")
            for s, f in zip(sigma_values, f_values):
                print(f"[ICESEE] sigma={s:.2e}, covariance_eq={f:.2e}")
        
        # Find a sign change
        for i in range(len(f_values) - 1):
            if f_values[i] * f_values[i + 1] < 0:
                a, b = sigma_values[i], sigma_values[i + 1]
                fa, fb = f_values[i], f_values[i + 1]
                break
        else:
            # Fallback: Estimate sigma based on rh
            warnings.warn("Could not find a bracketing interval. Using heuristic sigma based on rh.")
            sigma = 2 / rh  # Heuristic: sigma ~ 2/rh
            if verbose:
                print(f"[ICESEE] Fallback sigma: {sigma}")
    else:
        # Solve for sigma
        try:
            sigma = brentq(covariance_eq, a, b, rtol=1e-6)
            if verbose:
                print(f"[ICESEE] Solved sigma: {sigma}")
        except ValueError as e:
            warnings.warn(f"brentq failed: {str(e)}. Using heuristic sigma.")
            sigma = 2 / rh  # Fallback
            if verbose:
                print(f"[ICESEE] Fallback sigma: {sigma}")
    
    # Compute c from variance condition
    k2 = kx**2
    sum_exp = np.sum(np.exp(-2 * k2 / sigma**2))
    c2 = 1 / (dk * sum_exp)
    c = np.sqrt(c2)
    
    if verbose:
        print(f"[ICESEE] Computed c: {c}")
    
    # Compute amplitude
    A = c * np.sqrt(dk) * np.exp(-k2 / sigma**2)
    
    # Generate random phases with Hermitian symmetry
    phi = np.zeros(N_ext)
    I = np.arange(N_ext)
    I_conj = np.mod(-I, N_ext)
    self_conj_mask = (I == I_conj)  # Points where k=0 or k=pi
    mask_representative = (I <= I_conj)  # Choose half of the spectrum
    
    # Set phases: zero for self-conjugate points, random for representatives
    phi[mask_representative & ~self_conj_mask] = np.random.rand(np.sum(mask_representative & ~self_conj_mask))
    phi[~mask_representative] = (-phi[I_conj[~mask_representative]]) % 1
    
    # Fourier coefficients
    b_q = A * np.exp(2j * np.pi * phi)
    
    # Inverse FFT to get the field
    q_ext = np.real(np.fft.ifft(b_q) * N_ext)
    
    # Crop to original domain
    q = q_ext[:N]
    
    # Normalize to ensure unit variance
    q = q / np.std(q) * 1.0
    
    if verbose:
        print(f"[ICESEE] Field variance: {np.var(q)}")
        print(f"[ICESEE] Field mean: {np.mean(q)}")
    
    return q


def generate_pseudo_random_field_2D(N, M, Lx, Ly, rh, grid_extension=2, verbose=False):
    """
    Generate a 2D pseudo-random field with zero mean, unit variance, and specified covariance.
    
    Parameters:
    - N, M: Grid points in x and y directions
    - Lx, Ly: Physical domain sizes in x and y directions
    - rh: Decorrelation length for covariance
    - grid_extension: Factor to extend grid to avoid periodicity (default=2)
    - verbose: If True, print diagnostic information (default=False)
    
    Returns:
    - q: 2D array of shape (N, M) containing the random field
    """

    import numpy as np
    from scipy.optimize import brentq
    import warnings

    # Grid spacing
    dx = Lx / N
    dy = Ly / M
    
    # Validate parameters
    if rh < dx or rh < dy:
        warnings.warn(f"Decorrelation length rh={rh} is smaller than grid spacing (dx={dx}, dy={dy}). "
                      "Consider increasing rh to be at least dx or dy, or decreasing Lx, Ly, or increasing N, M.")
    
    # Extended grid to avoid periodicity
    N_ext = int(N * grid_extension)
    M_ext = int(M * grid_extension)
    
    # Wave numbers
    kx = np.fft.fftfreq(M_ext, d=dx) * 2 * np.pi
    ky = np.fft.fftfreq(N_ext, d=dy) * 2 * np.pi
    KY, KX = np.meshgrid(ky, kx, indexing='ij')
    
    # Delta k for Fourier summation
    dk = (2 * np.pi)**2 / (N_ext * M_ext * dx * dy)
    
    # Compute sigma by solving the covariance equation
    def covariance_eq(sigma):
        k2 = KX**2 + KY**2
        exp_term = np.exp(-2 * k2 / sigma**2)
        numerator = np.sum(exp_term * np.cos(KX * rh))
        denominator = np.sum(exp_term)
        return numerator / denominator - np.exp(-1)
    
    # Dynamically find a bracketing interval
    a, b = 1e-6, 100
    fa = covariance_eq(a)
    fb = covariance_eq(b)
    
    if verbose:
        print(f"[ICESEE] covariance_eq at sigma={a}: {fa}")
        print(f"[ICESEE] covariance_eq at sigma={b}: {fb}")
    
    # Try expanding the interval if signs are the same
    if fa * fb > 0:
        warnings.warn("Initial interval [1e-6, 100] does not bracket a root. Trying to find a new interval.")
        sigma_values = np.logspace(-6, 6, 25)  # Test a wider, finer range
        f_values = [covariance_eq(s) for s in sigma_values]
        
        if verbose:
            print("[ICESEE] Testing sigma values:")
            for s, f in zip(sigma_values, f_values):
                print(f"[ICESEE] sigma={s:.2e}, covariance_eq={f:.2e}")
        
        # Find a sign change
        for i in range(len(f_values) - 1):
            if f_values[i] * f_values[i + 1] < 0:
                a, b = sigma_values[i], sigma_values[i + 1]
                fa, fb = f_values[i], f_values[i + 1]
                break
        else:
            # Fallback: Estimate sigma based on rh
            warnings.warn("Could not find a bracketing interval. Using heuristic sigma based on rh.")
            sigma = 2 / rh  # Heuristic: sigma ~ 2/rh from Gaussian covariance approximation
            if verbose:
                print(f"[ICESEE] Fallback sigma: {sigma}")
    else:
        # Solve for sigma
        try:
            sigma = brentq(covariance_eq, a, b, rtol=1e-6)
            if verbose:
                print(f"[ICESEE] Solved sigma: {sigma}")
        except ValueError as e:
            warnings.warn(f"brentq failed: {str(e)}. Using heuristic sigma.")
            sigma = 2 / rh  # Fallback
            if verbose:
                print(f"[ICESEE] Fallback sigma: {sigma}")
    
    # Compute c from variance condition
    k2 = KX**2 + KY**2
    sum_exp = np.sum(np.exp(-2 * k2 / sigma**2))
    c2 = 1 / (dk * sum_exp)
    c = np.sqrt(c2)
    
    if verbose:
        print(f"[ICESEE] Computed c: {c}")
    
    # Compute amplitude
    A = c * np.sqrt(dk) * np.exp(-k2 / sigma**2)
    
    # Generate random phases with Hermitian symmetry
    phi = np.zeros((N_ext, M_ext))
    I, J = np.meshgrid(np.arange(N_ext), np.arange(M_ext), indexing='ij')
    I_conj = np.mod(-I, N_ext)
    J_conj = np.mod(-J, M_ext)
    self_conj_mask = (I == I_conj) & (J == J_conj)
    mask_representative = (I < I_conj) | ((I == I_conj) & (J <= J_conj))
    
    # Set phases: zero for self-conjugate points, random for representatives
    phi[mask_representative & ~self_conj_mask] = np.random.rand(np.sum(mask_representative & ~self_conj_mask))
    phi[~mask_representative] = (-phi[I_conj[~mask_representative], J_conj[~mask_representative]]) % 1
    
    # Fourier coefficients
    b_q = A * np.exp(2j * np.pi * phi)
    
    # Inverse FFT to get the field
    q_ext = np.real(np.fft.ifft2(b_q) * N_ext * M_ext)
    
    # Crop to original domain
    q = q_ext[:N, :M]
    
    # Normalize to ensure unit variance
    q = q / np.std(q) * 1.0
    
    if verbose:
        print(f"[ICESEE] Field variance: {np.var(q)}")
        print(f"[ICESEE] Field mean: {np.mean(q)}")
    
    return q

def sample_periodic_exp_cov(hdim: int, sigma2: float, Lx: float, rng=None):
    """
    Sample x ~ N(0, C) where C_ij = sigma2 * exp(-d(i,j)/Lx),
    d(i,j) = min(|i-j|, hdim-|i-j|) (periodic ring distance).

    Uses circulant diagonalization via FFT: C = F^* diag(lam) F.
    Returns a real sample of shape (hdim,).
    """
    if rng is None:
        rng = np.random.default_rng()

    n = int(hdim)
    if n <= 0:
        raise ValueError("hdim must be positive")
    if Lx <= 0:
        raise ValueError("Lx must be > 0")
    if sigma2 < 0:
        raise ValueError("sigma2 must be >= 0")

    # First row of the circulant covariance: c[k] = sigma2 * exp(-min(k, n-k)/Lx)
    k = np.arange(n, dtype=np.float64)
    d = np.minimum(k, n - k)
    c = sigma2 * np.exp(-d / Lx)

    # Eigenvalues of the circulant matrix are FFT of the first row
    lam = np.fft.rfft(c)  # real FFT -> length n//2 + 1, complex in general but should be real-ish
    lam = np.real(lam)

    # Numerical safety: tiny negatives can happen from roundoff
    lam[lam < 0] = 0.0

    # Sample in Fourier domain:
    # For a real spatial signal, rfft coefficients have special structure:
    #   - DC and Nyquist (if present) are real
    #   - others are complex with independent N(0,1) real/imag
    m = lam.shape[0]
    z = np.empty(m, dtype=np.complex128)

    # DC component (pure real)
    z[0] = rng.normal()

    # Nyquist component if n even (pure real)
    if n % 2 == 0:
        z[-1] = rng.normal()
        mid = m - 2
    else:
        mid = m - 1

    # Remaining positive frequencies (complex)
    if mid > 0:
        z[1:1+mid] = rng.normal(size=mid) + 1j * rng.normal(size=mid)

    # Scale by sqrt eigenvalues; rfft/irfft normalization:
    # numpy's irfft returns the time-domain signal with 1/n factor consistent with FFT conventions.
    # To get covariance C, scale by sqrt(lam * n).
    z *= np.sqrt(lam * n)

    x = np.fft.irfft(z, n=n)
    return x.astype(np.float64, copy=False)

def generate_enkf_field(ii_sig, Lx, hdim, num_vars, rh=None, grid_extension=2, verbose=False):
    """
    Generate a pseudo-random field for EnKF with specified DoF.

    Parameters:
    - Lx: Representative length scale (e.g., domain size in x)
    - hdim: Degrees of freedom per variable
    - num_vars: Number of variables
    - rh: Decorrelation length (float or dict with variable-specific values)
    - grid_extension: Factor to extend grid (default=2)
    - verbose: Print diagnostics (default=False)

    Returns:
    - q: Array of shape (hdim * num_vars, 1)
    """
    # N = hdim * num_vars
    if rh is None:
        rh = Lx / 10  # Default decorrelation length

    # Handle trivial case: no spatial dimension
    if hdim < 1e2:
        if verbose:
            print(f"[ICESEE] hdim={hdim} small — using FFT exp-cov sampling (no dense cov).")

        if isinstance(rh, (list, np.ndarray)):
            if ii_sig is None:
                # Separate fields for each variable
                q_total = []
                for i in range(num_vars):
                    # var_rh = rh.get(f'var{i+1}', Lx / 10)
                    var_rh = rh[i] if isinstance(rh, list) else rh
                    q_var = sample_periodic_exp_cov(hdim, var_rh, Lx)
                    q_total.append(q_var)
                return np.concatenate(q_total, axis=0)
            else:
                # we are in the for loop for perturbation update already
                q0 = sample_periodic_exp_cov(hdim, rh[ii_sig], Lx)
                return q0
        else:
            # Single field
            if ii_sig is None:
                q0 = sample_periodic_exp_cov(hdim * num_vars, rh, Lx)
            else:
                q0 = sample_periodic_exp_cov(hdim, rh, Lx)
            return q0
    else:
        # check if rh is a array
        if isinstance(rh, (list, np.ndarray)):
        
            if ii_sig is None:
                # Separate fields for each variable
                q_total = []
                for i in range(num_vars):
                    # var_rh = rh.get(f'var{i+1}', Lx / 10)
                    var_rh = rh[i] if isinstance(rh, list) else rh
                    q_var = generate_pseudo_random_field_1d(
                        N=hdim, Lx=Lx, rh=var_rh, grid_extension=grid_extension, verbose=verbose
                    )
                    q_total.append(q_var)
                return np.concatenate(q_total, axis=0)
            else:
                # we are in the for loop for perturbation update already
                q0 = generate_pseudo_random_field_1d(
                    N=hdim, Lx=Lx, rh=rh[ii_sig], grid_extension=grid_extension, verbose=verbose
                )
                return q0
        else:
            # Single field
            if ii_sig is None:
                q0 = generate_pseudo_random_field_1d(
                    N=hdim*num_vars, Lx=Lx, rh=rh, grid_extension=grid_extension, verbose=verbose
                )
            else:
                q0 = generate_pseudo_random_field_1d(
                    N=hdim, Lx=Lx, rh=rh, grid_extension=grid_extension, verbose=verbose
                )
            # print(f"[ICESEE] Field shape: {q0.shape}")
            return q0