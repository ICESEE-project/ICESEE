# ==============================================================================
# @des: This file contains run functions for any error generation
# supprorts: - "fft": old fast spectral method for uniform grids
#            -  "auto": choose automatically
#            -  "random_fields": use gstools to generate random fields with specified covariance
#            - "graph": sparse smoothing on topology/connectivity           
# @date: 2025-07-30
# @author: Brian Kyanjo
# ==============================================================================
    
# --- Imports ---
import numpy as np
import warnings
from scipy.optimize import brentq
import scipy.sparse as sp


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

def generate_pseudo_random_field_1d_(N, Lx, rh, grid_extension=2, verbose=False, **kwargs):
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


# -----> debuging
def generate_pseudo_random_field_1d(N=None, Lx=None, rh=None, grid_extension=2, verbose=False, **kwargs):
    """
    Generate a 1D pseudo-random field with zero mean, unit variance, and specified covariance.

    Backward compatible:
      - If called with only the original arguments, behavior remains the same.
      - Optional kwargs enable automatic handling for coords/connectivity/nonuniform grids.

    Original parameters
    -------------------
    N : int
        Number of grid points
    Lx : float
        Physical domain size
    rh : float
        Decorrelation length for covariance
    grid_extension : int
        Factor to extend grid to avoid periodicity
    verbose : bool
        Print diagnostics

    New optional kwargs
    -------------------
    method : {"auto", "fft", "graph"}, default="auto"
    coords : array_like, optional
        1D coordinates of length N for nonuniform grids.
    connectivity : sparse matrix or edge list, optional
        Graph connectivity for nonuniform/unstructured topology.
    seed : int, optional
        Random seed.
    num_passes : int or "auto", optional
        Graph smoothing passes.
    blend : float, optional
        Graph smoothing blend.
    k_neighbors : int, optional
        Number of neighbors for coords-based graph construction.

    Returns
    -------
    q : ndarray of shape (N,)
    """
    import numpy as np
    import warnings
    from scipy.optimize import brentq
    import scipy.sparse as sp

    try:
        from scipy.spatial import cKDTree
        _HAVE_KDTREE = True
    except Exception:
        _HAVE_KDTREE = False

    # ------------------------------------------------------------------
    # New optional controls. If none are provided, old behavior is used.
    # ------------------------------------------------------------------
    method = kwargs.get("method", "auto")
    coords = kwargs.get("coords", None)
    connectivity = kwargs.get("connectivity", None)
    seed = kwargs.get("seed", 42)
    num_passes = kwargs.get("num_passes", "auto")
    blend = kwargs.get("blend", 0.55)
    k_neighbors = kwargs.get("k_neighbors", 6)


    # if seed is not None:
    #     np.random.seed(int(seed))

    def _is_uniform_1d_coords(x, rtol=1e-6, atol=1e-12):
        x = np.asarray(x, dtype=float).ravel()
        if x.size < 3:
            return True
        dxs = np.diff(x)
        return np.allclose(dxs, dxs[0], rtol=rtol, atol=atol)

    def _build_chain_adjacency(n):
        rows, cols, data = [], [], []
        for i in range(n - 1):
            rows += [i, i + 1]
            cols += [i + 1, i]
            data += [1.0, 1.0]
        return sp.csr_matrix((data, (rows, cols)), shape=(n, n))

    def _build_connectivity_adjacency(conn, n):
        if sp.issparse(conn):
            return conn.tocsr().astype(float)

        conn = np.asarray(conn, dtype=int)
        if conn.ndim != 2 or conn.shape[1] != 2:
            raise ValueError("connectivity must be a sparse matrix or edge list of shape (E, 2).")

        i = conn[:, 0]
        j = conn[:, 1]
        rows = np.concatenate([i, j])
        cols = np.concatenate([j, i])
        data = np.ones(len(rows), dtype=float)
        return sp.csr_matrix((data, (rows, cols)), shape=(n, n))

    def _build_knn_adjacency(x, k):
        x = np.asarray(x, dtype=float)
        if x.ndim == 1:
            x = x[:, None]

        n = x.shape[0]
        if not _HAVE_KDTREE:
            raise ImportError("scipy.spatial.cKDTree is required for coords-based graph mode.")

        tree = cKDTree(x)
        kq = min(k + 1, n)
        dists, inds = tree.query(x, k=kq)

        rows, cols, data = [], [], []

        valid = dists[:, 1:].ravel()
        valid = valid[valid > 0]
        eps = np.median(valid) if valid.size else 1.0
        eps = max(eps, 1e-12)

        for i in range(n):
            for dist, j in zip(np.atleast_1d(dists[i])[1:], np.atleast_1d(inds[i])[1:]):
                if i == j:
                    continue
                w = np.exp(-(dist / eps) ** 2)
                rows.append(i)
                cols.append(j)
                data.append(w)

        W = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
        W = 0.5 * (W + W.T)
        return W

    def _graph_field_1d(n, x=None, conn=None):
        if conn is not None:
            W = _build_connectivity_adjacency(conn, n)
            mode_used = "connectivity"
        elif x is not None:
            x = np.asarray(x)
            if len(x) != n:
                raise ValueError(f"coords length {len(x)} does not match N={n}")
            W = _build_knn_adjacency(x, k_neighbors)
            mode_used = "coords"
        else:
            W = _build_chain_adjacency(n)
            mode_used = "chain"

        deg = np.asarray(W.sum(axis=1)).ravel()
        deg_safe = np.where(deg > 0, deg, 1.0)
        P = sp.diags(1.0 / deg_safe) @ W

        if num_passes == "auto":
            passes = int(np.clip(round(0.12 * np.sqrt(max(n, 4))), 3, 25))
        else:
            passes = int(num_passes)

        q = np.random.randn(n)
        q -= np.mean(q)

        for _ in range(passes):
            q = (1.0 - blend) * q + blend * (P @ q)

        q -= np.mean(q)
        std = np.std(q)
        if std > 0:
            q /= std

        if verbose:
            print(f"[ICESEE] graph mode = {mode_used}")
            print(f"[ICESEE] graph passes = {passes}")
            print(f"[ICESEE] graph blend = {blend}")
            print(f"[ICESEE] graph var = {np.var(q)}")
            print(f"[ICESEE] graph mean = {np.mean(q)}")

        return np.asarray(q).reshape(n,)

    # ------------------------------------------------------------------
    # AUTO selection. If nothing new is passed, fall through to old FFT.
    # ------------------------------------------------------------------
    if method == "auto":
        if connectivity is not None:
            method = "graph"
        elif coords is not None:
            coords_arr = np.asarray(coords).ravel()
            if coords_arr.size != N:
                raise ValueError(f"coords size {coords_arr.size} does not match N={N}")
            if _is_uniform_1d_coords(coords_arr):
                method = "fft"
                if Lx is None:
                    Lx = float(coords_arr.max() - coords_arr.min()) if N > 1 else 1.0
                if rh is None:
                    rh = Lx / 10.0
            else:
                method = "graph"
        else:
            method = "fft"

    if method == "graph":
        return _graph_field_1d(N, x=coords, conn=connectivity)

    # ------------------------------------------------------------------
    # Original FFT code path below: preserved as-is in spirit.
    # ------------------------------------------------------------------
    # Grid spacing
    if Lx is None:
        Lx = float(N)
    dx = Lx / N

    if rh is None:
        rh = Lx / 10.0

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
        sigma_values = np.logspace(-6, 6, 25)
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
            warnings.warn("Could not find a bracketing interval. Using heuristic sigma based on rh.")
            sigma = 2 / rh
            if verbose:
                print(f"[ICESEE] Fallback sigma: {sigma}")
    else:
        try:
            sigma = brentq(covariance_eq, a, b, rtol=1e-6)
            if verbose:
                print(f"[ICESEE] Solved sigma: {sigma}")
        except ValueError as e:
            warnings.warn(f"brentq failed: {str(e)}. Using heuristic sigma.")
            sigma = 2 / rh
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
    self_conj_mask = (I == I_conj)
    mask_representative = (I <= I_conj)

    phi[mask_representative & ~self_conj_mask] = np.random.rand(np.sum(mask_representative & ~self_conj_mask))
    phi[~mask_representative] = (-phi[I_conj[~mask_representative]]) % 1

    # Fourier coefficients
    b_q = A * np.exp(2j * np.pi * phi)

    # Inverse FFT to get the field
    q_ext = np.real(np.fft.ifft(b_q) * N_ext)

    # Crop to original domain
    q = q_ext[:N]

    # Normalize to ensure unit variance
    q = q - np.mean(q)
    std = np.std(q)
    if std > 0:
        q = q / std

    if verbose:
        print(f"[ICESEE] Field variance: {np.var(q)}")
        print(f"[ICESEE] Field mean: {np.mean(q)}")

    return q


def generate_pseudo_random_field_2D(nx=None, ny=None, Lx=None, Ly=None, rh=None,
                                    grid_extension=2, verbose=False, **kwargs):
    """
    Generate a 2D pseudo-random field with zero mean, unit variance, and specified covariance.

    Backward compatible:
      - If called with only the original arguments, behavior remains FFT-based.
      - Optional kwargs enable automatic handling for coords/connectivity/nonuniform grids.

    Original parameters
    -------------------
    nx, ny : int
        Number of grid points in x and y.
    Lx, Ly : float
        Physical domain sizes in x and y.
    rh : float
        Decorrelation length for covariance.
    grid_extension : int
        Factor to extend grid to avoid periodicity.
    verbose : bool
        Print diagnostics.

    New optional kwargs
    -------------------
    method : {"auto", "fft", "graph"}, default="auto"
    coords : array_like, optional
        Coordinates of shape (nx*ny, 2) or (ny, nx, 2) for nonuniform grids.
    connectivity : sparse matrix or edge list, optional
        Graph connectivity for nonuniform/unstructured topology.
    seed : int, optional
        Random seed.
    num_passes : int or "auto", optional
        Graph smoothing passes.
    blend : float, optional
        Graph smoothing blend.
    k_neighbors : int, optional
        Number of neighbors for coords-based graph construction.

    Returns
    -------
    q : ndarray of shape (ny, nx)
    """
    import numpy as np
    import warnings
    from scipy.optimize import brentq
    import scipy.sparse as sp

    try:
        from scipy.spatial import cKDTree
        _HAVE_KDTREE = True
    except Exception:
        _HAVE_KDTREE = False

    # ------------------------------------------------------------------
    # New optional controls. If none are provided, old behavior is used.
    # ------------------------------------------------------------------
    method = kwargs.get("method", "auto")
    coords = kwargs.get("coords", None)
    connectivity = kwargs.get("connectivity", None)
    seed = kwargs.get("seed", None)
    num_passes = kwargs.get("num_passes", "auto")
    blend = kwargs.get("blend", 0.55)
    k_neighbors = kwargs.get("k_neighbors", 6)

    if seed is not None:
        np.random.seed(seed)

    def _reshape_coords_2d(xy, nx_, ny_):
        xy = np.asarray(xy, dtype=float)
        if xy.ndim == 3:
            if xy.shape != (ny_, nx_, 2):
                raise ValueError(f"coords shape {xy.shape} must be (ny, nx, 2)=({ny_}, {nx_}, 2)")
            return xy.reshape(ny_ * nx_, 2)
        elif xy.ndim == 2:
            if xy.shape != (nx_ * ny_, 2):
                raise ValueError(f"coords shape {xy.shape} must be (nx*ny, 2)=({nx_ * ny_}, 2)")
            return xy
        else:
            raise ValueError("coords must have shape (ny, nx, 2) or (nx*ny, 2)")

    def _is_uniform_2d_coords(xy, nx_, ny_, rtol=1e-6, atol=1e-12):
        xy = _reshape_coords_2d(xy, nx_, ny_)
        X = xy[:, 0].reshape(ny_, nx_)
        Y = xy[:, 1].reshape(ny_, nx_)

        # x should vary regularly across columns, y across rows
        dx_rows = np.diff(X, axis=1)
        dy_cols = np.diff(Y, axis=0)

        uniform_x = True if dx_rows.size == 0 else np.allclose(dx_rows, dx_rows[0, 0], rtol=rtol, atol=atol)
        uniform_y = True if dy_cols.size == 0 else np.allclose(dy_cols, dy_cols[0, 0], rtol=rtol, atol=atol)

        # also check rectilinearity: x constant down columns, y constant across rows
        rect_x = np.allclose(X, X[0:1, :], rtol=rtol, atol=atol)
        rect_y = np.allclose(Y, Y[:, 0:1], rtol=rtol, atol=atol)

        return uniform_x and uniform_y and rect_x and rect_y

    def _build_grid_adjacency_2d(nx_, ny_):
        rows, cols, data = [], [], []

        def idx(j, i):
            return j * nx_ + i

        for j in range(ny_):
            for i in range(nx_):
                p = idx(j, i)

                if i + 1 < nx_:
                    q = idx(j, i + 1)
                    rows += [p, q]
                    cols += [q, p]
                    data += [1.0, 1.0]

                if j + 1 < ny_:
                    q = idx(j + 1, i)
                    rows += [p, q]
                    cols += [q, p]
                    data += [1.0, 1.0]

        return sp.csr_matrix((data, (rows, cols)), shape=(nx_ * ny_, nx_ * ny_))

    def _build_connectivity_adjacency(conn, n):
        if sp.issparse(conn):
            return conn.tocsr().astype(float)

        conn = np.asarray(conn, dtype=int)
        if conn.ndim != 2 or conn.shape[1] != 2:
            raise ValueError("connectivity must be a sparse matrix or edge list of shape (E, 2).")

        i = conn[:, 0]
        j = conn[:, 1]
        rows = np.concatenate([i, j])
        cols = np.concatenate([j, i])
        data = np.ones(len(rows), dtype=float)
        return sp.csr_matrix((data, (rows, cols)), shape=(n, n))

    def _build_knn_adjacency(xy, k):
        xy = np.asarray(xy, dtype=float)
        n = xy.shape[0]

        if not _HAVE_KDTREE:
            raise ImportError("scipy.spatial.cKDTree is required for coords-based graph mode.")

        tree = cKDTree(xy)
        kq = min(k + 1, n)
        dists, inds = tree.query(xy, k=kq)

        rows, cols, data = [], [], []

        valid = dists[:, 1:].ravel()
        valid = valid[valid > 0]
        eps = np.median(valid) if valid.size else 1.0
        eps = max(eps, 1e-12)

        for i in range(n):
            for dist, j in zip(np.atleast_1d(dists[i])[1:], np.atleast_1d(inds[i])[1:]):
                if i == j:
                    continue
                w = np.exp(-(dist / eps) ** 2)
                rows.append(i)
                cols.append(j)
                data.append(w)

        W = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
        W = 0.5 * (W + W.T)
        return W

    def _graph_field_2d(nx_, ny_, xy=None, conn=None):
        n = nx_ * ny_

        if conn is not None:
            W = _build_connectivity_adjacency(conn, n)
            mode_used = "connectivity"
        elif xy is not None:
            xy = _reshape_coords_2d(xy, nx_, ny_)
            W = _build_knn_adjacency(xy, k_neighbors)
            mode_used = "coords"
        else:
            W = _build_grid_adjacency_2d(nx_, ny_)
            mode_used = "grid"

        deg = np.asarray(W.sum(axis=1)).ravel()
        deg_safe = np.where(deg > 0, deg, 1.0)
        P = sp.diags(1.0 / deg_safe) @ W

        if num_passes == "auto":
            passes = int(np.clip(round(0.12 * np.sqrt(max(n, 4))), 3, 25))
        else:
            passes = int(num_passes)

        q = np.random.randn(n)
        q -= np.mean(q)

        for _ in range(passes):
            q = (1.0 - blend) * q + blend * (P @ q)

        q -= np.mean(q)
        std = np.std(q)
        if std > 0:
            q /= std

        q = np.asarray(q).reshape(ny_, nx_)

        if verbose:
            print(f"[ICESEE] graph mode = {mode_used}")
            print(f"[ICESEE] graph passes = {passes}")
            print(f"[ICESEE] graph blend = {blend}")
            print(f"[ICESEE] graph var = {np.var(q)}")
            print(f"[ICESEE] graph mean = {np.mean(q)}")

        return q

    # ------------------------------------------------------------------
    # AUTO selection. If nothing new is passed, fall through to old FFT.
    # ------------------------------------------------------------------
    if method == "auto":
        if connectivity is not None:
            method = "graph"
        elif coords is not None:
            if _is_uniform_2d_coords(coords, nx, ny):
                method = "fft"
                xy = _reshape_coords_2d(coords, nx, ny)
                X = xy[:, 0].reshape(ny, nx)
                Y = xy[:, 1].reshape(ny, nx)
                if Lx is None:
                    Lx = float(X.max() - X.min()) if nx > 1 else 1.0
                if Ly is None:
                    Ly = float(Y.max() - Y.min()) if ny > 1 else 1.0
                if rh is None:
                    rh = min(Lx, Ly) / 10.0
            else:
                method = "graph"
        else:
            method = "fft"

    if method == "graph":
        return _graph_field_2d(nx, ny, xy=coords, conn=connectivity)

    # ------------------------------------------------------------------
    # Original FFT-style 2D code path
    # ------------------------------------------------------------------
    if Lx is None:
        Lx = float(nx)
    if Ly is None:
        Ly = float(ny)

    dx = Lx / nx
    dy = Ly / ny

    if rh is None:
        rh = min(Lx, Ly) / 10.0

    if rh < min(dx, dy):
        warnings.warn(
            f"Decorrelation length rh={rh} is smaller than grid spacing min(dx,dy)={min(dx,dy)}. "
            "Consider increasing rh, decreasing Lx/Ly, or increasing nx/ny."
        )

    nx_ext = int(nx * grid_extension)
    ny_ext = int(ny * grid_extension)

    kx = np.fft.fftfreq(nx_ext, d=dx) * 2.0 * np.pi
    ky = np.fft.fftfreq(ny_ext, d=dy) * 2.0 * np.pi
    KX, KY = np.meshgrid(kx, ky, indexing="xy")

    k2 = KX**2 + KY**2

    dkx = 2.0 * np.pi / (nx_ext * dx)
    dky = 2.0 * np.pi / (ny_ext * dy)
    dk = dkx * dky

    def covariance_eq(sigma):
        exp_term = np.exp(-2.0 * k2 / sigma**2)
        numerator = np.sum(exp_term * np.cos(KX * rh))
        denominator = np.sum(exp_term)
        return numerator / denominator - np.exp(-1.0)

    a, b = 1e-6, 100
    fa = covariance_eq(a)
    fb = covariance_eq(b)

    if verbose:
        print(f"[ICESEE] covariance_eq at sigma={a}: {fa}")
        print(f"[ICESEE] covariance_eq at sigma={b}: {fb}")

    if fa * fb > 0:
        warnings.warn("Initial interval [1e-6, 100] does not bracket a root. Trying to find a new interval.")
        sigma_values = np.logspace(-6, 6, 25)
        f_values = [covariance_eq(s) for s in sigma_values]

        if verbose:
            print("[ICESEE] Testing sigma values:")
            for s, f in zip(sigma_values, f_values):
                print(f"[ICESEE] sigma={s:.2e}, covariance_eq={f:.2e}")

        for i in range(len(f_values) - 1):
            if f_values[i] * f_values[i + 1] < 0:
                a, b = sigma_values[i], sigma_values[i + 1]
                fa, fb = f_values[i], f_values[i + 1]
                break
        else:
            warnings.warn("Could not find a bracketing interval. Using heuristic sigma based on rh.")
            sigma = 2 / rh
            if verbose:
                print(f"[ICESEE] Fallback sigma: {sigma}")
    else:
        try:
            sigma = brentq(covariance_eq, a, b, rtol=1e-6)
            if verbose:
                print(f"[ICESEE] Solved sigma: {sigma}")
        except ValueError as e:
            warnings.warn(f"brentq failed: {str(e)}. Using heuristic sigma.")
            sigma = 2 / rh
            if verbose:
                print(f"[ICESEE] Fallback sigma: {sigma}")

    sum_exp = np.sum(np.exp(-2.0 * k2 / sigma**2))
    c2 = 1.0 / (dk * sum_exp)
    c = np.sqrt(c2)

    if verbose:
        print(f"[ICESEE] Computed c: {c}")

    A = c * np.sqrt(dk) * np.exp(-k2 / sigma**2)

    # 2D Hermitian-symmetric random phases
    phi = np.zeros((ny_ext, nx_ext))
    done = np.zeros((ny_ext, nx_ext), dtype=bool)

    for j in range(ny_ext):
        for i in range(nx_ext):
            jc = (-j) % ny_ext
            ic = (-i) % nx_ext

            if done[j, i]:
                continue

            if (j == jc) and (i == ic):
                phi[j, i] = 0.0
                done[j, i] = True
            else:
                r = np.random.rand()
                phi[j, i] = r
                phi[jc, ic] = (-r) % 1.0
                done[j, i] = True
                done[jc, ic] = True

    b_q = A * np.exp(2j * np.pi * phi)

    q_ext = np.real(np.fft.ifft2(b_q) * (nx_ext * ny_ext))
    q = q_ext[:ny, :nx]

    q = q - np.mean(q)
    std = np.std(q)
    if std > 0:
        q = q / std

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

# def generate_enkf_field(ii_sig, Lx, hdim, num_vars, rh=None, grid_extension=2, verbose=False, field_kwargs=None):
def generate_enkf_field(**kwargs):
    """
    Generate a pseudo-random field for EnKF with specified DoF.

    Parameters expected in kwargs
    -----------------------------
    ii_sig : int or None
        Variable index when generating one variable at a time.
    Lx : float or None
        Representative length scale.
    hdim : int
        Degrees of freedom per variable.
    num_vars : int
        Number of variables.
    rh : float, list, ndarray, or None
        Decorrelation length(s).
    grid_extension : int, optional
        FFT grid extension factor.
    verbose : bool, optional
        Print diagnostics.

    Additional kwargs are passed through to generate_pseudo_random_field_1d,
    e.g.:
        method, coords, connectivity, seed, coords_by_var,
        connectivity_by_var, num_passes, blend, k_neighbors, ...

    Returns
    -------
    q : ndarray
        Array of shape (hdim * num_vars,) or (hdim,)
    """
    import numpy as np

    # ------------------------------------------------------------
    # unpack core kwargs with defaults
    # ------------------------------------------------------------
    ii_sig = kwargs.get("ii_sig", None)
    Lx = kwargs.get("Lx_dim", None)
    hdim = kwargs.get("noise_dim", None)
    num_vars = kwargs.get("num_vars", None)
    rh = kwargs.get("rh", None)
    grid_extension = kwargs.get("grid_extension", 2)
    verbose = kwargs.get("verbose", False)
    kwargs['method'] = kwargs.get("method", "fft")

    if Lx == 1:
        Lx = None

    if hdim is None:
        raise ValueError("generate_enkf_field requires 'hdim'.")
    if num_vars is None:
        raise ValueError("generate_enkf_field requires 'num_vars'.")

    if rh is None:
        if Lx is not None:
            rh = Lx / 10.0
        else:
            rh = max(float(hdim) / 10.0, 1.0)

    # ------------------------------------------------------------
    # kwargs to pass down to generate_pseudo_random_field_1d
    # remove keys already consumed here
    # ------------------------------------------------------------
    passthrough_kwargs = dict(kwargs)
    for key in ["ii_sig", "Lx", "hdim", "num_vars", "rh", "grid_extension", "verbose"]:
        passthrough_kwargs.pop(key, None)

    def _local_kwargs(var_index=None):
        local_kwargs = dict(passthrough_kwargs)

        if var_index is not None:
            if "coords_by_var" in local_kwargs and "coords" not in local_kwargs:
                local_kwargs["coords"] = local_kwargs["coords_by_var"][var_index]
            if "connectivity_by_var" in local_kwargs and "connectivity" not in local_kwargs:
                local_kwargs["connectivity"] = local_kwargs["connectivity_by_var"][var_index]

        # do not forward these containers further down
        local_kwargs.pop("coords_by_var", None)
        local_kwargs.pop("connectivity_by_var", None)

        return local_kwargs

    # ------------------------------------------------------------
    # Handle trivial case: no spatial dimension
    # preserve existing behavior
    # ------------------------------------------------------------
    if hdim < 1e2:
        if verbose:
            print(f"[ICESEE] hdim={hdim} small — using FFT exp-cov sampling (no dense cov).")

        if isinstance(rh, (list, np.ndarray)):
            if ii_sig is None:
                q_total = []
                for i in range(num_vars):
                    var_rh = rh[i]
                    q_var = sample_periodic_exp_cov(hdim, var_rh, Lx)
                    q_total.append(q_var)
                return np.concatenate(q_total, axis=0)
            else:
                return sample_periodic_exp_cov(hdim, rh[ii_sig], Lx)
        else:
            if ii_sig is None:
                return sample_periodic_exp_cov(hdim * num_vars, rh, Lx)
            else:
                return sample_periodic_exp_cov(hdim, rh, Lx)

    # ------------------------------------------------------------
    # Main branch
    # preserve old output shapes exactly
    # ------------------------------------------------------------
    if isinstance(rh, (list, np.ndarray)):

        if ii_sig is None:
            q_total = []
            for i in range(num_vars):
                var_rh = rh[i]

                q_var = generate_pseudo_random_field_1d(
                    N=hdim,
                    Lx=Lx,
                    rh=var_rh,
                    grid_extension=grid_extension,
                    verbose=verbose,
                    **_local_kwargs(var_index=i)
                )
                q_total.append(q_var)

            return np.concatenate(q_total, axis=0)

        else:
            q0 = generate_pseudo_random_field_1d(
                N=hdim,
                Lx=Lx,
                rh=rh[ii_sig],
                grid_extension=grid_extension,
                verbose=verbose,
                **_local_kwargs(var_index=ii_sig)
            )
            return q0

    else:
        if ii_sig is None:
            q0 = generate_pseudo_random_field_1d(
                N=hdim * num_vars,
                Lx=Lx,
                rh=rh,
                grid_extension=grid_extension,
                verbose=verbose,
                **_local_kwargs(var_index=None)
            )
        else:
            q0 = generate_pseudo_random_field_1d(
                N=hdim,
                Lx=Lx,
                rh=rh,
                grid_extension=grid_extension,
                verbose=verbose,
                **_local_kwargs(var_index=ii_sig)
            )

        return q0