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
import functools

from ICESEE.src.utils.random_streams import initialization_seed


def _quantize_coordinate_key(coord, tolerance):
    """Deterministic, decomposition-independent integer key for one
    physical coordinate.

    Quantizes to ``tolerance`` (in the coordinate array's own units --
    typically meters for Icepack) so that harmless floating-point
    differences -- e.g. the identical physical mesh assembled by a
    different Firedrake/PETSc parallel decomposition, which can differ in
    the last few ULPs from a different parallel assembly/summation order
    -- still produce the same key. The caller is responsible for
    verifying this tolerance is well inside the mesh's actual minimum DOF
    spacing (see ``check_coordinate_uniqueness`` in
    ``src/tests/parallel_mpi/_coordinate_aware_compare.py``); this
    function does not itself detect a too-loose tolerance collapsing two
    distinct physical DOFs onto the same key.
    """
    scale = 1.0 / float(tolerance)
    return tuple(
        int(np.round(float(c) * scale)) & 0xFFFFFFFF
        for c in np.atleast_1d(coord)
    )


def coordinate_keyed_white_noise(coords, seed, coord_tolerance=1e-6):
    """Zero-mean, unit-variance white noise keyed by (seed, physical
    coordinate) instead of array position.

    ``np.random.randn(n)`` (the previous implementation of graph random
    fields' initial white-noise vector, before this fix) assigns value
    ``i`` to whichever DOF happens to occupy array position ``i`` --
    decomposition-dependent under a partitioned Firedrake mesh, since the
    same physical DOF occupies a different array position under a
    different rank count/decomposition. This instead assigns the SAME
    value to the SAME physical coordinate regardless of what array
    position it occupies -- (base_seed, ensemble member, variable,
    initialization-vs-process-noise namespace, [timestep for process
    noise]) are all already folded into ``seed`` by the caller before it
    reaches here (see ``src/utils/random_streams.py``'s
    ``initialization_seed``/``process_noise_seed``, and
    ``generate_initial_member_increment``/``add_member_process_noise``,
    which compute that seed and either pass it explicitly or seed
    NumPy's legacy global state with it before calling down into this
    module) -- so init-vs-process-noise separation and member/variable/
    timestep independence are inherited automatically, not re-derived
    here. Never keys on MPI rank, local array index, or any other
    decomposition-dependent numbering.

    Only used for the graph method's coords-based branch (physical
    coordinates are its own documented, intended notion of DOF identity
    -- see ``_graph_field_1d``'s "coords" mode). The chain/connectivity-
    only branches (no physical coordinates available) are unchanged.
    """
    coords = np.atleast_2d(coords)
    n = coords.shape[0]
    seed_int = int(seed) & 0xFFFFFFFF
    out = np.empty(n, dtype=np.float64)
    for i in range(n):
        key_ints = [
            seed_int,
            *_quantize_coordinate_key(coords[i], coord_tolerance),
            0x6EA9D0F,
        ]
        per_dof_seed = np.random.SeedSequence(key_ints).generate_state(
            1, dtype=np.uint32
        )[0]
        out[i] = np.random.default_rng(per_dof_seed).standard_normal()
    return out


@functools.lru_cache(maxsize=256)
def _cached_fft_calibration(N_ext, dx, rh):
    """
    Deterministic part of generate_pseudo_random_field_1d's FFT branch:
    wavenumbers, dk, and the amplitude array A. Identical to the
    original inline logic, just cached on (N_ext, dx, rh) since none of
    it depends on randomness.

    Cache-key note: N_ext is an int (exact), dx and rh are floats. As
    long as a given run passes the same literal Lx/rh/N values each call
    (the normal case - these come from config, not from a fresh
    computation each time), float equality holds and the cache hits
    reliably. A miss just means falling back to recomputing - never an
    incorrect result, only a missed speedup.
    """
    kx = np.fft.fftfreq(N_ext, d=dx) * 2 * np.pi
    dk = 2 * np.pi / (N_ext * dx)

    def covariance_eq(sigma):
        k2 = kx**2
        exp_term = np.exp(-2 * k2 / sigma**2)
        return np.sum(exp_term * np.cos(kx * rh)) / np.sum(exp_term) - np.exp(-1)

    a, b = 1e-6, 100
    fa, fb = covariance_eq(a), covariance_eq(b)
    sigma = None

    if fa * fb > 0:
        warnings.warn(
            "Initial interval [1e-6, 100] does not bracket a root. "
            "Trying to find a new interval."
        )
        sigma_values = np.logspace(-6, 6, 25)
        f_values = [covariance_eq(s) for s in sigma_values]
        for i in range(len(f_values) - 1):
            if f_values[i] * f_values[i + 1] < 0:
                a, b = sigma_values[i], sigma_values[i + 1]
                sigma = brentq(covariance_eq, a, b, rtol=1e-6)
                break
        if sigma is None:
            warnings.warn(
                "Could not find a bracketing interval. Using heuristic sigma based on rh."
            )
            sigma = 2 / rh
    else:
        try:
            sigma = brentq(covariance_eq, a, b, rtol=1e-6)
        except ValueError as e:
            warnings.warn(f"brentq failed: {str(e)}. Using heuristic sigma.")
            sigma = 2 / rh

    k2 = kx**2
    sum_exp = np.sum(np.exp(-2 * k2 / sigma**2))
    c = np.sqrt(1.0 / (dk * sum_exp))
    A = c * np.sqrt(dk) * np.exp(-k2 / sigma**2)

    return kx, A

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

def generate_pseudo_random_field_1d_(N, Lx, rh, grid_extension=2, verbose=False, **icesee_kwargs):
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

    kx, A = _cached_fft_calibration(N_ext, dx, rh)

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
def generate_pseudo_random_field_1d(N=None, Lx=None, rh=None, grid_extension=2, verbose=False, **icesee_kwargs):
    """
    Generate a 1D pseudo-random field with zero mean, unit variance, and specified covariance.

    Backward compatible:
      - If called with only the original arguments, behavior remains the same.
      - Optional icesee_kwargs enable automatic handling for coords/connectivity/nonuniform grids.

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

    New optional icesee_kwargs
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
    method = str(
        icesee_kwargs.get(
            "random_field_method",
            icesee_kwargs.get("enkf_field_method", icesee_kwargs.get("method", "auto")),
        )
    ).strip().lower()
    if method not in {"auto", "fft", "graph"}:
        raise ValueError(
            f"Unsupported random-field method {method!r}; "
            "expected 'auto', 'fft', or 'graph'"
        )
    coords = icesee_kwargs.get("coords", None)
    connectivity = icesee_kwargs.get("connectivity", None)
    seed = icesee_kwargs.get("seed", 42)
    num_passes = icesee_kwargs.get("num_passes", "auto")
    blend = icesee_kwargs.get("blend", 0.55)
    k_neighbors = icesee_kwargs.get("k_neighbors", 6)


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

    def _build_knn_adjacency(x, k, coord_tolerance=1e-6):
        x = np.asarray(x, dtype=float)
        if x.ndim == 1:
            x = x[:, None]

        n = x.shape[0]
        if not _HAVE_KDTREE:
            raise ImportError("scipy.spatial.cKDTree is required for coords-based graph mode.")

        tree = cKDTree(x)
        # Query a margin of extra candidates beyond k+1 (self + k
        # neighbors): cKDTree's own tie-breaking for points at EXACTLY
        # the same distance (common on a regular/symmetric mesh) is not
        # guaranteed independent of the order `x` is presented in --
        # decomposition-dependent under a partitioned Firedrake mesh,
        # where the same physical point set can arrive in a different
        # array order. Querying extra candidates and re-selecting the
        # final k ourselves, by a canonical (distance, quantized-
        # coordinate) key that depends only on physical position, makes
        # the selection deterministic regardless of input order. A
        # margin of 8 comfortably covers realistic tie multiplicities
        # (e.g. up to 8-fold-degenerate distances on a regular 2D/3D
        # grid) without materially changing the query cost for a typical
        # k (6 by default).
        tie_margin = 8
        kq = min(k + 1 + tie_margin, n)
        dists, inds = tree.query(x, k=kq)
        dists = np.atleast_2d(dists)
        inds = np.atleast_2d(inds)

        def _canonical_key(idx):
            return _quantize_coordinate_key(x[idx], coord_tolerance)

        per_point_chosen = []
        selected_dists = []
        for i in range(n):
            candidates = [
                (float(dist), int(j))
                for dist, j in zip(dists[i], inds[i])
                if int(j) != i
            ]
            # Stable, decomposition-independent ordering: primarily by
            # distance, then by the candidate's own canonical physical-
            # coordinate key (never by its array index j) -- so an exact
            # distance tie is always broken the same way for the same
            # physical points, regardless of which array position they
            # occupy.
            candidates.sort(key=lambda pair: (pair[0], _canonical_key(pair[1])))
            chosen = candidates[:k]
            per_point_chosen.append(chosen)
            selected_dists.extend(dist for dist, _ in chosen)

        eps = float(np.median(selected_dists)) if selected_dists else 1.0
        eps = max(eps, 1e-12)

        rows, cols, data = [], [], []
        for i, chosen in enumerate(per_point_chosen):
            for dist, j in chosen:
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
            W = _build_knn_adjacency(
                x, k_neighbors,
                coord_tolerance=float(icesee_kwargs.get("coord_key_tolerance", 1e-6)),
            )
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

        if mode_used == "coords":
            # Physical coordinates are graph's own documented notion of
            # DOF identity (unlike "connectivity"/"chain", which have no
            # coordinates to key by and are unchanged here): key the
            # initial white-noise draw to (seed, physical coordinate)
            # rather than array position, so the SAME physical DOF gets
            # the SAME value regardless of which rank/array-index it
            # occupies under a given Firedrake decomposition -- see
            # coordinate_keyed_white_noise's own docstring.
            coord_tolerance = float(icesee_kwargs.get("coord_key_tolerance", 1e-6))
            q = coordinate_keyed_white_noise(x, seed, coord_tolerance=coord_tolerance)
        else:
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
            coords_arr = np.asarray(coords, dtype=float)
            if coords_arr.ndim == 1:
                coords_arr = coords_arr[:, None]
            if coords_arr.ndim != 2 or coords_arr.shape[0] != N:
                raise ValueError(
                    f"coords shape {coords_arr.shape} does not have N={N} rows"
                )
            # A multi-dimensional coordinate array describes a physical mesh,
            # not a scalar uniform line, and must use graph smoothing.
            if coords_arr.shape[1] == 1 and _is_uniform_1d_coords(coords_arr[:, 0]):
                method = "fft"
                if Lx is None:
                    Lx = float(coords_arr[:, 0].max() - coords_arr[:, 0].min()) if N > 1 else 1.0
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
    kx, A = _cached_fft_calibration(N_ext, dx, rh)

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
                                    grid_extension=2, verbose=False, **icesee_kwargs):
    """
    Generate a 2D pseudo-random field with zero mean, unit variance, and specified covariance.

    Backward compatible:
      - If called with only the original arguments, behavior remains FFT-based.
      - Optional icesee_kwargs enable automatic handling for coords/connectivity/nonuniform grids.

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

    New optional icesee_kwargs
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
    method = str(
        icesee_kwargs.get(
            "random_field_method",
            icesee_kwargs.get("enkf_field_method", icesee_kwargs.get("method", "auto")),
        )
    ).strip().lower()
    if method not in {"auto", "fft", "graph"}:
        raise ValueError(
            f"Unsupported random-field method {method!r}; "
            "expected 'auto', 'fft', or 'graph'"
        )
    coords = icesee_kwargs.get("coords", None)
    connectivity = icesee_kwargs.get("connectivity", None)
    seed = icesee_kwargs.get("seed", None)
    num_passes = icesee_kwargs.get("num_passes", "auto")
    blend = icesee_kwargs.get("blend", 0.55)
    k_neighbors = icesee_kwargs.get("k_neighbors", 6)

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
    if Lx is None:
        Lx = float(n)

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
def generate_enkf_field(**icesee_kwargs):
    """
    Generate a pseudo-random field for EnKF with specified DoF.

    Parameters expected in icesee_kwargs
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
        FFT decorrelation length(s). The public YAML key ``length_scale`` is
        accepted as an alias. Neither setting is used by graph mode.
    grid_extension : int, optional
        FFT grid extension factor.
    verbose : bool, optional
        Print diagnostics.

    Additional runtime-context entries are passed to generate_pseudo_random_field_1d,
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
    # unpack core icesee_kwargs with defaults
    # ------------------------------------------------------------
    ii_sig = icesee_kwargs.get("ii_sig", None)
    Lx = icesee_kwargs.get("Lx_dim", None)
    hdim = icesee_kwargs.get("noise_dim", None)
    num_vars = icesee_kwargs.get("num_vars", None)
    # ``length_scale`` is the public YAML spelling used by the application
    # configurations. Keep ``rh`` as the low-level/API spelling, but make the
    # two names genuine aliases so configured FFT scales reach the generator.
    rh = icesee_kwargs.get("rh", None)
    if rh is None:
        configured_rh = icesee_kwargs.get("length_scale", icesee_kwargs.get("len_scale", None))
        if configured_rh is not None and np.asarray(configured_rh).size:
            rh = configured_rh
    grid_extension = icesee_kwargs.get("grid_extension", 2)
    verbose = icesee_kwargs.get("verbose", False)
    method = str(
        icesee_kwargs.get(
            "random_field_method",
            icesee_kwargs.get("enkf_field_method", icesee_kwargs.get("method", "fft")),
        )
    ).strip().lower()
    if method not in {"fft", "graph"}:
        raise ValueError(
            f"Unsupported random-field method {method!r}; expected 'fft' or 'graph'"
        )
    icesee_kwargs["method"] = method

    # Use the caller's member-keyed generator in every sampling branch.  The
    # small-state fast path formerly called ``default_rng()`` implicitly inside
    # ``sample_periodic_exp_cov``; that made otherwise identical execution
    # modes diverge during ensemble initialization.  Falling back to ``seed``
    # keeps direct API calls reproducible as well.
    rng = icesee_kwargs.get("rng")
    if rng is None:
        # ``seed`` is also a historical model/config parameter and can be a
        # YAML float.  The DA-wide base seed is the stable fallback shared by
        # execution modes; accept integral numeric spellings for compatibility.
        seed_value = icesee_kwargs.get(
            "base_seed", icesee_kwargs.get("seed", 42)
        )
        rng = np.random.default_rng(int(seed_value))

    # The run drivers pack registered application coordinates once into
    # icesee_kwargs immediately before ensemble initialization.  Retain a lazy
    # lookup here as well for tests and other direct generator callers.
    if icesee_kwargs.get("coords") is None and icesee_kwargs.get("mesh_coords") is not None:
        icesee_kwargs["coords"] = icesee_kwargs["mesh_coords"]
    if (
        method == "graph"
        and icesee_kwargs.get("coords") is None
        and icesee_kwargs.get("connectivity") is None
        and icesee_kwargs.get("coords_by_var") is None
        and icesee_kwargs.get("connectivity_by_var") is None
    ):
        from ICESEE.src.utils.localization import get_mesh_coordinates

        coords = get_mesh_coordinates(icesee_kwargs)
        if coords is None:
            raise ValueError(
                "random_field_method='graph' requires registered mesh "
                "coordinates, explicit coords, or explicit connectivity"
            )
        icesee_kwargs["coords"] = coords

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
    # icesee_kwargs to pass down to generate_pseudo_random_field_1d
    # remove keys already consumed here
    # ------------------------------------------------------------
    passthrough_kwargs = dict(icesee_kwargs)
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
    if hdim < 1e2 and method != "graph":
        if verbose:
            print(f"[ICESEE] hdim={hdim} small — using FFT exp-cov sampling (no dense cov).")

        if isinstance(rh, (list, np.ndarray)):
            if ii_sig is None:
                q_total = []
                for i in range(num_vars):
                    var_rh = rh[i]
                    q_var = sample_periodic_exp_cov(hdim, var_rh, Lx, rng=rng)
                    q_total.append(q_var)
                return np.concatenate(q_total, axis=0)
            else:
                return sample_periodic_exp_cov(hdim, rh[ii_sig], Lx, rng=rng)
        else:
            if ii_sig is None:
                return sample_periodic_exp_cov(
                    hdim * num_vars, rh, Lx, rng=rng
                )
            else:
                return sample_periodic_exp_cov(hdim, rh, Lx, rng=rng)

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
            if method == "graph":
                # Each variable is defined on the same physical mesh.  Do not
                # concatenate variables first: that would incorrectly require
                # num_vars*hdim distinct coordinates and create graph edges
                # between unrelated state/parameter blocks.
                q0 = np.concatenate([
                    generate_pseudo_random_field_1d(
                        N=hdim,
                        Lx=Lx,
                        rh=rh,
                        grid_extension=grid_extension,
                        verbose=verbose,
                        **_local_kwargs(var_index=i)
                    )
                    for i in range(num_vars)
                ])
            else:
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


def resolve_variable_block_sizes(vec_inputs, vector_size, scalar_inputs=None, var_nd=None):
    """Determine each state/parameter variable's block length within a flat
    state vector of total length ``vector_size``.

    Every variable in ``vec_inputs`` (in order) gets exactly one of:
      - 1, if its name is listed in ``scalar_inputs``;
      - the explicit size configured in ``var_nd[name]``, if present;
      - otherwise, the common "regular field" size, inferred as the
        remaining length after every scalar/``var_nd`` block is subtracted,
        split evenly across the remaining ("regular") variables.

    This generalizes (and is the single source of truth for) the
    ``var_nd = {var: (1 if var in scalar_inputs else variable_size) for var
    in vec_inputs}`` convention two call sites already build independently
    (``applications/lorenz_model/lorenz_utils/mode3_runner.py`` and
    ``applications/issm_model/examples/basal_friction_variation/run_da_issm.py``)
    -- both only ever use it for scalar/uniform-field layouts; this
    function also honors a ``var_nd`` size other than 1, which the
    configuration key's own name implies should be possible even though no
    current application exercises it (see
    ``generate_initial_member_increment``'s docstring for how an
    unresolvable layout fails).

    Raises ``ValueError`` if the declared layout cannot exactly cover
    ``vector_size`` (rather than silently truncating, overlapping blocks,
    or padding with zeros) -- a mismatched size declaration is a
    configuration error, not something to guess through.
    """
    scalar_inputs = set(scalar_inputs or [])
    var_nd = dict(var_nd or {})
    vector_size = int(vector_size)

    block_sizes = []
    regular_slots = []
    for i, name in enumerate(vec_inputs):
        if name in scalar_inputs:
            block_sizes.append(1)
        elif name in var_nd:
            size = int(var_nd[name])
            if size <= 0:
                raise ValueError(
                    f"var_nd[{name!r}] must be a positive integer, got {size}"
                )
            block_sizes.append(size)
        else:
            block_sizes.append(None)
            regular_slots.append(i)

    explicit_total = sum(size for size in block_sizes if size is not None)
    if regular_slots:
        remaining = vector_size - explicit_total
        if remaining <= 0 or remaining % len(regular_slots):
            raise ValueError(
                "Cannot resolve a common field size for the "
                f"{len(regular_slots)} non-scalar, non-var_nd variable(s) in "
                f"vec_inputs={list(vec_inputs)!r}: vector_size={vector_size}, "
                f"explicit (scalar/var_nd) blocks already total {explicit_total}, "
                f"leaving {remaining}, which does not split evenly across "
                f"{len(regular_slots)} variable(s)."
            )
        common_hdim = remaining // len(regular_slots)
        for i in regular_slots:
            block_sizes[i] = common_hdim
    elif explicit_total != vector_size:
        raise ValueError(
            f"Declared scalar/var_nd blocks total {explicit_total}, but "
            f"vector_size={vector_size}: every variable in vec_inputs is "
            "scalar- or var_nd-sized, and the declared sizes do not exactly "
            "cover the state vector."
        )

    if sum(block_sizes) != vector_size:
        raise ValueError(
            f"Resolved block sizes {block_sizes} sum to {sum(block_sizes)}, "
            f"not vector_size={vector_size}."
        )
    return block_sizes


def generate_initial_member_increment(
    hdim, icesee_kwargs, ensemble_id, vector_size=None
):
    """Generate the canonical ICESEE initial-ensemble perturbation.

    Each state/parameter block gets an independent, member-keyed spatial
    field -- seeded by ``(base_seed, ensemble_id, variable_index + 1)`` via
    ``initialization_seed``, never by a shared/unreseeded generator -- and
    is scaled by the corresponding ``sig_Q``. This is the single source of
    truth for initial-ensemble generation, shared by every execution mode
    (0/1/2) so that an application's initial ensemble does not depend on
    which runner produced it. Historically this lived only inside mode 2's
    MPI-specific runner, which is why modes 0/1 could not use it and instead
    kept an older, separately broken implementation (one shared RNG object
    reused across every member with no reseed -- collapsing ensemble spread
    -- and the raw, ``length_scale``-only-scaled field added directly
    instead of being scaled by ``sig_Q``). Relocating it here (model/MPI-
    agnostic) lets every mode share one correct implementation instead of
    maintaining parallel ones.

    The ensemble-wide ``initial_spread_factor`` is deliberately *not*
    applied here; callers apply that factor later about the realized
    ensemble mean.

    Returns ``(scaled_increment, raw_increment)``, both 1-D arrays of
    length ``vector_size`` (default ``total_state_param_vars * hdim``):
    ``scaled_increment`` is what callers should add to the initial state:
    ``raw_increment`` (unscaled by ``sig_Q``) is exposed for callers that
    need the underlying field itself (e.g. diagnostics).

    Heterogeneous layouts: when ``icesee_kwargs['vec_inputs']`` is present
    together with a non-empty ``scalar_inputs`` and/or a ``var_nd`` dict,
    each variable's block length is resolved via
    ``resolve_variable_block_sizes`` instead of assuming every variable
    shares ``hdim`` -- see that function's docstring. A configured-scalar
    variable (e.g. flowline_1d's ``xg``, a single grounding-line position,
    or a bare state/observation count -- there is no plausible spatial
    decorrelation concept for a single number) draws its perturbation from
    the same generator as a regular field, degenerated to one node: this
    matches the one piece of historical intent available (the original,
    since-removed scalar branch in ``src/EnKF/_ensemble_initialization.py``
    also called the field generator at ``noise_dim=1`` for its scalar
    case), and empirically produces a well-defined, nonzero-variance,
    ``sig_Q``-scaled single value (see
    ``src/tests/test_generate_initial_member_increment.py``). When neither
    ``scalar_inputs`` nor ``var_nd`` is configured, behavior is completely
    unchanged from before (bit-identical): every variable uses ``hdim``.
    """
    configured_vars = int(icesee_kwargs["total_state_param_vars"])
    if vector_size is None:
        vector_size = configured_vars * int(hdim)
    vector_size = int(vector_size)

    vec_inputs = icesee_kwargs.get("vec_inputs")
    scalar_inputs = icesee_kwargs.get("scalar_inputs") or []
    var_nd = icesee_kwargs.get("var_nd")

    if vec_inputs and (scalar_inputs or var_nd):
        block_sizes = resolve_variable_block_sizes(
            list(vec_inputs)[:configured_vars],
            vector_size,
            scalar_inputs=scalar_inputs,
            var_nd=var_nd,
        )
    else:
        # Unchanged from before heterogeneous-layout support: every
        # variable shares the caller-supplied hdim.
        if vector_size % int(hdim):
            raise ValueError(
                "Initial state-vector size must be an integer number of variable "
                f"blocks: vector_size={vector_size}, hdim={hdim}."
            )
        block_sizes = [int(hdim)] * (vector_size // int(hdim))

    nvars = len(block_sizes)
    sig_q = list(icesee_kwargs.get("sig_Q", []))
    lx = float(icesee_kwargs.get("Lx", 1.0))
    ly = float(icesee_kwargs.get("Ly", 1.0))
    blocks = []
    raw_blocks = []
    for variable_index in range(nvars):
        this_block_size = int(block_sizes[variable_index])
        seed = initialization_seed(
            icesee_kwargs.get("base_seed", 42),
            ensemble_id,
            variable_index + 1,
        )
        field_kwargs = dict(icesee_kwargs)
        field_kwargs.update(
            {
                "ens_id": int(ensemble_id),
                "ii_sig": variable_index,
                "seed": seed,
                "rank_seed": seed,
                "rng": np.random.default_rng(seed),
                "Lx_dim": np.sqrt(lx * ly),
                "noise_dim": this_block_size,
                "num_vars": configured_vars,
            }
        )
        old_state = np.random.get_state()
        try:
            np.random.seed(seed)
            field = np.asarray(
                generate_enkf_field(**field_kwargs), dtype=np.float64
            ).reshape(-1)
        finally:
            np.random.set_state(old_state)
        if field.size != this_block_size:
            raise ValueError(
                "Initial random field has the wrong block size: "
                f"variable {variable_index} produced {field.size}, "
                f"expected {this_block_size}."
            )
        sigma = float(sig_q[variable_index]) if variable_index < len(sig_q) else 0.0
        raw_blocks.append(field)
        blocks.append(sigma * field)
    return np.concatenate(blocks), np.concatenate(raw_blocks)
