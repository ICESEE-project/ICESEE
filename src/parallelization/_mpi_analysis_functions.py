# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-03-06
# @description: computes the X5 matrix for the EnKF
#               - the new formulation is based on the paper by Geir Evensen: The Ensemble Kalman Filter: Theoretical Formulation And Practical Implementation
#               - this formulation supports our need for mpi parallelization and no need for localizations
# =============================================================================

import gc
import os
import copy
import h5py
import numpy as np
import bigmpi4py as BM
from scipy.stats import multivariate_normal, beta
from mpi4py import MPI

# seed the random number generator
# np.random.seed(0)
from ICESEE.src.utils.localization import (
    active_observation_std,
    apply_local_patches,
    compute_local_patches_X5,
    compute_X5_from_matrices,
    generated_observation_terms,
    stochastic_observation_terms,
    restore_frozen_analysis_vars,
)

from ICESEE.src.utils.inference_plugin import (
    apply_bed_update_gate_local,
)

from ICESEE.src.parallelization._parallel_i_o import (
    partition_rows,
    compute_HAprime_Eta_Dprime_partitioned,
    write_analysis_partitioned,
    parallel_write_ensemble_scattered
)


def stabilize_analysis_increments(
    analysis_vec,
    forecast_vec,
    global_rows,
    vec_inputs,
    hdim,
    icesee_kwargs,
):
    """Damp and cap analysis increments block-by-block.

    The operation is applied to increments, rather than absolute state values,
    so physically large but valid velocities are not clipped.  Defaults are a
    no-op; reviewer experiments opt in through the YAML configuration.
    """
    default_relaxation = float(
        icesee_kwargs.get(
            "analysis_relaxation_factor",
            1.0,
        )
    )
    relaxation_by_var = icesee_kwargs.get(
        "analysis_relaxation_factors",
        {},
    ) or {}
    increment_limits = icesee_kwargs.get(
        "analysis_increment_limits",
        {},
    ) or {}

    if not 0.0 < default_relaxation <= 1.0:
        raise ValueError("analysis_relaxation_factor must be in (0, 1]")

    relaxation_by_var = {
        str(key).lower(): float(value)
        for key, value in relaxation_by_var.items()
    }
    increment_limits = {
        str(key).lower(): float(value)
        for key, value in increment_limits.items()
    }

    stabilized = np.array(analysis_vec, copy=True)
    for block_index, key in enumerate(vec_inputs):
        start = block_index * hdim
        end = start + hdim
        local_mask = (global_rows >= start) & (global_rows < end)
        if not np.any(local_mask):
            continue

        key_lower = str(key).lower()
        relaxation = relaxation_by_var.get(key_lower, default_relaxation)
        if not 0.0 < relaxation <= 1.0:
            raise ValueError(
                f"analysis relaxation for {key!r} must be in (0, 1]"
            )

        increment = relaxation * (
            analysis_vec[local_mask, :] - forecast_vec[local_mask, :]
        )
        limit = increment_limits.get(key_lower)
        if limit is not None:
            if not np.isfinite(limit) or limit <= 0.0:
                raise ValueError(
                    f"analysis increment limit for {key!r} must be positive"
                )
            increment = np.clip(increment, -limit, limit)

        stabilized[local_mask, :] = forecast_vec[local_mask, :] + increment

    if not np.all(np.isfinite(stabilized)):
        raise FloatingPointError("non-finite values produced by EnKF analysis")
    return stabilized


def apply_analysis_controls_local(
    analysis_vec,
    forecast_vec,
    global_rows,
    icesee_kwargs,
):
    """Apply the row-local post-analysis policy shared by modes 1 and 2.

    Keeping these operations in one function is part of the execution-mode
    parity contract: a full-parallel run must not silently use different
    relaxation, inflation, frozen-variable, or bed-update semantics merely
    because its state matrix is streamed from disk.
    """
    vec_inputs = list(icesee_kwargs.get("vec_inputs", []))
    if not vec_inputs:
        raise ValueError("vec_inputs is empty during EnKF analysis")

    nd = int(icesee_kwargs["nd"])
    # The current analysis controls operate on the common ICESEE block layout
    # (one mesh-sized block per configured state/parameter).  Reject an
    # incompatible layout rather than applying controls to the wrong rows.
    if nd % len(vec_inputs) != 0:
        raise ValueError(
            "analysis controls require equal-sized state/parameter blocks: "
            f"nd={nd}, blocks={len(vec_inputs)}"
        )
    hdim = nd // len(vec_inputs)

    controlled = stabilize_analysis_increments(
        analysis_vec=analysis_vec,
        forecast_vec=forecast_vec,
        global_rows=global_rows,
        vec_inputs=vec_inputs,
        hdim=hdim,
        icesee_kwargs=icesee_kwargs,
    )

    state_inflation = icesee_kwargs.get(
        "state_inflation_factor", icesee_kwargs.get("inflation_factor", 1.0)
    )
    param_inflation = icesee_kwargs.get(
        "param_inflation_factor", icesee_kwargs.get("inflation_factor", 1.0)
    )
    bed_inflation = icesee_kwargs.get("bed_inflation_factor", param_inflation)

    inflation_vec = np.full(global_rows.size, float(param_inflation))
    num_state_vars = int(icesee_kwargs.get("num_state_vars", 0))
    for block_index, key in enumerate(vec_inputs):
        start = block_index * hdim
        end = start + hdim
        local_mask = (global_rows >= start) & (global_rows < end)
        if not np.any(local_mask):
            continue
        key_lower = str(key).lower()
        if block_index < num_state_vars:
            inflation_vec[local_mask] = float(state_inflation)
        elif key_lower in {
            "bed", "bedrock", "bedtopography", "bed_topography",
            "bed_elevation",
        }:
            inflation_vec[local_mask] = float(bed_inflation)

    local_mean = np.mean(controlled, axis=1, keepdims=True)
    controlled = local_mean + inflation_vec[:, None] * (controlled - local_mean)
    controlled = restore_frozen_analysis_vars(
        controlled,
        forecast_vec,
        global_rows,
        vec_inputs,
        hdim,
        icesee_kwargs.get("frozen_analysis_vars", []),
    )
    controlled = apply_bed_update_gate_local(
        analysis_vec=controlled,
        forecast_vec=forecast_vec,
        global_rows=global_rows,
        vec_inputs=vec_inputs,
        hdim=hdim,
        icesee_kwargs=icesee_kwargs,
    )
    if not np.all(np.isfinite(controlled)):
        raise FloatingPointError("non-finite values produced by analysis controls")
    return controlled


# ============================ EnKF functions ============================
def EnKF_X5(k_obs, ensemble_vec, Nens, hu_obs, icesee_kwargs, UtilsFunctions):
    comm_world = icesee_kwargs.get("comm_world")

    if icesee_kwargs.get("partitioned_io_flag", False):
        if str(icesee_kwargs.get("bed_update_domain", "all")).lower() != "all":
            raise ValueError(
                "a restricted bed_update_domain currently requires "
                "partitioned_io_flag=False"
            )
        # ============================================== NEW BRANCH
        input_file = icesee_kwargs.get("_forecast_h5_path")
        timestep_forecast = icesee_kwargs.get("k")

        icesee_kwargs["hu_obs_loaded"] = hu_obs
        icesee_kwargs["km"] = k_obs

        U = UtilsFunctions(icesee_kwargs=icesee_kwargs, ensemble=None)
        obs_indices = U.JObs_indices(icesee_kwargs["nd"])

        y_full = hu_obs[:, k_obs].copy()
        y_full[np.isnan(y_full)] = 0.0
        d = y_full[obs_indices]

        error_mode = str(
            icesee_kwargs.get("enkf_observation_error_mode", "stochastic_R")
        ).lower()
        sigma = (
            active_observation_std(icesee_kwargs, k_obs, obs_indices)
            if error_mode == "stochastic_r"
            else None
        )
        obs_seed = int(icesee_kwargs.get("base_seed", 42)) + 1000003 * (
            int(k_obs) + 1
        )

        HAprime, Eta, Dprime = compute_HAprime_Eta_Dprime_partitioned(
            input_file,
            timestep_forecast,
            obs_indices,
            d,
            Nens,
            comm_world,
            sigma=sigma,
            seed=obs_seed,
            error_mode=error_mode,
        )

        X5 = compute_X5_from_matrices(HAprime, Eta, Dprime, Nens)

        local_patches = None
        if icesee_kwargs.get("local_analysis", False):
            local_patches = compute_local_patches_X5(
                vec_inputs=icesee_kwargs.get("vec_inputs", []),
                hdim=icesee_kwargs["nd"] // icesee_kwargs["total_state_param_vars"],
                HAprime=HAprime, Eta=Eta, Dprime=Dprime,
                Nens=Nens, obs_indices=obs_indices, icesee_kwargs=icesee_kwargs
            )

        icesee_kwargs["_analysis_h5_path"] = input_file
        icesee_kwargs["_analysis_timestep_forecast"] = timestep_forecast

        return X5, local_patches
        # ============================================== END NEW

    # ---- EXISTING in-memory path, unchanged below ----
    icesee_kwargs["hu_obs_loaded"] = hu_obs
    icesee_kwargs["km"] = k_obs

    U = UtilsFunctions(icesee_kwargs=icesee_kwargs, ensemble=ensemble_vec)
    obs_indices = U.JObs_indices(ensemble_vec.shape[0])

    y_full = hu_obs[:, k_obs].copy()
    y_full[np.isnan(y_full)] = 0.0
    d = y_full[obs_indices]

    error_mode = str(
        icesee_kwargs.get("enkf_observation_error_mode", "stochastic_R")
    ).lower()
    HA = ensemble_vec[obs_indices, :]
    if error_mode == "stochastic_r":
        sigma = active_observation_std(icesee_kwargs, k_obs, obs_indices)
        obs_seed = int(icesee_kwargs.get("base_seed", 42)) + 1000003 * (
            int(k_obs) + 1
        )
        HAprime, Eta, Dprime = stochastic_observation_terms(
            HA, d, sigma, obs_seed
        )
    elif error_mode == "legacy_prior_anomalies":
        HAprime = HA - np.mean(HA, axis=1, keepdims=True)
        Eta = HAprime.copy()
        Dprime = d.reshape(-1, 1) - HA
    elif error_mode == "generated_r":
        HAprime, Eta, Dprime = generated_observation_terms(
            HA, d, icesee_kwargs, obs_indices, ensemble_vec.shape[0], k_obs
        )
    else:
        raise ValueError(
            "enkf_observation_error_mode must be 'stochastic_R', "
            "'legacy_prior_anomalies', or 'generated_R'"
        )

    m = d.shape[0]
    nrmin = min(m, Nens)

    HAprime_eta = HAprime + Eta
    U_svd, sig, _ = np.linalg.svd(HAprime_eta, full_matrices=False)
    sig = sig ** 2

    sigsum = np.sum(sig[:nrmin])
    sigsum1 = 0.0
    for i in range(nrmin):
        if sigsum1 / sigsum < 0.999:
            sigsum1 += sig[i]
            sig[i] = 1.0 / sig[i]
        else:
            sig[i:nrmin] = 0.0
            break

    X1 = np.empty((nrmin, m))
    for j in range(m):
        for i in range(nrmin):
            X1[i, j] = sig[i] * U_svd[j, i]

    X2 = np.dot(X1, Dprime)
    X3 = np.dot(U_svd, X2)
    X4 = np.dot(HAprime.T, X3)

    X5 = X4 + np.eye(Nens)
    if np.sum(X5, axis=0).all() != 1.0:
        print(f"[ICESEE] Sum of each X5 column is not 1.0: {np.sum(X5, axis=0)}")

    local_patches = None
    if icesee_kwargs.get("local_analysis", False):
        local_patches = compute_local_patches_X5(
            vec_inputs=icesee_kwargs.get("vec_inputs", []),
            hdim=ensemble_vec.shape[0] // icesee_kwargs["total_state_param_vars"],
            HAprime=HAprime, Eta=Eta, Dprime=Dprime,
            Nens=Nens, obs_indices=obs_indices, icesee_kwargs=icesee_kwargs
        )

    del X2, X3, X4, U_svd, HAprime
    gc.collect()

    return X5, local_patches

def analysis_enkf_update(
    k, ens_mean, ensemble_vec, shape_ens, X5,
    time_analysis_mean_generation, time_analysis_file_writing,
    local_patches, UtilsFunctions, icesee_kwargs, smb_scale
):
    import gc
    import numpy as np
    from mpi4py import MPI

    comm_world = icesee_kwargs.get("comm_world")

    if icesee_kwargs.get("partitioned_io_flag", False):
        # ============================================== NEW BRANCH
        h5_path = icesee_kwargs.get("_analysis_h5_path")
        timestep_forecast = icesee_kwargs.get("_analysis_timestep_forecast")

        t0 = MPI.Wtime()
        write_analysis_partitioned(k, X5, local_patches, h5_path, timestep_forecast, icesee_kwargs, comm_world)
        time_analysis_file_writing += MPI.Wtime() - t0

        return time_analysis_mean_generation, time_analysis_file_writing
        # ============================================== END NEW

    # ---- EXISTING Scatterv-based path, unchanged below ----
    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()

    X5 = BM.bcast(X5, comm=comm_world)
    # EnKF_X5 constructs grouped local-analysis transforms on rank zero in
    # the in-memory mode.  Every rank owns a different slab of state rows, so
    # each rank must receive the same patch registry before applying its local
    # transforms.  Broadcasting only the global X5 silently made non-root
    # ranks fall back to the global analysis for localized variables.
    local_patches = comm_world.bcast(local_patches, root=0)
    time_analysis_mean_generation = BM.bcast(time_analysis_mean_generation, comm=comm_world)
    icesee_kwargs["X5"] = X5

    if rank_world != 0:
        ensemble_vec = None

    nd_total, nens = shape_ens

    counts_rows = np.array(
        [nd_total // size_world + (1 if r < nd_total % size_world else 0)
         for r in range(size_world)], dtype=np.int32
    )
    offsets_rows = np.insert(np.cumsum(counts_rows), 0, 0)[:-1].astype(np.int32)
    local_nd = int(counts_rows[rank_world])
    row_offset = int(offsets_rows[rank_world])
    global_rows = row_offset + np.arange(local_nd)

    counts = (counts_rows * nens).astype(np.int32)
    displs = (offsets_rows * nens).astype(np.int32)

    scatter_ensemble = np.empty((local_nd, nens), dtype=np.float64)
    if rank_world == 0 and np.shape(ensemble_vec) != (nd_total, nens):
        raise ValueError(
            "analysis shape metadata does not match the active ensemble: "
            f"shape_ens={(nd_total, nens)}, ensemble={np.shape(ensemble_vec)}"
        )
    sendbuf = np.ascontiguousarray(ensemble_vec) if rank_world == 0 else None
    comm_world.Scatterv([sendbuf, counts, displs, MPI.DOUBLE], scatter_ensemble, root=0)

    analysis_vec = np.dot(scatter_ensemble, X5)

    if icesee_kwargs.get("local_analysis", False):
        analysis_vec = apply_local_patches(analysis_vec, scatter_ensemble, global_rows, local_patches)

    t0 = MPI.Wtime()
    analysis_vec = apply_analysis_controls_local(
        analysis_vec=analysis_vec,
        forecast_vec=scatter_ensemble,
        global_rows=global_rows,
        icesee_kwargs=icesee_kwargs,
    )
    time_analysis_mean_generation += MPI.Wtime() - t0

    t0 = MPI.Wtime()
    parallel_write_ensemble_scattered(
        k + 1,
        ens_mean,
        analysis_vec,
        comm_world,
        icesee_kwargs,
        forecast_chunk=scatter_ensemble,
    )
    time_analysis_file_writing += MPI.Wtime() - t0

    del scatter_ensemble, analysis_vec
    gc.collect()

    return time_analysis_mean_generation, time_analysis_file_writing


# ============================ EnKF functions ============================

# ============================ DEnKF functions ============================
def DEnKF_X5(k,ensemble_vec, Cov_obs, Nens, d, icesee_kwargs,UtilsFunctions):
    """
    Function to compute the X5 matrix for the DEnKF
        - ensemble_vec: ensemble matrix of size (ndxNens)
        - Cov_obs: observation covariance matrix
        - Nens: ensemble size
        - d: observation vector
    """
    comm_world = icesee_kwargs.get("comm_world")
    H = UtilsFunctions(icesee_kwargs=icesee_kwargs, ensemble=ensemble_vec).JObs_fun(ensemble_vec.shape[0]) # mxNens, observation operator

    # -- get ensemble pertubations
    ensemble_perturbations = ensemble_vec - np.mean(ensemble_vec, axis=1).reshape(-1,1)

    # ----parallelize this step
    A_anomaly = np.zeros_like(ensemble_vec) # mxNens, ensemble pertubations
    Eta = np.dot(H, ensemble_perturbations) # mxNens, ensemble pertubations
    # D   = np.zeros_like(Eta) # mxNens #virtual observations
    ens_mean = np.mean(ensemble_vec, axis=1)
    # Eta = np.zeros((d.shape[0], Nens)) # mxNens
    HA  = np.zeros_like(Eta)
    ha = np.zeros_like(Eta)
    for ens in range(Nens):
        A_anomaly[:,ens] = ensemble_vec[:,ens] - ens_mean
        # D[:,ens] = d + Eta[:,ens]
        # HA[:,ens] = np.dot(H, ensemble_vec[:,ens])
        HA[:,ens] = np.dot(H, A_anomaly[:,ens])
        # ha[:,ens] = np.dot(H, ensemble_vec[:,ens])
    # # ---------------------------------------

    # # --- compute the innovations D` = D-HA
    # Dprime = D - HA # mxNens

    # --- compute HAbar
    # HAbar = np.mean(HA, axis=1) # mx1
    # --- compute HAprime
    # HAprime = HA - HAbar.reshape(-1,1) # mxNens (requires H to be linear)

    # Aprime = ensemble_vec@(np.eye(Nens) - one_N) # mxNens
    one_N = np.ones((Nens,Nens))/Nens
    HAprime=HA@(np.eye(Nens) - one_N) # mxNens

    # get the min(m,Nens)
    m_obs = d.shape[0]
    nrmin = min(m_obs, Nens)

    # --- compute HA' + eta
    HAprime_eta = HAprime + Eta

    # --- compute the SVD of HA' + eta
    U, sig, _ = np.linalg.svd(HAprime_eta, full_matrices=False)

    # --- convert s to eigenvalues
    sig = sig**2
    # for i in range(nrmin):
    #     sig[i] = sig[i]**2

    # ---compute the number of significant eigenvalues
    sigsum = np.sum(sig[:nrmin])  # Compute total sum of the first `nrmin` eigenvalues
    sigsum1 = 0.0
    nrsigma = 0

    for i in range(nrmin):
        if sigsum1 / sigsum < 0.999:
            nrsigma += 1
            sigsum1 += sig[i]
            sig[i] = 1.0 / sig[i]  # Inverse of eigenvalue
        else:
            sig[i:nrmin] = 0.0  # Set remaining eigenvalues to 0
            break  # Exit the loop

    # compute X1 = sig*UT #Nens x m_obs
    X1 = np.empty((nrmin, m_obs))
    for j in range(m_obs):
        for i in range(nrmin):
            X1[i,j] =sig[i]*U[j,i]

    # compute X2 = X1*Dprime # Nens x Nens
    # X2 = np.dot(X1, Dprime)
    X2 = np.dot(X1, HA) # Nens x Nens  #TODO  or np.dot(X1, HA)???
    # del Cov_obs, sig, X1, Dprime; gc.collect()

    # --get wprime
    wprime = d - np.dot(H, ens_mean)
    X2prime = np.dot(X1, wprime) # Nens x Nens

    # print(f"[ICESEE] Rank: {rank_world} X2 shape: {X2.shape}")
    #  compute X3 = U*X2 # m_obs x Nens
    X3 = np.dot(U, X2)
    X3prime = np.dot(U, X2prime) # m_obs x Nens

    # print(f"[ICESEE] Rank: {rank_world} X3 shape: {X3.shape}")
    # compute X4 = (HAprime.T)*X3 # Nens x Nens
    X4 = np.dot(HAprime.T, X3)
    X4prime = np.dot(HAprime.T, X3prime) # Nens x Nens
    del X2, X3, U, HAprime; gc.collect()

    # print(f"[ICESEE] Rank: {rank_world} X4 shape: {X4.shape}")
    # compute X5 = X4 + I
    # X5 = X4 + np.eye(Nens)
    X5 = 0.5*(2*np.eye(Nens) + np.dot(one_N, X4) - X4) #TODO check this
    # X5 = 0.5*(2*np.eye(Nens) - X4)
    X5prime = one_N + np.dot((np.eye(Nens) - one_N),X4prime) #TODO check this
    # X5prime = (one_N - X4prime)
    X5 =  (np.eye(Nens) - (0.5*(np.dot(np.eye(Nens) - one_N, X4)))) + one_N + np.dot((np.eye(Nens) - one_N),X4prime)
    # X5 = 0.5*(2*np.eye(Nens) - X4)
    # sum of each column of X5 should be 1
    if np.sum(X5, axis=0).all() != 1.0:
        print(f"[ICESEE] Sum of each X5 column is not 1.0: {np.sum(X5, axis=0)}")
    # print(f"[ICESEE] Rank: {comm_world.Get_rank()} X5 sum: {np.sum(X5, axis=0)}")
    del X4; gc.collect()

    # ===local computation
    if icesee_kwargs.get("local_analysis",False):
        analysis_vec_ij = np.empty_like(ensemble_vec)
        AssertionError("Local analysis is not implemented yet for DEnKF")
    else:
        analysis_vec_ij = None


    return X5, X5prime

def analysis_Denkf_update(k,ens_mean,ensemble_vec, shape_ens, X5, UtilsFunctions,icesee_kwargs,smb_scale):
    """
    Function to perform the analysis update using the EnKF
        - broadcast X5 to all processors
        - initialize an empty ensemble vector for the rest of the processors
        - scatter ensemble_vec to all processors
        - do the ensemble analysis update: A_j = Fj*X5
        - gather from all processors
    """


    if icesee_kwargs.get("local_analysis",False):
        pass
    else:
        comm_world = icesee_kwargs.get("comm_world")
        # get the rank and size of the world communicator
        rank_world = comm_world.Get_rank()
        # broadcast X5 to all processors
        X5 = BM.bcast(X5, comm=comm_world)
        # ens_mean = BM.bcast(ens_mean, comm=comm_world)
        # X5_diff = BM.bcast(X5_diff, comm=comm_world)

        # initialize the an empty ensemble vector for the rest of the processors
        if rank_world != 0:
            ensemble_vec = np.empty(shape_ens, dtype=np.float64)

        # --- scatter ensemble_vec to all processors ---
        scatter_ensemble = BM.scatter(ensemble_vec, comm_world)
        # -* instead of using scattter from root, if the ensemble vec doesn't fit in memory then
        # with h5py.File("icesee_ensemble_data.h5", 'r', driver='mpio', comm=comm_world) as f:
        #     scatter_ensemble = f['ensemble']
        #     total_rows = scatter_ensemble.shape[0]

        #     # calculate rows per rank
        #     rows_per_rank = total_rows // comm_world.Get_size()
        #     # remainder = total_rows % comm_world.Get_size()
        #     start_row = rank_world * rows_per_rank
        #     end_row = start_row + rows_per_rank if rank_world != comm_world.Get_size()-1 else total_rows

        #     # Each rank reads its chunk from the dataset
        #     scatter_ensemble = scatter_ensemble[start_row:end_row, :, k]
        # do the ensemble analysis update: A_j = Fj*X5
        analysis_vec = np.dot(scatter_ensemble, X5)
        # ens_mean_ = np.dot(scatter_ensemble, X5prime)

        # print(f"[ICESEE] Rank: {rank_world} analysis_vec shape: {analysis_vec.shape}, ens_mean shape: {ens_mean.shape}")

        # comm_world.Barrier()
        # analysis_vec = analysis_vec + ens_mean

        # ens_mean = np.mean(analysis_vec, axis=1)

        ndim = analysis_vec.shape[0] // icesee_kwargs["total_state_param_vars"]
        state_block_size = ndim*icesee_kwargs["num_state_vars"]

        # analysis_vec[state_block_size:,:] /= 10
        # analysis_vec[state_block_size:,:] *= (smb_scale)  # Scale SMB after analysis
        # icesee_kwargs['inflation_factor'] = 1.1
        # analysis_vec = UtilsFunctions(icesee_kwargs,  analysis_vec).inflate_ensemble(in_place=True)
        # ---> multiplicative inflation
        mean_params = np.mean(analysis_vec[state_block_size:,:], axis=1)
        mean_vars = np.mean(analysis_vec[:ndim,:], axis=1)
        #  compute parturbations
        pertubations = analysis_vec[state_block_size:,:] - mean_params.reshape(-1,1)
        pertubations_vars = analysis_vec[:ndim,:] - mean_vars.reshape(-1,1)
        # apply the inflation factor
        inflated_pertubations = pertubations * icesee_kwargs['inflation_factor']
        # inflated_pertubations_vars = pertubations_vars * icesee_kwargs['inflation_factor']

        # update the analysis vector
        analysis_vec[state_block_size:,:] = mean_params.reshape(-1,1) + inflated_pertubations
        # analysis_vec[:ndim,:] = mean_vars.reshape(-1,1) + inflated_pertubations_vars


        # check for negative thicknes and set to 1e-3 if vec_input contains h
        for i, var in enumerate(icesee_kwargs.get("vec_inputs",[])):
            if var == "h" or var == "thickness" or var == "ice_thickness" or var == "Thickness":
                start = i * ndim
                end = start + ndim
                analysis_vec[start:end, :] = np.maximum(analysis_vec[start:end, :], 1e-2)

        # # ISSM *------
        # di = 0.8930
        # rho_ice = 917.0
        # rho_sw = 1028.0
        # ocean_levelset = analysis_vec[:ndim,:] + analysis_vec[state_block_size:ndim,:]/di
        # # Floating ice (ocean_levelset < 0) find the indices
        # pos = np.where(ocean_levelset < 0)
        # thickness_floating = analysis_vec[:ndim,:]
        # surface = analysis_vec[ndim:2*ndim,:]
        # surface[pos] = thickness_floating[pos]* (rho_sw - rho_ice)/rho_sw
        # analysis_vec[ndim:2*ndim,:] = surface

        # # read base data from h5file and compute the mean base from all ensembles




        # *---------

        # dynamical model for parameters: from https://doi.org/10.1002/qj.3257
        # obs_index = icesee_kwargs.get("obs_index")
        # # #  check if k equals to the first observation index
        # # print(f"[ICESEE] Rank: {rank_world} km: {km} obs_index: {obs_index}")
        # if  (k+1 == obs_index[0]):
        # #     print(f"[ICESEE] [Debug] Rank: {rank_world} k: {km} obs_index: {obs_index}")
        #     params_analysis_0 = analysis_vec[state_block_size:, :]

        # # size of parameters
        # param_size = analysis_vec.shape[0] - state_block_size
        # alpha = np.ones(param_size)*2.0
        # beta_param = alpha
        # def compute_f_params(alpha, beta_param):
        #     mean_x = alpha/(alpha+beta_param)
        #     a = 1.0
        #     b = -a*mean_x
        #     return a,b

        # def update_theta(alpha, beta_param):
        #     # theta_f_t = np.zeros_like(theta_prev)
        #     f_x_ti = np.zeros((param_size,analysis_vec.shape[1]))
        #     for i in range(analysis_vec.shape[1]):
        #         a,b = compute_f_params(alpha[i], beta_param[i])
        #         x_ti = beta.rvs(alpha[i], beta_param[i])

        #         f_x_ti[:,i] = a*x_ti + b

        #         # theta_f_t[:,i] = theta_prev[:,i] + f_x_ti
        #     # return theta_f_t
        #     return f_x_ti

        # analysis_vec[state_block_size:,:] = params_analysis_0 +  update_theta(alpha, beta_param)

        # # X = beta.rvs(alpha, beta_param,param_size)
        # # linear_bijective_function = lambda x,a: 2*a*(x - 0.5) #zero mean
        # # analysis_vec[state_block_size:,:] = params_analysis_0 + linear_bijective_function(X,a=0.1)

        # params_analysis_0 = analysis_vec[state_block_size:, :]


        # gather from all processors
        # ensemble_vec = BM.allgather(analysis_vec, comm_world)
        parallel_write_ensemble_scattered(
            k + 1,
            ens_mean,
            analysis_vec,
            comm_world,
            icesee_kwargs,
            forecast_chunk=scatter_ensemble,
        )

        # clean the memory
        del scatter_ensemble, analysis_vec; gc.collect()


# ============================ EnSRF functions ============================


# ============================ EnTKF functions ============================


# ============================ Other functions ============================
