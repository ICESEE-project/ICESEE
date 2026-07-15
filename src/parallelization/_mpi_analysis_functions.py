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
from ICESEE.src.utils.localization import apply_local_patches, compute_local_patches_X5, compute_X5_from_matrices, apply_smb_physics_correction

from ICESEE.src.parallelization._parallel_i_o import (
    partition_rows,
    compute_HAprime_Eta_Dprime_partitioned,
    write_analysis_partitioned,
    parallel_write_ensemble_scattered
)

# ============================ EnKF functions ============================ 
def EnKF_X5(k_obs, ensemble_vec, Nens, hu_obs, model_kwargs, UtilsFunctions):
    params = model_kwargs.get("params")
    comm_world = model_kwargs.get("comm_world")

    if model_kwargs.get("partitioned_io_flag", False):
        # ============================================== NEW BRANCH
        input_file = model_kwargs.get("_forecast_h5_path")
        timestep_forecast = model_kwargs.get("k")

        model_kwargs["hu_obs_loaded"] = hu_obs
        model_kwargs["km"] = k_obs

        U = UtilsFunctions(params=params, model_kwargs=model_kwargs, ensemble=None)
        obs_indices = U.JObs_indices(params["nd"])

        y_full = hu_obs[:, k_obs].copy()
        y_full[np.isnan(y_full)] = 0.0
        d = y_full[obs_indices]

        HAprime, Eta, Dprime = compute_HAprime_Eta_Dprime_partitioned(
            input_file, timestep_forecast, obs_indices, d, Nens, comm_world
        )

        X5 = compute_X5_from_matrices(HAprime, Eta, Dprime, Nens)

        local_patches = None
        if model_kwargs.get("local_analysis", False):
            local_patches = compute_local_patches_X5(
                vec_inputs=model_kwargs.get("vec_inputs", []),
                hdim=params["nd"] // params["total_state_param_vars"],
                HAprime=HAprime, Eta=Eta, Dprime=Dprime,
                Nens=Nens, obs_indices=obs_indices, model_kwargs=model_kwargs
            )

        model_kwargs["_analysis_h5_path"] = input_file
        model_kwargs["_analysis_timestep_forecast"] = timestep_forecast

        return X5, local_patches
        # ============================================== END NEW

    # ---- EXISTING in-memory path, unchanged below ----
    generate_enkf_field = model_kwargs.get("generate_enkf_field", False)
    model_kwargs["hu_obs_loaded"] = hu_obs
    model_kwargs["km"] = k_obs

    U = UtilsFunctions(params=params, model_kwargs=model_kwargs, ensemble=ensemble_vec)
    obs_indices = U.JObs_indices(ensemble_vec.shape[0])

    y_full = hu_obs[:, k_obs].copy()
    y_full[np.isnan(y_full)] = 0.0
    d = y_full[obs_indices]

    use_ensemble_pertubations = model_kwargs.get("use_ensemble_pertubations", True)
    ensemble_mean = np.mean(ensemble_vec, axis=1).reshape(-1, 1)

    if use_ensemble_pertubations:
        ensemble_perturbations = ensemble_vec - ensemble_mean
        Eta = ensemble_perturbations[obs_indices, :] 
    else:
        if model_kwargs["joint_estimation"] or params["localization_flag"]:
            hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
        else:
            hdim = ensemble_vec.shape[0] // params["num_state_vars"]

        Lx, Ly = model_kwargs.get("Lx"), model_kwargs.get("Ly")

        if model_kwargs.get("inversion_flag", False):
            friction_idx = int(model_kwargs.get("friction_idx", -1))
            sig_obs_filtered = [params["sig_obs"][i] for i in range(len(params["sig_obs"])) if i != friction_idx]

        _eta = []
        for ens in range(Nens):
            noise_all = []
            for ii, sig in enumerate(sig_obs_filtered if model_kwargs.get("inversion_flag", False) else params["sig_obs"]):
                model_kwargs.update({"ii_sig": ii, "Lx_dim": np.sqrt(Lx * Ly), "noise_dim": hdim, "num_vars": params["total_state_param_vars"]})
                W = generate_enkf_field(**model_kwargs)
                noise_all.append(sig * W)
            noise_ = np.concatenate(noise_all, axis=0)
            _eta.append(noise_)
        _eta = np.array(_eta).T
        _eta -= np.mean(_eta, axis=1).reshape(-1, 1)
        Eta = _eta[obs_indices, :]
        
    Dprime = d.reshape(-1, 1) - ensemble_vec[obs_indices, :]

    one_N = np.ones((Nens, Nens)) / Nens
    Aprime = np.dot(ensemble_vec, (np.eye(Nens) - one_N))
    HAprime = Aprime[obs_indices, :]

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
    if model_kwargs.get("local_analysis", False):
        local_patches = compute_local_patches_X5(
            vec_inputs=model_kwargs.get("vec_inputs", []),
            hdim=ensemble_vec.shape[0] // params["total_state_param_vars"],
            HAprime=HAprime, Eta=Eta, Dprime=Dprime,
            Nens=Nens, obs_indices=obs_indices, model_kwargs=model_kwargs
        )

    del X2, X3, X4, U_svd, HAprime
    gc.collect()

    return X5, local_patches

def analysis_enkf_update(
    k, ens_mean, ensemble_vec, shape_ens, X5,
    time_analysis_mean_generation, time_analysis_file_writing,
    local_patches, UtilsFunctions, model_kwargs, smb_scale
):
    import gc
    import numpy as np
    from mpi4py import MPI

    params = model_kwargs.get("params")
    comm_world = model_kwargs.get("comm_world")

    if model_kwargs.get("partitioned_io_flag", False):
        # ============================================== NEW BRANCH
        h5_path = model_kwargs.get("_analysis_h5_path")
        timestep_forecast = model_kwargs.get("_analysis_timestep_forecast")

        t0 = MPI.Wtime()
        write_analysis_partitioned(k, X5, local_patches, h5_path, timestep_forecast, model_kwargs, comm_world)
        time_analysis_file_writing += MPI.Wtime() - t0

        return time_analysis_mean_generation, time_analysis_file_writing
        # ============================================== END NEW

    # ---- EXISTING Scatterv-based path, unchanged below ----
    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()

    X5 = BM.bcast(X5, comm=comm_world)
    time_analysis_mean_generation = BM.bcast(time_analysis_mean_generation, comm=comm_world)
    model_kwargs["X5"] = X5

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
    sendbuf = np.ascontiguousarray(ensemble_vec) if rank_world == 0 else None
    comm_world.Scatterv([sendbuf, counts, displs, MPI.DOUBLE], scatter_ensemble, root=0)

    analysis_vec = np.dot(scatter_ensemble, X5)

    if model_kwargs.get("local_analysis", False):
        analysis_vec = apply_local_patches(analysis_vec, scatter_ensemble, global_rows, local_patches)

    t0 = MPI.Wtime()
    state_inflation = params.get("state_inflation_factor", params.get("inflation_factor", 1.0))
    param_inflation = params.get("param_inflation_factor", params.get("inflation_factor", 1.0))
    bed_inflation = params.get("bed_inflation_factor", param_inflation)

    vec_inputs = model_kwargs.get("vec_inputs", [])
    nblocks = len(vec_inputs)
    hdim = model_kwargs.get("nd", params["nd"]) // nblocks

    inflation_vec = np.ones(local_nd) * param_inflation
    for ii, key in enumerate(vec_inputs):
        start, end = ii * hdim, ii * hdim + hdim
        local_mask = (global_rows >= start) & (global_rows < end)
        if not np.any(local_mask):
            continue
        key_l = key.lower()
        if ii < params["num_state_vars"]:
            inflation_vec[local_mask] = state_inflation
        elif key_l in ["bed", "bedrock", "bedtopography", "bed_topography"]:
            inflation_vec[local_mask] = bed_inflation
        else:
            inflation_vec[local_mask] = param_inflation

    local_mean = np.mean(analysis_vec, axis=1, keepdims=True)
    local_pert = analysis_vec - local_mean
    analysis_vec = local_mean + inflation_vec[:, None] * local_pert
    time_analysis_mean_generation += MPI.Wtime() - t0

    km = model_kwargs.get("km", None)
    bed_snap_cols = set(model_kwargs.get("bed_snap_cols", []))
    if km is not None and int(km) not in bed_snap_cols:
        for ii, key in enumerate(vec_inputs):
            if key.lower() not in ["bed", "bedrock", "bedtopography", "bed_topography"]:
                continue
            start, end = ii * hdim, ii * hdim + hdim
            local_mask = (global_rows >= start) & (global_rows < end)
            if np.any(local_mask):
                analysis_vec[local_mask, :] = scatter_ensemble[local_mask, :]

    t0 = MPI.Wtime()
    parallel_write_ensemble_scattered(k + 1, ens_mean, params, analysis_vec, comm_world, model_kwargs)
    time_analysis_file_writing += MPI.Wtime() - t0

    del scatter_ensemble, analysis_vec
    gc.collect()

    return time_analysis_mean_generation, time_analysis_file_writing


# ============================ EnKF functions ============================

# ============================ DEnKF functions ============================
def DEnKF_X5(k,ensemble_vec, Cov_obs, Nens, d, model_kwargs,UtilsFunctions):
    """
    Function to compute the X5 matrix for the DEnKF
        - ensemble_vec: ensemble matrix of size (ndxNens)
        - Cov_obs: observation covariance matrix
        - Nens: ensemble size
        - d: observation vector
    """
    params = model_kwargs.get("params")
    comm_world = model_kwargs.get("comm_world")
    H = UtilsFunctions(params =params, model_kwargs=model_kwargs,ensemble= ensemble_vec).JObs_fun(ensemble_vec.shape[0]) # mxNens, observation operator

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
    if model_kwargs.get("local_analysis",False):
        analysis_vec_ij = np.empty_like(ensemble_vec)
        AssertionError("Local analysis is not implemented yet for DEnKF")
    else:
        analysis_vec_ij = None
        

    return X5, X5prime

def analysis_Denkf_update(k,ens_mean,ensemble_vec, shape_ens, X5, UtilsFunctions,model_kwargs,smb_scale):
    """
    Function to perform the analysis update using the EnKF
        - broadcast X5 to all processors
        - initialize an empty ensemble vector for the rest of the processors
        - scatter ensemble_vec to all processors
        - do the ensemble analysis update: A_j = Fj*X5
        - gather from all processors
    """
    
    
    if model_kwargs.get("local_analysis",False):
        pass
    else:
        params = model_kwargs.get("params")
        comm_world = model_kwargs.get("comm_world")
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

        ndim = analysis_vec.shape[0] // params["total_state_param_vars"]
        state_block_size = ndim*params["num_state_vars"]

        # analysis_vec[state_block_size:,:] /= 10
        # analysis_vec[state_block_size:,:] *= (smb_scale)  # Scale SMB after analysis
        # params['inflation_factor'] = 1.1
        # analysis_vec = UtilsFunctions(params,  analysis_vec).inflate_ensemble(in_place=True)
        # ---> multiplicative inflation
        mean_params = np.mean(analysis_vec[state_block_size:,:], axis=1)
        mean_vars = np.mean(analysis_vec[:ndim,:], axis=1)
        #  compute parturbations
        pertubations = analysis_vec[state_block_size:,:] - mean_params.reshape(-1,1)
        pertubations_vars = analysis_vec[:ndim,:] - mean_vars.reshape(-1,1)
        # apply the inflation factor
        inflated_pertubations = pertubations * params['inflation_factor']
        # inflated_pertubations_vars = pertubations_vars * params['inflation_factor']

        # update the analysis vector
        analysis_vec[state_block_size:,:] = mean_params.reshape(-1,1) + inflated_pertubations
        # analysis_vec[:ndim,:] = mean_vars.reshape(-1,1) + inflated_pertubations_vars


        # check for negative thicknes and set to 1e-3 if vec_input contains h
        for i, var in enumerate(model_kwargs.get("vec_inputs",[])):
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
        # obs_index = model_kwargs.get("obs_index")
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
        parallel_write_ensemble_scattered(k+1,ens_mean, params,analysis_vec, comm_world,model_kwargs)

        # clean the memory
        del scatter_ensemble, analysis_vec; gc.collect()


# ============================ EnSRF functions ============================


# ============================ EnTKF functions ============================


# ============================ Other functions ============================

