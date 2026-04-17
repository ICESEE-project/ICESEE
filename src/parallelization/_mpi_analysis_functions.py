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

from ICESEE.src.parallelization._parallel_i_o import parallel_write_full_ensemble_from_root, \
                                                parallel_write_ensemble_scattered
from ICESEE.src.utils.tools import icesee_get_index, get_grid_dimensions

# from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager
# rank_seed, rng = ParallelManager().initialize_seed(MPI.COMM_WORLD)

# seed the random number generator
# np.random.seed(0)

# ============================ EnKF functions ============================ 
# def EnKF_X5(Cov_obs, Nens, D, HA, Eta, d): 
def EnKF_X5(km,ensemble_vec, Nens, hu_obs, model_kwargs,UtilsFunctions):
    """
    Function to compute the X5 matrix for the EnKF
        - ensemble_vec: ensemble matrix of size (ndxNens)
        - Cov_obs: observation covariance matrix
        - Nens: ensemble size
        - hu_obs: observations vector
    """
    params = model_kwargs.get("params")
    comm_world = model_kwargs.get("comm_world")
    generate_enkf_field = model_kwargs.get("generate_enkf_field", False)
    rng = model_kwargs.get("rng", np.random.default_rng())
    rank_seed = model_kwargs.get("rank_seed", 0)
    Cov_obs = None
    model_kwargs["hu_obs_loaded"] = hu_obs

    # np.random.seed(rank_seed)

    # H = UtilsFunctions(params =params, model_kwargs=model_kwargs,ensemble= ensemble_vec).JObs_fun(ensemble_vec.shape[0]) # mxNens, observation operator
    U = UtilsFunctions(params=params, model_kwargs=model_kwargs, ensemble=ensemble_vec)

    H = U.JObs_fun(ensemble_vec.shape[0])          # build once
    k_obs = model_kwargs.get("km")
    # d = U.Obs_fun(hu_obs[:, k_obs], H=H, km=km)       # reuse
    y_full = hu_obs[:, k_obs].copy()
    y_full[np.isnan(y_full)] = 0.0
    # d = U.Obs_fun(hu_obs[:, k_obs], H=H, km=k_obs)       # reuse
    d = U.Obs_fun(y_full, H=H, km=k_obs)       # reuse
    # if comm_world.Get_rank() == 0:
    #     print(f"[ICESEE] Rank: {comm_world.Get_rank()} EnKF_X5 at time step km: {km} any NaN in H: {np.isnan(H).any()} d: {np.isnan(d).any()}")
    # -- get ensemble pertubations
    use_ensemble_pertubations = model_kwargs.get("use_ensemble_pertubations", True)
    ensemble_mean = np.mean(ensemble_vec, axis=1).reshape(-1,1)
    if use_ensemble_pertubations:
        ensemble_perturbations = ensemble_vec - ensemble_mean # ensure mean is zero
        Eta = np.dot(H, ensemble_perturbations) # mxNens, ensemble pertubations
    else: #or use ensembles of perturbations
        # generate ensemble of perturbations # mxNens o---->
        if model_kwargs["joint_estimation"] or params["localization_flag"]:
            hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
        else:
            hdim = ensemble_vec.shape[0] // params["num_state_vars"]

        Lx, Ly = model_kwargs.get("Lx"), model_kwargs.get("Ly")
        len_scale =  model_kwargs.get("length_scale")  # Length scale for localization
        alpha = model_kwargs.get("alpha", 1.0)  # Mixing parameter for noise generation
        rho = model_kwargs.get("rho", 1.0)  # Correlation coefficient for noise generation
        dt = model_kwargs.get("dt", 1.0)  # Time step size
        noise = model_kwargs.get("noise", None)  # Noise vector, should be provided
        
        # ensure mean of noise is zero
        # ensure mean of noise is zero
        if model_kwargs.get("inversion_flag", False):
            # exclude friction index from params['sig_obs']
            friction_idx = int(model_kwargs.get("friction_idx", -1))  # int always
            sig_obs_filtered = [
                params["sig_obs"][i]
                for i in range(len(params["sig_obs"]))
                if i != friction_idx
            ]

        _eta = []
        for ens in range(Nens):
            noise_all = []
            for ii, sig in enumerate(sig_obs_filtered if model_kwargs.get("inversion_flag", False) else params["sig_obs"]):
                model_kwargs.update({"ii_sig": ii, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":params["total_state_param_vars"]})
                W = generate_enkf_field(**model_kwargs)
                noise_all.append(sig * W)

            noise_ = np.concatenate(noise_all, axis=0)
            _eta.append(noise_)
        _eta= np.array(_eta).T  # Convert to shape (nd, Nens)

        _eta -= np.mean(_eta, axis=1).reshape(-1, 1)  # Ensure mean is zero
        Eta = np.dot(H, _eta)  # mxNens, ensemble perturbations
        # o--->

    HAbar = np.dot(H, ensemble_mean)
    Dprime = d.reshape(-1, 1) - HAbar  # mxNens
    HAprime = copy.deepcopy(Eta)  # mxNens (requires H to be linear)

    # get the min(m,Nens)
    m = d.shape[0]
    nrmin = min(m, Nens)

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
    X1 = np.empty((nrmin, m))
    for j in range(m):
        for i in range(nrmin):
            X1[i,j] =sig[i]*U[j,i]
    
    # compute X2 = X1*Dprime # Nens x Nens
    X2 = np.dot(X1, Dprime)
    # del Cov_obs, sig, X1, Dprime; gc.collect()
    
    # print(f"[ICESEE] Rank: {rank_world} X2 shape: {X2.shape}")
    #  compute X3 = U*X2 # m_obs x Nens
    X3 = np.dot(U, X2)

    # print(f"[ICESEE] Rank: {rank_world} X3 shape: {X3.shape}")
    # compute X4 = (HAprime.T)*X3 # Nens x Nens
    X4 = np.dot(HAprime.T, X3)
    del X2, X3, U, HAprime; gc.collect()
    
    # print(f"[ICESEE] Rank: {rank_world} X4 shape: {X4.shape}")
    # compute X5 = X4 + I
    X5 = X4 + np.eye(Nens)
    # sum of each column of X5 should be 1
    if np.sum(X5, axis=0).all() != 1.0:
        print(f"[ICESEE] Sum of each X5 column is not 1.0: {np.sum(X5, axis=0)}")
    # print(f"[ICESEE] Rank: {comm_world.Get_rank()} X5 sum: {np.sum(X5, axis=0)}")
    del X4; gc.collect()

    # ===local computation
    if model_kwargs.get("local_analysis",False):
        nx, ny = model_kwargs.get("nx"), model_kwargs.get("ny")
        from scipy.spatial import distance
        # for each grid point
        h = UtilsFunctions(params =params, model_kwargs=model_kwargs,ensemble= ensemble_vec).Obs_fun 
        # d = UtilsFunctions(params =params, model_kwargs=model_kwargs,ensemble= ensemble_vec).Obs_fun(hu_obs[:,km])
        analysis_vec_ij = np.empty_like(ensemble_vec)
        dim = ensemble_vec.shape[0]//params["total_state_param_vars"]
        mx, my = get_grid_dimensions(nx, ny, dim)
        yg, xg = np.unravel_index(np.arange(dim), (my, mx))
        lscale = 20
        for ij in range(dim):
            # reference point xg[ij], yg[ij]
            # dist = np.sqrt((xg[ij] - xg)**2 + (yg[ij] - yg)**2)
            dist = distance.cdist(np.array([[xg[ij], yg[ij]]]), np.column_stack((xg, yg)))[0]
            # nearest_indices = np.argsort(dist)[:lscale]
            nearest_indices = dist < np.abs(lscale)
            # ensemble_vec_ij = ensemble_vec[nearest_indices,:]

            Eta_local = np.zeros(Nens)
            D_local   = np.zeros_like(Eta_local)
            HA_local  = np.zeros_like(D_local)
            for ens in range(Nens):
                for var in range(params["total_state_param_vars"]):
                    idx = var*dim + ij
                    # nearrest observations indices 
                    idx_obs_loc = var*dim + nearest_indices
                    print(f"[ICESEE] nearest_indices: {nearest_indices} idx_obs_loc: {idx_obs_loc}")
                    # d_loc = d[idx]
                    d_loc = d[idx_obs_loc]
                    # Cov_obs_loc = Cov_obs[idx,idx]
                    Cov_obs_loc = Cov_obs[idx_obs_loc,idx_obs_loc]
                    # mean = np.zeros(1)
                    # Eta_local[ens] = np.random.multivariate_normal(mean, cov=Cov_obs_loc)
                    Eta_local[ens] = np.random.normal(0, np.sqrt(Cov_obs_loc))
                    D_local[ens] = d_loc + Eta_local[ens]
                    # HA_local[ens] = h(ensemble_vec[idx,ens])
                    HA_local[ens] = UtilsFunctions(params, ensemble_vec[idx_obs_loc,ens]).Obs_fun(ensemble_vec[idx_obs_loc,ens])

            Dprime_local = D_local - HA_local
            HAbar_local = np.mean(HA_local)
            HAprime_local = HA_local - HAbar_local
            m_obs_local = d_loc.shape[0]
            nrmin_local = min(m_obs_local, Nens)
            HAprime_eta_local = HAprime_local + Eta_local
            U_local, sig_local, _ = np.linalg.svd(HAprime_eta_local, full_matrices=False)
            sig_local = sig_local**2
            sigsum_local = np.sum(sig_local[:nrmin_local])
            sigsum1_local = 0.0
            nrsigma_local = 0
            for i in range(nrmin_local):
                if sigsum1_local / sigsum_local < 0.999:
                    nrsigma_local += 1
                    sigsum1_local += sig_local[i]
                    sig_local[i] = 1.0 / sig_local[i]
                else:
                    sig_local[i:nrmin_local] = 0.0
                    break
            X1_local = np.empty((nrmin_local, m_obs_local))
            for j in range(m_obs_local):
                for i in range(nrmin_local):
                    X1_local[i,j] = sig_local[i]*U_local[j,i]
            X2_local = np.dot(X1_local, Dprime_local)
            X3_local = np.dot(U_local, X2_local)
            X4_local = np.dot(HAprime_local.T, X3_local)
            X5_local = X4_local + np.eye(Nens)

            # compute the diff
            X5_diff = X5_local - X5

            # compute analysis vector
            for var in range(params["total_state_param_vars"]):
                # idx = var*dim + ij
                idx = var*dim + nearest_indices
                analysis_vec_ij[ij,:] = np.dot(ensemble_vec[idx,:], X5) + np.dot(ensemble_vec[idx,:], X5_diff)
        
    else:
        analysis_vec_ij = None
        

    return X5, analysis_vec_ij

def analysis_enkf_update(k,ens_mean,ensemble_vec, shape_ens, X5,time_analysis_mean_generation,time_analysis_file_writing, analysis_vec_ij,UtilsFunctions,model_kwargs,smb_scale):
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
        time_analysis_mean_generation = BM.bcast(time_analysis_mean_generation, comm=comm_world)
        model_kwargs['X5'] = X5
        # X5_diff = BM.bcast(X5_diff, comm=comm_world)

        # initialize the an empty ensemble vector for the rest of the processors
        if rank_world != 0:
            ensemble_vec = np.empty(shape_ens, dtype=np.float64)

        # --- scatter ensemble_vec to all processors ---
        scatter_ensemble = BM.scatter(ensemble_vec, comm_world)
    
        # do the ensemble analysis update: A_j = Fj*X5 
        analysis_vec = np.dot(scatter_ensemble, X5)

        # ndim = analysis_vec.shape[0] // params["total_state_param_vars"]
        ndim = model_kwargs.get("nd", params["nd"])//len(model_kwargs.get("all_observed",[]))
        state_block_size = ndim*params["num_state_vars"]
        
        # ---> multiplicative inflation
        time_analysis_mean_generation1  = MPI.Wtime() 
        mean_params = np.mean(analysis_vec[state_block_size:,:], axis=1)
        time_analysis_mean_generation1 = MPI.Wtime() - time_analysis_mean_generation1
        time_analysis_mean_generation += time_analysis_mean_generation1

        #  compute parturbations
        pertubations = analysis_vec[state_block_size:,:] - mean_params.reshape(-1,1)
        # apply the inflation factor
        inflated_pertubations = pertubations * params['inflation_factor']

        # update the analysis vector
        analysis_vec[state_block_size:,:] = mean_params.reshape(-1,1) + inflated_pertubations
        # only inflate bed topography if it is directly observed
        observed_params = model_kwargs.get("observed_params", [])

        # gather from all processors
        # ensemble_vec = BM.allgather(analysis_vec, comm_world)
        _time_analysis_file_writing = MPI.Wtime()
        parallel_write_ensemble_scattered(k+1,ens_mean, params,analysis_vec, comm_world,model_kwargs)
        time_analysis_file_writing += MPI.Wtime() - _time_analysis_file_writing

        # clean the memory
        del scatter_ensemble, analysis_vec; gc.collect()

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

