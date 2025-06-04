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

# seed the random number generator
np.random.seed(0)

def parallel_write_ensemble_scattered(timestep, ensemble_mean, params, ensemble_chunk, comm, model_kwargs, output_file="icesee_ensemble_data.h5"):
    """
    Write ensemble data in parallel using h5py and MPI
    ensemble_chunk: local data on each rank with shape (local_nd, Nens)
    """
    # MPI setup
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Get local chunk dimensions
    local_nd = ensemble_chunk.shape[0]  # rows per rank
    Nens = ensemble_chunk.shape[1]      # ensemble members (same for all ranks)

    # Gather the number of rows from each rank
    local_nd_array = comm.gather(local_nd, root=0)
    # local_nd_array = BM.gather(local_nd, comm)
    
    if rank == 0:
        nd_total = sum(local_nd_array)
    else:
        nd_total = None
    nd_total = BM.bcast(nd_total, comm)

    # Calculate offsets for each rank
    offset = comm.scan(local_nd) - local_nd  # Exclusive scan gives starting position

    # check if _modelrun_datasets exists
    # _modelrun_datasets = f"_modelrun_datasets"
    # if rank == 0 and not os.path.exists(_modelrun_datasets):
    #     # cretate the directory
    #     os.makedirs(_modelrun_datasets, exist_ok=True)
    
    # comm.barrier() # wait for all processes to reach this point
    output_file = os.path.join(model_kwargs.get('data_path'), output_file)

    # Open file in parallel mode
    if timestep == 0.0:
        with h5py.File(output_file, 'w', driver='mpio', comm=comm) as f:
            # Create dataset with total dimensions
            dset = f.create_dataset('ensemble', (nd_total, Nens, model_kwargs.get('nt', params['nt']) +1), dtype=ensemble_chunk.dtype)
            
            # Each rank writes its chunk
            dset[offset:offset + local_nd, :,0] = ensemble_chunk

            # ens_mean 
            ens_mean = f.create_dataset('ensemble_mean', (local_nd, model_kwargs.get('nt', params['nt']) +1), dtype=ensemble_chunk.dtype)
            if rank == 0:
                ens_mean[:,0] = ensemble_mean

            DEnKF_flag = model_kwargs.get("DEnKF_flag",False)
            if DEnKF_flag:
                # dset = dset + ensemble_mean
                comm.barrier() # wait for all processes to reach this point
                if rank == 0:
                    ensemble_mean = np.mean(dset[:, :, 0], axis=1)
                    dset[:,:, 0] = dset[:, :, 0] + ensemble_mean[:,np.newaxis]
    else:
        with h5py.File(output_file, 'a', driver='mpio', comm=comm) as f:
            dset = f['ensemble']
            # dset[offset:offset + local_nd, :,timestep] = ensemble_chunk

            # ================
            if False: #TODO: test tomorrow
                # # extract bounds for the parameters
                # bounds = model_kwargs["bounds"]
                # # Function f: Linear, bijective mapping from [0,1] to [l_theta - theta^a, u_theta - theta^a]
                # def f(x, theta_a_i):
                #     # Map x (from Beta[0,1]) to the range [l_theta - theta_a_i, u_theta - theta_a_i]
                #     for i, vars in enumerate(model_kwargs["params_vec"]):
                #         param_bound = bounds[i]
                #         l_theta, u_theta = param_bound[0], param_bound[1]
                #         l_theta = np.ones((theta_a_i.shape[0],1))*l_theta
                #         u_theta = np.ones((theta_a_i.shape[0],1))*u_theta
                #         lower = l_theta - theta_a_i
                #         upper = u_theta - theta_a_i
                #         x = x*(upper - lower) + lower
                #     return x
                
                # ndim = ensemble_chunk.shape[0] // params["total_state_param_vars"]
                # state_block_size = ndim*params["num_state_vars"]
                # param_size = ensemble_chunk.shape[0] - state_block_size

                # alpha_t, beta_t = 2.0, 2.0  # Beta distribution parameters
                # X_t = beta.rvs(alpha_t, beta_t, ensemble_chunk.shape[1])
                # pertubations = np.array([f(X_t[i], ensemble_chunk[state_block_size:,i]) for i in range( ensemble_chunk.shape[1])])
                # prev_data = dset[offset:offset + local_nd, :, timestep-1]
                # ensemble_chunk[state_block_size:,:] = prev_data[state_block_size:,:] + pertubations

                # # ensure parameters stay within bounds
                # for i, vars in enumerate(model_kwargs["params_vec"]):
                #     param_bound = bounds[i]
                #     l_theta, u_theta = param_bound[0], param_bound[1]
                #     ensemble_chunk[state_block_size+i,:] = np.clip(ensemble_chunk[state_block_size+i,:], l_theta, u_theta)
                
                # dset[offset:offset + local_nd, :,timestep] = ensemble_chunk

                # ----------
                prev_data = dset[offset:offset + local_nd, :, timestep-1]
                n_params = len(model_kwargs["params_vec"])
                # Extract bounds
                bounds = model_kwargs["bounds"]

                # Fixed function f
                def func(x, theta_a_i, l_theta, u_theta):
                    # lower = l_theta - theta_a_i
                    # upper = u_theta - theta_a_i
                    # return lower + x * (upper - lower)
                    # scale = u_theta - l_theta
                    # return (x-0.5)*scale 
                    scale = u_theta - l_theta
                    current_spread = np.std(theta_a_i, axis=1, keepdims=True)
                    adaptive_scale = np.maximum(scale, current_spread * 2.0)  # Boost spread
                    return (x - 0.5) * adaptive_scale

                # Dimensions
                ndim = ensemble_chunk.shape[0] // params["total_state_param_vars"]  # e.g., 50 / 4 = 12
                state_block_size = ndim * params["num_state_vars"]  # e.g., 12 * 1 = 12
                param_size = ensemble_chunk.shape[0] - state_block_size  # e.g., 50 - 12 = 38

                # Perturbations
                alpha_t, beta_t = 2.0, 2.0
                X_t = beta.rvs(alpha_t, beta_t, size=(param_size, ensemble_chunk.shape[1]))
                perturbations = np.zeros((param_size, ensemble_chunk.shape[1]))
                # print(bounds)
                # print(bounds[0][0], bounds[0][1])
                # bounds = np.array([0.5, 1.9])
                for i in range(n_params):
                    l_theta, u_theta = bounds[i]
                    # l_theta, u_theta = bounds
                    idx_start = i * ndim
                    idx_end = (i + 1) * ndim
                    # param_block = ensemble_chunk[state_block_size + idx_start:state_block_size + idx_end, :]
                    param_block = prev_data[state_block_size + idx_start:state_block_size + idx_end, :]
                    perturbations[idx_start:idx_end, :] = func(X_t[idx_start:idx_end, :], param_block, l_theta, u_theta)

                # Update ensemble
                # prev_data = dset[offset:offset + local_nd, :, timestep-1]
                ensemble_chunk[state_block_size:, :] = prev_data[state_block_size:, :] + perturbations

                # Enforce bounds
                for i in range(n_params):
                    # l_theta, u_theta = bounds[i]
                    # l_theta, u_theta = bounds
                    idx_start = state_block_size + i * ndim
                    idx_end = idx_start + ndim
                    # ensemble_chunk[idx_start:idx_end, :] = np.clip(ensemble_chunk[idx_start:idx_end, :], l_theta, u_theta)

                # Write to dataset
                dset[offset:offset + local_nd, :, timestep] = ensemble_chunk


            # =================
            if False:
                ndim = ensemble_chunk.shape[0] // params["total_state_param_vars"]
                state_block_size = ndim*params["num_state_vars"]
                param_size = ensemble_chunk.shape[0] - state_block_size
                alpha = np.ones(param_size)*2.0
                beta_param = alpha
                def compute_f_params(alpha, beta_param):
                    mean_x = alpha/(alpha+beta_param)
                    a = 1.0
                    b = -a*mean_x
                    return a,b

                def update_theta(alpha, beta_param):
                    # theta_f_t = np.zeros_like(theta_prev)
                    f_x_ti = np.zeros((param_size,ensemble_chunk.shape[1]))
                    for i in range(ensemble_chunk.shape[1]):
                        a,b = compute_f_params(alpha[i], beta_param[i])
                        x_ti = beta.rvs(alpha[i], beta_param[i])
                        
                        f_x_ti[:,i] = a*x_ti + b

                        # theta_f_t[:,i] = theta_prev[:,i] + f_x_ti
                    # return theta_f_t
                    return f_x_ti
                
                # Update ensemble_chunk before writing
                if state_block_size < ensemble_chunk.shape[0]:
                    prev_data = dset[offset:offset + local_nd, :, timestep-1]
                    ensemble_chunk[state_block_size:,:] = prev_data[state_block_size:,:] + update_theta(alpha, beta_param)

                if False:
                     # ----------
                    n_params = len(model_kwargs["params_vec"])
                    # Extract bounds
                    # bounds = model_kwargs["bounds"]
                    bounds = np.array([0.2, 1.3])
                    # Enforce bounds
                    for i in range(n_params):
                        # l_theta, u_theta = bounds[i]
                        l_theta, u_theta = bounds
                        idx_start = state_block_size + i * ndim
                        idx_end = idx_start + ndim
                        ensemble_chunk[idx_start:idx_end, :] = np.clip(ensemble_chunk[idx_start:idx_end, :], l_theta, u_theta)

            # ensemble_chunk[state_block_size:,:] =  dset[offset:offset + local_nd, :,timestep-1] + update_theta(alpha, beta_param)
            dset[offset:offset + local_nd, :,timestep] = ensemble_chunk


            # ================

            if rank == 0:
                ens_mean = f['ensemble_mean']
                ens_mean[:,timestep] = ensemble_mean

            DEnKF_flag = model_kwargs.get("DEnKF_flag",False)
            if DEnKF_flag:
                comm.barrier() # wait for all processes to reach this point
                if rank == 0:
                    # dset = dset + ensemble_mean
                    ensemble_mean = np.mean(dset[:, :, timestep], axis=1)
                    dset[:,:, timestep] = dset[:, :, timestep] #- ensemble_mean[:,np.newaxis]

    comm.Barrier()

def parallel_write_data_from_root_2D(full_ensemble=None, comm=None, data_name=None, output_file="preliminary_data.h5"):
    """
    Write ensemble data in parallel where full matrix exists on rank 0
    full_ensemble: complete matrix on rank 0 with shape (nd, Nens)
    """
    # MPI setup
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Get dimensions on root and broadcast
    if rank == 0:
        nd = full_ensemble.shape[0]
        Nens = full_ensemble.shape[1]
        dtype = full_ensemble.dtype
    else:
        nd = None
        Nens = None
        dtype = None
    
    nd = comm.bcast(nd, root=0)
    Nens = comm.bcast(Nens, root=0)
    dtype = comm.bcast(dtype, root=0)

    # Calculate local chunk sizes
    local_nd = nd // size  # Base size per rank
    remainder = nd % size  # Extra rows to distribute
    
    # Determine local size and offset for each rank
    if rank < remainder:
        local_nd += 1  # Distribute remainder to first few ranks
    offset = rank * (nd // size) + min(rank, remainder)

    # Scatter the data (only if rank 0 has it)
    if rank == 0:
        chunks = np.array_split(full_ensemble, size, axis=0)
    else:
        chunks = None
    
    local_chunk = BM.scatter(chunks, comm)
    
    # comm.barrier() # wait for all processes to reach this point
    output_file = os.path.join("_modelrun_datasets", output_file)

    # Open file in parallel mode
    with h5py.File(output_file, 'w', driver='mpio', comm=comm) as f:
        # Create dataset with total dimensions
        dset = f.create_dataset(data_name, (nd, Nens), dtype=dtype)
        
        # Each rank writes its chunk
        dset[offset:offset + local_nd, :] = local_chunk

def parallel_write_vector_from_root(full_ensemble=None, comm=None, data_shape=None, data_name=None, output_file="icesee_ensemble_data.h5"):
    """
    Append ensemble data in parallel where the full matrix exists on rank 0.
    Each call appends a new time step, resulting in a dataset of shape (nd, Nens, nt).

    full_ensemble: complete matrix on rank 0 with shape (nd, Nens)
    comm: MPI communicator
    output_file: Name of the output HDF5 file
    """
    # MPI setup
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Get dimensions on root and broadcast
    if rank == 0:
        nd, Nens = full_ensemble.shape
        dtype = full_ensemble.dtype
    else:
        nd, Nens, dtype = None, None, None
    
    nd = comm.bcast(nd, root=0)
    Nens = comm.bcast(Nens, root=0)
    dtype = comm.bcast(dtype, root=0)

    # Calculate local chunk sizes
    local_nd = nd // size
    remainder = nd % size

    if rank < remainder:
        local_nd += 1
    offset = rank * (nd // size) + min(rank, remainder)

    # Scatter the data from rank 0
    if rank == 0:
        chunks = np.array_split(full_ensemble, size, axis=0)
    else:
        chunks = None
    
    local_chunk = BM.scatter(chunks, comm)

    # Define output file path
    output_file = os.path.join("_modelrun_datasets", output_file)

    # Open file in parallel mode
    with h5py.File(output_file, 'w', driver='mpio', comm=comm) as f:
        # Create dataset with total dimensions
        #  data_shape should be a tuple
        dset = f.create_dataset(data_name, data_shape, dtype=dtype)
        
        # Each rank writes its chunk
        dset[offset:offset + local_nd, :,0] = local_chunk

    comm.Barrier()


    
def parallel_write_full_ensemble_from_root(timestep, ensemble_mean, model_kwargs,full_ensemble=None, comm=None, output_file="icesee_ensemble_data.h5"):
    """
    Append ensemble data in parallel where the full matrix exists on rank 0.
    Each call appends a new time step, resulting in a dataset of shape (nd, Nens, nt).

    full_ensemble: complete matrix on rank 0 with shape (nd, Nens)
    comm: MPI communicator
    output_file: Name of the output HDF5 file
    """
    params = model_kwargs.get("params")

    # MPI setup
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Get dimensions on root and broadcast
    if rank == 0:
        nd, Nens = full_ensemble.shape
        dtype = full_ensemble.dtype
    else:
        nd, Nens, dtype = None, None, None
    
    nd = comm.bcast(nd, root=0)
    Nens = comm.bcast(Nens, root=0)
    dtype = comm.bcast(dtype, root=0)

    # Calculate local chunk sizes
    local_nd = nd // size
    remainder = nd % size

    if rank < remainder:
        local_nd += 1
    offset = rank * (nd // size) + min(rank, remainder)

    # Scatter the data from rank 0
    if rank == 0:
        chunks = np.array_split(full_ensemble, size, axis=0)
    else:
        chunks = None
    
    local_chunk = BM.scatter(chunks, comm)

    # Define output file path
    output_file = os.path.join(params.get('data_path'), output_file)

    # Open file in parallel mode
    if timestep == 0:
        with h5py.File(output_file, 'w', driver='mpio', comm=comm) as f:
            # Create dataset with total dimensions
            dset = f.create_dataset('ensemble', (nd, Nens, model_kwargs.get('nt', params['nt'])+1), dtype=dtype)
            
            # Each rank writes its chunk
            dset[offset:offset + local_nd, :,0] = local_chunk

            # ens_mean 
            ens_mean = f.create_dataset('ensemble_mean', (nd, model_kwargs.get('nt', params['nt'])+1), dtype=dtype)
            if rank == 0:
                ens_mean[:,0] = ensemble_mean
    else:
        with h5py.File(output_file, 'a', driver='mpio', comm=comm) as f:
            dset = f['ensemble']
            dset[offset:offset + local_nd, :,timestep] = local_chunk

            if rank == 0:
                ens_mean = f['ensemble_mean']
                ens_mean[:,timestep] = ensemble_mean
    comm.Barrier()

def get_grid_dimensions(nx, ny, ndim):
    """
    Calculate grid dimensions mx and my based on physical dimensions and total points.
    
    Parameters:
    nx (int): Number of elements in x-direction
    ny (int): Number of elements in y-direction
    ndim (int): Total number of grid points (mx * my)
    
    Returns:
    tuple: (mx, my) - number of grid points in x and y directions
    """
    # Calculate aspect ratio from physical dimensions
    alpha = nx / ny
    
    # Initial estimate based on aspect ratio and ndim
    # mx/my = alpha and mx*my = ndim
    # mx = sqrt(ndim * alpha), my = sqrt(ndim / alpha)
    mx = np.sqrt(ndim * alpha)
    my = np.sqrt(ndim / alpha)
    
    # Initial rounding
    if mx - int(mx) > 0.5:
        mx = int(np.ceil(mx))
        my = int(np.floor(my))
    elif my - int(my) > 0.5:
        my = int(np.ceil(my))
        mx = int(np.floor(mx))
    else:
        mx, my = int(mx), int(my)
    
    # Quick adjustment to reach ndim
    current_product = mx * my
    if current_product != ndim:
        # Calculate scale factor
        scale = np.sqrt(ndim / current_product)
        mx = int(round(mx * scale))
        my = int(round(my * scale))
        
        # Fast fine-tuning with minimal iterations
        product = mx * my
        if product < ndim:
            while product < ndim:
                if mx/my < alpha:
                    mx += 1
                else:
                    my += 1
                product = mx * my
        elif product > ndim:
            while product > ndim:
                if mx/my > alpha:
                    mx -= 1
                else:
                    my -= 1
                product = mx * my
    
    return mx, my
# ============================ EnKF functions ============================ 
# def EnKF_X5(Cov_obs, Nens, D, HA, Eta, d): 
def EnKF_X5(k,ensemble_vec, Cov_obs, Nens, d, model_kwargs,UtilsFunctions):
    """
    Function to compute the X5 matrix for the EnKF
        - ensemble_vec: ensemble matrix of size (ndxNens)
        - Cov_obs: observation covariance matrix
        - Nens: ensemble size
        - d: observation vector
    """
    params = model_kwargs.get("params")
    comm_world = model_kwargs.get("comm_world")
    H = UtilsFunctions(params, ensemble_vec).JObs_fun(ensemble_vec.shape[0]) # mxNens, observation operator

    # -- get ensemble pertubations
    ensemble_perturbations = ensemble_vec - np.mean(ensemble_vec, axis=1).reshape(-1,1)
    
    # ----parallelize this step
    Eta = np.zeros((d.shape[0], Nens)) # mxNens, ensemble pertubations
    Eta = np.dot(H, ensemble_perturbations) # mxNens, ensemble pertubations
    D   = np.zeros_like(Eta) # mxNens #virtual observations
    HA  = np.zeros_like(D)
    for ens in range(Nens):
        # Eta[:,ens] = np.random.multivariate_normal(mean=np.zeros(d.shape[0]), cov=Cov_obs) 
        D[:,ens] = d + Eta[:,ens]
        HA[:,ens] = np.dot(H, ensemble_vec[:,ens])
    # ---------------------------------------

    # --- compute the innovations D` = D-HA
    Dprime = D - HA # mxNens

    # --- compute HAbar
    HAbar = np.mean(HA, axis=1) # mx1
    # --- compute HAprime
    # HAprime = HA - HAbar.reshape(-1,1) # mxNens (requires H to be linear)
    
    # Aprime = ensemble_vec@(np.eye(Nens) - one_N) # mxNens
    one_N = np.ones((Nens,Nens))/Nens
    HAprime= HA@(np.eye(Nens) - one_N) # mxNens

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
    X2 = np.dot(X1, Dprime)
    # del Cov_obs, sig, X1, Dprime; gc.collect()
    
    # print(f"Rank: {rank_world} X2 shape: {X2.shape}")
    #  compute X3 = U*X2 # m_obs x Nens
    X3 = np.dot(U, X2)

    # print(f"Rank: {rank_world} X3 shape: {X3.shape}")
    # compute X4 = (HAprime.T)*X3 # Nens x Nens
    X4 = np.dot(HAprime.T, X3)
    del X2, X3, U, HAprime; gc.collect()
    
    # print(f"Rank: {rank_world} X4 shape: {X4.shape}")
    # compute X5 = X4 + I
    X5 = X4 + np.eye(Nens)
    # sum of each column of X5 should be 1
    if np.sum(X5, axis=0).all() != 1.0:
        print(f"Sum of each X5 column is not 1.0: {np.sum(X5, axis=0)}")
    # print(f"Rank: {comm_world.Get_rank()} X5 sum: {np.sum(X5, axis=0)}")
    del X4; gc.collect()

    # ===local computation
    if model_kwargs.get("local_analysis",False):
        nx, ny = model_kwargs.get("nx"), model_kwargs.get("ny")
        from scipy.spatial import distance
        # for each grid point
        h = UtilsFunctions(params, ensemble_vec).Obs_fun 
        # d = UtilsFunctions(params, ensemble_vec).Obs_fun(hu_obs[:,km])
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
                    print(f"nearest_indices: {nearest_indices} idx_obs_loc: {idx_obs_loc}")
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

def analysis_enkf_update(k,ens_mean,ensemble_vec, shape_ens, X5, analysis_vec_ij,UtilsFunctions,model_kwargs,smb_scale):
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

        ndim = analysis_vec.shape[0] // params["total_state_param_vars"]
        state_block_size = ndim*params["num_state_vars"]
        # analysis_vec[state_block_size:,:] /= 10
        # analysis_vec[state_block_size:,:] *= (smb_scale)  # Scale SMB after analysis
        # params['inflation_factor'] = 1.1
        # analysis_vec = UtilsFunctions(params,  analysis_vec).inflate_ensemble(in_place=True)
        # ---> multiplicative inflation
        mean_params = np.mean(analysis_vec[state_block_size:,:], axis=1)
        #  compute parturbations
        pertubations = analysis_vec[state_block_size:,:] - mean_params.reshape(-1,1)
        # apply the inflation factor
        inflated_pertubations = pertubations * params['inflation_factor']

        # update the analysis vector
        analysis_vec[state_block_size:,:] = mean_params.reshape(-1,1) + inflated_pertubations


        # check for negative thicknes and set to 1e-3 if vec_input contains h
        for i, var in enumerate(model_kwargs.get("vec_inputs",[])):
            if var == "h":
                start = i * ndim
                end = start + ndim
                analysis_vec[start:end, :] = np.maximum(analysis_vec[start:end, :], 1e-2)

        # dynamical model for parameters: from https://doi.org/10.1002/qj.3257
        # obs_index = model_kwargs.get("obs_index")
        # # #  check if k equals to the first observation index
        # # print(f"Rank: {rank_world} km: {km} obs_index: {obs_index}")
        # if  (k+1 == obs_index[0]):
        # #     print(f"[Debug] Rank: {rank_world} k: {km} obs_index: {obs_index}")
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
    H = UtilsFunctions(params, ensemble_vec).JObs_fun(ensemble_vec.shape[0]) # mxNens, observation operator

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
    
    # print(f"Rank: {rank_world} X2 shape: {X2.shape}")
    #  compute X3 = U*X2 # m_obs x Nens
    X3 = np.dot(U, X2)
    X3prime = np.dot(U, X2prime) # m_obs x Nens

    # print(f"Rank: {rank_world} X3 shape: {X3.shape}")
    # compute X4 = (HAprime.T)*X3 # Nens x Nens
    X4 = np.dot(HAprime.T, X3)
    X4prime = np.dot(HAprime.T, X3prime) # Nens x Nens
    del X2, X3, U, HAprime; gc.collect()
    
    # print(f"Rank: {rank_world} X4 shape: {X4.shape}")
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
        print(f"Sum of each X5 column is not 1.0: {np.sum(X5, axis=0)}")
    # print(f"Rank: {comm_world.Get_rank()} X5 sum: {np.sum(X5, axis=0)}")
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

        # print(f"Rank: {rank_world} analysis_vec shape: {analysis_vec.shape}, ens_mean shape: {ens_mean.shape}")

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
        #  compute parturbations
        pertubations = analysis_vec[state_block_size:,:] - mean_params.reshape(-1,1)
        # apply the inflation factor
        inflated_pertubations = pertubations * params['inflation_factor']

        # update the analysis vector
        analysis_vec[state_block_size:,:] = mean_params.reshape(-1,1) + inflated_pertubations


        # check for negative thicknes and set to 1e-3 if vec_input contains h
        for i, var in enumerate(model_kwargs.get("vec_inputs",[])):
            if var == "h":
                start = i * ndim
                end = start + ndim
                analysis_vec[start:end, :] = np.maximum(analysis_vec[start:end, :], 1e-2)

        # dynamical model for parameters: from https://doi.org/10.1002/qj.3257
        # obs_index = model_kwargs.get("obs_index")
        # # #  check if k equals to the first observation index
        # # print(f"Rank: {rank_world} km: {km} obs_index: {obs_index}")
        # if  (k+1 == obs_index[0]):
        # #     print(f"[Debug] Rank: {rank_world} k: {km} obs_index: {obs_index}")
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

def gather_and_broadcast_data_default_run(updated_state, subcomm, sub_rank, comm_world, rank_world, params):
    """
    Gathers, processes, and broadcasts ensemble data across MPI processes.

    Parameters:
    - updated_state: dict, contains state variables to be gathered
    - subcomm: MPI communicator for subgroups
    - sub_rank: int, rank within the subcommunicator
    - comm_world: MPI communicator for all processes
    - rank_world: int, rank within the world communicator
    - params: dict, contains necessary parameters like "total_state_param_vars"
    - BM: object with a `bcast` method for broadcasting data

    Returns:
    - ensemble_vec: The processed and broadcasted ensemble data
    """

    # Step 1: Gather data from all sub-ranks
    global_data = {key: subcomm.gather(data, root=0) for key, data in updated_state.items()}

    # Step 2: Process on sub_rank 0
    if sub_rank == 0:
        for key in global_data:
            global_data[key] = np.hstack(global_data[key])

        # Stack all variables into a single array
        stacked = np.hstack([global_data[key] for key in updated_state.keys()])
        shape_ = np.array(stacked.shape, dtype=np.int32)
    else:
        shape_ = np.empty(2, dtype=np.int32)

    # Step 3: Broadcast the shape to all processors
    shape_ = comm_world.bcast(shape_, root=0)

    # Step 4: Prepare the stacked array for non-root sub-ranks
    if sub_rank != 0:
        stacked = np.empty(shape_, dtype=np.float64)

    # Step 5: Gather the stacked arrays from all sub-ranks
    all_ens = comm_world.gather(stacked if sub_rank == 0 else None, root=0)

    # Step 6: Final processing on world rank 0
    if rank_world == 0:
        all_ens = [arr for arr in all_ens if isinstance(arr, np.ndarray)]
        ensemble_vec = np.column_stack(all_ens)
        hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
    else:
        ensemble_vec = np.empty((shape_[0], params["Nens"]), dtype=np.float64)

    # Step 7: Broadcast the final ensemble vector
    # ensemble_vec = BM.bcast(ensemble_vec, comm_world)

    return ensemble_vec, shape_

