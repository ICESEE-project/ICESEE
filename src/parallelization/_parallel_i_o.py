# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-03-06
# @description: This module contains functions for parallel I/O operations in the ICESEE model.
# =============================================================================

import gc
import os
import copy
import h5py
import numpy as np
import bigmpi4py as BM
from scipy.stats import multivariate_normal, beta
from mpi4py import MPI
from ICESEE.src.utils.tools import icesee_get_index

from ICESEE.src.utils.inference_plugin import (
    apply_bed_observation_anchor_global,
    apply_bed_domain_gate_global,
    apply_global_inference_hook,
)


# def parallel_write_ensemble_scattered(timestep, ensemble_mean, params, ensemble_chunk, comm, model_kwargs, output_file="icesee_ensemble_data.h5"):
#     """
#     Write ensemble data in parallel using h5py and MPI
#     ensemble_chunk: local data on each rank with shape (local_nd, Nens)
#     """
#     # MPI setup
#     rank = comm.Get_rank()
#     size = comm.Get_size()

#     # Get local chunk dimensions
#     local_nd = ensemble_chunk.shape[0]  # rows per rank
#     Nens = ensemble_chunk.shape[1]      # ensemble members (same for all ranks)

#     # Gather the number of rows from each rank
#     local_nd_array = comm.gather(local_nd, root=0)
#     # local_nd_array = BM.gather(local_nd, comm)
    
#     if rank == 0:
#         nd_total = sum(local_nd_array)
#     else:
#         nd_total = None
#     nd_total = BM.bcast(nd_total, comm)

#     # Calculate offsets for each rank
#     offset = comm.scan(local_nd) - local_nd  # Exclusive scan gives starting position

#     # check if _modelrun_datasets exists
#     # _modelrun_datasets = f"_modelrun_datasets"
#     # if rank == 0 and not os.path.exists(_modelrun_datasets):
#     #     # cretate the directory
#     #     os.makedirs(_modelrun_datasets, exist_ok=True)
    
#     # comm.barrier() # wait for all processes to reach this point
#     output_file = os.path.join(model_kwargs.get('data_path'), output_file)

#     # Open file in parallel mode
#     if timestep == 0.0:
#         with h5py.File(output_file, 'w', driver='mpio', comm=comm) as f:
#             # Create dataset with total dimensions
#             dset = f.create_dataset('ensemble', (nd_total, Nens, model_kwargs.get('nt', params['nt']) +1), dtype=ensemble_chunk.dtype)
            
#             # Each rank writes its chunk
#             dset[offset:offset + local_nd, :,0] = ensemble_chunk

#             # ens_mean 
#             ens_mean = f.create_dataset('ensemble_mean', (local_nd, model_kwargs.get('nt', params['nt']) +1), dtype=ensemble_chunk.dtype)
#             if rank == 0:
#                 ens_mean[:,0] = ensemble_mean

#             DEnKF_flag = model_kwargs.get("DEnKF_flag",False)
#             if DEnKF_flag:
#                 # dset = dset + ensemble_mean
#                 comm.barrier() # wait for all processes to reach this point
#                 if rank == 0:
#                     ensemble_mean = np.mean(dset[:, :, 0], axis=1)
#                     dset[:,:, 0] = dset[:, :, 0] + ensemble_mean[:,np.newaxis]
#     else:
#         with h5py.File(output_file, 'a', driver='mpio', comm=comm) as f:
#             dset = f['ensemble']
#             # dset[offset:offset + local_nd, :,timestep] = ensemble_chunk

#             # ================
#             if False: #TODO: test tomorrow
#                 # # extract bounds for the parameters
#                 # bounds = model_kwargs["bounds"]
#                 # # Function f: Linear, bijective mapping from [0,1] to [l_theta - theta^a, u_theta - theta^a]
#                 # def f(x, theta_a_i):
#                 #     # Map x (from Beta[0,1]) to the range [l_theta - theta_a_i, u_theta - theta_a_i]
#                 #     for i, vars in enumerate(model_kwargs["params_vec"]):
#                 #         param_bound = bounds[i]
#                 #         l_theta, u_theta = param_bound[0], param_bound[1]
#                 #         l_theta = np.ones((theta_a_i.shape[0],1))*l_theta
#                 #         u_theta = np.ones((theta_a_i.shape[0],1))*u_theta
#                 #         lower = l_theta - theta_a_i
#                 #         upper = u_theta - theta_a_i
#                 #         x = x*(upper - lower) + lower
#                 #     return x
                
#                 # ndim = ensemble_chunk.shape[0] // params["total_state_param_vars"]
#                 # state_block_size = ndim*params["num_state_vars"]
#                 # param_size = ensemble_chunk.shape[0] - state_block_size

#                 # alpha_t, beta_t = 2.0, 2.0  # Beta distribution parameters
#                 # X_t = beta.rvs(alpha_t, beta_t, ensemble_chunk.shape[1])
#                 # pertubations = np.array([f(X_t[i], ensemble_chunk[state_block_size:,i]) for i in range( ensemble_chunk.shape[1])])
#                 # prev_data = dset[offset:offset + local_nd, :, timestep-1]
#                 # ensemble_chunk[state_block_size:,:] = prev_data[state_block_size:,:] + pertubations

#                 # # ensure parameters stay within bounds
#                 # for i, vars in enumerate(model_kwargs["params_vec"]):
#                 #     param_bound = bounds[i]
#                 #     l_theta, u_theta = param_bound[0], param_bound[1]
#                 #     ensemble_chunk[state_block_size+i,:] = np.clip(ensemble_chunk[state_block_size+i,:], l_theta, u_theta)
                
#                 # dset[offset:offset + local_nd, :,timestep] = ensemble_chunk

#                 # ----------
#                 prev_data = dset[offset:offset + local_nd, :, timestep-1]
#                 n_params = len(model_kwargs["params_vec"])
#                 # Extract bounds
#                 bounds = model_kwargs["bounds"]

#                 # Fixed function f
#                 def func(x, theta_a_i, l_theta, u_theta):
#                     # lower = l_theta - theta_a_i
#                     # upper = u_theta - theta_a_i
#                     # return lower + x * (upper - lower)
#                     # scale = u_theta - l_theta
#                     # return (x-0.5)*scale 
#                     scale = u_theta - l_theta
#                     current_spread = np.std(theta_a_i, axis=1, keepdims=True)
#                     adaptive_scale = np.maximum(scale, current_spread * 2.0)  # Boost spread
#                     return (x - 0.5) * adaptive_scale

#                 # Dimensions
#                 ndim = ensemble_chunk.shape[0] // params["total_state_param_vars"]  # e.g., 50 / 4 = 12
#                 state_block_size = ndim * params["num_state_vars"]  # e.g., 12 * 1 = 12
#                 param_size = ensemble_chunk.shape[0] - state_block_size  # e.g., 50 - 12 = 38

#                 # Perturbations
#                 alpha_t, beta_t = 2.0, 2.0
#                 X_t = beta.rvs(alpha_t, beta_t, size=(param_size, ensemble_chunk.shape[1]))
#                 perturbations = np.zeros((param_size, ensemble_chunk.shape[1]))
#                 # print(bounds)
#                 # print(bounds[0][0], bounds[0][1])
#                 # bounds = np.array([0.5, 1.9])
#                 for i in range(n_params):
#                     l_theta, u_theta = bounds[i]
#                     # l_theta, u_theta = bounds
#                     idx_start = i * ndim
#                     idx_end = (i + 1) * ndim
#                     # param_block = ensemble_chunk[state_block_size + idx_start:state_block_size + idx_end, :]
#                     param_block = prev_data[state_block_size + idx_start:state_block_size + idx_end, :]
#                     perturbations[idx_start:idx_end, :] = func(X_t[idx_start:idx_end, :], param_block, l_theta, u_theta)

#                 # Update ensemble
#                 # prev_data = dset[offset:offset + local_nd, :, timestep-1]
#                 ensemble_chunk[state_block_size:, :] = prev_data[state_block_size:, :] + perturbations

#                 # Enforce bounds
#                 for i in range(n_params):
#                     # l_theta, u_theta = bounds[i]
#                     # l_theta, u_theta = bounds
#                     idx_start = state_block_size + i * ndim
#                     idx_end = idx_start + ndim
#                     # ensemble_chunk[idx_start:idx_end, :] = np.clip(ensemble_chunk[idx_start:idx_end, :], l_theta, u_theta)

#                 # Write to dataset
#                 dset[offset:offset + local_nd, :, timestep] = ensemble_chunk


#             # =================
#             if False:
#                 ndim = ensemble_chunk.shape[0] // params["total_state_param_vars"]
#                 state_block_size = ndim*params["num_state_vars"]
#                 param_size = ensemble_chunk.shape[0] - state_block_size
#                 alpha = np.ones(param_size)*2.0
#                 beta_param = alpha
#                 def compute_f_params(alpha, beta_param):
#                     mean_x = alpha/(alpha+beta_param)
#                     a = 1.0
#                     b = -a*mean_x
#                     return a,b

#                 def update_theta(alpha, beta_param):
#                     # theta_f_t = np.zeros_like(theta_prev)
#                     f_x_ti = np.zeros((param_size,ensemble_chunk.shape[1]))
#                     for i in range(ensemble_chunk.shape[1]):
#                         a,b = compute_f_params(alpha[i], beta_param[i])
#                         x_ti = beta.rvs(alpha[i], beta_param[i])
                        
#                         f_x_ti[:,i] = a*x_ti + b

#                         # theta_f_t[:,i] = theta_prev[:,i] + f_x_ti
#                     # return theta_f_t
#                     return f_x_ti
                
#                 # Update ensemble_chunk before writing
#                 if state_block_size < ensemble_chunk.shape[0]:
#                     prev_data = dset[offset:offset + local_nd, :, timestep-1]
#                     ensemble_chunk[state_block_size:,:] = prev_data[state_block_size:,:] + update_theta(alpha, beta_param)

#                 if False:
#                      # ----------
#                     n_params = len(model_kwargs["params_vec"])
#                     # Extract bounds
#                     # bounds = model_kwargs["bounds"]
#                     bounds = np.array([0.2, 1.3])
#                     # Enforce bounds
#                     for i in range(n_params):
#                         # l_theta, u_theta = bounds[i]
#                         l_theta, u_theta = bounds
#                         idx_start = state_block_size + i * ndim
#                         idx_end = idx_start + ndim
#                         ensemble_chunk[idx_start:idx_end, :] = np.clip(ensemble_chunk[idx_start:idx_end, :], l_theta, u_theta)

#             # ensemble_chunk[state_block_size:,:] =  dset[offset:offset + local_nd, :,timestep-1] + update_theta(alpha, beta_param)
#             dset[offset:offset + local_nd, :,timestep] = ensemble_chunk


#             # ================

#             if rank == 0:
#                 ens_mean = f['ensemble_mean']
#                 ens_mean[:,timestep] = ensemble_mean

#             DEnKF_flag = model_kwargs.get("DEnKF_flag",False)
#             if DEnKF_flag:
#                 comm.barrier() # wait for all processes to reach this point
#                 if rank == 0:
#                     # dset = dset + ensemble_mean
#                     ensemble_mean = np.mean(dset[:, :, timestep], axis=1)
#                     dset[:,:, timestep] = dset[:, :, timestep] #- ensemble_mean[:,np.newaxis]

#     comm.Barrier()

def parallel_write_ensemble_scattered(
    timestep,
    ensemble_mean,
    params,
    ensemble_chunk,
    comm,
    model_kwargs,
    output_file="icesee_ensemble_data.h5",
    forecast_chunk=None,
):
    """
    Write ensemble data using h5py and MPI, with only rank 0 writing to the dataset.
    Optimized for large datasets and many processes using MPI_Gatherv, without parallel I/O.

    This version replaces the single-mean inversion with a member-wise parallel inversion:
      - rank 0 gathers the full ensemble
      - all ranks participate in inversion
      - each rank processes a subset of ensemble members
      - updated members are gathered back to rank 0
      - rank 0 writes the final result

    Parameters
    ----------
    timestep : int
        Current timestep index.
    ensemble_mean : ndarray
        Mean ensemble at current timestep.
    params : dict
        Global parameters dictionary.
    ensemble_chunk : ndarray
        Local data on each rank with shape (local_nd, Nens).
    comm : MPI.Comm
        MPI communicator.
    model_kwargs : dict
        Model-specific arguments.
    output_file : str
        Output HDF5 filename.
    """
    import os
    import gc
    import copy
    import h5py
    import numpy as np
    from mpi4py import MPI

    # ---------------- MPI setup ----------------
    rank = comm.Get_rank()
    size = comm.Get_size()

    # ---------------- Local dimensions ----------------
    local_nd = ensemble_chunk.shape[0]
    Nens = ensemble_chunk.shape[1]

    # ---------------- Gather local row counts ----------------
    local_nd_array = comm.gather(local_nd, root=0)

    if rank == 0:
        nd_total = sum(local_nd_array)
        counts = [n * Nens for n in local_nd_array]
        displs = [sum(counts[:i]) for i in range(size)]
        recvbuf = np.empty((nd_total, Nens), dtype=ensemble_chunk.dtype)
    else:
        nd_total = None
        counts = None
        displs = None
        recvbuf = None

    nd_total = comm.bcast(nd_total, root=0)

    # ---------------- Gather scattered state to rank 0 ----------------
    mpi_dtype = MPI.DOUBLE
    if ensemble_chunk.dtype == np.float32:
        mpi_dtype = MPI.FLOAT
    elif ensemble_chunk.dtype == np.int32:
        mpi_dtype = MPI.INT
    elif ensemble_chunk.dtype == np.int64:
        mpi_dtype = MPI.LONG

    comm.Gatherv(ensemble_chunk, [recvbuf, counts, displs, mpi_dtype], root=0)

    forecast_recvbuf = None
    if (
        timestep not in (0, 0.0)
        and str(model_kwargs.get("bed_update_domain", "all")).lower()
        != "all"
    ):
        if forecast_chunk is None:
            raise ValueError(
                "grounded-only bed updates require the forecast ensemble chunk"
            )
        if rank == 0:
            forecast_recvbuf = np.empty(
                (nd_total, Nens), dtype=forecast_chunk.dtype
            )
        comm.Gatherv(
            forecast_chunk,
            [forecast_recvbuf, counts, displs, mpi_dtype],
            root=0,
        )

    output_file = os.path.join(model_kwargs.get("data_path"), output_file)

    # ============================================================
    # timestep = 0: only initialize file
    # ============================================================
    if timestep == 0 or timestep == 0.0:
        if rank == 0:
            with h5py.File(output_file, "w") as f:
                dset = f.create_dataset(
                    "ensemble",
                    (nd_total, Nens, model_kwargs.get("nt", params["nt"]) + 1),
                    dtype=ensemble_chunk.dtype
                )
                ens_mean = f.create_dataset(
                    "ensemble_mean",
                    (nd_total, model_kwargs.get("nt", params["nt"]) + 1),
                    dtype=ensemble_chunk.dtype
                )

                dset[:, :, 0] = recvbuf
                ens_mean[:, 0] = ensemble_mean

                if model_kwargs.get("DEnKF_flag", False):
                    mean0 = np.mean(dset[:, :, 0], axis=1)
                    dset[:, :, 0] += mean0[:, np.newaxis]

        comm.Barrier()
        return

    # ============================================================
    # timestep > 0
    # ============================================================

    # ------------------------------------------------------------
    # Step 1: rank 0 prepares recvbuf (bed relax + ISSM fixes)
    # ------------------------------------------------------------
    if rank == 0:
        with h5py.File(output_file, "a") as f:
            dset = f["ensemble"]
            ens_mean_ds = f["ensemble_mean"]

            # ---------------- index lookup for current recvbuf layout ----------------
            vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)

            thickness_idx = 0
            surface_idx = 0
            bed_idx = 0

            for ii, vec in enumerate(model_kwargs.get("vec_inputs", [])):
                vec_l = vec.lower()
                if vec_l in ["thickness", "ice_thickness", "h"]:
                    thickness_idx = indx_map[vec]
                if vec_l in ["surface", "ice_surface", "s"]:
                    surface_idx = indx_map[vec]
                if vec_l in ["bed", "bedrock", "base", "bedtopography"]:
                    bed_idx = indx_map[vec]

            # ---------------- bed relaxation ----------------
            if False: #TODO: for test runs
                for ii, vec in enumerate(model_kwargs.get("vec_inputs", [])):
                    vec_l = vec.lower()
                    if vec_l in ["bed", "bedrock", "base", "bedtopography"]:
                        bed_prior = dset[indx_map[vec], :, timestep - 1]
                        bed_now = recvbuf[indx_map[vec], :]

                        thickness = recvbuf[thickness_idx, :]
                        di = 0.8930
                        ocean_levelset = thickness + (recvbuf[bed_idx, :] / di)

                        dt = model_kwargs.get("dt", params["dt"])
                        t = timestep * dt

                        do_bed_snap = False
                        for bed_snap in model_kwargs.get("bed_obs_snapshot", []):
                            if np.isclose(t, bed_snap, rtol=0, atol=1e-12):
                                do_bed_snap = True
                                break
                        
                        do_bed_snap = False #TODO: for test runs
                        if do_bed_snap:
                            eta = 1.0
                            rho = model_kwargs.get("rho", 1.0)
                            sigma = 1e-3
                            X5 = model_kwargs.get("X5", None)
                            beta_t = model_kwargs.get("initial_bed_bias", 0.0015)

                            if X5 is not None:
                                for i in range(X5.shape[0]):
                                    for j in range(X5.shape[0]):
                                        beta_t *= X5[j, i]

                            for i_sig, sig in enumerate(params["sig_Q"]):
                                if i_sig == ii:
                                    sigma = sig

                            relaxation_factor = (eta + beta_t) * dt + np.sqrt(dt) * sigma * rho
                            if relaxation_factor > 1.5:
                                relaxation_factor = np.sqrt(dt) * sigma * rho
                            relaxation_factor = min(relaxation_factor, 0.5)

                            recvbuf[indx_map[vec], :] = bed_prior + relaxation_factor * (bed_now - bed_prior)
                        else:
                            relaxation_factor = model_kwargs.get("bed_relaxation_factor", 0.05)
                            recvbuf[indx_map[vec], :] = bed_prior + relaxation_factor * (bed_now - bed_prior)

            vec_inputs = model_kwargs.get("vec_inputs", [])

            if vec_inputs:
                hdim = recvbuf.shape[0] // len(vec_inputs)
                model_dt = model_kwargs.get("dt", params["dt"])
                current_model_time = timestep * model_dt

                recvbuf = apply_bed_domain_gate_global(
                    analysis_vec=recvbuf,
                    forecast_vec=forecast_recvbuf,
                    vec_inputs=vec_inputs,
                    hdim=hdim,
                    model_kwargs=model_kwargs,
                )
                recvbuf = apply_bed_observation_anchor_global(
                    analysis_vec=recvbuf,
                    vec_inputs=vec_inputs,
                    hdim=hdim,
                    model_kwargs=model_kwargs,
                    stage="pre",
                )
                if bed_idx is not None and forecast_recvbuf is not None:
                    model_kwargs["_bed_forecast_reference"] = np.asarray(
                        forecast_recvbuf[bed_idx, :], dtype=float
                    ).copy()

                recvbuf = apply_global_inference_hook(
                    analysis_vec=recvbuf,
                    vec_inputs=vec_inputs,
                    hdim=hdim,
                    params=params,
                    model_kwargs=model_kwargs,
                    timestep=timestep,
                    model_time=current_model_time,
                    stage="pre_geometry",  # bed only
                )

                # Spatial regularization can spread an increment across the
                # graph. Reapply the support gate so no smoothed correction
                # leaks beneath unobserved floating ice.
                recvbuf = apply_bed_domain_gate_global(
                    analysis_vec=recvbuf,
                    forecast_vec=forecast_recvbuf,
                    vec_inputs=vec_inputs,
                    hdim=hdim,
                    model_kwargs=model_kwargs,
                )
                # Re-anchor the actual survey nodes after graph smoothing. The
                # preceding gate still prevents inferred corrections beneath
                # unobserved floating ice.
                recvbuf = apply_bed_observation_anchor_global(
                    analysis_vec=recvbuf,
                    vec_inputs=vec_inputs,
                    hdim=hdim,
                    model_kwargs=model_kwargs,
                    stage="post",
                )
                if (
                    model_kwargs.get("physics_bed_inference", False)
                    and bed_idx is not None
                ):
                    model_kwargs["_bed_previous_applied"] = np.asarray(
                        recvbuf[bed_idx, :], dtype=float
                    ).copy()

            # ---------------- ISSM physical fixes ----------------
            if True:
                if model_kwargs.get("model_name", "").lower() == "issm":
                    di = 0.8930
                    rho_ice = 917.0
                    rho_sw = 1028.0

                    thickness = recvbuf[thickness_idx, :]
                    surface = recvbuf[surface_idx, :]
                    bed = recvbuf[bed_idx, :]

                    projection_mode = str(
                        model_kwargs.get(
                            "geometry_projection_mode", "preserve_thickness"
                        )
                    ).lower()
                    if projection_mode not in {
                        "preserve_thickness",
                        "preserve_surface",
                    }:
                        raise ValueError(
                            "geometry_projection_mode must be "
                            "'preserve_thickness' or 'preserve_surface'"
                        )

                    if projection_mode == "preserve_surface":
                        # Surface is the observed geometry variable in this
                        # profile. Diagnose H from S and b while retaining the
                        # forecast grounding classification for this projection.
                        grounded_before = thickness + (bed / di) >= 0.0
                        floating_before = ~grounded_before
                        thickness[grounded_before] = (
                            surface[grounded_before] - bed[grounded_before]
                        )
                        thickness[floating_before] = (
                            surface[floating_before]
                            * rho_sw
                            / (rho_sw - rho_ice)
                        )

                    pos = np.where(thickness < 1)
                    thickness[pos] = 1.0

                    ocean_levelset = thickness + (bed / di)

                    # floating ice
                    pos_float = np.where(ocean_levelset < 0)
                    surface[pos_float] = thickness[pos_float] * ((rho_sw - rho_ice) / rho_sw)

                    recvbuf[surface_idx, :] = surface
                    base = surface - thickness

                    pos_grounded = np.where(ocean_levelset > 0)
                    base[pos_grounded] = bed[pos_grounded]

                    recvbuf[surface_idx, :] = base + thickness
                    recvbuf[thickness_idx, :] = thickness

                    del thickness, surface, bed, ocean_levelset, pos, pos_float, base, pos_grounded
                    gc.collect()

            recvbuf = apply_global_inference_hook(
                analysis_vec=recvbuf,
                vec_inputs=vec_inputs,
                hdim=hdim,
                params=params,
                model_kwargs=model_kwargs,
                timestep=timestep,
                model_time=current_model_time,
                stage="post_geometry",  # SMB only
            )

            # Don't write yet if inversion is enabled. We'll do inversion first.
            if not model_kwargs.get("inversion_flag", False):
                dset[:, :, timestep] = recvbuf
                ens_mean_ds[:, timestep] = np.mean(recvbuf, axis=1)

                if model_kwargs.get("DEnKF_flag", False):
                    mean_now = np.mean(dset[:, :, timestep], axis=1)
                    dset[:, :, timestep] += mean_now[:, np.newaxis]

    # If no inversion, we are done after rank 0 writes
    if not model_kwargs.get("inversion_flag", False):
        comm.Barrier()
        return

    # ------------------------------------------------------------
    # Step 2: rank 0 prepares full inversion state
    # ------------------------------------------------------------
    if rank == 0:
        inv_kwargs_root = dict(model_kwargs)
        inv_kwargs_root["vec_inputs"] = copy.deepcopy(model_kwargs.get("vec_inputs_old", []))
        inv_kwargs_root["nd"] = model_kwargs.get("nd_old", None)

        vecs_inv, indx_map_inv, dim_per_proc_inv = icesee_get_index(**inv_kwargs_root)

        with h5py.File(
            f'{model_kwargs.get("data_path")}/ensemble_before_analysis_step_{timestep:04d}.h5',
            "a"
        ) as f_before:
            data_before = f_before["ensemble_before_analysis"]
            data_before_arr = data_before[:, :].copy()

        hdim = data_before_arr.shape[0] // len(inv_kwargs_root.get("vec_inputs", []))

        # inject updated analysis state from recvbuf into the full inversion state
        for ii, key in enumerate(model_kwargs.get("vec_inputs_new", [])):
            start = ii * hdim
            end = start + hdim
            data_before_arr[indx_map_inv[key], :] = recvbuf[start:end, :]

        full_shape = np.array(data_before_arr.shape, dtype=np.int32)
    else:
        inv_kwargs_root = None
        indx_map_inv = None
        data_before_arr = None
        full_shape = np.empty(2, dtype=np.int32)

    # broadcast shape
    comm.Bcast(full_shape, root=0)
    nd_full, nens_full = int(full_shape[0]), int(full_shape[1])

    # broadcast full state to all ranks
    if rank != 0:
        data_before_arr = np.empty((nd_full, nens_full), dtype=np.float64)
    comm.Bcast(data_before_arr, root=0)

    # ------------------------------------------------------------
    # Step 3: all ranks do member-wise inversion in parallel
    # ------------------------------------------------------------
    inv_kwargs = dict(model_kwargs)
    inv_kwargs["vec_inputs"] = copy.deepcopy(model_kwargs.get("vec_inputs_old", []))
    inv_kwargs["nd"] = model_kwargs.get("nd_old", None)

    # Each rank should treat each member independently
    # Using COMM_SELF is the safest non-breaking choice here
    inv_kwargs["comm"] = MPI.COMM_SELF

    vecs_inv, indx_map_inv, dim_per_proc_inv = icesee_get_index(**inv_kwargs)
    model_module = inv_kwargs.get("model_module", None)

    local_updates = []

    for ens_id in range(rank, Nens, size):
        member = data_before_arr[:, ens_id].copy()

        inv_kwargs_member = dict(inv_kwargs)
        inv_kwargs_member["ens_id"] = ens_id

        data = model_module.inverse_step_single(ensemble=member, **inv_kwargs_member)

        update_dict = {"ens_id": ens_id}

        for key, value in data.items():
            key_l = key.lower()

            if key_l in ["coefficient", "friction", "friction_coefficient", "fcoef", "frictioncoefficient"]:
                update_dict["friction_idx"] = indx_map_inv[key]
                update_dict["friction_val"] = np.asarray(value).copy()

            # elif key_l in ["vx", "velocity_x", "vel_x", "v_x"]:
            #     update_dict["vx_idx"] = indx_map_inv[key]
            #     update_dict["vx_val"] = np.asarray(value).copy()

            # elif key_l in ["vy", "velocity_y", "vel_y", "v_y"]:
            #     update_dict["vy_idx"] = indx_map_inv[key]
            #     update_dict["vy_val"] = np.asarray(value).copy()

        local_updates.append(update_dict)

        del member, data, inv_kwargs_member
        gc.collect()

    # Gather all updates to rank 0
    gathered_updates = comm.gather(local_updates, root=0)

    # ------------------------------------------------------------
    # Step 4: rank 0 applies inversion updates and writes to HDF5
    # ------------------------------------------------------------
    if rank == 0:
        for rank_updates in gathered_updates:
            for upd in rank_updates:
                ens_id = upd["ens_id"]

                if "friction_idx" in upd and "friction_val" in upd:
                    data_before_arr[upd["friction_idx"], ens_id] = upd["friction_val"]

                # if "vx_idx" in upd and "vx_val" in upd:
                    # data_before_arr[upd["vx_idx"], ens_id] = upd["vx_val"]

                # if "vy_idx" in upd and "vy_val" in upd:
                    # data_before_arr[upd["vy_idx"], ens_id] = upd["vy_val"]

        with h5py.File(output_file, "a") as f:
            dset = f["ensemble"]
            ens_mean_ds = f["ensemble_mean"]

            dset[:, :, timestep] = data_before_arr
            ens_mean_ds[:, timestep] = np.mean(data_before_arr, axis=1)

            if model_kwargs.get("DEnKF_flag", False):
                mean_now = np.mean(dset[:, :, timestep], axis=1)
                dset[:, :, timestep] += mean_now[:, np.newaxis]

        del data_before_arr
        gc.collect()

    comm.Barrier()
    
# # ---- Will uncomment above after fixing parallel i/o issues on the cluster ----
def parallel_write_ensemble_scattered_rank_0(timestep, ensemble_mean, params, ensemble_chunk, comm, model_kwargs, output_file="icesee_ensemble_data.h5"):
    """
    Write ensemble data using h5py and MPI, with only rank 0 writing to the dataset.
    Optimized for large datasets and many processes using MPI_Gatherv, without parallel I/O.
    ensemble_chunk: local data on each rank with shape (local_nd, Nens)
    """
    import numpy as np
    from mpi4py import MPI

    # MPI setup
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Get local chunk dimensions
    local_nd = ensemble_chunk.shape[0]  # rows per rank
    Nens = ensemble_chunk.shape[1]      # ensemble members (same for all ranks)

    # Gather the number of rows from each rank to rank 0
    local_nd_array = comm.gather(local_nd, root=0)
    
    if rank == 0:
        nd_total = sum(local_nd_array)
        # Prepare arrays for Gatherv
        counts = [n * Nens for n in local_nd_array]
        displs = [sum(counts[:i]) for i in range(size)]
        recvbuf = np.empty((nd_total, Nens), dtype=ensemble_chunk.dtype)
    else:
        nd_total = None
        counts = None
        displs = None
        recvbuf = None

    nd_total = comm.bcast(nd_total, root=0)

    # Gather ensemble chunks to rank 0
    comm.Gatherv(ensemble_chunk, [recvbuf, counts, displs, MPI.DOUBLE], root=0)

    output_file = os.path.join(model_kwargs.get('data_path'), output_file)

    if rank == 0:
        if timestep == 0.0:
            with h5py.File(output_file, 'w') as f:
                dset = f.create_dataset('ensemble', (nd_total, Nens, model_kwargs.get('nt', params['nt']) + 1), dtype=ensemble_chunk.dtype)
                ens_mean = f.create_dataset('ensemble_mean', (nd_total, model_kwargs.get('nt', params['nt']) + 1), dtype=ensemble_chunk.dtype)
                
                # Write gathered data
                dset[:, :, 0] = recvbuf
                ens_mean[:, 0] = ensemble_mean

                if model_kwargs.get("DEnKF_flag", False):
                    ensemble_mean = np.mean(dset[:, :, 0], axis=1)
                    dset[:, :, 0] += ensemble_mean[:, np.newaxis]
        else:
            with h5py.File(output_file, 'a') as f:
                dset = f['ensemble']
                ens_mean = f['ensemble_mean']

                # apply a relaxation factor to bed
                vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
                thickness_idx = 0; surface_idx = 0; bed_idx = 0
                for ii, vec in enumerate(model_kwargs.get("vec_inputs", [])):
                    # dumy variables to hold indices
                    # thickness_idx = 0; surface_idx = 0; bed_idx = 0
                    if vec.lower() in ["thickness","ice_thickness","h","Thickness"]:
                        thickness_idx = indx_map[vec]
                    # else:
                    #     thickness_idx = 0
                    if vec.lower() in ["surface","ice_surface","s","Surface"]:
                        surface_idx = indx_map[vec]
                    # else:
                    #     surface_idx = 0
                    if vec.lower() in ["Bed", "bed","bedrock","base","bedtopography"]:
                        bed_idx = indx_map[vec]
                    # else:
                    #     bed_idx = 0

                for ii, vec in enumerate(model_kwargs.get("vec_inputs", [])):
                    if vec.lower() in ["Bed", "bed","bedrock","base","bedtopography"]:
                        bed_prior = dset[indx_map[vec], :, timestep-1]
                        bed_now = recvbuf[indx_map[vec], :]

                        thickness = recvbuf[thickness_idx,:]
                        di = 0.8930
                        ocean_levelset = thickness + (recvbuf[bed_idx,:]/di)
                        pos_bed = np.where(ocean_levelset > 0) # no floating ice

                        # only update bed if observation is available
                        # for bed_snaps in model_kwargs.get("bed_obs_snapshot", []):
                        dt = model_kwargs.get("dt", params['dt'])
                            # if timestep*dt == bed_snaps:
                        t = timestep * dt
                        do_bed_snap = False
                        for bed_snap in model_kwargs.get("bed_obs_snapshot", []):
                            if np.isclose(t, bed_snap, rtol=0, atol=1e-12):
                                do_bed_snap = True
                                break
                            else:
                                do_bed_snap = False
                        if do_bed_snap:
                            # if True:
                            # relaxation_factor = 0.05
                            # relaxation_factor = (eta + beta_t)*dt + sqrt(dt)*sigma*rho*bed_error
                            eta = 1.0
                            rho = model_kwargs.get("rho", 1.0) 
                            sigma = 1e-3
                            X5 = model_kwargs.get("X5", None)
                            # beta_k = beta_0 \prod_{i=1}^{Nens} (X5_i)
                            beta_t = model_kwargs.get("initial_bed_bias", 0.0015)
                            for i in range(X5.shape[0]):
                                for j in range(X5.shape[0]):
                                    beta_t *= X5[j,i]
                            for i, sig in enumerate(params["sig_Q"]):
                                if i == ii:
                                    sigma = sig
                            relaxation_factor = (eta+beta_t)*dt + np.sqrt(dt)*sigma*rho
                            # put cap on relaxation factor to avoid instability 
                            if relaxation_factor > 1.5:
                                relaxation_factor = np.sqrt(dt)*sigma*rho
                            relaxation_factor = min(relaxation_factor, 0.5)
                            recvbuf[indx_map[vec], :] = bed_prior + relaxation_factor * (bed_now - bed_prior)

                            # update bed bias
                            # model_kwargs["initial_bed_bias"] = beta_t  
                        else:
                            # do_bed_snap = False
                            relaxation_factor = model_kwargs.get("bed_relaxation_factor", 0.05)
                            recvbuf[indx_map[vec], :] = bed_prior + relaxation_factor * (bed_now - bed_prior)
                            # recvbuf[pos_bed, :] = bed_prior[pos_bed]

                    
                # check for negative thickness
                # ISSM *------
                if model_kwargs.get("model_name", "").lower() == "issm":
                    di = 0.8930
                    rho_ice = 917.0
                    rho_sw = 1028.0

                    thickness = recvbuf[thickness_idx, :]
                    surface   = recvbuf[surface_idx, :]
                    bed       = recvbuf[bed_idx, :]

                    projection_mode = str(
                        model_kwargs.get(
                            "geometry_projection_mode", "preserve_thickness"
                        )
                    ).lower()
                    if projection_mode not in {
                        "preserve_thickness",
                        "preserve_surface",
                    }:
                        raise ValueError(
                            "geometry_projection_mode must be "
                            "'preserve_thickness' or 'preserve_surface'"
                        )

                    if projection_mode == "preserve_surface":
                        grounded_before = thickness + bed / di >= 0.0
                        floating_before = ~grounded_before
                        thickness[grounded_before] = (
                            surface[grounded_before] - bed[grounded_before]
                        )
                        thickness[floating_before] = (
                            surface[floating_before]
                            * rho_sw
                            / (rho_sw - rho_ice)
                        )

                    thickness[thickness < 1.0] = 1.0

                    ocean_levelset = thickness + bed / di

                    pos_float = np.where(ocean_levelset < 0)
                    surface[pos_float] = thickness[pos_float] * ((rho_sw - rho_ice) / rho_sw)

                    base = surface - thickness

                    pos_base = np.where(base < bed)
                    base[pos_base] = bed[pos_base]

                    pos_grounded = np.where(ocean_levelset > 0)
                    base[pos_grounded] = bed[pos_grounded]

                    surface = base + thickness

                    recvbuf[thickness_idx, :] = thickness
                    recvbuf[surface_idx, :] = surface

                    # -------*ISSM
                    del thickness, surface, bed, ocean_levelset, pos, base, pos_base, pos_grounded
                    gc.collect()

                    # get velcity and friction from inversion -------------
                    # mean_now = np.mean(recvbuf, axis=1)
                    if model_kwargs.get("inversion_flag", False):
                        # update model_kwargs for inversion
                        model_kwargs["vec_inputs"] = copy.deepcopy(model_kwargs.get("vec_inputs_old", []))
                        model_kwargs["nd"] = model_kwargs.get("nd_old", None)
                        vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
                        with h5py.File(f'{model_kwargs.get("data_path")}/ensemble_before_analysis_step_{timestep:04d}.h5', 'a') as f_before:
                            data_before = f_before['ensemble_before_analysis']
                            data_before_arr = data_before[:, :].copy()

                            hdim = data_before_arr.shape[0] // len(model_kwargs.get("vec_inputs", []))
                        
                            # update data_before with recvbuf
                            for ii, key in enumerate (model_kwargs.get("vec_inputs_new", [])):
                                # print('key:\n', key)
                                start = ii*hdim
                                end = start + hdim      
                                data_before_arr[indx_map[key], :] = recvbuf[start:end, :]
                                # data_before[indx_map[key], :] = recvbuf[indx_map[key], :]
                            #TODO: test this part
                            # recvbuf[:, :] =  dset[:, :, timestep-1]
                            # recompute the mean_now after updating friction
                            # mean_now = np.mean(recvbuf, axis=1)
                            mean_now = np.mean(data_before_arr, axis=1)
                            # call inverse step to get velocity fields, and new friction
                            model_module   = model_kwargs.get("model_module", None)
                            data = model_module.inverse_step_single(ensemble=mean_now, **model_kwargs)
                            for key, value in data.items():
                                # if key.lower() in ["vx","velocity_x","vel_x","v_x"]:
                                #     anomaly = recvbuf[indx_map[key], :] - value[:, np.newaxis]
                                #     velocity_x_prior = dset[indx_map[key], :, timestep-1]
                                #     recvbuf[indx_map[key], :] = value[:, np.newaxis]
                                # if key.lower() in ["vy","velocity_y","vel_y","v_y"]:
                                #     velocity_y_prior = dset[indx_map[key], :, timestep-1]
                                #     anomaly = recvbuf[indx_map[key], :] - value[:, np.newaxis]
                                    # recvbuf[indx_map[key], :] = value[:, np.newaxis]
                                if key.lower() in ["coefficient","friction","friction_coefficient", 'fcoef', "frictioncoefficient"]:
                                    data_before_arr[indx_map[key], :] = value[:, np.newaxis]
                            del mean_now
                            gc.collect()

                        dset[:, :, timestep] = data_before_arr
                        ens_mean[:, timestep] = np.mean(data_before_arr, axis=1)
                    
                    else:
                        dset[:, :, timestep] = recvbuf
                        ens_mean[:, timestep] = np.mean(recvbuf, axis=1)

                    # dset[:, :, timestep] = data_before
                    # ens_mean[:, timestep] = np.mean(data_before, axis=1)
                    del bed_prior, bed_now
                    gc.collect()
                else:
                    dset[:, :, timestep] = recvbuf
                    ens_mean[:, timestep] = np.mean(recvbuf, axis=1)

                if model_kwargs.get("DEnKF_flag", False):
                    ensemble_mean = np.mean(dset[:, :, timestep], axis=1)
                    dset[:, :, timestep] += ensemble_mean[:, np.newaxis]

    comm.Barrier()

def parallel_write_data_from_root_2D(
    full_ensemble=None,
    comm=None,
    data_name=None,
    output_file="preliminary_data.h5",
    data_path="_modelrun_datasets",
):
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
    if rank == 0:
        os.makedirs(data_path, exist_ok=True)
    comm.Barrier()
    output_file = os.path.join(data_path, output_file)

    # Open file in parallel mode
    # with h5py.File(output_file, 'w', driver='mpio', comm=comm) as f:
    with h5py.File(output_file, 'w') as f:
        # Create dataset with total dimensions
        dset = f.create_dataset(data_name, (nd, Nens), dtype=dtype)
        
        # Each rank writes its chunk
        dset[offset:offset + local_nd, :] = local_chunk

def parallel_write_vector_from_root(
    full_ensemble=None,
    comm=None,
    data_shape=None,
    data_name=None,
    output_file="icesee_ensemble_data.h5",
    data_path="_modelrun_datasets",
):
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
    if rank == 0:
        os.makedirs(data_path, exist_ok=True)
    comm.Barrier()
    output_file = os.path.join(data_path, output_file)

    # Open file in parallel mode
    # with h5py.File(output_file, 'w', driver='mpio', comm=comm) as f:
    with h5py.File(output_file, 'w') as f:
        # Create dataset with total dimensions
        #  data_shape should be a tuple
        dset = f.create_dataset(data_name, data_shape, dtype=dtype)
        
        # Each rank writes its chunk
        dset[offset:offset + local_nd, :,0] = local_chunk

    comm.Barrier()

# ---- Will uncomment above after fixing parallel i/o issues on the cluster ----
def parallel_write_full_ensemble_from_root(timestep, ensemble_mean, model_kwargs, full_ensemble=None, comm=None, output_file="icesee_ensemble_data.h5"):
    """
    Append ensemble data where the full matrix exists on rank 0, with only rank 0 writing to the dataset.
    Optimized for large datasets and many processes without parallel I/O.
    Each call appends a new time step, resulting in a dataset of shape (nd, Nens, nt).
    full_ensemble: complete matrix on rank 0 with shape (nd, Nens)
    comm: MPI communicator
    output_file: Name of the output HDF5 file
    """
    import numpy as np
    import h5py

    params = model_kwargs.get("params")

    # MPI setup
    rank = comm.Get_rank()

    # Get dimensions on root and broadcast
    if rank == 0:
        nd, Nens = full_ensemble.shape
        dtype = full_ensemble.dtype
    else:
        nd, Nens, dtype = None, None, None
    
    nd = comm.bcast(nd, root=0)
    Nens = comm.bcast(Nens, root=0)
    dtype = comm.bcast(dtype, root=0)

    # Define output file path
    output_file = os.path.join(params.get('data_path'), output_file)

    # Only rank 0 writes to the file
    if rank == 0:
        if timestep == 0:
            with h5py.File(output_file, 'w') as f:
                # Create dataset with total dimensions
                dset = f.create_dataset('ensemble', (nd, Nens, model_kwargs.get('nt', params['nt']) + 1), dtype=dtype)
                # Write full ensemble
                dset[:, :, 0] = full_ensemble[:,:Nens]

                # Create and write ensemble mean
                ens_mean = f.create_dataset('ensemble_mean', (nd, model_kwargs.get('nt', params['nt']) + 1), dtype=dtype)
                ens_mean[:, 0] = ensemble_mean
                # ens_mean[:, 0] = full_ensemble[:, 0]

                if model_kwargs.get("DEnKF_flag", False):
                    ensemble_mean = np.mean(dset[:, :, 0], axis=1)
                    dset[:, :, 0] += ensemble_mean[:, np.newaxis]
        else:
            with h5py.File(output_file, 'a') as f:
                dset = f['ensemble']
                # Write full ensemble for current timestep
                dset[:, :, timestep] = full_ensemble[:,:Nens]

                ens_mean = f['ensemble_mean']
                ens_mean[:, timestep] = ensemble_mean

                if model_kwargs.get("DEnKF_flag", False):
                    ensemble_mean = np.mean(dset[:, :, timestep], axis=1)
                    dset[:, :, timestep] += ensemble_mean[:, np.newaxis]

    comm.Barrier()

def parallel_write_full_ensemble_from_root_full_parallel_run(timestep, ensemble_mean, model_kwargs, full_ensemble=None, comm=None, output_file="icesee_ensemble_data.h5"):
        """
        Append ensemble data where the full matrix exists on rank 0, with only rank 0 writing to the dataset.
        Optimized for large datasets and many processes without parallel I/O.
        Each call appends a new time step, resulting in a dataset of shape (nd, Nens, nt).
        full_ensemble: complete matrix on rank 0 with shape (nd, Nens)
        comm: MPI communicator
        output_file: Name of the output HDF5 file
        """
        import numpy as np
        import h5py

        params = model_kwargs.get("params")

        # MPI setup
        rank = comm.Get_rank()

        # Get dimensions on root and broadcast
        if rank == 0:
            nd, Nens = full_ensemble.shape
            dtype = full_ensemble.dtype
        else:
            nd, Nens, dtype = None, None, None
        
        nd = comm.bcast(nd, root=0)
        Nens = comm.bcast(Nens, root=0)
        dtype = comm.bcast(dtype, root=0)

        # Define output file path
        output_file = os.path.join(params.get('data_path'), output_file)

        # Only rank 0 writes to the file
        if rank == 0:
            if timestep == 0:
                with h5py.File(output_file, 'w') as f:
                    # Create dataset with total dimensions
                    # dset = f.create_dataset('ensemble', (nd, Nens, model_kwargs.get('nt', params['nt']) + 1), dtype=dtype)
                    chunk_size = (min(5000, nd), 1)
                    dset = f.create_dataset('ensemble', (nd, Nens), dtype=dtype, chunks=chunk_size, compression="gzip", compression_opts=9)
                    # Write full ensemble
                    dset[:, :, 0] = full_ensemble

                    # Create and write ensemble mean
                    ens_mean = f.create_dataset('ensemble_mean', (nd, model_kwargs.get('nt', params['nt']) + 1), dtype=dtype)
                    ens_mean[:, 0] = ensemble_mean
                    # ens_mean[:, 0] = full_ensemble[:, 0]

                    if model_kwargs.get("DEnKF_flag", False):
                        ensemble_mean = np.mean(dset[:, :, 0], axis=1)
                        dset[:, :, 0] += ensemble_mean[:, np.newaxis]
            else:
                with h5py.File(output_file, 'a') as f:
                    dset = f['ensemble']
                    # Write full ensemble for current timestep
                    dset[:, :, timestep] = full_ensemble

                    ens_mean = f['ensemble_mean']
                    ens_mean[:, timestep] = ensemble_mean

                    if model_kwargs.get("DEnKF_flag", False):
                        ensemble_mean = np.mean(dset[:, :, timestep], axis=1)
                        dset[:, :, timestep] += ensemble_mean[:, np.newaxis]

        comm.Barrier()

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


def compute_HA_block(h5_path, timestep, obs_indices, Nens, chunk_members=1):
    """
    Evensen (2003) Sec 5.4 block algorithm: build HA (m_obs, Nens) by
    reading one ensemble member (or small chunk) from disk at a time,
    projecting into observation space (row-select via obs_indices, since
    H is a one-hot selector — see H_indices), and discarding the
    full-length (nd,) column immediately. Never holds the full
    (nd, Nens) ensemble in memory — only the much smaller (m_obs, Nens)
    result accumulates.

    h5_path   : path to the ensemble HDF5 file (the 'ensemble' dataset,
                shape (nd, Nens, nt+1), as written by
                parallel_write_ensemble_scattered)
    timestep  : which time slice to read (the forecast state before this
                analysis step)
    obs_indices : (m_obs,) global row indices, from H_indices/JObs_indices
    chunk_members : how many members to read per HDF5 call (1 = literal
                    Evensen "one member at a time"; higher trades memory
                    for fewer I/O calls)
    """
    import h5py
    import numpy as np

    with h5py.File(h5_path, "r") as f:
        dset = f["ensemble"]
        m_obs = obs_indices.size
        HA = np.empty((m_obs, Nens), dtype=np.float64)

        for start in range(0, Nens, chunk_members):
            end = min(start + chunk_members, Nens)
            member_chunk = dset[:, start:end, timestep]        # (nd, chunk)
            HA[:, start:end] = member_chunk[obs_indices, :]    # (m_obs, chunk)
            del member_chunk

    return HA


def compute_HAprime_block(h5_path, timestep, obs_indices, Nens, chunk_members=1):
    """
    Block version of HAprime = H @ (ensemble - mean) = HA - HAbar.
    Since H is linear (one-hot selector), H @ mean == mean(HA, axis=1),
    so this is exact, not an approximation.
    """
    HA = compute_HA_block(h5_path, timestep, obs_indices, Nens, chunk_members)
    HAbar = HA.mean(axis=1, keepdims=True)
    return HA - HAbar, HA  # return HA too — needed for Dprime below


def compute_Dprime_block(h5_path, timestep, obs_indices, d, Nens, chunk_members=1):
    """
    Block version of Dprime = d.reshape(-1,1) - H @ ensemble_vec
                              = d.reshape(-1,1) - HA.
    Reuses the same streamed HA rather than re-reading the file.
    """
    HA = compute_HA_block(h5_path, timestep, obs_indices, Nens, chunk_members)
    return d.reshape(-1, 1) - HA

def partition_rows(nd, size, rank):
    """Contiguous row partition, consistent with the rest of the partitioned path."""
    q, r = divmod(nd, size)
    start = rank * q + min(rank, r)
    stop = start + q + (1 if rank < r else 0)
    return start, stop


def write_ensemble_member_direct(h5_path, timestep, ens_id, member_vec, nd, Nens, nt, comm):
    """Each subcomm leader writes its own completed member directly —
    no Gatherv, no full-ensemble array anywhere. Ranks with no completed
    member this call (ens_id is None) still participate in the
    collective open/write with an empty selection."""
    import h5py
    import h5py.h5p as h5p, h5py.h5s as h5s, h5py.h5fd as h5fd
    import numpy as np

    with h5py.File(h5_path, "a", driver="mpio", comm=comm) as f:
        if "ensemble" not in f:
            f.create_dataset("ensemble", (nd, Nens, nt + 1), dtype="f8")
        dset = f["ensemble"]

        dxpl = h5p.create(h5p.DATASET_XFER)
        dxpl.set_dxpl_mpio(h5fd.MPIO_COLLECTIVE)

        file_space = dset.id.get_space()
        if ens_id is not None:
            file_space.select_hyperslab((0, ens_id, timestep), (nd, 1, 1))
            mem_space = h5s.create_simple((nd,))
            buf = np.ascontiguousarray(member_vec, dtype=np.float64)
        else:
            mem_space = h5s.create_simple((0,))
            file_space.select_none()
            buf = np.empty((0,), dtype=np.float64)

        dset.id.write(mem_space, file_space, buf, dxpl=dxpl)

    comm.Barrier()


def compute_and_apply_inflation_partitioned(h5_path, nd, Nens, alpha, comm, timestep=0):
    """
    Partial-parallel-native version of the mean+inflate-at-init step.
    Two-pass, row-partitioned, no full ensemble anywhere:
      Pass 1: each rank sums its own row-block across all Nens members
              (Allreduce not needed here — each rank already has ONLY
              its own rows; no cross-rank sum required since row-blocks
              are disjoint).
      Pass 2: each rank re-reads its row-block, applies
              mean + alpha*(member - mean), writes back.
    This mirrors compute_forecast_mean_chunked_v2's SPIRIT (running sum,
    then per-rank write) but is self-contained here — no dependency on
    EnKF_fully_parallel_IO.
    """
    import h5py
    import numpy as np

    rank = comm.Get_rank()
    size = comm.Get_size()
    row_start, row_stop = partition_rows(nd, size, rank)
    local_nd = row_stop - row_start

    with h5py.File(h5_path, "a", driver="mpio", comm=comm) as f:
        dset = f["ensemble"]

        if local_nd > 0:
            block = dset[row_start:row_stop, :, timestep]        # (local_nd, Nens)
            mean_local = block.mean(axis=1, keepdims=True)
            inflated = mean_local + alpha * (block - mean_local)
            dset[row_start:row_stop, :, timestep] = inflated

    comm.Barrier()

def compute_HAprime_Eta_Dprime_partitioned(
    h5_path,
    timestep,
    obs_indices,
    d,
    Nens,
    comm,
    sigma,
    seed,
    error_mode="stochastic_r",
):
    """Row-partitioned, Allreduce-based HAprime/Eta/Dprime. No rank ever
    holds the full (nd, Nens) ensemble."""
    import h5py
    from mpi4py import MPI
    import numpy as np

    rank = comm.Get_rank()
    size = comm.Get_size()

    with h5py.File(h5_path, "r", driver="mpio", comm=comm) as f:
        dset = f["ensemble"]
        nd = dset.shape[0]
        row_start, row_stop = partition_rows(nd, size, rank)
        local_nd = row_stop - row_start
        States_local = dset[row_start:row_stop, :, timestep] if local_nd > 0 else np.empty((0, Nens))

    m_obs = obs_indices.size
    HA_local = np.zeros((m_obs, Nens), dtype=np.float64)

    if local_nd > 0:
        in_range = (obs_indices >= row_start) & (obs_indices < row_stop)
        rows_here = np.nonzero(in_range)[0]
        local_rows = obs_indices[in_range] - row_start
        if rows_here.size > 0:
            HA_local[rows_here, :] = States_local[local_rows, :]

    HA = np.empty_like(HA_local)
    comm.Allreduce(HA_local, HA, op=MPI.SUM)

    error_mode = str(error_mode).lower()
    if error_mode == "stochastic_r":
        from ICESEE.src.utils.localization import stochastic_observation_terms

        HAprime, Eta, Dprime = stochastic_observation_terms(
            HA, d, sigma, seed
        )
    elif error_mode == "legacy_prior_anomalies":
        HAprime = HA - HA.mean(axis=1, keepdims=True)
        Eta = HAprime.copy()
        Dprime = d.reshape(-1, 1) - HA
    else:
        raise ValueError(
            "enkf_observation_error_mode must be 'stochastic_R' or "
            "'legacy_prior_anomalies'"
        )

    return HAprime, Eta, Dprime


def write_analysis_partitioned(k, X5, local_patches, h5_path, timestep_forecast, model_kwargs, comm):
    """Row-partitioned analysis write: inflation + bed-freeze + local
    patches applied per-rank on its own row-block only. No full
    ensemble anywhere."""
    import h5py
    import numpy as np
    from ICESEE.src.utils.localization import (
        apply_local_patches,
        restore_frozen_analysis_vars,
    )

    params = model_kwargs.get("params")
    rank = comm.Get_rank()
    size = comm.Get_size()

    with h5py.File(h5_path, "a", driver="mpio", comm=comm) as f:
        dset = f["ensemble"]
        nd, Nens, _ = dset.shape
        row_start, row_stop = partition_rows(nd, size, rank)
        local_nd = row_stop - row_start
        global_rows = row_start + np.arange(local_nd)

        if local_nd > 0:
            prior_local = dset[row_start:row_stop, :, timestep_forecast]
            analysis_local = prior_local @ X5

            if model_kwargs.get("local_analysis", False):
                analysis_local = apply_local_patches(analysis_local, prior_local, global_rows, local_patches)

            vec_inputs = model_kwargs.get("vec_inputs", [])
            hdim = nd // len(vec_inputs)
            state_inflation = params.get("state_inflation_factor", params.get("inflation_factor", 1.0))
            param_inflation = params.get("param_inflation_factor", params.get("inflation_factor", 1.0))
            bed_inflation = params.get("bed_inflation_factor", param_inflation)

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

            local_mean = np.mean(analysis_local, axis=1, keepdims=True)
            local_pert = analysis_local - local_mean
            analysis_local = local_mean + inflation_vec[:, None] * local_pert
            analysis_local = restore_frozen_analysis_vars(
                analysis_local,
                prior_local,
                global_rows,
                vec_inputs,
                hdim,
                model_kwargs.get("frozen_analysis_vars", []),
            )

            km = model_kwargs.get("km", None)
            bed_snap_cols = set(model_kwargs.get("bed_snap_cols", []))
            if km is not None and int(km) not in bed_snap_cols:
                for ii, key in enumerate(vec_inputs):
                    if key.lower() not in ["bed", "bedrock", "bedtopography", "bed_topography"]:
                        continue
                    start, end = ii * hdim, ii * hdim + hdim
                    local_mask = (global_rows >= start) & (global_rows < end)
                    if np.any(local_mask):
                        analysis_local[local_mask, :] = prior_local[local_mask, :]

            dset[row_start:row_stop, :, k + 1] = analysis_local

    comm.Barrier()


def write_ensemble_member_direct_h5(dset, timestep, ens_id, member_vec, nd, comm):
    """
    Write ONE member's column using an ALREADY-OPEN dataset handle.
    No file open/close here — caller manages the file lifecycle.
    Still a collective call (every rank must call this together, even
    with ens_id=None/member_vec=None) since the dataset uses mpio.
    """
    import h5py.h5p as h5p, h5py.h5s as h5s, h5py.h5fd as h5fd
    import numpy as np

    dxpl = h5p.create(h5p.DATASET_XFER)
    dxpl.set_dxpl_mpio(h5fd.MPIO_COLLECTIVE)

    file_space = dset.id.get_space()
    if ens_id is not None:
        file_space.select_hyperslab((0, ens_id, timestep), (nd, 1, 1))
        mem_space = h5s.create_simple((nd,))
        buf = np.ascontiguousarray(member_vec, dtype=np.float64)
    else:
        mem_space = h5s.create_simple((0,))
        file_space.select_none()
        buf = np.empty((0,), dtype=np.float64)

    dset.id.write(mem_space, file_space, buf, dxpl=dxpl)


def open_ensemble_file(h5_path, nd, Nens, nt, comm, mode="a"):
    """
    Open the shared ensemble file ONCE, creating the dataset if needed.
    Caller is responsible for closing it (use as a context manager or
    call .close() explicitly) — reused across many collective writes
    within a timestep/round-loop to avoid repeated mpio file-open cost.
    """
    import h5py
    f = h5py.File(h5_path, mode, driver="mpio", comm=comm)
    if "ensemble" not in f:
        f.create_dataset("ensemble", (nd, Nens, nt + 1), dtype="f8")
    return f
