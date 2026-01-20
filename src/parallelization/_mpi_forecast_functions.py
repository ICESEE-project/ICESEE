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

from ICESEE.src.utils.tools import icesee_get_index
# from ICESEE.src.run_model_da._parallel_i_o import parallel_write_full_ensemble_from_root
                                               
from ICESEE.src.run_model_da._error_generation import generate_enkf_field
from ICESEE.src.utils.utils import UtilsFunctions

from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager
# rank_seed, rng = ParallelManager().initialize_seed(MPI.COMM_WORLD)

def parallel_forecast_step_default_run(**model_kwargs):
    """
    Process ensemble data in parallel using MPI, with scalable communication via MPI_Gatherv.
    Only rank 0 in each subcommunicator and the global communicator gathers and processes results.
    """
    import numpy as np
    import h5py
    from mpi4py import MPI

    # --- cases 2 & 3 ---
    # case 2: Form batches of sub-communicators and distribute resources among them
    #          - only works for Nens >= size_world
    # case 3: form Nens sub-communicators and distribute resources among them
    #          - only works for size_world > Nens
    #          - even distribution and load balancing leading to performance improvement
    #          - best for size_world/Nens is a whole number

    # Get necessary parameters
    params                    = model_kwargs.get("params")
    Nens                      = params.get("Nens")
    rounds                    = model_kwargs.get("rounds")
    color                     = model_kwargs.get("color", 0)                    # Default color for the subcommunicator
    subcomm_size_min              = model_kwargs.get("subcomm_size_min", 1)            # Size of the subcommunicator
    subcomm                   = model_kwargs.get("subcomm", None)              # Subcommunicator for local ensemble processing
    comm_world                = model_kwargs.get("comm_world", MPI.COMM_WORLD) # Global communicator
    rank_world                = comm_world.Get_rank()                          # Global rank
    sub_rank                  = subcomm.Get_rank() if subcomm else 0          # Rank within the subcommunicator
    _modelrun_datasets        = model_kwargs.get("_modelrun_datasets", "_modelrun_datasets")
    alpha                     = model_kwargs.get("alpha", 0.0)                 # Default alpha value
    rho                       = model_kwargs.get("rho", 1.0)                   # Default rho value
    dt                        = model_kwargs.get("dt", 1.0)                    # Default time step
    Lx                        = model_kwargs.get("Lx", 1.0)                    # Default domain length in x-direction
    Ly                        = model_kwargs.get("Ly", 1.0)                    # Default domain length in y-direction
    len_scale                 = model_kwargs.get("len_scale", 1.0)             # Default length scale for noise generation
    model_module              = model_kwargs.get("model_module", None)         # Module containing the forecast step function
    k                         = model_kwargs.get("k", 0)                      # Time step index, default to 0 if not provided
    noise                     = model_kwargs.get("noise", None)                   # Noise vector, if provided
    rng                       = model_kwargs.get("rng", np.random.default_rng())  # Random number generator, default to numpy's default RNG
    rank_seed = model_kwargs.get("rank_seed", 0)

    # np.random.seed(rank_seed)

    size_world = comm_world.Get_size()  # Total number of ranks in the global communicator

    # get timing variables from model_kwargs 
    time_forecast_ensemble_generation = model_kwargs.get("time_forecast_ensemble_generation", 0.0)
    time_forecast_noise_generation = model_kwargs.get("time_forecast_noise_generation", 0.0)
    time_forecast_ensemble_mean_generation = model_kwargs.get("time_forecast_ensemble_mean_generation", 0.0)
    time_forecast_file_writing = model_kwargs.get("time_forecast_file_writing", 0.0)
                         
    # --- case 2: Form batches of sub-communicators and distribute resources among them ---
    if Nens >= size_world:
        # ensemble_vec, shape_ens = parallel_forecast_step_default_run(**model_kwargs)
        
        # store results for each round      
        ens_list = []
        for round_id in range(rounds):
            ensemble_id = color + round_id * subcomm_size_min  # Global ensemble index
            model_kwargs.update({'ens_id': ensemble_id, 'comm': subcomm})

            if ensemble_id < Nens:  # Only process valid ensembles
                # print(f"[ICESEE] Rank {rank_world} processing ensemble {ensemble_id} in round {round_id + 1}/{rounds}")

                # Ensure all ranks in the subcommunicator are synchronized before running
                subcomm.Barrier()
                ens = ensemble_id
                # ---- read from file ----
                input_file = f"{_modelrun_datasets}/icesee_ensemble_data.h5"
                # with h5py.File(input_file, "r", driver="mpio", comm=subcomm) as f:
                # with h5py.File(input_file, "r", driver="mpio", comm=comm_world) as f:
                with h5py.File(input_file, "r") as f:
                    ensemble_vec = f["ensemble"][:,ens,k]
                
                # print(f"[ICESEE] Finished reading ensemble {ens} ensemble shape: {ensemble_vec.shape} norm {np.linalg.norm(ensemble_vec)}")
                # ---- end of read from file ----

                # Call the forecast step function
                # hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
                # print(f"[ICESEE] Rank: {rank_world}, min ensemble: {np.min(ensemble_vec[:hdim])}, max ensemble: {np.max(ensemble_vec[:hdim])}")

                updated_state = model_module.forecast_step_single(ensemble=ensemble_vec,**model_kwargs)

                #fetch the updated state
                vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
                for key,value in updated_state.items():
                    ensemble_vec[indx_map[key]] = value

                #  add time evolution noise to the ensemble
                # if k == 0:
                #     # model noise, q0
                #     q0 = np.random.multivariate_normal(np.zeros(nd), Q_err)

                # # squence of white noise drawn from  a smooth pseudorandm fields,w0,
                # # w0 = np.random.normal(0, 1, nd) #TODO: look into this
                # # w0 = np.random.multivariate_normal(np.zeros(nd), np.eye(nd))
                # # q0 = alpha * q0 + np.sqrt(1 - alpha**2) * w0
                # q0 = np.random.multivariate_normal(np.zeros(nd), Q_err)
                # Q_err = Q_err[:state_block_size,:state_block_size]
                # q0 = multivariate_normal.rvs(np.zeros(state_block_size), Q_err)
                # q0 = np.sqrt(model_kwargs.get("dt",params["dt"]))*multivariate_normal.rvs(np.zeros(state_block_size), Q_err)

                # if k+1 <= max(model_kwargs["obs_index"]):
                #     ensemble_vec[:state_block_size] = ensemble_vec[:state_block_size] + q0[:state_block_size]
                # else:
                    #  create a guassian noise with zero mean and variance = 1
                    # q0 = np.random.normal(0, 1, state_block_size)
                    # ensemble_vec[:state_block_size] = ensemble_vec[:state_block_size] + q0[:state_block_size]
                #---------------------------------------------------------------
                if model_kwargs["joint_estimation"] or params["localization_flag"]:
                    hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
                else:
                    hdim = ensemble_vec.shape[0] // params["num_state_vars"]
                state_block_size = hdim * params["num_state_vars"]

                # --- time forecast noise generation ---
                _time_forecast_noise_generation = MPI.Wtime()
                # if k == 0:
                #     # noise = compute_noise_random_fields(ens, hdim, pos, gs_model, params["total_state_param_vars"], L_C)
                #     N_size = params["total_state_param_vars"] * hdim
                #     # noise = generate_pseudo_random_field_1d(N_size,np.sqrt(Lx*Ly), len_scale, verbose=0)
                #     noise = generate_enkf_field(None,np.sqrt(Lx*Ly),hdim, params["total_state_param_vars"], rh=len_scale, rng=rng, verbose=False)

                # noise = noise / np.max(np.abs(noise))
                # if k+1 <= max(model_kwargs["obs_index"]):
                    # W = np.random.normal(0, 1, state_block_size)

                # ======
                noise_all = []
                q0 = []
                for ii, sig in enumerate(params["sig_Q"]):
                    if ii <=params["num_state_vars"]:
                        # W = np.random.normal(0, 1, hdim)
                        # W = generate_pseudo_random_field_1d(hdim,np.sqrt(Lx*Ly), len_scale, verbose=0)
                        model_kwargs.update({"ii_sig": ii, "hdim":hdim, "num_vars":params["total_state_param_vars"]})
                        # W = generate_enkf_field(**model_kwargs)
                        W = generate_enkf_field(ii,np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)
                        noise_ = alpha*noise[ii*hdim:(ii+1)*hdim] + np.sqrt(1 - alpha**2)*W
                        q0.append(noise_)

                        Z = np.sqrt(dt)*sig*rho*noise_
                        noise_all.append(Z)
                noise_ = np.concatenate(noise_all, axis=0)
                ensemble_vec[:state_block_size] = ensemble_vec[:state_block_size] + noise_[:state_block_size]
                # ensemble_vec[hdim:state_block_size] = ensemble_vec[hdim:state_block_size] + noise_[hdim:state_block_size]

                # for ii, key in enumerate(model_kwargs['observed_vars']):
                #         if ii < params["num_state_vars"]:
                #             ensemble_vec[indx_map[key]] += noise_[indx_map[key]]
                noise = np.concatenate(q0, axis=0)
                model_kwargs.update({"noise": noise})  # save the noise to the model_kwargs dictionary

                # clean up memory
                del noise_all, q0, noise_, W
                time_forecast_noise_generation += MPI.Wtime() - _time_forecast_noise_generation

                # =====
                # pack
                
                
                # mean_x = np.mean(ensemble_vec[:state_block_size], axis=1)[:,np.newaxis]
                # ensemble_vec[:state_block_size] = ensemble_vec[:state_block_size] - mean_x
                # Ensure all ranks in the subcommunicator are synchronized before moving on
                # subcomm.Barrier()

                # Gather results within each subcommunicator
                # gathered_ensemble = subcomm.gather(ensemble_vec[:], root=0)
                # replace with GatherV
                local_shape = ensemble_vec.size
                local_shapes = subcomm.gather(local_shape, root=0)

                if sub_rank == 0:
                    total_size = sum(local_shapes)
                    counts = local_shapes
                    displs = [sum(counts[:i]) for i in range(subcomm.Get_size())]
                    gathered_ensemble = np.empty(total_size, dtype=ensemble_vec.dtype)
                else:
                    gathered_ensemble = None
                    counts = None
                    displs = None

                subcomm.Gatherv(ensemble_vec[:], [gathered_ensemble, counts, displs, MPI.DOUBLE], root=0)

                # Ensure only rank = 0 in each subcommunicator gathers the results
                if sub_rank == 0:
                        gathered_ensemble = np.hstack(gathered_ensemble)

                ens_list.append(gathered_ensemble if sub_rank == 0 else None)

            # Gather results from all subcommunicators
            gathered_ensemble_global = ParallelManager().gather_data(comm_world, ens_list, root=0)
            # gathered_ensemble_global = comm_world.gather(ens_list, root=0)
            # replace with GatherV
            # local_list_shape = np.array(ens_list).size
            # local_list_shapes = comm_world.gather(local_list_shape, root=0)
            # if rank_world == 0:
            #     total_size_global = sum(local_list_shapes)
            #     counts_global = local_list_shapes
            #     displs_global = [sum(counts_global[:i]) for i in range(size_world)]
            #     gathered_ensemble_global = np.empty(total_size_global, dtype=np.float64)
            # else:
            #     gathered_ensemble_global = None
            #     counts_global = None
            #     displs_global = None
            # comm_world.Gatherv(ens_list, [gathered_ensemble_global, counts_global, displs_global, MPI.DOUBLE], root=0)

            #  free up memory
    
        # subcomm.Barrier()
        del ens_list; del gathered_ensemble; gc.collect()
        if rank_world == 0:
            ensemble_vec = [arr for sublist in gathered_ensemble_global for arr in sublist if arr is not None]
            ensemble_vec = np.column_stack(ensemble_vec) 

            # cap ensemble_vec to (nd,nens) dimensions
            ensemble_vec = ensemble_vec[:, :Nens]
            
            # get the shape of the ensemble
            shape_ens = np.array(ensemble_vec.shape, dtype=np.int32)
        else:
            shape_ens = np.empty(2, dtype=np.int32)

        # broadcast the shape to all processors
        shape_ens = comm_world.bcast(shape_ens, root=0)

    # --- case 3: Form Nens sub-communicators and distribute resources among them ---
    elif Nens < size_world:
        # Ensure all ranks in subcomm are in sync 
        subcomm.Barrier()
        ens = color # each subcomm has a unique color
        model_kwargs.update({'ens_id': ens, 'comm': subcomm})

        # ---- read from file ----
        input_file = f"{_modelrun_datasets}/icesee_ensemble_data.h5"
        with h5py.File(input_file, "r", driver="mpio", comm=subcomm) as f:
            ensemble_vec = f["ensemble"][:,ens,k]
        # ---- end of read from file ----

        # Call the forecast step fucntion- Each subcomm runs the function indepadently
        updated_state = model_module.forecast_step_single(ensemble=ensemble_vec, **model_kwargs)

        # ensemble_vec = gather_and_broadcast_data_default_run(updated_state, subcomm, sub_rank, comm_world, rank_world, params)
        # ensemble_vec = BM.bcast(ensemble_vec, comm_world)
        subcomm.Barrier() #*---
        global_data = {key: subcomm.gather(data, root=0) for key, data in updated_state.items()}

        # Step 2: Process on sub_rank 0
        key_list = list(global_data.keys())
        state_keys = key_list[:params["num_state_vars"]] # Get the state variables to add noise
        if sub_rank == 0:
            # for key in global_data:
            for key in key_list:
                global_data[key] = np.hstack(global_data[key])
                if model_kwargs["joint_estimation"] or params["localization_flag"]:
                    hdim = global_data[key].shape[0] // params["total_state_param_vars"]
                else:
                    hdim = global_data[key].shape[0] // params["num_state_vars"]
                state_block_size = hdim * params["num_state_vars"]  # Compute the state block size
                # Add process noise to the ensembles variables only
                # if key in state_keys:
                #     Q_err = Q_err[:state_block_size, :state_block_size]
                #     q0 = multivariate_normal.rvs(np.zeros(state_block_size), Q_err)
                #     # q0 = np.sqrt(model_kwargs.get("dt",params["dt"]))*multivariate_normal.rvs(np.zeros(state_block_size), Q_err)
                #     global_data[key][:state_block_size] = global_data[key][:state_block_size] + q0[:state_block_size]

                # use pseudorandom fields 
                _time_forecast_noise_generation = MPI.Wtime()
                # if k == 0:
                #     N_size = params["total_state_param_vars"] * hdim
                #     noise = generate_enkf_field(ens,np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, rng=rng, verbose=False)

                noise_all = []
                q0 = []
                for ii, sig in enumerate(params["sig_Q"]):
                    if ii <=params["num_state_vars"]:
                        # W = np.random.normal(0, 1, hdim)
                        # W = generate_pseudo_random_field_1d(hdim,np.sqrt(Lx*Ly), len_scale, verbose=0)
                        model_kwargs.update({"ii_sig": ii, "hdim":hdim, "num_vars":params["total_state_param_vars"]})
                        # W = generate_enkf_field(**model_kwargs)
                        W = generate_enkf_field(ii,np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)
                        noise_ = alpha*noise[ii*hdim:(ii+1)*hdim] + np.sqrt(1 - alpha**2)*W
                        q0.append(noise_)

                        Z = np.sqrt(dt)*sig*rho*noise_
                        noise_all.append(Z)
                noise_ = np.concatenate(noise_all, axis=0)
                global_data[key][:state_block_size] = global_data[key][:state_block_size] + noise_[:state_block_size]
                noise = np.concatenate(q0, axis=0)
                
                del noise_all, q0, noise_, W
                time_forecast_noise_generation += MPI.Wtime() - _time_forecast_noise_generation
                
            # Stack all variables into a single array
            stacked = np.hstack([global_data[key] for key in updated_state.keys()])
            shape_ = np.array(stacked.shape, dtype=np.int32)
            
            # *- compute the mean on each color


            # *- Each color writes each ensemble to the h5 file
            # with h5py.File(input_file, "a", driver="mpio", comm=subcomm) as f:
            #     dset = f['ensemble']
            #     dset[:,ens:ens+1,k+1] = stacked

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

            # add some noise to the ensemble
            # if model_kwargs["joint_estimation"] or params["localization_flag"]:
            #     hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
            # else:
            #     hdim = ensemble_vec.shape[0] // params["num_state_vars"]
            # state_block_size = hdim * params["num_state_vars"]  # Compute the state block size
            # Q_err = Q_err[:state_block_size, :state_block_size]
            # q0 = multivariate_normal.rvs(np.zeros(state_block_size), Q_err)
            # ensemble_vec[:state_block_size, :] = ensemble_vec[:state_block_size, :] + q0[:state_block_size,np.newaxis]

            # hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
            shape_ens = np.array(ensemble_vec.shape, dtype=np.int32)
        else:
            shape_ens = np.empty(2, dtype=np.int32)
            ensemble_vec = np.empty((shape_[0], params["Nens"]), dtype=np.float64)

        # boradcast shape to all processors
        shape_ens = comm_world.bcast(shape_ens, root=0)

        # broadcast the ensemble to all processors
        # ensemble_vec = comm_world.bcast(ensemble_vec, root=0)

    # --- compute the mean
    _time_forecast_ensemble_mean_generation = MPI.Wtime()
    ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], Nens, comm_world, root=0)
    time_forecast_ensemble_mean_generation += MPI.Wtime() - _time_forecast_ensemble_mean_generation

    # update model_kwargs with timing variables and other parameters
    model_kwargs.update({
        "time_forecast_ensemble_generation": time_forecast_ensemble_generation,
        "time_forecast_noise_generation": time_forecast_noise_generation,
        "time_forecast_ensemble_mean_generation": time_forecast_ensemble_mean_generation,
        "time_forecast_file_writing": time_forecast_file_writing,
        "shape_ens": shape_ens,
        "noise": noise,
    })

    return model_kwargs, ensemble_vec, shape_ens, ens_mean



def parallel_forecast_step_default_full_parallel_run(**model_kwargs):
    """
    Process ensemble data in parallel using MPI, with scalable communication via MPI_Gatherv.
    Only rank 0 in each subcommunicator and the global communicator gathers and processes results.
    """
    import numpy as np
    import h5py
    from mpi4py import MPI

    # --- cases 2 & 3 ---
    # case 2: Form batches of sub-communicators and distribute resources among them
    #          - only works for Nens >= size_world
    # case 3: form Nens sub-communicators and distribute resources among them
    #          - only works for size_world > Nens
    #          - even distribution and load balancing leading to performance improvement
    #          - best for size_world/Nens is a whole number

    # Get necessary parameters
    params                    = model_kwargs.get("params")
    Nens                      = params.get("Nens")
    rounds                    = model_kwargs.get("rounds")
    color                     = model_kwargs.get("color", 0)                    # Default color for the subcommunicator
    subcomm_size_min              = model_kwargs.get("subcomm_size_min", 1)            # Size of the subcommunicator
    subcomm                   = model_kwargs.get("subcomm", None)              # Subcommunicator for local ensemble processing
    comm_world                = model_kwargs.get("comm_world", MPI.COMM_WORLD) # Global communicator
    rank_world                = comm_world.Get_rank()                          # Global rank
    sub_rank                  = subcomm.Get_rank() if subcomm else 0          # Rank within the subcommunicator
    _modelrun_datasets        = model_kwargs.get("_modelrun_datasets", "_modelrun_datasets")
    alpha                     = model_kwargs.get("alpha", 0.0)                 # Default alpha value
    rho                       = model_kwargs.get("rho", 1.0)                   # Default rho value
    dt                        = model_kwargs.get("dt", 1.0)                    # Default time step
    Lx                        = model_kwargs.get("Lx", 1.0)                    # Default domain length in x-direction
    Ly                        = model_kwargs.get("Ly", 1.0)                    # Default domain length in y-direction
    len_scale                 = model_kwargs.get("len_scale", 1.0)             # Default length scale for noise generation
    model_module              = model_kwargs.get("model_module", None)         # Module containing the forecast step function
    k                         = model_kwargs.get("k", 0)                      # Time step index, default to 0 if not provided
    noise                     = model_kwargs.get("noise", None)                   # Noise vector, if provided
    rng                       = model_kwargs.get("rng", np.random.default_rng())  # Random number generator, default to numpy's default RNG
    rank_seed = model_kwargs.get("rank_seed", 0)
    enkf_parallel_io = model_kwargs.get("enkf_parallel_io", None)

    # np.random.seed(rank_seed)

    size_world = comm_world.Get_size()  # Total number of ranks in the global communicator

    # get timing variables from model_kwargs 
    time_forecast_ensemble_generation = model_kwargs.get("time_forecast_ensemble_generation", 0.0)
    time_forecast_noise_generation = model_kwargs.get("time_forecast_noise_generation", 0.0)
    time_forecast_ensemble_mean_generation = model_kwargs.get("time_forecast_ensemble_mean_generation", 0.0)
    time_forecast_file_writing = model_kwargs.get("time_forecast_file_writing", 0.0)
                         
    # --- case 2: Form batches of sub-communicators and distribute resources among them ---
    if Nens >= size_world:
        # ensemble_vec, shape_ens = parallel_forecast_step_default_run(**model_kwargs)
        nd = model_kwargs.get("nd", params.get("nd"))
        nt = model_kwargs.get("nt", params.get("nt"))

        # store results for each round      
        ens_list = []
        for round_id in range(rounds):
            ensemble_id = color + round_id * subcomm_size_min  # Global ensemble index
            model_kwargs.update({'ens_id': ensemble_id, 'comm': subcomm})

            if ensemble_id < Nens:  # Only process valid ensembles
                # print(f"[ICESEE] Rank {rank_world} processing ensemble {ensemble_id} in round {round_id + 1}/{rounds}")

                # Ensure all ranks in the subcommunicator are synchronized before running
                subcomm.Barrier()
                ens = ensemble_id
                _local_time_forecast_file_writing_0 = MPI.Wtime()
                # ---- read from file ----
                ensemble_vec = enkf_parallel_io.read_forecast(k, ens)
                # ---- end of read from file ----
                time_forecast_file_writing_0 = MPI.Wtime() - _local_time_forecast_file_writing_0

                # Call the forecast step function
                # hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
                # print(f"[ICESEE] Rank: {rank_world}, min ensemble: {np.min(ensemble_vec[:hdim])}, max ensemble: {np.max(ensemble_vec[:hdim])}")

                updated_state = model_module.forecast_step_single(ensemble=ensemble_vec,**model_kwargs)

                #fetch the updated state
                vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
                for key,value in updated_state.items():
                    ensemble_vec[indx_map[key]] = value

                #  add time evolution noise to the ensemble
                # if k == 0:
                #     # model noise, q0
                #     q0 = np.random.multivariate_normal(np.zeros(nd), Q_err)

                # # squence of white noise drawn from  a smooth pseudorandm fields,w0,
                # # w0 = np.random.normal(0, 1, nd) #TODO: look into this
                # # w0 = np.random.multivariate_normal(np.zeros(nd), np.eye(nd))
                # # q0 = alpha * q0 + np.sqrt(1 - alpha**2) * w0
                # q0 = np.random.multivariate_normal(np.zeros(nd), Q_err)
                # Q_err = Q_err[:state_block_size,:state_block_size]
                # q0 = multivariate_normal.rvs(np.zeros(state_block_size), Q_err)
                # q0 = np.sqrt(model_kwargs.get("dt",params["dt"]))*multivariate_normal.rvs(np.zeros(state_block_size), Q_err)

                # if k+1 <= max(model_kwargs["obs_index"]):
                #     ensemble_vec[:state_block_size] = ensemble_vec[:state_block_size] + q0[:state_block_size]
                # else:
                    #  create a guassian noise with zero mean and variance = 1
                    # q0 = np.random.normal(0, 1, state_block_size)
                    # ensemble_vec[:state_block_size] = ensemble_vec[:state_block_size] + q0[:state_block_size]
                #---------------------------------------------------------------
                if model_kwargs["joint_estimation"] or params["localization_flag"]:
                    hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
                else:
                    hdim = ensemble_vec.shape[0] // params["num_state_vars"]
                state_block_size = hdim * params["num_state_vars"]

                # --- time forecast noise generation ---
                _time_forecast_noise_generation = MPI.Wtime()
                # if k == 0:
                #     # noise = compute_noise_random_fields(ens, hdim, pos, gs_model, params["total_state_param_vars"], L_C)
                #     N_size = params["total_state_param_vars"] * hdim
                #     # noise = generate_pseudo_random_field_1d(N_size,np.sqrt(Lx*Ly), len_scale, verbose=0)
                #     noise = generate_enkf_field(None,np.sqrt(Lx*Ly),hdim, params["total_state_param_vars"], rh=len_scale, rng=rng, verbose=False)

                # noise = noise / np.max(np.abs(noise))
                # if k+1 <= max(model_kwargs["obs_index"]):
                    # W = np.random.normal(0, 1, state_block_size)

                # ======
                noise_all = []
                q0 = []
                for ii, sig in enumerate(params["sig_Q"]):
                    if ii <=params["num_state_vars"]:
                        # W = np.random.normal(0, 1, hdim)
                        # W = generate_pseudo_random_field_1d(hdim,np.sqrt(Lx*Ly), len_scale, verbose=0)
                        model_kwargs.update({"ii_sig": ii, "hdim":hdim, "num_vars":params["total_state_param_vars"]})
                        # W = generate_enkf_field(**model_kwargs)
                        W = generate_enkf_field(ii,np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)
                        noise_ = alpha*noise[ii*hdim:(ii+1)*hdim] + np.sqrt(1 - alpha**2)*W
                        q0.append(noise_)

                        Z = np.sqrt(dt)*sig*rho*noise_
                        noise_all.append(Z)
                noise_ = np.concatenate(noise_all, axis=0)
                ensemble_vec[:state_block_size] = ensemble_vec[:state_block_size] + noise_[:state_block_size]
                noise = np.concatenate(q0, axis=0)
                model_kwargs.update({"noise": noise})  # save the noise to the model_kwargs dictionary
                
                # clean up memory
                del noise_all, q0, noise_, W
                time_forecast_noise_generation += MPI.Wtime() - _time_forecast_noise_generation

                #  time forecast file writing
                _time_forecast_file_writing = MPI.Wtime()
                # enkf_parallel_io.write_forecast(k + 1 if k < nt - 1 else k, ensemble_vec, ens)
                # ensemble_vec_block = 
                enkf_parallel_io.write_forecast(k + 1 if k < nt - 1 else k, ensemble_vec, ens)
                time_forecast_file_writing += MPI.Wtime() - _time_forecast_file_writing + time_forecast_file_writing_0

                shape_ens = np.array(ensemble_vec.shape, dtype=np.int32)

    # --- case 3: Form Nens sub-communicators and distribute resources among them ---
    elif Nens < size_world:
        # Ensure all ranks in subcomm are in sync 
        subcomm.Barrier()
        ens = color # each subcomm has a unique color
        model_kwargs.update({'ens_id': ens, 'comm': subcomm})

        # ---- read from file ----
        input_file = f"{_modelrun_datasets}/icesee_ensemble_data.h5"
        with h5py.File(input_file, "r", driver="mpio", comm=subcomm) as f:
            ensemble_vec = f["ensemble"][:,ens,k]
        # ---- end of read from file ----

        # Call the forecast step fucntion- Each subcomm runs the function indepadently
        updated_state = model_module.forecast_step_single(ensemble=ensemble_vec, **model_kwargs)

        # ensemble_vec = gather_and_broadcast_data_default_run(updated_state, subcomm, sub_rank, comm_world, rank_world, params)
        # ensemble_vec = BM.bcast(ensemble_vec, comm_world)
        subcomm.Barrier() #*---
        global_data = {key: subcomm.gather(data, root=0) for key, data in updated_state.items()}

        # Step 2: Process on sub_rank 0
        key_list = list(global_data.keys())
        state_keys = key_list[:params["num_state_vars"]] # Get the state variables to add noise
        if sub_rank == 0:
            # for key in global_data:
            for key in key_list:
                global_data[key] = np.hstack(global_data[key])
                if model_kwargs["joint_estimation"] or params["localization_flag"]:
                    hdim = global_data[key].shape[0] // params["total_state_param_vars"]
                else:
                    hdim = global_data[key].shape[0] // params["num_state_vars"]
                state_block_size = hdim * params["num_state_vars"]  # Compute the state block size
                # Add process noise to the ensembles variables only
                # if key in state_keys:
                #     Q_err = Q_err[:state_block_size, :state_block_size]
                #     q0 = multivariate_normal.rvs(np.zeros(state_block_size), Q_err)
                #     # q0 = np.sqrt(model_kwargs.get("dt",params["dt"]))*multivariate_normal.rvs(np.zeros(state_block_size), Q_err)
                #     global_data[key][:state_block_size] = global_data[key][:state_block_size] + q0[:state_block_size]

                # use pseudorandom fields 
                _time_forecast_noise_generation = MPI.Wtime()
                # if k == 0:
                #     N_size = params["total_state_param_vars"] * hdim
                #     noise = generate_enkf_field(ens,np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, rng=rng, verbose=False)

                noise_all = []
                q0 = []
                for ii, sig in enumerate(params["sig_Q"]):
                    if ii <=params["num_state_vars"]:
                        # W = np.random.normal(0, 1, hdim)
                        # W = generate_pseudo_random_field_1d(hdim,np.sqrt(Lx*Ly), len_scale, verbose=0)
                        model_kwargs.update({"ii_sig": ii, "hdim":hdim, "num_vars":params["total_state_param_vars"]})
                        # W = generate_enkf_field(**model_kwargs)
                        W = generate_enkf_field(ii,np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)
                        noise_ = alpha*noise[ii*hdim:(ii+1)*hdim] + np.sqrt(1 - alpha**2)*W
                        q0.append(noise_)

                        Z = np.sqrt(dt)*sig*rho*noise_
                        noise_all.append(Z)
                noise_ = np.concatenate(noise_all, axis=0)
                global_data[key][:state_block_size] = global_data[key][:state_block_size] + noise_[:state_block_size]
                noise = np.concatenate(q0, axis=0)
                
                del noise_all, q0, noise_, W
                time_forecast_noise_generation += MPI.Wtime() - _time_forecast_noise_generation
                
            # Stack all variables into a single array
            stacked = np.hstack([global_data[key] for key in updated_state.keys()])
            shape_ = np.array(stacked.shape, dtype=np.int32)
            
            # *- compute the mean on each color


            # *- Each color writes each ensemble to the h5 file
            # with h5py.File(input_file, "a", driver="mpio", comm=subcomm) as f:
            #     dset = f['ensemble']
            #     dset[:,ens:ens+1,k+1] = stacked

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

            # add some noise to the ensemble
            # if model_kwargs["joint_estimation"] or params["localization_flag"]:
            #     hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
            # else:
            #     hdim = ensemble_vec.shape[0] // params["num_state_vars"]
            # state_block_size = hdim * params["num_state_vars"]  # Compute the state block size
            # Q_err = Q_err[:state_block_size, :state_block_size]
            # q0 = multivariate_normal.rvs(np.zeros(state_block_size), Q_err)
            # ensemble_vec[:state_block_size, :] = ensemble_vec[:state_block_size, :] + q0[:state_block_size,np.newaxis]

            # hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
            shape_ens = np.array(ensemble_vec.shape, dtype=np.int32)
        else:
            shape_ens = np.empty(2, dtype=np.int32)
            ensemble_vec = np.empty((shape_[0], params["Nens"]), dtype=np.float64)

        # boradcast shape to all processors
        shape_ens = comm_world.bcast(shape_ens, root=0)

        # broadcast the ensemble to all processors
        # ensemble_vec = comm_world.bcast(ensemble_vec, root=0)

    # --- compute the mean
    _time_forecast_ensemble_mean_generation = MPI.Wtime()
    # enkf_parallel_io.compute_forecast_mean_chunked(k + 1 if k < nt - 1 else k)
     # only compute the mean only when we have to observe (sine we can gnerate the during post-processing)
    km = model_kwargs.get("km", 0)
    k = model_kwargs.get("k", 0)
    tobserve = model_kwargs.get("tobserve")
    m_obs = model_kwargs.get("m_obs", params["number_obs_instants"])
    obs_index = model_kwargs["obs_index"]
    if (km < params["number_obs_instants"]) and (k == obs_index[km]):
        enkf_parallel_io.compute_forecast_mean_chunked_v2(k + 1 if k < nt - 1 else k, flag='initial')
        # enkf_parallel_io.compute_forecast_mean_chunked_v2(k, flag='initial')
    time_forecast_ensemble_mean_generation += MPI.Wtime() - _time_forecast_ensemble_mean_generation

    # update model_kwargs with timing variables and other parameters
    model_kwargs.update({
        "time_forecast_ensemble_generation": time_forecast_ensemble_generation,
        "time_forecast_noise_generation": time_forecast_noise_generation,
        "time_forecast_ensemble_mean_generation": time_forecast_ensemble_mean_generation,
        "time_forecast_file_writing": time_forecast_file_writing,
        "shape_ens": shape_ens,
        "noise": noise,
    })

    return model_kwargs
    # return model_kwargs, ensemble_vec, shape_ens, enkf_parallel_io.ensemble_mean



def parallel_forecast_step_even_distribution_run(**model_kwargs):
    """ Parallel run of the forecast step for each ensemble member.
        This function is designed to be used in a distributed environment where each rank processes a single ensemble member.
        It assumes that the number of ensemble members (Nens) is divisible by the size of the world communicator (size_world).
    """

    # unpack model_kwargs
    params = model_kwargs.get("params")
    model_module = model_kwargs.get("model_module")
    comm_world = model_kwargs.get("comm_world", MPI.COMM_WORLD)
    rank_world = comm_world.Get_rank()
    Nens = params.get("Nens", 1)  # Number of ensemble members
    nd = params.get("nd", 1)  # Dimension of the state vector
    Q_err = model_kwargs.get("Q_err", np.eye(nd))  # Error covariance matrix
    state_block_size = model_kwargs.get("state_block_size", nd)  # Size of the state block
    size_world = comm_world.Get_size()  # Total number of ranks in the world communicato
    ensemble_vec = model_kwargs.get("ensemble_vec", np.zeros((nd, Nens), dtype=np.float64))  # Initialize ensemble vector
    ensemble_vec_mean = model_kwargs.get("ensemble_vec_mean", np.zeros((nd, params.get("nt", params["nt"]) + 1), dtype=np.float64))  # Initialize ensemble mean vector
    shape_ens = np.array(ensemble_vec.shape, dtype=np.int32)  # Shape of the ensemble vector
    ensemble_local = model_kwargs.get("ensemble_local", np.zeros((nd, Nens), dtype=np.float64))  # Local ensemble vector
    k = model_kwargs.get("k", 0)  # Time step index, default to 0 if not provided

    # check if Nens is divisible by size_world and greater or equal to size_world
    if Nens >= size_world and Nens % size_world == 0:
        for ens in range(ensemble_local.shape[1]):
            ensemble_local[:, ens] = model_module.forecast_step_single(ensemble=ensemble_local, **model_kwargs)
            # q0 = np.random.multivariate_normal(np.zeros(nd), Q_err)
            Q_err = Q_err[:state_block_size,:state_block_size]
            q0 = multivariate_normal.rvs(np.zeros(state_block_size), Q_err)
            ensemble_local[:state_block_size,ens] = ensemble_local[:state_block_size,ens] + q0[:state_block_size]

        # --- compute the ensemble mean ---
        ensemble_vec_mean[:,k+1] = ParallelManager().compute_mean_from_local_matrix(ensemble_local, comm_world)

        # --- gather all local ensembles from all processors to root---
        gathered_ensemble = ParallelManager().gather_data(comm_world, ensemble_local, root=0)
        if rank_world == 0:
            ensemble_vec = np.hstack(gathered_ensemble)
        else:
            ensemble_vec = np.empty((nd, Nens), dtype=np.float64)

    return ensemble_vec, ensemble_vec_mean, shape_ens

def parallel_forecast_step_squential_run(**model_kwargs):
    """ Squential run of the forecast step for each ensemble member.
        This function is designed to be used in a distributed environment where each rank processes a single ensemble member.
        #TODO: still under development, not fully tested.
    """

    # unpack model_kwargs
    params = model_kwargs.get("params")
    model_module = model_kwargs.get("model_module")
    comm_world = model_kwargs.get("comm_world", MPI.COMM_WORLD)
    rank_world = comm_world.Get_rank()
    Nens = params.get("Nens", 1)  # Number of ensemble members
    nd = params.get("nd", 1)  # Dimension of the state vector
    Q_err = model_kwargs.get("Q_err", np.eye(nd))  # Error covariance matrix
    state_block_size = model_kwargs.get("state_block_size", nd)  # Size of the state block
    ensemble_vec = model_kwargs.get("ensemble_vec", np.zeros((nd, Nens), dtype=np.float64))  # Initialize ensemble vector
    ensemble_vec_mean = model_kwargs.get("ensemble_vec_mean", np.zeros((nd, params.get("nt", params["nt"]) + 1), dtype=np.float64))  # Initialize ensemble mean vector
    shape_ens = np.array(ensemble_vec.shape, dtype=np.int32)  # Shape of the ensemble vector
    ensemble_local = model_kwargs.get("ensemble_local", np.zeros((nd, Nens), dtype=np.float64))  # Local ensemble vector
    k = model_kwargs.get("k", 0)  # Time step index, default to 0 if not provided

    ensemble_col_stack = []
    for ens in range(Nens):
        comm_world.Barrier() # make sure all processors are in sync
        ensemble_vec[:,ens] = model_module.forecast_step_single(ens=ens, ensemble=ensemble_vec, nd=nd,  **model_kwargs)
        q0 = np.random.multivariate_normal(np.zeros(nd), Q_err)
        ensemble_vec[:state_block_size,ens] = ensemble_vec[:state_block_size,ens] + q0[:state_block_size]
        comm_world.Barrier() # make sure all processors reach this point before moving on
        
        # gather the ensemble from all processors to rank 0
        gathered_ensemble = ParallelManager().gather_data(comm_world, ensemble_vec, root=0)
        if rank_world == 0:
            # print(f"[ICESEE] [Rank {rank_world}] Gathered shapes: {[arr.shape for arr in ens_all]}")
            ensemble_stack = np.hstack(gathered_ensemble)
            # print(f"[ICESEE] Ensemble stack shape: {ensemble_stack.shape}")
            ensemble_col_stack.append(ensemble_stack)
    
    # transpose the ensemble column
    if rank_world == 0:
        ens_T = np.array(ensemble_col_stack).T
        print(f"[ICESEE] Ensemble column shape: {ens_T.shape}")
        shape_ens = np.array(ens_T.shape, dtype=np.int32) # send shape info
    else:
        shape_ens = np.empty(2, dtype=np.int32)
    exit()
    # broadcast the shape to all processors
    comm_world.Bcast([shape_ens, MPI.INT], root=0)

    if rank_world != 0:
        # if k == 0:
        ens_T = np.empty(shape_ens, dtype=np.float64)

    # broadcast the ensemble to all processors
    comm_world.Bcast([ens_T, MPI.DOUBLE], root=0)
    # print(f"[ICESEE] Rank: {rank_world}, Ensemble shape: {ens_T.shape}")

    # compute the ensemble mean
    # if k == 0: # only do this at the first time step
    #     # gather from all processors ensemble_vec_mean[:,k+1]
    #     gathered_ensemble_vec_mean = comm_world.allgather(ensemble_vec_mean[:,k])
    #     if rank_world == 0:
    #         # print(f"[ICESEE] Ensemble mean shape: {[arr.shape for arr in gathered_ensemble_vec_mean]}")
    #         stack_ensemble_vec_mean = np.hstack(gathered_ensemble_vec_mean)
    #         ensemble_vec_mean = np.empty((shape_ens[0],model_kwargs.get("nt",params["nt"])+1), dtype=np.float64)
    #         ensemble_vec_mean[:,k] = np.mean(stack_ensemble_vec_mean, axis=1)
    #     else: 
    #         ensemble_vec_mean = np.empty((shape_ens[0],model_kwargs.get("nt",params["nt"])), dtype=np.float64)
        
    #     # broadcast the ensemble mean to all processors
    #     comm_world.Bcast([ensemble_vec_mean, MPI.DOUBLE], root=0)
    #     print(f"[ICESEE] Rank: {rank_world}, Ensemble mean shape: {ensemble_vec_mean.shape}") 

    ensemble_vec_mean[:,k+1] = np.mean(ens_T[:nd,:], axis=1)
    # ensemble_vec_mean[:,k+1] = ParallelManager().compute_mean(ens_T[:nd,:], comm_world)

    
    return model_kwargs, ensemble_vec, ensemble_vec_mean, shape_ens