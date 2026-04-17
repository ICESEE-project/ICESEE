# ==============================================================================
# @des: This file contains run functions for the ICESEE model to initialize the ensemble.
# @date: 2025-07-30
# @author: Brian Kyanjo
# ==============================================================================

# --- import necessary libraries ---
import numpy as np
import h5py
import gc
import zarr
import os
from mpi4py import MPI

from ICESEE.src.utils.tools import icesee_get_index, env_flag
from ICESEE.src.run_model_da._error_generation import compute_Q_err_random_fields, \
                              compute_noise_random_fields, \
                              generate_pseudo_random_field_1d, \
                              generate_pseudo_random_field_2D, \
                              generate_enkf_field
from ICESEE.src.parallelization._parallel_i_o import parallel_write_full_ensemble_from_root,parallel_write_full_ensemble_from_root_full_parallel_run

from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager
# rank_seed, rng = ParallelManager().initialize_seed(MPI.COMM_WORLD)

def ensemble_initialization(**model_kwargs):
    """Initialize the ensemble for the ICESEE model.
    """
    
    # unpack model_kwargs
    params         = model_kwargs.get("params")
    model_module   = model_kwargs.get("model_module", None)
    comm_world     = model_kwargs.get("comm_world", MPI.COMM_WORLD)
    subcomm        = model_kwargs.get("subcomm", None)
    color          = model_kwargs.get("color", 0)
    pos            = model_kwargs.get("pos", None)
    gs_model       = model_kwargs.get("gs_model", None)
    L_C           = model_kwargs.get("L_C", None)
    Lx             = model_kwargs.get("Lx", 1.0)
    Ly             = model_kwargs.get("Ly", 1.0)
    len_scale      = model_kwargs.get("len_scale", 1.0)
    Q_rho          = model_kwargs.get("Q_rho", 1.0)
    model_nprocs   = model_kwargs.get("model_nprocs", 1)
    total_cores    = model_kwargs.get("total_cores", 1)
    base_total_procs = model_kwargs.get("base_total_procs", 1)
    rounds         = model_kwargs.get("rounds", 1)
    subcomm_size_min   = model_kwargs.get("subcomm_size_min", 1)
    rng           = model_kwargs.get("rng", np.random.default_rng())
    rank_seed = model_kwargs.get("rank_seed", 0)
    alpha = model_kwargs.get("initial_spread_factor")


    sub_rank  = subcomm.Get_rank()
    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()

    time_init_noise_generation = 0.0
    time_init_file_writing     = 0.0
    time_init_ensemble_mean_computation = 0.0

    if params["even_distribution"] or (params["default_run"] and size_world <= params["Nens"]):
        if params["default_run"] and size_world <= params["Nens"] and not (model_kwargs.get("sequential_ensemble_initialization", False)):
        # if False:
            if rank_world == 0:
                print("[ICESEE] Initializing the ensemble ...")

            # model_kwargs.update({'ens_id': rank_world})
            Nens = params["Nens"]
            model_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})

            model_kwargs.update({"statevec_ens":np.zeros([params["nd"], params["Nens"]])})

            # get the ensemble matrix   
            vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
            ensemble_vec = np.zeros_like(model_kwargs["statevec_ens"])

            if model_kwargs["joint_estimation"] or params["localization_flag"]:
                    hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
            else:
                hdim = ensemble_vec.shape[0] // params["num_state_vars"]
            state_block_size = hdim * params["num_state_vars"]
            
            if sub_rank == 0:
                ens_list_init = []
            else:
                ens_list_init = []

            for round_id in range(rounds):
                ensemble_id = color + (round_id * subcomm_size_min)
                model_kwargs.update({'ens_id': ensemble_id})

                if ensemble_id < Nens:
                    # Synchronize the ensemble initialization
                    subcomm.Barrier()
                    ens = ensemble_id

                    # Call the model to initialize the ensemble
                    data = model_module.initialize_ensemble(ens, **model_kwargs)
                    for key, value in data.items():
                        ensemble_vec[indx_map[key], ens] = value

                    # Add process noise in-place to avoid temporary array
                    _time_init_noise_generation = MPI.Wtime()
                    model_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":params["total_state_param_vars"]})
                    noise = generate_enkf_field(**model_kwargs)
                    time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                    ensemble_vec[:, ens] += alpha*noise

                    # for ii, sig in enumerate(params["sig_Q"]):
                    #     if ii <=params["num_state_vars"]:
                    #         start_idx = ii * hdim
                    #         end_idx = start_idx + hdim
                    #         # make sure the noise has zero mean
                    #         noise[start_idx:end_idx] *=sig
                    # #         noise[start_idx:end_idx] -= np.mean(noise[start_idx:end_idx])
                    # #         inflated_noise = noise[start_idx:end_idx] / np.max(np.abs(noise[start_idx:end_idx]))
                    #         ensemble_vec[start_idx:end_idx, ens] += noise[start_idx:end_idx]

                    # for ii, sig in enumerate(params["sig_Q"]):
                    #     # if ii <=params["num_state_vars"]:
                    #     start_idx = ii * hdim
                    #     end_idx = start_idx + hdim
                    #     ensemble_vec[start_idx:end_idx, ens] += noise[start_idx:end_idx] * sig

                    del noise  # Free memory immediately

                    # Gather ensemble data efficiently
                    gathered_ensemble = subcomm.gather(ensemble_vec[:, ens], root=0)
                    
                    if sub_rank == 0:
                        # Use np.concatenate with pre-allocated array to avoid np.hstack memory overhead
                        gathered_ensemble = np.concatenate(gathered_ensemble, axis=0)
                        ens_list_init.append(gathered_ensemble)
                    
                    del gathered_ensemble  # Free memory after appending

                # subcomm.Barrier()
                # Gather all ensembles from all subcommunicators
                gathered_ensemble_global = ParallelManager().gather_data(comm_world, ens_list_init, root=0)

            # Final processing on rank 0
            # subcomm.Barrier()
            del ens_list_init; gc.collect()
            if rank_world == 0:
                # Flatten and filter None values
                ensemble_vec = [arr for sublist in gathered_ensemble_global for arr in sublist if arr is not None]
                # Pre-allocate final array to avoid memory spikes during np.column_stack
                final_shape = (len(ensemble_vec[0]), len(ensemble_vec)) if ensemble_vec else (0, 0)
                ensemble_vec_final = np.empty(final_shape, dtype=ensemble_vec[0].dtype)
                for i, arr in enumerate(ensemble_vec):
                    ensemble_vec_final[:, i] = arr
                shape_ens = np.array(ensemble_vec_final.shape, dtype=np.int32)
                ensemble_vec = ensemble_vec_final  # Replace ensemble_vec with final array
                del ensemble_vec_final  # Free memory
            else:
                shape_ens = np.empty(2, dtype=np.int32)

            # Broadcast the shape of the ensemble
            shape_ens = comm_world.bcast(shape_ens, root=0)

        else:
            if rank_world == 0:
                print("[ICESEE] Initializing the ensemble ...")
                model_kwargs.update({'ens_id': rank_world})
                if params["even_distribution"]:
                    model_kwargs.update({'rank': rank_world, 'color': color, 'comm': comm_world})
                else:
                    model_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})

                model_kwargs.update({"statevec_ens":np.zeros([params["nd"], params["Nens"]])})
                
                # get the ensemble matrix   
                vecs, indx_map, dim_per_proc = icesee_get_index(model_kwargs["statevec_ens"], **model_kwargs)
                ensemble_vec = np.zeros_like(model_kwargs["statevec_ens"])

                if model_kwargs["joint_estimation"] or params["localization_flag"]:
                        hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
                else:
                    hdim = ensemble_vec.shape[0] // params["num_state_vars"]
                state_block_size = hdim * params["num_state_vars"]

                # # --- get the process noise ---
                # pos, gs_model, L_C = compute_Q_err_random_fields(hdim, params["total_state_param_vars"], params["sig_Q"], Q_rho, len_scale)

                # process_noise = []
                for ens in range(params["Nens"]):
                    # model_kwargs.update({"ens_id": ens})
                    data = model_module.initialize_ensemble(ens,**model_kwargs)
                
                    # iterate over the data and update the ensemble
                    for key, value in data.items():
                        ensemble_vec[indx_map[key],ens] = value

                    # --->
                    # noise = compute_noise_random_fields(ens, hdim, pos, gs_model, params["total_state_param_vars"], L_C)
                    # ensemble_vec[:,ens] += noise
                    #----->
                    _time_init_noise_generation = MPI.Wtime()
                    N_size = params["total_state_param_vars"] * hdim
                    # noise = generate_pseudo_random_field_1d(N_size,np.sqrt(Lx*Ly), len_scale, verbose=True)
                    model_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":params["total_state_param_vars"]})
                    noise = generate_enkf_field(**model_kwargs)
                    time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation

                    # lets inflate the noise to increase the spread
                    model_kwargs['observed_vars_params'] = (model_kwargs['observed_vars'] + model_kwargs['observed_params'])
                    # exclude bed variables from observed variables
                    all_observed = model_kwargs['observed_vars_params']
                    model_kwargs['all_observed'] = all_observed; params['all_observed'] = all_observed
                    nd_new = len(all_observed) * hdim
                    model_kwargs['nd_observed'] = nd_new
                    # for ii, key in enumerate(all_observed):
                    #     # if ii < params["num_state_vars"]:
                    #     ensemble_vec[indx_map[key], ens] += alpha * noise[indx_map[key]]
                        # if ii < params["num_state_vars"]:
                        #     start_idx = ii * hdim
                        #     end_idx = start_idx + hdim
                        #     ensemble_vec[start_idx:end_idx, ens] += alpha * noise[start_idx:end_idx]
                            
                    ensemble_vec[:,ens] += alpha*noise
                    # for ii, sig in enumerate(params["sig_Q"]):
                    #     if ii <=params["num_state_vars"]:
                    #         start_idx = ii * hdim
                    #         end_idx = start_idx + hdim
                    #         ensemble_vec[start_idx:end_idx, ens] += noise[start_idx:end_idx] * sig
                    
                shape_ens = np.array(ensemble_vec.shape,dtype=np.int32)
                
    
            else:
                ensemble_vec = np.empty((params["nd"],params["Nens"]),dtype=np.float64)
                shape_ens = np.empty(2,dtype=np.int32)
                # pos, gs_model, L_C

        comm_world.Barrier()

        # now reset the model_nprocs
        if rank_world == 0:
            diff = total_cores - base_total_procs 
            if diff >= 0:
                # split the diff amaongest all processors
                min_model_nprocs = max(model_nprocs-1, 1) 
                if model_kwargs.get('ICESEE_PERFORMANCE_TEST') or env_flag("ICESEE_PERFORMANCE_TEST", default=False):
                    model_nprocs = model_nprocs
                else:
                    model_nprocs = max(min_model_nprocs, model_nprocs + (diff // size_world))
            else:
                model_nprocs = model_nprocs

        model_nprocs = comm_world.bcast(model_nprocs, root=0)
        model_kwargs.update({'model_nprocs': model_nprocs})

        if params["even_distribution"]:
            # Bcast the ensemble
            comm_world.Bcast(ensemble_vec, root=0)
            ensemble_bg = np.empty((params["nd"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)
            ensemble_vec_mean = np.empty((params["nd"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)
            ensemble_vec_full = np.empty((params["nd"],params["Nens"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)
            ensemble_vec_mean[:,0] = np.mean(ensemble_vec, axis=1)
            ensemble_vec_full[:,:,0] = ensemble_vec
            ensemble_bg[:,0] = ensemble_vec_mean[:,0]
        else:
            # broadcast the shape of the ensemble
            shape_ens = comm_world.bcast(shape_ens, root=0)
            # write the ensemble to the file

            # -- time ensemble mean computation ---
            _time_init_ensemble_mean_computation = MPI.Wtime()
            ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], params['Nens'], comm_world, root=0)
            time_init_ensemble_mean_computation += MPI.Wtime() - _time_init_ensemble_mean_computation

            # ---time file writing ---
            _time_init_file_writing = MPI.Wtime()
            parallel_write_full_ensemble_from_root(0, ens_mean, model_kwargs,ensemble_vec,comm_world)
            time_init_file_writing += MPI.Wtime() - _time_init_file_writing

        # comm_world.Bcast(ensemble_vec, root=0)
        # hdim = params["nd"] // params["total_state_param_vars"]
        # print(f"[ICESEE] Rank: {rank_world}, min ensemble: {np.min(ensemble_vec[hdim,:])}, max ensemble: {np.max(ensemble_vec[hdim,:])}")
        # exit()
    else:
        if rank_world == 0:
            print("[ICESEE] Initializing the ensemble ...")
        
        if params["default_run"] and size_world > params["Nens"]:
            # debug
            sub_shape = model_kwargs['dim_list'][sub_rank]
            model_kwargs.update({"statevec_ens":np.zeros((sub_shape, params["Nens"]))})

            model_kwargs.update({"ens_id": color, "rank": sub_rank, "color": color, "comm": subcomm})

            # ensemble_vec, shape_ens  = model_module.initialize_ensemble_debug(color,**model_kwargs)
            # ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], params['Nens'], comm_world, root=0)
            # parallel_write_full_ensemble_from_root(0, ens_mean, params,ensemble_vec,comm_world)
            # -----------------------------------------------------

            ens = color
            # model_kwargs.update({"statevec_ens":np.zeros((model_kwargs['global_shape'], params["Nens"]))})
            initialilaized_state = model_module.initialize_ensemble(ens,**model_kwargs)
            # ensemble_vec, shape_ens = gather_and_broadcast_data_default_run(initialilaized_state, subcomm, sub_rank, comm_world, rank_world, params)
            # ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], params['Nens'], comm_world, root=0)
            # parallel_write_full_ensemble_from_root(0, ens_mean, params,ensemble_vec,comm_world)
            # ensemble_vec = BM.bcast(ensemble_vec, comm_world)
            
            initial_data = {key: subcomm.gather(value, root=0) for key, value in initialilaized_state.items()}
            key_list = list(initial_data.keys())
            state_keys = key_list[:params["num_state_vars"]]
            if sub_rank == 0:
                # for key in initial_data:
                for key in key_list:
                    initial_data[key] = np.hstack(initial_data[key])
                    if model_kwargs["joint_estimation"] or params["localization_flag"]:
                        hdim = initial_data[key].shape[0] // params["total_state_param_vars"]
                    else:
                        hdim = initial_data[key].shape[0] // params["num_state_vars"]
                    state_block_size = hdim*params["num_state_vars"]
                    full_block_size = hdim*params["total_state_param_vars"]
                    # if key in state_keys:
                        # noise = np.random.normal(0, 0.1, state_block_size)
                        # Q_err = np.eye(state_block_size) * params["sig_Q"] ** 2
                        # Q_err = np.eye(state_block_size) * 0.01 ** 2
                    if model_kwargs.get("random_fields",False):
                        Q_err = np.zeros((full_block_size,full_block_size))
                        for i, sig in enumerate(params["sig_Q"]):
                            start_idx = i *hdim
                            end_idx = start_idx + hdim
                            Q_err[start_idx:end_idx,start_idx:end_idx] = np.eye(hdim) * sig ** 2

                        # noise = multivariate_normal.rvs(mean=np.zeros(state_block_size), cov=Q_err)
                        _time_init_noise_generation = MPI.Wtime()
                        noise = compute_noise_random_fields(ens, hdim, pos, gs_model, params["total_state_param_vars"], L_C)
                        time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                        # initial_data[key][:state_block_size] += noise[:state_block_size]
                        # noise = noise / np.max(np.abs(noise))
                        initial_data[key] += noise
                    else:
                        N_size = params["total_state_param_vars"] * hdim
                        _time_init_noise_generation = MPI.Wtime()
                        model_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":params["total_state_param_vars"]})
                        noise = generate_enkf_field(**model_kwargs)
                        time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                        initial_data[key] += noise
                        # for ii, sig in enumerate(params["sig_Q"]):
                        #     start_idx = ii *hdim
                        #     end_idx = start_idx + hdim
                        #     initial_data[key][start_idx:end_idx] += noise[start_idx:end_idx]*sig
                    
                # stack all variables together into a single array
                stacked = np.hstack([initial_data[key] for key in initialilaized_state.keys()])
                shape_ens = np.array(stacked.shape,dtype=np.int32)
            else:
                shape_ens = np.empty(2,dtype=np.int32)

            # broadcast the shape of the initialized ensemble
            shape_ens = comm_world.bcast(shape_ens, root=0)

            if sub_rank != 0:
                stacked = np.empty(shape_ens,dtype=np.float64)

            all_init = comm_world.gather(stacked if sub_rank == 0 else None, root=0)

            if rank_world == 0:
                all_init = [arr for arr in all_init if isinstance(arr, np.ndarray)]
                ensemble_vec = np.column_stack(all_init)
                # print(f"[ICESEE] Shape of the ensemble: {ensemble_vec.shape}")
            else:
                ensemble_vec = np.empty((model_kwargs["global_shape"],params["Nens"]),dtype=np.float64)
            
            time_init_ensemble_mean_computation = MPI.Wtime()
            ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], params['Nens'], comm_world, root=0)
            time_init_ensemble_mean_computation = MPI.Wtime() - time_init_ensemble_mean_computation

            _time_init_file_writing = MPI.Wtime()
            parallel_write_full_ensemble_from_root(0, ens_mean, model_kwargs,ensemble_vec,comm_world)
            time_init_file_writing += MPI.Wtime() - _time_init_file_writing
            
        elif params["sequential_run"]:
            comm_world.Barrier()
            sub_shape = model_kwargs['dim_list'][rank_world]
            model_kwargs.update({"statevec_ens":np.zeros([model_kwargs["global_shape"], params["Nens"]]),
                                "statevec_ens_mean":np.zeros([model_kwargs["global_shape"], model_kwargs.get("nt",params["nt"]) + 1]),
                                "statevec_ens_full":np.zeros([model_kwargs["global_shape"], params["Nens"], model_kwargs.get("nt",params["nt"]) + 1]),
                                "statevec_bg":np.zeros([model_kwargs["global_shape"], model_kwargs.get("nt",params["nt"]) + 1])})
            ensemble_bg, ensemble_vec, ensemble_vec_mean, ensemble_vec_full = model_module.initialize_ensemble(**model_kwargs)

            # gather from every rank to rank 0
            gathered_ensemble = comm_world.gather(ensemble_vec[:sub_shape,:], root=0)
            if rank_world == 0:
                ensemble_vec = np.vstack(gathered_ensemble)
                print(f"[ICESEE] Shape of the ensemble: {ensemble_vec.shape}")
                ensemble_vec_mean[:,0] = np.mean(ensemble_vec, axis=1)
                ensemble_vec_full[:,:,0] = ensemble_vec
            else:
                ensemble_vec = np.empty((model_kwargs["global_shape"],params["Nens"]),dtype=np.float64)
                ensemble_vec_mean = np.empty((model_kwargs["global_shape"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)
                ensemble_vec_full = np.empty((model_kwargs["global_shape"],params["Nens"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)

            # else:
            #     ensemble_bg = np.empty((model_kwargs["global_shape"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)
            #     ensemble_vec = np.empty((model_kwargs["global_shape"],params["Nens"]),dtype=np.float64)
            #     ensemble_vec_mean = np.empty((model_kwargs["global_shape"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)
            #     ensemble_vec_full = np.empty((model_kwargs["global_shape"],params["Nens"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)

            # # Bcast the ensemble
            # comm_world.Bcast(ensemble_bg, root=0)
            comm_world.Bcast(ensemble_vec, root=0)
            comm_world.Bcast(ensemble_vec_mean, root=0)
            comm_world.Bcast(ensemble_vec_full, root=0)

            # hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
            # print(f"[ICESEE] rank: {rank_world}, subrank: {sub_rank}, min ensemble: {np.min(ensemble_vec[hdim,:])}, max ensemble: {np.max(ensemble_vec[hdim,:])}")

    if params.get("default_run", False):
        return model_kwargs, ensemble_vec, time_init_noise_generation, \
               time_init_ensemble_mean_computation, time_init_file_writing, \
                shape_ens, None, None, None
    else:
        return model_kwargs, ensemble_vec, time_init_noise_generation, \
               time_init_ensemble_mean_computation, time_init_file_writing, \
                shape_ens,ensemble_bg,  ensemble_vec_mean, ensemble_vec_full
    


def ensemble_initialization_full_parallel_run(**model_kwargs):
    """Initialize the ensemble for the ICESEE model.
    """
    
    # unpack model_kwargs
    params         = model_kwargs.get("params")
    model_module   = model_kwargs.get("model_module", None)
    comm_world     = model_kwargs.get("comm_world", MPI.COMM_WORLD)
    subcomm        = model_kwargs.get("subcomm", None)
    color          = model_kwargs.get("color", 0)
    pos            = model_kwargs.get("pos", None)
    gs_model       = model_kwargs.get("gs_model", None)
    L_C           = model_kwargs.get("L_C", None)
    Lx             = model_kwargs.get("Lx", 1.0)
    Ly             = model_kwargs.get("Ly", 1.0)
    len_scale      = model_kwargs.get("len_scale", 1.0)
    Q_rho          = model_kwargs.get("Q_rho", 1.0)
    model_nprocs   = model_kwargs.get("model_nprocs", 1)
    total_cores    = model_kwargs.get("total_cores", 1)
    base_total_procs = model_kwargs.get("base_total_procs", 1)
    rounds         = model_kwargs.get("rounds", 1)
    subcomm_size_min   = model_kwargs.get("subcomm_size_min", 1)
    rng           = model_kwargs.get("rng", np.random.default_rng())
    rank_seed = model_kwargs.get("rank_seed", 0)
    data_path = model_kwargs.get("data_path", "_modeldatasets")
    enkf_parallel_io = model_kwargs.get("enkf_parallel_io", None)
    alpha       = model_kwargs.get("initial_spread_factor")

    sub_rank     = subcomm.Get_rank()
    rank_world   = comm_world.Get_rank()
    size_world   = comm_world.Get_size()

    time_init_noise_generation = 0.0
    time_init_file_writing     = 0.0
    time_init_ensemble_mean_computation = 0.0

    if params["even_distribution"] or (params["default_run"] and size_world <= params["Nens"]):
        if params["default_run"] and size_world <= params["Nens"] and not (model_kwargs.get("sequential_ensemble_initialization", False)):
        # if False:
            if rank_world == 0:
                print("[ICESEE] Initializing the ensemble ...")

            # model_kwargs.update({'ens_id': rank_world})
            Nens = params["Nens"]
            nd = model_kwargs.get("nd", params["nd"])
            model_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})

            model_kwargs.update({"statevec_ens":np.zeros([params["nd"], params["Nens"]])})

            # get the ensemble matrix   
            vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
            # ensemble_vec = np.zeros_like(model_kwargs["statevec_ens"])
            # store=f"{data_path}/ensemble_initialization_{color}.zarr"
            # chunk_size = (min(nd, 1000), 1)
            # ensemble_vec = zarr.create_array(store=store, shape=(params["nd"], params["Nens"]), chunks=chunk_size, dtype=np.float64, overwrite=True)
            ensemble_vec = np.zeros(nd, dtype=np.float64)

            if model_kwargs["joint_estimation"] or params["localization_flag"]:
                    hdim = nd // params["total_state_param_vars"]
            else:
                hdim = nd // params["num_state_vars"]
            state_block_size = hdim * params["num_state_vars"]

            for round_id in range(rounds):
                ensemble_id = color + (round_id * subcomm_size_min)
                model_kwargs.update({'ens_id': ensemble_id})

                if ensemble_id < Nens:
                    # Synchronize the ensemble initialization
                    # subcomm.Barrier()
                    # comm_world.Barrier()
                    ens = ensemble_id

                    # Call the model to initialize the ensemble
                    data = model_module.initialize_ensemble(ens, **model_kwargs)
                    for key, value in data.items():
                        # ensemble_vec[indx_map[key], ens] = value
                        ensemble_vec[indx_map[key]] = value

                    # Add process noise in-place to avoid temporary array
                    _time_init_noise_generation = MPI.Wtime()
                    model_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":params["total_state_param_vars"]})
                    noise = generate_enkf_field(**model_kwargs)
                    time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                    # ensemble_vec[:,ens] += noise
                    ensemble_vec += alpha*noise
                      
                    _time_init_file_writing = MPI.Wtime()    
                    enkf_parallel_io.write_forecast(0, ensemble_vec, ensemble_id)
                    # enkf_parallel_io.datasets[0][:, ens] = ensemble_vec
                    time_init_file_writing += MPI.Wtime() - _time_init_file_writing
         
        else:
            if rank_world == 0:
                print("[ICESEE] Initializing the ensemble ...")
                model_kwargs.update({'ens_id': rank_world})
                if params["even_distribution"]:
                    model_kwargs.update({'rank': rank_world, 'color': color, 'comm': comm_world})
                else:
                    model_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})

                nd = model_kwargs.get("nd", params["nd"])
                
                # get the ensemble matrix   
                vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
                # ensemble_vec = np.zeros_like(model_kwargs["statevec_ens"])
                store=f"{data_path}/ensemble_initialization.zarr"
                chunk_size = (min(nd, 1000), 1)
                ensemble_vec = zarr.create_array(store=store, shape=(params["nd"], params["Nens"]), chunks=chunk_size, dtype=np.float64, overwrite=True)

                if model_kwargs["joint_estimation"] or params["localization_flag"]:
                        hdim = nd // params["total_state_param_vars"]
                else:
                    hdim = nd // params["num_state_vars"]
                state_block_size = hdim * params["num_state_vars"]

                # # --- get the process noise ---
                # pos, gs_model, L_C = compute_Q_err_random_fields(hdim, params["total_state_param_vars"], params["sig_Q"], Q_rho, len_scale)

                # process_noise = []
                for ens in range(params["Nens"]):
                    # model_kwargs.update({"ens_id": ens})
                    data = model_module.initialize_ensemble(ens,**model_kwargs)
                
                    # iterate over the data and update the ensemble
                    for key, value in data.items():
                        ensemble_vec[indx_map[key],ens] = value

                    # --->
                    # noise = compute_noise_random_fields(ens, hdim, pos, gs_model, params["total_state_param_vars"], L_C)
                    # ensemble_vec[:,ens] += noise
                    #----->
                    _time_init_noise_generation = MPI.Wtime()
                    N_size = params["total_state_param_vars"] * hdim
                    # noise = generate_pseudo_random_field_1d(N_size,np.sqrt(Lx*Ly), len_scale, verbose=True)
                    model_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":params["total_state_param_vars"]})
                    noise = generate_enkf_field(**model_kwargs)
                    time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                    ensemble_vec[:,ens] += alpha*noise
                    # for ii, sig in enumerate(params["sig_Q"]):
                    #     if ii <=params["num_state_vars"]:
                    #         start_idx = ii * hdim
                    #         end_idx = start_idx + hdim
                    #         ensemble_vec[start_idx:end_idx, ens] += noise[start_idx:end_idx] * sig

                    # enkf_parallel_io.write_forecast(0, ensemble_vec[:,ens], ens)
                    # enkf_parallel_io.datasets[0][:, ens] = ensemble_vec[:,ens]
                    # print(f"[ICESEE] Rank {rank_world}: Ensemble initialization completed and written to disk with norm {np.linalg.norm(ensemble_vec)}")
                
                shape_ens = np.array(ensemble_vec.shape,dtype=np.int32)
                # print(f"[ICESEE] Rank {rank_world}: Ensemble initialization completed for all members.")
                
    
            else:
                ensemble_vec = np.empty((params["nd"],params["Nens"]),dtype=np.float64)
                shape_ens = np.empty(2,dtype=np.int32)
                # pos, gs_model, L_C

            _time_init_file_writing = MPI.Wtime()    
            # scatter  enkf_parallel_io.nd_local_world of the ensemble to all processors
            localshape = enkf_parallel_io.nd_local_world
            all_local_shapes = comm_world.gather(localshape)
            if rank_world == 0:
                counts_rows = np.array(all_local_shapes) 
                displacement_rows = np.insert(np.cumsum(counts_rows), 0, 0)[0:-1]
                counts_rows = counts_rows * params["Nens"]
                displacement_rows = displacement_rows * params["Nens"]
            else:
                counts_rows = None
                displacement_rows = None
            
            local_ensemble = np.empty((localshape, params["Nens"]), dtype=np.float64)
            comm_world.Scatterv([ensemble_vec, counts_rows, displacement_rows, MPI.DOUBLE], local_ensemble, root=0)
            enkf_parallel_io.datasets[0][localshape, :] = local_ensemble
            time_init_file_writing += MPI.Wtime() - _time_init_file_writing

        comm_world.Barrier()
        _time_init_ensemble_mean_computation = MPI.Wtime()
        # enkf_parallel_io.compute_forecast_mean_chunked(0)
        enkf_parallel_io.compute_forecast_mean_chunked_v2(k=0,flag="initial")
        # ens_mean = enkf_parallel_io.compute_forecast_mean(0)
        # ens_mean = .datasets[0][:, :].mean(axis=1)
        time_init_ensemble_mean_computation += MPI.Wtime() - _time_init_ensemble_mean_computation

        # now reset the model_nprocs
        if rank_world == 0:
            diff = total_cores - base_total_procs 
            if diff >= 0:
                # split the diff amaongest all processors
                min_model_nprocs = max(model_nprocs-1, 1) 
                if model_kwargs.get('ICESEE_PERFORMANCE_TEST') or env_flag("ICESEE_PERFORMANCE_TEST", default=False):
                    model_nprocs = params.get("model_nprocs", 1)
                else:
                    model_nprocs = max(min_model_nprocs, model_nprocs + (diff // size_world))
            else:
                model_nprocs = model_nprocs

        model_nprocs = comm_world.bcast(model_nprocs, root=0)
        model_kwargs.update({'model_nprocs': model_nprocs})

    else:
        if rank_world == 0:
            print("[ICESEE] Initializing the ensemble ...")
        
        if params["default_run"] and size_world > params["Nens"]:
            # debug
            sub_shape = model_kwargs['dim_list'][sub_rank]
            model_kwargs.update({"statevec_ens":np.zeros((sub_shape, params["Nens"]))})

            model_kwargs.update({"ens_id": color, "rank": sub_rank, "color": color, "comm": subcomm})

            # ensemble_vec, shape_ens  = model_module.initialize_ensemble_debug(color,**model_kwargs)
            # ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], params['Nens'], comm_world, root=0)
            # parallel_write_full_ensemble_from_root(0, ens_mean, params,ensemble_vec,comm_world)
            # -----------------------------------------------------

            ens = color
            # model_kwargs.update({"statevec_ens":np.zeros((model_kwargs['global_shape'], params["Nens"]))})
            initialilaized_state = model_module.initialize_ensemble(ens,**model_kwargs)
            # ensemble_vec, shape_ens = gather_and_broadcast_data_default_run(initialilaized_state, subcomm, sub_rank, comm_world, rank_world, params)
            # ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], params['Nens'], comm_world, root=0)
            # parallel_write_full_ensemble_from_root(0, ens_mean, params,ensemble_vec,comm_world)
            # ensemble_vec = BM.bcast(ensemble_vec, comm_world)
            
            initial_data = {key: subcomm.gather(value, root=0) for key, value in initialilaized_state.items()}
            key_list = list(initial_data.keys())
            state_keys = key_list[:params["num_state_vars"]]
            if sub_rank == 0:
                # for key in initial_data:
                for key in key_list:
                    initial_data[key] = np.hstack(initial_data[key])
                    if model_kwargs["joint_estimation"] or params["localization_flag"]:
                        hdim = initial_data[key].shape[0] // params["total_state_param_vars"]
                    else:
                        hdim = initial_data[key].shape[0] // params["num_state_vars"]
                    state_block_size = hdim*params["num_state_vars"]
                    full_block_size = hdim*params["total_state_param_vars"]
                    # if key in state_keys:
                        # noise = np.random.normal(0, 0.1, state_block_size)
                        # Q_err = np.eye(state_block_size) * params["sig_Q"] ** 2
                        # Q_err = np.eye(state_block_size) * 0.01 ** 2
                    if model_kwargs.get("random_fields",False):
                        Q_err = np.zeros((full_block_size,full_block_size))
                        for i, sig in enumerate(params["sig_Q"]):
                            start_idx = i *hdim
                            end_idx = start_idx + hdim
                            Q_err[start_idx:end_idx,start_idx:end_idx] = np.eye(hdim) * sig ** 2

                        # noise = multivariate_normal.rvs(mean=np.zeros(state_block_size), cov=Q_err)
                        _time_init_noise_generation = MPI.Wtime()
                        noise = compute_noise_random_fields(ens, hdim, pos, gs_model, params["total_state_param_vars"], L_C)
                        time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                        # initial_data[key][:state_block_size] += noise[:state_block_size]
                        # noise = noise / np.max(np.abs(noise))
                        initial_data[key] += alpha*noise
                    else:
                        N_size = params["total_state_param_vars"] * hdim
                        _time_init_noise_generation = MPI.Wtime()
                        model_kwargs.update({"ii_sig": None, "hdim":hdim, "num_vars":params["total_state_param_vars"]})
                        noise = generate_enkf_field(**model_kwargs)
                        time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                        initial_data[key] += noise
                        # for ii, sig in enumerate(params["sig_Q"]):
                        #     start_idx = ii *hdim
                        #     end_idx = start_idx + hdim
                        #     initial_data[key][start_idx:end_idx] += noise[start_idx:end_idx]*sig
                    
                # stack all variables together into a single array
                stacked = np.hstack([initial_data[key] for key in initialilaized_state.keys()])
                shape_ens = np.array(stacked.shape,dtype=np.int32)
            else:
                shape_ens = np.empty(2,dtype=np.int32)

            # broadcast the shape of the initialized ensemble
            shape_ens = comm_world.bcast(shape_ens, root=0)

            if sub_rank != 0:
                stacked = np.empty(shape_ens,dtype=np.float64)

            all_init = comm_world.gather(stacked if sub_rank == 0 else None, root=0)

            if rank_world == 0:
                all_init = [arr for arr in all_init if isinstance(arr, np.ndarray)]
                ensemble_vec = np.column_stack(all_init)
                # print(f"[ICESEE] Shape of the ensemble: {ensemble_vec.shape}")
            else:
                ensemble_vec = np.empty((model_kwargs["global_shape"],params["Nens"]),dtype=np.float64)
            
            time_init_ensemble_mean_computation = MPI.Wtime()
            ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], params['Nens'], comm_world, root=0)
            time_init_ensemble_mean_computation += MPI.Wtime() - _time_init_ensemble_mean_computation

            _time_init_file_writing = MPI.Wtime()
            parallel_write_full_ensemble_from_root(0, ens_mean, model_kwargs,ensemble_vec,comm_world)
            time_init_file_writing += MPI.Wtime() - _time_init_file_writing
            
        elif params["sequential_run"]:
            comm_world.Barrier()
            sub_shape = model_kwargs['dim_list'][rank_world]
            model_kwargs.update({"statevec_ens":np.zeros([model_kwargs["global_shape"], params["Nens"]]),
                                "statevec_ens_mean":np.zeros([model_kwargs["global_shape"], model_kwargs.get("nt",params["nt"]) + 1]),
                                "statevec_ens_full":np.zeros([model_kwargs["global_shape"], params["Nens"], model_kwargs.get("nt",params["nt"]) + 1]),
                                "statevec_bg":np.zeros([model_kwargs["global_shape"], model_kwargs.get("nt",params["nt"]) + 1])})
            ensemble_bg, ensemble_vec, ensemble_vec_mean, ensemble_vec_full = model_module.initialize_ensemble(**model_kwargs)

            # gather from every rank to rank 0
            gathered_ensemble = comm_world.gather(ensemble_vec[:sub_shape,:], root=0)
            if rank_world == 0:
                ensemble_vec = np.vstack(gathered_ensemble)
                print(f"[ICESEE] Shape of the ensemble: {ensemble_vec.shape}")
                ensemble_vec_mean[:,0] = np.mean(ensemble_vec, axis=1)
                ensemble_vec_full[:,:,0] = ensemble_vec
            else:
                ensemble_vec = np.empty((model_kwargs["global_shape"],params["Nens"]),dtype=np.float64)
                ensemble_vec_mean = np.empty((model_kwargs["global_shape"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)
                ensemble_vec_full = np.empty((model_kwargs["global_shape"],params["Nens"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)

            # else:
            #     ensemble_bg = np.empty((model_kwargs["global_shape"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)
            #     ensemble_vec = np.empty((model_kwargs["global_shape"],params["Nens"]),dtype=np.float64)
            #     ensemble_vec_mean = np.empty((model_kwargs["global_shape"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)
            #     ensemble_vec_full = np.empty((model_kwargs["global_shape"],params["Nens"],model_kwargs.get("nt",params["nt"])+1),dtype=np.float64)

            # # Bcast the ensemble
            # comm_world.Bcast(ensemble_bg, root=0)
            comm_world.Bcast(ensemble_vec, root=0)
            comm_world.Bcast(ensemble_vec_mean, root=0)
            comm_world.Bcast(ensemble_vec_full, root=0)

            # hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
            # print(f"[ICESEE] rank: {rank_world}, subrank: {sub_rank}, min ensemble: {np.min(ensemble_vec[hdim,:])}, max ensemble: {np.max(ensemble_vec[hdim,:])}")

    if params.get("default_run", False):
        return model_kwargs, None, time_init_noise_generation, \
               time_init_ensemble_mean_computation,time_init_file_writing, \
                None, None, None, None
    else:
        return model_kwargs, ensemble_vec, time_init_noise_generation, \
               time_init_ensemble_mean_computation, time_init_file_writing, \
                shape_ens,ensemble_bg,  ensemble_vec_mean, ensemble_vec_full