# ==============================================================================
# @des: This file contains run functions for any model with data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2024-11-4
# @author: Brian Kyanjo
# ==============================================================================
    
# ==== Imports ========================================================
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

# ==== ICESEE utility imports ========================================
from ICESEE.src.utils import tools, utils                                     # utility functions for the model 
from ICESEE.src.utils.utils import UtilsFunctions
from ICESEE.src.EnKF.python_enkf.EnKF import EnsembleKalmanFilter as EnKF     # Ensemble Kalman Filter
from ICESEE.applications.supported_models import SupportedModels              # supported models for data assimilation routine
from ICESEE.src.utils.tools import icesee_get_index, display_timing, \
                                save_all_data
from ICESEE.src.run_model_da._error_generation import compute_Q_err_random_fields, \
                              compute_noise_random_fields, \
                              generate_pseudo_random_field_1d, \
                              generate_pseudo_random_field_2D, \
                              generate_enkf_field
from ICESEE.src.run_model_da._localization_functions import gaspari_cohn, localization 

# ======================== Run model with EnKF ========================
def icesee_model_data_assimilation(**model_kwargs): 
    """ General function to run any kind of model with the Ensemble Kalman Filter """

    # --- unpack the data assimilation arguments
    filter_type       = model_kwargs.get("filter_type", "EnKF")      # filter type
    model             = model_kwargs.get("model_name",None)          # model name
    parallel_flag     = model_kwargs.get("parallel_flag",False)      # parallel flag
    params            = model_kwargs.get("params",None)              # parameters
    Q_err             = model_kwargs.get("Q_err",None)               # process noise
    commandlinerun    = model_kwargs.get("commandlinerun",None)      # run through the terminal
    Lx, Ly            = model_kwargs.get("Lx",1.0), model_kwargs.get("Ly",1.0)
    nx, ny            = model_kwargs.get("nx",0.2), model_kwargs.get("ny",0.2)
    b_in, b_out       = model_kwargs.get("b_in",0.0), model_kwargs.get("b_out",0.0) 

    # --- call the ICESEE mpi parallel manager ---
    if re.match(r"\AMPI_model\Z", parallel_flag, re.IGNORECASE):
        from mpi4py import MPI
        from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager
        from ICESEE.src.run_model_da._mpi_analysis_functions import analysis_enkf_update, EnKF_X5, DEnKF_X5, \
                                            gather_and_broadcast_data_default_run, \
                                            parallel_write_full_ensemble_from_root, \
                                            parallel_write_data_from_root_2D, \
                                            parallel_write_vector_from_root, \
                                            analysis_Denkf_update, \
                                            parallel_forecast_step_default_run
        
        from ICESEE.src.run_model_da._mpi_generate_true_wrong_state import generate_true_wrong_state
        from ICESEE.src.run_model_da._mpi_generate_synthetic_observations import generate_synthetic_observations
        
        # start the timer
        global_start_time = MPI.Wtime()

        # --- icesee mpi parallel manager ---------------------------------------------------
        # --- ensemble load distribution --
        rounds, color, sub_rank, sub_size, subcomm, subcomm_size, rank_world, size_world, comm_world, start, stop = ParallelManager().icesee_mpi_ens_distribution(params)
        model_kwargs.update({'size_world': size_world, 'comm_world': comm_world})

        # --- call curently supported model Class
        model_module = SupportedModels(model=model,comm=comm_world,verbose=params.get('verbose')).call_model()


        # pack the global communicator, the subcommunicator and other important parameters
        model_kwargs.update({"comm_world": comm_world, "subcomm": subcomm,
                             "rank_world": rank_world, "sub_rank": sub_rank,
                             "size_world": size_world, "sub_size": sub_size,
                             "rounds": rounds, "color": color,
                             "start": start, "stop": stop,
                             "subcomm_size": subcomm_size,
                             "model_module": model_module})

        # --- check if the modelrun dataset directory is present ---
        _modelrun_datasets = model_kwargs.get("data_path",None)
        if rank_world == 0 and not os.path.exists(_modelrun_datasets):
            # cretate the directory
            os.makedirs(_modelrun_datasets, exist_ok=True)

        comm_world.Barrier()
        # --- file_names
        _true_nurged   = f'{ _modelrun_datasets}/true_nurged_states.h5'
        _synthetic_obs = f'{ _modelrun_datasets}/synthetic_obs.h5'

        # --update model_kwargs with the file names
        model_kwargs.update({"true_nurged_file": _true_nurged, "synthetic_obs_file": _synthetic_obs})

        # --- initialize seed for reproducibility ---
        ParallelManager().initialize_seed(comm_world, base_seed=0)

        # fetch model nprocs
        model_nprocs = params.get("model_nprocs", 1)

        # set modeel_nprocs adaptively
        total_cores = os.cpu_count()
        base_total_procs = size_world + (size_world * model_nprocs)  # MPI + MATLAB processes
        diff = total_cores - base_total_procs  # Available or deficit cores

        # Dynamic process allocation
        if rank_world == 0:
            # Prioritize rank 0: Allocate extra cores or handle deficit
            if diff >= 0:
                # Extra cores available: give rank 0 up to 2x model_nprocs or more
                extra_procs = min(diff, model_nprocs * 2)  # Cap at 2x base for safety
                effective_model_nprocs = model_nprocs + extra_procs
            else:
                # Deficit: Maintain base model_nprocs or slightly reduce
                effective_model_nprocs = max(1, model_nprocs + (diff // size_world))
        else:
            # Other ranks: Minimize MATLAB processes, ensure at least 1
            if diff >= 0:
                effective_model_nprocs = model_nprocs
            else:
                effective_model_nprocs = max(1, model_nprocs + (diff // size_world))

        # Ensure total processes don’t exceed cores
        total_matlab_procs = effective_model_nprocs if rank_world == 0 else effective_model_nprocs * (size_world - 1)
        total_procs = size_world + total_matlab_procs
        if total_procs > total_cores:
            # Scale down proportionally
            scale_factor = total_cores / total_procs
            min_model_nprocs = max(model_nprocs-1, 1)  # Ensure at least 1 process
            effective_model_nprocs = max(min_model_nprocs, np.floor(effective_model_nprocs * scale_factor))

        # update model_kwargs with the effective model_nprocs
        model_kwargs.update({'model_nprocs': effective_model_nprocs})

        # --- Generate True and Nurged States -------------------------------------------------------------------
        # -- time generation of true state ----
        time_generation_true_and_wrong_state = MPI.Wtime()
        # call the generate_true_wrong_state function
        model_kwargs = generate_true_wrong_state(**model_kwargs)
        # --- time generation of true state and nurged state ---
        time_generation_true_and_wrong_state = MPI.Wtime() - time_generation_true_and_wrong_state

        comm_world.Barrier()

        # --- Generate the Synthetic ObservationsObservations ---------------------------------------------------
        # --- time generation of synthetic observations ---
        time_generation_synthetic_obs = MPI.Wtime()
        # call the generate_synthetic_observations function
        model_kwargs =  generate_synthetic_observations(**model_kwargs)
        # --- time generation of synthetic observations ---
        time_generation_synthetic_obs = MPI.Wtime() - time_generation_synthetic_obs

        # --- Initialize the ensemble -------------------------------------------------------------------------
        comm_world.Barrier()
        Q_rho     = model_kwargs.get("Q_rho")
        len_scale = model_kwargs.get("length_scale")
        hdim  = params["nd"] // params["total_state_param_vars"]
        # --- get the process noise --->
        # pos, gs_model, L_C = compute_Q_err_random_fields(hdim, params["total_state_param_vars"], params["sig_Q"], Q_rho, len_scale)
        
        # -- time ensemble initialization ---
        time_ensemble_initialization = MPI.Wtime()
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
                vecs, indx_map, dim_per_proc = icesee_get_index(model_kwargs["statevec_ens"], **model_kwargs)
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
                    ensemble_id = color + (round_id * subcomm_size)
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
                        time_init_noise_generation = MPI.Wtime()
                        noise = generate_enkf_field(None, np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)
                        time_init_noise_generation = MPI.Wtime() - time_init_noise_generation
                        ensemble_vec[:, ens] += noise
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

                    # ensemble = ensemble_vec
                    # for ens in range(params["Nens"]):
                        # add process noise
                        # if model_kwargs["joint_estimation"] or params["localization_flag"]:
                        #     hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
                        # else:
                        #     hdim = ensemble_vec.shape[0] // params["num_state_vars"]
                        # state_block_size = hdim * params["num_state_vars"]

                        # --->
                        # noise = compute_noise_random_fields(ens, hdim, pos, gs_model, params["total_state_param_vars"], L_C)
                        # ensemble_vec[:,ens] += noise
                        #----->
                        time_init_noise_generation = MPI.Wtime()
                        N_size = params["total_state_param_vars"] * hdim
                        # noise = generate_pseudo_random_field_1d(N_size,np.sqrt(Lx*Ly), len_scale, verbose=True)
                        noise = generate_enkf_field(None,np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)
                        time_init_noise_generation = MPI.Wtime() - time_init_noise_generation
                        ensemble_vec[:,ens] += noise
                        # for ii, sig in enumerate(params["sig_Q"]):
                        #     if ii <=params["num_state_vars"]:
                        #         start_idx = ii * hdim
                        #         end_idx = start_idx + hdim
                        #         ensemble_vec[start_idx:end_idx, ens] += noise[start_idx:end_idx] * sig


                        # -----------------------------
                        # full_block_size = hdim * params["total_state_param_vars"]
                        # Q_err = np.zeros((full_block_size,full_block_size))
                        # for i, sig in enumerate(params["sig_Q"]):
                        #     start_idx = i *hdim
                        #     end_idx = start_idx + hdim
                        #     Q_err[start_idx:end_idx,start_idx:end_idx] = np.eye(hdim) * sig ** 2

                        # # print(f"[ICESEE] [Q_err] Q_err shape: {Q_err.shape}, Q_err: {Q_err[:10,:10]}")


                        # # noise = multivariate_normal.rvs(mean=np.zeros(state_block_size), cov=Q_err[:state_block_size,:state_block_size])
                        # noise = multivariate_normal.rvs(mean=np.zeros(full_block_size), cov=Q_err)
                        # ensemble_vec[:,ens] += noise

                        # # print(f"[ICESEE] [Debug] Ensemble vector shape: {ensemble_vec.shape}, noise shape: {noise.shape}")
                        # ------------------------------

                        # add a spread to the smb
                        # if model_kwargs["joint_estimation"] or params["localization_flag"]:
                            # ensemble_vec[state_block_size:,ens] = ensemble_vec[state_block_size:,ens] + np.diag(Q_err[state_block_size:,state_block_size:])

                        
                        
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
                time_init_ensemble_mean_computation = MPI.Wtime()
                ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], params['Nens'], comm_world, root=0)
                time_init_ensemble_mean_computation = MPI.Wtime() - time_init_ensemble_mean_computation

                # ---time file writing ---
                time_init_file_writing = MPI.Wtime()
                parallel_write_full_ensemble_from_root(0, ens_mean, model_kwargs,ensemble_vec,comm_world)
                time_init_file_writing = MPI.Wtime() - time_init_file_writing

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
                            time_init_noise_generation = MPI.Wtime()
                            noise = compute_noise_random_fields(ens, hdim, pos, gs_model, params["total_state_param_vars"], L_C)
                            time_init_noise_generation = MPI.Wtime() - time_init_noise_generation
                            # initial_data[key][:state_block_size] += noise[:state_block_size]
                            # noise = noise / np.max(np.abs(noise))
                            initial_data[key] += noise
                        else:
                            N_size = params["total_state_param_vars"] * hdim
                            time_init_noise_generation = MPI.Wtime()
                            noise = generate_enkf_field(None,np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)
                            time_init_noise_generation = MPI.Wtime() - time_init_noise_generation
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

                time_init_file_writing = MPI.Wtime()
                parallel_write_full_ensemble_from_root(0, ens_mean, model_kwargs,ensemble_vec,comm_world)
                time_init_file_writing = MPI.Wtime() - time_init_file_writing
                
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

        # --- time ensemble initialization ---
        time_ensemble_initialization = MPI.Wtime() - time_ensemble_initialization

        # --- get the ensemble size
        nd, Nens = ensemble_vec.shape

        if params["even_distribution"]:
            ensemble_local = copy.deepcopy(ensemble_vec[:,start:stop])
                
        # --- row vector load distribution ---   
        # local_rows, start_row, end_row = ParallelManager().icesee_mpi_row_distribution(ensemble_vec, params)
        comm_world.Barrier()
        parallel_manager = None # debugging flag for now
        
    else:
        parallel_manager = None

        # --- call curently supported model Class
        model_module = SupportedModels(model=model,verbose=params.get('verbose')).call_model()

        # --- get the ensemble size
        # nd, Nens = ensemble_vec.shape
        nd = params["nd"]
        Nens = params["Nens"]
        size_world = 1
        rank_world = 0
        sub_rank = 0
        color = 0

        _modelrun_datasets = model_kwargs.get("data_path",None)
        if rank_world == 0 and not os.path.exists(_modelrun_datasets):
            # cretate the directory
            os.makedirs(_modelrun_datasets, exist_ok=True)

        # comm_world.Barrier()
        # --- file_names
        _true_nurged   = f'{ _modelrun_datasets}/true_nurged_states.h5'
        _synthetic_obs = f'{ _modelrun_datasets}/synthetic_obs.h5'

        # -- generate the true and nurged states



        # --- generate synthetic observations


        # --- initialize the ensemble ---
        if params["even_distribution"] or (params["default_run"] and size_world <= params["Nens"]):
            if rank_world == 0:
                print("[ICESEE] Initializing the ensemble ...")
                model_kwargs.update({"statevec_ens":np.zeros([params["nd"], params["Nens"]])})
                
                # get the ensemble matrix   
                vecs, indx_map, dim_per_proc = icesee_get_index(model_kwargs["statevec_ens"], **model_kwargs)
                ensemble_vec = np.zeros_like(model_kwargs["statevec_ens"])

                if model_kwargs["joint_estimation"] or params["localization_flag"]:
                    hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
                else:
                    hdim = ensemble_vec.shape[0] // params["num_state_vars"]
                state_block_size = hdim * params["num_state_vars"]

                for ens in range(params["Nens"]):
                    # model_kwargs.update({"ens_id": ens})
                    data = model_module.initialize_ensemble(ens,**model_kwargs)
                
                    # iterate over the data and update the ensemble
                    for key, value in data.items():
                        ensemble_vec[indx_map[key],ens] = value

                    N_size = params["total_state_param_vars"] * hdim
                    noise = generate_enkf_field(None,np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)

                    # for ii, sig in enumerate(params["sig_Q"]):
                    #     start_idx = ii *hdim
                    #     end_idx = start_idx + hdim
                    #     ensemble_vec[start_idx:end_idx,ens] += noise[start_idx:end_idx]*sig

                    ensemble_vec[:,ens] += noise

                shape_ens = np.array(ensemble_vec.shape,dtype=np.int32)

            else:
                ensemble_vec = np.empty((params["nd"],params["Nens"]),dtype=np.float64)

            shape_ens = comm_world.bcast(shape_ens, root=0)

            ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], params['Nens'], comm_world, root=0)
            parallel_write_full_ensemble_from_root(0, ens_mean, model_kwargs,ensemble_vec,comm_world)

    # --- hdim based on nd or global_shape ---
    if params["even_distribution"] or (params["default_run"] and size_world <= params["Nens"]):
        if model_kwargs["joint_estimation"] or params["localization_flag"]:
            hdim = nd // params["total_state_param_vars"]
        else:
            hdim = nd // params["num_state_vars"]
    else:   
        if model_kwargs["joint_estimation"] or params["localization_flag"]:
            hdim = model_kwargs["global_shape"] // params["total_state_param_vars"]
        else:
            hdim = model_kwargs["global_shape"] // params["num_state_vars"]
        
    state_block_size = hdim * params["num_state_vars"]

    # --- compute the process noise covariance matrix ---
    # check if scalar or matrix
    if isinstance(params["sig_Q"], float):
        nd = hdim*params["total_state_param_vars"]
        # params["nd"] = nd
        Q_err = np.eye(nd) * params["sig_Q"] ** 2
    else:
        nd = hdim*params["total_state_param_vars"]
        # params["nd"] = nd
        # Q_err = np.diag(params["sig_Q"] ** 2)
        Q_err = np.zeros((nd,nd))
        for i, sig in enumerate(params["sig_Q"]):
            start_idx = i *hdim
            end_idx = start_idx + hdim
            Q_err[start_idx:end_idx,start_idx:end_idx] = np.eye(hdim) * sig ** 2

        # with h5py.File(_synthetic_obs, 'r') as f:
        #     error_R = f['error_R'][:]
        #     Cov_obs = np.cov(error_R)
        #  --- get the observation noise ---
        # pos_obs, gs_model_obs, L_C_obs = compute_Q_err_random_fields(hdim, params["total_state_param_vars"], params["sig_obs"], Q_rho, len_scale) #TODO:will start from here
    
    # save the process noise to the model_kwargs dictionary
    model_kwargs.update({"Q_err": Q_err})
    

    # --- Define filter flags
    EnKF_flag   = re.match(r"\AEnKF\Z", filter_type, re.IGNORECASE)
    DEnKF_flag  = re.match(r"\ADEnKF\Z", filter_type, re.IGNORECASE)
    EnRSKF_flag = re.match(r"\AEnRSKF\Z", filter_type, re.IGNORECASE)
    EnTKF_flag  = re.match(r"\AEnTKF\Z", filter_type, re.IGNORECASE)

    # get the grid points
    if params.get("localization_flag", False):
        #  for both localization and joint estimation
        # - apply Gaspari-Cohn localization to only state variables [h,u,v] in [h,u,v,smb]
        # - for parameters eg. smb and others, don't apply localization
        # if model_kwargs["joint_estimation"]:
        # get state variables indices
        num_state_vars = params["num_state_vars"]
        num_params = params["num_param_vars"]
        # get the the inital smb
        # smb_init = ensemble_vec[num_state_vars*hdim:,:]
        inflation_factor = params["inflation_factor"] #TODO: store this, for localization debuging
        
        if True:
            # --- call the localization function (with adaptive localization) ---
            state_size = params["total_state_param_vars"]*hdim
            adaptive_localization = False   
            if not adaptive_localization:
                x_points = np.linspace(0, model_kwargs["Lx"], model_kwargs["nx"]+1)
                y_points = np.linspace(0, model_kwargs["Ly"], model_kwargs["ny"]+1)
                grid_x, grid_y = np.meshgrid(x_points, y_points)

                grid_points = np.vstack((grid_x.ravel(), grid_y.ravel())).T

                # Adjust grid if n_points != nx * ny (interpolating for 425 points)
                n_points = hdim
                missing_rows = n_points - grid_points.shape[0]
                if missing_rows > 0:
                    last_row = grid_points[-1]  # Get the last available row
                    extrapolated_rows = np.tile(last_row, (missing_rows, 1))  # Repeat last row
                    grid_points = np.vstack([grid_points, extrapolated_rows])  # Append extrapolated rows

                dist_matrix = distance_matrix(grid_points, grid_points) 

                # Normalize distance matrix
                L = 2654
                r_matrix = dist_matrix / L
            else:
                loc_matrix = localization(Lx,Ly,nx, ny, hdim, params["total_state_param_vars"], Nens, state_size)
    

    # --- Initialize the EnKF class ---
    EnKFclass = EnKF(parameters=params, parallel_manager=parallel_manager, parallel_flag = parallel_flag)

    # tqdm progress bar
    # Initialize progress bar on the root process
    if rank_world == 0:
        nt = model_kwargs.get("nt", params["nt"])
        print(f"[ICESEE] Launching {model} with data assimilation using the {filter_type} filter across {size_world} MPI ranks.")
        pbar = tqdm(
            total=nt,
            desc=f"[ICESEE] Assimilation progress ({size_world} ranks)",
            position=0,
            leave=True,
            dynamic_ncols=True
        )

    # ==== Time loop =======================================================================================
    # comm_world.Barrier()
    # --- timing intializations
    time_forecast_step = 0.0
    time_analysis_step = 0.0
    time_forecast_noise_generation = 0.0
    time_forecast_file_writing = 0.0
    time_analysis_file_writing = 0.0


    # specified decorrelation length scale, tau,
    min_tau = 200
    max_tau = 500
    dt  = model_kwargs.get("dt",params["dt"])
    tau = max(max_tau,max(min_tau, dt))

    # tau = max(model_kwargs.get("dt",params["dt"]),10)
    alpha = 1 - dt/tau
    # make sure  0=<alpha<1
    if alpha <= 0 or alpha > 1:
        alpha = 0.5
    n = model_kwargs.get("nt",params["nt"])
    # rho = np.sqrt((1-alpha**2)/(dt*(n - 2*alpha - n*alpha**2 + 2*alpha**(n+1))))
    rho = np.sqrt((1/dt)*((1-alpha)**2)*(1/(n - (2*alpha) - (n*alpha**2) + (2*alpha**(n+1)))))
    params_analysis_0 = np.zeros((2, Nens))
    km = 0
    for k in range(model_kwargs.get("nt",params["nt"])):

        model_kwargs.update({"k": k, "km":km, "alpha": alpha, "rho": rho, "tau": tau, "dt": dt,"n": n})
        model_kwargs.update({"generate_enkf_field": generate_enkf_field}) #save the function to generate the enkf field

        # background step
        # ensemble_bg = model_module.background_step(k,ensemble_bg, hdim, **model_kwargs)

        # save a copy of initial ensemble
        # ensemble_init = ensemble_vec.copy()

        if re.match(r"\AMPI_model\Z", parallel_flag, re.IGNORECASE):                                   
            
            # -- time forecast step ---
            _time_forecast_step = MPI.Wtime()

            # load all needed parameters and variables into model_kwargs
            model_kwargs.update({"rounds":rounds, 
                                "color": color, 
                                "subcomm_size": subcomm_size, 
                                "subcomm": subcomm, 
                                "comm_world": comm_world, 
                                "rank_world": rank_world, 
                                "sub_rank": sub_rank, 
                                "_modelrun_datasets": _modelrun_datasets,
                                "alpha": alpha, 
                                "rho": rho, 
                                "dt": dt, 
                                "Lx": Lx, 
                                "Ly": Ly, 
                                "len_scale": len_scale,
                                "model_module": model_module,
                                "UtilsFunctions": UtilsFunctions,
                                "generate_enkf_field": generate_enkf_field,
                                "icesee_get_index": icesee_get_index})
            
            # === Four approaches of forecast step mpi parallelization ===
            # --- case 1: Each forecast runs squentially using all available processors
            if params.get("sequential_run", False):
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

                # Analysis step
                obs_index = model_kwargs["obs_index"]
                if (km < params["number_obs_instants"]) and (k+1 == obs_index[km]):
                #     local_ensemble_centered = ensemble_local -  np.mean(ensemble_local, axis=1).reshape(-1,1)  # Center data
                    if EnKF_flag or DEnKF_flag:
                        diff = ens_T[:nd,:] - ensemble_vec_mean[:,k+1].reshape(-1,1)
                        Cov_model = diff @ diff.T / (Nens - 1)
                    elif EnRSKF_flag or EnTKF_flag:
                        diff = ens_T[:nd,:] - ensemble_vec_mean[:,k+1].reshape(-1,1)
                        Cov_model = diff / (Nens - 1)
                
                    # localization
                    if params.get("localization_flag", False):
                        # try with the localization matrix
                        cutoff_distance = 6000

                        # rho = np.zeros_like(Cov_model)
                        rho = np.ones_like(Cov_model)
                        # for j in range(Cov_model.shape[0]):
                        #     for i in range(Cov_model.shape[1]):
                        #         rad_x = np.abs(X[j] - X[i])
                        #         rad_y = np.abs(Y[j] - Y[i])
                        #         rad = np.sqrt(rad_x**2 + rad_y**2)
                        #         rad = rad/cutoff_distance
                        #         rho[j,i] = gaspari_cohn(rad)

                        Cov_model = rho * Cov_model

                        # if EnKF_flag or DEnKF_flag:
                    analysis  = EnKF(Observation_vec=  UtilsFunctions(params, ensemble_vec).Obs_fun(hu_obs[:nd,km]), 
                                            Cov_obs=params["sig_obs"][k+1]**2 * np.eye(2*params["number_obs_instants"]+1), \
                                            Cov_model= Cov_model, \
                                            Observation_function=UtilsFunctions(params, ensemble_vec).Obs_fun, \
                                            Obs_Jacobian=UtilsFunctions(params, ensemble_vec).JObs_fun, \
                                            parameters=  params,\
                                            parallel_flag=   parallel_flag)
                    # compute the analysis ensemble
                    if EnKF_flag:
                        ens_T[:nd,:], Cov_model = analysis.EnKF_Analysis(ens_T[:nd,:])
                    elif DEnKF_flag:
                        ens_T[:nd,:], Cov_model = analysis.DEnKF_Analysis(ens_T[:nd,:])
                    elif EnRSKF_flag:
                        ens_T[:nd,:], Cov_model = analysis.EnRSKF_Analysis(ens_T[:nd,:])
                    elif EnTKF_flag:
                        ens_T[:nd,:], Cov_model = analysis.EnTKF_Analysis(ens_T[:nd,:])
                    else:
                        raise ValueError("Filter type not supported")
                    
                    # update the ensemble mean
                    ensemble_vec_mean[:,k+1] = np.mean(ens_T[:nd,:], axis=1)

                    # update observation index
                    km += 1

                    # inflate the ensemble
                    ens_T = UtilsFunctions(params, ens_T[:nd,:]).inflate_ensemble(in_place=True)
                            
                    # ensemble_vec = copy.deepcopy(ens_T[:nd,:])
                    ensemble_vec = ens_T[:,:]
            
                # save the ensemble
                ensemble_vec_full[:,:,k+1] = ensemble_vec[:nd,:]
                
                # before exiting the time loop, we have to gather data from all processors
                if k == model_kwargs.get("nt",params["nt"]) - 1:
                    # we are interested in ensemble_vec_full, ensemble_vec_mean, ensemble_bg
                    gathered_ens_vec_mean = comm_world.allgather(ensemble_vec_mean)
                    gathered_ens_vec_full = comm_world.allgather(ensemble_vec_full)
                    if rank_world == 0:
                        # print(f"[ICESEE] Ensemble mean shape: {[arr.shape for arr in gathered_ens_vec_mean]}")
                        ensemble_vec_mean = np.vstack(gathered_ens_vec_mean)
                        ensemble_vec_full = np.vstack(gathered_ens_vec_full)
                        print(f"[ICESEE] Ensemble mean shape: {ensemble_vec_mean.shape}")
                    else:
                        ensemble_vec_mean = np.empty((shape_ens[0],model_kwargs.get("nt",params["nt"])+1), dtype=np.float64)
                        ensemble_vec_full = np.empty((shape_ens[0],Nens,model_kwargs.get("nt",params["nt"])+1), dtype=np.float64)

            # ------------------------------------------------- end of case 1 -------------------------------------------------

            # --- cases 2 & 3 ---
            # case 2: Form batches of sub-communicators and distribute resources among them
            #          - only works for Nens >= size_world
            # case 3: form Nens sub-communicators and distribute resources among them
            #          - only works for size_world > Nens
            #          - even distribution and load balancing leading to performance improvement
            #          - best for size_world/Nens is a whole number

            
            if params["default_run"]:
                # --- case 2: Form batches of sub-communicators and distribute resources among them ---
                if Nens >= size_world:
                    # ensemble_vec, shape_ens = parallel_forecast_step_default_run(**model_kwargs)
                   
                    # store results for each round      
                    ens_list = []
                    for round_id in range(rounds):
                        ensemble_id = color + round_id * subcomm_size  # Global ensemble index
                        model_kwargs.update({'ens_id': ensemble_id, 'comm': subcomm})

                        if ensemble_id < Nens:  # Only process valid ensembles
                            # print(f"[ICESEE] Rank {rank_world} processing ensemble {ensemble_id} in round {round_id + 1}/{rounds}")

                            # Ensure all ranks in the subcommunicator are synchronized before running
                            subcomm.Barrier()
                            ens = ensemble_id
                            # ---- read from file ----
                            input_file = f"{_modelrun_datasets}/icesee_ensemble_data.h5"
                            # with h5py.File(input_file, "r", driver="mpio", comm=subcomm) as f:
                            with h5py.File(input_file, "r") as f:
                                ensemble_vec = f["ensemble"][:,ens,k]
                            # ---- end of read from file ----

                            # Call the forecast step function
                            # hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
                            # print(f"[ICESEE] Rank: {rank_world}, min ensemble: {np.min(ensemble_vec[:hdim])}, max ensemble: {np.max(ensemble_vec[:hdim])}")

                            updated_state = model_module.forecast_step_single(ensemble=ensemble_vec,**model_kwargs)

                            #fetch the updated state
                            vecs, indx_map, dim_per_proc = icesee_get_index(ensemble_vec, **model_kwargs)
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
                            if k == 0:
                                # noise = compute_noise_random_fields(ens, hdim, pos, gs_model, params["total_state_param_vars"], L_C)
                                N_size = params["total_state_param_vars"] * hdim
                                # noise = generate_pseudo_random_field_1d(N_size,np.sqrt(Lx*Ly), len_scale, verbose=0)
                                noise = generate_enkf_field(None,np.sqrt(Lx*Ly),hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)

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
                            if k == 0:
                                N_size = params["total_state_param_vars"] * hdim
                                noise = generate_enkf_field(ens,np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)

                            noise_all = []
                            q0 = []
                            for ii, sig in enumerate(params["sig_Q"]):
                                if ii <=params["num_state_vars"]:
                                    # W = np.random.normal(0, 1, hdim)
                                    # W = generate_pseudo_random_field_1d(hdim,np.sqrt(Lx*Ly), len_scale, verbose=0)
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
                time_forecast_ensemble_generation = MPI.Wtime()
                ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], Nens, comm_world, root=0)
                time_forecast_ensemble_generation = MPI.Wtime() - time_forecast_ensemble_generation

                # --- end time forecast step
                time_forecast_step += MPI.Wtime() - _time_forecast_step

                # ===== Global analysis step =====
                if model_kwargs.get('global_analysis', True) or model_kwargs.get('local_analysis', False):
                    # -- time global analysis step ---
                    _time_analysis_step = MPI.Wtime()
                    obs_index = model_kwargs["obs_index"]
                    if (km < params["number_obs_instants"]) and (k+1 == obs_index[km]):
                        # *- parallelize the getting Eta, D and HA steps
                        # if Nens >= size_world:
                        #     with h5py.File(input_file, "r", driver="mpio", comm=comm_world) as f:
                        #         hu_obs = f["hu_obs"][:]
                        #         local_Nens = Nens // size_world
                        #         remainder = Nens % size_world
                        #         start = rank_world*local_Nens + min(rank_world, remainder)
                        #         if rank_world < remainder:
                        #             local_Nens += 1
                        #         stop = start + local_Nens
                        #         ensemble_vec = f["ensemble"][:,start:stop,k+1]
                        #         d = UtilsFunctions(params, ensemble_vec).Obs_fun(hu_obs[:,km])
                        #         Cov_obs = params["sig_obs"][k+1]**2 * np.eye(2*params["number_obs_instants"]+1)
                        #         Eta = np.zeros((d.shape[0], local_Nens))
                        #         D = np.zeros_like(Eta)
                        #         HA = np.zeros_like(Eta)
                        #         for i, ens in enumerate(range(start, stop)):
                        #             Eta[:,i] = np.random.multivariate_normal(mean=np.zeros(d.shape[0]), cov=Cov_obs)
                        #             D[:,i] = d + Eta[:,i]
                        #             HA[:,i] = UtilsFunctions(params, ensemble_vec[:,i]).Obs_fun(ensemble_vec[:,i])

                        #         # gather the results
                        #         gathered_Eta = comm_world.gather(Eta, root=0)
                        #         gathered_D = comm_world.gather(D, root=0)
                        #         gathered_HA = comm_world.gather(HA, root=0)
                        #         if rank_world == 0:
                        #             Eta = np.hstack(gathered_Eta)
                        #             D = np.hstack(gathered_D)
                        #             HA = np.hstack(gathered_HA)
                        #         else:
                        #             Eta = np.zeros((d.shape[0], Nens))
                        #             D = np.zeros_like(Eta)
                        #             HA = np.zeros_like(Eta)
                        # else:
                        #     if rank_world < Nens:
                        #         with h5py.File(input_file, "r", driver="mpio", comm=comm_world) as f:
                        #             hu_obs = f["hu_obs"][:]
                        #             ensemble_vec = f["ensemble"][:,rank_world,k+1]
                        #             d = UtilsFunctions(params, ensemble_vec).Obs_fun(hu_obs[:,km])
                        #             Cov_obs = params["sig_obs"][k+1]**2 * np.eye(2*params["number_obs_instants"]+1)
                        #             Eta = np.random.multivariate_normal(mean=np.zeros(d.shape[0]), cov=Cov_obs)
                        #             D = d + Eta
                        #             HA = UtilsFunctions(params, ensemble_vec).Obs_fun(ensemble_vec)

                        #             # gather the results
                        #             gathered_Eta = comm_world.gather(Eta, root=0)
                        #             gathered_D = comm_world.gather(D, root=0)
                        #             gathered_HA = comm_world.gather(HA, root=0)
                        #             if rank_world == 0:
                        #                 Eta = np.hstack(gatheEta)
                        #                 D = np.hstack(gathered_D)
                        #                 HA = np.hstack(gathered_HA)
                        #             else:
                        #                 Eta = np.zeros((d.shape[0], Nens))
                        #                 D = np.zeros_like(Eta)
                        #                 HA = np.zeros_like(Eta)

                        #     else:
                        #         global_shape = model_kwargs["global_shape"]
                        #         Eta = np.zeros((global_shape, Nens))
                        #         D = np.zeros_like(Eta)
                        #         HA = np.zeros_like(Eta)

                        # comm_world.Barrier()
                        if rank_world == 0:

                            ndim = ensemble_vec.shape[0]//params["total_state_param_vars"]  
                            state_block_size = ndim*params["num_state_vars"]
                        
                            # -------------
                            # H = UtilsFunctions(params, ensemble_vec).JObs_fun(ensemble_vec.shape[0]) 
                            # h = UtilsFunctions(params, ensemble_vec).Obs_fun # observation operator

                            # compute the observation covariance matrix
                            # Cov_obs = params["sig_obs"][k+1]**2 * np.eye(2*params["number_obs_instants"]+1)
                            # Cov_obs = error_R[:,k+1]**2 * np.eye(2*params["number_obs_instants"]+1)

                            # --- vector of measurements
                            with h5py.File(_synthetic_obs, 'r') as f:
                                hu_obs  = f['hu_obs'][:]
                                error_R = f['R'][:]
                                # Cov_obs = np.cov(error_R)
                                Cov_obs = np.zeros(error_R.shape)

                            d = UtilsFunctions(params, ensemble_vec).Obs_fun(hu_obs[:,km])
                            model_kwargs.update({"error_R": error_R}) # store the error covariance matrix
                            #  -------------

                            # get parameter
                            # parameter_estimated = ensemble_vec[state_block_size:,:]
                            eta = 0.0 # trend term
                            beta = np.ones(nd)
                            # ensemble_vec[state_block_size:,:] = ensemble_vec[state_block_size:,:] + (eta + beta)*model_kwargs.get("dt",params["dt"]) + np.sqrt(model_kwargs.get("dt",params["dt"])) * alpha*rho*q0[state_block_size:]

                            if EnKF_flag:
                                # compute the X5 matrix
                                X5,analysis_vec_ij = EnKF_X5(k,ensemble_vec, Cov_obs, Nens, d, model_kwargs,UtilsFunctions)
                                # X5 = EnKF_X5(Cov_obs, Nens, D, HA, Eta, d)
                                y_i = np.sum(X5, axis=1)
                                # ensemble_vec_mean[:,k+1] = (1/Nens)*(ensemble_vec @ y_i.reshape(-1,1)).ravel()
                                time_analysis_mean_generation = MPI.Wtime()
                                ens_mean = (1/Nens)*(ensemble_vec @ y_i.reshape(-1,1)).ravel()
                                time_analysis_mean_generation = MPI.Wtime() - time_analysis_mean_generation

                            elif DEnKF_flag:
                                # compute the X5 matrix
                                X5,X5prime = DEnKF_X5(k,ensemble_vec, Cov_obs, Nens, d, model_kwargs,UtilsFunctions)
                                # y_i = np.sum(X5, axis=1)
                                # ens_mean = (1/Nens)*(ensemble_vec @ y_i.reshape(-1,1)).ravel()
                                # H = UtilsFunctions(params, ensemble_vec).JObs_fun(ensemble_vec.shape[0])
                                # Cov_model = np.cov(ensemble_vec)
                                # ens_mean = np.mean(ensemble_vec, axis=1)
                                # diff = (ensemble_vec -np.tile(ens_mean.reshape(-1,1),Nens) )
                                # Cov_model = 1/(Nens-1) * (diff @ diff.T)
                                # epsilon = 1e-6
                                # inv_matrix = np.linalg.pinv(H @ Cov_model @ H.T + Cov_obs + epsilon * np.eye(Cov_obs.shape[0]))
                                # KalGain = Cov_model @ H.T @ inv_matrix
                                # X5prime = KalGain@(d - np.dot(H, ens_mean))
                                # ens_mean = ens_mean + X5prime
                                # print(f"[ICESEE] X5prime shape: {X5prime.shape}")
                                analysis_vec_ij = None
                        else:
                            X5 = np.empty((Nens, Nens))
                            time_analysis_mean_generation = 0.0
                            analysis_vec_ij = None
                            smb_scale = 0.0
                            if DEnKF_flag:
                                ens_mean = np.empty((nd, 1))

                        if model_kwargs.get('local_analysis', False):
                            shape_ens = ensemble_vec.shape
                            ens_mean = ParallelManager().compute_mean_matrix_from_root(analysis_vec_ij, shape_ens[0], params['Nens'], comm_world, root=0)
                            parallel_write_full_ensemble_from_root(k+1,ens_mean, model_kwargs,analysis_vec_ij,comm_world)
                        
                        # smb_scale = comm_world.bcast(smb_scale, root=0)
                        smb_scale = 1.0

                        with h5py.File(_synthetic_obs, 'r', driver='mpio', comm=comm_world) as f:
                            hu_obs  = f['hu_obs'][:]

                        # fetch the upper and lower bounds for every paramerter from observed data
                        ndim = hu_obs.shape[0]//params["total_state_param_vars"]
                        state_block_size = ndim*params["num_state_vars"]
                        bounds = []
                        for i, var in enumerate(model_kwargs["params_vec"]):
                            bound_idx = (params["num_state_vars"] + i) * ndim
                            bound_idx_end = bound_idx + ndim

                            param_slice = hu_obs[bound_idx:bound_idx_end, km]
                            param_min = np.min(param_slice)
                            param_max = np.max(param_slice)

                            bounds.append(np.array([param_min, param_max]))

                        # pack the bunds into model_kwargs
                        model_kwargs.update({"bounds": bounds})
                            

                        # call the analysis update function
                        if EnKF_flag:
                            time_analysis_mean_generation, time_analysis_file_writing = analysis_enkf_update(k,ens_mean,ensemble_vec, \
                                                                                                             shape_ens, X5, time_analysis_mean_generation, \
                                                                                                                time_analysis_file_writing, analysis_vec_ij,\
                                                                                                            UtilsFunctions,model_kwargs,smb_scale)
                        elif DEnKF_flag:
                            model_kwargs.update({"DEnKF_flag": True})
                            analysis_Denkf_update(k,ens_mean,ensemble_vec, shape_ens, X5,UtilsFunctions,model_kwargs,smb_scale)
                            # analysis_enkf_update(k,ens_mean,ensemble_vec, shape_ens, X5, analysis_vec_ij,UtilsFunctions,model_kwargs,smb_scale)
                    
                        # update the observation index
                        km += 1
                        # hu_obs[state_block_size:,:] *= smb_scale
                        del hu_obs
                        gc.collect()
                        
                        # --- end time analysis step ---
                        time_analysis_step += MPI.Wtime() - _time_analysis_step

                    else: 
                        # if Nens < size_world:
                        
                        _time_forecast_file_writing = MPI.Wtime()

                        parallel_write_full_ensemble_from_root(k+1,ens_mean, model_kwargs,ensemble_vec,comm_world)

                        # --time forecast file writing ---
                        _time_forecast_file_writing = MPI.Wtime() - _time_forecast_file_writing
                        time_forecast_file_writing += _time_forecast_file_writing
                        time_forecast_step = time_forecast_step + _time_forecast_file_writing
                        del ensemble_vec; gc.collect()
                            # parallel_write_full_ensemble_from_root(ensemble_vec,ensemble_vec_full,comm_world,k)
                    

                # ======= Local analyais step =======
                if model_kwargs.get('local_analysis', False):
                    # --- compute the local X5 for each horizontal grid point ---
                    pass


            # -------------------------------------------------- end of cases 2 & 3 --------------------------------------------

            # --- case 4: Evenly distribute ensemble members among processors 
            #         - each processor runs a subset of ensemble members
            #         - best for size_world/Nens is a whole number and Nens >= size_world
            #         - size_world = 2^n where n is an integer
            if params["even_distribution"]:
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

                    # Analysis step
                    obs_index = model_kwargs["obs_index"]
                    if (km < params["number_obs_instants"]) and (k+1 == obs_index[km]):
                
                        if rank_world == 0:
    
                            H = UtilsFunctions(params, ensemble_vec).JObs_fun(ensemble_vec.shape[0])
                            h = UtilsFunctions(params, ensemble_vec).Obs_fun # observation operator

                            # compute the observation covariance matrix
                            Cov_obs = params["sig_obs"][k+1]**2 * np.eye(2*params["number_obs_instants"]+1)

                            # --- vector of measurements
                            d = UtilsFunctions(params, ensemble_vec).Obs_fun(hu_obs[:,km])

                            if EnKF_flag:
                                # compute the X5 matrix
                                X5 = EnKF_X5(ensemble_vec, Cov_obs, Nens, h, d)
                                y_i = np.sum(X5, axis=1)
                                ensemble_vec_mean[:,k+1] = (1/Nens)*(ensemble_vec @ y_i.reshape(-1,1)).ravel()
                                
                        else:
                            X5 = np.empty((Nens, Nens))

                        # clean the memory
                        del ensemble_local, gathered_ensemble; gc.collect()

                        # call the analysis update function
                        shape_ens = ensemble_vec.shape # get the shape of the ensemble
                        ensemble_vec = analysis_enkf_update(ensemble_vec, shape_ens, X5, comm_world)

                        # update the ensemble with observations instants
                        km += 1

                        # inflate the ensemble
                        # params["inflation_factor"] = inflation_factor
                        ensemble_vec = UtilsFunctions(params, ensemble_vec).inflate_ensemble(in_place=True)
                        # ensemble_vec = UtilsFunctions(params, ensemble_vec)._inflate_ensemble()
                    
                        # update the local ensemble
                        ensemble_local = copy.deepcopy(ensemble_vec[:,start:stop])

                    # Save the ensemble
                    if rank_world == 0:
                        ensemble_vec_full[:,:,k+1] = ensemble_vec
                    else:
                        ensemble_vec_full = np.empty((nd, Nens, model_kwargs.get("nt",params["nt"])+1), dtype=np.float64)

                    # free up memory
                    del ensemble_vec; gc.collect()
                else:
                    raise ValueError("Nens must be divisible by size_world and greater or equal to size_world. size_world must be a power of 2")

            # -------------------------------------------------- end of case 4 -------------------------------------------------    
        #  ====== Serial run ======
        else:
            input_file = f"{_modelrun_datasets}/icesee_ensemble_data.h5"
            with h5py.File(input_file, "r") as f:
                ensemble_vec = f["ensemble"][:,:,k]

            ensemble_vec = EnKFclass.forecast_step(ensemble_vec, \
                                               model_module.forecast_step_single, \
                                                Q_err, **model_kwargs)


            #  compute the ensemble mean
            ensemble_vec_mean[:,k+1] = np.mean(ensemble_vec, axis=1)

            # Analysis step
            obs_index = model_kwargs["obs_index"]
            if (km < params["number_obs_instants"]) and (k+1 == obs_index[km]):

                # Compute the model covariance
                diff = ensemble_vec - np.tile(ensemble_vec_mean[:,k+1].reshape(-1,1),Nens)
                if EnKF_flag or DEnKF_flag:
                    Cov_model = 1/(Nens-1) * diff @ diff.T
                elif EnRSKF_flag or EnTKF_flag:
                    Cov_model = 1/(Nens-1) * diff 

                # --- localization ---
                if params["localization_flag"]:
                    if not adaptive_localization:
                        # call the gahpari-cohn localization function
                        loc_matrix_spatial = gaspari_cohn(r_matrix)

                        # expand to full state space
                        loc_matrix = np.empty_like(Cov_model)
                        for var_i in range(params["total_state_param_vars"]):
                            for var_j in range(params["total_state_param_vars"]):
                                start_i, start_j = var_i * hdim, var_j * hdim
                                loc_matrix[start_i:start_i+hdim, start_j:start_j+hdim] = loc_matrix_spatial
                        
                        # apply the localization matrix
                        # Cov_model = loc_matrix * Cov_model
                        
                    Cov_model = loc_matrix * Cov_model

                    # inflate the top-left (smb h) and bottom-right (h smb) blocks of the covariance matrix 
                    state_block_size = num_state_vars*hdim
                    h_smb_block = Cov_model[:hdim,state_block_size:]
                    smb_h_block = Cov_model[state_block_size:,:hdim]

                    # apply the inflation factor
                    params["inflation_factor"] = 1.2
                    smb_h_block = UtilsFunctions(params, smb_h_block).inflate_ensemble(in_place=True)
                    h_smb_block = UtilsFunctions(params, h_smb_block).inflate_ensemble(in_place=True)

                    # update the covariance matrix
                    Cov_model[:hdim,state_block_size:] = h_smb_block
                    Cov_model[state_block_size:,:hdim] = smb_h_block

                # check if params["sig_obs"] is a scalar
                if isinstance(params["sig_obs"], (int, float)):
                    params["sig_obs"] = np.ones(model_kwargs.get("nt",params["nt"])+1) * params["sig_obs"]


                # Call the EnKF class for the analysis step
                analysis  = EnKF(Observation_vec=  UtilsFunctions(params, ensemble_vec).Obs_fun(hu_obs[:,km]), 
                                Cov_obs=params["sig_obs"][k+1]**2 * np.eye(2*params["number_obs_instants"]+1), \
                                Cov_model= Cov_model, \
                                Observation_function=UtilsFunctions(params, ensemble_vec).Obs_fun, \
                                Obs_Jacobian=UtilsFunctions(params, ensemble_vec).JObs_fun, \
                                parameters=  params,\
                                parallel_flag=   parallel_flag)
                
                # Compute the analysis ensemble
                if EnKF_flag:
                    ensemble_vec, Cov_model = analysis.EnKF_Analysis(ensemble_vec)
                elif DEnKF_flag:
                    ensemble_vec, Cov_model = analysis.DEnKF_Analysis(ensemble_vec)
                elif EnRSKF_flag:
                    ensemble_vec, Cov_model = analysis.EnRSKF_Analysis(ensemble_vec)
                elif EnTKF_flag:
                    ensemble_vec, Cov_model = analysis.EnTKF_Analysis(ensemble_vec)
                else:
                    raise ValueError("Filter type not supported")

                ensemble_vec_mean[:,k+1] = np.mean(ensemble_vec, axis=1)
                
                # update the ensemble with observations instants
                km += 1

                # inflate the ensemble
                ensemble_vec = UtilsFunctions(params, ensemble_vec).inflate_ensemble(in_place=True)
                # ensemble_vec = UtilsFunctions(params, ensemble_vec)._inflate_ensemble()
            
                # ensemble_vec_mean[:,k+1] = np.mean(ensemble_vec, axis=1)

            # Save the ensemble
            ensemble_vec_full[:,:,k+1] = ensemble_vec

        # update the progress bar
        if rank_world == 0:
            pbar.update(1)

    # close the progress bar
    if rank_world == 0:
        pbar.close()
    # comm_world.Barrier()

    # ====== load data to be written to file ======
    # print("[ICESEE] Saving data ...")
    if params["even_distribution"]:
        save_all_data(
            enkf_params=model_kwargs['enkf_params'],
            nofilter=True,
            t=model_kwargs["t"], b_io=np.array([b_in,b_out]),
            Lxy=np.array([Lx,Ly]),nxy=np.array([nx,ny]),
            ensemble_true_state=ensemble_true_state,
            ensemble_nurged_state=ensemble_nurged_state, 
            obs_max_time=np.array([params["obs_max_time"]]),
            obs_index=model_kwargs["obs_index"],
            w=hu_obs,
            run_mode= np.array([params["execution_flag"]])
        )

        # --- Save final data ---
        save_all_data(
            enkf_params=model_kwargs['enkf_params'],
            ensemble_vec_full=ensemble_vec_full,
            ensemble_vec_mean=ensemble_vec_mean,
            ensemble_bg=ensemble_bg
        )
    else:
        save_all_data(
            enkf_params=model_kwargs['enkf_params'],
            nofilter=True,
            t=model_kwargs["t"], b_io=np.array([b_in,b_out]),
            Lxy=np.array([Lx,Ly]),nxy=np.array([nx,ny]),
            # ensemble_true_state=ensemble_true_state,
            # ensemble_nurged_state=ensemble_nurged_state, 
            obs_max_time=np.array([params["obs_max_time"]]),
            obs_index=model_kwargs["obs_index"],
            # w=hu_obs,
            run_mode= np.array([params["execution_flag"]])
        )

    # ─────────────────────────────────────────────────────────────
    #  End Timer and Aggregate Elapsed Time Across Processors
    # ─────────────────────────────────────────────────────────────
    # # --- global elapsed time ---
    # global_end_time = MPI.Wtime()
    # global_elapsed_time = global_end_time - global_start_time
    # # Reduce elapsed time across all processors (sum across ranks)
    # total_elapsed_time = comm_world.allreduce(global_elapsed_time, op=MPI.SUM)
    # total_wall_time = comm_world.allreduce(global_elapsed_time, op=MPI.MAX)

    # Display elapsed time on rank 0
    # comm_world.Barrier()
    # if rank_world == 0:
    #     display_timing(total_elapsed_time, total_wall_time)
    # else:
    #     None

    # ─────────────────────────────────────────────────────────────
    #  End Timer and Aggregate Elapsed Time Across Processors
    # ─────────────────────────────────────────────────────────────
    # --total elapsed time
    global_end_time = MPI.Wtime()
    global_elapsed_time = global_end_time - global_start_time
    # Reduce elapsed time across all processors (sum across ranks)
    total_elapsed_time = comm_world.allreduce(global_elapsed_time, op=MPI.SUM)
    total_wall_time = comm_world.allreduce(global_elapsed_time, op=MPI.MAX)

    # -- timing true and wrong state generation
    true_wrong_time = comm_world.allreduce(time_generation_true_and_wrong_state, op=MPI.MAX)

    # -- timing ensemble initialization
    ensemble_init_time = comm_world.allreduce(time_ensemble_initialization, op=MPI.MAX)

    # -- timing forecast step
    forecast_step_time = comm_world.allreduce(time_forecast_step, op=MPI.MAX)

    # -- timing forecast noise generation
    forecast_noise_time = comm_world.allreduce(time_forecast_noise_generation, op=MPI.MAX)

    # -- timing analysis step
    analysis_step_time = comm_world.allreduce(time_analysis_step, op=MPI.MAX)

    # -- total assimilation time = ensemble init + forecast step + analysis step
    assimilation_time = ensemble_init_time + forecast_step_time + analysis_step_time

    # --- time forecast file writing ---
    forecast_file_time = comm_world.allreduce(time_forecast_file_writing, op=MPI.MAX)

    # --- time analysis file writing ---
    analysis_file_time = comm_world.allreduce(time_analysis_file_writing, op=MPI.MAX)

    # total file writing time initialization file writing + forecast file writing + analysis file writing
    init_file_time = comm_world.allreduce(time_init_file_writing, op=MPI.MAX)
    total_file_time = init_file_time + forecast_file_time + analysis_file_time

    comm_world.Barrier()
    if rank_world == 0:
        display_timing(
            computational_time=total_elapsed_time,
            wallclock_time=total_wall_time,
            true_wrong_time=true_wrong_time,
            assimilation_time=assimilation_time,
            forecast_step_time=forecast_step_time,
            analysis_step_time=analysis_step_time,
            ensemble_init_time=ensemble_init_time,
            init_file_time=init_file_time,
            forecast_file_time=forecast_file_time,
            analysis_file_time=analysis_file_time,
            total_file_time=total_file_time,
            forecast_noise_time=forecast_noise_time, comm=comm_world
        )
    else:
        None

