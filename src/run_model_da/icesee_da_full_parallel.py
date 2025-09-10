# ==============================================================================
# @des: This script is for the fully parallelized ICESEE model-data assimilation
#       - It uses only the default MPI parallelization strategy.
#       - Uses parallel I/O batch I/O via both h5py and Zarr. (see EnKF_parallel_io.py)
# @date: 2025-09-8
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
import zarr
import numpy as np
from tqdm import tqdm 
import bigmpi4py as BM # BigMPI for large data transfer and communication
from mpi4py import MPI

# ==== ICESEE utility imports ========================================
from ICESEE.src.utils import tools, utils                                     # utility functions for the model 
from ICESEE.src.utils.utils import UtilsFunctions
from ICESEE.applications.supported_models import SupportedModels              # supported models for data assimilation routine
from ICESEE.src.utils.tools import icesee_get_index, display_timing_default,display_timing_verbose, save_all_data
from ICESEE.src.run_model_da._error_generation import compute_Q_err_random_fields, \
                              compute_noise_random_fields, \
                              generate_pseudo_random_field_1d, \
                              generate_pseudo_random_field_2D, \
                              generate_enkf_field

# --- call the ICESEE mpi parallel manager ---
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager
from ICESEE.src.parallelization._mpi_forecast_functions import parallel_forecast_step_default_full_parallel_run
from ICESEE.src.parallelization._mpi_generate_true_wrong_state import generate_true_wrong_state
from ICESEE.src.parallelization._mpi_ensemble_intialization import ensemble_initialization_full_parallel_run
from ICESEE.src.parallelization.EnKF_parallel_io import EnKF_fully_parallel_IO

# ======================== Run model with EnKF ========================
def icesee_model_data_assimilation_full_parallel(**model_kwargs): 
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
    data_path         = model_kwargs.get("data_path","_modelrun_datasets")      # data path


    # start the timer
    global_start_time = MPI.Wtime()

    # --- icesee mpi parallel manager ---------------------------------------------------
    # --- ensemble load distribution --
    rounds, color, sub_rank, sub_size, subcomm, subcomm_size_min, rank_world, size_world, comm_world, start, stop = ParallelManager().icesee_mpi_ens_distribution(params)
    model_kwargs.update({'size_world': size_world, 'comm_world': comm_world})

    # --- call curently supported model Class
    model_module = SupportedModels(model=model,comm=comm_world,verbose=params.get('verbose')).call_model()
    # pack the global communicator, the subcommunicator and other important parameters
    model_kwargs.update({"comm_world": comm_world, "subcomm": subcomm,
                            "rank_world": rank_world, "sub_rank": sub_rank,
                            "size_world": size_world, "sub_size": sub_size,
                            "rounds": rounds, "color": color,
                            "start": start, "stop": stop,
                            "subcomm_size_min": subcomm_size_min,
                            "model_module": model_module})

    # pack the global communicator and the subcommunicator
    model_kwargs.update({"comm_world": comm_world, "subcomm": subcomm})

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

    # --- intialize EnKF I/O handler class ---
    nd = model_kwargs.get("nd",params["nd"])
    nt = model_kwargs.get("nt",params["nt"])
    Nens = model_kwargs.get("Nens",params["Nens"])
    
    batch_size = model_kwargs.get("batch_size",100)
    # serial_file_creation = model_kwargs.get("serial_file_creation",True)
    serial_file_creation = True
    enkf_parallel_io = EnKF_fully_parallel_IO('icesee_enkf', nd, Nens, nt, subcomm, comm_world, \
                                             params, serial_file_creation, base_path=_modelrun_datasets, \
                                             batch_size=batch_size)
    # Update model_kwargs with the EnKF I/O handler
    model_kwargs.update({"enkf_parallel_io": enkf_parallel_io})

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
        effective_model_nprocs = max(1, np.floor(effective_model_nprocs * scale_factor)) 

    # update model_kwargs with the effective model_nprocs
    model_kwargs.update({'model_nprocs': effective_model_nprocs,
                            "total_cores": total_cores,
                            "base_total_procs": base_total_procs,
                        })

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
    # model_kwargs =  generate_synthetic_observations(**model_kwargs)
    synthetic_obs_zarr_path=f"{_modelrun_datasets}/synthetic_observations.zarr"
    error_R_zarr_path=f"{_modelrun_datasets}/error_R.zarr"
    model_kwargs.update({'synthetic_obs_zarr_path': synthetic_obs_zarr_path, 'error_R_zarr_path': error_R_zarr_path})
    tobserve, m_obs = enkf_parallel_io._create_synthetic_observations(**model_kwargs)
    model_kwargs.update({'tobserve': tobserve, 'm_obs': m_obs})
    # --- time generation of synthetic observations ---
    time_generation_synthetic_obs = MPI.Wtime() - time_generation_synthetic_obs

    comm_world.Barrier()
    #  --- generate the H file
    rank = comm_world.Get_rank()
    if rank == 0:
        print("[ICESEE] Generating H matrix and saving to Zarr...")
        H_matrix_zarr_path = "output/H_matrix.zarr"
        model_kwargs.update({'H_matrix_zarr_path': H_matrix_zarr_path})
        enkf_parallel_io.H_matrix(**model_kwargs)
                
    # --- Initialize the ensemble ---------------------------------------------------
    comm_world.Barrier()
    Q_rho     = model_kwargs.get("Q_rho")
    len_scale = model_kwargs.get("length_scale")
    hdim  = params["nd"] // params["total_state_param_vars"]
    model_kwargs.update({"hdim": hdim, "Q_rho": Q_rho, "len_scale": len_scale})

        # --- get the process noise --->
    if params.get("use_random_fields", False):
        pos, gs_model, L_C = compute_Q_err_random_fields(hdim, params["total_state_param_vars"], params["sig_Q"], Q_rho, len_scale)
        model_kwargs.update({"pos": pos, "gs_model": gs_model, "L_C": L_C})
    
    # -- time ensemble initialization ---
    time_ensemble_initialization = MPI.Wtime()
    # call the ensemble_initialization function
    model_kwargs, ensemble_vec, time_init_noise_generation, \
    time_init_ensemble_mean_computation, time_init_file_writing, \
    shape_ens,ensemble_bg,  ensemble_vec_mean, ensemble_vec_full = ensemble_initialization_full_parallel_run(**model_kwargs)
    # --- time ensemble initialization ---
    time_ensemble_initialization = MPI.Wtime() - time_ensemble_initialization

    # get updated model_nprocs
    model_nprocs = model_kwargs.get("model_nprocs", 1)

    # --- Define filter flags
    EnKF_flag   = re.match(r"\AEnKF\Z", filter_type, re.IGNORECASE)
    DEnKF_flag  = re.match(r"\ADEnKF\Z", filter_type, re.IGNORECASE)
    EnRSKF_flag = re.match(r"\AEnRSKF\Z", filter_type, re.IGNORECASE)
    EnTKF_flag  = re.match(r"\AEnTKF\Z", filter_type, re.IGNORECASE)

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
    # --- timing intializations
    time_forecast_step = 0.0
    time_analysis_step = 0.0
    time_forecast_noise_generation = 0.0
    time_forecast_file_writing = 0.0
    time_analysis_file_writing = 0.0
    time_forecast_ensemble_mean_generation = 0.0

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

    #--- generate inital noise
    if params.get("use_random_fields", False):
        # with h5py.File(_synthetic_obs, 'r') as f:
        #     error_R = f['error_R'][:]
        #     Cov_obs = np.cov(error_R)
        #  --- get the observation noise ---
        pos_obs, gs_model_obs, L_C_obs = compute_Q_err_random_fields(hdim, params["total_state_param_vars"], params["sig_obs"], Q_rho, len_scale)
    else:
        N_size = params["total_state_param_vars"] * hdim
        # noise = generate_pseudo_random_field_1d(N_size,np.sqrt(Lx*Ly), len_scale, verbose=0)
        model_kwargs.update({"ii_sig": None, "hdim":hdim, "num_vars":params["total_state_param_vars"]})
        # noise = generate_enkf_field(**model_kwargs)
        noise = generate_enkf_field(None, np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)
    
    # synchronize all processes before starting the time loop
    comm_world.Barrier()

    for k in range(model_kwargs.get("nt",params["nt"])):

        model_kwargs.update({"k": k, "km":km, "alpha": alpha, "rho": rho, "tau": tau, "dt": dt,"n": n})
        model_kwargs.update({"generate_enkf_field": generate_enkf_field}) #save the function to generate the enkf field

        if re.match(r"\AMPI_model\Z", parallel_flag, re.IGNORECASE):      
            # -- time forecast step ---
            _time_forecast_step = MPI.Wtime()

            # get the state block size
            ndim = nd//params["total_state_param_vars"]
            state_block_size = ndim*params["num_state_vars"]

            # load all needed parameters and variables into model_kwargs
            model_kwargs.update({"_modelrun_datasets": _modelrun_datasets,
                                "alpha": alpha, 
                                "rho": rho, 
                                "dt": dt, 
                                "Lx": Lx, 
                                "Ly": Ly, 
                                "len_scale": len_scale,
                                "model_module": model_module,
                                "time_forecast_step": time_forecast_step,
                                "time_analysis_step": time_analysis_step,
                                "time_forecast_noise_generation": time_forecast_noise_generation,
                                "time_forecast_file_writing": time_forecast_file_writing,
                                "time_analysis_file_writing": time_analysis_file_writing,
                                "time_forecast_ensemble_mean_generation": time_forecast_ensemble_mean_generation,
                                "state_block_size": state_block_size, "noise": noise, "rng": None, "rank_seed": None,})
            
            if params["default_run"]:
                # call the parallel_forecast_step_default_run function
                model_kwargs = parallel_forecast_step_default_full_parallel_run(**model_kwargs)
                time_forecast_step = model_kwargs.get("time_forecast_step", 0.0)
                time_forecast_noise_generation = model_kwargs.get("time_forecast_noise_generation", 0.0)
                time_forecast_file_writing = model_kwargs.get("time_forecast_file_writing", 0.0)
                time_forecast_ensemble_mean_generation = model_kwargs.get("time_forecast_ensemble_mean_generation", 0.0)
            
                comm_world.Barrier()
                # --- end time forecast step
                time_forecast_step += MPI.Wtime() - _time_forecast_step

                # ===== Global analysis step =====
                if model_kwargs.get('global_analysis', True) or model_kwargs.get('local_analysis', False):
                    # -- time global analysis step ---
                    _time_analysis_step = MPI.Wtime()
    #                 obs_index = model_kwargs["obs_index"]
    #                 if (km < params["number_obs_instants"]) and (k+1 == obs_index[km]):
    #                     #

    #                     # comm_world.Barrier()
    #                     if rank_world == 0:

    #                         ndim = ensemble_vec.shape[0]//params["total_state_param_vars"]  
    #                         state_block_size = ndim*params["num_state_vars"]
                        
    #                         # -------------
    #                         # H = UtilsFunctions(params, ensemble_vec).JObs_fun(ensemble_vec.shape[0]) 
    #                         # h = UtilsFunctions(params, ensemble_vec).Obs_fun # observation operator

    #                         # compute the observation covariance matrix
    #                         # Cov_obs = params["sig_obs"][k+1]**2 * np.eye(2*params["number_obs_instants"]+1)
    #                         # Cov_obs = error_R[:,k+1]**2 * np.eye(2*params["number_obs_instants"]+1)

    #                         # --- vector of measurements
    #                         with h5py.File(_synthetic_obs, 'r') as f:
    #                             hu_obs  = f['hu_obs'][:]
    #                             error_R = f['R'][:]
    #                             # Cov_obs = np.cov(error_R)
    #                             Cov_obs = np.zeros(error_R.shape)

    #                         d = UtilsFunctions(params, ensemble_vec).Obs_fun(hu_obs[:,km])
    #                         model_kwargs.update({"error_R": error_R}) # store the error covariance matrix
    #                         #  -------------

    #                         # get parameter
    #                         # parameter_estimated = ensemble_vec[state_block_size:,:]
    #                         eta = 0.0 # trend term
    #                         beta = np.ones(nd)
    #                         # ensemble_vec[state_block_size:,:] = ensemble_vec[state_block_size:,:] + (eta + beta)*model_kwargs.get("dt",params["dt"]) + np.sqrt(model_kwargs.get("dt",params["dt"])) * alpha*rho*q0[state_block_size:]

    #                         if EnKF_flag:
    #                             # compute the X5 matrix
    #                             X5,analysis_vec_ij = EnKF_X5(k,ensemble_vec, Cov_obs, Nens, d, model_kwargs,UtilsFunctions)
    #                             # X5 = EnKF_X5(Cov_obs, Nens, D, HA, Eta, d)
    #                             y_i = np.sum(X5, axis=1)
    #                             # ensemble_vec_mean[:,k+1] = (1/Nens)*(ensemble_vec @ y_i.reshape(-1,1)).ravel()
    #                             time_analysis_mean_generation = MPI.Wtime()
    #                             ens_mean = (1/Nens)*(ensemble_vec @ y_i.reshape(-1,1)).ravel()
    #                             time_analysis_mean_generation = MPI.Wtime() - time_analysis_mean_generation

    #                         elif DEnKF_flag:
    #                             # compute the X5 matrix
    #                             X5,X5prime = DEnKF_X5(k,ensemble_vec, Cov_obs, Nens, d, model_kwargs,UtilsFunctions)
    #                             # y_i = np.sum(X5, axis=1)
    #                             # ens_mean = (1/Nens)*(ensemble_vec @ y_i.reshape(-1,1)).ravel()
    #                             # H = UtilsFunctions(params, ensemble_vec).JObs_fun(ensemble_vec.shape[0])
    #                             # Cov_model = np.cov(ensemble_vec)
    #                             # ens_mean = np.mean(ensemble_vec, axis=1)
    #                             # diff = (ensemble_vec -np.tile(ens_mean.reshape(-1,1),Nens) )
    #                             # Cov_model = 1/(Nens-1) * (diff @ diff.T)
    #                             # epsilon = 1e-6
    #                             # inv_matrix = np.linalg.pinv(H @ Cov_model @ H.T + Cov_obs + epsilon * np.eye(Cov_obs.shape[0]))
    #                             # KalGain = Cov_model @ H.T @ inv_matrix
    #                             # X5prime = KalGain@(d - np.dot(H, ens_mean))
    #                             # ens_mean = ens_mean + X5prime
    #                             # print(f"[ICESEE] X5prime shape: {X5prime.shape}")
    #                             analysis_vec_ij = None
    #                     else:
    #                         X5 = np.empty((Nens, Nens))
    #                         time_analysis_mean_generation = 0.0
    #                         analysis_vec_ij = None
    #                         smb_scale = 0.0
    #                         if DEnKF_flag:
    #                             ens_mean = np.empty((nd, 1))

    #                     if model_kwargs.get('local_analysis', False):
    #                         shape_ens = ensemble_vec.shape
    #                         ens_mean = ParallelManager().compute_mean_matrix_from_root(analysis_vec_ij, shape_ens[0], params['Nens'], comm_world, root=0)
    #                         parallel_write_full_ensemble_from_root(k+1,ens_mean, model_kwargs,analysis_vec_ij,comm_world)
                        
    #                     # smb_scale = comm_world.bcast(smb_scale, root=0)
    #                     smb_scale = 1.0

    #                     with h5py.File(_synthetic_obs, 'r', driver='mpio', comm=comm_world) as f:
    #                         hu_obs  = f['hu_obs'][:]

    #                     # fetch the upper and lower bounds for every paramerter from observed data
    #                     ndim = hu_obs.shape[0]//params["total_state_param_vars"]
    #                     state_block_size = ndim*params["num_state_vars"]
    #                     bounds = []
    #                     for i, var in enumerate(model_kwargs["params_vec"]):
    #                         bound_idx = (params["num_state_vars"] + i) * ndim
    #                         bound_idx_end = bound_idx + ndim

    #                         param_slice = hu_obs[bound_idx:bound_idx_end, km]
    #                         param_min = np.min(param_slice)
    #                         param_max = np.max(param_slice)

    #                         bounds.append(np.array([param_min, param_max]))

    #                     # pack the bunds into model_kwargs
    #                     model_kwargs.update({"bounds": bounds})
                            

    #                     # call the analysis update function
    #                     if EnKF_flag:
    #                         time_analysis_mean_generation, time_analysis_file_writing = analysis_enkf_update(k,ens_mean,ensemble_vec, \
    #                                                                                                          shape_ens, X5, time_analysis_mean_generation, \
    #                                                                                                             time_analysis_file_writing, analysis_vec_ij,\
    #                                                                                                         UtilsFunctions,model_kwargs,smb_scale)
    #                     elif DEnKF_flag:
    #                         model_kwargs.update({"DEnKF_flag": True})
    #                         analysis_Denkf_update(k,ens_mean,ensemble_vec, shape_ens, X5,UtilsFunctions,model_kwargs,smb_scale)
    #                         # analysis_enkf_update(k,ens_mean,ensemble_vec, shape_ens, X5, analysis_vec_ij,UtilsFunctions,model_kwargs,smb_scale)
                    
    #                     # update the observation index
    #                     km += 1
    #                     # hu_obs[state_block_size:,:] *= smb_scale
    #                     del hu_obs
    #                     gc.collect()
                        
    #                     # --- end time analysis step ---
    #                     time_analysis_step += MPI.Wtime() - _time_analysis_step

    #                 else: 
    #                     # if Nens < size_world:
                        
    #                     _time_forecast_file_writing = MPI.Wtime()

    #                     parallel_write_full_ensemble_from_root(k+1,ens_mean, model_kwargs,ensemble_vec,comm_world)

    #                     # --time forecast file writing ---
    #                     _time_forecast_file_writing = MPI.Wtime() - _time_forecast_file_writing
    #                     time_forecast_file_writing += _time_forecast_file_writing
    #                     time_forecast_step = time_forecast_step + _time_forecast_file_writing
    #                     del ensemble_vec; gc.collect()
    #                         # parallel_write_full_ensemble_from_root(ensemble_vec,ensemble_vec_full,comm_world,k)
                    

    #             # ======= Local analyais step =======
    #             if model_kwargs.get('local_analysis', False):
    #                 # --- compute the local X5 for each horizontal grid point ---
    #                 pass

    #         # -------------------------------------------------- end of cases 2 & 3 --------------------------------------------

    #     # update the progress bar
    #     if rank_world == 0:
    #         pbar.update(1)

    # # close the progress bar
    # if rank_world == 0:
    #     pbar.close()
    # # comm_world.Barrier()

    # # ====== load data to be written to file ======
    # # print("[ICESEE] Saving data ...")
    # save_all_data(
    #     enkf_params=model_kwargs['enkf_params'],
    #     nofilter=True,
    #     t=model_kwargs["t"], b_io=np.array([b_in,b_out]),
    #     Lxy=np.array([Lx,Ly]),nxy=np.array([nx,ny]),
    #     # ensemble_true_state=ensemble_true_state,
    #     # ensemble_nurged_state=ensemble_nurged_state, 
    #     obs_max_time=np.array([params["obs_max_time"]]),
    #     obs_index=model_kwargs["obs_index"],
    #     # w=hu_obs,
    #     run_mode= np.array([params["execution_flag"]])
    # )
    enkf_parallel_io.close()

    # # ─────────────────────────────────────────────────────────────
    # #  End Timer and Aggregate Elapsed Time Across Processors
    # # ─────────────────────────────────────────────────────────────
    # # --total elapsed time
    # global_end_time = MPI.Wtime()
    # global_elapsed_time = global_end_time - global_start_time
    # # Reduce elapsed time across all processors (sum across ranks)
    # total_elapsed_time = comm_world.allreduce(global_elapsed_time, op=MPI.SUM)
    # total_wall_time = comm_world.allreduce(global_elapsed_time, op=MPI.MAX)

    # # -- timing true and wrong state generation
    # true_wrong_time = comm_world.allreduce(time_generation_true_and_wrong_state, op=MPI.MAX)

    # # -- timing ensemble initialization
    # ensemble_init_time = comm_world.allreduce(time_ensemble_initialization, op=MPI.MAX)

    # # -- timing forecast step
    # forecast_step_time = comm_world.allreduce(time_forecast_step, op=MPI.MAX)

    # # -- timing forecast noise generation
    # forecast_noise_time = comm_world.allreduce(time_forecast_noise_generation, op=MPI.MAX)

    # # -- timing analysis step
    # analysis_step_time = comm_world.allreduce(time_analysis_step, op=MPI.MAX)

    # # -- total assimilation time = ensemble init + forecast step + analysis step
    # assimilation_time = ensemble_init_time + forecast_step_time + analysis_step_time

    # # --- time forecast file writing ---
    # forecast_file_time = comm_world.allreduce(time_forecast_file_writing, op=MPI.MAX)

    # # --- time analysis file writing ---
    # analysis_file_time = comm_world.allreduce(time_analysis_file_writing, op=MPI.MAX)

    # # total file writing time initialization file writing + forecast file writing + analysis file writing
    # init_file_time = comm_world.allreduce(time_init_file_writing, op=MPI.MAX)
    # total_file_time = init_file_time + forecast_file_time + analysis_file_time

    # # Display elapsed time on rank 0
    # comm_world.Barrier()
    # if rank_world == 0:
    #     verbose = model_kwargs.get("verbose", False)
    #     # if verbose:
    #     if True:
    #          display_timing_verbose(
    #         computational_time=total_elapsed_time,
    #         wallclock_time=total_wall_time,
    #         true_wrong_time=true_wrong_time,
    #         assimilation_time=assimilation_time,
    #         forecast_step_time=forecast_step_time,
    #         analysis_step_time=analysis_step_time,
    #         ensemble_init_time=ensemble_init_time,
    #         init_file_time=init_file_time,
    #         forecast_file_time=forecast_file_time,
    #         analysis_file_time=analysis_file_time,
    #         total_file_time=total_file_time,
    #         forecast_noise_time=forecast_noise_time, comm=comm_world
    #     )
    #     else:
    #         display_timing_default(total_elapsed_time, total_wall_time)
    # else:
    #     None


