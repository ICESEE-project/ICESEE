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
from ICESEE.src.utils.tools import icesee_get_index, display_timing_default,display_timing_verbose, save_all_data
from ICESEE.src.EnKF._localization_inflation import  LocalizationInflationUtils
from ICESEE.src.EnKF._generate_synthetic_observations import generate_synthetic_observations
from ICESEE.src.EnKF._generate_true_wrong_state import generate_true_wrong_state
from ICESEE.src.EnKF._ensemble_initialization import ensemble_initialization
from ICESEE.src.run_model_da._error_generation import compute_Q_err_random_fields, \
                              compute_noise_random_fields, \
                              generate_pseudo_random_field_1d, \
                              generate_pseudo_random_field_2D, \
                              generate_enkf_field

# ======================== Run model with EnKF ========================
def icesee_model_data_assimilation_serial(**model_kwargs): 
    """ General function to run any kind of model with the Ensemble Kalman Filter """

    # print(f"model_kwargs: {model_kwargs}")

    # --- unpack the data assimilation arguments
    filter_type       = model_kwargs.get("filter_type", "EnKF")      # filter type
    model             = model_kwargs.get("model_name",None)          # model name
    parallel_flag     = model_kwargs.get("parallel_flag",False)      # parallel flag
    params            = model_kwargs.get("params",None)              # parameters
    Q_err             = model_kwargs.get("Q_err",None)               # process noise
    commandlinerun    = model_kwargs.get("commandlinerun",None)      # run through the terminal
    Lx, Ly            = model_kwargs.get("Lx",1.0), model_kwargs.get("Ly",1.0)
    nx, ny            = model_kwargs.get("nx",1), model_kwargs.get("ny",1)
    b_in, b_out       = model_kwargs.get("b_in",0.0), model_kwargs.get("b_out",0.0) 

    #  call the gaspari function and localization function
    localization = LocalizationInflationUtils(params).localization
    gaspari_cohn = LocalizationInflationUtils(params).gaspari_cohn

    parallel_manager = None

    # --- call curently supported model Class
    model_module = SupportedModels(model=model,verbose=params.get('verbose')).call_model()

    # --- get the ensemble size
    # nd, Nens = ensemble_vec.shape
    nd = model_kwargs.get('nd', params.get("nd", 3))
    Nens = params["Nens"]
    size_world = 1
    rank_world = 0
    sub_rank = 0
    color = 0
    sub_size = 1
    rounds = 1
    start = 0
    stop = 0
    subcomm_size_min = 1
    comm_world = None
    subcomm = None

    # pack the global communicator, the subcommunicator and other important parameters
    model_kwargs.update({"comm_world": comm_world, "subcomm": subcomm,
                            "rank_world": rank_world, "sub_rank": sub_rank,
                            "size_world": size_world, "sub_size": sub_size,
                            "rounds": rounds, "color": color,
                            "start": start, "stop": stop,
                            "subcomm_size_min": subcomm_size_min, "dim_list": [nd],
                            "model_module": model_module})
    
    # observed variables
    model_kwargs['observed_vars_params'] = (model_kwargs['observed_vars'] + model_kwargs['observed_params'])
    # exclude bed variables from observed variables
    all_observed = model_kwargs['observed_vars_params']
    model_kwargs['all_observed'] = all_observed; params['all_observed'] = all_observed

    _modelrun_datasets = model_kwargs.get("data_path",None)
    if rank_world == 0 and not os.path.exists(_modelrun_datasets):
        # cretate the directory
        os.makedirs(_modelrun_datasets, exist_ok=True)

    # --- file_names
    _true_nurged   = f'{ _modelrun_datasets}/true_nurged_states.h5'
    _synthetic_obs = f'{ _modelrun_datasets}/synthetic_obs.h5'

    # --update model_kwargs with the file names
    model_kwargs.update({"true_nurged_file": _true_nurged, "synthetic_obs_file": _synthetic_obs})

     # fetch model nprocs
    model_nprocs = params.get("model_nprocs", 1)

    # set modeel_nprocs adaptively
    if model_kwargs.get('ICESEE_PERFORMANCE_TEST') or os.environ.get("ICESEE_PERFORMANCE_TEST"):
        total_cores = size_world * model_nprocs
    else:
        # Get total cores from SLURM environment (more reliable than os.cpu_count())
        try:
            total_cores = int(os.environ.get("SLURM_NTASKS", os.cpu_count()))
            slurm_nodes = int(os.environ.get("SLURM_JOB_NUM_NODES", 1))
        except ValueError:
            total_cores = os.cpu_count()  # Fallback if not in SLURM
            slurm_nodes = 1

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
                         'total_cores': total_cores,
                         'base_total_procs': base_total_procs,
                         'model_module': model_module,
                         'vec_inputs_old': model_kwargs.get('vec_inputs', params.get('vec_inputs', None)),
                        })

    #  --- Generate True and Nurged States -------------------------------------------------------------------
    # -- time generation of true state ----
    time_generation_true_and_wrong_state = time.time()
    # call the generate_true_wrong_state function
    model_kwargs = generate_true_wrong_state(**model_kwargs)
    # --- time generation of true state and nurged state ---
    time_generation_true_and_wrong_state = time.time() - time_generation_true_and_wrong_state


    # --- Generate the Synthetic ObservationsObservations ---------------------------------------------------
    # --- time generation of synthetic observations ---
    time_generation_synthetic_obs = time.time()
    # call the generate_synthetic_observations function
    model_kwargs =  generate_synthetic_observations(**model_kwargs)
    # --- time generation of synthetic observations ---
    time_generation_synthetic_obs = time.time() - time_generation_synthetic_obs

   # --- Initialize the ensemble ---------------------------------------------------
    Q_rho     = model_kwargs.get("Q_rho")
    len_scale = model_kwargs.get("length_scale")
    hdim  = params["nd"] // params["total_state_param_vars"]
    model_kwargs.update({"hdim": hdim, "Q_rho": Q_rho, "len_scale": len_scale})

        # --- get the process noise --->
    if params.get("use_random_fields", False):
        pos, gs_model, L_C = compute_Q_err_random_fields(hdim, params["total_state_param_vars"], params["sig_Q"], Q_rho, len_scale)
        model_kwargs.update({"pos": pos, "gs_model": gs_model, "L_C": L_C})
    
    # -- time ensemble initialization ---
    time_ensemble_initialization = time.time()
    # call the ensemble_initialization function
    if model_kwargs.get("initialize_ensemble", True):
        model_kwargs, ensemble_vec, time_init_noise_generation, \
        time_init_ensemble_mean_computation, time_init_file_writing, \
        shape_ens,ensemble_bg,  ensemble_vec_mean, ensemble_vec_full = ensemble_initialization(**model_kwargs)
    else:
        time_init_noise_generation = 0.0
        time_init_ensemble_mean_computation = 0.0
        time_init_file_writing = 0.0
        with h5py.File(_modelrun_datasets + '/icesee_ensemble_data.h5', 'r') as f:
            ensemble_vec = f['ensemble'][:, :, 0]
    # --- time ensemble initialization ---
    time_ensemble_initialization = time.time() - time_ensemble_initialization
    
    # --- get the ensemble size
    nd, Nens = ensemble_vec.shape
    module_nprocs = model_kwargs.get("model_nprocs", 1)
   
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
        # for i, sig in enumerate(params["sig_Q"]):
        #     start_idx = i *hdim
        #     end_idx = start_idx + hdim
            # Q_err[start_idx:end_idx,start_idx:end_idx] = np.eye(hdim) * sig ** 2

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
                x_points = np.linspace(0, Lx, nx+1)
                y_points = np.linspace(0, Ly, ny+1)
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
        model_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":params["total_state_param_vars"]})

        if (len(model_kwargs.get("scalar_inputs", [])) > 0) or (model_kwargs.get("var_nd", None) is not None):
                noise_1 = generate_enkf_field(**model_kwargs)
                ndim = 1 if len(model_kwargs.get("scalar_inputs", [])) > 0 else (model_kwargs["var_nd"][model_kwargs["scalar_inputs"][0]])
                model_kwargs.update({ "noise_dim": ndim})
                noise_2 = generate_enkf_field(**model_kwargs)
                # concatenate noise_1 and noise_2 
                noise = np.concatenate((noise_1, noise_2))[:-1]

        else:
            noise = generate_enkf_field(**model_kwargs)

    for k in range(model_kwargs.get("nt",params["nt"])):

        model_kwargs.update({"k": k, "km":km, "alpha": alpha, "rho": rho, "tau": tau, "dt": dt,"n": n, "noise": noise})
        model_kwargs.update({"generate_enkf_field": generate_enkf_field}) #save the function to generate the enkf field

        input_file = f"{_modelrun_datasets}/icesee_ensemble_data.h5"
        with h5py.File(input_file, "r") as f:
            ensemble_vec = f["ensemble"][:,:,k]

        # -------------- Forecast step
        ensemble_vec = EnKFclass.forecast_step(ensemble_vec, \
                                            model_module.forecast_step_single, \
                                             **model_kwargs)
       
        #  compute the ensemble mean
        with h5py.File(input_file, "a") as f:
            f["ensemble"][:,:,k+1] = ensemble_vec
            f["ensemble_mean"][:,k+1] = np.mean(ensemble_vec, axis=1)
       
        # -------------- Analysis step
        # generate Observation schedule
        obs_t, ind_m, m_obs = UtilsFunctions(params, ensemble_vec).generate_observation_schedule(**model_kwargs)
        params.update({"number_obs_instants": m_obs, 'obs_index': ind_m})
        model_kwargs.update({"m_obs": m_obs, "obs_t": obs_t, "obs_index": ind_m})
        model_kwargs.update({"params": params})
        if (km < m_obs) and (k == ind_m[km]):

            # read the ensemble mean from file
            with h5py.File(input_file, "r") as f:
                ensemble_vec_mean = f["ensemble_mean"][:,k+1]

            # Compute the model covariance
            diff = ensemble_vec - np.tile(ensemble_vec_mean.reshape(-1,1),Nens)
            if EnKF_flag or DEnKF_flag:
                Cov_model = 1/(Nens-1) * diff @ diff.T
            elif EnRSKF_flag or EnTKF_flag:
                Cov_model = 1/(Nens-1) * diff 

            # --- localization ---
            if params["localization_flag"]:
                if not adaptive_localization:
                    if hasattr(model_module, "localization_function") and callable(model_module.localization_function):
                        loc_matrix = model_module.localization_function(**model_kwargs)
                    else:
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

                # # inflate the top-left (smb h) and bottom-right (h smb) blocks of the covariance matrix 
                # state_block_size = num_state_vars*hdim
                # h_smb_block = Cov_model[:hdim,state_block_size:]
                # smb_h_block = Cov_model[state_block_size:,:hdim]

                # # apply the inflation factor
                # params["inflation_factor"] = 1.2
                # smb_h_block = UtilsFunctions(params, smb_h_block).inflate_ensemble(in_place=True)
                # h_smb_block = UtilsFunctions(params, h_smb_block).inflate_ensemble(in_place=True)

                # # update the covariance matrix
                # Cov_model[:hdim,state_block_size:] = h_smb_block
                # Cov_model[state_block_size:,:hdim] = smb_h_block

            # check if params["sig_obs"] is a scalar
            if isinstance(params["sig_obs"], (int, float)):
                params["sig_obs"] = np.ones(model_kwargs.get("nt",params["nt"])+1) * params["sig_obs"]

            # get synthetic observations if they exist
            with h5py.File(_synthetic_obs, 'r') as f:
                hu_obs  = f['hu_obs'][:]
      

            # print(f"R matrix: {R.shape},  cov shape: {Cov_model.shape}")
            # choose between user-defined functions or default functions
            # Obs = model_module.Obs_fun(hu_obs[:,km]) or UtilsFunctions(params, ensemble_vec).Obs_fun(hu_obs[:,km])
            # JObs = model_module.JObs_fun(nd) or UtilsFunctions(params, ensemble_vec).JObs_fun
            # analysis  = EnKF(Observation_vec=  UtilsFunctions(params, ensemble_vec).Obs_fun(hu_obs[:,km]), 
            #                 Cov_obs=R, \
            #                 Cov_model= Cov_model, \
            #                 Observation_function=UtilsFunctions(params, ensemble_vec).Obs_fun, \
            #                 Obs_Jacobian=UtilsFunctions(params, ensemble_vec).JObs_fun, \
            #                 parameters=  params,\
            #                 parallel_flag=   parallel_flag)

            # Create default functions object once
            utils = UtilsFunctions(params, ensemble_vec)

            if hasattr(model_module, "Cov_Obs_fun") and callable(model_module.Cov_Obs_fun):
                R = model_module.Cov_Obs_fun(sig_obs=params["sig_obs"][0],  nd=nd, kwargs=model_kwargs)
                # print(f"R.shape from model module: {R.shape}")
            else:
                R = np.eye(m_obs)
                for ii, sig in enumerate(params["sig_obs"]):
                    start_idx = ii *2
                    end_idx = start_idx +2
                    R[start_idx:end_idx,start_idx:end_idx] = np.eye(2) * sig ** 2

            # --- Select Obs_fun ---
            if hasattr(model_module, "Obs_fun") and callable(model_module.Obs_fun):
                Obs_fun = model_module.Obs_fun
            else:
                Obs_fun = utils.Obs_fun

            # --- Select JObs_fun ---
            if hasattr(model_module, "JObs_fun") and callable(model_module.JObs_fun):
                JObs_fun = model_module.JObs_fun
            else:
                JObs_fun = utils.JObs_fun

            # --- Evaluate observation vector for this time step ---
            Obs_vec = Obs_fun(hu_obs[:, km])

            # --- EnKF Analysis ---
            analysis = EnKF(
                Observation_vec      = Obs_vec,
                Cov_obs              = R,
                Cov_model            = Cov_model,
                Observation_function = Obs_fun,   # pass function handle
                Obs_Jacobian         = JObs_fun,  # pass function handle
                parameters           = params,
                parallel_flag        = parallel_flag
            )

            # Compute the analysis ensemble
            if EnKF_flag:
                ensemble_vec = analysis.EnKF_Analysis(ensemble_vec)
            elif DEnKF_flag:
                ensemble_vec = analysis.DEnKF_Analysis(ensemble_vec)
            elif EnRSKF_flag:
                ensemble_vec = analysis.EnRSKF_Analysis(ensemble_vec)
            elif EnTKF_flag:
                ensemble_vec = analysis.EnTKF_Analysis(ensemble_vec)
            else:
                raise ValueError("Filter type not supported")
            
            # model updates after analysis if any
            if hasattr(model_module, "post_analysis_update") and callable(model_module.post_analysis_update):
                model_kwargs, ensemble_vec = model_module.post_analysis_update(ensemble_vec, **model_kwargs)
            else:
                pass # no post analysis updates present

            # save the ensemble mean after analysis
            with h5py.File(input_file, "a") as f:
                f["ensemble_mean"][:,k+1] = np.mean(ensemble_vec, axis=1)
            
            # update the ensemble with observations instants
            km += 1

            # inflate the ensemble
            # ensemble_vec = UtilsFunctions(params, ensemble_vec).inflate_ensemble(in_place=True)
            # ensemble_vec = LocalizationInflationUtils(params, ensemble_vec).inflate_ensemble(in_place=True)
            # ensemble_vec = UtilsFunctions(params, ensemble_vec)._inflate_ensemble()
           
        # Save the ensemble
        with h5py.File(input_file, "a") as f:
            f["ensemble"][:,:,k+1] = ensemble_vec

        # update the progress bar
        if rank_world == 0:
            pbar.update(1)

    # close the progress bar
    if rank_world == 0:
        pbar.close()
    # comm_world.Barrier()

    # ====== load data to be written to file ======
    # print("[ICESEE] Saving data ...")
    # if params["even_distribution"]:
    #     save_all_data(
    #         enkf_params=model_kwargs['enkf_params'],
    #         nofilter=True,
    #         t=model_kwargs["t"], b_io=np.array([b_in,b_out]),
    #         Lxy=np.array([Lx,Ly]),nxy=np.array([nx,ny]),
    #         ensemble_true_state=ensemble_true_state,
    #         ensemble_nurged_state=ensemble_nurged_state, 
    #         obs_max_time=np.array([params["obs_max_time"]]),
    #         obs_index=model_kwargs["obs_index"],
    #         w=hu_obs,
    #         run_mode= np.array([params["execution_flag"]])
    #     )

    #     # --- Save final data ---
    #     save_all_data(
    #         enkf_params=model_kwargs['enkf_params'],
    #         ensemble_vec_full=ensemble_vec_full,
    #         ensemble_vec_mean=ensemble_vec_mean,
    #         ensemble_bg=ensemble_bg
    #     )
    # else:
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
    # --total elapsed time
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

    # Display elapsed time on rank 0
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


