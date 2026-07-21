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
from ICESEE.src.utils.tools import icesee_get_index, display_timing_default,display_timing_verbose, save_all_data, \
                                   load_bed_masks_from_h5
from ICESEE.src.run_model_da._error_generation import compute_Q_err_random_fields, \
                              compute_noise_random_fields, \
                              generate_enkf_field

from ICESEE.src.utils.inference_plugin import (
    reset_inference_plugin_state,
)

# ======================== Run model with EnKF ========================
def icesee_model_data_assimilation_partial_parallel(**model_kwargs): 
    """ General function to run any kind of model with the Ensemble Kalman Filter """

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
    resume_timestep   = 0

    # --- call the ICESEE mpi parallel manager ---
    if re.match(r"\AMPI_model\Z", parallel_flag, re.IGNORECASE):
        from mpi4py import MPI
        from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager
        from ICESEE.src.parallelization._mpi_analysis_functions import analysis_enkf_update, EnKF_X5, DEnKF_X5, \
                                             analysis_Denkf_update                                    
        from ICESEE.src.parallelization._mpi_forecast_functions import parallel_forecast_step_default_run, \
                                                                    parallel_forecast_step_squential_run, \
                                                                    parallel_forecast_step_even_distribution_run

        from ICESEE.src.parallelization._parallel_i_o import parallel_write_data_from_root_2D, \
                                            parallel_write_vector_from_root, parallel_write_full_ensemble_from_root, \
                                            parallel_write_ensemble_scattered, gather_and_broadcast_data_default_run
        from ICESEE.src.parallelization._mpi_generate_true_wrong_state import generate_true_wrong_state
        from ICESEE.src.parallelization._mpi_generate_synthetic_observations import generate_synthetic_observations
        from ICESEE.src.parallelization._mpi_ensemble_intialization import ensemble_initialization

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
                             "model_module": model_module,
                             'vec_inputs_old': model_kwargs.get('vec_inputs', params.get('vec_inputs', None)),})

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
        model_kwargs =  generate_synthetic_observations(**model_kwargs)
        # --- time generation of synthetic observations ---
        time_generation_synthetic_obs = MPI.Wtime() - time_generation_synthetic_obs
                    
        # --- Initialize the ensemble ---------------------------------------------------
        comm_world.Barrier()

        # end here if user choses to only generate the synthetic observations and the true and nurged states
        if model_kwargs.get("generate_true_wrong_state_only", False):
            if rank_world == 0:
                print("[ICESEE] Data generation complete. Exiting as per user request.")
            exit(0)

        Q_rho     = model_kwargs.get("Q_rho")
        len_scale = model_kwargs.get("length_scale")
        hdim  = params["nd"] // params["total_state_param_vars"]
        model_kwargs.update({"hdim": hdim, "Q_rho": Q_rho, "len_scale": len_scale})

        # The partial-parallel driver historically always restarted at zero.
        # For long ISSM reviewer experiments, allow an explicit restart from a
        # known-good ensemble slice.  This preserves the existing history and
        # rewinds the member exchange files before the next forecast.
        resume_timestep = int(model_kwargs.get("resume_from_timestep", 0) or 0)
        resume_ensemble = None
        if resume_timestep:
            if not params.get("default_run", False):
                raise ValueError(
                    "resume_from_timestep currently requires default_run"
                )
            if rank_world == 0:
                ensemble_history = os.path.join(
                    _modelrun_datasets, "icesee_ensemble_data.h5"
                )
                if not os.path.exists(ensemble_history):
                    raise FileNotFoundError(
                        f"Cannot resume: {ensemble_history} does not exist"
                    )
                with h5py.File(ensemble_history, "r") as f:
                    dset = f["ensemble"]
                    if not 0 < resume_timestep < dset.shape[2]:
                        raise ValueError(
                            "resume_from_timestep must select an existing "
                            f"positive slice; got {resume_timestep} for "
                            f"history length {dset.shape[2]}"
                        )
                    resume_ensemble = np.asarray(
                        dset[:, :, resume_timestep], dtype=float
                    )
                    if not np.all(np.isfinite(resume_ensemble)):
                        raise ValueError(
                            f"Ensemble slice {resume_timestep} is not finite"
                        )
                print(
                    "[ICESEE] Resuming partial run from stored ensemble "
                    f"timestep {resume_timestep}."
                )
            resume_shape = comm_world.bcast(
                None if rank_world else resume_ensemble.shape, root=0
            )
        else:
            resume_shape = None

         # --- get the process noise --->
        if params.get("use_random_fields", False):
            pos, gs_model, L_C = compute_Q_err_random_fields(hdim, params["total_state_param_vars"], params["sig_Q"], Q_rho, len_scale)
            model_kwargs.update({"pos": pos, "gs_model": gs_model, "L_C": L_C})
        
        # -- time ensemble initialization ---
        time_ensemble_initialization = MPI.Wtime()
        if resume_timestep:
            if rank_world == 0:
                ensemble_vec = resume_ensemble
                ensemble_vec_mean = np.mean(ensemble_vec, axis=1)
                ensemble_vec_full = ensemble_vec
                ensemble_bg = ensemble_vec.copy()

                # Rewind the per-member ISSM exchange files to the same known-
                # good slice.  The model MAT files retain boundary conditions;
                # run_model replaces all six coupled fields from these files.
                vec_inputs = model_kwargs.get(
                    "vec_inputs", params.get("vec_inputs", [])
                )
                if len(vec_inputs) * hdim != ensemble_vec.shape[0]:
                    raise ValueError(
                        "Cannot resume: vec_inputs do not match ensemble rows"
                    )
                # ISSM ensemble IDs in this driver are normally zero based
                # (ensemble_output_0 ... ensemble_output_39).  Retain support
                # for older one-based datasets without assuming either form.
                zero_based = os.path.exists(
                    os.path.join(_modelrun_datasets, "ensemble_output_0.h5")
                )
                member_id_offset = 0 if zero_based else 1
                for ens in range(ensemble_vec.shape[1]):
                    member_id = ens + member_id_offset
                    member_file = os.path.join(
                        _modelrun_datasets, f"ensemble_output_{member_id}.h5"
                    )
                    if not os.path.exists(member_file):
                        raise FileNotFoundError(
                            f"Cannot resume: {member_file} does not exist"
                        )
                    with h5py.File(member_file, "a") as f:
                        for block, key in enumerate(vec_inputs):
                            values = ensemble_vec[
                                block * hdim : (block + 1) * hdim, ens
                            ]
                            if key not in f:
                                raise KeyError(
                                    f"Cannot resume: /{key} missing from "
                                    f"{member_file}"
                                )
                            f[key][...] = values.reshape(f[key].shape)
                shape_ens = np.asarray(ensemble_vec.shape, dtype=np.int32)
            else:
                ensemble_vec = np.empty(resume_shape, dtype=float)
                ensemble_vec_mean = None
                ensemble_vec_full = None
                ensemble_bg = None
                shape_ens = None
            shape_ens = comm_world.bcast(shape_ens, root=0)
            comm_world.Barrier()
            time_init_noise_generation = 0.0
            time_init_ensemble_mean_computation = 0.0
            time_init_file_writing = 0.0
        elif model_kwargs.get("initialize_ensemble", True):
            # call the ensemble_initialization function
            model_kwargs, ensemble_vec, time_init_noise_generation, \
            time_init_ensemble_mean_computation, time_init_file_writing, \
            shape_ens,ensemble_bg,  ensemble_vec_mean, ensemble_vec_full = ensemble_initialization(**model_kwargs)
        else:
            time_init_noise_generation = 0.0
            time_init_ensemble_mean_computation = 0.0
            time_init_file_writing = 0.0
        # --- time ensemble initialization ---
        time_ensemble_initialization = MPI.Wtime() - time_ensemble_initialization
        
        # --- get the ensemble size
        nd = model_kwargs.get("nd", params["nd"])
        Nens = model_kwargs.get("Nens", params["Nens"])
        module_nprocs = model_kwargs.get("model_nprocs", 1)

        if params["even_distribution"]:
            ensemble_local = copy.deepcopy(ensemble_vec[:,start:stop])
                
        comm_world.Barrier()
        parallel_manager = None # debugging flag for now
        
    else:
        parallel_manager = None
        model_module = SupportedModels(model=model,verbose=params.get('verbose')).call_model()

        nd = params["nd"]
        Nens = params["Nens"]
        size_world = 1
        rank_world = 0
        sub_rank = 0
        color = 0

        _modelrun_datasets = model_kwargs.get("data_path",None)
        if rank_world == 0 and not os.path.exists(_modelrun_datasets):
            os.makedirs(_modelrun_datasets, exist_ok=True)

        _true_nurged   = f'{ _modelrun_datasets}/true_nurged_states.h5'
        _synthetic_obs = f'{ _modelrun_datasets}/synthetic_obs.h5'

        if params["even_distribution"] or (params["default_run"] and size_world <= params["Nens"]):
            if rank_world == 0:
                print("[ICESEE] Initializing the ensemble ...")
                model_kwargs.update({"statevec_ens":np.zeros([params["nd"], params["Nens"]])})
                
                vecs, indx_map, dim_per_proc = icesee_get_index(model_kwargs["statevec_ens"], **model_kwargs)
                ensemble_vec = np.zeros_like(model_kwargs["statevec_ens"])

                if model_kwargs["joint_estimation"] or params["localization_flag"]:
                    hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
                else:
                    hdim = ensemble_vec.shape[0] // params["num_state_vars"]
                state_block_size = hdim * params["num_state_vars"]

                for ens in range(params["Nens"]):
                    data = model_module.initialize_ensemble(ens,**model_kwargs)
                    for key, value in data.items():
                        ensemble_vec[indx_map[key],ens] = value

                    N_size = params["total_state_param_vars"] * hdim
                    model_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":params["total_state_param_vars"]})
                    noise = generate_enkf_field(**model_kwargs)
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
    if isinstance(params["sig_Q"], float):
        nd = hdim*params["total_state_param_vars"]
        Q_err = np.eye(nd) * params["sig_Q"] ** 2
    else:
        nd = hdim*params["total_state_param_vars"]
        Q_err = np.zeros((nd,nd))
    
    model_kwargs.update({"Q_err": Q_err})

    # --- Define filter flags
    EnKF_flag   = re.match(r"\AEnKF\Z", filter_type, re.IGNORECASE)
    DEnKF_flag  = re.match(r"\ADEnKF\Z", filter_type, re.IGNORECASE)
    EnRSKF_flag = re.match(r"\AEnRSKF\Z", filter_type, re.IGNORECASE)
    EnTKF_flag  = re.match(r"\AEnTKF\Z", filter_type, re.IGNORECASE)
    
    if rank_world == 0:
            nt = model_kwargs.get("nt", params["nt"])
            print(f"[ICESEE] Launching {model} with data assimilation using the {filter_type} filter across {size_world*(params['model_nprocs']+1)} MPI ranks.")
            pbar = tqdm(
                total=nt,
                initial=resume_timestep,
                desc=f"[ICESEE] Assimilation progress ({size_world*(params['model_nprocs']+1)} ranks)",
                position=0,
                leave=True,
                dynamic_ncols=True
        )

    # ==== Time loop =======================================================================================
    time_forecast_step = 0.0
    time_analysis_step = 0.0
    time_forecast_noise_generation = 0.0
    time_forecast_file_writing = 0.0
    time_analysis_file_writing = 0.0
    time_forecast_ensemble_mean_generation = 0.0
    time_analysis_ensemble_mean_generation = 0.0

    min_tau = 200
    max_tau = 500
    dt  = model_kwargs.get("dt",params["dt"])
    tau = max(max_tau,max(min_tau, dt))

    alpha = 1 - dt/tau
    if alpha <= 0 or alpha > 1:
        alpha = 0.5

    n = model_kwargs.get("nt",params["nt"])
    rho = np.sqrt((1/dt)*((1-alpha)**2)*(1/(n - (2*alpha) - (n*alpha**2) + (2*alpha**(n+1)))))
    params_analysis_0 = np.zeros((2, Nens))
    if resume_timestep:
        obs_index = np.asarray(model_kwargs.get("obs_index", []), dtype=int)
        km = int(np.searchsorted(obs_index, resume_timestep, side="left"))
    else:
        km = 0

    if params.get("use_random_fields", False):
        pos_obs, gs_model_obs, L_C_obs = compute_Q_err_random_fields(hdim, params["total_state_param_vars"], params["sig_obs"], Q_rho, len_scale)
    else:
        N_size = params["total_state_param_vars"] * hdim
        model_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":params["total_state_param_vars"]})
        noise = generate_enkf_field(**model_kwargs)
        model_kwargs.update({"noise": noise})

    # Reset private inference memory once at the beginning of a fresh run.
    reset_inference_plugin_state(model_kwargs)

    for k in range(resume_timestep, model_kwargs.get("nt",params["nt"])):

        model_kwargs.update({"k": k, "km":km, "alpha": alpha, "rho": rho, "tau": tau, "dt": dt,"n": n})
        model_kwargs.update({"generate_enkf_field": generate_enkf_field})

        if re.match(r"\AMPI_model\Z", parallel_flag, re.IGNORECASE):      
            _time_forecast_step = MPI.Wtime()

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
                                "state_block_size": state_block_size, "rng": None, "rank_seed": None,
                                "noise": noise,
                                }) 
            
            if not params.get("default_run", False):
                model_kwargs.update({"ensemble_vec": ensemble_vec,
                                "ensemble_vec_mean": ensemble_vec_mean,
                                "ensemble_vec_full": ensemble_vec_full,
                                "hdim": hdim,
                                "Nens": Nens,
                                "ensemble_local": ensemble_local if params.get("even_distribution", False) else None,
                                })                             
            
            if params["default_run"]:
                model_kwargs, ensemble_vec, shape_ens,ens_mean = parallel_forecast_step_default_run(**model_kwargs)
                time_forecast_step = model_kwargs.get("time_forecast_step", 0.0)
                time_forecast_noise_generation = model_kwargs.get("time_forecast_noise_generation", 0.0)
                time_forecast_file_writing = model_kwargs.get("time_forecast_file_writing", 0.0)
                time_forecast_ensemble_mean_generation = model_kwargs.get("time_forecast_ensemble_mean_generation", 0.0)

                time_forecast_step += MPI.Wtime() - _time_forecast_step

                # ===== Global analysis step =====
                if model_kwargs.get('global_analysis', True) or model_kwargs.get('local_analysis', False):
                   
                    obs_index = model_kwargs["obs_index"]
                    if (km < params["number_obs_instants"]) and (k == obs_index[km]):
                        _time_analysis_step = MPI.Wtime()
                        model_kwargs.update({"km": km})
                        inversion_enabled = bool(model_kwargs.get(
                            "inversion_enabled",
                            model_kwargs.get("inversion_flag", False),
                        ))
                        inversion_start_time = float(
                            model_kwargs.get("inversion_start_time", 0.0)
                        )
                        cycle_time = float(np.asarray(model_kwargs["t"])[k])
                        inversion_flag = (
                            inversion_enabled
                            and cycle_time + 1.0e-12 >= inversion_start_time
                        )
                        # Downstream analysis writers use inversion_flag to
                        # choose the five-block state/inversion workflow.
                        # Recompute it from inversion_enabled every cycle so
                        # an intentionally delayed first inversion does not
                        # disable all subsequent inversions.
                        model_kwargs["inversion_flag"] = inversion_flag
                        if rank_world == 0 and inversion_enabled and not inversion_flag:
                            print(
                                "[ICESEE] Deferring friction inversion at "
                                f"t={cycle_time:g} yr; configured start is "
                                f"{inversion_start_time:g} yr."
                            )
                        nd_old = model_kwargs.get("nd", nd)
                        model_kwargs.update({"nd_old": nd_old})

                        partitioned_io = model_kwargs.get("partitioned_io_flag", False)  # NEW

                        if inversion_flag:
                            # ---- EXISTING inversion-state-reduction logic, fully unchanged ----
                            if rank_world == 0:
                                with h5py.File(f'{_modelrun_datasets}/ensemble_before_analysis_step_{k+1:04d}.h5', 'w') as f:
                                    f.create_dataset("ensemble_before_analysis", data=ensemble_vec)
                                nd = ensemble_vec.shape[0]
                                Nens = ensemble_vec.shape[1]
                                model_kwargs.update({"nd_old": nd}); params.update({"nd_old": nd})
                                friction_idx = model_kwargs.get("friction_idx")
                                excluded_indices = [friction_idx]
                                vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
                                hdim = nd // (params["total_state_param_vars"] )
                                nd_new = hdim * (params["total_state_param_vars"] - 1)
                                ensemble_vec_reduced = np.zeros((nd_new, Nens))
                                for ii, key in enumerate(model_kwargs['vec_inputs']):
                                    if ii not in excluded_indices:
                                        start = ii * hdim
                                        end = start + hdim
                                        ensemble_vec_reduced[start:end, :] = ensemble_vec[indx_map[key], :]
                                ensemble_vec = copy.deepcopy(ensemble_vec_reduced)
                                model_kwargs.update({"nd": nd_new}); params["nd"] = nd_new

                                with h5py.File(_synthetic_obs, 'r') as f:
                                    _hu_obs  = f['hu_obs'][:]
                                    error_R = f['R'][:]
                                    bed_mask_map_static, bed_mask_map_cols, bed_snap_cols, obs_model_to_col = load_bed_masks_from_h5(f)
                                    Cov_obs = np.zeros(error_R.shape)

                                hu_obs = np.zeros((nd_new, _hu_obs.shape[1]))
                                for ii, key in enumerate(model_kwargs['vec_inputs']):
                                    if ii not in excluded_indices:
                                        start = ii * hdim
                                        end = start + hdim
                                        hu_obs[start:end, :] = _hu_obs[indx_map[key], :]

                                model_kwargs.update({
                                        "bed_mask_map_static": bed_mask_map_static,
                                        "bed_mask_map_cols": bed_mask_map_cols,
                                        "bed_snap_cols": bed_snap_cols,
                                        "obs_model_to_col": obs_model_to_col,
                                    })

                                model_kwargs['vec_inputs_new'] = [key for ii, key in enumerate(model_kwargs['vec_inputs']) if ii not in excluded_indices]
                                model_kwargs.update({'vec_inputs': model_kwargs['vec_inputs_new']})
                                vec_inputs = model_kwargs['vec_inputs_new']
                                params.update({'vec_inputs': vec_inputs})
                                nd = nd_new
                                params.update({"nd": nd, "total_state_param_vars": len(vec_inputs)})
                                model_kwargs.update({'excluded_indices': excluded_indices})
                            else:
                                nd_new = 0
                                vec_inputs = None
                                hu_obs = None
                                bed_mask_map = None
                                excluded_indices = None
                                bed_mask_map_static, bed_mask_map_cols, bed_snap_cols, obs_model_to_col = None, None, None, None

                            excluded_indices = comm_world.bcast(excluded_indices, root=0)
                            model_kwargs.update({'excluded_indices': excluded_indices})
                            hu_obs = comm_world.bcast(hu_obs, root=0)
                            bed_mask_map_static = comm_world.bcast(bed_mask_map_static, root=0)
                            bed_mask_map_cols = comm_world.bcast(bed_mask_map_cols, root=0)
                            bed_snap_cols = comm_world.bcast(bed_snap_cols, root=0)
                            obs_model_to_col = comm_world.bcast(obs_model_to_col, root=0)
                            model_kwargs.update({
                                        "bed_mask_map_static": bed_mask_map_static,
                                        "bed_mask_map_cols": bed_mask_map_cols,
                                        "bed_snap_cols": bed_snap_cols,
                                        "obs_model_to_col": obs_model_to_col,
                                    })
                            nd_new = comm_world.bcast(nd_new, root=0)
                            model_kwargs.update({'nd': nd_new}); params.update({'nd': nd_new})
                            # The inversion workflow removes the friction block
                            # before EnKF analysis. Keep the MPI Scatterv shape
                            # synchronized with that reduced five-block vector;
                            # retaining the full six-block shape over-allocates
                            # rows and corrupts hdim in the analysis writer.
                            shape_ens = np.array([nd_new, Nens], dtype=np.int32)
                            vec_inputs = comm_world.bcast(vec_inputs, root=0)
                            model_kwargs.update({'vec_inputs': vec_inputs})
                            params.update({'vec_inputs': vec_inputs})
                            params.update({"total_state_param_vars": len(vec_inputs)})
                            model_kwargs.update({"params": params})
                            
                        else:
                            with h5py.File(_synthetic_obs, 'r') as f:
                                hu_obs = f['hu_obs'][:]
                                error_R = f['R'][:]
                                bed_mask_map_static, bed_mask_map_cols, bed_snap_cols, obs_model_to_col = load_bed_masks_from_h5(f)
                                mask = ~np.isnan(hu_obs[:, km])
                                model_kwargs.update({
                                                "bed_mask_map_static": bed_mask_map_static,
                                                "bed_mask_map_cols": bed_mask_map_cols,
                                                "bed_snap_cols": bed_snap_cols,
                                                "obs_model_to_col": obs_model_to_col,
                                            })
                            
                        model_kwargs['observed_vars_params'] = (model_kwargs['observed_vars'] + model_kwargs['observed_params'])
                        all_observed = model_kwargs['observed_vars_params']
                        nd_new = len(all_observed)* hdim
                        model_kwargs.update({'all_observed': all_observed}); params.update({'all_observed': all_observed})

                        # ============================================== FIXED SECTION
                        if partitioned_io:
                            # Collective — every rank participates (mpio + Allreduce inside)
                            if EnKF_flag:
                                model_kwargs.update({"error_R": error_R})
                                X5, analysis_vec_ij = EnKF_X5(km, None, Nens, hu_obs, model_kwargs, UtilsFunctions)
                                time_analysis_mean_generation = 0.0
                                ens_mean = None
                            elif DEnKF_flag:
                                raise NotImplementedError(
                                    "partitioned_io_flag is not yet supported for DEnKF."
                                )
                            smb_scale = 1.0
                        else:
                            # ---- EXISTING behavior, unchanged ----
                            if rank_world == 0:
                                ndim = ensemble_vec.shape[0]//params["total_state_param_vars"]  
                                state_block_size = ndim*params["num_state_vars"]
                                model_kwargs.update({"error_R": error_R})

                                eta = 0.0
                                beta = np.ones(nd)

                                if EnKF_flag:
                                    X5, analysis_vec_ij = EnKF_X5(km, ensemble_vec, Nens, hu_obs, model_kwargs, UtilsFunctions)
                                    y_i = np.sum(X5, axis=1)
                                    time_analysis_mean_generation = MPI.Wtime()
                                    ens_mean = (1/Nens)*(ensemble_vec @ y_i.reshape(-1,1)).ravel()
                                    time_analysis_mean_generation = MPI.Wtime() - time_analysis_mean_generation

                                elif DEnKF_flag:
                                    X5, X5prime = DEnKF_X5(km, ensemble_vec, Cov_obs, Nens, model_kwargs, UtilsFunctions)
                                    analysis_vec_ij = None
                            else:
                                X5 = np.empty((Nens, Nens))
                                time_analysis_mean_generation = 0.0
                                analysis_vec_ij = None
                                smb_scale = 0.0
                                if DEnKF_flag:
                                    ens_mean = np.empty((nd, 1))

                            smb_scale = 1.0
                        # ============================================== END FIXED SECTION

                        # ---- Refresh hu_obs / bed masks for the next cycle — unconditional, unchanged ----
                        vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
                        with h5py.File(_synthetic_obs, 'r') as f:
                            hu_obs = f['hu_obs'][:]
                            bed_mask_map_static, bed_mask_map_cols, bed_snap_cols, obs_model_to_col = load_bed_masks_from_h5(f)

                        model_kwargs.update({
                                        "bed_mask_map_static": bed_mask_map_static,
                                        "bed_mask_map_cols": bed_mask_map_cols,
                                        "bed_snap_cols": bed_snap_cols,
                                        "obs_model_to_col": obs_model_to_col,
                                    })

                        # fetch the upper and lower bounds for every paramerter from observed data
                        ndim = hu_obs.shape[0]//params["total_state_param_vars"]
                        state_block_size = ndim*params["num_state_vars"]
                        bounds = []
                        for i, var in enumerate(model_kwargs["params_vec"]):
                            bound_idx = (params["num_state_vars"] + i) * ndim
                            bound_idx_end = bound_idx + ndim

                            param_slice = hu_obs[bound_idx:bound_idx_end, km]
                            param_slice = param_slice[~np.isnan(param_slice)]

                            if param_slice.size > 0:
                                param_min = np.min(param_slice)
                                param_max = np.max(param_slice)
                            else:
                                param_min = None
                                param_max = None

                            bounds.append(np.array([param_min, param_max]))

                        model_kwargs.update({"bounds": bounds})

                        # ============================================== FIXED SECTION
                        # call the analysis update function
                        if EnKF_flag:
                            time_analysis_mean_generation, time_analysis_file_writing = analysis_enkf_update(
                                k, ens_mean, ensemble_vec, shape_ens, X5,
                                time_analysis_mean_generation, time_analysis_file_writing,
                                analysis_vec_ij, UtilsFunctions, model_kwargs, smb_scale
                            )
                        elif DEnKF_flag:
                            model_kwargs.update({"DEnKF_flag": True})
                            analysis_Denkf_update(k, ens_mean, ensemble_vec, shape_ens, X5, UtilsFunctions, model_kwargs, smb_scale)
                        # ============================================== END FIXED SECTION
                    
                        km += 1
                        del hu_obs
                        gc.collect()
                        
                        time_analysis_step += MPI.Wtime() - _time_analysis_step

                        if inversion_flag:
                            nd = model_kwargs.get('nd_old')
                            model_kwargs.update({'nd': nd}); params.update({'nd': nd})
                            model_kwargs.update({'vec_inputs': model_kwargs.get('vec_inputs_old')})
                            params.update({'vec_inputs': model_kwargs.get('vec_inputs_old')})
                            params["total_state_param_vars"] = len(model_kwargs.get('vec_inputs_old'))
                            model_kwargs.update({"params": params})

                    else: 
                        partitioned_io = model_kwargs.get("partitioned_io_flag", False)  # NEW

                        if not partitioned_io:
                            # ---- EXISTING behavior, unchanged ----
                            _time_forecast_file_writing = MPI.Wtime()
                            parallel_write_full_ensemble_from_root(k+1,ens_mean, model_kwargs,ensemble_vec,comm_world)
                            _time_forecast_file_writing = MPI.Wtime() - _time_forecast_file_writing
                            time_forecast_file_writing += _time_forecast_file_writing
                            time_forecast_step = time_forecast_step + _time_forecast_file_writing
                            del ensemble_vec; gc.collect()
                        # else: partitioned path — the forecast step already wrote this
                        # timestep's ensemble directly via write_ensemble_member_direct;
                        # nothing further to do here.

                # ======= Local analyais step =======
                if model_kwargs.get('local_analysis', False):
                    pass

        if rank_world == 0:
            pbar.update(1)

    if rank_world == 0:
        pbar.close()

    save_all_data(
            enkf_params=model_kwargs['enkf_params'],
            nofilter=True,
            t=model_kwargs["t"], b_io=np.array([b_in,b_out]),
            Lxy=np.array([Lx,Ly]),nxy=np.array([nx,ny]),
            obs_max_time=np.array([params["obs_max_time"]]),
            obs_index=model_kwargs["obs_index"],
            run_mode= np.array([params["execution_flag"]])
        )

    global_end_time = MPI.Wtime()
    global_elapsed_time = global_end_time - global_start_time
    total_elapsed_time = comm_world.allreduce(global_elapsed_time, op=MPI.SUM)
    total_wall_time = comm_world.allreduce(global_elapsed_time, op=MPI.MAX)
    true_wrong_time = comm_world.allreduce(time_generation_true_and_wrong_state, op=MPI.MAX)
    ensemble_init_time = comm_world.allreduce(time_ensemble_initialization, op=MPI.MAX)
    forecast_step_time = comm_world.allreduce(time_forecast_step, op=MPI.MAX)
    forecast_noise_time = comm_world.allreduce(time_forecast_noise_generation, op=MPI.MAX)
    analysis_step_time = comm_world.allreduce(time_analysis_step, op=MPI.MAX)
    assimilation_time = ensemble_init_time + forecast_step_time + analysis_step_time
    forecast_file_time = comm_world.allreduce(time_forecast_file_writing, op=MPI.MAX)
    analysis_file_time = comm_world.allreduce(time_analysis_file_writing, op=MPI.MAX)
    init_file_time = comm_world.allreduce(time_init_file_writing, op=MPI.MAX)
    total_file_time = init_file_time + forecast_file_time + analysis_file_time
    time_analysis_ensemble_mean = comm_world.allreduce(time_analysis_ensemble_mean_generation, op=MPI.MAX)
    time_forecast_ensemble_mean= comm_world.allreduce(time_forecast_ensemble_mean_generation, op=MPI.MAX)
    time_init_ensemble_mean = comm_world.allreduce(time_init_ensemble_mean_computation, op=MPI.MAX)

    comm_world.Barrier()
    if rank_world == 0:
        verbose = model_kwargs.get("verbose", False)
        if True:
             display_timing_verbose(
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
            forecast_noise_time=forecast_noise_time, 
            time_init_ensemble_mean_computation=time_init_ensemble_mean,
            time_forecast_ensemble_mean_computation=time_forecast_ensemble_mean,
            time_analysis_ensemble_mean_computation=time_analysis_ensemble_mean,
            comm=comm_world
        )
        else:
            display_timing_default(total_elapsed_time, total_wall_time)
    else:
        None
