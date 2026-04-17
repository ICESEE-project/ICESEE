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
import shutil
import traceback
import numpy as np
from tqdm import tqdm 
import bigmpi4py as BM # BigMPI for large data transfer and communication
from mpi4py import MPI
import json, glob, tempfile

CKPT_DIRNAME = "_checkpoints"
CKPT_BASENAME = "icesee_ckpt.json"
FNAME_PATTERN = r'icesee_enkf_ens_(\d+)\.h5$'  # matches ..._0000.h5, ..._12.h5, etc.

# ==== ICESEE utility imports ========================================
from ICESEE.src.utils import tools, utils                                     # utility functions for the model 
from ICESEE.src.utils.utils import UtilsFunctions
from ICESEE.applications.supported_models import SupportedModels              # supported models for data assimilation routine
from ICESEE.src.utils.tools import icesee_get_index, display_timing_default,display_timing_verbose, \
                                    save_all_data, finalize_stack, _extract_time_from_name, _sorted_step_files,\
                                    _last_completed_step, _ckpt_path, _atomic_write_json, save_checkpoint, load_checkpoint, \
                                    compute_km_from_tobserve, step_already_done, reseed_for_step, icesee_fingerprint, h5_has_dataset_with_shape, \
                                    h5_attr_equals, mark_h5_with_fingerprint, env_flag
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
    nx, ny            = model_kwargs.get("nx",1), model_kwargs.get("ny",1)
    b_in, b_out       = model_kwargs.get("b_in",0.0), model_kwargs.get("b_out",0.0) 
    data_path         = model_kwargs.get("data_path","_modelrun_datasets")      # data path
    restart_enabled   = model_kwargs.get("restart_enabled", True)   # turn on/off restart
    force_fresh_start = model_kwargs.get("force_fresh_start", False) # ignore old files/ckpt
    checkpoint_every  = model_kwargs.get("checkpoint_every", 1)     # write ckpt every N steps
    base_seed         = model_kwargs.get("base_seed", 0)             # for reproducible reseed
    nd               = model_kwargs.get("nd",params["nd"])      # model dimension
    Nens             = model_kwargs.get("Nens",params["Nens"])  # number of ensemble members
    nt               = model_kwargs.get("nt",params["nt"])      # number of time steps


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
    
    model_kwargs['observed_vars_params'] = (model_kwargs['observed_vars'] + model_kwargs['observed_params'])
    all_observed = model_kwargs['observed_vars_params']

    model_kwargs.update({'all_observed': all_observed}); params.update({'all_observed': all_observed})
    model_kwargs.update({'params': params})

    # pack the global communicator and the subcommunicator
    model_kwargs.update({"comm_world": comm_world, "subcomm": subcomm})

    # --- check if the modelrun dataset directory is present ---
    _modelrun_datasets = model_kwargs.get("data_path",None)
    if rank_world == 0 and not os.path.exists(_modelrun_datasets):
        # cretate the directory
        os.makedirs(_modelrun_datasets, exist_ok=True)
        
    # Ensure checkpoint dir exists
    if rank_world == 0:
        os.makedirs(os.path.join(_modelrun_datasets, CKPT_DIRNAME), exist_ok=True)
    comm_world.Barrier()

    ckpt = None
    if restart_enabled and not force_fresh_start:
        if rank_world == 0:
            ckpt = load_checkpoint(_modelrun_datasets)
        ckpt = comm_world.bcast(ckpt, root=0)

    comm_world.Barrier()
    # --- file_names
    _true_nurged   = f'{ _modelrun_datasets}/true_nurged_states.h5'
    _synthetic_obs = f'{ _modelrun_datasets}/synthetic_obs.h5'

    # --update model_kwargs with the file names
    model_kwargs.update({"true_nurged_file": _true_nurged, "synthetic_obs_file": _synthetic_obs})

    # Build a reproducibility fingerprint from current config
    fp = icesee_fingerprint({
        "model_name": model,
        "nd": nd,
        "nt": nt,
        "Nens": Nens,
        "base_seed": base_seed,
    })

    # Should we reuse prior artifacts?
    reuse_allowed = restart_enabled and not force_fresh_start

    # --- initialize seed for reproducibility ---
    ParallelManager().initialize_seed(comm_world, base_seed=base_seed)

    # --- intialize EnKF I/O handler class ---
    time_file_io_initialization = MPI.Wtime()
    batch_size = model_kwargs.get("batch_size", nt if nt <= 100 else max(1, (nt + 9) // 10))
    serial_file_creation = model_kwargs.get("serial_file_creation",True)
    h5_file_compression = model_kwargs.get("h5_file_compression",None)
    h5_file_compression_level = model_kwargs.get("h5_file_compression_level",4)
    h5_file_chunk_size = model_kwargs.get("h5_file_chunk_size",1000)
    enkf_parallel_io = EnKF_fully_parallel_IO('icesee_enkf', nd, Nens, nt, subcomm, comm_world, \
                                             params, serial_file_creation, base_path=_modelrun_datasets, \
                                             batch_size=batch_size, h5_file_compression=h5_file_compression, \
                                             h5_file_compression_level=h5_file_compression_level, \
                                             h5_file_chunk_size=h5_file_chunk_size)
    # Update model_kwargs with the EnKF I/O handler
    model_kwargs.update({"enkf_parallel_io": enkf_parallel_io})
    time_file_io_initialization = MPI.Wtime() - time_file_io_initialization

    try:

        # fetch model nprocs
        model_nprocs = params.get("model_nprocs", 1)

        # set modeel_nprocs adaptively
        # total_cores = os.cpu_count()
        if model_kwargs.get('ICESEE_PERFORMANCE_TEST') or env_flag("ICESEE_PERFORMANCE_TEST", default=False):
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
            effective_model_nprocs = int(max(1, np.floor(effective_model_nprocs * scale_factor)))

        # update model_kwargs with the effective model_nprocs
        model_kwargs.update({'model_nprocs': effective_model_nprocs,
                                "total_cores": total_cores,
                                "base_total_procs": base_total_procs,
                            })

        # --- Generate True and Nurged States -------------------------------------------------------------------
        # ---- Generate True and Nurged States (skip on restart if valid) ----
        time_generation_true_and_wrong_state = MPI.Wtime()

        need_true = True
        true_reason = None

        if reuse_allowed and os.path.exists(_true_nurged):
            if h5_has_dataset_with_shape(_true_nurged, "true_state", (nd, nt+1)) \
            and h5_has_dataset_with_shape(_true_nurged, "nurged_state", (nd, nt+1)):
                if h5_attr_equals(_true_nurged, "icesee_fingerprint", fp) \
                or model_kwargs.get("allow_reuse_without_fingerprint", False):
                    need_true = False
                else:
                    true_reason = "fingerprint mismatch/absent"
            else:
                true_reason = "dataset(s) missing or shape mismatch"
        else:
            true_reason = "file missing/fresh start"

        if not need_true and rank_world == 0:
            print("[ICESEE][RESTART] Using existing true/nurged states.")
        elif need_true and rank_world == 0:
            print(f"[ICESEE][RESTART] Regenerating true/nurged states ({true_reason}).")

        if need_true:
            model_kwargs = generate_true_wrong_state(**model_kwargs)
            if rank_world == 0 and os.path.exists(_true_nurged):
                mark_h5_with_fingerprint(_true_nurged, value=fp, extra={
                    "dataset_name_true": "true_state",
                    "dataset_name_nurged": "nurged_state",
                })

        time_generation_true_and_wrong_state = MPI.Wtime() - time_generation_true_and_wrong_state
        comm_world.Barrier()

        # --- Generate the Synthetic ObservationsObservations ---------------------------------------------------
        # ---- Synthetic observations (skip if present & matching) ----
        time_generation_synthetic_obs = MPI.Wtime()

        # tobserve = None
        # m_obs = None

        syn_ok = False
        if reuse_allowed and os.path.exists(_synthetic_obs):
            syn_ok = True
            # try:
            #     with h5py.File(_synthetic_obs, "r") as f:
            #         # must exist
            #         if "tobserve" in f and "m_obs" in f:
            #             tobserve = f["tobserve"][...]           # shape (m_obs,)
            #             m_obs = int(f["m_obs"][0])
            #             # basic sanity: increasing times, within [1..nt]
            #             if tobserve.ndim == 1 and m_obs == len(tobserve) and m_obs >= 0:
            #                 if np.all(np.diff(tobserve) >= 0) and int(tobserve[-1]) <= int(nt):
            #                     # fingerprint (strict) or allow fallback if absent
            #                     if h5_attr_equals(_synthetic_obs, "icesee_fingerprint", fp) \
            #                     or model_kwargs.get("allow_reuse_without_fingerprint", False):
            #                         syn_ok = True
            # except Exception:
            #     syn_ok = False

        if model_kwargs.get("generate_synthetic_obs", True) and not syn_ok:
            synthetic_obs_zarr_path = f"{_modelrun_datasets}/synthetic_observations.zarr"
            error_R_zarr_path = f"{_modelrun_datasets}/error_R.zarr"
            model_kwargs.update({'synthetic_obs_zarr_path': synthetic_obs_zarr_path, 'error_R_zarr_path': error_R_zarr_path})
            tobserve, m_obs = enkf_parallel_io._create_synthetic_observations(**model_kwargs)
            # if rank_world == 0:
            #     with h5py.File(_synthetic_obs, "w") as f:
            #         f.create_dataset("tobserve", data=np.asarray(tobserve, dtype=np.int64))
            #         f.create_dataset("m_obs", data=np.asarray([m_obs], dtype=np.int64))
            #         f.attrs["icesee_fingerprint"] = fp
        else:
            if rank_world == 0:
                print("[ICESEE][RESTART] Using existing synthetic observations.")
            _, tobserve, m_obs = enkf_parallel_io.generate_observation_schedule(**model_kwargs)

        model_kwargs.update({"tobserve": tobserve, "m_obs": m_obs})
        time_generation_synthetic_obs = MPI.Wtime() - time_generation_synthetic_obs
        comm_world.Barrier()

        # ----- Decide restart step (supports explicit override) -----
        k_start = 0

        # User override (highest priority)
        k_start_override = model_kwargs.get("k_start_override", None)
        if k_start_override is not None and not force_fresh_start:
            if not (0 <= int(k_start_override) <= int(nt)):
                raise ValueError(f"k_start_override={k_start_override} out of range [0, {nt}]")
            k_start = int(k_start_override)
        else:
            # Normal restart logic
            if restart_enabled and not force_fresh_start:
                if ckpt is not None and "last_done_k" in ckpt:
                    k_start = int(ckpt["last_done_k"] + 1)
                else:
                    last_k = _last_completed_step(_modelrun_datasets)
                    if last_k is not None:
                        k_start = int(last_k + 1)

        # Clamp
        k_start = min(max(0, k_start), nt)

        # (Optional) safety: ensure files before k_start exist contiguously
        if model_kwargs.get("enforce_contiguous_history", True) and k_start > 0:
            missing = []
            for kk in range(k_start):
                if not step_already_done(_modelrun_datasets, kk):
                    missing.append(kk)
            if missing:
                raise RuntimeError(
                    f"Cannot start at k={k_start}: missing completed steps {missing[:10]}{' ...' if len(missing)>10 else ''}. "
                    "Either lower k_start_override or disable enforce_contiguous_history."
                )

        # (Optional) destructive but safe: truncate files AFTER k_start-1
        if model_kwargs.get("truncate_after_k_start", False):
            # delete any step files >= k_start
            for fname in _sorted_step_files(_modelrun_datasets):
                kk = _extract_time_from_name(fname)
                if kk >= k_start:
                    try: os.remove(fname)
                    except Exception: pass
            # clear checkpoint so we honor the override on subsequent restarts
            if rank_world == 0:
                ckpt_path = _ckpt_path(_modelrun_datasets)
                if os.path.exists(ckpt_path):
                    try: os.remove(ckpt_path)
                    except Exception: pass
        comm_world.Barrier()

        # Recompute km consistent with your (k+1 == tobserve[km]) condition
        # print(f"\n[ICESEE] Starting at k={k_start} (nt={nt}) on rank {rank_world}.\n")
        km = compute_km_from_tobserve(np.asarray(tobserve), k_start, m_obs)
        model_kwargs.update({"km": km})

        # If we’re resuming, let the user know
        if rank_world == 0:
            if k_start > 0:
                print(f"[ICESEE][RESTART] Resuming at k={k_start} (nt={nt}); km={km}")
            else:
                print(f"[ICESEE] Fresh start at k=0 (nt={nt})")

        comm_world.Barrier()
        #  --- generate the H file
        # ---- H matrix (skip if present & matching) ----
        H_matrix_zarr_path = f"{_modelrun_datasets}/H_matrix.zarr"
        # if rank_world == 0:
        #     need_H = True
        #     meta_h5 = f"{_modelrun_datasets}/H_matrix_meta.h5"
        #     if reuse_allowed and os.path.isdir(H_matrix_zarr_path) and tools.h5_attr_equals(meta_h5, "icesee_fingerprint", fp):
        #         need_H = False
        #         print("[ICESEE][RESTART] Using existing H matrix.")
        #     if need_H:
        #         print("[ICESEE] Generating H matrix and saving to Zarr...")
        #         model_kwargs.update({'H_matrix_zarr_path': H_matrix_zarr_path})
        #         enkf_parallel_io.H_matrix(**model_kwargs)
        #         # record a tiny meta tag
        #         with h5py.File(meta_h5, "w") as f:
        #             f.attrs["icesee_fingerprint"] = fp
            
        model_kwargs.update({'H_matrix_zarr_path': H_matrix_zarr_path})
        enkf_parallel_io.H_matrix(**model_kwargs)
        # comm_world.Barrier()
                    
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
        time_ensemble_initialization = MPI.Wtime()

        init_ok = reuse_allowed and tools.h5_has_dataset_with_shape(
            os.path.join(_modelrun_datasets, "icesee_enkf_ens_0000.h5"),
            "states", (nd, Nens)
        )
        if init_ok:
            if rank_world == 0:
                print("[ICESEE][RESTART] Skipping ensemble initialization (found step 0).")
            
            # load only dictionary essentials
            if size_world <= params["Nens"]:
                model_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})
                dim_list = comm_world.allgather(model_kwargs.get("nd", params["nd"]))
                model_kwargs.update({"global_shape": model_kwargs.get("nd", params["nd"]), "dim_list": dim_list})
            else:
                model_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})
                model_kwargs.update({'ens_id': color}) # Nens = color
                # gather all the vector dimensions from all processors
                dim_list = subcomm.allgather(model_kwargs.get("nd", params["nd"]))
                global_shape = sum(dim_list)
                model_kwargs.update({"global_shape": global_shape, "dim_list": dim_list})
            
            time_init_file_writing = 0.0
            time_init_noise_generation = 0.0
            time_init_ensemble_mean_computation = 0.0
        else:
            # call the ensemble_initialization function
            if model_kwargs.get("initialize_ensemble", True):
                model_kwargs, ensemble_vec, time_init_noise_generation, \
                time_init_ensemble_mean_computation, time_init_file_writing, \
                shape_ens,ensemble_bg,  ensemble_vec_mean, ensemble_vec_full = ensemble_initialization_full_parallel_run(**model_kwargs)
            else:
                # If ensemble initialization is disabled, set default values
                ensemble_vec = None
                time_init_noise_generation = 0.0
                time_init_ensemble_mean_computation = 0.0
                time_init_file_writing = 0.0
                shape_ens = (params["nd"], params["nd"])
                ensemble_bg = np.zeros(shape_ens)
                ensemble_vec_mean = np.zeros((params["nd"], 1))
                ensemble_vec_full = np.zeros(shape_ens)

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
            print(f"[ICESEE] Launching {model} with data assimilation using the {filter_type} filter across {size_world*(params['model_nprocs']+1)} MPI ranks.")
            pbar = tqdm(
                total=nt,
                desc=f"[ICESEE] Assimilation progress ({size_world*(params['model_nprocs']+1)} ranks)",
                position=0,
                leave=True,
                dynamic_ncols=True,
                initial=k_start   # <-- start from resumed step
            )

        # synchronize all processes before starting the time loop
        comm_world.Barrier()

        # ==== Time loop =======================================================================================
        # --- timing intializations
        time_forecast_step = 0.0
        time_analysis_step = 0.0
        time_forecast_noise_generation = 0.0
        time_forecast_file_writing = 0.0
        time_analysis_file_writing = 0.0
        time_forecast_ensemble_mean_generation = 0.0
        time_analysis_ensemble_mean_generation = 0.0

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
        # km = 0

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
            noise = generate_enkf_field(**model_kwargs)
        
        # synchronize all processes before starting the time loop
        comm_world.Barrier()

        # for k in range(model_kwargs.get("nt",params["nt"])):
        for k in range(k_start, model_kwargs.get("nt",params["nt"])): # resume from k_start

            # Idempotency guard: if output for this k exists, skip safely
            if (not force_fresh_start) and step_already_done(_modelrun_datasets, k):
                if (km < m_obs) and (k+1 == model_kwargs.get("tobserve")[km]):
                    km += 1
                    model_kwargs.update({"km": km})
                if rank_world == 0:
                    pbar.update(1)
                continue

            # Deterministic reseed per step (optional but recommended)
            rank_seed = reseed_for_step(base_seed, rank_world, k)
            model_kwargs.update({"rank_seed": rank_seed})

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
                                    "km": km,
                                    "k": k,
                                    "len_scale": len_scale,
                                    "model_module": model_module,
                                    "time_forecast_step": time_forecast_step,
                                    "time_analysis_step": time_analysis_step,
                                    "time_forecast_noise_generation": time_forecast_noise_generation,
                                    "time_forecast_file_writing": time_forecast_file_writing,
                                    "time_analysis_file_writing": time_analysis_file_writing,
                                    "time_forecast_ensemble_mean_generation": time_forecast_ensemble_mean_generation,
                                    "time_analysis_ensemble_mean_generation": time_analysis_ensemble_mean_generation,
                                    "state_block_size": state_block_size, "noise": noise, "rng": None, "rank_seed": None,})
                
                if params["default_run"]:
                    # call the parallel_forecast_step_default_run function
                    model_kwargs = parallel_forecast_step_default_full_parallel_run(**model_kwargs)
                    time_forecast_step = model_kwargs.get("time_forecast_step", 0.0)
                    time_forecast_noise_generation = model_kwargs.get("time_forecast_noise_generation", 0.0)
                    time_forecast_file_writing = model_kwargs.get("time_forecast_file_writing", 0.0)
                    time_forecast_ensemble_mean_generation = model_kwargs.get("time_forecast_ensemble_mean_generation", 0.0)
                
                    comm_world.Barrier()
                    # print(f"[ICESEE] Rank {rank_world}, completed time step {k+1}/{params['nt']} with forecast time {time_forecast_step:.2f}s.")
                    # --- end time forecast step
                    time_forecast_step += MPI.Wtime() - _time_forecast_step

                    # ===== Global analysis step =====
                    if model_kwargs.get('global_analysis', True) or model_kwargs.get('local_analysis', False):
                        
                        tobserve = model_kwargs.get("tobserve")
                        m_obs = model_kwargs.get("m_obs", params["number_obs_instants"])
                        # if (km < m_obs) and (k+1 == tobserve[km]):
                        # if (km < m_obs) and (k == tobserve[km]):
                        obs_index = model_kwargs["obs_index"]
                        if (km < params["number_obs_instants"]) and (k == obs_index[km]):
                            # -- time global analysis step ---
                            _time_analysis_step = MPI.Wtime()
                            model_kwargs.update({'km': km, 'k': k})

                            inversion_flag = model_kwargs.get("inversion_flag", False)
                            nd_old = model_kwargs.get("nd", nd)
                            model_kwargs.update({"nd_old": nd_old})
        
                            # call the analysis update function
                            if EnKF_flag:
                                model_kwargs = enkf_parallel_io.compute_analysis_update(**model_kwargs)
                                time_analysis_ensemble_mean_generation = model_kwargs.get("time_analysis_ensemble_mean_generation", 0.0)
                                time_analysis_file_writing = model_kwargs.get("time_analysis_file_writing", 0.0)
                            
                            # update the observation index
                            km += 1
        #                    
                            # --- end time analysis step ---
                            time_analysis_step += MPI.Wtime() - _time_analysis_step
                    
                    # Step k fully completer; checkpoint if needed
                    if restart_enabled and (k % checkpoint_every == 0 or k == nt - 1):
                        # Build a minimal state; only rank 0 writes
                        if rank_world == 0:
                            ck = {
                                "last_done_k": k,
                                "km": int(km),
                                "nt": int(nt),
                                "nd": int(nd),
                                "nens": int(Nens),
                                "dataset_dir": os.path.abspath(_modelrun_datasets),
                                "timestamp": time.time(),
                                "base_seed": int(base_seed),
                            }
                            try:
                                save_checkpoint(_modelrun_datasets, **ck)
                            except Exception as e:
                                print(f"[ICESEE][WARN] Failed to save checkpoint at k={k}: {e}")

            # update the progress bar
            if rank_world == 0:
                pbar.update(1)

        # close the progress bar
        if rank_world == 0:
            pbar.close()
        comm_world.Barrier()
        time_file_io_closing = MPI.Wtime()
        enkf_parallel_io.close()
        # --- Build the Virtual Dataset view for the entire run ---
        # if rank_world == 0:
        #     print("[ICESEE] Building unified Virtual Dataset...")
        # enkf_parallel_io.create_virtual_dataset()
        time_file_io_initialization += MPI.Wtime() - time_file_io_closing

        # comm_world.Barrier()  
        # # ====== load data to be written to file ======
        if rank_world == 0:
            print("[ICESEE] Saving data ...")
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

        # ───────── Collective finalize (safe across ranks) ─────────
        # comm_world.Barrier()  # ensure all finished compute before finalize
        # print(f"[ICESEE] Rank {rank_world} entering finalize.")
        t0_final = MPI.Wtime()
        finalize_ok = True
        finalize_err = ""
        
        if rank_world == 0:
            try:
                # --- create the ensemble dataset ---
                if model_kwargs.get("create_ensemble_dataset", True):
                    print("[ICESEE] Creating ensemble dataset...")
                    # Option A: no-copy, instant
                    # out_vds = finalize_stack(_modelrun_datasets, mode="vds", dset_name="states")
                    # print("VDS ready:", out_vds)

                    # Option B: portable single file
                    out_h5 = finalize_stack(_modelrun_datasets, mode="h5", dset_name="states",
                                            allow_missing=False, compression="gzip", compression_opts=4)
                    print("Materialized file:", out_h5)
                # --- remove all .zarr files ---
                cleanup_intermediates = model_kwargs.get("cleanup_intermediates", True)
                if cleanup_intermediates:
                    for item in os.listdir(_modelrun_datasets):
                        if item.endswith(".zarr"):
                            item_path = os.path.join(_modelrun_datasets, item)
                            if os.path.isdir(item_path):
                                shutil.rmtree(item_path, ignore_errors=True)
                                print(f"[ICESEE] Removed {item_path}")

            except Exception as e:
                finalize_ok = False
                finalize_err = f"{type(e).__name__}: {e}"

        # Broadcast finalize status to all ranks so nobody hangs at a barrier
        finalize_ok = comm_world.bcast(finalize_ok, root=0)
        finalize_err = comm_world.bcast(finalize_err, root=0)

        if not finalize_ok:
            # Raise collectively so all ranks exit the same way
            raise RuntimeError(f"[ICESEE][FINALIZE] Root finalize failed: {finalize_err}")
        # print(f"[ICESEE] Rank {rank_world} finalize successful.")
        comm_world.Barrier()  # all ranks leave finalize together
        # print(f"[ICESEE] Rank {rank_world} passed finalize barrier.\n")
        time_final_file_writing = MPI.Wtime() - t0_final
        # ──────── end collective finalize ────────
        # ─────────────────────────────────────────────────────────────
        #  End Timer and Aggregate Elapsed Time Across Processors
        # ─────────────────────────────────────────────────────────────
        # ── Collective timing reductions (exception-safe) ──
        timing_ok = True
        timing_err = ""

        try:
            # --total elapsed time
            global_end_time = MPI.Wtime()
            global_elapsed_time = global_end_time - global_start_time

            # Reduce elapsed time across all processors (sum across ranks)
            # print(f"\n[ICESEE] Rank {rank_world} starting elapsed time reduction.")
            total_elapsed_time = comm_world.allreduce(global_elapsed_time, op=MPI.SUM)
            # print(f"[ICESEE] Rank {rank_world} finished elapsed time reduction.\n")

            # print(f"\n[ICESEE] Rank {rank_world} starting wall time reduction.")  
            total_wall_time = comm_world.allreduce(global_elapsed_time, op=MPI.MAX)
            # print(f"[ICESEE] Rank {rank_world} finished wall time reduction.\n")

            # -- timing true and wrong state generation
            # print(f"\n[ICESEE] Rank {rank_world} starting true/wrong state time reduction.")
            true_wrong_time = comm_world.allreduce(time_generation_true_and_wrong_state, op=MPI.MAX)
            # print(f"[ICESEE] Rank {rank_world} finished true/wrong state time reduction.\n")

            # -- timing ensemble initialization
            # print(f"\n[ICESEE] Rank {rank_world} starting ensemble initialization time reduction.") 
            ensemble_init_time = comm_world.allreduce(time_ensemble_initialization, op=MPI.MAX)
            # print(f"[ICESEE] Rank {rank_world} finished ensemble initialization time reduction.\n")

            # -- timing forecast step
            # print(f"\n[ICESEE] Rank {rank_world} starting forecast step time reduction.")
            forecast_step_time = comm_world.allreduce(time_forecast_step, op=MPI.MAX)
            # print(f"[ICESEE] Rank {rank_world} finished forecast step time reduction.\n")

            # -- timing forecast noise generation
            # print(f"\n[ICESEE] Rank {rank_world} starting forecast noise generation time reduction.")
            forecast_noise_time = comm_world.allreduce(time_forecast_noise_generation, op=MPI.MAX)
            # print(f"[ICESEE] Rank {rank_world} finished forecast noise generation time reduction.\n")

            # -- timing analysis step
            # print(f"\n[ICESEE] Rank {rank_world} starting analysis step time reduction.")
            analysis_step_time = comm_world.allreduce(time_analysis_step, op=MPI.MAX)
            # print(f"[ICESEE] Rank {rank_world} finished analysis step time reduction.\n")

            # -- total assimilation time = ensemble init + forecast step + analysis step
            assimilation_time = ensemble_init_time + forecast_step_time + analysis_step_time

            # --- time forecast file writing ---
            # print(f"\n[ICESEE] Rank {rank_world} starting forecast file writing time reduction.")
            forecast_file_time = comm_world.allreduce(time_forecast_file_writing, op=MPI.MAX)
            # print(f"[ICESEE] Rank {rank_world} finished forecast file writing time reduction.\n")

            # --- time analysis file writing ---
            # print(f"\n[ICESEE] Rank {rank_world} starting analysis file writing time reduction.")
            analysis_file_time = comm_world.allreduce(time_analysis_file_writing, op=MPI.MAX)
            # print(f"[ICESEE] Rank {rank_world} finished analysis file writing time reduction.\n")

            # total file writing time initialization file writing + forecast file writing + analysis file writing
            # print(f"\n[ICESEE] Rank {rank_world} starting initialization file writing time reduction.")
            init_file_time = comm_world.allreduce(time_init_file_writing, op=MPI.MAX)      
            # print(f"[ICESEE] Rank {rank_world} finished initialization file writing time reduction.\n")
            total_file_time = init_file_time + forecast_file_time + analysis_file_time + time_final_file_writing + time_file_io_initialization

            time_analysis_ensemble_mean = comm_world.allreduce(time_analysis_ensemble_mean_generation, op=MPI.MAX)
            time_forecast_ensemble_mean= comm_world.allreduce(time_forecast_ensemble_mean_generation, op=MPI.MAX)
            time_init_ensemble_mean = comm_world.allreduce(time_init_ensemble_mean_computation, op=MPI.MAX)
            # print(f"[ICESEE] Rank {rank_world}
        except Exception as e:
            timing_ok = False
            timing_err = f"{type(e).__name__}: {e}"
            tb_str = "".join(traceback.format_exception(*sys.exc_info()))
            print(f"Traceback details:\n{tb_str}")

        # Broadcast timing status to all ranks so nobody hangs at a barrier
        timing_ok = comm_world.bcast(timing_ok, root=0)
        timing_err = comm_world.bcast(timing_err, root=0)
        if not timing_ok:
            # Raise collectively so all ranks exit the same way
            raise RuntimeError(f"[ICESEE][TIMING] Collective timing reduction failed: {timing_err}")
        # ── end collective timing reductions ──

        # Display elapsed time on rank 0
        # print(f"[ICESEE] Rank {rank_world} finished in {global_elapsed_time:.2f}s (wall {global_end_time - global_start_time:.2f}s).")
        comm_world.Barrier()
        # print(f"[ICESEE] Rank {rank_world} passed timing barrier.")
        if rank_world == 0:
            verbose = model_kwargs.get("verbose", False)
            # if verbose:
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
                comm=comm_world, model_nprocs=params['model_nprocs']
            )
            else:
                display_timing_default(total_elapsed_time, total_wall_time)
        else:
            None

    except Exception as e:
        # Handle exceptions and print error messages
        comm_world.Barrier()
        if rank_world == 0:
            print(f"[ICESEE] An error occurred on in icesee_model_data_assimilation_full_parallel: {str(e)}")
            print(f"[ICESEE] You can restart from the previous checkpoint if enabled.")
        
        try:
            # Try to salvage k and km if present
            cur_k = model_kwargs.get("k", None)
            cur_km = model_kwargs.get("km", None)
            if restart_enabled and rank_world == 0 and cur_k is not None:
                ck = {
                    "last_done_k": max(int(cur_k) - 1, -1),  # last fully done; conservative
                    "km": int(cur_km) if cur_km is not None else None,
                    "nt": int(model_kwargs.get("nt", params["nt"])),
                    "nd": int(model_kwargs.get("nd", params["nd"])),
                    "nens": int(model_kwargs.get("Nens", params["Nens"])),
                    "dataset_dir": os.path.abspath(_modelrun_datasets),
                    "timestamp": time.time(),
                    "base_seed": int(base_seed),
                    "crash_message": str(e),
                }
                save_checkpoint(_modelrun_datasets, **ck)
                print(f"[ICESEE][RESTART] Checkpoint saved after error; you can restart safely.")
        except Exception as _ckerr:
            if rank_world == 0:
                print(f"[ICESEE][WARN] Could not save crash checkpoint: {_ckerr}")

        # close the EnKF I/O handler
        enkf_parallel_io.close()
        comm_world.Barrier()
        if model_kwargs.get("create_ensemble_dataset", True):
            if rank_world == 0:
                print("[ICESEE] Creating ensemble dataset...")
                # Option A: no-copy, instant
                out_vds = finalize_stack("_modelrun_datasets", mode="vds", dset_name="states")
                print("VDS ready:", out_vds)
                # Option B: portable single file
                # out_h5 = finalize_stack("_modelrun_datasets", mode="h5", dset_name="states",
                #                         allow_missing=False, compression="gzip", compression_opts=4)
                
            comm_world.Barrier()
        tb_str = "".join(traceback.format_exception(*sys.exc_info()))
        print(f"Traceback details:\n{tb_str}")
        # comm_world.Abort(1)  # Abort all processes in the communicator


