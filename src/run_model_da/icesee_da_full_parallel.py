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
from ICESEE.src.utils.icesee_context import normalize_execution_mode, normalize_icesee_kwargs
from ICESEE.src.utils.localization import prepare_random_field_coordinates

# --- call the ICESEE mpi parallel manager ---
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager
from ICESEE.src.parallelization._mpi_forecast_functions import parallel_forecast_step_default_full_parallel_run
from ICESEE.src.parallelization._mpi_generate_true_wrong_state import generate_true_wrong_state
from ICESEE.src.parallelization._mpi_generate_synthetic_observations import synchronize_observation_schedule
from ICESEE.src.parallelization._mpi_ensemble_intialization import ensemble_initialization_full_parallel_run
from ICESEE.src.parallelization.EnKF_parallel_io import EnKF_fully_parallel_IO

# ======================== Run model with EnKF ========================
def icesee_model_data_assimilation_full_parallel(**icesee_kwargs):
    """ General function to run any kind of model with the Ensemble Kalman Filter """

    icesee_kwargs = normalize_icesee_kwargs(icesee_kwargs)

    # --- unpack the data assimilation arguments
    filter_type       = icesee_kwargs.get("filter_type", "EnKF")      # filter type
    model             = icesee_kwargs.get("model_name",None)          # model name
    execution_mode    = normalize_execution_mode(icesee_kwargs, expected=2)
    Q_err             = icesee_kwargs.get("Q_err",None)               # process noise
    commandlinerun    = icesee_kwargs.get("commandlinerun",None)      # run through the terminal
    Lx, Ly            = icesee_kwargs.get("Lx",1.0), icesee_kwargs.get("Ly",1.0)
    nx, ny            = icesee_kwargs.get("nx",1), icesee_kwargs.get("ny",1)
    b_in, b_out       = icesee_kwargs.get("b_in",0.0), icesee_kwargs.get("b_out",0.0)
    data_path         = icesee_kwargs.get("data_path","_modelrun_datasets")      # data path
    restart_enabled   = icesee_kwargs.get("restart_enabled", True)   # turn on/off restart
    force_fresh_start = icesee_kwargs.get("force_fresh_start", False) # ignore old files/ckpt
    checkpoint_every  = icesee_kwargs.get("checkpoint_every", 1)     # write ckpt every N steps
    base_seed         = int(icesee_kwargs.get("base_seed", 42))       # for reproducible reseed
    nd               = icesee_kwargs["nd"]      # model dimension
    Nens             = icesee_kwargs["Nens"]    # number of ensemble members
    nt               = icesee_kwargs["nt"]      # number of time steps


    # start the timer
    global_start_time = MPI.Wtime()

    # --- icesee mpi parallel manager ---------------------------------------------------
    # --- ensemble load distribution --
    rounds, color, sub_rank, sub_size, subcomm, subcomm_size_min, rank_world, size_world, comm_world, start, stop = ParallelManager().icesee_mpi_ens_distribution(icesee_kwargs)
    icesee_kwargs.update({'size_world': size_world, 'comm_world': comm_world})

    # --- call curently supported model Class
    model_module = SupportedModels(model=model,comm=comm_world,verbose=icesee_kwargs.get('verbose')).call_model()
    # pack the global communicator, the subcommunicator and other important parameters
    icesee_kwargs.update({"comm_world": comm_world, "subcomm": subcomm,
                            "rank_world": rank_world, "sub_rank": sub_rank,
                            "size_world": size_world, "sub_size": sub_size,
                            "rounds": rounds, "color": color,
                            "start": start, "stop": stop,
                            "subcomm_size_min": subcomm_size_min,
                            "model_module": model_module,
                            'vec_inputs_old': icesee_kwargs.get('vec_inputs')})

    icesee_kwargs['observed_vars_params'] = (icesee_kwargs['observed_vars'] + icesee_kwargs['observed_params'])
    all_observed = icesee_kwargs['observed_vars_params']

    icesee_kwargs.update({'all_observed': all_observed}); icesee_kwargs.update({'all_observed': all_observed})

    # pack the global communicator and the subcommunicator
    icesee_kwargs.update({"comm_world": comm_world, "subcomm": subcomm})

    # --- check if the modelrun dataset directory is present ---
    _modelrun_datasets = icesee_kwargs.get("data_path") or "_modelrun_datasets"
    os.environ["ICESEE_RESULTS_DIR"] = str(_modelrun_datasets)
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

    # --update icesee_kwargs with the file names
    icesee_kwargs.update({"true_nurged_file": _true_nurged, "synthetic_obs_file": _synthetic_obs})

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
    # Timestep files contain the complete ensemble.  Creating a large batch
    # up front can reserve tens or hundreds of GiB before those steps run.
    # Keep mode 2 lazy by default; users can raise this only when metadata
    # latency dominates and storage is known to be sufficient.
    batch_size = max(1, int(icesee_kwargs.get("batch_size", 2)))
    requested_history_mode = str(
        icesee_kwargs.get("ensemble_history_mode", "auto")
    ).strip().lower()
    if requested_history_mode not in {"auto", "full", "rolling"}:
        raise ValueError(
            "ensemble_history_mode must be 'auto', 'full', or 'rolling'."
        )
    serial_file_creation = icesee_kwargs.get("serial_file_creation",True)
    h5_file_compression = icesee_kwargs.get("h5_file_compression",None)
    h5_file_compression_level = icesee_kwargs.get("h5_file_compression_level",4)
    h5_file_chunk_size = icesee_kwargs.get("h5_file_chunk_size",1000)
    enkf_parallel_io = EnKF_fully_parallel_IO('icesee_enkf_ens', nd, Nens, nt, subcomm, comm_world, \
                                             icesee_kwargs, serial_file_creation, base_path=_modelrun_datasets, \
                                             batch_size=batch_size, h5_file_compression=h5_file_compression, \
                                             h5_file_compression_level=h5_file_compression_level, \
                                             h5_file_chunk_size=h5_file_chunk_size)
    history_mode = enkf_parallel_io.history_mode
    icesee_kwargs["ensemble_history_mode"] = history_mode
    # Update icesee_kwargs with the EnKF I/O handler
    icesee_kwargs.update({"enkf_parallel_io": enkf_parallel_io})
    time_file_io_initialization = MPI.Wtime() - time_file_io_initialization

    try:

        # fetch model nprocs
        model_nprocs = icesee_kwargs.get("model_nprocs", 1)

        # set modeel_nprocs adaptively
        # total_cores = os.cpu_count()
        if icesee_kwargs.get('ICESEE_PERFORMANCE_TEST') or env_flag("ICESEE_PERFORMANCE_TEST", default=False):
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

        # update icesee_kwargs with the effective model_nprocs
        icesee_kwargs.update({'model_nprocs': effective_model_nprocs,
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
                or icesee_kwargs.get("allow_reuse_without_fingerprint", False):
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
            icesee_kwargs = generate_true_wrong_state(**icesee_kwargs)
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
            #                     or icesee_kwargs.get("allow_reuse_without_fingerprint", False):
            #                         syn_ok = True
            # except Exception:
            #     syn_ok = False

        if icesee_kwargs.get("generate_synthetic_obs", True) and not syn_ok:
            synthetic_obs_zarr_path = f"{_modelrun_datasets}/synthetic_observations.zarr"
            error_R_zarr_path = f"{_modelrun_datasets}/error_R.zarr"
            icesee_kwargs.update({'synthetic_obs_zarr_path': synthetic_obs_zarr_path, 'error_R_zarr_path': error_R_zarr_path})
            tobserve, m_obs = enkf_parallel_io._create_synthetic_observations(**icesee_kwargs)
            # if rank_world == 0:
            #     with h5py.File(_synthetic_obs, "w") as f:
            #         f.create_dataset("tobserve", data=np.asarray(tobserve, dtype=np.int64))
            #         f.create_dataset("m_obs", data=np.asarray([m_obs], dtype=np.int64))
            #         f.attrs["icesee_fingerprint"] = fp
        else:
            if rank_world == 0:
                print("[ICESEE][RESTART] Using existing synthetic observations.")
            _, tobserve, m_obs = enkf_parallel_io.generate_observation_schedule(**icesee_kwargs)

        # Treat the schedule stored beside the generated observation columns as
        # authoritative, then broadcast it to all ranks.  This is shared with
        # execution mode 1 and prevents analysis from occurring at stale steps.
        icesee_kwargs = synchronize_observation_schedule(
            icesee_kwargs, obs_file=_synthetic_obs
        )
        tobserve = icesee_kwargs["tobserve"]
        m_obs = icesee_kwargs["m_obs"]
        time_generation_synthetic_obs = MPI.Wtime() - time_generation_synthetic_obs
        comm_world.Barrier()

        # ----- Decide restart step (supports explicit override) -----
        k_start = 0

        # User override (highest priority)
        k_start_override = icesee_kwargs.get("k_start_override", None)
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
                        # Shard t is the input state for model cycle t.  In
                        # the absence of a checkpoint, resume *from* that
                        # state rather than advancing it a second time.
                        k_start = int(last_k)

        # Clamp
        k_start = min(max(0, k_start), nt)

        # (Optional) safety: ensure files before k_start exist contiguously
        if (history_mode == "full"
                and icesee_kwargs.get("enforce_contiguous_history", True)
                and k_start > 0):
            missing = []
            for kk in range(k_start):
                if not step_already_done(_modelrun_datasets, kk):
                    missing.append(kk)
            if missing:
                raise RuntimeError(
                    f"Cannot start at k={k_start}: missing completed steps {missing[:10]}{' ...' if len(missing)>10 else ''}. "
                    "Either lower k_start_override or disable enforce_contiguous_history."
                )

        # (Optional) remove outputs newer than the requested restart state.
        # The shard at ``k_start`` is the input state for the replayed cycle
        # and must be retained.
        if icesee_kwargs.get("truncate_after_k_start", False):
            # Delete step files strictly newer than k_start.
            for fname in _sorted_step_files(_modelrun_datasets):
                kk = _extract_time_from_name(fname)
                if kk > k_start:
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
        icesee_kwargs.update({"km": km})

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
        #         icesee_kwargs.update({'H_matrix_zarr_path': H_matrix_zarr_path})
        #         enkf_parallel_io.H_matrix(**icesee_kwargs)
        #         # record a tiny meta tag
        #         with h5py.File(meta_h5, "w") as f:
        #             f.attrs["icesee_fingerprint"] = fp

        icesee_kwargs.update({'H_matrix_zarr_path': H_matrix_zarr_path})
        enkf_parallel_io.H_matrix(**icesee_kwargs)
        # comm_world.Barrier()

        # --- Initialize the ensemble ---------------------------------------------------
        Q_rho     = icesee_kwargs.get("Q_rho")
        len_scale = icesee_kwargs.get("length_scale")
        hdim  = icesee_kwargs["nd"] // icesee_kwargs["total_state_param_vars"]
        icesee_kwargs.update({"hdim": hdim, "Q_rho": Q_rho, "len_scale": len_scale})
        prepare_random_field_coordinates(icesee_kwargs, expected_nodes=hdim)

            # --- get the process noise --->
        if icesee_kwargs.get("use_random_fields", False):
            pos, gs_model, L_C = compute_Q_err_random_fields(hdim, icesee_kwargs["total_state_param_vars"], icesee_kwargs["sig_Q"], Q_rho, len_scale)
            icesee_kwargs.update({"pos": pos, "gs_model": gs_model, "L_C": L_C})

        # -- time ensemble initialization ---
        time_ensemble_initialization = MPI.Wtime()

        restart_state_t = int(k_start) if k_start > 0 else 0
        run_already_complete = k_start >= nt
        init_ok = run_already_complete or (reuse_allowed and tools.h5_has_dataset_with_shape(
            os.path.join(
                _modelrun_datasets,
                f"icesee_enkf_ens_{restart_state_t:04d}.h5",
            ),
            "states", (nd, Nens)
        ))
        if init_ok:
            if rank_world == 0:
                if run_already_complete:
                    print("[ICESEE][RESTART] Run is already complete.")
                else:
                    print(
                        "[ICESEE][RESTART] Skipping ensemble initialization "
                        f"(found restart state {restart_state_t})."
                    )

            # load only dictionary essentials
            if size_world <= icesee_kwargs["Nens"]:
                icesee_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})
                dim_list = comm_world.allgather(icesee_kwargs.get("nd", icesee_kwargs["nd"]))
                icesee_kwargs.update({"global_shape": icesee_kwargs.get("nd", icesee_kwargs["nd"]), "dim_list": dim_list})
            else:
                icesee_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})
                icesee_kwargs.update({'ens_id': color}) # Nens = color
                # gather all the vector dimensions from all processors
                dim_list = subcomm.allgather(icesee_kwargs.get("nd", icesee_kwargs["nd"]))
                global_shape = sum(dim_list)
                icesee_kwargs.update({"global_shape": global_shape, "dim_list": dim_list})

            time_init_file_writing = 0.0
            time_init_noise_generation = 0.0
            time_init_ensemble_mean_computation = 0.0
        else:
            # call the ensemble_initialization function
            if icesee_kwargs.get("initialize_ensemble", True):
                icesee_kwargs, ensemble_vec, time_init_noise_generation, \
                time_init_ensemble_mean_computation, time_init_file_writing, \
                shape_ens,ensemble_bg,  ensemble_vec_mean, ensemble_vec_full = ensemble_initialization_full_parallel_run(**icesee_kwargs)
            else:
                # If ensemble initialization is disabled, set default values
                ensemble_vec = None
                time_init_noise_generation = 0.0
                time_init_ensemble_mean_computation = 0.0
                time_init_file_writing = 0.0
                shape_ens = (icesee_kwargs["nd"], icesee_kwargs["nd"])
                ensemble_bg = np.zeros(shape_ens)
                ensemble_vec_mean = np.zeros((icesee_kwargs["nd"], 1))
                ensemble_vec_full = np.zeros(shape_ens)

        # --- time ensemble initialization ---
        time_ensemble_initialization = MPI.Wtime() - time_ensemble_initialization

        # get updated model_nprocs
        model_nprocs = icesee_kwargs.get("model_nprocs", 1)

        # --- Define filter flags
        EnKF_flag   = re.match(r"\AEnKF\Z", filter_type, re.IGNORECASE)
        DEnKF_flag  = re.match(r"\ADEnKF\Z", filter_type, re.IGNORECASE)
        EnRSKF_flag = re.match(r"\AEnRSKF\Z", filter_type, re.IGNORECASE)
        EnTKF_flag  = re.match(r"\AEnTKF\Z", filter_type, re.IGNORECASE)

        # tqdm progress bar
        # Initialize progress bar on the root process
        if rank_world == 0:
            nt = icesee_kwargs.get("nt", icesee_kwargs["nt"])
            print(f"[ICESEE] Launching {model} with data assimilation using the {filter_type} filter across {size_world*(icesee_kwargs['model_nprocs']+1)} MPI ranks.")
            pbar = tqdm(
                total=nt,
                desc=f"[ICESEE] Assimilation progress ({size_world*(icesee_kwargs['model_nprocs']+1)} ranks)",
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
        dt  = icesee_kwargs.get("dt",icesee_kwargs["dt"])
        tau = max(max_tau,max(min_tau, dt))

        # tau = max(icesee_kwargs.get("dt",icesee_kwargs["dt"]),10)
        alpha = 1 - dt/tau
        # make sure  0=<alpha<1
        if alpha <= 0 or alpha > 1:
            alpha = 0.5

        n = icesee_kwargs.get("nt",icesee_kwargs["nt"])
        # rho = np.sqrt((1-alpha**2)/(dt*(n - 2*alpha - n*alpha**2 + 2*alpha**(n+1))))
        rho = np.sqrt((1/dt)*((1-alpha)**2)*(1/(n - (2*alpha) - (n*alpha**2) + (2*alpha**(n+1)))))
        params_analysis_0 = np.zeros((2, Nens))
        # km = 0

        #--- generate inital noise
        if icesee_kwargs.get("use_random_fields", False):
            # with h5py.File(_synthetic_obs, 'r') as f:
            #     error_R = f['error_R'][:]
            #     Cov_obs = np.cov(error_R)
            #  --- get the observation noise ---
            pos_obs, gs_model_obs, L_C_obs = compute_Q_err_random_fields(hdim, icesee_kwargs["total_state_param_vars"], icesee_kwargs["sig_obs"], Q_rho, len_scale)
        else:
            N_size = icesee_kwargs["total_state_param_vars"] * hdim
            # noise = generate_pseudo_random_field_1d(N_size,np.sqrt(Lx*Ly), len_scale, verbose=0)
            icesee_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":icesee_kwargs["total_state_param_vars"]})
            noise = generate_enkf_field(**icesee_kwargs)

        # synchronize all processes before starting the time loop
        comm_world.Barrier()

        # for k in range(icesee_kwargs.get("nt",icesee_kwargs["nt"])):
        for k in range(k_start, icesee_kwargs.get("nt",icesee_kwargs["nt"])): # resume from k_start

            # Deterministic reseed per step (optional but recommended)
            rank_seed = reseed_for_step(base_seed, rank_world, k)
            icesee_kwargs.update({"rank_seed": rank_seed})

            icesee_kwargs.update({"k": k, "km":km, "alpha": alpha, "rho": rho, "tau": tau, "dt": dt,"n": n})
            icesee_kwargs.update({"generate_enkf_field": generate_enkf_field}) #save the function to generate the enkf field

            if execution_mode == 2:
                # -- time forecast step ---
                _time_forecast_step = MPI.Wtime()

                # get the state block size
                ndim = nd//icesee_kwargs["total_state_param_vars"]
                state_block_size = ndim*icesee_kwargs["num_state_vars"]

                # load all needed parameters and variables into icesee_kwargs
                icesee_kwargs.update({"_modelrun_datasets": _modelrun_datasets,
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

                if icesee_kwargs["default_run"]:
                    # call the parallel_forecast_step_default_run function
                    icesee_kwargs = parallel_forecast_step_default_full_parallel_run(**icesee_kwargs)
                    time_forecast_step = icesee_kwargs.get("time_forecast_step", 0.0)
                    time_forecast_noise_generation = icesee_kwargs.get("time_forecast_noise_generation", 0.0)
                    time_forecast_file_writing = icesee_kwargs.get("time_forecast_file_writing", 0.0)
                    time_forecast_ensemble_mean_generation = icesee_kwargs.get("time_forecast_ensemble_mean_generation", 0.0)

                    comm_world.Barrier()
                    # print(f"[ICESEE] Rank {rank_world}, completed time step {k+1}/{icesee_kwargs['nt']} with forecast time {time_forecast_step:.2f}s.")
                    # --- end time forecast step
                    time_forecast_step += MPI.Wtime() - _time_forecast_step

                    # ===== Global analysis step =====
                    if icesee_kwargs.get('global_analysis', True) or icesee_kwargs.get('local_analysis', False):

                        tobserve = icesee_kwargs.get("tobserve")
                        m_obs = icesee_kwargs.get("m_obs", icesee_kwargs["number_obs_instants"])
                        # if (km < m_obs) and (k+1 == tobserve[km]):
                        # if (km < m_obs) and (k == tobserve[km]):
                        obs_index = icesee_kwargs["obs_index"]
                        if (km < icesee_kwargs["number_obs_instants"]) and (k == obs_index[km]):
                            # -- time global analysis step ---
                            _time_analysis_step = MPI.Wtime()
                            icesee_kwargs.update({'km': km, 'k': k})

                            inversion_enabled = bool(icesee_kwargs.get(
                                "inversion_enabled",
                                icesee_kwargs.get("inversion_flag", False),
                            ))
                            inversion_start_time = float(
                                icesee_kwargs.get("inversion_start_time", 0.0)
                            )
                            cycle_time = float(np.asarray(icesee_kwargs["t"])[k])
                            inversion_flag = (
                                inversion_enabled
                                and cycle_time + 1.0e-12 >= inversion_start_time
                            )
                            icesee_kwargs["inversion_flag"] = inversion_flag
                            if rank_world == 0 and inversion_enabled and not inversion_flag:
                                print(
                                    "[ICESEE] Deferring friction inversion at "
                                    f"t={cycle_time:g} yr; configured start is "
                                    f"{inversion_start_time:g} yr."
                                )
                            nd_old = icesee_kwargs.get("nd", nd)
                            icesee_kwargs.update({"nd_old": nd_old})

                            # call the analysis update function
                            if EnKF_flag:
                                icesee_kwargs = enkf_parallel_io.compute_analysis_update(**icesee_kwargs)
                                time_analysis_ensemble_mean_generation = icesee_kwargs.get("time_analysis_ensemble_mean_generation", 0.0)
                                time_analysis_file_writing = icesee_kwargs.get("time_analysis_file_writing", 0.0)

                            # update the observation index
                            km += 1
        #
                            # --- end time analysis step ---
                            time_analysis_step += MPI.Wtime() - _time_analysis_step

                    # Step k is now fully complete.  Like execution mode 1,
                    # mode 2 retains the initial state at shard 0 and writes
                    # every advanced state to shard k + 1, including the
                    # terminal state at shard nt.
                    write_k = k + 1

                    # Checkpoint before pruning so the retained rolling shard
                    # and last_done_k always describe the same restart point.
                    if restart_enabled and (
                        history_mode == "rolling"
                        or k % checkpoint_every == 0
                        or k == nt - 1
                    ):
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
                                "ensemble_state_t": int(write_k),
                            }
                            try:
                                save_checkpoint(_modelrun_datasets, **ck)
                            except Exception as e:
                                print(f"[ICESEE][WARN] Failed to save checkpoint at k={k}: {e}")

                    # Bound disk usage independently of run length.  This is
                    # collective because active HDF5 handles must be closed on
                    # every rank before an old shard is unlinked.
                    enkf_parallel_io.prune_history(write_k)

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
            icesee_kwargs,
            data_path=_modelrun_datasets,
            nofilter=True,
            t=icesee_kwargs["t"],
            b_io=np.array([b_in, b_out]),
            Lxy=np.array([Lx, Ly]),
            nxy=np.array([nx, ny]),
            obs_max_time=np.array([icesee_kwargs["obs_max_time"]]),
            obs_index=icesee_kwargs["obs_index"],
            run_mode=np.array([icesee_kwargs["execution_flag"]]),
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
                if icesee_kwargs.get("create_ensemble_dataset", True):
                    finalize_mode = str(icesee_kwargs.get(
                        "ensemble_finalize_mode", "vds"
                    )).lower()
                    if history_mode == "rolling":
                        print(
                            "[ICESEE] Rolling ensemble history enabled; "
                            "skipping full-history ensemble finalization."
                        )
                    elif finalize_mode not in {"none", "off", "false"}:
                        print(f"[ICESEE] Creating {finalize_mode} ensemble view...")
                        out_h5 = finalize_stack(
                            _modelrun_datasets, mode=finalize_mode,
                            dset_name="states",
                            row_chunk_size=int(icesee_kwargs.get(
                                "finalize_row_chunk_size", 16384
                            )),
                        )
                        print("Ensemble view ready:", out_h5)
                # --- remove all .zarr files ---
                cleanup_intermediates = icesee_kwargs.get("cleanup_intermediates", True)
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
            verbose = icesee_kwargs.get("verbose", False)
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
                comm=comm_world, model_nprocs=icesee_kwargs['model_nprocs']
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
            cur_k = icesee_kwargs.get("k", None)
            cur_km = icesee_kwargs.get("km", None)
            if (restart_enabled and history_mode != "rolling"
                    and rank_world == 0 and cur_k is not None):
                ck = {
                    "last_done_k": max(int(cur_k) - 1, -1),  # last fully done; conservative
                    "km": int(cur_km) if cur_km is not None else None,
                    "nt": int(icesee_kwargs.get("nt", icesee_kwargs["nt"])),
                    "nd": int(icesee_kwargs.get("nd", icesee_kwargs["nd"])),
                    "nens": int(icesee_kwargs.get("Nens", icesee_kwargs["Nens"])),
                    "dataset_dir": os.path.abspath(_modelrun_datasets),
                    "timestamp": time.time(),
                    "base_seed": int(base_seed),
                    "crash_message": str(e),
                }
                save_checkpoint(_modelrun_datasets, **ck)
                print(f"[ICESEE][RESTART] Checkpoint saved after error; you can restart safely.")
            elif restart_enabled and history_mode == "rolling" and rank_world == 0:
                print(
                    "[ICESEE][RESTART] Retaining the last completed rolling "
                    "checkpoint after the error."
                )
        except Exception as _ckerr:
            if rank_world == 0:
                print(f"[ICESEE][WARN] Could not save crash checkpoint: {_ckerr}")

        # close the EnKF I/O handler
        enkf_parallel_io.close()
        comm_world.Barrier()
        if icesee_kwargs.get("create_ensemble_dataset", True):
            if rank_world == 0:
                finalize_mode = str(icesee_kwargs.get(
                    "ensemble_finalize_mode", "vds"
                )).lower()
                if history_mode == "rolling":
                    print(
                        "[ICESEE] Rolling ensemble history enabled; no "
                        "full-history recovery view can be created."
                    )
                elif finalize_mode not in {"none", "off", "false"}:
                    print(f"[ICESEE] Creating {finalize_mode} recovery view...")
                    out_vds = finalize_stack(
                        _modelrun_datasets, mode=finalize_mode,
                        dset_name="states",
                        row_chunk_size=int(icesee_kwargs.get(
                            "finalize_row_chunk_size", 16384
                        )),
                    )
                    print("Recovery view ready:", out_vds)

            comm_world.Barrier()
        tb_str = "".join(traceback.format_exception(*sys.exc_info()))
        print(f"Traceback details:\n{tb_str}")
        # comm_world.Abort(1)  # Abort all processes in the communicator
