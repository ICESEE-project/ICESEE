# ==============================================================================
# @des: This file contains run functions for the ICESEE model to initialize the ensemble. Serial version
# @date: 2025-07-30
# @author: Brian Kyanjo
# ==============================================================================

# --- import necessary libraries ---
import numpy as np
import h5py
import gc
import zarr
import os
import time

from ICESEE.src.utils.tools import icesee_get_index, env_flag
from ICESEE.src.run_model_da._error_generation import compute_Q_err_random_fields, \
                              compute_noise_random_fields, \
                              generate_pseudo_random_field_1d, \
                              generate_pseudo_random_field_2D, \
                              generate_enkf_field

def ensemble_initialization(**model_kwargs):
    """Initialize the ensemble for the ICESEE model.
    """
    
    # unpack model_kwargs
    params         = model_kwargs.get("params")
    model_module   = model_kwargs.get("model_module", None)
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
    rng           = model_kwargs.get("rng", np.random.default_rng())
    rank_seed = model_kwargs.get("rank_seed", 0)

    rank_world = 0
    size_world = 1

    time_init_noise_generation = 0.0
    time_init_file_writing     = 0.0
    time_init_ensemble_mean_computation = 0.0

    nd = model_kwargs.get("nd", params["nd"])
    Nens = model_kwargs.get("Nens", params["Nens"])

    if rank_world == 0:
        print("[ICESEE] Initializing the ensemble ...")
        model_kwargs.update({'ens_id': rank_world})

        model_kwargs.update({"statevec_ens":np.zeros([nd, Nens])})
        
        # get the ensemble matrix   
        vecs, indx_map, dim_per_proc = icesee_get_index(model_kwargs["statevec_ens"], **model_kwargs)
        ensemble_vec = np.zeros_like(model_kwargs["statevec_ens"])

        hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
        
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
            _time_init_noise_generation = time.time()
            N_size = params["total_state_param_vars"] * hdim
            # noise = generate_pseudo_random_field_1d(N_size,np.sqrt(Lx*Ly), len_scale, verbose=True)
            model_kwargs.update({"ii_sig": None, "hdim":hdim, "num_vars":params["total_state_param_vars"]})
            # noise = generate_enkf_field(**model_kwargs)

            if (len(model_kwargs.get("scalar_inputs", [])) > 0) or (model_kwargs.get("var_nd", None) is not None):
                noise_1 = generate_enkf_field(None, np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)
                ndim = 1 if len(model_kwargs.get("scalar_inputs", [])) > 0 else (model_kwargs["var_nd"][model_kwargs["scalar_inputs"][0]])
                noise_2 = generate_enkf_field(None, np.sqrt(Lx*Ly), ndim, params["total_state_param_vars"], rh=len_scale, verbose=False)
                # concatenate noise_1 and noise_2 
                noise = np.concatenate((noise_1, noise_2))[:-1]

            else:
                noise = generate_enkf_field(None, np.sqrt(Lx*Ly), hdim, params["total_state_param_vars"], rh=len_scale, verbose=False)
            
            time_init_noise_generation += time.time() - _time_init_noise_generation
            # print(f"\nensemble_vec[:,{ens}]: {ensemble_vec[:,ens]} noise: {noise}, hdim: {hdim} Lx: {Lx}, Ly: {Ly}, len_scale: {len_scale}, total_params: {params['total_state_param_vars']}\n")
            ensemble_vec[:,ens] += noise
            # for ii, sig in enumerate(params["sig_Q"]):
            #     if ii <=params["num_state_vars"]:
            #         start_idx = ii * hdim
            #         end_idx = start_idx + hdim
            #         ensemble_vec[start_idx:end_idx, ens] += noise[start_idx:end_idx] * sig
            # print(f"\nensemble_vec[:,{ens}]: {ensemble_vec[:,ens]}\n")
        shape_ens = np.array(ensemble_vec.shape,dtype=np.int32)
    
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

    model_kwargs.update({'model_nprocs': model_nprocs})

    # -- time ensemble mean computation ---
    _time_init_ensemble_mean_computation = time.time()
    ens_mean = np.mean(ensemble_vec, axis=1)
    time_init_ensemble_mean_computation += time.time() - _time_init_ensemble_mean_computation

    # ---time file writing ---
    _time_init_file_writing = time.time()
    # serial write from root
    output_file = os.path.join(params.get('data_path'), "icesee_ensemble_data.h5")
    nd, Nens = ensemble_vec.shape
    with h5py.File(output_file, 'w') as f:
        # Create dataset with total dimensions
        dset = f.create_dataset('ensemble', (nd, Nens, model_kwargs.get('nt', params['nt']) + 1), dtype='f8')
        # Write full ensemble
        dset[:, :, 0] = ensemble_vec[:,:Nens]

        # Create and write ensemble mean
        ensemble_mean = f.create_dataset('ensemble_mean', (nd, model_kwargs.get('nt', params['nt']) + 1), dtype='f8')
        ensemble_mean[:, 0] = ens_mean

    time_init_file_writing += time.time() - _time_init_file_writing

    if params.get("default_run", False):
        return model_kwargs, ensemble_vec, time_init_noise_generation, \
               time_init_ensemble_mean_computation, time_init_file_writing, \
                shape_ens, None, None, None
