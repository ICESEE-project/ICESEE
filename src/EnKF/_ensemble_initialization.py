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

def ensemble_initialization(**icesee_kwargs):
    """Initialize the ensemble for the ICESEE model.
    """

    # unpack icesee_kwargs
    model_module   = icesee_kwargs.get("model_module", None)
    pos            = icesee_kwargs.get("pos", None)
    gs_model       = icesee_kwargs.get("gs_model", None)
    L_C           = icesee_kwargs.get("L_C", None)
    Lx             = icesee_kwargs.get("Lx", 1.0)
    Ly             = icesee_kwargs.get("Ly", 1.0)
    len_scale      = icesee_kwargs.get("len_scale", 1.0)
    Q_rho          = icesee_kwargs.get("Q_rho", 1.0)
    model_nprocs   = icesee_kwargs.get("model_nprocs", 1)
    total_cores    = icesee_kwargs.get("total_cores", 1)
    base_total_procs = icesee_kwargs.get("base_total_procs", 1)
    rng           = icesee_kwargs.get("rng", np.random.default_rng())
    rank_seed = icesee_kwargs.get("rank_seed", 0)

    rank_world = 0
    size_world = 1

    time_init_noise_generation = 0.0
    time_init_file_writing     = 0.0
    time_init_ensemble_mean_computation = 0.0

    nd = icesee_kwargs.get("nd", icesee_kwargs["nd"])
    Nens = icesee_kwargs.get("Nens", icesee_kwargs["Nens"])

    if rank_world == 0:
        print("[ICESEE] Initializing the ensemble ...")
        icesee_kwargs.update({'ens_id': rank_world})

        icesee_kwargs.update({"statevec_ens":np.zeros([nd, Nens])})

        # get the ensemble matrix
        vecs, indx_map, dim_per_proc = icesee_get_index(icesee_kwargs["statevec_ens"], **icesee_kwargs)
        ensemble_vec = np.zeros_like(icesee_kwargs["statevec_ens"])

        hdim = ensemble_vec.shape[0] // icesee_kwargs["total_state_param_vars"]

        state_block_size = hdim * icesee_kwargs["num_state_vars"]

        # # --- get the process noise ---
        # pos, gs_model, L_C = compute_Q_err_random_fields(hdim, icesee_kwargs["total_state_param_vars"], icesee_kwargs["sig_Q"], Q_rho, len_scale)

        # process_noise = []
        for ens in range(icesee_kwargs["Nens"]):
            # icesee_kwargs.update({"ens_id": ens})
            data = model_module.initialize_ensemble(ens,**icesee_kwargs)

            # iterate over the data and update the ensemble
            for key, value in data.items():
                ensemble_vec[indx_map[key],ens] = value

            # --->
            # noise = compute_noise_random_fields(ens, hdim, pos, gs_model, icesee_kwargs["total_state_param_vars"], L_C)
            # ensemble_vec[:,ens] += noise
            #----->
            _time_init_noise_generation = time.time()
            N_size = icesee_kwargs["total_state_param_vars"] * hdim
            # noise = generate_pseudo_random_field_1d(N_size,np.sqrt(Lx*Ly), len_scale, verbose=True)
            icesee_kwargs.update({"ii_sig": None, "hdim":hdim, "num_vars":icesee_kwargs["total_state_param_vars"]})
            # noise = generate_enkf_field(**icesee_kwargs)

            if (len(icesee_kwargs.get("scalar_inputs", [])) > 0) or (icesee_kwargs.get("var_nd", None) is not None):
                icesee_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim})
                noise_1 = generate_enkf_field(**icesee_kwargs)
                ndim = 1 if len(icesee_kwargs.get("scalar_inputs", [])) > 0 else (icesee_kwargs["var_nd"][icesee_kwargs["scalar_inputs"][0]])
                icesee_kwargs.update({ "noise_dim": ndim})
                noise_2 = generate_enkf_field(**icesee_kwargs)
                # concatenate noise_1 and noise_2
                noise = np.concatenate((noise_1, noise_2))[:-1]

            else:
                icesee_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim})
                noise = generate_enkf_field(**icesee_kwargs)

            time_init_noise_generation += time.time() - _time_init_noise_generation
            # print(f"\nensemble_vec[:,{ens}]: {ensemble_vec[:,ens]} noise: {noise}, hdim: {hdim} Lx: {Lx}, Ly: {Ly}, len_scale: {len_scale}, total_params: {icesee_kwargs['total_state_param_vars']}\n")
            ensemble_vec[:,ens] += noise
            # for ii, sig in enumerate(icesee_kwargs["sig_Q"]):
            #     if ii <=icesee_kwargs["num_state_vars"]:
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
            if icesee_kwargs.get('ICESEE_PERFORMANCE_TEST') or env_flag("ICESEE_PERFORMANCE_TEST", default=False):
                model_nprocs = model_nprocs
            else:
                model_nprocs = max(min_model_nprocs, model_nprocs + (diff // size_world))
        else:
            model_nprocs = model_nprocs

    icesee_kwargs.update({'model_nprocs': model_nprocs})

    # -- time ensemble mean computation ---
    _time_init_ensemble_mean_computation = time.time()
    ens_mean = np.mean(ensemble_vec, axis=1)
    time_init_ensemble_mean_computation += time.time() - _time_init_ensemble_mean_computation

    # ---time file writing ---
    _time_init_file_writing = time.time()
    # serial write from root
    output_file = os.path.join(icesee_kwargs.get('data_path'), "icesee_ensemble_data.h5")
    nd, Nens = ensemble_vec.shape
    with h5py.File(output_file, 'w') as f:
        # Create dataset with total dimensions
        dset = f.create_dataset('ensemble', (nd, Nens, icesee_kwargs.get('nt', icesee_kwargs['nt']) + 1), dtype='f8')
        # Write full ensemble
        dset[:, :, 0] = ensemble_vec[:,:Nens]

        # Create and write ensemble mean
        ensemble_mean = f.create_dataset('ensemble_mean', (nd, icesee_kwargs.get('nt', icesee_kwargs['nt']) + 1), dtype='f8')
        ensemble_mean[:, 0] = ens_mean

    time_init_file_writing += time.time() - _time_init_file_writing

    if icesee_kwargs.get("default_run", False):
        return icesee_kwargs, ensemble_vec, time_init_noise_generation, \
               time_init_ensemble_mean_computation, time_init_file_writing, \
                shape_ens, None, None, None
