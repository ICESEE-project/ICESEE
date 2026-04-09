# ==============================================================================
# @des: This file contains run functions for the ICESEE model to generate true and nurged states.
# @date: 2025-07-30
# @author: Brian Kyanjo
# ==============================================================================

# --- import necessary libraries ---
import numpy as np
import h5py
import gc
from mpi4py import MPI

from ICESEE.src.utils.utils import UtilsFunctions
from ICESEE.src.utils.tools import icesee_get_index


def generate_synthetic_observations(**model_kwargs):
    """Generate synthetic observations for the ICESEE model.
    """

    # unpack model_kwargs
    params         = model_kwargs.get("params", {})
    model_module   = model_kwargs.get("model_module", None)
    comm_world     = model_kwargs.get("comm_world", MPI.COMM_WORLD)
    _synthetic_obs = model_kwargs.get("synthetic_obs_file")
    _true_nurged   = model_kwargs.get("true_nurged_file")
    color          = model_kwargs.get("color", 0)
    subcomm        = model_kwargs.get("subcomm", None)
    sub_rank       = model_kwargs.get("sub_rank", 0)
    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()


    if model_kwargs.get("generate_synthetic_obs", True):   
        if params["even_distribution"] or (params["default_run"] and size_world <= params["Nens"]):
            if rank_world == 0:
                # --- Synthetic Observations ---
                print("[ICESEE] Generating synthetic observations ...")
                with h5py.File(_true_nurged, "r") as f:
                    ensemble_true_state = f['true_state'][:]

                utils_funs = UtilsFunctions(params, ensemble_true_state)
                model_kwargs.update({"statevec_true": ensemble_true_state})
                hu_obs, error_R, bed_masks, kwargs = utils_funs._create_synthetic_observations(**model_kwargs)

                # check if the best_mask_map is generated
                model_kwargs.update({"bed_mask_map": bed_masks})
                
                # observe or don't observe parameters.
                vecs, indx_map,_ = icesee_get_index(hu_obs, **model_kwargs)
                # check if model_kwargs['observe_params'] is empty
                if len(model_kwargs['observed_params']) == 0:
                    for key in model_kwargs['params_vec']:
                        hu_obs[indx_map[key],:] = 0.0
                        error_R[:,indx_map[key]] = 0.0
                else: 
                    for key in model_kwargs['params_vec']:
                        if key not in model_kwargs['observed_params']:
                            hu_obs[indx_map[key],:] = 0.0
                            error_R[:,indx_map[key]] = 0.0

                # -- write data to file
                with h5py.File(_synthetic_obs, "w") as f:
                    f.create_dataset("hu_obs", data=hu_obs)
                    f.create_dataset("R", data=error_R)

                    # ---- bed masks ----
                    g_masks = f.create_group("bed_masks")

                    g_static = g_masks.create_group("static")
                    for key, mask in bed_masks["static"].items():
                        g_static.create_dataset(
                            key,
                            data=np.asarray(mask, dtype=np.uint8),
                            compression="gzip",
                            compression_opts=4,
                        )

                    g_cols = g_masks.create_group("cols")
                    for key, mask_cols in bed_masks["cols"].items():
                        g_cols.create_dataset(
                            key,
                            data=np.asarray(mask_cols, dtype=np.uint8),
                            compression="gzip",
                            compression_opts=4,
                        )

                    # ---- metadata needed to rebuild H consistently ----
                    f.create_dataset("bed_snap_cols", data=np.asarray(kwargs["bed_snap_cols"], dtype=int))
                    f.create_dataset("ind_m", data=np.asarray(kwargs["ind_m"], dtype=int))
                    f.create_dataset("obs_t", data=np.asarray(kwargs["obs_t"], dtype=float))

                    # obs_model_to_col is a dict -> store as parallel arrays
                    m = kwargs.get("obs_model_to_col", {})
                    keys = np.asarray(list(m.keys()), dtype=int)
                    vals = np.asarray([m[k] for k in keys], dtype=int)
                    f.create_dataset("obs_model_to_col_keys", data=keys)
                    f.create_dataset("obs_model_to_col_vals", data=vals)

                # --- clear memory
                del hu_obs
                del error_R
                gc.collect()

            else:
                pass
                # hu_obs = np.empty((params["nd"],params["number_obs_instants"]),dtype=np.float64)
                # error_R = np.empty((params["number_obs_instants"], params["nd"]),dtype=np.float64)

            if params["even_distribution"]:
                # Bcast the observations
                comm_world.Bcast(hu_obs, root=0)
            else:
                pass
                # hu_obs = comm_world.bcast(hu_obs, root=0)
                # error_R = comm_world.bcast(error_R, root=0)
                # *--- write observations to file ---
                # parallel_write_data_from_root_2D(full_ensemble=hu_obs, comm=comm_world, data_name='hu_obs', output_file="icesee_ensemble_data.h5")
        else:
            # --- Synthetic Observations ---
            if rank_world == 0:
                print("[ICESEE] Generating synthetic observations ...")

            if params["default_run"] and size_world > params["Nens"]:
                subcomm.Barrier()
                # comm_world.Bcast(hu_obs, root=0)
                if sub_rank == 0:
                    utils_funs = UtilsFunctions(params, ensemble_true_state)
                    model_kwargs.update({"statevec_true": ensemble_true_state})
                    hu_obs, error_R, bed_mask_map, kwargs = utils_funs._create_synthetic_observations(**model_kwargs)
                    model_kwargs.update({"bed_mask_map": bed_mask_map})

                    # observe or don't observe parameters.
                    vecs, indx_map,_ = icesee_get_index(hu_obs, **model_kwargs)
                    # check if model_kwargs['observe_params'] is empty
                    if len(model_kwargs['observed_params']) == 0:
                        for key in model_kwargs['params_vec']:
                            hu_obs[indx_map[key],:] = 0.0 
                            error_R[:,indx_map[key]] = 0.0
                    else: 
                        for key in model_kwargs['params_vec']:
                            if key not in model_kwargs['observed_params']:
                                hu_obs[indx_map[key],:] = 0.0
                                error_R[:,indx_map[key]] = 0.0

                    shape_ = np.array(hu_obs.shape,dtype=np.int32)
                    shape_R = np.array(error_R.shape,dtype=np.int32)

                    # write data to the file
                    with h5py.File(_synthetic_obs, 'w', driver='mpio', comm=subcomm) as f:
                        f.create_dataset("hu_obs", data=hu_obs)
                        f.create_dataset("R", data=error_R)
                else:
                    shape_ = np.empty(2,dtype=np.int32)
                    shape_R = np.empty(2,dtype=np.int32)

                subcomm.Bcast(shape_, root=0)
                subcomm.Bcast(shape_R, root=0)
                if sub_rank != 0:
                    hu_obs = np.empty(shape_,dtype=np.float64)
                    error_R = np.empty(shape_R,dtype=np.float64)


                # bcast the synthetic observations
                # subcomm.Bcast(hu_obs, root=0)
                # subcomm.Bcast(error_R, root=0)
                #- write observations to file
                # parallel_write_data_from_root_2D(full_ensemble=hu_obs, comm=subcomm, data_name='hu_obs', output_file="icesee_ensemble_data.h5")


                # broadcast to the global communicator
                # comm_world.Bcast(hu_obs, root=0)
                # print(f"[ICESEE] rank {rank_world} Shape of the observations: {hu_obs.shape}")
                # exit()    
            elif params["sequential_run"]:
                comm_world.Barrier()
                # g_shape = model_kwargs['dim_list'][rank_world]
                # utils_funs = UtilsFunctions(params, ensemble_true_state)
                # model_kwargs.update({"statevec_true": ensemble_true_state})
                # hu_obs = utils_funs._create_synthetic_observations(**model_kwargs)
                # # gather from every rank to rank 0
                # gathered_obs = comm_world.gather(hu_obs[:g_shape,:], root=0)
                # if rank_world == 0:
                #     print(f"[ICESEE] {[arr.shape for arr in gathered_obs]}")
                #     hu_obs = np.vstack(gathered_obs)
                # else:
                #     hu_obs = np.empty((model_kwargs["global_shape"],params["number_obs_instants"]),dtype=np.float64)
                
                # comm_world.Bcast(hu_obs, root=0)
                if rank_world == 0:
                    utils_funs = UtilsFunctions(params, ensemble_true_state)
                    model_kwargs.update({"statevec_true": ensemble_true_state})
                    hu_obs, error_R, bed_mask_map, kwargs = utils_funs._create_synthetic_observations(**model_kwargs)
                    model_kwargs.update({"bed_mask_map": bed_mask_map})
                    shape_ = np.array(hu_obs.shape,dtype=np.int32)
                    shape_R = np.array(error_R.shape,dtype=np.int32)
                else:
                    shape_ = np.empty(2,dtype=np.int32)
                    shape_R = np.empty(2,dtype=np.int32)

                comm_world.Bcast(shape_, root=0)
                comm_world.Bcast(shape_R, root=0)

                if rank_world != 0:
                    hu_obs = np.empty(shape_,dtype=np.float64)
                    error_R = np.empty(shape_R,dtype=np.float64)

                # bcast the synthetic observations
                comm_world.Bcast(hu_obs, root=0)
                comm_world.Bcast(error_R, root=0)

    return model_kwargs