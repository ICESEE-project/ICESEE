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


def synchronize_observation_schedule(icesee_kwargs, obs_file=None):
    """Load and broadcast the canonical observation schedule.

    Synthetic observations are written by one MPI rank, while every rank
    evaluates the forecast/analysis schedule.  Returning the rank-local
    ``icesee_kwargs`` from the generator can therefore leave non-writer ranks
    with the schedule computed during configuration parsing.  The observation
    artifact is the authoritative source because its columns were generated
    with exactly these indices.

    This helper deliberately mutates and returns ``icesee_kwargs`` so both
    execution modes use identical schedule metadata.
    """
    comm = icesee_kwargs.get("comm_world", MPI.COMM_WORLD)
    rank = comm.Get_rank()
    path = obs_file or icesee_kwargs.get("synthetic_obs_file")

    payload = None
    if rank == 0:
        if not path or not h5py.is_hdf5(path):
            raise FileNotFoundError(
                f"Cannot synchronize observation schedule: '{path}' is not a valid HDF5 file"
            )
        with h5py.File(path, "r") as f:
            if "obs_index" in f:
                obs_index = np.asarray(f["obs_index"][:], dtype=np.int64)
            elif "ind_m" in f:
                obs_index = np.asarray(f["ind_m"][:], dtype=np.int64)
            else:
                raise KeyError(f"'{path}' contains neither /obs_index nor /ind_m")

            if "obs_t" in f:
                obs_t = np.asarray(f["obs_t"][:], dtype=float)
            else:
                # Older compact files did not retain obs_t.  Recover it from
                # the configured model time vector without changing indices.
                model_t = np.asarray(icesee_kwargs.get("t", []), dtype=float)
                if model_t.size and np.all((obs_index >= 0) & (obs_index < model_t.size)):
                    obs_t = model_t[obs_index]
                else:
                    dt = float(icesee_kwargs.get("dt", 1.0))
                    obs_t = obs_index.astype(float) * dt
        if obs_t.size != obs_index.size:
            raise ValueError(
                f"Observation schedule mismatch in '{path}': "
                f"{obs_index.size} indices but {obs_t.size} times"
            )
        payload = (obs_index, obs_t)

    obs_index, obs_t = comm.bcast(payload, root=0)
    obs_index = np.asarray(obs_index, dtype=np.int64)
    obs_t = np.asarray(obs_t, dtype=float)
    m_obs = int(obs_index.size)
    mapping = {int(step): int(col) for col, step in enumerate(obs_index)}
    icesee_kwargs.update({
        "ind_m": obs_index,
        "obs_index": obs_index,
        "obs_t": obs_t,
        "tobserve": obs_index,
        "m_obs": m_obs,
        "number_obs_instants": m_obs,
        "obs_model_to_col": mapping,
    })
    return icesee_kwargs


def generate_synthetic_observations(**icesee_kwargs):
    """Generate synthetic observations for the ICESEE model.
    """

    # unpack icesee_kwargs
    model_module   = icesee_kwargs.get("model_module", None)
    comm_world     = icesee_kwargs.get("comm_world", MPI.COMM_WORLD)
    _synthetic_obs = icesee_kwargs.get("synthetic_obs_file")
    _true_nurged   = icesee_kwargs.get("true_nurged_file")
    color          = icesee_kwargs.get("color", 0)
    subcomm        = icesee_kwargs.get("subcomm", None)
    sub_rank       = icesee_kwargs.get("sub_rank", 0)
    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()


    if icesee_kwargs.get("generate_synthetic_obs", True):
        if icesee_kwargs["even_distribution"] or (icesee_kwargs["default_run"] and size_world <= icesee_kwargs["Nens"]):
            if rank_world == 0:
                # --- Synthetic Observations ---
                print("[ICESEE] Generating synthetic observations ...")
                with h5py.File(_true_nurged, "r") as f:
                    ensemble_true_state = f['true_state'][:]

                utils_funs = UtilsFunctions(
                    icesee_kwargs=icesee_kwargs,
                    ensemble=ensemble_true_state
                )
                icesee_kwargs.update({"statevec_true": ensemble_true_state})
                hu_obs, error_R, bed_masks, icesee_kwargs = utils_funs._create_synthetic_observations(**icesee_kwargs)

                # check if the best_mask_map is generated
                icesee_kwargs.update({"bed_mask_map": bed_masks})

                # observe or don't observe parameters.
                vecs, indx_map,_ = icesee_get_index(hu_obs, **icesee_kwargs)
                # check if icesee_kwargs['observe_params'] is empty
                if len(icesee_kwargs['observed_params']) == 0:
                    for key in icesee_kwargs['params_vec']:
                        hu_obs[indx_map[key],:] = 0.0
                        error_R[:,indx_map[key]] = 0.0
                else:
                    for key in icesee_kwargs['params_vec']:
                        if key not in icesee_kwargs['observed_params']:
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
                    f.create_dataset("bed_snap_cols", data=np.asarray(icesee_kwargs["bed_snap_cols"], dtype=int))
                    obs_index = np.asarray(icesee_kwargs["ind_m"], dtype=int)
                    obs_t = np.asarray(icesee_kwargs["obs_t"], dtype=float)
                    f.create_dataset("ind_m", data=obs_index)
                    f.create_dataset("obs_t", data=obs_t)
                    # Plotting-facing aliases keep all observation data and
                    # metadata together in synthetic_obs.h5.  The original
                    # names remain for backward compatibility with existing
                    # analysis and restart code.
                    f.create_dataset("obs_index", data=obs_index)
                    f.create_dataset(
                        "obs_max_time",
                        data=np.asarray(
                            [np.max(obs_t) if obs_t.size else 0.0], dtype=float
                        ),
                    )

                    # obs_model_to_col is a dict -> store as parallel arrays
                    m = icesee_kwargs.get("obs_model_to_col", {})
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
                # hu_obs = np.empty((icesee_kwargs["nd"],icesee_kwargs["number_obs_instants"]),dtype=np.float64)
                # error_R = np.empty((icesee_kwargs["number_obs_instants"], icesee_kwargs["nd"]),dtype=np.float64)

            if icesee_kwargs["even_distribution"]:
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

            if icesee_kwargs["default_run"] and size_world > icesee_kwargs["Nens"]:
                subcomm.Barrier()
                # comm_world.Bcast(hu_obs, root=0)
                if sub_rank == 0:
                    utils_funs = UtilsFunctions(
                        icesee_kwargs=icesee_kwargs,
                        ensemble=ensemble_true_state
                    )
                    icesee_kwargs.update({"statevec_true": ensemble_true_state})
                    hu_obs, error_R, bed_mask_map, icesee_kwargs = utils_funs._create_synthetic_observations(**icesee_kwargs)
                    icesee_kwargs.update({"bed_mask_map": bed_mask_map})

                    # observe or don't observe parameters.
                    vecs, indx_map,_ = icesee_get_index(hu_obs, **icesee_kwargs)
                    # check if icesee_kwargs['observe_params'] is empty
                    if len(icesee_kwargs['observed_params']) == 0:
                        for key in icesee_kwargs['params_vec']:
                            hu_obs[indx_map[key],:] = 0.0
                            error_R[:,indx_map[key]] = 0.0
                    else:
                        for key in icesee_kwargs['params_vec']:
                            if key not in icesee_kwargs['observed_params']:
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
            elif icesee_kwargs["sequential_run"]:
                comm_world.Barrier()
                # g_shape = icesee_kwargs['dim_list'][rank_world]
                # utils_funs = UtilsFunctions(icesee_kwargs, ensemble_true_state)
                # icesee_kwargs.update({"statevec_true": ensemble_true_state})
                # hu_obs = utils_funs._create_synthetic_observations(**icesee_kwargs)
                # # gather from every rank to rank 0
                # gathered_obs = comm_world.gather(hu_obs[:g_shape,:], root=0)
                # if rank_world == 0:
                #     print(f"[ICESEE] {[arr.shape for arr in gathered_obs]}")
                #     hu_obs = np.vstack(gathered_obs)
                # else:
                #     hu_obs = np.empty((icesee_kwargs["global_shape"],icesee_kwargs["number_obs_instants"]),dtype=np.float64)

                # comm_world.Bcast(hu_obs, root=0)
                if rank_world == 0:
                    utils_funs = UtilsFunctions(
                        icesee_kwargs=icesee_kwargs,
                        ensemble=ensemble_true_state
                    )
                    icesee_kwargs.update({"statevec_true": ensemble_true_state})
                    hu_obs, error_R, bed_mask_map, icesee_kwargs = utils_funs._create_synthetic_observations(**icesee_kwargs)
                    icesee_kwargs.update({"bed_mask_map": bed_mask_map})
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

    return icesee_kwargs
