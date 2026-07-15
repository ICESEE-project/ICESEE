# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-03-06
# @description: computes the X5 matrix for the EnKF
#               - the new formulation is based on the paper by Geir Evensen: The Ensemble Kalman Filter: Theoretical Formulation And Practical Implementation
#               - this formulation supports our need for mpi parallelization and no need for localizations
# =============================================================================

import gc
import os
import copy
import h5py
import numpy as np
import bigmpi4py as BM
from scipy.stats import multivariate_normal, beta
from mpi4py import MPI

from ICESEE.src.utils.tools import icesee_get_index
# from ICESEE.src.run_model_da._parallel_i_o import parallel_write_full_ensemble_from_root
                                               
from ICESEE.src.run_model_da._error_generation import generate_enkf_field
from ICESEE.src.utils.utils import UtilsFunctions

from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager
# rank_seed, rng = ParallelManager().initialize_seed(MPI.COMM_WORLD)
from ICESEE.src.parallelization._parallel_i_o import (write_ensemble_member_direct, write_ensemble_member_direct_h5, open_ensemble_file)

def parallel_forecast_step_default_run(**model_kwargs):
    import h5py
    import numpy as np
    from mpi4py import MPI

    params = model_kwargs.get("params")
    Nens = params["Nens"]
    nd = model_kwargs.get("nd", params.get("nd"))

    rounds = model_kwargs.get("rounds")
    color = model_kwargs.get("color", 0)
    subcomm_size_min = model_kwargs.get("subcomm_size_min", 1)

    subcomm = model_kwargs.get("subcomm", MPI.COMM_SELF)
    comm_world = model_kwargs.get("comm_world", MPI.COMM_WORLD)

    rank_world = comm_world.Get_rank()
    sub_rank = subcomm.Get_rank()

    model_module = model_kwargs["model_module"]
    k = model_kwargs.get("k", 0)

    _modelrun_datasets = model_kwargs.get("_modelrun_datasets", "_modelrun_datasets")
    input_file = f"{_modelrun_datasets}/icesee_ensemble_data.h5"

    partitioned_io = model_kwargs.get("partitioned_io_flag", False)  # NEW

    alpha = model_kwargs.get("alpha", 0.0)
    rho = model_kwargs.get("rho", 1.0)
    dt = model_kwargs.get("dt", 1.0)
    Lx = model_kwargs.get("Lx", 1.0)
    Ly = model_kwargs.get("Ly", 1.0)

    noise = model_kwargs.get("noise", None)
    if noise is None:
        raise ValueError("model_kwargs must contain `noise`.")

    time_forecast_ensemble_generation = model_kwargs.get("time_forecast_ensemble_generation", 0.0)
    time_forecast_noise_generation = model_kwargs.get("time_forecast_noise_generation", 0.0)
    time_forecast_ensemble_mean_generation = model_kwargs.get("time_forecast_ensemble_mean_generation", 0.0)
    time_forecast_file_writing = model_kwargs.get("time_forecast_file_writing", 0.0)

    root_comm = comm_world.Split(
        color=0 if sub_rank == 0 else MPI.UNDEFINED,
        key=rank_world,
    )
    is_root_leader = root_comm is not MPI.COMM_NULL
    root_rank = root_comm.Get_rank() if is_root_leader else None
    root_is_world_root = is_root_leader and rank_world == 0

    vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)

    def _add_process_noise(ensemble_vec, ens_id, local_kwargs):
        nonlocal noise
        nonlocal time_forecast_noise_generation

        if model_kwargs["joint_estimation"] or params["localization_flag"]:
            hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
        else:
            hdim = ensemble_vec.shape[0] // params["num_state_vars"]

        state_block_size = hdim * params["num_state_vars"]
        _t_noise = MPI.Wtime()

        noise_all = []
        q0 = []
        for ii, sig in enumerate(params["sig_Q"]):
            if ii < params["num_state_vars"]:
                noise_kwargs = dict(local_kwargs)
                noise_kwargs.update({"ii_sig": ii, "Lx_dim": np.sqrt(Lx * Ly), "noise_dim": hdim, "num_vars": params["total_state_param_vars"]})
                W = generate_enkf_field(**noise_kwargs)
                prev_noise = noise[ii * hdim: (ii + 1) * hdim]
                noise_i = alpha * prev_noise + np.sqrt(1.0 - alpha**2) * W
                q0.append(noise_i)
                noise_all.append(np.sqrt(dt) * sig * rho * noise_i)

        if noise_all:
            noise_update = np.concatenate(noise_all, axis=0)
            ensemble_vec[:state_block_size] += noise_update[:state_block_size]
            noise = np.concatenate(q0, axis=0)

        time_forecast_noise_generation += MPI.Wtime() - _t_noise
        return ensemble_vec

    def _gather_root_vectors(local_vec, ens_id):
        if not is_root_leader:
            return None, None
        active = local_vec is not None
        send_id = np.array([ens_id if active else -1], dtype=np.int64)
        recv_ids = None
        if root_rank == 0:
            recv_ids = np.empty(root_comm.Get_size(), dtype=np.int64)
        root_comm.Gather(send_id, recv_ids, root=0)

        send_vec = np.asarray(local_vec, dtype=np.float64) if active else np.empty(nd, dtype=np.float64)
        recv_mat = None
        if root_rank == 0:
            recv_mat = np.empty((root_comm.Get_size(), nd), dtype=np.float64)
        root_comm.Gather(send_vec, recv_mat, root=0)
        return recv_mat, recv_ids

    # ---- allocate full ensemble ONLY when not partitioned ----          # NEW
    if not partitioned_io:
        ensemble_vec = np.empty((nd, Nens), dtype=np.float64)
    else:
        ensemble_vec = None                                              # NEW

    _t_forecast = MPI.Wtime()

    # NEW — open the file ONCE for the whole round loop when partitioned
    ens_file = None
    ens_dset = None
    if partitioned_io:
        ens_file = open_ensemble_file(input_file, nd, Nens, model_kwargs.get("nt", params["nt"]), comm_world)
        ens_dset = ens_file["ensemble"]

    if Nens >= comm_world.Get_size():
        for round_id in range(rounds):
            ens_id = color + round_id * subcomm_size_min
            local_vec = None

            if ens_id < Nens:
                local_kwargs = dict(model_kwargs)
                local_kwargs.update({"ens_id": ens_id, "comm": subcomm})
                subcomm.Barrier()

                with h5py.File(input_file, "r") as f:
                    ensemble_local = f["ensemble"][:, ens_id, k].astype(np.float64, copy=True)

                updated_state = model_module.forecast_step_single(ensemble=ensemble_local, **local_kwargs)
                for key, value in updated_state.items():
                    ensemble_local[indx_map[key]] = value

                obs_index = model_kwargs["obs_index"]
                km = model_kwargs.get("km")
                if (km < params["number_obs_instants"]) and (k == obs_index[km]):
                    ensemble_local = _add_process_noise(ensemble_local, ens_id, local_kwargs)

                if sub_rank == 0:
                    local_vec = ensemble_local

            # ============================================== NEW BRANCH
            if partitioned_io:
                            member_to_write = local_vec if sub_rank == 0 else None
                            this_ens_id = ens_id if (sub_rank == 0 and ens_id < Nens) else None
                            write_ensemble_member_direct_h5(ens_dset, k + 1, this_ens_id, member_to_write, nd, comm_world)
            else:
                # ---- EXISTING behavior, unchanged ----
                recv_mat, recv_ids = _gather_root_vectors(local_vec, ens_id)
                if root_is_world_root and recv_mat is not None:
                    for row, eid in enumerate(recv_ids):
                        eid = int(eid)
                        if 0 <= eid < Nens:
                            ensemble_vec[:, eid] = recv_mat[row, :]
            # ============================================== END NEW

    else:
        ens_id = color
        local_vec = None

        if ens_id < Nens:
            local_kwargs = dict(model_kwargs)
            local_kwargs.update({"ens_id": ens_id, "comm": subcomm})
            subcomm.Barrier()

            with h5py.File(input_file, "r") as f:
                ensemble_local = f["ensemble"][:, ens_id, k].astype(np.float64, copy=True)

            updated_state = model_module.forecast_step_single(ensemble=ensemble_local, **local_kwargs)
            for key, value in updated_state.items():
                ensemble_local[indx_map[key]] = value

            ensemble_local = _add_process_noise(ensemble_local, ens_id, local_kwargs)

            if sub_rank == 0:
                local_vec = ensemble_local

        # ============================================== NEW BRANCH
        if partitioned_io:
            member_to_write = local_vec if sub_rank == 0 else None
            this_ens_id = ens_id if (sub_rank == 0 and ens_id < Nens) else None
            write_ensemble_member_direct_h5(ens_dset, k + 1, this_ens_id, member_to_write, nd, comm_world)
        else:
            # ---- EXISTING behavior, unchanged ----
            recv_mat, recv_ids = _gather_root_vectors(local_vec, ens_id)
            if root_is_world_root and recv_mat is not None:
                for row, eid in enumerate(recv_ids):
                    eid = int(eid)
                    if 0 <= eid < Nens:
                        ensemble_vec[:, eid] = recv_mat[row, :]
        # ============================================== END NEW

    # NEW — single Barrier + close after all rounds, instead of per-round
    if partitioned_io:
        comm_world.Barrier()
        ens_file.close()

    time_forecast_ensemble_generation += MPI.Wtime() - _t_forecast

    if rank_world == 0:
        shape_ens = np.array((nd, Nens), dtype=np.int32)
    else:
        shape_ens = np.empty(2, dtype=np.int32)
    shape_ens = comm_world.bcast(shape_ens, root=0)

    _t_mean = MPI.Wtime()
    if not partitioned_io:
        # ---- EXISTING behavior, unchanged ----
        ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], Nens, comm_world, root=0)
    else:
        ens_mean = None  # NEW — not needed downstream in the partitioned path
    time_forecast_ensemble_mean_generation += MPI.Wtime() - _t_mean

    model_kwargs.update({
        "time_forecast_ensemble_generation": time_forecast_ensemble_generation,
        "time_forecast_noise_generation": time_forecast_noise_generation,
        "time_forecast_ensemble_mean_generation": time_forecast_ensemble_mean_generation,
        "time_forecast_file_writing": time_forecast_file_writing,
        "shape_ens": shape_ens,
        "noise": noise,
        "_forecast_h5_path": input_file,   # NEW — consumed by EnKF_X5's partitioned branch
    })

    if is_root_leader:
        root_comm.Free()

    return model_kwargs, ensemble_vec, shape_ens, ens_mean

def parallel_forecast_step_default_full_parallel_run(**model_kwargs):
    """
    Full-parallel forecast step for ICESEE.

    This version is file-backed:
    - each subcommunicator advances one ensemble member at a time;
    - only sub_rank == 0 writes the completed ensemble vector;
    - no global ensemble matrix is gathered in memory;
    - both Nens >= size_world and Nens < size_world are handled consistently.
    """

    import numpy as np
    from mpi4py import MPI

    params = model_kwargs.get("params")
    Nens = params["Nens"]

    comm_world = model_kwargs.get("comm_world", MPI.COMM_WORLD)
    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()

    subcomm = model_kwargs.get("subcomm", MPI.COMM_SELF)
    sub_rank = subcomm.Get_rank()

    color = model_kwargs.get("color", 0)
    rounds = model_kwargs.get("rounds", 1)
    subcomm_size_min = model_kwargs.get("subcomm_size_min", 1)

    model_module = model_kwargs.get("model_module")
    enkf_parallel_io = model_kwargs.get("enkf_parallel_io")

    if enkf_parallel_io is None:
        raise ValueError("parallel_forecast_step_default_full_parallel_run requires enkf_parallel_io.")

    nd = model_kwargs.get("nd", params.get("nd"))
    nt = model_kwargs.get("nt", params.get("nt"))
    k = model_kwargs.get("k", 0)

    alpha = model_kwargs.get("alpha", 0.0)
    rho = model_kwargs.get("rho", 1.0)
    dt = model_kwargs.get("dt", 1.0)
    Lx = model_kwargs.get("Lx", 1.0)
    Ly = model_kwargs.get("Ly", 1.0)

    noise = model_kwargs.get("noise", None)
    if noise is None:
        raise ValueError("model_kwargs must contain `noise`.")

    time_forecast_ensemble_generation = model_kwargs.get(
        "time_forecast_ensemble_generation", 0.0
    )
    time_forecast_noise_generation = model_kwargs.get(
        "time_forecast_noise_generation", 0.0
    )
    time_forecast_ensemble_mean_generation = model_kwargs.get(
        "time_forecast_ensemble_mean_generation", 0.0
    )
    time_forecast_file_writing = model_kwargs.get(
        "time_forecast_file_writing", 0.0
    )

    vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)

    def _add_process_noise(ensemble_vec, ens_id, local_kwargs):
        nonlocal noise
        nonlocal time_forecast_noise_generation

        if model_kwargs["joint_estimation"] or params["localization_flag"]:
            hdim = ensemble_vec.shape[0] // params["total_state_param_vars"]
        else:
            hdim = ensemble_vec.shape[0] // params["num_state_vars"]

        state_block_size = hdim * params["num_state_vars"]

        _t_noise = MPI.Wtime()

        noise_all = []
        q0 = []

        for ii, sig in enumerate(params["sig_Q"]):
            # IMPORTANT: strict <, not <=
            if ii < params["num_state_vars"]:
                noise_kwargs = dict(local_kwargs)
                noise_kwargs.update(
                    {
                        "ens_id": ens_id,
                        "ii_sig": ii,
                        "Lx_dim": np.sqrt(Lx * Ly),
                        "noise_dim": hdim,
                        "num_vars": params["total_state_param_vars"],
                    }
                )

                W = generate_enkf_field(**noise_kwargs)

                prev_noise = noise[ii * hdim : (ii + 1) * hdim]
                noise_i = alpha * prev_noise + np.sqrt(1.0 - alpha**2) * W

                q0.append(noise_i)
                noise_all.append(np.sqrt(dt) * sig * rho * noise_i)

        if noise_all:
            noise_update = np.concatenate(noise_all, axis=0)
            ensemble_vec[:state_block_size] += noise_update[:state_block_size]
            noise = np.concatenate(q0, axis=0)

        time_forecast_noise_generation += MPI.Wtime() - _t_noise

        return ensemble_vec

    def _run_one_ensemble(ens_id):
        nonlocal time_forecast_file_writing

        if ens_id < 0 or ens_id >= Nens:
            return

        local_kwargs = dict(model_kwargs)
        local_kwargs.update(
            {
                "ens_id": ens_id,
                "comm": subcomm,
            }
        )

        subcomm.Barrier()

        _t_read = MPI.Wtime()
        ensemble_vec = enkf_parallel_io.read_forecast(k, ens_id).astype(
            np.float64,
            copy=True,
        )
        time_read = MPI.Wtime() - _t_read

        updated_state = model_module.forecast_step_single(
            ensemble=ensemble_vec,
            **local_kwargs,
        )

        for key, value in updated_state.items():
            ensemble_vec[indx_map[key]] = value

        ensemble_vec = _add_process_noise(
            ensemble_vec,
            ens_id,
            local_kwargs,
        )

        # To avoid duplicate writes, only the subcommunicator root writes
        # the completed ensemble member.
        if sub_rank == 0:
            _t_write = MPI.Wtime()
            write_k = k + 1 if k < nt - 1 else k
            enkf_parallel_io.write_forecast(write_k, ensemble_vec, ens_id)
            time_forecast_file_writing += MPI.Wtime() - _t_write + time_read

        subcomm.Barrier()

    _t_forecast = MPI.Wtime()

    # ------------------------------------------------------------------
    # Case 2: Nens >= size_world
    # each subcomm processes multiple ensemble members over rounds.
    # ------------------------------------------------------------------
    if Nens >= size_world:
        for round_id in range(rounds):
            ens_id = color + round_id * subcomm_size_min
            _run_one_ensemble(ens_id)

    # ------------------------------------------------------------------
    # Case 3: Nens < size_world
    # one subcommunicator per ensemble member.
    # ------------------------------------------------------------------
    else:
        ens_id = color
        _run_one_ensemble(ens_id)

    time_forecast_ensemble_generation += MPI.Wtime() - _t_forecast

    # Shape is known globally; no need to gather ensemble.
    shape_ens = np.array((nd, Nens), dtype=np.int32)
    shape_ens = comm_world.bcast(shape_ens, root=0)

    # ------------------------------------------------------------------
    # Compute forecast mean only at observation time.
    # ------------------------------------------------------------------
    _t_mean = MPI.Wtime()

    km = model_kwargs.get("km", 0)
    obs_index = model_kwargs["obs_index"]

    if (km < params["number_obs_instants"]) and (k == obs_index[km]):
        write_k = k + 1 if k < nt - 1 else k
        enkf_parallel_io.compute_forecast_mean_chunked_v2(
            write_k,
            flag="initial",
        )

    time_forecast_ensemble_mean_generation += MPI.Wtime() - _t_mean

    model_kwargs.update(
        {
            "time_forecast_ensemble_generation": time_forecast_ensemble_generation,
            "time_forecast_noise_generation": time_forecast_noise_generation,
            "time_forecast_ensemble_mean_generation": time_forecast_ensemble_mean_generation,
            "time_forecast_file_writing": time_forecast_file_writing,
            "shape_ens": shape_ens,
            "noise": noise,
        }
    )

    return model_kwargs



def parallel_forecast_step_even_distribution_run(**model_kwargs):
    """ Parallel run of the forecast step for each ensemble member.
        This function is designed to be used in a distributed environment where each rank processes a single ensemble member.
        It assumes that the number of ensemble members (Nens) is divisible by the size of the world communicator (size_world).
    """

    # unpack model_kwargs
    params = model_kwargs.get("params")
    model_module = model_kwargs.get("model_module")
    comm_world = model_kwargs.get("comm_world", MPI.COMM_WORLD)
    rank_world = comm_world.Get_rank()
    Nens = params.get("Nens", 1)  # Number of ensemble members
    nd = params.get("nd", 1)  # Dimension of the state vector
    Q_err = model_kwargs.get("Q_err", np.eye(nd))  # Error covariance matrix
    state_block_size = model_kwargs.get("state_block_size", nd)  # Size of the state block
    size_world = comm_world.Get_size()  # Total number of ranks in the world communicato
    ensemble_vec = model_kwargs.get("ensemble_vec", np.zeros((nd, Nens), dtype=np.float64))  # Initialize ensemble vector
    ensemble_vec_mean = model_kwargs.get("ensemble_vec_mean", np.zeros((nd, params.get("nt", params["nt"]) + 1), dtype=np.float64))  # Initialize ensemble mean vector
    shape_ens = np.array(ensemble_vec.shape, dtype=np.int32)  # Shape of the ensemble vector
    ensemble_local = model_kwargs.get("ensemble_local", np.zeros((nd, Nens), dtype=np.float64))  # Local ensemble vector
    k = model_kwargs.get("k", 0)  # Time step index, default to 0 if not provided

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

    return ensemble_vec, ensemble_vec_mean, shape_ens

def parallel_forecast_step_squential_run(**model_kwargs):
    """ Squential run of the forecast step for each ensemble member.
        This function is designed to be used in a distributed environment where each rank processes a single ensemble member.
        #TODO: still under development, not fully tested.
    """

    # unpack model_kwargs
    params = model_kwargs.get("params")
    model_module = model_kwargs.get("model_module")
    comm_world = model_kwargs.get("comm_world", MPI.COMM_WORLD)
    rank_world = comm_world.Get_rank()
    Nens = params.get("Nens", 1)  # Number of ensemble members
    nd = params.get("nd", 1)  # Dimension of the state vector
    Q_err = model_kwargs.get("Q_err", np.eye(nd))  # Error covariance matrix
    state_block_size = model_kwargs.get("state_block_size", nd)  # Size of the state block
    ensemble_vec = model_kwargs.get("ensemble_vec", np.zeros((nd, Nens), dtype=np.float64))  # Initialize ensemble vector
    ensemble_vec_mean = model_kwargs.get("ensemble_vec_mean", np.zeros((nd, params.get("nt", params["nt"]) + 1), dtype=np.float64))  # Initialize ensemble mean vector
    shape_ens = np.array(ensemble_vec.shape, dtype=np.int32)  # Shape of the ensemble vector
    ensemble_local = model_kwargs.get("ensemble_local", np.zeros((nd, Nens), dtype=np.float64))  # Local ensemble vector
    k = model_kwargs.get("k", 0)  # Time step index, default to 0 if not provided

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

    
    return model_kwargs, ensemble_vec, ensemble_vec_mean, shape_ens