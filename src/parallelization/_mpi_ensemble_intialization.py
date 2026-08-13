# ==============================================================================
# @des: This file contains run functions for the ICESEE model to initialize the ensemble.
# @date: 2025-07-30
# @author: Brian Kyanjo
# ==============================================================================

# --- import necessary libraries ---
import numpy as np
import h5py
import gc
import zarr
import os
from mpi4py import MPI

from ICESEE.src.utils.tools import icesee_get_index, env_flag
from ICESEE.src.run_model_da._error_generation import compute_Q_err_random_fields, \
                              compute_noise_random_fields, \
                              generate_pseudo_random_field_1d, \
                              generate_pseudo_random_field_2D, \
                              generate_enkf_field

from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager
# rank_seed, rng = ParallelManager().initialize_seed(MPI.COMM_WORLD)

from ICESEE.src.parallelization._parallel_i_o import (
    parallel_write_full_ensemble_from_root,
    parallel_write_full_ensemble_from_root_full_parallel_run,
    write_ensemble_member_direct,
    compute_and_apply_inflation_partitioned,
)


def _assemble_initialized_members(gathered_rank_members, nens):
    """Reassemble variable-length MPI member lists in global member order."""
    members_by_id = {}
    for rank_members in gathered_rank_members:
        for member_id, member_vector in rank_members:
            member_id = int(member_id)
            if member_id in members_by_id:
                raise RuntimeError(
                    "Duplicate ensemble member produced during initialization: "
                    f"member {member_id}."
                )
            members_by_id[member_id] = np.asarray(member_vector).reshape(-1)

    expected_ids = set(range(nens))
    actual_ids = set(members_by_id)
    missing_ids = sorted(expected_ids - actual_ids)
    unexpected_ids = sorted(actual_ids - expected_ids)
    if missing_ids or unexpected_ids:
        raise RuntimeError(
            "MPI ensemble initialization did not produce the configured "
            "member set. "
            f"Missing ids={missing_ids}; unexpected ids={unexpected_ids}."
        )

    member_lengths = {members_by_id[i].size for i in range(nens)}
    if len(member_lengths) != 1:
        raise RuntimeError(
            "MPI ensemble initialization produced inconsistent state-vector "
            f"lengths: {sorted(member_lengths)}."
        )

    return np.column_stack([members_by_id[i] for i in range(nens)])

def ensemble_initialization(**icesee_kwargs):
    """Initialize the ensemble for the ICESEE model.
    """

    # unpack icesee_kwargs
    model_module   = icesee_kwargs.get("model_module", None)
    comm_world     = icesee_kwargs.get("comm_world", MPI.COMM_WORLD)
    subcomm        = icesee_kwargs.get("subcomm", None)
    color          = icesee_kwargs.get("color", 0)
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
    rounds         = icesee_kwargs.get("rounds", 1)
    subcomm_size_min   = icesee_kwargs.get("subcomm_size_min", 1)
    rng           = icesee_kwargs.get("rng", np.random.default_rng())
    rank_seed = icesee_kwargs.get("rank_seed", 0)
    alpha = icesee_kwargs.get("initial_spread_factor")

    partitioned_io = icesee_kwargs.get("partitioned_io_flag", False)  # NEW

    sub_rank  = subcomm.Get_rank()
    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()

    time_init_noise_generation = 0.0
    time_init_file_writing     = 0.0
    time_init_ensemble_mean_computation = 0.0

    observed_vars = icesee_kwargs.get("observed_vars", [])
    observed_params = icesee_kwargs.get("observed_params", [])

    all_observed = list(observed_vars) + list(observed_params)

    icesee_kwargs["observed_vars_params"] = all_observed
    icesee_kwargs["all_observed"] = all_observed
    icesee_kwargs["all_observed"] = all_observed
    icesee_kwargs["nd_observed"] = len(all_observed) * (icesee_kwargs["nd"] // icesee_kwargs["total_state_param_vars"])

    if icesee_kwargs["even_distribution"] or (icesee_kwargs["default_run"] and size_world <= icesee_kwargs["Nens"]):
        if icesee_kwargs["default_run"] and size_world <= icesee_kwargs["Nens"] and not (icesee_kwargs.get("sequential_ensemble_initialization", False)):
            if rank_world == 0:
                print("[ICESEE] Initializing the ensemble ...")

            Nens = icesee_kwargs["Nens"]
            icesee_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})
            icesee_kwargs.update({"statevec_ens":np.zeros([icesee_kwargs["nd"], icesee_kwargs["Nens"]])})

            vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)
            ensemble_vec = np.zeros_like(icesee_kwargs["statevec_ens"])

            if icesee_kwargs["joint_estimation"] or icesee_kwargs["localization_flag"]:
                hdim = ensemble_vec.shape[0] // icesee_kwargs["total_state_param_vars"]
            else:
                hdim = ensemble_vec.shape[0] // icesee_kwargs["num_state_vars"]
            state_block_size = hdim * icesee_kwargs["num_state_vars"]

            ens_list_init = []

            for round_id in range(rounds):
                ensemble_id = color + (round_id * subcomm_size_min)
                icesee_kwargs.update({'ens_id': ensemble_id})

                if ensemble_id < Nens:
                    subcomm.Barrier()
                    ens = ensemble_id

                    data = model_module.initialize_ensemble(ens, **icesee_kwargs)
                    for key, value in data.items():
                        ensemble_vec[indx_map[key], ens] = value

                    _time_init_noise_generation = MPI.Wtime()
                    icesee_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":icesee_kwargs["total_state_param_vars"]})
                    noise = generate_enkf_field(**icesee_kwargs)
                    time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation

                    noise_all = []
                    q0 = []
                    dt = icesee_kwargs.get("dt", 1.0)
                    min_tau = 200
                    max_tau = 500
                    dt  = icesee_kwargs.get("dt",icesee_kwargs["dt"])
                    tau = max(max_tau,max(min_tau, dt))

                    alpha_ar1 = 1 - dt/tau
                    if alpha_ar1 <= 0 or alpha_ar1 > 1:
                        alpha_ar1 = 0.5

                    n = icesee_kwargs.get("nt",icesee_kwargs["nt"])
                    rho = np.sqrt((1/dt)*((1-alpha_ar1)**2)*(1/(n - (2*alpha_ar1) - (n*alpha_ar1**2) + (2*alpha_ar1**(n+1)))))
                    for ii, sig in enumerate(icesee_kwargs["sig_Q"]):
                        if ii <= icesee_kwargs["total_state_param_vars"]:
                            icesee_kwargs.update({"ii_sig": ii, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":icesee_kwargs["total_state_param_vars"]})
                            W = generate_enkf_field(**icesee_kwargs)
                            noise_ = alpha_ar1*noise[ii*hdim:(ii+1)*hdim] + np.sqrt(1 - alpha_ar1**2)*W
                            q0.append(noise_)
                            Z = np.sqrt(dt)*sig*rho*noise_
                            noise_all.append(Z)
                    noise_ = np.concatenate(noise_all, axis=0)

                    for ii, sig in enumerate(icesee_kwargs["sig_Q"]):
                        if ii <= icesee_kwargs["total_state_param_vars"]:
                            start_idx = ii * hdim
                            end_idx = start_idx + hdim
                            ensemble_vec[start_idx:end_idx, ens] += noise[start_idx:end_idx] * sig

                    icesee_kwargs.update({"noise": noise})
                    del noise

                    # ============================================== NEW BRANCH
                    if partitioned_io:
                        member_to_write = ensemble_vec[:, ens] if sub_rank == 0 else None
                        this_ens_id = ensemble_id if sub_rank == 0 else None
                        write_ensemble_member_direct(
                            f"{icesee_kwargs.get('data_path')}/icesee_ensemble_data.h5",
                            0, this_ens_id, member_to_write,
                            icesee_kwargs["nd"], Nens, icesee_kwargs.get("nt", icesee_kwargs["nt"]),
                            comm_world
                        )
                    else:
                        # ---- EXISTING behavior, unchanged ----
                        gathered_ensemble = subcomm.gather(ensemble_vec[:, ens], root=0)
                        if sub_rank == 0:
                            gathered_ensemble = np.concatenate(gathered_ensemble, axis=0)
                            # Retain the global member id.  The final MPI round
                            # can be only partially populated when Nens is not
                            # divisible by the number of model groups, so the
                            # result must not be reconstructed by list position.
                            ens_list_init.append((ensemble_id, gathered_ensemble))
                        del gathered_ensemble
                    # ============================================== END NEW

            if not partitioned_io:
                # Python-object gather deliberately supports different list
                # lengths on different ranks.  The former fixed-shape numeric
                # Gather used the root rank's list length for every rank; on a
                # partial final round that could manufacture padded members
                # (for example 64 columns for Nens=60).
                gathered_ensemble_global = comm_world.gather(ens_list_init, root=0)

            # ============================================== NEW BRANCH
            if partitioned_io:
                del ens_list_init; gc.collect()
                comm_world.Barrier()
                compute_and_apply_inflation_partitioned(
                    f"{icesee_kwargs.get('data_path')}/icesee_ensemble_data.h5",
                    icesee_kwargs["nd"], Nens, alpha, comm_world, timestep=0
                )
                ensemble_vec = None
                shape_ens = np.array([icesee_kwargs["nd"], Nens], dtype=np.int32)
            else:
                # ---- EXISTING reassembly + inflation, unchanged ----
                del ens_list_init; gc.collect()
                assembly_error = None
                if rank_world == 0:
                    try:
                        ensemble_vec_final = _assemble_initialized_members(
                            gathered_ensemble_global, Nens
                        )
                        shape_ens = np.array(
                            ensemble_vec_final.shape, dtype=np.int32
                        )
                        ensemble_vec = ensemble_vec_final

                        mean_params = np.mean(ensemble_vec, axis=1)
                        pertubations = ensemble_vec - mean_params.reshape(-1,1)
                        inflated_pertubations = pertubations * alpha
                        ensemble_vec = mean_params.reshape(-1,1) + inflated_pertubations
                        del ensemble_vec_final
                    except Exception as exc:
                        assembly_error = f"{type(exc).__name__}: {exc}"
                        shape_ens = np.empty(2, dtype=np.int32)
                else:
                    shape_ens = np.empty(2, dtype=np.int32)

                assembly_error = comm_world.bcast(assembly_error, root=0)
                if assembly_error is not None:
                    raise RuntimeError(
                        "Collective ensemble initialization failed: "
                        f"{assembly_error}"
                    )
                shape_ens = comm_world.bcast(shape_ens, root=0)
            # ============================================== END NEW

        else:
            # ---- EXISTING sequential-ensemble-initialization branch, fully unchanged ----
            if rank_world == 0:
                print("[ICESEE] Initializing the ensemble ...")
                icesee_kwargs.update({'ens_id': rank_world})
                if icesee_kwargs["even_distribution"]:
                    icesee_kwargs.update({'rank': rank_world, 'color': color, 'comm': comm_world})
                else:
                    icesee_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})

                icesee_kwargs.update({"statevec_ens":np.zeros([icesee_kwargs["nd"], icesee_kwargs["Nens"]])})
                vecs, indx_map, dim_per_proc = icesee_get_index(icesee_kwargs["statevec_ens"], **icesee_kwargs)
                ensemble_vec = np.zeros_like(icesee_kwargs["statevec_ens"])

                if icesee_kwargs["joint_estimation"] or icesee_kwargs["localization_flag"]:
                    hdim = ensemble_vec.shape[0] // icesee_kwargs["total_state_param_vars"]
                else:
                    hdim = ensemble_vec.shape[0] // icesee_kwargs["num_state_vars"]
                state_block_size = hdim * icesee_kwargs["num_state_vars"]

                for ens in range(icesee_kwargs["Nens"]):
                    data = model_module.initialize_ensemble(ens,**icesee_kwargs)
                    for key, value in data.items():
                        ensemble_vec[indx_map[key],ens] = value

                    _time_init_noise_generation = MPI.Wtime()
                    N_size = icesee_kwargs["total_state_param_vars"] * hdim
                    icesee_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":icesee_kwargs["total_state_param_vars"]})
                    noise = generate_enkf_field(**icesee_kwargs)
                    time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation

                    for ii, sig in enumerate(icesee_kwargs["sig_Q"]):
                        if ii <= icesee_kwargs["total_state_param_vars"]:
                            start_idx = ii * hdim
                            end_idx = start_idx + hdim
                            ensemble_vec[start_idx:end_idx, ens] += noise[start_idx:end_idx] * sig

                    icesee_kwargs.update({"noise": noise})

                shape_ens = np.array(ensemble_vec.shape,dtype=np.int32)

                mean_params = np.mean(ensemble_vec, axis=1)
                pertubations = ensemble_vec - mean_params.reshape(-1,1)
                inflated_pertubations = pertubations * alpha
                ensemble_vec = mean_params.reshape(-1,1) + inflated_pertubations

            else:
                ensemble_vec = np.empty((icesee_kwargs["nd"],icesee_kwargs["Nens"]),dtype=np.float64)
                shape_ens = np.empty(2,dtype=np.int32)

        comm_world.Barrier()

        # now reset the model_nprocs
        if rank_world == 0:
            diff = total_cores - base_total_procs
            if diff >= 0:
                min_model_nprocs = max(model_nprocs-1, 1)
                if icesee_kwargs.get('ICESEE_PERFORMANCE_TEST') or env_flag("ICESEE_PERFORMANCE_TEST", default=False):
                    model_nprocs = model_nprocs
                else:
                    model_nprocs = max(min_model_nprocs, model_nprocs + (diff // size_world))
            else:
                model_nprocs = model_nprocs

        model_nprocs = comm_world.bcast(model_nprocs, root=0)
        icesee_kwargs.update({'model_nprocs': model_nprocs})

        if icesee_kwargs["even_distribution"]:
            comm_world.Bcast(ensemble_vec, root=0)
            ensemble_bg = np.empty((icesee_kwargs["nd"],icesee_kwargs.get("nt",icesee_kwargs["nt"])+1),dtype=np.float64)
            ensemble_vec_mean = np.empty((icesee_kwargs["nd"],icesee_kwargs.get("nt",icesee_kwargs["nt"])+1),dtype=np.float64)
            ensemble_vec_full = np.empty((icesee_kwargs["nd"],icesee_kwargs["Nens"],icesee_kwargs.get("nt",icesee_kwargs["nt"])+1),dtype=np.float64)
            ensemble_vec_mean[:,0] = np.mean(ensemble_vec, axis=1)
            ensemble_vec_full[:,:,0] = ensemble_vec
            ensemble_bg[:,0] = ensemble_vec_mean[:,0]
        else:
            # ============================================== NEW BRANCH
            if partitioned_io:
                # ensemble already written directly, mean not needed here
                # in the same form as parallel_write_full_ensemble_from_root
                ens_mean = None
            else:
                # ---- EXISTING behavior, unchanged ----
                shape_ens = comm_world.bcast(shape_ens, root=0)
                _time_init_ensemble_mean_computation = MPI.Wtime()
                ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], icesee_kwargs['Nens'], comm_world, root=0)
                time_init_ensemble_mean_computation += MPI.Wtime() - _time_init_ensemble_mean_computation

                _time_init_file_writing = MPI.Wtime()
                parallel_write_full_ensemble_from_root(0, ens_mean, icesee_kwargs,ensemble_vec,comm_world)
                time_init_file_writing += MPI.Wtime() - _time_init_file_writing
            # ============================================== END NEW

    else:
        # ---- EXISTING size_world > Nens branch, fully unchanged ----
        if rank_world == 0:
            print("[ICESEE] Initializing the ensemble ...")

        if icesee_kwargs["default_run"] and size_world > icesee_kwargs["Nens"]:
            sub_shape = icesee_kwargs['dim_list'][sub_rank]
            icesee_kwargs.update({"statevec_ens":np.zeros((sub_shape, icesee_kwargs["Nens"]))})
            icesee_kwargs.update({"ens_id": color, "rank": sub_rank, "color": color, "comm": subcomm})

            ens = color
            initialilaized_state = model_module.initialize_ensemble(ens,**icesee_kwargs)

            initial_data = {key: subcomm.gather(value, root=0) for key, value in initialilaized_state.items()}
            key_list = list(initial_data.keys())
            state_keys = key_list[:icesee_kwargs["num_state_vars"]]
            if sub_rank == 0:
                for key in key_list:
                    initial_data[key] = np.hstack(initial_data[key])
                    if icesee_kwargs["joint_estimation"] or icesee_kwargs["localization_flag"]:
                        hdim = initial_data[key].shape[0] // icesee_kwargs["total_state_param_vars"]
                    else:
                        hdim = initial_data[key].shape[0] // icesee_kwargs["num_state_vars"]
                    state_block_size = hdim*icesee_kwargs["num_state_vars"]
                    full_block_size = hdim*icesee_kwargs["total_state_param_vars"]
                    if icesee_kwargs.get("random_fields",False):
                        Q_err = np.zeros((full_block_size,full_block_size))
                        for i, sig in enumerate(icesee_kwargs["sig_Q"]):
                            start_idx = i *hdim
                            end_idx = start_idx + hdim
                            Q_err[start_idx:end_idx,start_idx:end_idx] = np.eye(hdim) * sig ** 2
                        _time_init_noise_generation = MPI.Wtime()
                        noise = compute_noise_random_fields(ens, hdim, pos, gs_model, icesee_kwargs["total_state_param_vars"], L_C)
                        time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                        initial_data[key] += noise
                    else:
                        N_size = icesee_kwargs["total_state_param_vars"] * hdim
                        _time_init_noise_generation = MPI.Wtime()
                        icesee_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":icesee_kwargs["total_state_param_vars"]})
                        noise = generate_enkf_field(**icesee_kwargs)
                        time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                        initial_data[key] += noise

                stacked = np.hstack([initial_data[key] for key in initialilaized_state.keys()])
                shape_ens = np.array(stacked.shape,dtype=np.int32)
            else:
                shape_ens = np.empty(2,dtype=np.int32)

            shape_ens = comm_world.bcast(shape_ens, root=0)

            if sub_rank != 0:
                stacked = np.empty(shape_ens,dtype=np.float64)

            all_init = comm_world.gather(stacked if sub_rank == 0 else None, root=0)

            if rank_world == 0:
                all_init = [arr for arr in all_init if isinstance(arr, np.ndarray)]
                ensemble_vec = np.column_stack(all_init)
            else:
                ensemble_vec = np.empty((icesee_kwargs["global_shape"],icesee_kwargs["Nens"]),dtype=np.float64)

            time_init_ensemble_mean_computation = MPI.Wtime()
            ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], icesee_kwargs['Nens'], comm_world, root=0)
            time_init_ensemble_mean_computation = MPI.Wtime() - time_init_ensemble_mean_computation

            _time_init_file_writing = MPI.Wtime()
            parallel_write_full_ensemble_from_root(0, ens_mean, icesee_kwargs,ensemble_vec,comm_world)
            time_init_file_writing += MPI.Wtime() - _time_init_file_writing

        elif icesee_kwargs["sequential_run"]:
            comm_world.Barrier()
            sub_shape = icesee_kwargs['dim_list'][rank_world]
            icesee_kwargs.update({"statevec_ens":np.zeros([icesee_kwargs["global_shape"], icesee_kwargs["Nens"]]),
                                "statevec_ens_mean":np.zeros([icesee_kwargs["global_shape"], icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1]),
                                "statevec_ens_full":np.zeros([icesee_kwargs["global_shape"], icesee_kwargs["Nens"], icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1]),
                                "statevec_bg":np.zeros([icesee_kwargs["global_shape"], icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1])})
            ensemble_bg, ensemble_vec, ensemble_vec_mean, ensemble_vec_full = model_module.initialize_ensemble(**icesee_kwargs)

            gathered_ensemble = comm_world.gather(ensemble_vec[:sub_shape,:], root=0)
            if rank_world == 0:
                ensemble_vec = np.vstack(gathered_ensemble)
                ensemble_vec_mean[:,0] = np.mean(ensemble_vec, axis=1)
                ensemble_vec_full[:,:,0] = ensemble_vec
            else:
                ensemble_vec = np.empty((icesee_kwargs["global_shape"],icesee_kwargs["Nens"]),dtype=np.float64)
                ensemble_vec_mean = np.empty((icesee_kwargs["global_shape"],icesee_kwargs.get("nt",icesee_kwargs["nt"])+1),dtype=np.float64)
                ensemble_vec_full = np.empty((icesee_kwargs["global_shape"],icesee_kwargs["Nens"],icesee_kwargs.get("nt",icesee_kwargs["nt"])+1),dtype=np.float64)

            comm_world.Bcast(ensemble_vec, root=0)
            comm_world.Bcast(ensemble_vec_mean, root=0)
            comm_world.Bcast(ensemble_vec_full, root=0)

    if icesee_kwargs.get("default_run", False):
        return icesee_kwargs, ensemble_vec, time_init_noise_generation, \
               time_init_ensemble_mean_computation, time_init_file_writing, \
                shape_ens, None, None, None
    else:
        return icesee_kwargs, ensemble_vec, time_init_noise_generation, \
               time_init_ensemble_mean_computation, time_init_file_writing, \
                shape_ens,ensemble_bg,  ensemble_vec_mean, ensemble_vec_full


def ensemble_initialization_full_parallel_run(**icesee_kwargs):
    """Initialize the ensemble for the ICESEE model.
    """

    # unpack icesee_kwargs
    model_module   = icesee_kwargs.get("model_module", None)
    comm_world     = icesee_kwargs.get("comm_world", MPI.COMM_WORLD)
    subcomm        = icesee_kwargs.get("subcomm", None)
    color          = icesee_kwargs.get("color", 0)
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
    rounds         = icesee_kwargs.get("rounds", 1)
    subcomm_size_min   = icesee_kwargs.get("subcomm_size_min", 1)
    rng           = icesee_kwargs.get("rng", np.random.default_rng())
    rank_seed = icesee_kwargs.get("rank_seed", 0)
    data_path = icesee_kwargs.get("data_path", "_modeldatasets")
    enkf_parallel_io = icesee_kwargs.get("enkf_parallel_io", None)
    alpha       = icesee_kwargs.get("initial_spread_factor")

    sub_rank     = subcomm.Get_rank()
    rank_world   = comm_world.Get_rank()
    size_world   = comm_world.Get_size()

    time_init_noise_generation = 0.0
    time_init_file_writing     = 0.0
    time_init_ensemble_mean_computation = 0.0

    observed_vars = icesee_kwargs.get("observed_vars", [])
    observed_params = icesee_kwargs.get("observed_params", [])

    all_observed = list(observed_vars) + list(observed_params)

    icesee_kwargs["observed_vars_params"] = all_observed
    icesee_kwargs["all_observed"] = all_observed
    icesee_kwargs["all_observed"] = all_observed
    icesee_kwargs["nd_observed"] = len(all_observed) * (icesee_kwargs["nd"] // icesee_kwargs["total_state_param_vars"])

    if icesee_kwargs["even_distribution"] or (icesee_kwargs["default_run"] and size_world <= icesee_kwargs["Nens"]):
        if icesee_kwargs["default_run"] and size_world <= icesee_kwargs["Nens"] and not (icesee_kwargs.get("sequential_ensemble_initialization", False)):
        # if False:
            if rank_world == 0:
                print("[ICESEE] Initializing the ensemble ...")

            # icesee_kwargs.update({'ens_id': rank_world})
            Nens = icesee_kwargs["Nens"]
            nd = icesee_kwargs.get("nd", icesee_kwargs["nd"])
            icesee_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})

            icesee_kwargs.update({"statevec_ens":np.zeros([icesee_kwargs["nd"], icesee_kwargs["Nens"]])})

            # get the ensemble matrix
            vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)
            # ensemble_vec = np.zeros_like(icesee_kwargs["statevec_ens"])
            # store=f"{data_path}/ensemble_initialization_{color}.zarr"
            # chunk_size = (min(nd, 1000), 1)
            # ensemble_vec = zarr.create_array(store=store, shape=(icesee_kwargs["nd"], icesee_kwargs["Nens"]), chunks=chunk_size, dtype=np.float64, overwrite=True)
            ensemble_vec = np.zeros(nd, dtype=np.float64)

            if icesee_kwargs["joint_estimation"] or icesee_kwargs["localization_flag"]:
                    hdim = nd // icesee_kwargs["total_state_param_vars"]
            else:
                hdim = nd // icesee_kwargs["num_state_vars"]
            state_block_size = hdim * icesee_kwargs["num_state_vars"]

            for round_id in range(rounds):
                ensemble_id = color + (round_id * subcomm_size_min)
                icesee_kwargs.update({'ens_id': ensemble_id})

                if ensemble_id < Nens:
                    # Synchronize the ensemble initialization
                    # subcomm.Barrier()
                    # comm_world.Barrier()
                    ens = ensemble_id

                    # Call the model to initialize the ensemble
                    data = model_module.initialize_ensemble(ens, **icesee_kwargs)
                    for key, value in data.items():
                        # ensemble_vec[indx_map[key], ens] = value
                        ensemble_vec[indx_map[key]] = value

                    # Add process noise in-place to avoid temporary array
                    _time_init_noise_generation = MPI.Wtime()
                    icesee_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":icesee_kwargs["total_state_param_vars"]})
                    noise = generate_enkf_field(**icesee_kwargs)
                    time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                    # ensemble_vec[:,ens] += noise
                    ensemble_vec += alpha*noise

                    _time_init_file_writing = MPI.Wtime()
                    enkf_parallel_io.write_forecast(0, ensemble_vec, ensemble_id)
                    # enkf_parallel_io.datasets[0][:, ens] = ensemble_vec
                    time_init_file_writing += MPI.Wtime() - _time_init_file_writing

        else:
            if rank_world == 0:
                print("[ICESEE] Initializing the ensemble ...")
                icesee_kwargs.update({'ens_id': rank_world})
                if icesee_kwargs["even_distribution"]:
                    icesee_kwargs.update({'rank': rank_world, 'color': color, 'comm': comm_world})
                else:
                    icesee_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})

                nd = icesee_kwargs.get("nd", icesee_kwargs["nd"])

                # get the ensemble matrix
                vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)
                # ensemble_vec = np.zeros_like(icesee_kwargs["statevec_ens"])
                store=f"{data_path}/ensemble_initialization.zarr"
                chunk_size = (min(nd, 1000), 1)
                ensemble_vec = zarr.create_array(store=store, shape=(icesee_kwargs["nd"], icesee_kwargs["Nens"]), chunks=chunk_size, dtype=np.float64, overwrite=True)

                if icesee_kwargs["joint_estimation"] or icesee_kwargs["localization_flag"]:
                        hdim = nd // icesee_kwargs["total_state_param_vars"]
                else:
                    hdim = nd // icesee_kwargs["num_state_vars"]
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
                    _time_init_noise_generation = MPI.Wtime()
                    N_size = icesee_kwargs["total_state_param_vars"] * hdim
                    # noise = generate_pseudo_random_field_1d(N_size,np.sqrt(Lx*Ly), len_scale, verbose=True)
                    icesee_kwargs.update({"ii_sig": None, "Lx_dim": np.sqrt(Lx*Ly), "noise_dim": hdim, "num_vars":icesee_kwargs["total_state_param_vars"]})
                    noise = generate_enkf_field(**icesee_kwargs)
                    time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                    ensemble_vec[:,ens] += alpha*noise
                    # for ii, sig in enumerate(icesee_kwargs["sig_Q"]):
                    #     if ii <=icesee_kwargs["num_state_vars"]:
                    #         start_idx = ii * hdim
                    #         end_idx = start_idx + hdim
                    #         ensemble_vec[start_idx:end_idx, ens] += noise[start_idx:end_idx] * sig

                    # enkf_parallel_io.write_forecast(0, ensemble_vec[:,ens], ens)
                    # enkf_parallel_io.datasets[0][:, ens] = ensemble_vec[:,ens]
                    # print(f"[ICESEE] Rank {rank_world}: Ensemble initialization completed and written to disk with norm {np.linalg.norm(ensemble_vec)}")

                shape_ens = np.array(ensemble_vec.shape,dtype=np.int32)
                # print(f"[ICESEE] Rank {rank_world}: Ensemble initialization completed for all members.")


            else:
                ensemble_vec = np.empty((icesee_kwargs["nd"],icesee_kwargs["Nens"]),dtype=np.float64)
                shape_ens = np.empty(2,dtype=np.int32)
                # pos, gs_model, L_C

            _time_init_file_writing = MPI.Wtime()
            # scatter  enkf_parallel_io.nd_local_world of the ensemble to all processors
            localshape = enkf_parallel_io.nd_local_world
            all_local_shapes = comm_world.gather(localshape)
            if rank_world == 0:
                counts_rows = np.array(all_local_shapes)
                displacement_rows = np.insert(np.cumsum(counts_rows), 0, 0)[0:-1]
                counts_rows = counts_rows * icesee_kwargs["Nens"]
                displacement_rows = displacement_rows * icesee_kwargs["Nens"]
            else:
                counts_rows = None
                displacement_rows = None

            local_ensemble = np.empty((localshape, icesee_kwargs["Nens"]), dtype=np.float64)
            comm_world.Scatterv([ensemble_vec, counts_rows, displacement_rows, MPI.DOUBLE], local_ensemble, root=0)
            enkf_parallel_io.datasets[0][localshape, :] = local_ensemble
            time_init_file_writing += MPI.Wtime() - _time_init_file_writing

        comm_world.Barrier()
        _time_init_ensemble_mean_computation = MPI.Wtime()
        # enkf_parallel_io.compute_forecast_mean_chunked(0)
        enkf_parallel_io.compute_forecast_mean_chunked_v2(k=0,flag="initial")
        # ens_mean = enkf_parallel_io.compute_forecast_mean(0)
        # ens_mean = .datasets[0][:, :].mean(axis=1)
        time_init_ensemble_mean_computation += MPI.Wtime() - _time_init_ensemble_mean_computation

        # now reset the model_nprocs
        if rank_world == 0:
            diff = total_cores - base_total_procs
            if diff >= 0:
                # split the diff amaongest all processors
                min_model_nprocs = max(model_nprocs-1, 1)
                if icesee_kwargs.get('ICESEE_PERFORMANCE_TEST') or env_flag("ICESEE_PERFORMANCE_TEST", default=False):
                    model_nprocs = icesee_kwargs.get("model_nprocs", 1)
                else:
                    model_nprocs = max(min_model_nprocs, model_nprocs + (diff // size_world))
            else:
                model_nprocs = model_nprocs

        model_nprocs = comm_world.bcast(model_nprocs, root=0)
        icesee_kwargs.update({'model_nprocs': model_nprocs})

    else:
        if rank_world == 0:
            print("[ICESEE] Initializing the ensemble ...")

        if icesee_kwargs["default_run"] and size_world > icesee_kwargs["Nens"]:
            # debug
            sub_shape = icesee_kwargs['dim_list'][sub_rank]
            icesee_kwargs.update({"statevec_ens":np.zeros((sub_shape, icesee_kwargs["Nens"]))})

            icesee_kwargs.update({"ens_id": color, "rank": sub_rank, "color": color, "comm": subcomm})

            # ensemble_vec, shape_ens  = model_module.initialize_ensemble_debug(color,**icesee_kwargs)
            # ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], icesee_kwargs['Nens'], comm_world, root=0)
            # parallel_write_full_ensemble_from_root(0, ens_mean, icesee_kwargs,ensemble_vec,comm_world)
            # -----------------------------------------------------

            ens = color
            # icesee_kwargs.update({"statevec_ens":np.zeros((icesee_kwargs['global_shape'], icesee_kwargs["Nens"]))})
            initialilaized_state = model_module.initialize_ensemble(ens,**icesee_kwargs)
            # ensemble_vec, shape_ens = gather_and_broadcast_data_default_run(initialilaized_state, subcomm, sub_rank, comm_world, rank_world, icesee_kwargs)
            # ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], icesee_kwargs['Nens'], comm_world, root=0)
            # parallel_write_full_ensemble_from_root(0, ens_mean, icesee_kwargs,ensemble_vec,comm_world)
            # ensemble_vec = BM.bcast(ensemble_vec, comm_world)

            initial_data = {key: subcomm.gather(value, root=0) for key, value in initialilaized_state.items()}
            key_list = list(initial_data.keys())
            state_keys = key_list[:icesee_kwargs["num_state_vars"]]
            if sub_rank == 0:
                # for key in initial_data:
                for key in key_list:
                    initial_data[key] = np.hstack(initial_data[key])
                    if icesee_kwargs["joint_estimation"] or icesee_kwargs["localization_flag"]:
                        hdim = initial_data[key].shape[0] // icesee_kwargs["total_state_param_vars"]
                    else:
                        hdim = initial_data[key].shape[0] // icesee_kwargs["num_state_vars"]
                    state_block_size = hdim*icesee_kwargs["num_state_vars"]
                    full_block_size = hdim*icesee_kwargs["total_state_param_vars"]
                    # if key in state_keys:
                        # noise = np.random.normal(0, 0.1, state_block_size)
                        # Q_err = np.eye(state_block_size) * icesee_kwargs["sig_Q"] ** 2
                        # Q_err = np.eye(state_block_size) * 0.01 ** 2
                    if icesee_kwargs.get("random_fields",False):
                        Q_err = np.zeros((full_block_size,full_block_size))
                        for i, sig in enumerate(icesee_kwargs["sig_Q"]):
                            start_idx = i *hdim
                            end_idx = start_idx + hdim
                            Q_err[start_idx:end_idx,start_idx:end_idx] = np.eye(hdim) * sig ** 2

                        # noise = multivariate_normal.rvs(mean=np.zeros(state_block_size), cov=Q_err)
                        _time_init_noise_generation = MPI.Wtime()
                        noise = compute_noise_random_fields(ens, hdim, pos, gs_model, icesee_kwargs["total_state_param_vars"], L_C)
                        time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                        # initial_data[key][:state_block_size] += noise[:state_block_size]
                        # noise = noise / np.max(np.abs(noise))
                        initial_data[key] += alpha*noise
                    else:
                        N_size = icesee_kwargs["total_state_param_vars"] * hdim
                        _time_init_noise_generation = MPI.Wtime()
                        icesee_kwargs.update({"ii_sig": None, "hdim":hdim, "num_vars":icesee_kwargs["total_state_param_vars"]})
                        noise = generate_enkf_field(**icesee_kwargs)
                        time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                        initial_data[key] += noise
                        # for ii, sig in enumerate(icesee_kwargs["sig_Q"]):
                        #     start_idx = ii *hdim
                        #     end_idx = start_idx + hdim
                        #     initial_data[key][start_idx:end_idx] += noise[start_idx:end_idx]*sig

                # stack all variables together into a single array
                stacked = np.hstack([initial_data[key] for key in initialilaized_state.keys()])
                shape_ens = np.array(stacked.shape,dtype=np.int32)
            else:
                shape_ens = np.empty(2,dtype=np.int32)

            # broadcast the shape of the initialized ensemble
            shape_ens = comm_world.bcast(shape_ens, root=0)

            if sub_rank != 0:
                stacked = np.empty(shape_ens,dtype=np.float64)

            all_init = comm_world.gather(stacked if sub_rank == 0 else None, root=0)

            if rank_world == 0:
                all_init = [arr for arr in all_init if isinstance(arr, np.ndarray)]
                ensemble_vec = np.column_stack(all_init)
                # print(f"[ICESEE] Shape of the ensemble: {ensemble_vec.shape}")
            else:
                ensemble_vec = np.empty((icesee_kwargs["global_shape"],icesee_kwargs["Nens"]),dtype=np.float64)

            time_init_ensemble_mean_computation = MPI.Wtime()
            ens_mean = ParallelManager().compute_mean_matrix_from_root(ensemble_vec, shape_ens[0], icesee_kwargs['Nens'], comm_world, root=0)
            time_init_ensemble_mean_computation += MPI.Wtime() - _time_init_ensemble_mean_computation

            _time_init_file_writing = MPI.Wtime()
            parallel_write_full_ensemble_from_root(0, ens_mean, icesee_kwargs,ensemble_vec,comm_world)
            time_init_file_writing += MPI.Wtime() - _time_init_file_writing

        elif icesee_kwargs["sequential_run"]:
            comm_world.Barrier()
            sub_shape = icesee_kwargs['dim_list'][rank_world]
            icesee_kwargs.update({"statevec_ens":np.zeros([icesee_kwargs["global_shape"], icesee_kwargs["Nens"]]),
                                "statevec_ens_mean":np.zeros([icesee_kwargs["global_shape"], icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1]),
                                "statevec_ens_full":np.zeros([icesee_kwargs["global_shape"], icesee_kwargs["Nens"], icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1]),
                                "statevec_bg":np.zeros([icesee_kwargs["global_shape"], icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1])})
            ensemble_bg, ensemble_vec, ensemble_vec_mean, ensemble_vec_full = model_module.initialize_ensemble(**icesee_kwargs)

            # gather from every rank to rank 0
            gathered_ensemble = comm_world.gather(ensemble_vec[:sub_shape,:], root=0)
            if rank_world == 0:
                ensemble_vec = np.vstack(gathered_ensemble)
                print(f"[ICESEE] Shape of the ensemble: {ensemble_vec.shape}")
                ensemble_vec_mean[:,0] = np.mean(ensemble_vec, axis=1)
                ensemble_vec_full[:,:,0] = ensemble_vec
            else:
                ensemble_vec = np.empty((icesee_kwargs["global_shape"],icesee_kwargs["Nens"]),dtype=np.float64)
                ensemble_vec_mean = np.empty((icesee_kwargs["global_shape"],icesee_kwargs.get("nt",icesee_kwargs["nt"])+1),dtype=np.float64)
                ensemble_vec_full = np.empty((icesee_kwargs["global_shape"],icesee_kwargs["Nens"],icesee_kwargs.get("nt",icesee_kwargs["nt"])+1),dtype=np.float64)

            # else:
            #     ensemble_bg = np.empty((icesee_kwargs["global_shape"],icesee_kwargs.get("nt",icesee_kwargs["nt"])+1),dtype=np.float64)
            #     ensemble_vec = np.empty((icesee_kwargs["global_shape"],icesee_kwargs["Nens"]),dtype=np.float64)
            #     ensemble_vec_mean = np.empty((icesee_kwargs["global_shape"],icesee_kwargs.get("nt",icesee_kwargs["nt"])+1),dtype=np.float64)
            #     ensemble_vec_full = np.empty((icesee_kwargs["global_shape"],icesee_kwargs["Nens"],icesee_kwargs.get("nt",icesee_kwargs["nt"])+1),dtype=np.float64)

            # # Bcast the ensemble
            # comm_world.Bcast(ensemble_bg, root=0)
            comm_world.Bcast(ensemble_vec, root=0)
            comm_world.Bcast(ensemble_vec_mean, root=0)
            comm_world.Bcast(ensemble_vec_full, root=0)

            # hdim = ensemble_vec.shape[0] // icesee_kwargs["total_state_param_vars"]
            # print(f"[ICESEE] rank: {rank_world}, subrank: {sub_rank}, min ensemble: {np.min(ensemble_vec[hdim,:])}, max ensemble: {np.max(ensemble_vec[hdim,:])}")

    if icesee_kwargs.get("default_run", False):
        return icesee_kwargs, None, time_init_noise_generation, \
               time_init_ensemble_mean_computation,time_init_file_writing, \
                None, None, None, None
    else:
        return icesee_kwargs, ensemble_vec, time_init_noise_generation, \
               time_init_ensemble_mean_computation, time_init_file_writing, \
                shape_ens,ensemble_bg,  ensemble_vec_mean, ensemble_vec_full
