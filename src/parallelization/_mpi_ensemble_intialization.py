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
from contextlib import contextmanager
from mpi4py import MPI

from ICESEE.src.utils.tools import icesee_get_index, env_flag
from ICESEE.src.utils.state_ownership import (
    resolve_state_ownership,
    combine_member_state,
    ensure_state_array,
)
from ICESEE.src.utils.random_streams import initialization_seed
from ICESEE.src.run_model_da._error_generation import compute_Q_err_random_fields, \
                              compute_noise_random_fields, \
                              generate_pseudo_random_field_1d, \
                              generate_pseudo_random_field_2D, \
                              generate_enkf_field, \
                              generate_initial_member_increment

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


class _EnsembleShapeProxy:
    """Compatibility view of an ensemble matrix without allocating it.

    Older application adapters inspect ``statevec_ens.shape`` during member
    initialization even though they never read or write the matrix.  Mode 2
    intentionally has no in-memory ``nd x Nens`` array, so expose only the
    legacy shape contract while keeping memory use O(1).
    """

    __slots__ = ("shape",)

    def __init__(self, nd, nens):
        self.shape = (int(nd), int(nens))


@contextmanager
def _member_initialization_context(icesee_kwargs, ensemble_id, stream=0):
    """Make application initialization independent of MPI member ordering.

    Several model adapters still consume NumPy's legacy global generator,
    while newer adapters consume ``rng`` or ``rank_seed`` from the runtime
    dictionary.  Supplying all three interfaces gives execution modes 1 and 2
    the same member stream without requiring application-specific changes.
    """
    seed = initialization_seed(
        icesee_kwargs.get("base_seed", 42), ensemble_id, stream
    )
    old_state = np.random.get_state()
    np.random.seed(seed)
    local_kwargs = dict(icesee_kwargs)
    local_kwargs.update(
        {
            "ens_id": int(ensemble_id),
            "seed": seed,
            "rank_seed": seed,
            "rng": np.random.default_rng(seed),
        }
    )
    if (
        "statevec_ens" not in local_kwargs
        and "nd" in local_kwargs
        and "Nens" in local_kwargs
    ):
        local_kwargs["statevec_ens"] = _EnsembleShapeProxy(
            local_kwargs["nd"], local_kwargs["Nens"]
        )
    try:
        yield local_kwargs
    finally:
        np.random.set_state(old_state)


# generate_initial_member_increment now lives in
# ICESEE.src.run_model_da._error_generation (imported above) instead of
# being defined here. It used to be mode-2-only, which is why modes 0/1
# (src/EnKF/_ensemble_initialization.py) could not use it and instead kept
# an older, separately-broken inline initial-noise implementation
# (unreseeded shared RNG across members, no sig_Q scaling). Relocating it
# to a model-agnostic, MPI-free module lets every execution mode share one
# correctly seeded/scaled implementation -- see _error_generation.py's
# docstring on the function itself, and CHANGELOG.md for the historical
# reconciliation this addressed.


def apply_initial_spread_filebacked(enkf_parallel_io, alpha):
    """Apply ``mean + alpha * anomaly`` without materializing nd-by-Nens.

    Every world rank owns a disjoint row slab.  Member columns are read in
    small batches, their mean is accumulated in float64, and the same slab is
    rewritten in bounded member batches.  This is the file-backed equivalent
    of mode 1's dense initial-spread operation.
    """
    if alpha is None:
        alpha = 1.0
    alpha = float(alpha)
    io = enkf_parallel_io
    io._ensure_batch(0)
    batch_index = 0 - io.current_batch_start
    dataset = io.datasets[batch_index]
    start, stop = io.nd_start_world, io.nd_end_world
    local_rows = max(0, stop - start)
    member_batch = max(
        1,
        int(io.icesee_kwargs.get("initialization_member_chunk_size", 1)),
    )
    if local_rows:
        local_sum = np.zeros(local_rows, dtype=np.float64)
        for first in range(0, io.nens, member_batch):
            last = min(io.nens, first + member_batch)
            local_sum += np.asarray(
                dataset[start:stop, first:last], dtype=np.float64
            ).sum(axis=1)
        local_mean = local_sum / float(io.nens)
        for first in range(0, io.nens, member_batch):
            last = min(io.nens, first + member_batch)
            block = np.asarray(
                dataset[start:stop, first:last], dtype=np.float64
            )
            block = local_mean[:, None] + alpha * (
                block - local_mean[:, None]
            )
            dataset[start:stop, first:last] = block.astype(
                io.storage_dtype, copy=False
            )
    io.mpi_comm.Barrier()

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

                    with _member_initialization_context(
                        icesee_kwargs, ens
                    ) as member_kwargs:
                        data = model_module.initialize_ensemble(
                            ens, **member_kwargs
                        )
                    for key, value in data.items():
                        ensemble_vec[indx_map[key], ens] = value

                    _time_init_noise_generation = MPI.Wtime()
                    increment, noise = generate_initial_member_increment(
                        hdim, icesee_kwargs, ens, ensemble_vec.shape[0]
                    )
                    time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                    ensemble_vec[:, ens] += increment

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
                    with _member_initialization_context(
                        icesee_kwargs, ens
                    ) as member_kwargs:
                        data = model_module.initialize_ensemble(
                            ens, **member_kwargs
                        )
                    for key, value in data.items():
                        ensemble_vec[indx_map[key],ens] = value

                    _time_init_noise_generation = MPI.Wtime()
                    increment, noise = generate_initial_member_increment(
                        hdim,
                        icesee_kwargs,
                        ens,
                        ensemble_vec.shape[0],
                    )
                    time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                    ensemble_vec[:, ens] += increment

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

    # A spare rank under the hierarchical resource plan (resource_plan.py)
    # is handed MPI.COMM_NULL for subcomm (Get_rank() raises MPI_ERR_COMM
    # on it); icesee_kwargs["sub_rank"] was already safely computed as
    # None for that case by icesee_mpi_ens_distribution and is preferred
    # here over recomputing from subcomm directly.
    if "sub_rank" in icesee_kwargs:
        sub_rank = icesee_kwargs["sub_rank"]
    elif subcomm is not None and subcomm != MPI.COMM_NULL:
        sub_rank = subcomm.Get_rank()
    else:
        sub_rank = None
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

    # ranks_per_model == 1 (resolved by the centralized resource plan --
    # resource_plan.py), not a raw size_world/Nens comparison, is the
    # actual condition for "no ICESEE-managed multi-rank model
    # communicator". See the matching change and comment in
    # _mpi_generate_true_wrong_state.py. Only this function
    # (ensemble_initialization_full_parallel_run, mode 2) is affected --
    # ensemble_initialization (modes 0/1) is untouched.
    if icesee_kwargs["even_distribution"] or (icesee_kwargs["default_run"] and icesee_kwargs.get("ranks_per_model", 1) == 1):
        if icesee_kwargs["default_run"] and icesee_kwargs.get("ranks_per_model", 1) == 1 and not (icesee_kwargs.get("sequential_ensemble_initialization", False)):
        # if False:
            if rank_world == 0:
                print("[ICESEE] Initializing the ensemble ...")

            Nens = icesee_kwargs["Nens"]

            # Every world rank reaches this point identically. In the last
            # round `ensemble_id < Nens` below is false for some ranks and
            # true for others (Nens need not be a multiple of
            # subcomm_size_min), so `write_forecast` inside that guard runs
            # on a rank subset and must not be the one to open the shared
            # HDF5 batch window itself -- see
            # EnKF_fully_parallel_IO._require_batch.
            enkf_parallel_io._ensure_batch(0)

            # A spare rank (color is None -- resource_plan.py) belongs to
            # no model group and must not enter this round loop, nor touch
            # anything derived from `subcomm` -- it is MPI.COMM_NULL for a
            # spare rank (Stage 4A/resource_plan.py), and icesee_get_index
            # (called below to build indx_map for this model group) calls
            # collectives on whatever communicator it is handed. This is
            # only reachable via an explicit ranks_per_model == 1 request
            # combined with world_size > Nens; the legacy/auto policy never
            # produces a spare rank when ranks_per_model resolves to 1, so
            # no previously working configuration reached this branch with
            # a spare rank before Stage 4A's ResourcePlan made that
            # combination requestable.
            if color is not None:
                nd = icesee_kwargs.get("nd", icesee_kwargs["nd"])
                icesee_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})

                # get the ensemble matrix
                vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)
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
                        with _member_initialization_context(
                            icesee_kwargs, ens
                        ) as member_kwargs:
                            data = model_module.initialize_ensemble(
                                ens, **member_kwargs
                            )
                        for key, value in data.items():
                            # ensemble_vec[indx_map[key], ens] = value
                            ensemble_vec[indx_map[key]] = value

                        # Add process noise in-place to avoid temporary array
                        _time_init_noise_generation = MPI.Wtime()
                        increment, noise = generate_initial_member_increment(
                            hdim, icesee_kwargs, ens, ensemble_vec.shape[0]
                        )
                        time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                        ensemble_vec += increment

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
                # Sequential initialization is intentionally member-streamed:
                # rank 0 holds one state vector, writes it, and reuses the buffer.
                ensemble_vec = np.zeros(nd, dtype=np.float64)

                if icesee_kwargs["joint_estimation"] or icesee_kwargs["localization_flag"]:
                        hdim = nd // icesee_kwargs["total_state_param_vars"]
                else:
                    hdim = nd // icesee_kwargs["num_state_vars"]
            else:
                nd = None
                indx_map = None
                ensemble_vec = None
                hdim = None

            # The model is initialized sequentially on rank zero, but opening
            # and writing a mode-2 HDF5 shard is collective.  Every rank must
            # therefore enter the loop and participate in the row-slab
            # scatter/write for every member.
            for ens in range(icesee_kwargs["Nens"]):
                if rank_world == 0:
                    with _member_initialization_context(
                        icesee_kwargs, ens
                    ) as member_kwargs:
                        data = model_module.initialize_ensemble(
                            ens, **member_kwargs
                        )

                    # iterate over the data and update the ensemble
                    ensemble_vec.fill(0.0)
                    for key, value in data.items():
                        ensemble_vec[indx_map[key]] = value

                    # --->
                    # noise = compute_noise_random_fields(ens, hdim, pos, gs_model, icesee_kwargs["total_state_param_vars"], L_C)
                    # ensemble_vec[:,ens] += noise
                    #----->
                    _time_init_noise_generation = MPI.Wtime()
                    increment, noise = generate_initial_member_increment(
                        hdim, icesee_kwargs, ens, ensemble_vec.shape[0]
                    )
                    time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                    ensemble_vec += increment

                _time_init_file_writing = MPI.Wtime()
                enkf_parallel_io.write_forecast_from_root(
                    0, ensemble_vec, ens, root=0
                )
                time_init_file_writing += MPI.Wtime() - _time_init_file_writing

            if rank_world == 0:
                shape_ens = np.array((nd, icesee_kwargs["Nens"]), dtype=np.int32)
            else:
                shape_ens = np.empty(2,dtype=np.int32)

            shape_ens = comm_world.bcast(shape_ens, root=0)

        comm_world.Barrier()
        apply_initial_spread_filebacked(enkf_parallel_io, alpha)
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

        if icesee_kwargs["default_run"] and icesee_kwargs.get("ranks_per_model", 1) > 1:
            # Every world rank reaches this point identically (unlike the
            # `if sub_rank == 0:` write below), so this is the one place
            # that may safely open the shared HDF5 batch window ahead of
            # the subcommunicator-root-only `write_forecast(0, ...)` a few
            # lines down -- see EnKF_fully_parallel_IO._require_batch. Must
            # stay unconditional (not inside the `if color is not None:`
            # guard below) so a spare rank still joins this comm_world
            # collective.
            enkf_parallel_io._ensure_batch(0)

            # A spare rank (color is None -- resource_plan.py) belongs to
            # no model group and must not touch subcomm-level state at
            # all, but it still must reach the comm_world.Barrier() and
            # compute_forecast_mean_chunked_v2() below, which every active
            # rank in this branch also reaches -- so only this inner body
            # is guarded, not the branch itself.
            if color is not None:
                icesee_kwargs.update({"rank": sub_rank, "color": color, "comm": subcomm})

                # This branch used to hard-code `ens = color`, silently
                # initializing only the round-0 member of every model
                # group and leaving every later round's member at HDF5's
                # zero fill value (never written at all) -- reproduced via
                # a real Icepack P=4,Nens=4,ranks_per_model=2 (2 groups, 2
                # rounds) run: members 2 and 3 (round 1) fed an all-zero
                # state into the first forecast solve, which immediately
                # diverged with DIVERGED_FNORM_NAN. The forecast step
                # (_mpi_forecast_functions.py's `_run_one_ensemble`)
                # already loops every round correctly with exactly this
                # `color + round_id * subcomm_size_min` formula; ensemble
                # initialization must agree, or a later round's member
                # never gets a real initial condition to forecast from.
                for round_id in range(rounds):
                    ens = color + round_id * subcomm_size_min
                    if ens < 0 or ens >= icesee_kwargs["Nens"]:
                        continue
                    icesee_kwargs.update({"ens_id": ens})

                    # A dense np.zeros((global_shape, Nens)) placeholder was tried
                    # here before (see the commented-out line this replaced) but
                    # abandoned -- it would materialize a full dense ensemble
                    # matrix, defeating mode 2's bounded-memory design. The
                    # `Nens >= size_world` round-loop branch above solves the same
                    # "some adapters read icesee_kwargs['statevec_ens'].shape"
                    # problem with a O(1) shape-only proxy instead
                    # (_member_initialization_context); this branch was missing
                    # that wrapper entirely, so any adapter reading
                    # icesee_kwargs["statevec_ens"] (e.g. Lorenz-96's
                    # initialize_ensemble) raised a bare KeyError here.
                    with _member_initialization_context(icesee_kwargs, ens) as member_kwargs:
                        initialilaized_state = model_module.initialize_ensemble(ens, **member_kwargs)

                    # Replicated: every rank in the subcomm already computed the
                    # identical full member -- no communication needed. Distributed:
                    # gathers+concatenates each rank's disjoint slice, exactly as
                    # the previous unconditional gather+hstack did. See
                    # src/utils/state_ownership.py.
                    ownership = resolve_state_ownership(icesee_kwargs, subcomm)
                    initial_data = combine_member_state(
                        ownership, subcomm, initialilaized_state, root=0
                    )
                    if sub_rank == 0:
                        key_list = list(initial_data.keys())
                        # hdim (the per-variable block length) is a property of
                        # this member's total state size, not of any one variable's
                        # own value -- it does not depend on `key` and must be
                        # computed from global_shape, matching every other correct
                        # use of this quantity in this codebase (e.g.
                        # _mpi_generate_true_wrong_state.py's `hdim = nd //
                        # num_state_vars`).
                        if icesee_kwargs["joint_estimation"] or icesee_kwargs["localization_flag"]:
                            hdim = icesee_kwargs["global_shape"] // icesee_kwargs["total_state_param_vars"]
                        else:
                            hdim = icesee_kwargs["global_shape"] // icesee_kwargs["num_state_vars"]

                        for key in key_list:
                            # A variable's local block can be a single element
                            # (hdim == 1); some adapters (e.g. Lorenz-96) return a
                            # bare NumPy scalar rather than a one-element array for
                            # that case -- both are the same state component, just
                            # different NumPy representations. Normalize once, here,
                            # generically (not by model name) before any
                            # shape-derived logic or arithmetic below. See
                            # ensure_state_array's docstring.
                            initial_data[key] = ensure_state_array(initial_data[key])

                        # stack all variables together into a single array
                        stacked = np.hstack([initial_data[key] for key in initialilaized_state.keys()])
                        if stacked.size != icesee_kwargs["global_shape"]:
                            raise ValueError(
                                "Initialized member has the wrong state-vector size: "
                                f"member {ens} has {stacked.size}, expected "
                                f"{icesee_kwargs['global_shape']}."
                            )

                        # Decomposition-invariant initial perturbation: this used
                        # to call generate_enkf_field per variable directly here,
                        # with no rng/seed explicitly keyed to (base_seed, ens,
                        # variable_index) -- generate_enkf_field's own fallback
                        # (icesee_kwargs.get("rng"), else default_rng(base_seed))
                        # is not ensemble-member-aware, so this branch's
                        # perturbation could depend on however many times a
                        # freshly-seeded rng had already been consumed rather
                        # than on which member/variable it was for. The
                        # ranks_per_model == 1 branch above never had this
                        # problem because it already called
                        # generate_initial_member_increment -- the shared,
                        # model/MPI-agnostic, (base_seed, ensemble_id,
                        # variable_index + 1)-seeded generator every other
                        # execution mode uses (see that function's own
                        # docstring). Calling the same function here, on the
                        # already-assembled GLOBAL `stacked` vector (not a
                        # per-rank local piece), means one member's initial
                        # perturbation is now the same logical realization
                        # regardless of how many ranks/model groups spatially
                        # partition it -- decomposition-invariant by
                        # construction, not by coincidence.
                        _time_init_noise_generation = MPI.Wtime()
                        increment, _raw_increment = generate_initial_member_increment(
                            hdim, icesee_kwargs, ens, stacked.size
                        )
                        time_init_noise_generation += MPI.Wtime() - _time_init_noise_generation
                        stacked = stacked + increment

                        # Each subcommunicator root owns exactly one initialized member.
                        # Write that column directly into the mode-2 shard instead of
                        # gathering an nd-by-Nens matrix on world rank zero.
                        _time_init_file_writing = MPI.Wtime()
                        enkf_parallel_io.write_forecast(0, stacked, ens)
                        time_init_file_writing += MPI.Wtime() - _time_init_file_writing

            comm_world.Barrier()
            _time_init_ensemble_mean_computation = MPI.Wtime()
            enkf_parallel_io.compute_forecast_mean_chunked_v2(k=0, flag="initial")
            time_init_ensemble_mean_computation += (
                MPI.Wtime() - _time_init_ensemble_mean_computation
            )

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
