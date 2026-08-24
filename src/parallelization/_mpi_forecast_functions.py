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
import tempfile
import json
from contextlib import contextmanager
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


def process_noise_is_due(icesee_kwargs, k=None):
    """Return whether process noise is applied at this forecast step.

    ``observations`` is the historical mode-1 contract for the normal
    Nens >= MPI-size path.  Keeping the decision here prevents execution
    mode and MPI layout from silently changing the stochastic model.
    """
    schedule = str(
        icesee_kwargs.get("process_noise_schedule", "observations")
    ).strip().lower()
    if schedule in {"none", "off", "disabled"}:
        return False
    if schedule in {"every_step", "all", "forecast"}:
        return True
    if schedule not in {"observations", "observation", "analysis"}:
        raise ValueError(
            "process_noise_schedule must be 'observations', 'every_step', "
            f"or 'none'; got {schedule!r}"
        )

    if k is None:
        k = int(icesee_kwargs.get("k", 0))
    obs_index = np.asarray(icesee_kwargs.get("obs_index", []), dtype=np.int64)
    number_obs = int(
        icesee_kwargs.get("number_obs_instants", obs_index.size)
    )
    return bool(np.any(obs_index[:number_obs] == int(k)))


def _process_noise_seed(base_seed, timestep, ens_id, variable_index):
    """Stable member-keyed seed independent of rank and execution order."""
    words = np.random.SeedSequence(
        [
            int(base_seed) & 0xFFFFFFFF,
            int(timestep) & 0xFFFFFFFF,
            int(ens_id) & 0xFFFFFFFF,
            int(variable_index) & 0xFFFFFFFF,
            0x1CE5EE,
        ]
    ).generate_state(1, dtype=np.uint32)
    return int(words[0])


def _forecast_member_seed(base_seed, timestep, ens_id):
    """Stable application-forecast seed independent of MPI scheduling."""
    return int(
        np.random.SeedSequence(
            [
                int(base_seed) & 0xFFFFFFFF,
                int(timestep) & 0xFFFFFFFF,
                int(ens_id) & 0xFFFFFFFF,
                0xF04ECA57,
            ]
        ).generate_state(1, dtype=np.uint32)[0]
    )


@contextmanager
def forecast_member_context(icesee_kwargs, ens_id):
    """Expose one deterministic RNG stream to every model adapter.

    Legacy adapters consume NumPy's global generator, while newer adapters
    consume ``rng`` or ``rank_seed``.  A member/timestep keyed stream makes
    these interfaces identical in execution modes 1 and 2 and invariant to
    rank count, round ordering, and restart scheduling.
    """
    seed = _forecast_member_seed(
        icesee_kwargs.get("base_seed", 42),
        icesee_kwargs.get("k", 0),
        ens_id,
    )
    old_state = np.random.get_state()
    np.random.seed(seed)
    local_kwargs = dict(icesee_kwargs)
    local_kwargs.update(
        {
            "ens_id": int(ens_id),
            "seed": seed,
            "rank_seed": seed,
            "rng": np.random.default_rng(seed),
        }
    )
    try:
        yield local_kwargs
    finally:
        np.random.set_state(old_state)


def _process_noise_state_path(icesee_kwargs, ens_id):
    base = icesee_kwargs.get("_modelrun_datasets") or icesee_kwargs.get(
        "data_path", "_modelrun_datasets"
    )
    return os.path.join(
        str(base), "_process_noise_state", f"member_{int(ens_id):06d}.h5"
    )


def _read_member_process_noise(icesee_kwargs, ens_id, size, timestep):
    """Read the member's AR(1) state, resetting stale state on a fresh run."""
    path = _process_noise_state_path(icesee_kwargs, ens_id)
    if not os.path.exists(path):
        return np.zeros(int(size), dtype=np.float64)
    try:
        with h5py.File(path, "r") as handle:
            last_k = int(handle.attrs.get("last_timestep", -1))
            stored_seed = int(handle.attrs.get("base_seed", -1))
            values = np.asarray(handle["noise"], dtype=np.float64)
        # A state from the same or a later cycle belongs to an older fresh run.
        if (
            last_k >= int(timestep)
            or values.size != int(size)
            or stored_seed != int(icesee_kwargs.get("base_seed", 42))
        ):
            return np.zeros(int(size), dtype=np.float64)
        return values
    except (OSError, KeyError, ValueError):
        return np.zeros(int(size), dtype=np.float64)


def _write_member_process_noise(icesee_kwargs, ens_id, values, timestep):
    """Atomically checkpoint one member's AR(1) process-noise state."""
    path = _process_noise_state_path(icesee_kwargs, ens_id)
    directory = os.path.dirname(path)
    os.makedirs(directory, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".member_{int(ens_id):06d}_", suffix=".h5", dir=directory
    )
    os.close(fd)
    try:
        with h5py.File(temporary, "w") as handle:
            handle.create_dataset("noise", data=np.asarray(values, dtype=np.float64))
            handle.attrs["last_timestep"] = int(timestep)
            handle.attrs["base_seed"] = int(icesee_kwargs.get("base_seed", 42))
            handle.flush()
        os.replace(temporary, path)
    except Exception:
        try:
            os.remove(temporary)
        except OSError:
            pass
        raise


def add_member_process_noise(ensemble_vec, ens_id, icesee_kwargs):
    """Apply restartable member-specific AR(1) noise to state variables.

    Random fields are keyed by (base seed, timestep, ensemble member,
    variable), so changing MPI rank counts, subcommunicator sizes, or member
    scheduling cannot alter the realization.  The AR state is checkpointed per
    member rather than shared by members in execution order.
    """
    k = int(icesee_kwargs.get("k", 0))
    if not process_noise_is_due(icesee_kwargs, k):
        return ensemble_vec

    num_state_vars = int(icesee_kwargs["num_state_vars"])
    total_vars = int(icesee_kwargs["total_state_param_vars"])
    if icesee_kwargs["joint_estimation"] or icesee_kwargs["localization_flag"]:
        hdim = ensemble_vec.shape[0] // total_vars
    else:
        hdim = ensemble_vec.shape[0] // num_state_vars
    state_block_size = hdim * num_state_vars

    previous = _read_member_process_noise(
        icesee_kwargs, ens_id, state_block_size, k
    )
    alpha = float(icesee_kwargs.get("alpha", 0.0))
    rho = float(icesee_kwargs.get("rho", 1.0))
    dt = float(icesee_kwargs.get("dt", 1.0))
    Lx = float(icesee_kwargs.get("Lx", 1.0))
    Ly = float(icesee_kwargs.get("Ly", 1.0))
    sig_q = list(icesee_kwargs.get("sig_Q", []))
    base_seed = int(icesee_kwargs.get("base_seed", 42))

    updated_noise = np.empty(state_block_size, dtype=np.float64)
    noise_increment = np.empty(state_block_size, dtype=np.float64)
    old_rng_state = np.random.get_state()
    try:
        for ii in range(num_state_vars):
            start, stop = ii * hdim, (ii + 1) * hdim
            np.random.seed(_process_noise_seed(base_seed, k, ens_id, ii))
            noise_kwargs = dict(icesee_kwargs)
            noise_kwargs.update(
                {
                    "ens_id": int(ens_id),
                    "ii_sig": ii,
                    "seed": _process_noise_seed(base_seed, k, ens_id, ii),
                    "Lx_dim": np.sqrt(Lx * Ly),
                    "noise_dim": hdim,
                    "num_vars": total_vars,
                }
            )
            white = np.asarray(generate_enkf_field(**noise_kwargs), dtype=np.float64)
            white = white.reshape(-1)[:hdim]
            current = (
                alpha * previous[start:stop]
                + np.sqrt(max(0.0, 1.0 - alpha**2)) * white
            )
            updated_noise[start:stop] = current
            sigma = float(sig_q[ii]) if ii < len(sig_q) else 0.0
            noise_increment[start:stop] = np.sqrt(dt) * sigma * rho * current
    finally:
        np.random.set_state(old_rng_state)

    ensemble_vec[:state_block_size] += noise_increment
    sub_rank = int(icesee_kwargs.get("sub_rank", 0))
    if sub_rank == 0:
        _write_member_process_noise(
            icesee_kwargs, ens_id, updated_noise, k
        )
    return ensemble_vec


def advance_ensemble_member(
    ensemble_vec,
    ens_id,
    icesee_kwargs,
    indx_map,
    *,
    sub_rank=0,
):
    """Advance one member using the execution-mode-independent kernel.

    Execution modes 1 and 2 differ only in how member vectors are stored and
    transported.  The model call, state-vector packing, deterministic random
    stream, and process-noise update must remain identical.  Keeping those
    operations in one function prevents the two drivers from gradually
    developing different scientific results.
    """
    local_kwargs = dict(icesee_kwargs)
    local_kwargs.update(
        {
            "ens_id": int(ens_id),
            "sub_rank": int(sub_rank),
        }
    )
    member = np.asarray(ensemble_vec, dtype=np.float64).copy()
    member_before = member.copy() if local_kwargs.get("execution_parity_trace") else None

    with forecast_member_context(local_kwargs, ens_id) as member_kwargs:
        member_kwargs["sub_rank"] = int(sub_rank)
        updated_state = icesee_kwargs["model_module"].forecast_step_single(
            ensemble=member,
            **member_kwargs,
        )

    for key, value in updated_state.items():
        member[indx_map[key]] = value

    member_after_model = (
        member.copy() if local_kwargs.get("execution_parity_trace") else None
    )

    # ``add_member_process_noise`` owns the schedule decision.  Calling it
    # unconditionally here keeps the two execution modes on exactly the same
    # branch and avoids duplicated schedule logic in their I/O wrappers.
    member = add_member_process_noise(member, ens_id, member_kwargs)
    if local_kwargs.get("execution_parity_trace") and int(sub_rank) == 0:
        trace_root = os.path.join(
            str(
                local_kwargs.get("_modelrun_datasets")
                or local_kwargs.get("data_path", "_modelrun_datasets")
            ),
            "_execution_parity_trace",
        )
        os.makedirs(trace_root, exist_ok=True)
        timestep = int(local_kwargs.get("k", 0))
        trace_path = os.path.join(
            trace_root,
            f"member_{int(ens_id):06d}_step_{timestep:08d}.npz",
        )
        metadata = {
            "execution_mode": int(local_kwargs.get("execution_mode", -1)),
            "member": int(ens_id),
            "timestep": timestep,
            "dt": float(local_kwargs.get("dt", np.nan)),
            "sigma_96": float(local_kwargs.get("sigma_96", np.nan)),
            "beta_96": float(local_kwargs.get("beta_96", np.nan)),
            "rho_96": float(local_kwargs.get("rho_96", np.nan)),
            "process_noise_due": bool(
                process_noise_is_due(local_kwargs, timestep)
            ),
        }
        np.savez(
            trace_path,
            member_before=member_before,
            member_after_model=member_after_model,
            member_after_noise=member,
            metadata=np.asarray(json.dumps(metadata)),
        )
    return member

def parallel_forecast_step_default_run(**icesee_kwargs):
    import h5py
    import numpy as np
    from mpi4py import MPI

    Nens = icesee_kwargs["Nens"]
    nd = icesee_kwargs.get("nd")

    rounds = icesee_kwargs.get("rounds")
    color = icesee_kwargs.get("color", 0)
    subcomm_size_min = icesee_kwargs.get("subcomm_size_min", 1)

    subcomm = icesee_kwargs.get("subcomm", MPI.COMM_SELF)
    comm_world = icesee_kwargs.get("comm_world", MPI.COMM_WORLD)

    rank_world = comm_world.Get_rank()
    sub_rank = subcomm.Get_rank()

    model_module = icesee_kwargs["model_module"]
    k = icesee_kwargs.get("k", 0)

    _modelrun_datasets = icesee_kwargs.get("_modelrun_datasets", "_modelrun_datasets")
    input_file = f"{_modelrun_datasets}/icesee_ensemble_data.h5"

    partitioned_io = icesee_kwargs.get("partitioned_io_flag", False)  # NEW

    alpha = icesee_kwargs.get("alpha", 0.0)
    rho = icesee_kwargs.get("rho", 1.0)
    dt = icesee_kwargs.get("dt", 1.0)
    Lx = icesee_kwargs.get("Lx", 1.0)
    Ly = icesee_kwargs.get("Ly", 1.0)

    # Retained for backward-compatible return dictionaries only.  Process
    # noise is member-keyed and restartable; it must never depend on one
    # shared, rank-local ``noise`` array.
    noise = icesee_kwargs.get("noise", None)

    time_forecast_ensemble_generation = icesee_kwargs.get("time_forecast_ensemble_generation", 0.0)
    time_forecast_noise_generation = icesee_kwargs.get("time_forecast_noise_generation", 0.0)
    time_forecast_ensemble_mean_generation = icesee_kwargs.get("time_forecast_ensemble_mean_generation", 0.0)
    time_forecast_file_writing = icesee_kwargs.get("time_forecast_file_writing", 0.0)

    root_comm = comm_world.Split(
        color=0 if sub_rank == 0 else MPI.UNDEFINED,
        key=rank_world,
    )
    is_root_leader = root_comm is not MPI.COMM_NULL
    root_rank = root_comm.Get_rank() if is_root_leader else None
    root_is_world_root = is_root_leader and rank_world == 0

    vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)

    def _add_process_noise(ensemble_vec, ens_id, local_kwargs):
        nonlocal time_forecast_noise_generation
        _t_noise = MPI.Wtime()
        local_kwargs.update({"sub_rank": sub_rank})
        ensemble_vec = add_member_process_noise(
            ensemble_vec, ens_id, local_kwargs
        )
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
        ens_file = open_ensemble_file(input_file, nd, Nens, icesee_kwargs.get("nt", icesee_kwargs["nt"]), comm_world)
        ens_dset = ens_file["ensemble"]

    if Nens >= comm_world.Get_size():
        for round_id in range(rounds):
            ens_id = color + round_id * subcomm_size_min
            local_vec = None

            if ens_id < Nens:
                local_kwargs = dict(icesee_kwargs)
                local_kwargs.update({"ens_id": ens_id, "comm": subcomm})
                subcomm.Barrier()

                with h5py.File(input_file, "r") as f:
                    ensemble_local = f["ensemble"][:, ens_id, k].astype(np.float64, copy=True)

                _t_noise = MPI.Wtime()
                ensemble_local = advance_ensemble_member(
                    ensemble_local,
                    ens_id,
                    local_kwargs,
                    indx_map,
                    sub_rank=sub_rank,
                )
                if process_noise_is_due(icesee_kwargs, k):
                    time_forecast_noise_generation += MPI.Wtime() - _t_noise

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
            local_kwargs = dict(icesee_kwargs)
            local_kwargs.update({"ens_id": ens_id, "comm": subcomm})
            subcomm.Barrier()

            with h5py.File(input_file, "r") as f:
                ensemble_local = f["ensemble"][:, ens_id, k].astype(np.float64, copy=True)

            _t_noise = MPI.Wtime()
            ensemble_local = advance_ensemble_member(
                ensemble_local,
                ens_id,
                local_kwargs,
                indx_map,
                sub_rank=sub_rank,
            )
            if process_noise_is_due(icesee_kwargs, k):
                time_forecast_noise_generation += MPI.Wtime() - _t_noise

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

    icesee_kwargs.update({
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

    return icesee_kwargs, ensemble_vec, shape_ens, ens_mean

def parallel_forecast_step_default_full_parallel_run(**icesee_kwargs):
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

    Nens = icesee_kwargs["Nens"]

    comm_world = icesee_kwargs.get("comm_world", MPI.COMM_WORLD)
    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()

    subcomm = icesee_kwargs.get("subcomm", MPI.COMM_SELF)
    sub_rank = subcomm.Get_rank()

    color = icesee_kwargs.get("color", 0)
    rounds = icesee_kwargs.get("rounds", 1)
    subcomm_size_min = icesee_kwargs.get("subcomm_size_min", 1)

    model_module = icesee_kwargs.get("model_module")
    enkf_parallel_io = icesee_kwargs.get("enkf_parallel_io")

    if enkf_parallel_io is None:
        raise ValueError("parallel_forecast_step_default_full_parallel_run requires enkf_parallel_io.")

    nd = icesee_kwargs.get("nd")
    nt = icesee_kwargs.get("nt")
    k = icesee_kwargs.get("k", 0)

    alpha = icesee_kwargs.get("alpha", 0.0)
    rho = icesee_kwargs.get("rho", 1.0)
    dt = icesee_kwargs.get("dt", 1.0)
    Lx = icesee_kwargs.get("Lx", 1.0)
    Ly = icesee_kwargs.get("Ly", 1.0)

    # Backward-compatible metadata only; stochastic state is maintained per
    # member by ``add_member_process_noise``.
    noise = icesee_kwargs.get("noise", None)

    time_forecast_ensemble_generation = icesee_kwargs.get(
        "time_forecast_ensemble_generation", 0.0
    )
    time_forecast_noise_generation = icesee_kwargs.get(
        "time_forecast_noise_generation", 0.0
    )
    time_forecast_ensemble_mean_generation = icesee_kwargs.get(
        "time_forecast_ensemble_mean_generation", 0.0
    )
    time_forecast_file_writing = icesee_kwargs.get(
        "time_forecast_file_writing", 0.0
    )

    vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)

    def _add_process_noise(ensemble_vec, ens_id, local_kwargs):
        nonlocal time_forecast_noise_generation
        _t_noise = MPI.Wtime()
        local_kwargs.update({"sub_rank": sub_rank})
        ensemble_vec = add_member_process_noise(
            ensemble_vec, ens_id, local_kwargs
        )
        time_forecast_noise_generation += MPI.Wtime() - _t_noise
        return ensemble_vec

    def _run_one_ensemble(ens_id):
        nonlocal time_forecast_file_writing
        nonlocal time_forecast_noise_generation

        if ens_id < 0 or ens_id >= Nens:
            return

        local_kwargs = dict(icesee_kwargs)
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

        _t_noise = MPI.Wtime()
        ensemble_vec = advance_ensemble_member(
            ensemble_vec,
            ens_id,
            local_kwargs,
            indx_map,
            sub_rank=sub_rank,
        )
        if process_noise_is_due(icesee_kwargs, k):
            time_forecast_noise_generation += MPI.Wtime() - _t_noise

        # To avoid duplicate writes, only the subcommunicator root writes
        # the completed ensemble member.
        if sub_rank == 0:
            _t_write = MPI.Wtime()
            # History contains the initial condition (0) plus one state for
            # each of the ``nt`` model advances, exactly as in mode 1.
            write_k = k + 1
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

    km = icesee_kwargs.get("km", 0)
    obs_index = icesee_kwargs["obs_index"]

    if (km < icesee_kwargs["number_obs_instants"]) and (k == obs_index[km]):
        write_k = k + 1 if k < nt - 1 else k
        enkf_parallel_io.compute_forecast_mean_chunked_v2(
            write_k,
            flag="initial",
        )

    time_forecast_ensemble_mean_generation += MPI.Wtime() - _t_mean

    icesee_kwargs.update(
        {
            "time_forecast_ensemble_generation": time_forecast_ensemble_generation,
            "time_forecast_noise_generation": time_forecast_noise_generation,
            "time_forecast_ensemble_mean_generation": time_forecast_ensemble_mean_generation,
            "time_forecast_file_writing": time_forecast_file_writing,
            "shape_ens": shape_ens,
            "noise": noise,
        }
    )

    return icesee_kwargs



def parallel_forecast_step_even_distribution_run(**icesee_kwargs):
    """ Parallel run of the forecast step for each ensemble member.
        This function is designed to be used in a distributed environment where each rank processes a single ensemble member.
        It assumes that the number of ensemble members (Nens) is divisible by the size of the world communicator (size_world).
    """

    # unpack icesee_kwargs
    model_module = icesee_kwargs.get("model_module")
    comm_world = icesee_kwargs.get("comm_world", MPI.COMM_WORLD)
    rank_world = comm_world.Get_rank()
    Nens = icesee_kwargs.get("Nens", 1)  # Number of ensemble members
    nd = icesee_kwargs.get("nd", 1)  # Dimension of the state vector
    Q_err = icesee_kwargs.get("Q_err", np.eye(nd))  # Error covariance matrix
    state_block_size = icesee_kwargs.get("state_block_size", nd)  # Size of the state block
    size_world = comm_world.Get_size()  # Total number of ranks in the world communicato
    ensemble_vec = icesee_kwargs.get("ensemble_vec", np.zeros((nd, Nens), dtype=np.float64))  # Initialize ensemble vector
    ensemble_vec_mean = icesee_kwargs.get("ensemble_vec_mean", np.zeros((nd, icesee_kwargs.get("nt", icesee_kwargs["nt"]) + 1), dtype=np.float64))  # Initialize ensemble mean vector
    shape_ens = np.array(ensemble_vec.shape, dtype=np.int32)  # Shape of the ensemble vector
    ensemble_local = icesee_kwargs.get("ensemble_local", np.zeros((nd, Nens), dtype=np.float64))  # Local ensemble vector
    k = icesee_kwargs.get("k", 0)  # Time step index, default to 0 if not provided

    # check if Nens is divisible by size_world and greater or equal to size_world
    if Nens >= size_world and Nens % size_world == 0:
        for ens in range(ensemble_local.shape[1]):
            ensemble_local[:, ens] = model_module.forecast_step_single(ensemble=ensemble_local, **icesee_kwargs)
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

def parallel_forecast_step_squential_run(**icesee_kwargs):
    """ Squential run of the forecast step for each ensemble member.
        This function is designed to be used in a distributed environment where each rank processes a single ensemble member.
        #TODO: still under development, not fully tested.
    """

    # unpack icesee_kwargs
    model_module = icesee_kwargs.get("model_module")
    comm_world = icesee_kwargs.get("comm_world", MPI.COMM_WORLD)
    rank_world = comm_world.Get_rank()
    Nens = icesee_kwargs.get("Nens", 1)  # Number of ensemble members
    nd = icesee_kwargs.get("nd", 1)  # Dimension of the state vector
    Q_err = icesee_kwargs.get("Q_err", np.eye(nd))  # Error covariance matrix
    state_block_size = icesee_kwargs.get("state_block_size", nd)  # Size of the state block
    ensemble_vec = icesee_kwargs.get("ensemble_vec", np.zeros((nd, Nens), dtype=np.float64))  # Initialize ensemble vector
    ensemble_vec_mean = icesee_kwargs.get("ensemble_vec_mean", np.zeros((nd, icesee_kwargs.get("nt", icesee_kwargs["nt"]) + 1), dtype=np.float64))  # Initialize ensemble mean vector
    shape_ens = np.array(ensemble_vec.shape, dtype=np.int32)  # Shape of the ensemble vector
    ensemble_local = icesee_kwargs.get("ensemble_local", np.zeros((nd, Nens), dtype=np.float64))  # Local ensemble vector
    k = icesee_kwargs.get("k", 0)  # Time step index, default to 0 if not provided

    ensemble_col_stack = []
    for ens in range(Nens):
        comm_world.Barrier() # make sure all processors are in sync
        ensemble_vec[:,ens] = model_module.forecast_step_single(ens=ens, ensemble=ensemble_vec, nd=nd,  **icesee_kwargs)
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
    #         ensemble_vec_mean = np.empty((shape_ens[0],icesee_kwargs.get("nt",icesee_kwargs["nt"])+1), dtype=np.float64)
    #         ensemble_vec_mean[:,k] = np.mean(stack_ensemble_vec_mean, axis=1)
    #     else:
    #         ensemble_vec_mean = np.empty((shape_ens[0],icesee_kwargs.get("nt",icesee_kwargs["nt"])), dtype=np.float64)

    #     # broadcast the ensemble mean to all processors
    #     comm_world.Bcast([ensemble_vec_mean, MPI.DOUBLE], root=0)
    #     print(f"[ICESEE] Rank: {rank_world}, Ensemble mean shape: {ensemble_vec_mean.shape}")

    ensemble_vec_mean[:,k+1] = np.mean(ens_T[:nd,:], axis=1)
    # ensemble_vec_mean[:,k+1] = ParallelManager().compute_mean(ens_T[:nd,:], comm_world)


    return icesee_kwargs, ensemble_vec, ensemble_vec_mean, shape_ens
