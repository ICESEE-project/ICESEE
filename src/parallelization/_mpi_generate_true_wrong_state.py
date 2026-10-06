# ==============================================================================
# @des: This file contains run functions for the ICESEE model to generate true and nurged states.
# @date: 2025-07-30
# @author: Brian Kyanjo
# ==============================================================================

# --- import necessary libraries ---
import numpy as np
import h5py
import gc
import zarr
import os
import shutil
import time
from mpi4py import MPI


from ICESEE.src.utils.tools import icesee_get_index
from ICESEE.src.utils.performance import record_phase
from ICESEE.src.utils.state_ownership import (
    resolve_state_ownership,
    combine_member_state,
)


def _timed_phase(phase, generate, icesee_kwargs):
    """Call a model trajectory generator and record its time for the
    end-of-run performance summary."""
    start = time.perf_counter()
    result = generate(**icesee_kwargs)
    record_phase(phase, time.perf_counter() - start)
    return result


def generate_true_wrong_state(**icesee_kwargs):
    """"Generate true and nurged states for the ICESEE model.
    """
    # unpack icesee_kwargs
    model_module   = icesee_kwargs.get("model_module", None)
    comm_world     = icesee_kwargs.get("comm_world", MPI.COMM_WORLD)
    _true_nurged   = icesee_kwargs.get("true_nurged_file")
    color          = icesee_kwargs.get("color", 0)
    subcomm        = icesee_kwargs.get("subcomm", None)
    sub_rank       = icesee_kwargs.get("sub_rank", 0)
    data_path      = icesee_kwargs.get("data_path", "output/")
    chunk_size      = icesee_kwargs.get("chunk_size", 5000)
    icesee_path         = icesee_kwargs.get('icesee_path')


    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()


    # ranks_per_model == 1 (the centralized resource plan's resolved
    # value, not a raw size_world/Nens comparison -- see resource_plan.py)
    # is the actual condition for "no ICESEE-managed multi-rank model
    # communicator, rank 0 alone generates the single shared true/nurged
    # trajectory". size_world <= Nens is only equivalent to this for
    # configurations that never request an explicit ranks_per_model; an
    # explicit request (e.g. P=4,Nens=8,ranks_per_model=2) can have
    # ranks_per_model > 1 even though size_world <= Nens, which the old
    # condition could not express.
    if icesee_kwargs["even_distribution"] or (icesee_kwargs["default_run"] and icesee_kwargs.get("ranks_per_model", 1) == 1):
        if icesee_kwargs["even_distribution"]:
            icesee_kwargs.update({'rank': rank_world, 'color': color, 'comm': comm_world})
        else:
            icesee_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})

        dim_list = comm_world.allgather(icesee_kwargs.get("nd", icesee_kwargs["nd"]))
        # print(f"[ICESEE] Dim list: {dim_list}")
        # save model_nprocs before update if rank_world == 0
        # model_nprocs = icesee_kwargs.get("model_nprocs", 1)
        nd   = int(icesee_kwargs.get("nd", icesee_kwargs["nd"]))
        ntp1 = 1 if icesee_kwargs.get("initial_state_only", False) else int(
            icesee_kwargs.get("nt", icesee_kwargs["nt"]) + 1
        )


        if rank_world == 0:

            icesee_kwargs.update({'ens_id': rank_world})
            icesee_kwargs.update({"global_shape": icesee_kwargs.get("nd", icesee_kwargs["nd"]), "dim_list": dim_list})
            # icesee_kwargs.update({'model_nprocs': (model_nprocs * size_world) - size_world}) # update the model_nprocs to include all processors for the external model run
            # Define shape and dtype
            nd = icesee_kwargs.get("nd", icesee_kwargs["nd"])
            npt1 = icesee_kwargs.get("nt", icesee_kwargs["nt"]) + 1   # +1 as in your np.zeros

            if icesee_kwargs["joint_estimation"] or icesee_kwargs["localization_flag"]:
                hdim = nd // icesee_kwargs["total_state_param_vars"]
            else:
                hdim = nd // icesee_kwargs["num_state_vars"]

            chunk_size = (hdim, 1)  # row-wise chunks, 1 time slice per chunk
            # chunk_size = (nd,1)

            gen_true   = bool(icesee_kwargs.get("generate_true_state", True))
            gen_nurged = bool(icesee_kwargs.get("generate_nurged_state", True))

            if not gen_true and not os.path.exists(_true_nurged):
                raise FileNotFoundError(f"{_true_nurged} not found, but generation is disabled.")

            # If neither is requested, do nothing (assume existing file/datasets are already there)
            if not (gen_true or gen_nurged):
                if icesee_kwargs.get("verbose", False):
                    print(f"[ICESEE] true/nurged generation disabled — reusing existing: {_true_nurged}")
            else:
                # A fresh run must not append to a stale or partially written
                # initialization file.  In particular, an interrupted HDF5
                # create may leave only the user block/header on disk, which
                # exists but cannot be opened in append mode.
                force_fresh = bool(
                    icesee_kwargs.get(
                        "force_fresh_start",
                        icesee_kwargs.get("force_fresh_start", False),
                    )
                )
                mode = "w" if force_fresh or not os.path.exists(_true_nurged) else "a"
                with h5py.File(_true_nurged, mode) as f:
                    # Helper: create or replace dataset safely if shape mismatch
                    def require_dataset(name: str, shape, dtype="f8", chunks=None):
                        if name in f:
                            d = f[name]
                            if d.shape != tuple(shape) or d.dtype != np.dtype(dtype):
                                # Replace only this dataset, not the whole file
                                del f[name]
                                d = f.create_dataset(name, shape=shape, dtype=dtype, chunks=chunks)
                        else:
                            d = f.create_dataset(name, shape=shape, dtype=dtype, chunks=chunks)
                        return d

                    # ---------- TRUE STATE ----------
                    if gen_true:
                        print("[ICESEE] Generating true state ...")
                        d_true = require_dataset("true_state", shape=(nd, ntp1), dtype="f8", chunks=chunk_size)
                        icesee_kwargs["statevec_true"] = d_true  # write target

                        out_true = _timed_phase("truth_generation", model_module.generate_true_state, icesee_kwargs)

                        # If function returns data, write it (else assume in-place write)
                        if out_true is not None:
                            vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)
                            if isinstance(out_true, dict):
                                for key, value in out_true.items():
                                    d_true[indx_map[key], :] = value
                            else:
                                d_true[:, :] = out_true

                    # ---------- NURGED STATE ----------
                    if gen_nurged:
                        print("[ICESEE] Generating nurged state ...")
                        d_nurged = require_dataset("nurged_state", shape=(nd, ntp1), dtype="f8", chunks=chunk_size)
                        icesee_kwargs["statevec_nurged"] = d_nurged  # write target

                        out_nurged = _timed_phase("wrong_reference_generation", model_module.generate_nurged_state, icesee_kwargs)

                        # If function returns data, write it (else assume in-place write)
                        if out_nurged is not None:
                            vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)
                            if isinstance(out_nurged, dict):
                                for key, value in out_nurged.items():
                                    d_nurged[indx_map[key], :] = value
                            else:
                                d_nurged[:, :] = out_nurged

        else:
            pass

        comm_world.Barrier()

        # -- write both the true and nurged states to file --
        data_shape = (icesee_kwargs.get("nd", icesee_kwargs["nd"]), icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1)

        icesee_kwargs.update({"dim_list": dim_list})

        # update model_nprocs back to the original value before proceeding to the # next step
        # icesee_kwargs.update({'model_nprocs': model_nprocs})

    else:
        # --- Generate True and Nurged States ---

        if icesee_kwargs["default_run"] and icesee_kwargs.get("ranks_per_model", 1) > 1:
            # A spare rank (color is None, subcomm is MPI.COMM_NULL --
            # resource_plan.py) belongs to no model group and must never
            # touch subcomm-level collectives (allgather/gather below), but
            # it still must reach the same comm_world.Barrier() every
            # active rank in this branch reaches, so the branch condition
            # above intentionally does not exclude it -- only the
            # subcomm-dependent body does.
            updated_true_state = None
            if color is not None:
                icesee_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})
                icesee_kwargs.update({'ens_id': color}) # Nens = color
                # A member's model communicator (subcomm) may hold either a
                # replicated full state (e.g. Lorenz-96: every rank reports the
                # same nd) or a distributed partition (e.g. Icepack: nd is a
                # genuine per-rank dof count). Summing/concatenating a
                # replicated nd across the subcomm double(+)-counts it -- see
                # src/utils/state_ownership.py.
                ownership = resolve_state_ownership(icesee_kwargs, subcomm)
                global_shape = ownership.global_size
                if ownership.distribution == "distributed":
                    dim_list = subcomm.allgather(icesee_kwargs.get("nd", icesee_kwargs["nd"]))
                else:
                    dim_list = [ownership.local_size]
                icesee_kwargs.update({"global_shape": global_shape, "dim_list": dim_list})

                if rank_world == 0:
                    print(
                        f"[ICESEE][ownership] distribution={ownership.distribution} "
                        f"local_size={ownership.local_size} "
                        f"global_size={ownership.global_size} "
                        f"model_comm_size={subcomm.Get_size()}",
                        flush=True,
                    )

                # The true/nurged trajectory is a single reference the whole
                # ensemble is compared against, not a per-member quantity --
                # exactly as the Nens >= size_world branch above already treats
                # it (gated to rank_world == 0 only). Every rank previously
                # reached this point with ens_id = color and independently
                # wrote the same shared _true_nurged file, racing every other
                # color's subcomm to create the same dataset name (confirmed:
                # this raised "Unable to synchronously create dataset (name
                # already exists)" once the file-open collective-mismatch and
                # state-size bugs above no longer masked it). Only color 0's
                # subcomm generates and writes it here; every other color skips
                # straight to the barrier and reads the shared file later, like
                # color 0's own non-root ranks already do.
                if color == 0 and icesee_kwargs.get("generate_true_state", True):
                    if rank_world == 0:
                        print("[ICESEE] Generating true state ...  ")
                    # statevec_true = np.zeros([icesee_kwargs['dim_list'][sub_rank], icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1])
                    statevec_true = np.zeros([global_shape, icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1])
                    icesee_kwargs.update({"statevec_true": statevec_true})
                    # generate the true state
                    updated_true_state = _timed_phase("truth_generation", model_module.generate_true_state, icesee_kwargs)
                    # Replicated: every rank already computed the identical full
                    # state, so this is a no-op (no communication). Distributed:
                    # gathers+concatenates each rank's disjoint slice, exactly as
                    # before.
                    global_data = combine_member_state(
                        ownership, subcomm, updated_true_state, root=0
                    )

                    if sub_rank == 0:
                        # stack all variables together into a single array
                        stacked = np.vstack([global_data[key] for key in updated_true_state.keys()])
                        shape_ = np.array(stacked.shape,dtype=np.int32)
                        hdim = stacked.shape[0] // icesee_kwargs["total_state_param_vars"]
                        # print(f"[ICESEE] Shape of the true state: {stacked.shape} min ensemble true: {np.min(stacked[hdim,:])}, max ensemble true: {np.max(stacked[hdim,:])}")
                        if icesee_kwargs.get("generate_true_state"):
                            # ``stacked`` at this point exists only on sub_rank
                            # == 0 (combine_member_state already assembled it
                            # there); this write is a single-rank operation, not
                            # a subcomm-collective one. Opening with
                            # ``comm=subcomm`` was a pre-existing bug: h5py's
                            # mpio driver requires every rank in the given
                            # communicator to call File() collectively, but only
                            # sub_rank == 0 ever reaches this line -- whenever a
                            # member's subcomm had more than one rank, the other
                            # rank(s) never joined the open and this call hung
                            # forever. COMM_SELF makes the (already correct)
                            # single-rank intent explicit and collective-safe.
                            with h5py.File(_true_nurged, "w", driver='mpio', comm=MPI.COMM_SELF) as f:
                                f.create_dataset("true_state", data=stacked)
                    # The pre-existing shape_/hdim comm_world.bcast and the
                    # np.empty(shape_) placeholder it fed (previously here) were
                    # dead weight: neither is read by anything downstream (the
                    # actual array broadcast this was scaffolding for is the
                    # commented-out line below, never enabled), and now that
                    # this block only runs on color 0's subcomm, a comm_world
                    # collective here would be a genuine rank-subset mismatch
                    # against the other colors, which no longer enter this
                    # block at all. Removed rather than rescoped.
                    # ensemble_true_state = comm_world.bcast(stacked, root=0)

                if color == 0 and icesee_kwargs.get("generate_nurged_state", True):
                    if rank_world == 0:
                        print("[ICESEE] Generating nurged state ... ")
                    # statevec_nurged = np.zeros([icesee_kwargs['dim_list'][sub_rank], icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1])
                    statevec_nurged = np.zeros([global_shape, icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1])
                    icesee_kwargs.update({"statevec_nurged": statevec_nurged})
                    updated_nurged_state = _timed_phase("wrong_reference_generation", model_module.generate_nurged_state, icesee_kwargs)
                    # Same ownership-aware assembly as the true-state block
                    # above: a no-op for replicated state, gather+concatenate
                    # for distributed.
                    global_nurged_data = combine_member_state(
                        ownership, subcomm, updated_nurged_state, root=0
                    )

                    if sub_rank == 0:
                        stacked_nurged = np.vstack(
                            [global_nurged_data[key] for key in updated_nurged_state.keys()]
                        )
                        # COMM_SELF for the same reason as the true-state write
                        # above: this is a single-rank operation, and the
                        # previous ``comm=comm_world`` open additionally passed
                        # a dict (not an array) as ``data=``, which every rank
                        # would have hit as a TypeError had execution reached
                        # it without first hanging on the collective mismatch.
                        with h5py.File(_true_nurged, "a", driver='mpio', comm=MPI.COMM_SELF) as f:
                            f.create_dataset("nurged_state", data=stacked_nurged)
                    del updated_nurged_state

            comm_world.Barrier()
            # clean memory
            if updated_true_state is not None:
                del updated_true_state
            gc.collect()

            # exit()
        elif icesee_kwargs["sequential_run"]:
            # ``sequential_run`` has no subcomm; every world rank shares the
            # same top-level communicator, so ownership is resolved against
            # comm_world directly (see resolve_state_ownership's docstring
            # for why "distributed" -- summing a genuine local nd -- is the
            # correct default). This mode is CLI-only (``--sequential_run``)
            # and not exercised by any shipped example's params.yaml.
            ownership = resolve_state_ownership(icesee_kwargs, comm_world)
            global_shape = ownership.global_size
            if ownership.distribution == "distributed":
                dim_list = comm_world.allgather(icesee_kwargs.get("nd", icesee_kwargs["nd"]))
            else:
                dim_list = [ownership.local_size]
            icesee_kwargs.update({"global_shape": global_shape, "dim_list": dim_list})
            statevec_true = np.zeros([icesee_kwargs["global_shape"], icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1])
            icesee_kwargs.update({"statevec_true": statevec_true})
            # generate the true state
            ensemble_true_state = _timed_phase("truth_generation", model_module.generate_true_state, icesee_kwargs)

            # generate the nurged state
            statevec_nurged = np.zeros([icesee_kwargs["global_shape"], icesee_kwargs.get("nt",icesee_kwargs["nt"]) + 1])
            icesee_kwargs.update({"statevec_nurged": statevec_nurged})
            ensemble_nurged_state = _timed_phase("wrong_reference_generation", model_module.generate_nurged_state, icesee_kwargs)

    # return new and updated icesee_kwargs
    # icesee_kwargs.update({"dim_list": dim_list, "global_shape": global_shape})

    # Skipped phases are recorded with zero operations ("no events").
    if not icesee_kwargs.get("generate_true_state", True):
        record_phase("truth_generation", 0.0, operations=0)
    if not icesee_kwargs.get("generate_nurged_state", True):
        record_phase("wrong_reference_generation", 0.0, operations=0)
    return icesee_kwargs
