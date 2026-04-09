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
from mpi4py import MPI


from ICESEE.src.utils.tools import icesee_get_index

def generate_true_wrong_state(**model_kwargs):
    """"Generate true and nurged states for the ICESEE model.
    """

    # unpack model_kwargs
    params         = model_kwargs.get("params", {})
    model_module   = model_kwargs.get("model_module", None)
    comm_world     = model_kwargs.get("comm_world", MPI.COMM_WORLD)
    _true_nurged   = model_kwargs.get("true_nurged_file")
    color          = model_kwargs.get("color", 0)
    subcomm        = model_kwargs.get("subcomm", None)
    sub_rank       = model_kwargs.get("sub_rank", 0)
    data_path      = model_kwargs.get("data_path", "output/")
    chunk_size      = model_kwargs.get("chunk_size", 5000)
    icesee_path         = model_kwargs.get('icesee_path')


    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()


    if params["even_distribution"] or (params["default_run"] and size_world <= params["Nens"]):
        if params["even_distribution"]:
            model_kwargs.update({'rank': rank_world, 'color': color, 'comm': comm_world})
        else:
            model_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})

        dim_list = comm_world.allgather(model_kwargs.get("nd", params["nd"]))
        # print(f"[ICESEE] Dim list: {dim_list}")
        # save model_nprocs before update if rank_world == 0
        # model_nprocs = params.get("model_nprocs", 1)
        nd   = int(model_kwargs.get("nd", params["nd"]))
        ntp1 = int(model_kwargs.get("nt", params["nt"]) + 1)

        
        if rank_world == 0:
            
            model_kwargs.update({'ens_id': rank_world})
            model_kwargs.update({"global_shape": model_kwargs.get("nd", params["nd"]), "dim_list": dim_list})
            # model_kwargs.update({'model_nprocs': (model_nprocs * size_world) - size_world}) # update the model_nprocs to include all processors for the external model run
            # Define shape and dtype
            nd = model_kwargs.get("nd", params["nd"])
            npt1 = model_kwargs.get("nt", params["nt"]) + 1   # +1 as in your np.zeros

            if model_kwargs["joint_estimation"] or params["localization_flag"]:
                hdim = nd // params["total_state_param_vars"]
            else:
                hdim = nd // params["num_state_vars"]

            chunk_size = (hdim, 1)  # row-wise chunks, 1 time slice per chunk
            # chunk_size = (nd,1)

            gen_true   = bool(model_kwargs.get("generate_true_state", True))
            gen_nurged = bool(model_kwargs.get("generate_nurged_state", True))

            if not gen_true and not os.path.exists(_true_nurged):
                raise FileNotFoundError(f"{_true_nurged} not found, but generation is disabled.")

            # If neither is requested, do nothing (assume existing file/datasets are already there)
            if not (gen_true or gen_nurged):
                if model_kwargs.get("verbose", False):
                    print(f"[ICESEE] true/nurged generation disabled — reusing existing: {_true_nurged}")
            else:
                # Open existing file if present; otherwise create it.
                mode = "a" if os.path.exists(_true_nurged) else "w"
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
                        model_kwargs["statevec_true"] = d_true  # write target

                        out_true = model_module.generate_true_state(**model_kwargs)

                        # If function returns data, write it (else assume in-place write)
                        if out_true is not None:
                            vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
                            if isinstance(out_true, dict):
                                for key, value in out_true.items():
                                    d_true[indx_map[key], :] = value
                            else:
                                d_true[:, :] = out_true

                    # ---------- NURGED STATE ----------
                    if gen_nurged:
                        print("[ICESEE] Generating nurged state ...")
                        d_nurged = require_dataset("nurged_state", shape=(nd, ntp1), dtype="f8", chunks=chunk_size)
                        model_kwargs["statevec_nurged"] = d_nurged  # write target

                        out_nurged = model_module.generate_nurged_state(**model_kwargs)

                        # If function returns data, write it (else assume in-place write)
                        if out_nurged is not None:
                            vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
                            if isinstance(out_nurged, dict):
                                for key, value in out_nurged.items():
                                    d_nurged[indx_map[key], :] = value
                            else:
                                d_nurged[:, :] = out_nurged
            
        else:
            pass
            
        comm_world.Barrier()
        
        # -- write both the true and nurged states to file --
        data_shape = (model_kwargs.get("nd", params["nd"]), model_kwargs.get("nt",params["nt"]) + 1)
        
        model_kwargs.update({"dim_list": dim_list})

        # update model_nprocs back to the original value before proceeding to the # next step
        # model_kwargs.update({'model_nprocs': model_nprocs})

    else:
        # --- Generate True and Nurged States ---

        if params["default_run"] and size_world > params["Nens"]:
            model_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})
            model_kwargs.update({'ens_id': color}) # Nens = color
            # gather all the vector dimensions from all processors
            dim_list = subcomm.allgather(model_kwargs.get("nd", params["nd"]))
            global_shape = sum(dim_list)
            model_kwargs.update({"global_shape": global_shape, "dim_list": dim_list})

            if model_kwargs.get("generate_true_state", True):
                if rank_world == 0:
                    print("[ICESEE] Generating true state ...  ")
                # statevec_true = np.zeros([model_kwargs['dim_list'][sub_rank], model_kwargs.get("nt",params["nt"]) + 1])
                statevec_true = np.zeros([global_shape, model_kwargs.get("nt",params["nt"]) + 1])
                model_kwargs.update({"statevec_true": statevec_true})
                # generate the true state
                updated_true_state = model_module.generate_true_state(**model_kwargs)
                # ensemble_true_state = gather_and_broadcast_data_default_run(updated_true_state, subcomm, sub_rank, comm_world, rank_world, params)
                global_data = {key: subcomm.gather(data, root=0) for key, data in updated_true_state.items()}

                if sub_rank == 0:
                    for key in global_data:
                        # print(f"[ICESEE] Key: {key}, shape: {[arr.shape for arr in global_data[key]]}")
                        global_data[key] = np.vstack(global_data[key])

                    # stack all variables together into a single array
                    stacked = np.vstack([global_data[key] for key in updated_true_state.keys()])
                    shape_ = np.array(stacked.shape,dtype=np.int32)
                    hdim = stacked.shape[0] // params["total_state_param_vars"]
                    # print(f"[ICESEE] Shape of the true state: {stacked.shape} min ensemble true: {np.min(stacked[hdim,:])}, max ensemble true: {np.max(stacked[hdim,:])}")
                    if model_kwargs.get("generate_true_state"):
                        # write data to the file
                        with h5py.File(_true_nurged, "w", driver='mpio', comm=subcomm) as f:
                            f.create_dataset("true_state", data=stacked)
                        
                    hdim = stacked.shape[0] // params["total_state_param_vars"]

                else:
                    shape_ = np.empty(2,dtype=np.int32)
                    hdim = 0

                # broadcast the shape of the true state
                shape_ = comm_world.bcast(shape_, root=0)
                hdim   = comm_world.bcast(hdim, root=0)

                if sub_rank != 0:
                    stacked = np.empty(shape_,dtype=np.float64)
            

                # write data to the file instead for memory management



                # broadcast the true state
                # ensemble_true_state = comm_world.bcast(stacked, root=0)
                # hdim = ensemble_true_state.shape[0] // params["total_state_param_vars"]
            
            if model_kwargs.get("generate_nurged_state", True):
                if rank_world == 0:
                    print("[ICESEE] Generating nurged state ... ")
                # statevec_nurged = np.zeros([model_kwargs['dim_list'][sub_rank], model_kwargs.get("nt",params["nt"]) + 1])
                statevec_nurged = np.zeros([global_shape, model_kwargs.get("nt",params["nt"]) + 1])
                model_kwargs.update({"statevec_nurged": statevec_nurged})
                ensemble_nurged_state = model_module.generate_nurged_state(**model_kwargs)

                with h5py.File(_true_nurged, "a", driver='mpio', comm=comm_world) as f:
                    f.create_dataset("nurged_state", data=ensemble_nurged_state)
                del ensemble_nurged_state 

            comm_world.Barrier()
            # clean memory
            if model_kwargs.get("generate_true_state"):
                del updated_true_state
            gc.collect()

            # exit()
        elif params["sequential_run"]:
            # gather all the vector dimensions from all processors
            dim_list = comm_world.allgather(model_kwargs.get("nd", params["nd"]))
            global_shape = sum(dim_list)
            model_kwargs.update({"global_shape": global_shape, "dim_list": dim_list})
            statevec_true = np.zeros([model_kwargs["global_shape"], model_kwargs.get("nt",params["nt"]) + 1])
            model_kwargs.update({"statevec_true": statevec_true})
            # generate the true state
            ensemble_true_state = model_module.generate_true_state(**model_kwargs)

            # generate the nurged state
            statevec_nurged = np.zeros([model_kwargs["global_shape"], model_kwargs.get("nt",params["nt"]) + 1])
            model_kwargs.update({"statevec_nurged": statevec_nurged})
            ensemble_nurged_state = model_module.generate_nurged_state(**model_kwargs)

    # return new and updated model_kwargs
    # model_kwargs.update({"dim_list": dim_list, "global_shape": global_shape})

    return model_kwargs
        

def generate_true_wrong_state_full_parallel(**model_kwargs):
    """"Generate true and nurged states for the ICESEE model.
    """

    # unpack model_kwargs
    params         = model_kwargs.get("params", {})
    model_module   = model_kwargs.get("model_module", None)
    comm_world     = model_kwargs.get("comm_world", MPI.COMM_WORLD)
    _true_nurged   = model_kwargs.get("true_nurged_file")
    color          = model_kwargs.get("color", 0)
    subcomm        = model_kwargs.get("subcomm", None)
    sub_rank       = model_kwargs.get("sub_rank", 0)


    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()


    if params["even_distribution"] or (params["default_run"] and size_world <= params["Nens"]):
        if params["even_distribution"]:
            model_kwargs.update({'rank': rank_world, 'color': color, 'comm': comm_world})
        else:
            model_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})

        dim_list = comm_world.allgather(model_kwargs.get("nd", params["nd"]))
        # print(f"[ICESEE] Dim list: {dim_list}")
        # save model_nprocs before update if rank_world == 0
        # model_nprocs = params.get("model_nprocs", 1)

        chunk_size = (hdim, 1)  # row-wise chunks, 1 time slice per chunk
        # chunk_size = (nd,1)
        nd   = int(model_kwargs.get("nd", params["nd"]))
        ntp1 = int(model_kwargs.get("nt", params["nt"]) + 1)

        
        if rank_world == 0:
            
            model_kwargs.update({'ens_id': rank_world})
            model_kwargs.update({"global_shape": model_kwargs.get("nd", params["nd"]), "dim_list": dim_list})

            gen_true   = bool(model_kwargs.get("generate_true_state", True))
            gen_nurged = bool(model_kwargs.get("generate_nurged_state", True))

            if not gen_true and not os.path.exists(_true_nurged):
                raise FileNotFoundError(f"{_true_nurged} not found, but generation is disabled.")

            # If neither is requested, do nothing (assume existing file/datasets are already there)
            if not (gen_true or gen_nurged):
                if model_kwargs.get("verbose", False):
                    print(f"[ICESEE] true/nurged generation disabled — reusing existing: {_true_nurged}")
            else:
                # Open existing file if present; otherwise create it.
                mode = "a" if os.path.exists(_true_nurged) else "w"
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
                        model_kwargs["statevec_true"] = d_true  # write target

                        out_true = model_module.generate_true_state(**model_kwargs)

                        # If function returns data, write it (else assume in-place write)
                        if out_true is not None:
                            vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
                            if isinstance(out_true, dict):
                                for key, value in out_true.items():
                                    d_true[indx_map[key], :] = value
                            else:
                                d_true[:, :] = out_true

                    # ---------- NURGED STATE ----------
                    if gen_nurged:
                        print("[ICESEE] Generating nurged state ...")
                        d_nurged = require_dataset("nurged_state", shape=(nd, ntp1), dtype="f8", chunks=chunk_size)
                        model_kwargs["statevec_nurged"] = d_nurged  # write target

                        out_nurged = model_module.generate_nurged_state(**model_kwargs)

                        # If function returns data, write it (else assume in-place write)
                        if out_nurged is not None:
                            vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
                            if isinstance(out_nurged, dict):
                                for key, value in out_nurged.items():
                                    d_nurged[indx_map[key], :] = value
                            else:
                                d_nurged[:, :] = out_nurged

        else:
            pass
            
        comm_world.Barrier()
        
        # -- write both the true and nurged states to file --
        data_shape = (model_kwargs.get("nd", params["nd"]), model_kwargs.get("nt",params["nt"]) + 1)
        
        model_kwargs.update({"dim_list": dim_list})

        # update model_nprocs back to the original value before proceeding to the # next step
        # model_kwargs.update({'model_nprocs': model_nprocs})

    else:
        # --- Generate True and Nurged States ---

        if params["default_run"] and size_world > params["Nens"]:
            model_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})
            model_kwargs.update({'ens_id': color}) # Nens = color
            # gather all the vector dimensions from all processors
            dim_list = subcomm.allgather(model_kwargs.get("nd", params["nd"]))
            global_shape = sum(dim_list)
            model_kwargs.update({"global_shape": global_shape, "dim_list": dim_list})

            if model_kwargs.get("generate_true_state", True):
                if rank_world == 0:
                    print("[ICESEE] Generating true state ...  ")
                # statevec_true = np.zeros([model_kwargs['dim_list'][sub_rank], model_kwargs.get("nt",params["nt"]) + 1])
                statevec_true = np.zeros([global_shape, model_kwargs.get("nt",params["nt"]) + 1])
                model_kwargs.update({"statevec_true": statevec_true})
                # generate the true state
                updated_true_state = model_module.generate_true_state(**model_kwargs)
                # ensemble_true_state = gather_and_broadcast_data_default_run(updated_true_state, subcomm, sub_rank, comm_world, rank_world, params)
                global_data = {key: subcomm.gather(data, root=0) for key, data in updated_true_state.items()}

                if sub_rank == 0:
                    for key in global_data:
                        # print(f"[ICESEE] Key: {key}, shape: {[arr.shape for arr in global_data[key]]}")
                        global_data[key] = np.vstack(global_data[key])

                    # stack all variables together into a single array
                    stacked = np.vstack([global_data[key] for key in updated_true_state.keys()])
                    shape_ = np.array(stacked.shape,dtype=np.int32)
                    hdim = stacked.shape[0] // params["total_state_param_vars"]
                    # print(f"[ICESEE] Shape of the true state: {stacked.shape} min ensemble true: {np.min(stacked[hdim,:])}, max ensemble true: {np.max(stacked[hdim,:])}")
                    if model_kwargs.get("generate_true_state"):
                        # write data to the file
                        with h5py.File(_true_nurged, "w", driver='mpio', comm=subcomm) as f:
                            f.create_dataset("true_state", data=stacked)
                        
                    hdim = stacked.shape[0] // params["total_state_param_vars"]

                else:
                    shape_ = np.empty(2,dtype=np.int32)
                    hdim = 0

                # broadcast the shape of the true state
                shape_ = comm_world.bcast(shape_, root=0)
                hdim   = comm_world.bcast(hdim, root=0)

                if sub_rank != 0:
                    stacked = np.empty(shape_,dtype=np.float64)
            

                # write data to the file instead for memory management



                # broadcast the true state
                # ensemble_true_state = comm_world.bcast(stacked, root=0)
                # hdim = ensemble_true_state.shape[0] // params["total_state_param_vars"]
            
            if model_kwargs.get("generate_nurged_state", True):
                if rank_world == 0:
                    print("[ICESEE] Generating nurged state ... ")
                # statevec_nurged = np.zeros([model_kwargs['dim_list'][sub_rank], model_kwargs.get("nt",params["nt"]) + 1])
                statevec_nurged = np.zeros([global_shape, model_kwargs.get("nt",params["nt"]) + 1])
                model_kwargs.update({"statevec_nurged": statevec_nurged})
                ensemble_nurged_state = model_module.generate_nurged_state(**model_kwargs)

                with h5py.File(_true_nurged, "a", driver='mpio', comm=comm_world) as f:
                    f.create_dataset("nurged_state", data=ensemble_nurged_state)
                del ensemble_nurged_state 

            comm_world.Barrier()
            # clean memory
            if model_kwargs.get("generate_true_state"):
                del updated_true_state
            gc.collect()

            # exit()
        elif params["sequential_run"]:
            # gather all the vector dimensions from all processors
            dim_list = comm_world.allgather(model_kwargs.get("nd", params["nd"]))
            global_shape = sum(dim_list)
            model_kwargs.update({"global_shape": global_shape, "dim_list": dim_list})
            statevec_true = np.zeros([model_kwargs["global_shape"], model_kwargs.get("nt",params["nt"]) + 1])
            model_kwargs.update({"statevec_true": statevec_true})
            # generate the true state
            ensemble_true_state = model_module.generate_true_state(**model_kwargs)

            # generate the nurged state
            statevec_nurged = np.zeros([model_kwargs["global_shape"], model_kwargs.get("nt",params["nt"]) + 1])
            model_kwargs.update({"statevec_nurged": statevec_nurged})
            ensemble_nurged_state = model_module.generate_nurged_state(**model_kwargs)

    # return new and updated model_kwargs
    # model_kwargs.update({"dim_list": dim_list, "global_shape": global_shape})

    return model_kwargs
        