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


    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()


    if params["even_distribution"] or (params["default_run"] and size_world <= params["Nens"]):
        if params["even_distribution"]:
            model_kwargs.update({'rank': rank_world, 'color': color, 'comm': comm_world})
        else:
            model_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})

        dim_list = comm_world.allgather(params["nd"])
        # print(f"[ICESEE] Dim list: {dim_list}")
        # save model_nprocs before update if rank_world == 0
        # model_nprocs = params.get("model_nprocs", 1)

        
        if rank_world == 0:
            
            model_kwargs.update({'ens_id': rank_world})
            # model_kwargs.update({'model_nprocs': (model_nprocs * size_world) - size_world}) # update the model_nprocs to include all processors for the external model run

            if model_kwargs.get("generate_true_state", True):
                print("[ICESEE] Generating true state ...")
                # dim_list = np.tile(params["nd"],size_world) # all processors have the same dimension
                model_kwargs.update({"global_shape": params["nd"], "dim_list": dim_list})
                statevec_true = np.zeros([params["nd"], model_kwargs.get("nt",params["nt"]) + 1])
                model_kwargs.update({"statevec_true": statevec_true})
                updated_true_state = model_module.generate_true_state(**model_kwargs)

                # unpack the dictionaery
                vecs, indx_map, dim_per_proc = icesee_get_index(statevec_true, **model_kwargs)
                ensemble_true_state = np.zeros_like(statevec_true)
                for key, value in updated_true_state.items():
                    ensemble_true_state[indx_map[key], :] = value

            if model_kwargs.get("generate_nurged_state",True):
                print("[ICESEE] Generating nurged state ...")
                model_kwargs.update({"statevec_nurged": np.zeros([params["nd"], model_kwargs.get("nt",params["nt"]) + 1])})
                ensemble_nurged_state = model_module.generate_nurged_state(**model_kwargs)

            # Write data to file
            if model_kwargs.get("generate_true_state",True) or model_kwargs.get("generate_nurged_state",True):
                with h5py.File(_true_nurged, "w") as f:
                    if model_kwargs.get("generate_true_state"):
                        f.create_dataset("true_state", data=ensemble_true_state)
                    if model_kwargs.get("generate_nurged_state"):
                        f.create_dataset("nurged_state", data=ensemble_nurged_state)

            # clean memory 
            if model_kwargs.get("generate_true_state",True):
                del updated_true_state
            if model_kwargs.get("generate_nurged_state",True):
                del ensemble_nurged_state
            gc.collect()

        else:
            pass
            
        comm_world.Barrier()
        
        # -- write both the true and nurged states to file --
        data_shape = (params["nd"], model_kwargs.get("nt",params["nt"]) + 1)
        
        model_kwargs.update({"dim_list": dim_list})

        # update model_nprocs back to the original value before proceeding to the # next step
        # model_kwargs.update({'model_nprocs': model_nprocs})

    else:
        # --- Generate True and Nurged States ---

        if params["default_run"] and size_world > params["Nens"]:
            model_kwargs.update({'rank': sub_rank, 'color': color, 'comm': subcomm})
            model_kwargs.update({'ens_id': color}) # Nens = color
            # gather all the vector dimensions from all processors
            dim_list = subcomm.allgather(params["nd"])
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
            dim_list = comm_world.allgather(params["nd"])
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
        