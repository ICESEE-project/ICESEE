# ==============================================================================
# @des: This file contains run functions for the ICESEE model to generate true and nurged states. Serial version
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
import psutil
def rss_gb():
    return psutil.Process(os.getpid()).memory_info().rss / 1e9

from ICESEE.src.utils.tools import icesee_get_index

def generate_true_wrong_state(**model_kwargs):
    """"Generate true and nurged states for the ICESEE model.
    """

    # unpack model_kwargs
    params         = model_kwargs.get("params", {})
    model_module   = model_kwargs.get("model_module", None)
    _true_nurged   = model_kwargs.get("true_nurged_file")
    color          = model_kwargs.get("color", 0)
    subcomm        = model_kwargs.get("subcomm", None)
    sub_rank       = model_kwargs.get("sub_rank", 0)
    data_path      = model_kwargs.get("data_path", "output/")
    chunk_size      = model_kwargs.get("chunk_size", 5000)
    icesee_path         = model_kwargs.get('icesee_path')

    # Set up serial MPI parameters
    rank_world = 0
    size_world = 1
    
        
    dim_list = [model_kwargs.get("nd", params["nd"])] * size_world
        
    if rank_world == 0:
        
        model_kwargs.update({'ens_id': rank_world})
        # model_kwargs.update({'model_nprocs': (model_nprocs * size_world) - size_world}) # update the model_nprocs to include all processors for the external model run
        # Define shape and dtype
        nd = model_kwargs.get("nd", params["nd"])
        nt = model_kwargs.get("nt", params["nt"]) + 1   # +1 as in your np.zeros
        # print(f"[ICESEE] Generating true and nurged states with shape ({nd}, {nt}) ...")

        if model_kwargs["joint_estimation"] or params["localization_flag"]:
            hdim = nd // params["total_state_param_vars"]
        else:
            hdim = nd // params["num_state_vars"]

        chunk_size = (hdim, 1)  # row-wise chunks, 1 time slice per chunk
        # chunk_size = (nd,1)

        if model_kwargs.get("generate_true_state", True):
            print("[ICESEE] Generating true state (partial parallel) ...")
            # dim_list = np.tile(model_kwargs.get("nd", params["nd"]),size_world) # all processors have the same dimension
            model_kwargs.update({"global_shape": model_kwargs.get("nd", params["nd"]), "dim_list": dim_list})
            # statevec_true = np.zeros([model_kwargs.get("nd", params["nd"]), model_kwargs.get("nt",params["nt"]) + 1])

            # (Optional) remove old store if shape/chunks might have changed
            store_path = f"{icesee_path}/{data_path}/statevec_true.zarr"
            if os.path.exists(store_path):
                shutil.rmtree(store_path)

            statevec_true = zarr.zeros(
                shape=(nd, nt),
                chunks=chunk_size,   # row-wise chunks, 1 time slice per chunk
                dtype="f8", 
                store=store_path,
            )

            model_kwargs.update({"statevec_true": statevec_true})
            # updated_true_state = model_module.generate_true_state(**model_kwargs)
            # print("RSS after zarr init:", rss_gb())
            model_module.generate_true_state(**model_kwargs)
            # print("RSS after generate_true_state:", rss_gb())

        if model_kwargs.get("generate_nurged_state",True):
            print("[ICESEE] Generating nurged state ...")
            store_path = f"{icesee_path}/{data_path}/statevec_nurged.zarr"
            if os.path.exists(store_path):
                shutil.rmtree(store_path)

            statevec_nurged = zarr.zeros(
                shape=(nd, nt),
                chunks=chunk_size,   # row-wise chunks, 1 time slice per chunk
                dtype="f8", 
                store=store_path,
            )
            model_kwargs.update({"statevec_nurged": statevec_nurged})
            # print("RSS after zarr init nurged:", rss_gb())
            model_module.generate_nurged_state(**model_kwargs)
            # print("RSS after generate_nurged_state:", rss_gb())
            
        # Write data to file
        if model_kwargs.get("generate_true_state",True) or model_kwargs.get("generate_nurged_state",True):
            with h5py.File(_true_nurged, "w") as f:
                # print("RSS after opening h5:", rss_gb())
                if model_kwargs.get("generate_true_state"):
                    dtrue = f.create_dataset("true_state", shape=(nd, nt), dtype="f8", chunks=chunk_size)
                    # write in time slices (small memory)
                    for j in range(nt):
                        dtrue[:, j] = statevec_true[:, j]

                if model_kwargs.get("generate_nurged_state"):
                    dn = f.create_dataset("nurged_state", shape=(nd, nt), dtype="f8", chunks=chunk_size)
                    for j in range(nt):
                        dn[:, j] = statevec_nurged[:, j]
                # print("RSS after writing datasets:", rss_gb())

        # -- write both the true and nurged states to file --
        data_shape = (model_kwargs.get("nd", params["nd"]), model_kwargs.get("nt",params["nt"]) + 1)
        
        model_kwargs.update({"dim_list": dim_list})
        gc.collect()

        # update model_nprocs back to the original value before proceeding to the # next step
        # model_kwargs.update({'model_nprocs': model_nprocs})

    return model_kwargs