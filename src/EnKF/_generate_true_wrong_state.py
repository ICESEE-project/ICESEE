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
        nd   = int(model_kwargs.get("nd", params["nd"]))
        ntp1 = int(model_kwargs.get("nt", params["nt"]) + 1)

        with h5py.File(_true_nurged, "w") as f:
            d_true   = None
            d_nurged = None

            if model_kwargs.get("generate_true_state", True):
                print("[ICESEE] Generating true state (Serial mode) ...")
                d_true = f.create_dataset(
                    "true_state",
                    shape=(nd, ntp1),
                    dtype="f8",
                    chunks=chunk_size,
                )
                model_kwargs.update({"statevec_true": d_true})
                updated_state = model_module.generate_true_state(**model_kwargs)
                vecs, indx_map, dim_per_proc = icesee_get_index(**model_kwargs)
                for key, value in updated_state.items():
                    d_true[indx_map[key], :] = value

            if model_kwargs.get("generate_nurged_state", True):
                print("[ICESEE] Generating nurged state (Serial mode) ...")
                d_nurged = f.create_dataset(
                    "nurged_state",
                    shape=(nd, ntp1),
                    dtype="f8",
                    chunks=chunk_size,
                )
                model_kwargs.update({"statevec_nurged": d_nurged})
                d_nurged = model_module.generate_nurged_state(**model_kwargs)

            model_kwargs.update({"dim_list": dim_list})

    return model_kwargs