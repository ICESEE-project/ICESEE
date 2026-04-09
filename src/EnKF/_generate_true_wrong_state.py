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

            model_kwargs.update({"dim_list": dim_list})

    return model_kwargs