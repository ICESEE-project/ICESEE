# ==============================================================================
# @des: This file contains run functions for the ICESEE model to generate true and nurged states. Serial version
# @date: 2025-07-30
# @author: Brian Kyanjo
# ==============================================================================

# --- import necessary libraries ---
import numpy as np
import h5py
import gc

from ICESEE.src.utils.utils import UtilsFunctions
from ICESEE.src.utils.tools import icesee_get_index


def generate_synthetic_observations(**model_kwargs):
    """Generate synthetic observations for the ICESEE model.
    """

    # unpack model_kwargs
    params         = model_kwargs.get("params", {})
    model_module   = model_kwargs.get("model_module", None)
    _synthetic_obs = model_kwargs.get("synthetic_obs_file")
    _true_nurged   = model_kwargs.get("true_nurged_file")
    color          = model_kwargs.get("color", 0)
    subcomm        = model_kwargs.get("subcomm", None)
    sub_rank       = model_kwargs.get("sub_rank", 0)
    rank_world = 0
    size_world = 1

    if model_kwargs.get("generate_synthetic_obs", True):   
        if params["even_distribution"] or (params["default_run"] and size_world <= params["Nens"]):
            if rank_world == 0:
                # --- Synthetic Observations ---
                print("[ICESEE] Generating synthetic observations ...")
                with h5py.File(_true_nurged, "r") as f:
                    ensemble_true_state = f['true_state'][:]

                utils_funs = UtilsFunctions(params, ensemble_true_state)
                model_kwargs.update({"statevec_true": ensemble_true_state})
                hu_obs, error_R, model_kwargs['bed_mask_map'] = utils_funs._create_synthetic_observations(**model_kwargs)

                # observe or don't observe parameters.
                vecs, indx_map,_ = icesee_get_index(hu_obs, **model_kwargs)
                # check if model_kwargs['observe_params'] is empty
                if len(model_kwargs['observed_params']) == 0:
                    for key in model_kwargs['params_vec']:
                        hu_obs[indx_map[key],:] = 0.0
                        error_R[:,indx_map[key]] = 0.0
                else: 
                    for key in model_kwargs['params_vec']:
                        if key not in model_kwargs['observed_params']:
                            hu_obs[indx_map[key],:] = 0.0
                            error_R[:,indx_map[key]] = 0.0

                # -- write data to file
                with h5py.File(_synthetic_obs, 'w') as f:
                    f.create_dataset("hu_obs", data=hu_obs)
                    f.create_dataset("R", data=error_R)

                # --- clear memory
                del hu_obs
                del error_R
                gc.collect()

            else:
                pass

    return model_kwargs