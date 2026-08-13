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


def generate_synthetic_observations(**icesee_kwargs):
    """Generate synthetic observations for the ICESEE model.
    """

    # unpack icesee_kwargs
    model_module   = icesee_kwargs.get("model_module", None)
    _synthetic_obs = icesee_kwargs.get("synthetic_obs_file")
    _true_nurged   = icesee_kwargs.get("true_nurged_file")
    color          = icesee_kwargs.get("color", 0)
    subcomm        = icesee_kwargs.get("subcomm", None)
    sub_rank       = icesee_kwargs.get("sub_rank", 0)
    rank_world = 0
    size_world = 1

    if icesee_kwargs.get("generate_synthetic_obs", True):
        if icesee_kwargs["even_distribution"] or (icesee_kwargs["default_run"] and size_world <= icesee_kwargs["Nens"]):
            if rank_world == 0:
                # --- Synthetic Observations ---
                print("[ICESEE] Generating synthetic observations ...")
                with h5py.File(_true_nurged, "r") as f:
                    ensemble_true_state = f['true_state'][:]

                icesee_kwargs.update({"statevec_true": ensemble_true_state})
                utils_funs = UtilsFunctions(
                    icesee_kwargs=icesee_kwargs,
                    ensemble=ensemble_true_state,
                )
                hu_obs, error_R, icesee_kwargs['bed_mask_map'], icesee_kwargs = utils_funs._create_synthetic_observations(**icesee_kwargs)

                # observe or don't observe parameters.
                vecs, indx_map,_ = icesee_get_index(hu_obs, **icesee_kwargs)
                all_observed =  icesee_kwargs['all_observed']
                # check if  icesee_kwargs['all_observed'] is empty
                if len( icesee_kwargs['all_observed']) == 0:
                    for key in icesee_kwargs['vec_inputs']:
                        hu_obs[indx_map[key],:] = 0.0
                        error_R[:,indx_map[key]] = 0.0
                else:
                    for key in icesee_kwargs['vec_inputs']:
                        if key not in icesee_kwargs['all_observed']:
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

    return icesee_kwargs
