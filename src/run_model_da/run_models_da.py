# ==============================================================================
# @desc: Entry point for running ICESEE data assimilation (EnKF variants) with:
#   • mode="full"    → Fully parallelized workflow: parallel batch file creation
#                      (serial/parallel), parallel I/O via both h5py and Zarr,
#                      memory-optimized routines, and a fully parallel analysis step.
#                      Intended for large models and datasets.
#   • mode="partial" → Partially parallelized workflow: parallel I/O but only rank 0
#                      writes; the analysis vector X5 is computed on rank 0 and
#                      broadcast to all ranks. State vectors are scattered for local
#                      dot products (StateVec · X5) and gathered for the analysis update.
#                      Suited to small/medium models.
#   • mode="serial"  → Single-process execution. Useful for small/experimental models
#                      and testing serial EnKF variants (DEnKF, EnKF, EnTKF, and EnRSKF,
#                      with adaptive localization.
#
# Notes:
#   - Parallel I/O supports both HDF5 (h5py/mpi) and Zarr backends.
#   - Exactly one mode must be chosen per run.
#
# @date: 2024-11-04
# @author: Brian Kyanjo
# ==============================================================================

# ==== ICESEE utility imports ========================================
from ICESEE.src.run_model_da.icesee_da_serial import (
    icesee_model_data_assimilation_serial,  # (development in progress)
)
from ICESEE.src.run_model_da.icesee_da_full_parallel import (
    icesee_model_data_assimilation_full_parallel,
)
from ICESEE.src.run_model_da.icesee_da_partial_parallel import (
    icesee_model_data_assimilation_partial_parallel,
)

import traceback
import sys

# ======================== Run model with EnKF ========================
def icesee_model_data_assimilation(**model_kwargs):
    """
    Run ICESEE model-data assimilation using the Ensemble Kalman Filter (and variants).

    Parameters
    ----------
    **model_kwargs :
        Arbitrary keyword arguments passed through to the selected runner.
        Expected to contain:
          params : dict
              Must include:
                - 'mode': str
                    One of {"serial", "partial", "full"}.
                    Defaults to "partial" if not provided.

    Raises
    ------
    ValueError
        If 'mode' is not one of {"serial", "partial", "full"}.
    """

    params = model_kwargs.get("params", {})
    mode = params.get("mode", "partial")

    if isinstance(mode, dict):
        mode = next((k for k, v in mode.items() if v), "partial")


    if mode not in {"serial", "partial", "full"}:
        raise ValueError(
            f"Invalid mode '{mode}'. Must be one of: 'serial', 'partial', or 'full'."
        )

    # Dispatch map
    runners = {
        "serial": icesee_model_data_assimilation_serial,
        "partial": icesee_model_data_assimilation_partial_parallel,
        "full": icesee_model_data_assimilation_full_parallel,
    }

    runner = runners[mode]

    try:
        return runner(**model_kwargs)
    except Exception:
        print(f"[ICESEE] Error in {mode} run mode:")
        tb_str = "".join(traceback.format_exception(*sys.exc_info()))
        print(f"Traceback details:\n{tb_str}")
        # raise  # optionally re-raise so the calling code can handle it