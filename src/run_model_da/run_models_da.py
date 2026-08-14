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

import importlib
import sys, traceback
from ICESEE.src.utils.icesee_context import normalize_icesee_kwargs

_MODE_TO_TARGET = {
    "serial":  ("ICESEE.src.run_model_da.icesee_da_serial",
                "icesee_model_data_assimilation_serial"),
    "partial": ("ICESEE.src.run_model_da.icesee_da_partial_parallel",
                "icesee_model_data_assimilation_partial_parallel"),
    "full":    ("ICESEE.src.run_model_da.icesee_da_full_parallel",
                "icesee_model_data_assimilation_full_parallel"),
}

def _resolve_mode(icesee_kwargs) -> str:
    """Accept 'mode' as str or dict {'serial':1,'partial':0,'full':0}. Default 'partial'."""
    if not isinstance(icesee_kwargs, dict):
        return "partial"
    mode = icesee_kwargs.get("mode", "partial")
    if isinstance(mode, dict):
        mode = next((k for k, v in mode.items() if v), "partial")
    return mode

# ======================== Run model with EnKF ========================
def icesee_model_data_assimilation(**icesee_kwargs):
    """
    Run ICESEE model-data assimilation using the Ensemble Kalman Filter (and variants).

    Keyword arguments form the single flat ICESEE runtime context. ``mode``
    may be ``"serial"``, ``"partial"``, or ``"full"``; the legacy Boolean
    mode mapping is also accepted at this public boundary.
    """
    icesee_kwargs = normalize_icesee_kwargs(icesee_kwargs)
    mode = _resolve_mode(icesee_kwargs)

    # ``batch_size`` is an I/O-window control, not an ensemble/time batching
    # heuristic.  In particular, mode 2 must honor a value of one or two;
    # silently expanding ``1`` to a sizeable fraction of ``nt`` can reserve a
    # large ensemble history before the first forecast and defeats bounded
    # large-data execution.
    batch_size = icesee_kwargs.get("batch_size")
    if batch_size is None:
        icesee_kwargs["batch_size"] = 2 if mode == "full" else 1
    else:
        icesee_kwargs["batch_size"] = max(1, int(batch_size))

    if mode not in _MODE_TO_TARGET:
        raise ValueError(
            f"Invalid mode '{mode}'. Must be one of: 'serial', 'partial', or 'full'."
        )

    module_name, func_name = _MODE_TO_TARGET[mode]

    try:
        # Lazy import only the selected runner
        mod = importlib.import_module(module_name)
        runner = getattr(mod, func_name)
        return runner(**icesee_kwargs)
    except Exception:
        print(f"[ICESEE] Error in {mode} run mode:")
        tb_str = "".join(traceback.format_exception(*sys.exc_info()))
        print(f"Traceback details:\n{tb_str}")
