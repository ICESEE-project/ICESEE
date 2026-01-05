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

_MODE_TO_TARGET = {
    "serial":  ("ICESEE.src.run_model_da.icesee_da_serial",
                "icesee_model_data_assimilation_serial"),
    "partial": ("ICESEE.src.run_model_da.icesee_da_partial_parallel",
                "icesee_model_data_assimilation_partial_parallel"),
    "full":    ("ICESEE.src.run_model_da.icesee_da_full_parallel",
                "icesee_model_data_assimilation_full_parallel"),
}

def _resolve_mode(params) -> str:
    """Accept 'mode' as str or dict {'serial':1,'partial':0,'full':0}. Default 'partial'."""
    if not isinstance(params, dict):
        return "partial"
    mode = params.get("mode", "partial")
    if isinstance(mode, dict):
        mode = next((k for k, v in mode.items() if v), "partial")
    return mode

# ======================== Run model with EnKF ========================
def icesee_model_data_assimilation(**model_kwargs):
    """
    Run ICESEE model-data assimilation using the Ensemble Kalman Filter (and variants).

    model_kwargs expected:
      params: dict with
        - mode: "serial" | "partial" | "full"   OR
        - mode: dict like {"serial":1,"partial":0,"full":0}
    """
    params = model_kwargs.get("params", {}) or {}
    mode = _resolve_mode(params)

    nt = model_kwargs.get('nt', params['nt'])
    batch_size = model_kwargs.get('batch_size', params.get('batch_size', None))
    params['batch_size'] = batch_size if batch_size!=1 else (nt if nt <= 20 else max(1, (nt + 9) // 5))
    # params.update({'batch_size', nt if nt <= 100 else max(1, (nt + 9) // 10)})

    if mode not in _MODE_TO_TARGET:
        raise ValueError(
            f"Invalid mode '{mode}'. Must be one of: 'serial', 'partial', or 'full'."
        )

    module_name, func_name = _MODE_TO_TARGET[mode]

    try:
        # Lazy import only the selected runner
        mod = importlib.import_module(module_name)
        runner = getattr(mod, func_name)
        return runner(**model_kwargs)
    except Exception:
        print(f"[ICESEE] Error in {mode} run mode:")
        tb_str = "".join(traceback.format_exception(*sys.exc_info()))
        print(f"Traceback details:\n{tb_str}")
