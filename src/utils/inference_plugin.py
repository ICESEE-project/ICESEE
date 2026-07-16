"""Feature-flagged integration hooks for SMB and bed inference.

Copy this module together with ``stable_smb_inference.py`` and
``stable_bed_inference.py`` into the package containing the analysis/writer
code. With ``inference_plugin_enabled=False`` (the default), the local bed gate
reproduces the legacy freeze behavior and the global hook is a no-op.
"""

import numpy as np

from .stable_bed_inference import apply_bed_regularized_correction
from .stable_smb_inference import apply_smb_physics_correction


_LEGACY_BED_ALIASES = {
    "bed",
    "bedrock",
    "bedtopography",
    "bed_topography",
}

_BED_ALIASES = {
    *_LEGACY_BED_ALIASES,
    "bed_elevation",
}


def _local_block_mask(global_rows, block_index, hdim):
    start = block_index * hdim
    end = start + hdim
    return (global_rows >= start) & (global_rows < end)


def apply_bed_update_gate_local(
    analysis_vec,
    forecast_vec,
    global_rows,
    vec_inputs,
    hdim,
    model_kwargs,
):
    """Gate distributed bed analysis updates without changing other fields.

    Modes
    -----
    legacy:
        Reproduce the original behavior exactly: when ``km`` is defined, bed
        is restored from the forecast unless ``km`` is in ``bed_snap_cols``.
    snapshots:
        Explicit snapshot-only inference. A nonempty ``bed_snap_cols`` is
        required, preventing accidental permanent freezing.
    continuous:
        Retain EnKF cross-covariance updates from H/u/v every analysis cycle.

    If the plugin itself is disabled, ``legacy`` is forced regardless of the
    configured new mode.
    """
    enabled = bool(model_kwargs.get("inference_plugin_enabled", False))
    mode = str(model_kwargs.get("bed_update_mode", "legacy")).lower()
    if not enabled:
        mode = "legacy"

    if mode not in {"legacy", "snapshots", "continuous"}:
        raise ValueError(
            "bed_update_mode must be 'legacy', 'snapshots', or 'continuous'"
        )

    km = model_kwargs.get("km")
    snapshot_columns = {
        int(value) for value in model_kwargs.get("bed_snap_cols", [])
    }
    is_snapshot = km is not None and int(km) in snapshot_columns

    if mode == "legacy":
        freeze_bed = km is not None and not is_snapshot
    elif mode == "snapshots":
        if not snapshot_columns:
            raise ValueError(
                "bed_update_mode='snapshots' requires nonempty bed_snap_cols"
            )
        freeze_bed = not is_snapshot
    else:
        freeze_bed = False

    model_kwargs["_bed_update_active"] = not freeze_bed
    model_kwargs["_bed_is_snapshot"] = is_snapshot

    if not freeze_bed:
        return analysis_vec

    aliases = _LEGACY_BED_ALIASES if mode == "legacy" else _BED_ALIASES
    for block_index, key in enumerate(vec_inputs):
        if str(key).lower() not in aliases:
            continue
        local_mask = _local_block_mask(global_rows, block_index, hdim)
        if np.any(local_mask):
            analysis_vec[local_mask, :] = forecast_vec[local_mask, :]
    return analysis_vec


def apply_global_inference_hook(
    analysis_vec,
    vec_inputs,
    hdim,
    params,
    model_kwargs,
    timestep,
    model_time=None,
    stage="all",
):
    """Apply enabled inference components to the gathered global ensemble.

    The hook is deliberately a no-op unless ``inference_plugin_enabled=True``.
    Use ``stage='pre_geometry'`` for bed before model-native geometry fixes and
    ``stage='post_geometry'`` for SMB after those fixes. ``stage='all'`` is
    retained for backward compatibility but should not be used for ISSM.
    """
    if not model_kwargs.get("inference_plugin_enabled", False):
        return analysis_vec

    stage = str(stage).lower()
    if stage not in {"pre_geometry", "post_geometry", "all"}:
        raise ValueError(
            "stage must be 'pre_geometry', 'post_geometry', or 'all'"
        )

    if model_time is None:
        dt = float(model_kwargs.get("dt", params.get("dt", 1.0)))
        model_time = float(timestep) * dt

    if stage in {"pre_geometry", "all"} and model_kwargs.get(
        "physics_bed_inference", False
    ):
        # In snapshot mode, do not repeatedly post-process a frozen forecast.
        bed_active = bool(model_kwargs.get("_bed_update_active", True))
        if bed_active:
            analysis_vec = apply_bed_regularized_correction(
                analysis_vec,
                vec_inputs,
                hdim,
                model_kwargs,
                timestep=timestep,
                model_time=model_time,
            )

    if stage in {"post_geometry", "all"} and model_kwargs.get(
        "physics_smb_inference", False
    ):
        analysis_vec = apply_smb_physics_correction(
            analysis_vec,
            vec_inputs,
            hdim,
            model_kwargs,
            timestep=timestep,
            model_time=model_time,
        )

    return analysis_vec


def reset_inference_plugin_state(model_kwargs):
    """Remove only private runtime state created by the inference plugin."""
    keys = (
        "_bed_update_active",
        "_bed_is_snapshot",
        "_bed_initial_reference",
        "_bed_previous_applied",
        "_bed_regularization_cache",
        "_bed_inference_call_count",
        "_smb_inference_history",
        "_smb_regularized_previous",
        "_smb_spinup_reference",
        "_smb_regularization_cache",
        "_smb_inference_call_count",
    )
    for key in keys:
        model_kwargs.pop(key, None)
