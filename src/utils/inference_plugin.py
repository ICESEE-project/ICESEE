"""Feature-flagged integration hooks for SMB and bed inference.

Copy this module together with ``stable_smb_inference.py`` and
``stable_bed_inference.py`` into the package containing the analysis/writer
code. With ``inference_plugin_enabled=False`` (the default), the local bed gate
reproduces the legacy freeze behavior and the global hook is a no-op.
"""

import os

import h5py
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


def resolve_analysis_cycle_time(icesee_kwargs, k, km):
    """Return the scheduled assimilation time and its mapped model time.

    Observation schedules retain both the requested times (``obs_t``) and the
    nearest model indices (``obs_index``).  Controls such as
    ``inversion_start_time`` describe the observation schedule, so they must
    be evaluated against ``obs_t[km]`` rather than the slightly offset model
    grid time ``t[k]``.  The latter is returned as well for transparent logs.
    """
    model_time = float(np.asarray(icesee_kwargs["t"], dtype=float)[int(k)])
    observation_times = np.asarray(
        icesee_kwargs.get("obs_t", []), dtype=float
    ).reshape(-1)
    if 0 <= int(km) < observation_times.size:
        return float(observation_times[int(km)]), model_time
    return model_time, model_time


def _local_block_mask(global_rows, block_index, hdim):
    start = block_index * hdim
    end = start + hdim
    return (global_rows >= start) & (global_rows < end)


def resolve_bed_update_schedule(icesee_kwargs):
    """Return ``(active, is_snapshot, mode)`` for the current analysis cycle.

    This is intentionally a pure calculation.  Earlier code communicated the
    result from row-local analysis to global finalization through private keys
    in ``icesee_kwargs``.  That is unsafe because several execution-mode-2
    boundaries copy the dictionary; global finalization could consequently
    see a stale or missing ``_bed_update_active`` value and skip the direct bed
    anchor.  Both execution modes now resolve the same schedule independently
    from the public configuration and current observation column.
    """
    enabled = bool(icesee_kwargs.get("inference_plugin_enabled", False))
    mode = str(icesee_kwargs.get("bed_update_mode", "legacy")).lower()
    if not enabled:
        mode = "legacy"

    if mode not in {"legacy", "snapshots", "continuous"}:
        raise ValueError(
            "bed_update_mode must be 'legacy', 'snapshots', or 'continuous'"
        )

    km = icesee_kwargs.get("km")
    snapshot_columns = {
        int(value) for value in icesee_kwargs.get("bed_snap_cols", [])
    }
    is_snapshot = km is not None and int(km) in snapshot_columns

    if mode == "legacy":
        active = not (km is not None and not is_snapshot)
    elif mode == "snapshots":
        if not snapshot_columns:
            raise ValueError(
                "bed_update_mode='snapshots' requires nonempty bed_snap_cols"
            )
        active = is_snapshot
    else:
        active = True

    return bool(active), bool(is_snapshot), mode


def apply_bed_update_gate_local(
    analysis_vec,
    forecast_vec,
    global_rows,
    vec_inputs,
    hdim,
    icesee_kwargs,
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
    active, is_snapshot, mode = resolve_bed_update_schedule(icesee_kwargs)
    freeze_bed = not active

    icesee_kwargs["_bed_update_active"] = active
    icesee_kwargs["_bed_is_snapshot"] = is_snapshot

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


def apply_bed_domain_gate_global(
    analysis_vec,
    forecast_vec,
    vec_inputs,
    hdim,
    icesee_kwargs,
):
    """Restrict bed increments to the configured physical/support domain.

    ``grounded_only`` diagnoses grounding independently for each forecast
    member. ``observed_only`` is stricter: it uses the truth-supported bed mask
    stored with the synthetic observations and therefore cannot leak an update
    beneath unobserved floating ice through a misclassified ensemble member.
    """
    domain = str(icesee_kwargs.get("bed_update_domain", "all")).lower()
    if domain == "all":
        return analysis_vec
    if domain not in {"grounded_only", "observed_only"}:
        raise ValueError(
            "bed_update_domain must be 'all', 'grounded_only', or "
            "'observed_only'"
        )

    km = icesee_kwargs.get("km")
    snapshot_columns = {
        int(value) for value in icesee_kwargs.get("bed_snap_cols", [])
    }
    if km is None or int(km) not in snapshot_columns:
        return analysis_vec

    thickness_block = None
    bed_block = None
    for block, key in enumerate(vec_inputs):
        key_l = str(key).lower()
        if key_l in {"thickness", "ice_thickness", "h"}:
            thickness_block = block
        elif key_l in _BED_ALIASES:
            bed_block = block
    if bed_block is None or (domain == "grounded_only" and thickness_block is None):
        raise ValueError(
            f"{domain} bed updates require the necessary state blocks"
        )

    thickness_slice = None
    if thickness_block is not None:
        thickness_slice = slice(
            thickness_block * hdim, (thickness_block + 1) * hdim
        )
    bed_slice = slice(bed_block * hdim, (bed_block + 1) * hdim)
    forecast_bed = np.asarray(forecast_vec[bed_slice, :], dtype=float)

    bed_analysis = np.asarray(analysis_vec[bed_slice, :], dtype=float).copy()
    if domain == "grounded_only":
        forecast_thickness = np.asarray(
            forecast_vec[thickness_slice, :], dtype=float
        )
        density_ratio = float(icesee_kwargs.get("di", 0.8930))
        allowed = forecast_thickness + forecast_bed / density_ratio > 0.0
    else:
        masks_by_key = icesee_kwargs.get("bed_mask_map_cols", {})
        support = None
        for key in vec_inputs:
            if str(key).lower() not in _BED_ALIASES:
                continue
            mask_columns = masks_by_key.get(key)
            if mask_columns is None:
                # HDF5 keys and configured aliases may differ only by case.
                for stored_key, stored_value in masks_by_key.items():
                    if str(stored_key).lower() == str(key).lower():
                        mask_columns = stored_value
                        break
            if mask_columns is not None:
                mask_columns = np.asarray(mask_columns, dtype=bool)
                if (
                    mask_columns.ndim == 2
                    and mask_columns.shape[1] == hdim
                    and mask_columns.shape[0] != hdim
                ):
                    mask_columns = mask_columns.T
                if mask_columns.shape[0] != hdim or int(km) >= mask_columns.shape[1]:
                    raise ValueError(
                        "bed observation-support mask has an incompatible "
                        f"shape: mask={mask_columns.shape}, hdim={hdim}, "
                        f"obs_column={km}"
                    )
                support = mask_columns[:, int(km)]
                break
        if support is None:
            raise ValueError(
                "bed_update_domain='observed_only' requires bed_mask_map_cols"
            )
        allowed = support[:, None]

    bed_analysis[~np.broadcast_to(allowed, bed_analysis.shape)] = (
        forecast_bed[~np.broadcast_to(allowed, forecast_bed.shape)]
    )
    analysis_vec[bed_slice, :] = bed_analysis
    return analysis_vec


def apply_bed_observation_anchor_global(
    analysis_vec,
    vec_inputs,
    hdim,
    icesee_kwargs,
    stage="pre",
):
    """Anchor surveyed bed nodes directly to their active observations.

    The augmented EnKF still supplies spatial cross-covariances, while this
    scalar relaxation guarantees that a direct bed observation cannot increase
    its own innovation.  Unobserved nodes are untouched here and are handled by
    the regularized/localized increment plus the physical domain gate.
    """
    active, is_snapshot, _ = resolve_bed_update_schedule(icesee_kwargs)
    # Keep these keys only as diagnostics/backward-compatible plugin state;
    # correctness no longer depends on a mutation surviving a kwargs copy.
    icesee_kwargs["_bed_update_active"] = active
    icesee_kwargs["_bed_is_snapshot"] = is_snapshot
    if not active:
        return analysis_vec

    stage = str(stage).lower()
    if stage not in {"pre", "post"}:
        raise ValueError("bed observation anchor stage must be 'pre' or 'post'")
    factor_key = (
        "bed_observation_nudge_factor"
        if stage == "pre"
        else "bed_observation_post_anchor_factor"
    )
    factor = float(icesee_kwargs.get(factor_key, 0.0))
    if factor == 0.0:
        return analysis_vec
    if not 0.0 <= factor <= 1.0:
        raise ValueError(f"{factor_key} must be in [0, 1]")

    km = icesee_kwargs.get("km")
    if km is None:
        return analysis_vec

    bed_block = None
    for block_index, key in enumerate(vec_inputs):
        if str(key).lower() in _BED_ALIASES:
            bed_block = block_index
            break
    if bed_block is None:
        return analysis_vec

    bed_slice = slice(bed_block * hdim, (bed_block + 1) * hdim)
    observations = icesee_kwargs.get("hu_obs_loaded")
    if observations is not None:
        observations = np.asarray(observations, dtype=float)
        if observations.ndim != 2 or int(km) >= observations.shape[1]:
            raise ValueError("bed observation array has an incompatible shape")
        if bed_slice.stop > observations.shape[0]:
            raise ValueError("bed observation block is outside hu_obs_loaded")
        target = observations[bed_slice, int(km)]
    else:
        # Fully parallel mode stores observations in compact form to avoid an
        # O(full-state-size * observation-times) dense array.  Resolve only
        # the current bed column here so the shared finalization contract is
        # numerically identical to partial-parallel mode without sacrificing
        # the compact-storage design.
        obs_path = icesee_kwargs.get("synthetic_obs_file")
        if not obs_path:
            obs_path = os.path.join(
                str(icesee_kwargs.get("data_path", ".")), "synthetic_obs.h5"
            )
        if not os.path.exists(obs_path):
            raise FileNotFoundError(
                "Bed observation anchoring requires synthetic_obs.h5; "
                f"not found at {obs_path}"
            )
        with h5py.File(obs_path, "r") as handle:
            if not {"hu_obs_compact", "obs_indices"}.issubset(handle.keys()):
                raise ValueError(
                    "synthetic_obs.h5 contains neither hu_obs_loaded nor the "
                    "required compact observation datasets"
                )
            values_dset = handle["hu_obs_compact"]
            if values_dset.ndim != 2 or int(km) >= values_dset.shape[1]:
                raise ValueError(
                    "compact bed observation array has an incompatible shape"
                )
            global_rows = np.asarray(handle["obs_indices"], dtype=np.int64).reshape(-1)
            if global_rows.size != values_dset.shape[0]:
                raise ValueError(
                    "obs_indices and hu_obs_compact have incompatible shapes"
                )
            in_bed = (global_rows >= bed_slice.start) & (global_rows < bed_slice.stop)
            positions = np.flatnonzero(in_bed)
            target = np.full(hdim, np.nan, dtype=float)
            if positions.size:
                values = np.asarray(values_dset[positions, int(km)], dtype=float)
                if "obs_active" in handle:
                    active_dset = handle["obs_active"]
                    if active_dset.shape != values_dset.shape:
                        raise ValueError(
                            "obs_active and hu_obs_compact have incompatible shapes"
                        )
                    compact_active = np.asarray(
                        active_dset[positions, int(km)], dtype=bool
                    )
                else:
                    compact_active = np.isfinite(values)
                local_rows = global_rows[positions] - bed_slice.start
                target[local_rows[compact_active]] = values[compact_active]

    active = np.isfinite(target)
    if not np.any(active):
        return analysis_vec

    bed = np.asarray(analysis_vec[bed_slice, :], dtype=float).copy()
    mean_before = np.mean(bed[active, :], axis=1)
    bed[active, :] += factor * (target[active, None] - bed[active, :])
    if bool(icesee_kwargs.get("bed_observation_diagnostics", False)):
        mean_after = np.mean(bed[active, :], axis=1)
        rmse_before = np.sqrt(np.mean((mean_before - target[active]) ** 2))
        rmse_after = np.sqrt(np.mean((mean_after - target[active]) ** 2))
        print(
            "[ICESEE][bed_anchor] "
            f"stage={stage} obs_column={int(km)} n={int(np.sum(active))} "
            f"innovation_RMSE={rmse_before:.3f}->{rmse_after:.3f} m"
        )
    analysis_vec[bed_slice, :] = bed
    return analysis_vec


def apply_global_inference_hook(
    analysis_vec,
    vec_inputs,
    hdim,
    icesee_kwargs,
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
    if not icesee_kwargs.get("inference_plugin_enabled", False):
        return analysis_vec

    stage = str(stage).lower()
    if stage not in {"pre_geometry", "post_geometry", "all"}:
        raise ValueError(
            "stage must be 'pre_geometry', 'post_geometry', or 'all'"
        )

    if model_time is None:
        dt = float(icesee_kwargs.get("dt", 1.0))
        model_time = float(timestep) * dt

    if stage in {"pre_geometry", "all"} and icesee_kwargs.get(
        "physics_bed_inference", False
    ):
        # In snapshot mode, do not repeatedly post-process a frozen forecast.
        bed_active = bool(icesee_kwargs.get("_bed_update_active", True))
        if bed_active:
            analysis_vec = apply_bed_regularized_correction(
                analysis_vec,
                vec_inputs,
                hdim,
                icesee_kwargs,
                timestep=timestep,
                model_time=model_time,
            )

    if stage in {"post_geometry", "all"} and icesee_kwargs.get(
        "physics_smb_inference", False
    ):
        analysis_vec = apply_smb_physics_correction(
            analysis_vec,
            vec_inputs,
            hdim,
            icesee_kwargs,
            timestep=timestep,
            model_time=model_time,
        )

    return analysis_vec


def reset_inference_plugin_state(icesee_kwargs):
    """Remove only private runtime state created by the inference plugin."""
    keys = (
        "_bed_update_active",
        "_bed_is_snapshot",
        "_bed_initial_reference",
        "_bed_previous_applied",
        "_bed_forecast_reference",
        "_bed_regularization_cache",
        "_bed_inference_call_count",
        "_smb_inference_history",
        "_smb_regularized_previous",
        "_smb_spinup_reference",
        "_smb_regularization_cache",
        "_smb_inference_call_count",
    )
    for key in keys:
        icesee_kwargs.pop(key, None)
