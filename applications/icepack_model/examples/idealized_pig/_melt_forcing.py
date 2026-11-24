"""Scientific settings for prescribed or estimated PIG melt forcing."""
from __future__ import annotations

import hashlib
import json
import numpy as np

_MELT_PRIOR_NAMESPACE = 0xBA5A1


def melt_spatial_policy(options):
    policy = str(options.get("basal_melt_spatial_policy", "evolving_depth")).lower()
    if policy not in {"evolving_depth", "fixed_initial_depth"}:
        raise ValueError("basal_melt_spatial_policy must be evolving_depth or fixed_initial_depth")
    return policy


def bind_melt_reference(options, context):
    """Use the model group's initial geometry, never another MPI mesh's fields."""
    if melt_spatial_policy(options) == "fixed_initial_depth":
        return dict(options, h0=context.h0, s0=context.s0)
    return options


def melt_forecast_mode(options):
    mode = str(options.get("basal_melt_forecast", "prescribed_wrong")).strip().lower()
    if mode not in {"prescribed_wrong", "estimated_field", "estimated_depth_law"}:
        raise ValueError("basal_melt_forecast must be prescribed_wrong, estimated_field or estimated_depth_law")
    return mode


def _interval(options, key):
    value = options.get(key)
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(f"{key} must explicitly specify [lower, upper]")
    lower, upper = map(float, value)
    if not np.isfinite([lower, upper]).all() or lower < 0 or upper <= lower:
        raise ValueError(f"{key} requires finite 0 <= lower < upper")
    return lower, upper


def estimated_melt_settings(options):
    if melt_forecast_mode(options) == "prescribed_wrong":
        return None
    if "basal_melt_field" not in options.get("vec_inputs", ()):
        raise ValueError("estimated_field requires basal_melt_field in the assimilated state")
    return (_interval(options, "basal_melt_prior_scale_bounds"),
            _interval(options, "basal_melt_bounds"))


def initial_melt_scale(options, member_id):
    settings = estimated_melt_settings(options)
    if settings is None:
        return 1.0
    if int(member_id) < 0:
        raise ValueError("melt-prior member_id must be nonnegative")
    prior, _ = settings
    seed = np.random.SeedSequence([int(options.get("base_seed", 42)),
                                   int(member_id), _MELT_PRIOR_NAMESPACE])
    return float(np.random.default_rng(seed).uniform(*prior))


def bounded_melt_values(values, options):
    values = np.asarray(values, dtype=np.float64)
    if not np.isfinite(values).all():
        raise FloatingPointError("melt state contains non-finite values")
    settings = estimated_melt_settings(options)
    return values.copy() if settings is None else np.clip(values, *settings[1])


def initial_melt_values(nominal_values, options, member_id):
    return bounded_melt_values(
        np.asarray(nominal_values) * initial_melt_scale(options, member_id), options)


def melt_experiment_identity(options):
    """Separate changed scientific policies from compatible old checkpoints."""
    mode = melt_forecast_mode(options)
    estimated_melt_settings(options)
    error_mode = str(options.get("enkf_observation_error_mode", "legacy_prior_anomalies")).lower()
    initial = options.get("initial_perturbation_vars")
    bounds = options.get("analysis_bounds") or {}
    from ICESEE.src.utils.ensemble_spread import spread_relaxation_settings
    spread = spread_relaxation_settings(options, options.get("vec_inputs", ()))
    spatial_policy = melt_spatial_policy(options)
    if mode == "prescribed_wrong" and error_mode == "legacy_prior_anomalies" and initial is None and not bounds and analysis_state_policy(options) == "independent_fields" and options.get("analysis_enabled", True) and not any(spread.values()) and spatial_policy == "evolving_depth":
        return None
    keys = ("basal_melt_forecast", "basal_melt_prior_scale_bounds", "basal_melt_bounds",
            "base_seed", "initial_perturbation_vars", "analysis_bounds", "initial_spread_factor",
            "joint_estimation", "wrong_basal_melt_field", "sig_obs", "sig_Q",
            "inflation_factor", "state_inflation_factor", "param_inflation_factor",
            "local_analysis", "global_analysis", "localization_flag", "localization_radius")
    payload = {key:options.get(key) for key in keys}
    if spatial_policy != "evolving_depth":
        payload["basal_melt_spatial_policy"] = "fixed-initial-control-depth-v1"
    if any(spread.values()):
        payload["analysis_spread_policy"] = {"method": "rtps-v1-prior-direction-at-zero", "factors": spread}
    for key in ("enkf_gain_formulation", "basal_melt_amplitude_trend", "analysis_enabled", "basal_melt_noise_basis"):
        if key in options:
            payload[key] = options[key]
    payload.update(basal_melt_forecast=mode, enkf_observation_error_mode=error_mode)
    if analysis_state_policy(options) != "independent_fields":
        payload.update(pig_analysis_state_policy=analysis_state_policy(options),
                       hThresh=options.get("hThresh"),
                       frozen_analysis_vars=options.get("frozen_analysis_vars"))
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"),
                         default=lambda value:np.asarray(value).tolist())
    return "pig-forcing-v1-" + hashlib.sha256(encoded.encode()).hexdigest()[:20]


def analysis_state_policy(options):
    policy = str(options.get("pig_analysis_state_policy", "independent_fields")).lower()
    if policy not in {"independent_fields", "diagnostic_balance"}:
        raise ValueError("pig_analysis_state_policy must be independent_fields or diagnostic_balance")
    return policy


def configure_melt_experiment(options):
    from ICESEE.src.parallelization.forecast_recovery import ForecastRecoveryPolicy
    recovery = ForecastRecoveryPolicy.from_options(options)
    if recovery is not None and int(options.get("execution_mode", 3)) != 3:
        raise NotImplementedError("PIG forecast_recovery requires the native Mode-3 adapter")
    from ICESEE.src.parallelization.distributed_analysis import resolve_gain_formulation
    resolve_gain_formulation(options)
    melt_spatial_policy(options)
    from ICESEE.src.utils.ensemble_spread import spread_relaxation_settings
    spread_relaxation_settings(options, options.get("vec_inputs", ()))
    if options.get("analysis_stage_diagnostics", False):
        if int(options.get("execution_mode", 3)) != 3 or options.get("use_store_streaming_analysis", False):
            raise NotImplementedError("analysis_stage_diagnostics requires the persistent native Mode-3 path")
    noise_basis = options.get("basal_melt_noise_basis", "field")
    if noise_basis not in {"field", "depth_amplitude"}:
        raise ValueError("basal_melt_noise_basis must be field or depth_amplitude")
    if noise_basis == "depth_amplitude" and (
            melt_forecast_mode(options) != "estimated_depth_law"
            or options.get("process_noise_normalization") != "physical_ar1"
            or options.get("process_noise_spatial_policy", {}).get("basal_melt_field") != "global_scalar"):
        raise ValueError("depth_amplitude noise requires estimated_depth_law, physical_ar1 and global_scalar melt noise")
    limit = options.get("pig_max_thickness_projection")
    if limit is not None and (not np.isfinite(float(limit)) or float(limit) < 0):
        raise ValueError("pig_max_thickness_projection must be finite and nonnegative")
    if melt_forecast_mode(options) == "estimated_depth_law":
        if int(options.get("execution_mode", 3)) != 3:
            raise ValueError("estimated_depth_law currently requires execution_mode=3")
        raw_trend = options.get("basal_melt_amplitude_trend")
        trend = float(raw_trend) if raw_trend is not None else np.nan
        if not np.isfinite(trend):
            raise ValueError("estimated_depth_law requires an explicit basal_melt_amplitude_trend (m/year per year)")
        if analysis_state_policy(options) != "diagnostic_balance":
            raise ValueError("estimated_depth_law requires diagnostic_balance")
    if analysis_state_policy(options) == "diagnostic_balance":
        if int(options.get("execution_mode", 3)) != 3:
            raise ValueError("diagnostic_balance currently requires execution_mode=3")
        required = {"h", "u", "v", "s"}
        if not required.issubset(options.get("vec_inputs", ())):
            raise ValueError("diagnostic_balance requires h, u, v and s fields")
        if not np.isfinite(float(options.get("hThresh", np.nan))) or float(options["hThresh"]) <= 0:
            raise ValueError("diagnostic_balance requires the existing positive hThresh")
        options["frozen_analysis_vars"] = sorted(set(options.get("frozen_analysis_vars") or ()) | {"u", "v", "s"})
    settings = estimated_melt_settings(options)
    if settings is not None:
        options.setdefault("initial_perturbation_vars", [])
        bounds = dict(options.get("analysis_bounds") or {})
        matches = [key for key in bounds if str(key).lower()=="basal_melt_field"]
        if matches and tuple(map(float,bounds[matches[0]])) != settings[1]:
            raise ValueError("analysis_bounds for basal_melt_field disagree with basal_melt_bounds")
        bounds[matches[0] if matches else "basal_melt_field"] = list(settings[1])
        options["analysis_bounds"] = bounds
    options["scientific_experiment_id"] = melt_experiment_identity(options)
    return options


def melt_scenario(model_time):
    """Known control/warm chronology shared by the depth-law predictor."""
    return "warm" if any(a <= model_time < b for a, b in (
        (6, 15), (18, 20), (25, 27), (31, 40), (48, 50), (59, 64), (69, 76))) else "control"


def advance_melt_amplitude(numerator, denominator, dt, options):
    """Area-weighted projection followed by the explicitly configured drift.

    The numerator/denominator are global finite-element integrals. No truth
    values enter this predictor. A zero-support field cannot identify amplitude.
    """
    if not np.isfinite([numerator, denominator, dt]).all() or denominator <= 0 or dt < 0:
        raise ValueError("estimated depth-law amplitude is unidentifiable or has invalid time")
    amplitude = numerator / denominator + float(options["basal_melt_amplitude_trend"]) * dt
    return float(np.clip(amplitude, *estimated_melt_settings(options)[1]))


def resolve_stop_step(options):
    """Limit executed cycles without changing the configured forcing horizon."""
    nt, dt = int(options["nt"]), float(options["dt"])
    stop = options.get("run_stop_time")
    if stop is None:
        return nt
    stop = float(stop)
    index = stop / dt
    if (not np.isfinite(index) or stop <= 0 or index > nt
            or not np.isclose(index, round(index), rtol=0, atol=1e-8)):
        raise ValueError("run_stop_time must lie on the existing time grid within the forcing horizon")
    if int(options.get("checkpoint_every", 1)) != 1:
        raise ValueError("run_stop_time requires checkpoint_every=1 to preserve the stopping state")
    return int(round(index))


def validate_reference_config(options):
    """Check the declared frozen configuration before reusing reference files."""
    from pathlib import Path
    import yaml
    root = options.get("reference_data_path")
    if root is None:
        return
    config_path = options.get("reference_params_file")
    if config_path is None:
        raise ValueError("reference_data_path requires reference_params_file from the generating run")
    config = yaml.safe_load(Path(config_path).read_text())
    original = {}
    for section in config.values():
        if isinstance(section, dict):
            original.update(section)
    if melt_spatial_policy(original) != melt_spatial_policy(options):
        raise ValueError("reference configuration differs in basal_melt_spatial_policy; generate separate truth/observations")
    for key in ("num_years", "timesteps_per_year", "sig_obs", "obs_start_time", "obs_max_time", "freq_obs"):
        if key not in original:
            raise ValueError(f"reference_params_file is missing {key}")
    keys = ("num_years", "timesteps_per_year", "initFile", "paramsFile", "meshFile", "SMBFile",
            "water_to_ice", "water_density", "GLThresh", "uThresh", "hThresh", "bmr_increase_time",
            "sig_obs", "observed_vars", "observed_params", "obs_start_time", "obs_max_time", "freq_obs",
            "vec_inputs", "seed", "base_seed")
    for key in keys:
        if key in original and not np.array_equal(np.asarray(original[key]), np.asarray(options.get(key))):
            raise ValueError(f"reference configuration differs in {key}; regenerate references for that experiment")
    root = Path(root).resolve()
    from ICESEE.src.utils.reference_data import validate_reference_paths
    validate_reference_paths(options["data_path"], root)
    for name in ("true_nurged_states.h5", "synthetic_obs.h5"):
        if not (root / name).is_file():
            raise ValueError(f"missing cached reference file: {root / name}")


def diagnostic_fork_info(options):
    """Declare a new experiment from a closed pre-first-analysis checkpoint.

    This is not a compatible restart. The parent state and manifest identity
    are retained, but changed scientific policies begin at the fork boundary.
    """
    from pathlib import Path
    source = options.get("diagnostic_fork_checkpoint")
    if source is None:
        return None
    if options.get("mode3_export_ensemble_history", False):
        raise ValueError("diagnostic forks write checkpoint windows, not a complete history from time zero")
    validate_reference_config(options)
    import yaml
    parent_config = yaml.safe_load(Path(options["reference_params_file"]).read_text())
    parent_schedule = parent_config.get("enkf-parameters", {}).get("process_noise_schedule", "observations")
    if parent_schedule not in {"observations", "observation", "analysis"}:
        raise ValueError("diagnostic fork requires an observation-only parent noise schedule with no previous updates")
    path = Path(source).resolve()
    from ICESEE.src.utils.reference_data import validate_reference_paths
    validate_reference_paths(options["data_path"], path)
    raw = (path / "manifest.json").read_bytes()
    manifest = json.loads(raw)
    meta = manifest.get("metadata", {})
    next_step = int(meta.get("next_step", -1))
    events = np.asarray(options.get("obs_index", []), dtype=int)
    events = events[:int(options.get("number_obs_instants", events.size))]
    valid = (manifest.get("complete") is True
        and int(manifest.get("number_of_members", -1)) == int(options["Nens"])
        and int(meta.get("nt", -1)) == int(options["nt"])
        and int(meta.get("analysis_events_completed", -1)) == 0
        and next_step == int(manifest.get("timestep", -2)) + 1
        and next_step > 0 and events.size > 0 and next_step < int(events.min())
        and next_step < resolve_stop_step(options)
        and np.isclose(float(meta.get("model_time", np.nan)), next_step * float(options["dt"]), rtol=0, atol=1e-10))
    if not valid:
        raise ValueError("diagnostic fork requires a complete matching pre-first-analysis checkpoint and a later stopping time")
    return dict(path=str(path), parent_run_id=manifest.get("run_id"),
        manifest_sha256=hashlib.sha256(raw).hexdigest(), completed_step=next_step-1,
        next_step=next_step, model_time=float(meta["model_time"]))
