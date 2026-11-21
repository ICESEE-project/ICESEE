"""Scientific settings for prescribed or estimated PIG melt forcing."""
from __future__ import annotations

import hashlib
import json
import numpy as np

_MELT_PRIOR_NAMESPACE = 0xBA5A1


def melt_forecast_mode(options):
    mode = str(options.get("basal_melt_forecast", "prescribed_wrong")).strip().lower()
    if mode not in {"prescribed_wrong", "estimated_field"}:
        raise ValueError("basal_melt_forecast must be 'prescribed_wrong' or 'estimated_field'")
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
    if melt_forecast_mode(options) != "estimated_field":
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
    if mode == "prescribed_wrong" and error_mode == "legacy_prior_anomalies" and initial is None and not bounds:
        return None
    keys = ("basal_melt_forecast", "basal_melt_prior_scale_bounds", "basal_melt_bounds",
            "base_seed", "initial_perturbation_vars", "analysis_bounds", "initial_spread_factor",
            "joint_estimation", "wrong_basal_melt_field", "sig_obs", "sig_Q",
            "inflation_factor", "state_inflation_factor", "param_inflation_factor",
            "local_analysis", "global_analysis", "localization_flag", "localization_radius")
    payload = {key:options.get(key) for key in keys}
    payload.update(basal_melt_forecast=mode, enkf_observation_error_mode=error_mode)
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"),
                         default=lambda value:np.asarray(value).tolist())
    return "pig-forcing-v1-" + hashlib.sha256(encoded.encode()).hexdigest()[:20]


def configure_melt_experiment(options):
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
