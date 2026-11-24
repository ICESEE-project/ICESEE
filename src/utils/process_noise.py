"""Temporal scaling shared by the canonical process-noise paths."""
from __future__ import annotations

import numpy as np


def process_noise_normalization(options):
    """Resolve a noise policy; omitted settings retain historical experiments."""
    policy = str(options.get("process_noise_normalization", "legacy_run_length")).strip().lower()
    if policy not in {"legacy_run_length", "per_step", "physical_ar1"}:
        raise ValueError("process_noise_normalization must be legacy_run_length, per_step or physical_ar1")
    return policy


def process_noise_coefficients(options):
    """Return (alpha, rho, tau) for zero-start AR(1) forecast increments.

    The increment is sqrt(dt)*sig_Q*rho*z, with
    z = alpha*z_previous + sqrt(1-alpha**2)*white.
    per_step uses rho=1: extending nt cannot change an earlier increment.
    sig_Q remains a coefficient, not the marginal increment standard deviation.
    """
    dt = float(options["dt"])
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("process-noise dt must be finite and positive")
    policy = process_noise_normalization(options)
    if policy == "physical_ar1":
        tau = float(options.get("process_noise_correlation_time", np.nan))
        if not np.isfinite(tau) or tau < 0:
            raise ValueError("physical_ar1 requires a nonnegative process_noise_correlation_time in model time units")
        if options.get("process_noise_start") not in {"stationary", "zero"}:
            raise ValueError("physical_ar1 requires process_noise_start=stationary or zero")
        return (np.exp(-dt / tau) if tau else 0.0), 1.0 / np.sqrt(dt), tau
    tau = max(500.0, dt)
    alpha = 1.0 - dt / tau
    if alpha <= 0 or alpha > 1:
        alpha = 0.5
    if policy == "per_step":
        return alpha, 1.0, tau
    nt = int(options["nt"])
    if nt <= 0:
        raise ValueError("legacy process-noise nt must be positive")
    denominator = nt - 2*alpha - nt*alpha**2 + 2*alpha**(nt+1)
    if not np.isfinite(denominator) or denominator <= 0:
        raise ValueError("legacy process-noise normalization has a nonpositive denominator")
    rho = np.sqrt((1/dt) * (1-alpha)**2 * (1/denominator))
    return alpha, float(rho), tau


def process_noise_event_coefficients(options, timestep):
    """Use elapsed physical time between actual AR updates."""
    if process_noise_normalization(options) != "physical_ar1":
        return (float(options.get("alpha", 0.0)),
                float(options.get("rho", 1.0)))
    _, rho, tau = process_noise_coefficients(options)
    k = int(timestep)
    start = int(options.get("process_noise_start_step", 0))
    schedule = str(options.get("process_noise_schedule", "observations")).lower()
    if schedule in {"every_step", "all", "forecast"}:
        previous = k - 1 if k > start else None
    elif schedule in {"observations", "observation", "analysis"}:
        events = np.asarray(options.get("obs_index", []), dtype=np.int64)
        events = events[:int(options.get("number_obs_instants", events.size))] - 1
        earlier = events[(events < k) & (events >= start)]
        previous = int(earlier.max()) if earlier.size else None
    else:
        raise ValueError("physical noise update requires an active noise schedule")
    if previous is None and options["process_noise_start"] == "stationary":
        return 0.0, rho
    elapsed = (k - previous if previous is not None else k - start + 1) * float(options["dt"])
    return (np.exp(-elapsed / tau) if tau else 0.0), rho


def process_noise_variable_mask(options):
    if "process_noise_vars" not in options:
        return [True] * int(options["num_state_vars"])
    names = list(options["vec_inputs"])[:int(options["num_state_vars"])]
    if len(names) != int(options["num_state_vars"]):
        raise ValueError("process_noise_vars requires names for every state block")
    selected = options.get("process_noise_vars", names)
    if not isinstance(selected, (list, tuple)) or set(selected) - set(names):
        raise ValueError("process_noise_vars must name configured state blocks")
    return [name in selected for name in names]


def process_noise_spatial_mode(options, variable):
    """Select a field draw or one scalar draw shared across a member's block."""
    policies = options.get("process_noise_spatial_policy", {})
    if not isinstance(policies, dict):
        raise ValueError("process_noise_spatial_policy must map state names to field or global_scalar")
    if not policies:
        return "field"
    names = list(options["vec_inputs"])[:int(options["num_state_vars"])]
    if set(policies) - set(names) or any(mode not in {"field", "global_scalar"} for mode in policies.values()):
        raise ValueError("invalid state name or process-noise spatial policy")
    return policies.get(names[int(variable)], "field")
