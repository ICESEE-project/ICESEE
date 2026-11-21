"""Temporal scaling shared by the canonical process-noise paths."""
from __future__ import annotations

import numpy as np


def process_noise_normalization(options):
    """Resolve a noise policy; omitted settings retain historical experiments."""
    policy = str(options.get("process_noise_normalization", "legacy_run_length")).strip().lower()
    if policy not in {"legacy_run_length", "per_step"}:
        raise ValueError("process_noise_normalization must be 'legacy_run_length' or 'per_step'")
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
