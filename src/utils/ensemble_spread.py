"""Row-local relaxation to prior spread for ensemble analyses."""
from collections.abc import Mapping

import numpy as np


def spread_relaxation_settings(options, field_names):
    """Return validated RTPS coefficients; omitted settings leave the analysis unchanged."""
    names = tuple(map(str, field_names))
    default = options.get("analysis_spread_relaxation", 0.0)
    overrides = options.get("analysis_spread_relaxation_factors", {})
    if not isinstance(overrides, Mapping):
        raise ValueError("analysis_spread_relaxation_factors must map field names to coefficients")
    unknown = set(overrides) - set(names)
    if unknown:
        raise ValueError(f"unknown spread-relaxation fields: {sorted(unknown)}")

    def coefficient(value):
        if isinstance(value, (bool, np.bool_)):
            raise ValueError("spread-relaxation coefficients must be numbers in [0, 1]")
        try:
            value = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError("spread-relaxation coefficients must be numbers in [0, 1]") from exc
        if not np.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("spread-relaxation coefficients must be finite and in [0, 1]")
        return value

    default = coefficient(default)
    return {name: coefficient(overrides[name]) if name in overrides else default for name in names}


def relax_analysis_spread(analysis, forecast, global_rows, field_intervals, options):
    """Preserve each analysis row mean and relax its sample SD toward the prior.

    ``sigma_new = (1-alpha)*sigma_analysis + alpha*sigma_forecast``.
    The caller supplies a bounded row block containing every logical member.
    No random draws or communication occur here. Exactly collapsed analysis
    rows use the centered prior direction; a prior with no spread stays collapsed.
    Bounds and frozen-field policies must be applied afterwards by the caller.
    """
    settings = spread_relaxation_settings(options, field_intervals)
    if not any(settings.values()):
        return analysis
    analysis = np.asarray(analysis, dtype=np.float64)
    forecast = np.asarray(forecast, dtype=np.float64)
    rows = np.asarray(global_rows)
    if analysis.ndim != 2 or analysis.shape != forecast.shape or analysis.shape[1] < 2:
        raise ValueError("RTPS requires aligned row-by-member arrays with at least two members")
    if rows.shape != (analysis.shape[0],) or rows.dtype.kind not in "iu":
        raise ValueError("RTPS global rows must be an aligned integer vector")
    if not np.isfinite(analysis).all() or not np.isfinite(forecast).all():
        raise FloatingPointError("RTPS requires finite analysis and forecast members")
    coefficients = np.zeros(len(rows))
    coverage = np.zeros(len(rows), dtype=int)
    for name, (begin, end) in field_intervals.items():
        mask = (rows >= int(begin)) & (rows < int(end))
        coefficients[mask] = settings[name]
        coverage[mask] += 1
    if np.any(coverage != 1):
        raise ValueError("RTPS field intervals must cover every selected row exactly once")
    active = coefficients > 0
    result = analysis.copy()
    mean = analysis[active].mean(axis=1, keepdims=True)
    anomalies = analysis[active] - mean
    prior = forecast[active] - forecast[active].mean(axis=1, keepdims=True)
    # Recenter before normalization to limit roundoff for large physical means.
    anomalies -= anomalies.mean(axis=1, keepdims=True)
    prior -= prior.mean(axis=1, keepdims=True)
    sd_a = np.linalg.norm(anomalies, axis=1) / np.sqrt(analysis.shape[1] - 1)
    sd_f = np.linalg.norm(prior, axis=1) / np.sqrt(analysis.shape[1] - 1)
    alpha = coefficients[active]
    target = (1-alpha)*sd_a + alpha*sd_f
    direction = np.zeros_like(anomalies)
    present = sd_a > 0
    direction[present] = anomalies[present] / sd_a[present, None]
    fallback = ~present & (sd_f > 0)
    direction[fallback] = prior[fallback] / sd_f[fallback, None]
    updated = direction * target[:, None]
    updated -= updated.mean(axis=1, keepdims=True)
    result[active] = mean + updated
    if not np.isfinite(result).all():
        raise FloatingPointError("RTPS produced non-finite members")
    return result
