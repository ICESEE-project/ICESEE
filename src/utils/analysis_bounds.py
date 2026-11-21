"""Optional row-local variable bounds for an analyzed ensemble."""
from collections.abc import Mapping
import numpy as np


def bound_analysis_rows(values, global_rows, variable_intervals, bounds=None):
    if not bounds:
        return values
    if not isinstance(bounds, Mapping):
        raise ValueError("analysis_bounds must map variable names to [lower, upper]")
    intervals = {str(name).lower():tuple(interval) for name,interval in variable_intervals.items()}
    specifications = []
    for name,limits in bounds.items():
        key = str(name).lower()
        if key not in intervals or not isinstance(limits,(list,tuple)) or len(limits)!=2:
            raise ValueError(f"invalid analysis_bounds for {name!r}")
        lower,upper = map(float,limits)
        if not np.isfinite([lower,upper]).all() or lower>upper:
            raise ValueError(f"analysis_bounds for {name!r} require finite lower <= upper")
        specifications.append((intervals[key],lower,upper))
    rows = np.asarray(global_rows)
    array = np.asarray(values)
    if rows.ndim!=1 or array.ndim!=2 or array.shape[0]!=rows.size:
        raise ValueError("analysis bounds require one global row for each ensemble row")
    result = array.copy()
    for (start,stop),lower,upper in specifications:
        mask = (rows>=int(start)) & (rows<int(stop))
        result[mask] = np.clip(result[mask],lower,upper)
    return result
