"""Canonical runtime-context helpers for ICESEE.

ICESEE passes one flat mutable dictionary, ``icesee_kwargs``, through its
configuration, model, analysis, and I/O layers.  The normalizer retains a
narrow compatibility boundary for callers using the former mappings.
"""

from collections.abc import Mapping
from pathlib import Path
import warnings

import numpy as np


_LEGACY_CONTEXT_KEYS = ("params", "model_kwargs", "kwargs")


def normalize_icesee_kwargs(context=None, **values):
    """Return one flat ICESEE runtime context.

    Direct values take precedence over legacy nested mappings. The original
    dictionary is returned when no merge is needed so shared runtime state is
    preserved.
    """
    if context is None:
        context = {}
    elif not isinstance(context, Mapping):
        raise TypeError("icesee_kwargs must be a mapping")
    elif not isinstance(context, dict):
        context = dict(context)

    legacy = []
    for key in _LEGACY_CONTEXT_KEYS:
        nested = values.pop(key, None)
        if isinstance(nested, Mapping):
            legacy.append((key, nested))
    if isinstance(context, dict):
        for key in _LEGACY_CONTEXT_KEYS:
            nested = context.get(key)
            if isinstance(nested, Mapping):
                legacy.append((key, nested))

    if not legacy and not values:
        return context

    merged = {}
    for key, nested in legacy:
        warnings.warn(
            f"'{key}' is deprecated; pass one flat icesee_kwargs dictionary",
            DeprecationWarning,
            stacklevel=2,
        )
        merged.update(nested)
    merged.update({k: v for k, v in context.items() if k not in _LEGACY_CONTEXT_KEYS})
    merged.update(values)
    return merged


_UNSUPPORTED = object()


def _matlab_value(value):
    """Convert a runtime-context value to a scipy.io.savemat-safe value."""
    if value is None or callable(value):
        return _UNSUPPORTED
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, bytes, bool, int, float, complex, np.generic)):
        return value
    if isinstance(value, np.ndarray):
        if value.dtype != object:
            return value
        converted = [_matlab_value(item) for item in value.flat]
        if any(item is _UNSUPPORTED for item in converted):
            return _UNSUPPORTED
        return np.asarray(converted, dtype=object).reshape(value.shape)
    if isinstance(value, Mapping):
        return {
            key: item
            for key, nested in value.items()
            if (item := _matlab_value(nested)) is not _UNSUPPORTED
        }
    if isinstance(value, (list, tuple)):
        converted = [_matlab_value(item) for item in value]
        if any(item is _UNSUPPORTED for item in converted):
            return _UNSUPPORTED
        return converted
    return _UNSUPPORTED


def matlab_icesee_kwargs(icesee_kwargs):
    """Return the serializable view used at the ICESEE--MATLAB boundary.

    The live runtime context remains the single source of truth. Runtime-only
    objects such as MPI communicators, server handles, modules, and callables
    are excluded from the MAT-file handoff rather than being copied into a
    second configuration dictionary.
    """
    context = normalize_icesee_kwargs(icesee_kwargs)
    return {
        key: value
        for key, raw_value in context.items()
        if (value := _matlab_value(raw_value)) is not _UNSUPPORTED
    }
