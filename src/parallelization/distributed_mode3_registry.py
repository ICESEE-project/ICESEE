"""Per-application opt-in registry for ``execution_mode: 3``.

``execution_mode`` 0, 1, and 2 are available to every ICESEE application
unconditionally. Mode 3 is different: per
``docs/execution-mode-3-design.md``, it requires a model adapter that
exposes genuinely distributed state (no rank holding a complete member) plus
a production DA-cycle runner built on
``src/parallelization/distributed_native_cycle.py``. Most applications do not
have one yet, and a reference-only adapter (for example
``applications/lorenz_model/lorenz_utils/distributed_adapter.py``, whose own
docstring states it is "used only for state-only parity development") must
not be mistaken for production support.

This module makes mode-3 access application-specific rather than globally
on/off: an application registers itself here only once it has a real
production adapter and runner, and every other application continues to get
a clear, named error instead of either a silent global block or a
misleadingly "successful" run through unfinished machinery. See
``src/run_model_da/icesee_da_distributed.py`` for the runner that consults
this registry.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class Mode3Registration:
    """One application's production mode-3 entry point.

    ``runner`` is called exactly like the mode 0/1/2 runners
    (``runner(**icesee_kwargs)``); it is responsible for driving the
    application's full mode-3 DA cycle using
    ``src/parallelization/distributed_native_cycle.py`` and its own
    production-grade native adapter.
    """

    model_name: str
    runner: Callable[..., object]
    notes: str = ""


_REGISTRY: dict[str, Mode3Registration] = {}


def register_execution_mode_3(
    model_name: str,
    runner: Callable[..., object],
    *,
    notes: str = "",
) -> None:
    """Register ``model_name``'s production mode-3 runner.

    Intended to be called once, at application import time, by the
    application itself -- mirroring how ``SupportedModels`` already lets
    each application register its modes 0-2 model class. Re-registering the
    same ``model_name`` overwrites the previous entry rather than erroring,
    so an application module can be re-imported (for example under pytest)
    without needing extra guard logic.
    """

    if not model_name:
        raise ValueError("model_name must be a non-empty string")
    if not callable(runner):
        raise TypeError("runner must be callable")
    _REGISTRY[model_name] = Mode3Registration(
        model_name=model_name, runner=runner, notes=notes
    )


def get_execution_mode_3(model_name: str | None) -> Mode3Registration | None:
    """Return ``model_name``'s registration, or ``None`` if unsupported."""

    if not model_name:
        return None
    return _REGISTRY.get(model_name)


def registered_execution_mode_3_models() -> tuple[str, ...]:
    """Return the currently registered model names, sorted for display."""

    return tuple(sorted(_REGISTRY))
