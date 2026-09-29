# ==============================================================================
# @desc: Entry point for execution_mode 3 ("spatially distributed").
#
# Unlike modes 0-2, no single implementation here is model-agnostic in the
# same sense: mode 3 requires a production native adapter (genuinely
# distributed state; no rank holds a complete member) and a DA-cycle runner
# built on src/parallelization/distributed_native_cycle.py, per
# docs/execution-mode-3-design.md. This module does not implement that cycle
# itself -- it looks up the calling application's registration in
# src/parallelization/distributed_mode3_registry.py and delegates to it.
#
# Applications that have not registered support get a clear, named error
# here instead of ICESEE silently pretending mode 3 does not exist (the
# previous behavior, enforced by icesee_context.normalize_execution_mode
# rejecting any value other than 0, 1, or 2).
# ==============================================================================

from __future__ import annotations

import importlib
import importlib.util
import os
import re
import sys

from ICESEE.src.parallelization.distributed_mode3_registry import (
    get_execution_mode_3,
    registered_execution_mode_3_models,
)
from ICESEE.src.utils.icesee_context import normalize_execution_mode, normalize_icesee_kwargs


def _jit_import_mode3_runner(model_name):
    """Best-effort fallback registration (Stage 4D.3 production bootstrap fix).

    Registration is a module-level side effect of importing an
    application's own ``mode3_runner.py`` (e.g.
    ``applications/icepack_model/examples/idealized_pig/mode3_runner.py``,
    imported by that example's own ``run_da_icepack.py``). That makes it
    fragile to exactly the failure a real production run hit: any entry
    point, wrapper script, or stale checkout that reaches
    ``icesee_model_data_assimilation`` without having first imported (or
    having successfully imported -- an earlier, unrelated exception in
    that same script before the mode3_runner import line would silently
    prevent it too) that one specific line sees an empty registry, even
    though the real adapter exists in the repository and nothing about
    ``model_capabilities``/``supported_models`` is actually broken.

    This mirrors ``applications.supported_models.SupportedModels.
    call_model()``'s own existing, already-trusted dynamic-resolution
    convention exactly: ``<model>_model.examples.<cwd_dir_name>.
    mode3_runner``, with the same ``applications/`` directory on
    ``sys.path``. A missing module (``find_spec`` returns ``None``,  the
    normal case for every model that has no mode-3 support at all) is not
    an error -- only a genuine import-time exception inside an existing
    ``mode3_runner.py`` propagates, so a real bug there is still visible
    rather than silently masked.
    """
    if not model_name:
        return
    normalized_model = str(model_name).lower()
    if not re.match(r'^[a-zA-Z_][a-zA-Z0-9_]*$', normalized_model):
        return

    _here = os.path.dirname(os.path.abspath(__file__))
    _project_root = _here
    while not os.path.isdir(os.path.join(_project_root, "applications")):
        _parent = os.path.dirname(_project_root)
        if _parent == _project_root:
            return
        _project_root = _parent
    application_dir = os.path.join(_project_root, "applications")
    if application_dir not in sys.path:
        sys.path.insert(0, application_dir)

    dir_name = os.path.basename(os.getcwd())
    if not re.match(r'^[a-zA-Z_][a-zA-Z0-9_]*$', dir_name):
        return

    # Some examples' own internal packages (confirmed for idealized_pig's
    # "modelfunc": its __init__.py does a bare "from modelfunc.argusMesh
    # import argusMesh", not a relative or ICESEE.-qualified import)
    # assume the example's own directory is directly on sys.path -- which
    # normally happens for free only because Python auto-inserts a
    # directly-executed script's own directory at sys.path[0]. An
    # application entry point invoked any other way (imported as a
    # module, launched via a wrapper, `python -m`, ...) would never get
    # that for free, and mode3_runner.py's own import chain (which
    # reaches idealized_pig's _icepack_model.py -> modelfunc) would fail
    # with a plain ModuleNotFoundError for "modelfunc" instead of
    # registering -- reproduced directly while building this fallback.
    example_dir = os.path.join(
        application_dir, f"{normalized_model}_model", "examples", dir_name
    )
    if os.path.isdir(example_dir) and example_dir not in sys.path:
        sys.path.insert(0, example_dir)

    dynamic_module_path = f"{normalized_model}_model.examples.{dir_name}.mode3_runner"

    try:
        # find_spec raises ModuleNotFoundError (rather than returning
        # None) when the top-level package itself doesn't exist -- e.g.
        # an unrecognized model_name with no "<model>_model" package at
        # all, not just a model that has no mode3_runner in this example
        # directory. Both cases mean "no mode-3 support to discover
        # here", not an error.
        module_spec = importlib.util.find_spec(dynamic_module_path)
    except (ImportError, ModuleNotFoundError, ValueError):
        return
    if module_spec is None:
        return
    importlib.import_module(dynamic_module_path)


def icesee_model_data_assimilation_distributed(**icesee_kwargs):
    """Dispatch execution_mode 3 to the calling application's registration.

    Raises ``NotImplementedError`` naming the requested model and the
    currently registered models when no registration exists, rather than
    running unfinished machinery or failing with an unrelated error deep
    inside a partially wired path.
    """

    icesee_kwargs = normalize_icesee_kwargs(icesee_kwargs)
    normalize_execution_mode(icesee_kwargs, expected=3)
    model = icesee_kwargs.get("model_name")
    registration = get_execution_mode_3(model)
    if registration is None:
        # Self-heal before failing: see _jit_import_mode3_runner's own
        # docstring for why the registry can legitimately be empty here
        # even when a real adapter exists in the repository.
        _jit_import_mode3_runner(model)
        registration = get_execution_mode_3(model)
    if registration is None:
        supported = registered_execution_mode_3_models()
        raise NotImplementedError(
            f"execution_mode 3 is not yet available for model {model!r}. "
            f"Currently registered models: {list(supported) or 'none'}. "
            "Mode 3 requires a production native adapter and DA-cycle "
            "runner (see docs/execution-mode-3-design.md); modes 0, 1, and "
            "2 remain available for every application."
        )
    return registration.runner(**icesee_kwargs)
