# ==============================================================================
# @des: Level-1 tests for the top-level execution-mode dispatcher in
# src/run_model_da/run_models_da.py, covering the mode-3 recovery: modes
# 0/1/2 route unchanged, mode 3 routes to icesee_da_distributed when
# registered and fails with a clear, controlled error when it is not.
#
# run_models_da.py itself is safe to import directly (it only imports
# icesee_context, which has no CLI/import-time side effects); the actual
# mode 0/1/2 driver modules are imported lazily inside
# icesee_model_data_assimilation only once a mode is dispatched, so this
# file never needs the config/_utility_imports argv/data_path safety dance
# used elsewhere -- it dispatches only mode 3 (icesee_da_distributed.py has
# no such side effects either) end-to-end, and checks modes 0/1/2 at the
# resolution/target-table level without actually invoking their drivers.
# ==============================================================================
import contextlib
import importlib
import os
import sys
import tempfile

import pytest

from ICESEE.src.run_model_da.run_models_da import (
    _MODE_TO_TARGET,
    _resolve_mode,
    icesee_model_data_assimilation,
)
from ICESEE.src.parallelization.distributed_mode3_registry import (
    register_execution_mode_3,
)

# icesee_da_serial/partial/full import config/_utility_imports transitively
# (via EnKF.py), which parses sys.argv and touches data_path at import time
# -- see test_enkf_serial_process_noise.py for the same pattern. Only
# needed for the importability check below; icesee_da_distributed.py has
# no such side effects.
_lorenz96_params = os.path.join(
    os.path.dirname(__file__),
    "..", "..", "applications", "lorenz_model", "examples", "lorenz96", "params.yaml",
)


def _import_module_with_safe_argv(module_name):
    safe_data_path = tempfile.mkdtemp(prefix="icesee_run_models_da_dispatch_test_")
    saved_argv = sys.argv
    sys.argv = [sys.argv[0], "-F", _lorenz96_params, "--data_path", safe_data_path]
    try:
        return importlib.import_module(module_name)
    finally:
        sys.argv = saved_argv


@pytest.fixture(autouse=True)
def _clean_registry():
    from ICESEE.src.parallelization import distributed_mode3_registry as module

    before = set(module._REGISTRY)
    yield
    for name in set(module._REGISTRY) - before:
        del module._REGISTRY[name]


@pytest.mark.parametrize(
    "mode_number,expected_key",
    [(0, "serial"), (1, "partial"), (2, "full"), (3, "distributed")],
)
def test_resolve_mode_is_unchanged_for_0_1_2_and_adds_3(mode_number, expected_key):
    assert _resolve_mode({"execution_mode": mode_number}) == expected_key


@pytest.mark.parametrize(
    "key,module_name,func_name",
    [
        ("serial", "ICESEE.src.run_model_da.icesee_da_serial", "icesee_model_data_assimilation_serial"),
        ("partial", "ICESEE.src.run_model_da.icesee_da_partial_parallel", "icesee_model_data_assimilation_partial_parallel"),
        ("full", "ICESEE.src.run_model_da.icesee_da_full_parallel", "icesee_model_data_assimilation_full_parallel"),
        ("distributed", "ICESEE.src.run_model_da.icesee_da_distributed", "icesee_model_data_assimilation_distributed"),
    ],
)
def test_mode_to_target_entries_are_unchanged_and_importable(key, module_name, func_name):
    assert _MODE_TO_TARGET[key] == (module_name, func_name)
    mod = _import_module_with_safe_argv(module_name)
    assert callable(getattr(mod, func_name))


def test_mode_3_registered_dispatches_through_icesee_model_data_assimilation():
    received = {}

    def fake_runner(**kwargs):
        received.update(kwargs)
        return "mode-3-ran"

    register_execution_mode_3("dispatch-test-model", fake_runner)

    result = icesee_model_data_assimilation(
        execution_mode=3,
        model_name="dispatch-test-model",
        Nens=3,
        batch_size=1,
    )

    assert result == "mode-3-ran"
    assert received["model_name"] == "dispatch-test-model"
    assert received["execution_mode"] == 3


def test_mode_3_unregistered_fails_clearly_not_silently():
    with pytest.raises(NotImplementedError, match="dispatch-test-unregistered-model"):
        icesee_model_data_assimilation(
            execution_mode=3,
            model_name="dispatch-test-unregistered-model",
            Nens=3,
            batch_size=1,
        )


def test_execution_mode_4_is_rejected():
    with pytest.raises(ValueError):
        icesee_model_data_assimilation(execution_mode=4, model_name="anything")


# --- Stage 4D.3 production bootstrap fix: self-healing mode-3 registration ---
#
# A real production run reached "Execution mode 3: spatially distributed
# (application-specific)" and then failed with "Currently registered
# models: none" even though idealized_pig's own mode3_runner.py (which
# registers "icepack" as a module-level side effect of being imported)
# exists and is unmodified. Root cause could not be reproduced exactly in
# this sandbox (a separate, unrelated Firedrake/UFL API-version crash in
# idealized_pig's own run_da_icepack.py -- "Non-contiguous argument
# numbers in interpolate" -- blocked the real script from ever reaching
# the mode-3 dispatch call at all here, before mode3_runner.py's own
# import could even matter); every import-order/import-prefix hypothesis
# tested directly (manual import, full-bootstrap-matching import order,
# with/without the "ICESEE." prefix) registered correctly. Given
# registration's only failure mode that DOES reproduce is "some entry
# point reaches icesee_model_data_assimilation without mode3_runner.py's
# side-effect import having happened" -- exactly what these tests
# simulate by calling the dispatcher directly, cwd'd into the real
# example directory, without ever importing mode3_runner first -- the fix
# is a self-healing fallback (_jit_import_mode3_runner in
# icesee_da_distributed.py) rather than a one-off patch to a single
# import line, so it is robust to whichever exact entry point or stale-
# import scenario caused the original failure.
#
# These tests exercise the REAL application bootstrap convention (cwd +
# model_name -> <model>_model.examples.<dir>.mode3_runner, mirroring
# SupportedModels.call_model()'s own dynamic resolution), not a manual
# "import mode3_runner first" shortcut.

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _real_bootstrap_registers(model_name, example_dir):
    """cd into a real example directory and dispatch mode 3 directly,
    without ever having imported that example's mode3_runner.py --
    returns whatever exception the dispatch raised (None if it returned).
    """
    saved_cwd = os.getcwd()
    saved_argv = sys.argv
    # mode3_runner.py's import chain transitively reaches
    # config/_utility_imports.py, which parses sys.argv for --key=value
    # CLI overrides at import time (see test_icepack_mode3_runner.py's own
    # identical --data_path-only argv override for the same reason) --
    # pytest's own argv (test paths, -q, ...) is not a recognized ICESEE
    # key and would otherwise raise there, before mode3_runner.py's own
    # module-level register_execution_mode_3 call is ever reached.
    sys.argv = [sys.argv[0], "--data_path", tempfile.mkdtemp(
        prefix="icesee_mode3_bootstrap_test_"
    )]
    os.chdir(example_dir)
    try:
        from ICESEE.src.run_model_da.icesee_da_distributed import (
            icesee_model_data_assimilation_distributed,
        )
        try:
            icesee_model_data_assimilation_distributed(
                model_name=model_name, execution_mode=3,
            )
            return None
        except Exception as e:  # noqa: BLE001 -- deliberately broad, see below
            return e
    finally:
        os.chdir(saved_cwd)
        sys.argv = saved_argv


@contextlib.contextmanager
def _temporarily_unregistered(model_name):
    """Remove a model from the mode-3 registry for the duration of the
    `with` block, restoring whatever was there before on exit (even on
    failure) so this test never leaks a permanent registry change into
    whichever sibling test happens to run after it.

    Some sibling test files (test_icepack_mode3_runner.py, test_issm_
    mode3_runner.py) import their example's mode3_runner.py at their own
    module TOP LEVEL (collection time, outside any fixture) purely to
    exercise its runner in isolation -- that registration is never
    cleaned up and persists for the rest of the pytest session, so this
    file's own _clean_registry fixture (which only removes entries ADDED
    during one test) cannot guarantee these tests' "starts unregistered"
    precondition when run as part of the full suite. Reproduced directly:
    these tests pass in isolation but failed when the full suite ran
    test_icepack_mode3_runner.py first.
    """
    from ICESEE.src.parallelization import distributed_mode3_registry as module

    previous = module._REGISTRY.pop(model_name, None)
    try:
        yield
    finally:
        module._REGISTRY.pop(model_name, None)
        if previous is not None:
            module._REGISTRY[model_name] = previous


def test_icepack_real_bootstrap_self_heals_without_a_prior_mode3_runner_import():
    example_dir = os.path.join(
        _REPO_ROOT, "applications", "icepack_model", "examples", "idealized_pig",
    )
    from ICESEE.src.parallelization.distributed_mode3_registry import get_execution_mode_3

    with _temporarily_unregistered("icepack"):
        assert get_execution_mode_3("icepack") is None

        exc = _real_bootstrap_registers("icepack", example_dir)

        # The real point of this test: registration must have happened.
        # What the dispatch raised AFTER that (a downstream error from
        # calling the real runner with a deliberately minimal
        # icesee_kwargs, since a full 37GB idealized_pig run is not
        # appropriate for a unit test) is not itself asserted beyond "it
        # is not the registration failure".
        assert not isinstance(exc, NotImplementedError), (
            f"mode 3 dispatch still failed to register icepack via the real "
            f"bootstrap path: {exc!r}"
        )
        assert get_execution_mode_3("icepack") is not None


def test_issm_real_bootstrap_self_heals_without_a_prior_mode3_runner_import():
    example_dir = os.path.join(
        _REPO_ROOT, "applications", "issm_model", "examples", "ISMIP_Choi",
    )
    from ICESEE.src.parallelization.distributed_mode3_registry import get_execution_mode_3

    with _temporarily_unregistered("issm"):
        assert get_execution_mode_3("issm") is None

        exc = _real_bootstrap_registers("issm", example_dir)

        assert not isinstance(exc, NotImplementedError), (
            f"mode 3 dispatch still failed to register issm via the real "
            f"bootstrap path: {exc!r}"
        )
        assert get_execution_mode_3("issm") is not None


def test_genuinely_unsupported_model_still_fails_clearly_via_real_bootstrap():
    # A model with no "<model>_model" package at all must still raise
    # NotImplementedError, not some other exception leaking out of the
    # self-heal fallback's own import machinery (regression guard for the
    # find_spec-raises-ModuleNotFoundError-for-a-missing-top-level-package
    # edge case found and fixed while writing _jit_import_mode3_runner).
    example_dir = os.path.join(
        _REPO_ROOT, "applications", "icepack_model", "examples", "idealized_pig",
    )
    exc = _real_bootstrap_registers("totally_unsupported_model_xyz", example_dir)
    assert isinstance(exc, NotImplementedError)


def test_model_with_no_mode3_runner_in_this_example_dir_still_fails_clearly():
    # A real, supported model (icepack) run from a real example directory
    # that has no mode3_runner.py of its own (synthetic_ice_stream) must
    # still raise NotImplementedError, not silently pick up a different
    # example's registration or crash differently.
    example_dir = os.path.join(
        _REPO_ROOT, "applications", "icepack_model", "examples", "synthetic_ice_stream",
    )
    from ICESEE.src.parallelization.distributed_mode3_registry import get_execution_mode_3

    with _temporarily_unregistered("icepack"):
        assert get_execution_mode_3("icepack") is None
        exc = _real_bootstrap_registers("icepack", example_dir)
        assert isinstance(exc, NotImplementedError)
