# ==============================================================================
# @des: Tests for the historical CLI-override and data_path auto-clean
# functionality ported from test/ICESEE's config/_utility_imports.py into
# main during Pass 2 reconciliation.
#
# config/_utility_imports.py executes its entire CLI-parsing/YAML-loading
# pipeline at IMPORT TIME (it is module-level code, not a callable
# function), and that pipeline now also deletes and recreates
# icesee_kwargs['data_path'] as a side effect. To test it safely:
#   - the module is imported exactly once at collection time, with argv
#     pointed at Lorenz96's example params.yaml (a real, complete config)
#     and --data_path redirected to a disposable temp directory, so the
#     auto-clean guard never touches anything real;
#   - the pure helper functions (_parse_generic_cli_overrides,
#     _coerce_cli_override, apply_generic_cli_overrides) are then tested
#     directly against synthetic dicts -- they have no side effects of
#     their own;
#   - the few tests that need the full module pipeline under different
#     argv (execution_mode CLI precedence, the data_path guard itself) use
#     importlib.reload with a fresh --data_path/--execution_mode each time.
# ==============================================================================
import importlib
import os
import shutil
import sys
import tempfile

import numpy as np
import pytest

_lorenz96_params = os.path.join(
    os.path.dirname(__file__),
    "..", "..", "applications", "lorenz_model", "examples", "lorenz96", "params.yaml",
)


def _fresh_tmp_data_path():
    path = tempfile.mkdtemp(prefix="icesee_utility_imports_test_")
    shutil.rmtree(path)  # exercise the "missing directory" case by default
    return path


_initial_data_path = _fresh_tmp_data_path()
_saved_argv = sys.argv
sys.argv = [sys.argv[0], "-F", _lorenz96_params, "--data_path", _initial_data_path]
try:
    from ICESEE.config import _utility_imports as cli_module
finally:
    sys.argv = _saved_argv


def _reload_with_argv(argv_tail):
    saved = sys.argv
    sys.argv = [sys.argv[0], "-F", _lorenz96_params] + argv_tail
    try:
        importlib.reload(cli_module)
    finally:
        sys.argv = saved


# --- _coerce_cli_override: representative types --------------------------

def test_coerce_int():
    assert cli_module._coerce_cli_override("40", 30) == 40
    assert isinstance(cli_module._coerce_cli_override("40", 30), int)


def test_coerce_fractional_value_is_not_truncated_by_an_integer_default():
    # e.g. params.yaml obs_start_time: 10 (an int) with --obs_start_time=0.25
    assert cli_module._coerce_cli_override("0.25", 10) == 0.25
    assert cli_module._coerce_cli_override("2.0", 10) == 2
    assert isinstance(cli_module._coerce_cli_override("2.0", 10), int)


def test_coerce_float():
    assert cli_module._coerce_cli_override("1.5", 1.0) == 1.5


def test_coerce_bool_true_and_false():
    assert cli_module._coerce_cli_override("true", False) is True
    assert cli_module._coerce_cli_override("false", True) is False


def test_coerce_string():
    assert cli_module._coerce_cli_override("fft", "graph") == "fft"


def test_coerce_none():
    # No isinstance branch matches a None current_value, so the parsed
    # value is returned as-is -- this is the historical fallback behavior,
    # not a special-cased None handler.
    assert cli_module._coerce_cli_override("null", None) is None
    assert cli_module._coerce_cli_override("5", None) == 5


def test_coerce_list_into_ndarray_preserves_dtype():
    current = np.array([0.01, 0.01, 0.01], dtype=np.float64)
    result = cli_module._coerce_cli_override("[0.02, 0.03, 0.04]", current)
    assert isinstance(result, np.ndarray)
    assert result.dtype == current.dtype
    np.testing.assert_allclose(result, [0.02, 0.03, 0.04])


def test_coerce_malformed_yaml_falls_back_to_raw_string():
    # Historical behavior: an unparsable value degrades to the literal
    # string rather than raising, and (since the current value here is a
    # string) is returned unchanged by the final fallback branch.
    result = cli_module._coerce_cli_override("[1,2", "placeholder")
    assert result == "[1,2"


# --- _parse_generic_cli_overrides ----------------------------------------

def test_parse_key_equals_value_syntax():
    overrides = cli_module._parse_generic_cli_overrides(["--Nens=40", "--sig_Q=[0.02,0.02]"])
    assert overrides == {"Nens": "40", "sig_Q": "[0.02,0.02]"}


def test_parse_key_space_value_syntax():
    overrides = cli_module._parse_generic_cli_overrides(["--Nens", "40"])
    assert overrides == {"Nens": "40"}


def test_parse_bare_flag_defaults_to_true():
    overrides = cli_module._parse_generic_cli_overrides(["--verbose"])
    assert overrides == {"verbose": "true"}


# --- apply_generic_cli_overrides: policy ----------------------------------

def test_apply_overrides_existing_keys_only():
    icesee_kwargs = {"Nens": 30, "verbose": False}
    result = cli_module.apply_generic_cli_overrides(dict(icesee_kwargs), ["--Nens=40", "--verbose"])
    assert result == {
        "Nens": 40,
        "verbose": True,
        # The explicitly requested values are recorded for later checks.
        "cli_overrides": {"Nens": 40, "verbose": True},
    }


def test_apply_overrides_rejects_unknown_key():
    # Overrides apply only to keys that already exist in icesee_kwargs;
    # this policy (typo protection) is preserved unchanged from the
    # historical implementation, not silently relaxed to allow new keys.
    with pytest.raises(ValueError):
        cli_module.apply_generic_cli_overrides({"Nens": 30}, ["--nense=40"])


# --- execution_mode CLI precedence ----------------------------------------

def test_execution_mode_cli_flag_can_represent_mode_3():
    # The parser can REPRESENT mode 3 (it is a real, documented value that
    # normalize_execution_mode already accepts); this does not claim mode
    # 3 is production-ready for Lorenz96 -- that is a separate, per-model
    # registry decision made downstream, not by this loader.
    tmp = _fresh_tmp_data_path()
    try:
        _reload_with_argv(["--data_path", tmp, "--execution_mode", "3"])
        assert cli_module.icesee_kwargs["execution_mode"] == 3
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_execution_mode_falls_back_to_yaml_when_cli_omitted():
    tmp = _fresh_tmp_data_path()
    try:
        _reload_with_argv(["--data_path", tmp])
        # Lorenz96's example params.yaml default execution_mode (this repo's
        # own value, independent of the CLI-precedence behavior under test).
        assert cli_module.icesee_kwargs["execution_mode"] == 2
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_execution_mode_rejects_out_of_range_value():
    tmp = _fresh_tmp_data_path()
    try:
        with pytest.raises(SystemExit):
            _reload_with_argv(["--data_path", tmp, "--execution_mode", "4"])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
        # argparse's choices rejection calls sys.exit before rebuilding
        # icesee_kwargs; reload the module back to a good state for any
        # later test in this session.
        _reload_with_argv(["--data_path", _fresh_tmp_data_path()])


# --- data_path auto-clean guard -------------------------------------------

def test_data_path_stale_file_is_removed_and_directory_recreated():
    tmp = tempfile.mkdtemp(prefix="icesee_utility_imports_test_")
    stale_file = os.path.join(tmp, "stale_from_previous_run.h5")
    with open(stale_file, "w") as handle:
        handle.write("stale")
    try:
        _reload_with_argv(["--data_path", tmp])
        assert os.path.isdir(tmp)
        assert not os.path.exists(stale_file)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


@pytest.mark.parametrize(
    "resume_flag", [["--resume_from_checkpoint=True"], ["--resume_from_checkpoint"]]
)
def test_resume_from_checkpoint_keeps_existing_data_path(resume_flag):
    tmp = tempfile.mkdtemp(prefix="icesee_utility_imports_test_")
    checkpoint = os.path.join(tmp, "_mode3_state_history", "steps", "checkpoint_00000525")
    os.makedirs(checkpoint)
    with open(os.path.join(checkpoint, "manifest.json"), "w") as handle:
        handle.write("{}")
    try:
        _reload_with_argv(["--data_path", tmp] + resume_flag)
        assert os.path.isfile(os.path.join(checkpoint, "manifest.json"))
        assert cli_module.icesee_kwargs["resume_from_checkpoint"] is True
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
        _reload_with_argv(["--data_path", _fresh_tmp_data_path()])


def test_checkpoint_retention_defaults_are_bounded():
    kwargs = cli_module.icesee_kwargs
    assert kwargs["resume_from_checkpoint"] is False
    assert kwargs["checkpoint_every"] == 1
    assert kwargs["checkpoint_keep_last"] == 2
    assert kwargs["checkpoint_keep_analysis"] is False


def test_data_path_missing_directory_is_created_without_error():
    tmp = _fresh_tmp_data_path()  # already removed by the helper
    assert not os.path.exists(tmp)
    try:
        _reload_with_argv(["--data_path", tmp])
        assert os.path.isdir(tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


@pytest.mark.parametrize("unsafe_path", [os.getcwd(), os.path.expanduser("~"), "/"])
def test_data_path_refuses_unsafe_paths(unsafe_path):
    with pytest.raises(ValueError):
        _reload_with_argv(["--data_path", unsafe_path])
    # Restore a good state for any later test in this session.
    _reload_with_argv(["--data_path", _fresh_tmp_data_path()])


def test_data_path_empty_string_does_not_silently_resolve_to_cwd_or_delete_it():
    # An empty data_path resolves (os.path.abspath('')) to the current
    # working directory, which the unsafe-path guard must catch -- not a
    # special case, the same guard that protects against '.' explicitly.
    with pytest.raises(ValueError):
        _reload_with_argv(["--data_path", ""])
    _reload_with_argv(["--data_path", _fresh_tmp_data_path()])
