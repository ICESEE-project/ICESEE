# ==============================================================================
# @des: Regression test for the observation-schedule nt architecture fix
# (config/_utility_imports.py, 2026-09-28, second pass).
#
# Background: different applications interpret `timesteps_per_year`
# differently (Icepack/Lorenz96/ISSM's own run_da_*.py all self-derive
# nt=num_years/timesteps_per_year in their OWN code, well AFTER config
# load -- but config/_utility_imports.py computes its own one-time
# internal nt/t ONCE, at config-load time, purely to feed the generic
# generate_observation_schedule() call, before any application gets a
# chance to self-derive anything). A first pass "fixed" this generic
# nt by switching its formula from multiplication to division -- but
# that just swapped one hardcoded, application-blind assumption for
# another (confirmed: this changed ISSM's and synthetic_ice_stream's
# internal nt here too, since they never asked for that).
#
# The actual fix: config/_utility_imports.py never assumes a
# convention. An application's own params.yaml may declare its
# already-resolved `nt` directly (`modeling-parameters.nt`) -- an
# explicit, application-owned value, not a generic reinterpretation of
# `timesteps_per_year`. Without that key, the ORIGINAL historical
# formula (num_years * timesteps_per_year) is used, so every existing
# configuration that does not opt in behaves EXACTLY as it always has.
# No model-name conditionals anywhere in this file.
#
# These tests exercise the REAL config loader end to end for each
# application's REAL params.yaml (not a duplicated formula), confirming:
#   1. Idealized PIG (which opts in via its own `nt: 1640`) keeps its
#      correct 82-year/1640-step production schedule.
#   2. ISSM (both examples) and synthetic_ice_stream (which do NOT opt
#      in) get back their exact ORIGINAL, pre-2026-09-28 nt -- unchanged
#      by any of this reconciliation's work, despite touching zero ISSM
#      files.
#   3. The nt-selection mechanism itself is generic: an explicit `nt`
#      always wins, and its absence always falls back to the historical
#      formula, for a synthetic (non-real-application) configuration too.
# ==============================================================================
from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_IDEALIZED_PIG_DIR = _REPO_ROOT / "applications" / "icepack_model" / "examples" / "idealized_pig"
_ISMIP_CHOI_DIR = _REPO_ROOT / "applications" / "issm_model" / "examples" / "ISMIP_Choi"
_BASAL_FRICTION_DIR = _REPO_ROOT / "applications" / "issm_model" / "examples" / "basal_friction_variation"
_SYNTHETIC_ICE_STREAM_DIR = _REPO_ROOT / "applications" / "icepack_model" / "examples" / "synthetic_ice_stream"


def _load_real_kwargs(app_dir, tmp_path, **cli_overrides):
    """Load a REAL params.yaml, from a REAL application directory, through
    the REAL generic config loader -- not a duplicated formula. Runs
    in-process (config loading itself doesn't import Firedrake/ISSM),
    matching this suite's existing convention for config-loader-only
    tests (e.g. test_icepack_compact_initialization.py)."""

    _argv_backup = sys.argv[:]
    _cwd_backup = os.getcwd()
    argv = [sys.argv[0], "--data_path", str(tmp_path)]
    for key, value in cli_overrides.items():
        argv.append(f"--{key}={value}")
    sys.argv = argv
    os.chdir(app_dir)
    try:
        # Force the module-level config loader to re-run each call via
        # importlib.reload (in place) rather than sys.modules.pop + a
        # fresh import -- popping replaces the sys.modules entry with a
        # brand-new module object, which orphans any reference to the
        # old one that other test modules captured at collection time
        # (e.g. test_utility_imports_cli_overrides.py's module-level
        # `cli_module`), breaking their own later importlib.reload calls.
        # reload() re-executes the module body but preserves identity.
        module_name = "ICESEE.config._utility_imports"
        if module_name in sys.modules:
            module = importlib.reload(sys.modules[module_name])
        else:
            module = importlib.import_module(module_name)
        return dict(module.icesee_kwargs)
    finally:
        sys.argv = _argv_backup
        os.chdir(_cwd_backup)


# --- 1. Idealized PIG: production schedule preserved via its own explicit nt ---

def test_idealized_pig_production_keeps_nt_1640(tmp_path):
    kwargs = _load_real_kwargs(_IDEALIZED_PIG_DIR, tmp_path)
    assert kwargs["num_years"] == 82
    assert kwargs["timesteps_per_year"] == 0.05
    assert kwargs["nt"] == 1640, (
        "Idealized PIG declares nt explicitly in its own params.yaml "
        "(modeling-parameters.nt) -- this must never depend on the "
        "generic loader's fallback formula."
    )


def test_idealized_pig_production_observation_window_gives_31_events_at_years_10_to_40(tmp_path):
    kwargs = _load_real_kwargs(_IDEALIZED_PIG_DIR, tmp_path)
    assert kwargs["obs_start_time"] == 10
    assert kwargs["obs_max_time"] == 40
    assert kwargs["freq_obs"] == 1

    obs_index = np.asarray(kwargs["obs_index"])
    assert kwargs["number_obs_instants"] == 31
    assert len(obs_index) == 31

    # step = freq_obs / timesteps_per_year = 1 / 0.05 = 20 timesteps/year
    expected_index = np.arange(200, 801, 20)
    np.testing.assert_array_equal(obs_index, expected_index)

    t = np.asarray(kwargs["t"])
    observed_years = t[obs_index]
    np.testing.assert_allclose(observed_years, np.arange(10, 41, 1), atol=1e-9)


def test_idealized_pig_calibration_config_gives_nt_40_and_one_event_at_year_1(tmp_path):
    kwargs = _load_real_kwargs(
        _IDEALIZED_PIG_DIR, tmp_path,
        **{"force-params": str(_IDEALIZED_PIG_DIR / "pace" / "params_calibration.yaml")},
    )
    assert kwargs["nt"] == 40
    obs_index = np.asarray(kwargs["obs_index"])
    assert kwargs["number_obs_instants"] == 1
    np.testing.assert_array_equal(obs_index, [20])
    t = np.asarray(kwargs["t"])
    np.testing.assert_allclose(t[obs_index], [1.0], atol=1e-9)


# --- 2. ISSM/synthetic_ice_stream: NOT touched, NOT opted in, behavior unchanged ---
# (config-loading only -- no ISSM/Firedrake application code is imported or run;
# this reads real params.yaml files nobody has edited, exactly as before.)

def test_issm_ismip_choi_internal_nt_is_unchanged_from_before_this_reconciliation(tmp_path):
    kwargs = _load_real_kwargs(_ISMIP_CHOI_DIR, tmp_path)
    assert "nt" not in _read_modeling_section_keys(_ISMIP_CHOI_DIR), (
        "this test's premise is that ISSM has NOT opted in -- if this "
        "fails, someone added an `nt` key to ISSM's own params.yaml"
    )
    assert kwargs["num_years"] == 60
    assert kwargs["timesteps_per_year"] == 0.2
    # Original historical formula (num_years * timesteps_per_year) --
    # unaffected by this reconciliation. This is config/_utility_imports.py's
    # OWN internal nt (feeds only the generic observation-schedule call);
    # ISSM's REAL simulation nt is computed separately, in ISSM's own
    # run_da_issm.py (nt=(num_years-tinitial)/timesteps_per_year=300),
    # which this generic loader has never touched and does not affect.
    assert kwargs["nt"] == 12


def test_issm_basal_friction_variation_internal_nt_is_unchanged(tmp_path):
    kwargs = _load_real_kwargs(_BASAL_FRICTION_DIR, tmp_path)
    assert kwargs["num_years"] == 200
    assert kwargs["timesteps_per_year"] == 0.2
    assert kwargs["nt"] == 40  # int(200 * 0.2), original historical formula


def test_synthetic_ice_stream_internal_nt_is_unchanged(tmp_path):
    kwargs = _load_real_kwargs(_SYNTHETIC_ICE_STREAM_DIR, tmp_path)
    assert kwargs["num_years"] == 120
    assert kwargs["timesteps_per_year"] == 2
    assert kwargs["nt"] == 240  # int(120 * 2), original historical formula -- unaffected


def _read_modeling_section_keys(app_dir):
    import yaml
    with open(app_dir / "params.yaml") as f:
        return set(yaml.safe_load(f).get("modeling-parameters", {}).keys())


# --- 3. The mechanism itself is generic: explicit nt always wins, absence
#        always falls back to the historical formula. No model-name
#        conditionals -- this is proven directly against the loader's own
#        source, and against a synthetic (non-real-application) YAML. ---

def test_nt_selection_is_generic_not_a_model_name_conditional_in_source():
    source = (
        Path(__file__).resolve().parents[2] / "config" / "_utility_imports.py"
    ).read_text()
    assert "'nt' in _modeling_section" in source, (
        "nt must be selected by whether the application's own YAML "
        "declares it, not by inspecting which model/application is running"
    )
    for banned in ("model == \"icepack\"", "model=='icepack'", "if model ==", 'model_name ==', "== 'icepack'"):
        assert banned not in source.split("def apply_generic_cli_overrides")[0][:6000], (
            f"generic nt selection must not contain a model-name conditional near {banned!r}"
        )


def test_explicit_nt_key_wins_over_the_historical_formula(tmp_path):
    """A synthetic, non-real-application config: proves the opt-in
    mechanism itself, independent of any specific application."""

    import yaml
    app_dir = tmp_path / "synthetic_app"
    app_dir.mkdir()
    params = {
        "modeling-parameters": {
            "num_years": 10,
            "timesteps_per_year": 3,  # historical formula would give nt=30
            "nt": 999,  # explicit override must win instead
        },
        "enkf-parameters": {
            "Nens": 1, "freq_obs": 1, "obs_start_time": 0, "obs_max_time": 0,
            "num_state_vars": 1, "num_param_vars": 0, "vec_inputs": ["x"],
            "model_name": "synthetic", "model_solver": "none", "filter_type": "EnKF",
            "seed": 1, "data_path": str(tmp_path / "out"),
            "joint_estimation": 0, "parameter_estimation": False,
            "joint_estimated_params": [], "execution_mode": 0,
        },
    }
    with open(app_dir / "params.yaml", "w") as f:
        yaml.safe_dump(params, f)

    kwargs = _load_real_kwargs(app_dir, tmp_path / "out")
    assert kwargs["nt"] == 999


def test_absent_nt_key_falls_back_to_historical_formula(tmp_path):
    import yaml
    app_dir = tmp_path / "synthetic_app_no_nt"
    app_dir.mkdir()
    params = {
        "modeling-parameters": {
            "num_years": 10,
            "timesteps_per_year": 3,  # no explicit nt -> falls back to 10*3=30
        },
        "enkf-parameters": {
            "Nens": 1, "freq_obs": 1, "obs_start_time": 0, "obs_max_time": 0,
            "num_state_vars": 1, "num_param_vars": 0, "vec_inputs": ["x"],
            "model_name": "synthetic", "model_solver": "none", "filter_type": "EnKF",
            "seed": 1, "data_path": str(tmp_path / "out2"),
            "joint_estimation": 0, "parameter_estimation": False,
            "joint_estimated_params": [], "execution_mode": 0,
        },
    }
    with open(app_dir / "params.yaml", "w") as f:
        yaml.safe_dump(params, f)

    kwargs = _load_real_kwargs(app_dir, tmp_path / "out2")
    assert kwargs["nt"] == 30
