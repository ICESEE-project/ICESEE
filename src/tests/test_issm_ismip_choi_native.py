"""Unit tests for the ISSM ISMIP_Choi mode-3 whole-member native wiring.

No MATLAB, ISSM, or real filesystem I/O is exercised here: the wiring
functions in ``_issm_native.py`` are pure translation between
``WholeMemberState`` and the existing ``initialize_ensemble``/
``forecast_step_single``/``run_model_inverse`` calling convention, so those
three functions are monkeypatched with in-memory fakes -- exactly the same
strategy used for the Icepack native adapter in
``test_icepack_idealized_pig_native.py``.

Import setup mirrors that same file: ``config/_utility_imports.py`` (pulled
in transitively by ``_issm_model.py``/``_issm_enkf.py``) parses ``sys.argv``
with ``argparse`` and loads ``params.yaml`` relative to the current working
directory, both at import time. Satisfy those here so this test runs under
plain ``pytest`` from the repo root like every other test.
"""

from __future__ import annotations

import os
import sys
import tempfile
import types
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_ISMIP_CHOI_DIR = (
    _REPO_ROOT / "applications" / "issm_model" / "examples" / "ISMIP_Choi"
)
for _extra_path in (str(_REPO_ROOT), str(_REPO_ROOT.parent), str(_ISMIP_CHOI_DIR)):
    if _extra_path not in sys.path:
        sys.path.insert(0, _extra_path)

_ARGV_BACKUP = sys.argv[:]
_CWD_BACKUP = os.getcwd()
sys.argv = [sys.argv[0], "--data_path", tempfile.mkdtemp(prefix="icesee_test_data_path_")]
os.chdir(_ISMIP_CHOI_DIR)
try:
    from applications.issm_model.examples.ISMIP_Choi import _issm_native as native
finally:
    sys.argv = _ARGV_BACKUP
    os.chdir(_CWD_BACKUP)

from src.parallelization.distributed_whole_member_adapter import WholeMemberState


def _topology():
    return types.SimpleNamespace(spatial_ranks=1)


_VEC_INPUTS = ["Thickness", "Surface", "Vx", "Vy"]


def _base_kwargs(**overrides):
    kwargs = {"vec_inputs": list(_VEC_INPUTS), "t": np.array([0.0, 0.2, 0.4])}
    kwargs.update(overrides)
    return kwargs


# --- initialize_member -------------------------------------------------


def test_initialize_member_delegates_to_initialize_ensemble(monkeypatch):
    calls = []

    def fake_initialize_ensemble(ens, **icesee_kwargs):
        calls.append((ens, icesee_kwargs.get("vec_inputs")))
        return {
            "Thickness": np.array([1.0, 2.0]),
            "Surface": np.array([3.0, 4.0]),
            "Vx": np.array([5.0, 6.0]),
            "Vy": np.array([7.0, 8.0]),
        }

    monkeypatch.setattr(native, "initialize_ensemble", fake_initialize_ensemble)

    state = native.initialize_member(
        3, topology=_topology(), icesee_kwargs=_base_kwargs()
    )

    assert isinstance(state, WholeMemberState)
    assert calls == [(3, _VEC_INPUTS)]
    assert state.model_context == {"ens_id": 3}
    np.testing.assert_array_equal(
        state.registry.pack_owned(), [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    )


def test_initialize_member_requires_vec_inputs(monkeypatch):
    monkeypatch.setattr(
        native, "initialize_ensemble", lambda ens, **kw: {"Thickness": np.zeros(1)}
    )
    with pytest.raises(ValueError):
        native.initialize_member(0, topology=_topology(), icesee_kwargs={})


# --- forecast_member -----------------------------------------------------


def test_forecast_member_packs_ensemble_and_returns_updates(monkeypatch):
    state = WholeMemberState(
        variables={
            "Thickness": np.array([1.0, 2.0]),
            "Surface": np.array([3.0, 4.0]),
            "Vx": np.array([5.0, 6.0]),
            "Vy": np.array([7.0, 8.0]),
        },
        vec_inputs=_VEC_INPUTS,
        model_context={"ens_id": 2},
    )

    captured = {}

    def fake_forecast_step_single(ensemble=None, **icesee_kwargs):
        captured["ensemble"] = ensemble
        captured["k"] = icesee_kwargs.get("k")
        captured["ens_id"] = icesee_kwargs.get("ens_id")
        return {
            "Thickness": ensemble[0:2] + 100.0,
            "Surface": ensemble[2:4] + 100.0,
            "Vx": ensemble[4:6] + 100.0,
            "Vy": ensemble[6:8] + 100.0,
        }

    monkeypatch.setattr(
        native, "forecast_step_single", fake_forecast_step_single
    )

    result = native.forecast_member(
        state, 4, topology=_topology(), icesee_kwargs=_base_kwargs()
    )

    np.testing.assert_array_equal(
        captured["ensemble"], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    )
    assert captured["k"] == 4
    assert captured["ens_id"] == 2
    np.testing.assert_array_equal(result["Thickness"], [101.0, 102.0])
    np.testing.assert_array_equal(result["Vy"], [107.0, 108.0])


def test_forecast_member_result_applies_through_adapter(monkeypatch):
    def fake_initialize_ensemble(ens, **icesee_kwargs):
        return {
            "Thickness": np.array([1.0, 2.0]),
            "Surface": np.array([3.0, 4.0]),
            "Vx": np.array([5.0, 6.0]),
            "Vy": np.array([7.0, 8.0]),
        }

    def fake_forecast_step_single(ensemble=None, **icesee_kwargs):
        return {name: ensemble[i * 2:(i + 1) * 2] + 1.0 for i, name in enumerate(_VEC_INPUTS)}

    monkeypatch.setattr(native, "initialize_ensemble", fake_initialize_ensemble)
    monkeypatch.setattr(native, "forecast_step_single", fake_forecast_step_single)

    adapter = native.ISMIP_CHOI_NATIVE_ADAPTER
    member = adapter.initialize_native_member(
        0, topology=_topology(), icesee_kwargs=_base_kwargs()
    )
    adapter.forecast_native_member(
        member, 0, topology=_topology(), icesee_kwargs=_base_kwargs()
    )
    np.testing.assert_array_equal(
        member.fields.pack_owned(), [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    )


# --- inverse_member -------------------------------------------------------


def test_inverse_member_delegates_to_run_model_inverse(monkeypatch):
    state = WholeMemberState(
        variables={
            "Thickness": np.array([1.0, 2.0]),
            "Surface": np.array([3.0, 4.0]),
            "Vx": np.array([5.0, 6.0]),
            "Vy": np.array([7.0, 8.0]),
        },
        vec_inputs=_VEC_INPUTS,
        model_context={"ens_id": 5},
    )

    captured = {}

    def fake_run_model_inverse(ensemble, **icesee_kwargs):
        captured["ensemble"] = ensemble
        captured["ens_id"] = icesee_kwargs.get("ens_id")
        return {"Vx": ensemble[4:6] + 10.0, "Vy": ensemble[6:8] + 10.0}

    monkeypatch.setattr(native, "run_model_inverse", fake_run_model_inverse)

    result = native.inverse_member(
        state, 0, topology=_topology(), icesee_kwargs=_base_kwargs(km=0)
    )

    assert captured["ens_id"] == 5
    np.testing.assert_array_equal(result["Vx"], [15.0, 16.0])
    assert "Thickness" not in result


def test_inverse_member_none_result_is_passed_through(monkeypatch):
    state = WholeMemberState(
        variables={
            "Thickness": np.array([1.0]),
            "Surface": np.array([1.0]),
            "Vx": np.array([1.0]),
            "Vy": np.array([1.0]),
        },
        vec_inputs=_VEC_INPUTS,
        model_context={"ens_id": 0},
    )
    monkeypatch.setattr(
        native, "run_model_inverse", lambda ensemble, **kw: None
    )
    result = native.inverse_member(
        state, 0, topology=_topology(), icesee_kwargs=_base_kwargs(km=0)
    )
    assert result is None


def test_adapter_exposes_inversion_capability():
    assert hasattr(native.ISMIP_CHOI_NATIVE_ADAPTER, "inverse_native_member")
