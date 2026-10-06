# ==============================================================================
# @des: CI guard for Idealized PIG compact initialization, using a small
#       synthetic checkpoint instead of the ~39GB extended_beta1000yrs.h5.
#
#       The production initFile stores velocity/thickness/surface as a spin-up
#       TIME SERIES; initializeRun only needs the final snapshot (idx=20000).
#       tools/build_compact_initialization.py extracts that snapshot into a
#       small non-timestepped file, and initializeRun(compact_initialization=
#       True) reads it without idx=. This module checks, on a synthetic
#       checkpoint with several history indices holding DIFFERENT values:
#         * the compact builder selects exactly the final (idx=20000) state;
#         * the real initializeRun produces bit-identical initialized fields
#           from the full-history path and from the compact path;
#         * the compact path works with the full-history file absent;
#         * the compact file does not grow with the length of the history.
#       Only readSMB (GeoTIFF input) and the diagnostic solve are stubbed; the
#       checkpoint reads, initialState, and flotationHeight run for real.
#       test_icepack_compact_initialization.py keeps the equivalent check
#       against the real dataset for machines that have it.
# ==============================================================================
from __future__ import annotations

import importlib.util
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

firedrake = pytest.importorskip("firedrake")
pytest.importorskip("icepack")

_REPO_ROOT = Path(__file__).resolve().parents[2]
_IDEALIZED_PIG_DIR = _REPO_ROOT / "applications" / "icepack_model" / "examples" / "idealized_pig"
_TOOLS_DIR = _IDEALIZED_PIG_DIR / "tools"
for _extra_path in (str(_REPO_ROOT), str(_REPO_ROOT.parent), str(_IDEALIZED_PIG_DIR)):
    if _extra_path not in sys.path:
        sys.path.insert(0, _extra_path)

_FINAL_IDX = 20000
_HISTORY = (0, 1, 7, _FINAL_IDX)

_ARGV_BACKUP = sys.argv[:]
_CWD_BACKUP = os.getcwd()
try:
    # _icepack_model imports the config loader, which parses sys.argv and
    # auto-cleans data_path at import time: read idealized_pig's params.yaml
    # but point data_path at a throwaway directory so no run output is touched.
    sys.argv = [sys.argv[0], f"--data_path={tempfile.mkdtemp(prefix='icesee_compact_init_')}"]
    os.chdir(_IDEALIZED_PIG_DIR)
    from applications.icepack_model.examples.idealized_pig import _icepack_model as model
finally:
    sys.argv = _ARGV_BACKUP
    os.chdir(_CWD_BACKUP)


def _load_build_module():
    spec = importlib.util.spec_from_file_location(
        "build_compact_initialization", _TOOLS_DIR / "build_compact_initialization.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_full_history(path, history):
    """Write a small checkpoint laid out like extended_beta1000yrs.h5:
    timestepped velocity/thickness/surface at every index in ``history`` and
    static bed/grounded/floating/fluidity/extended_beta. Each index carries a
    distinct offset so reading the wrong index is detectable."""
    mesh = firedrake.UnitSquareMesh(3, 3)
    Q = firedrake.FunctionSpace(mesh, "CG", 1)
    V = firedrake.VectorFunctionSpace(mesh, "CG", 1)
    x, y = firedrake.SpatialCoordinate(mesh)

    with firedrake.CheckpointFile(str(path), "w") as out:
        out.save_mesh(mesh)
        velocity = firedrake.Function(V, name="velocity")
        thickness = firedrake.Function(Q, name="thickness")
        surface = firedrake.Function(Q, name="surface")
        for idx in history:
            velocity.interpolate(firedrake.as_vector((100.0 + idx + x, 3.0 * idx + y)))
            thickness.interpolate(500.0 + idx + 10.0 * x)
            surface.interpolate(50.0 + 0.5 * idx + y)
            out.save_function(velocity, idx=idx)
            out.save_function(thickness, idx=idx)
            out.save_function(surface, idx=idx)
        static = {
            "bed": -400.0 + 100.0 * x,
            "grounded": firedrake.conditional(x < 0.5, 1.0, 0.0),
            "floating": firedrake.conditional(x < 0.5, 0.0, 1.0),
            "fluidity": 20.0 + y,
            "extended_beta": 1000.0 + 5.0 * x * y,
        }
        for name, expr in static.items():
            out.save_function(firedrake.Function(Q, name=name).interpolate(expr))


class _RecordingSolver:
    """Stands in for icepack's FlowSolver: records exactly which fields
    initializeRun read from the checkpoint and returns the input velocity."""

    def __init__(self):
        self.calls = []

    def diagnostic_solve(self, **fields):
        self.calls.append(fields)
        return fields["velocity"].copy(deepcopy=True)


def _run_initialize(init_file, compact, monkeypatch):
    monkeypatch.setattr(model, "readSMB", lambda icesee_kwargs, Q: None)
    with firedrake.CheckpointFile(str(init_file), "r") as checkpoint:
        mesh = checkpoint.load_mesh()
    Q = firedrake.FunctionSpace(mesh, "CG", 1)
    V = firedrake.VectorFunctionSpace(mesh, "CG", 1)
    solver = _RecordingSolver()
    kwargs = {
        "initFile": str(init_file),
        "compact_initialization": compact,
        "uThresh": 300.0,
        "comm": mesh.comm,
    }
    result = model.initializeRun(kwargs, solver, mesh, Q, V)
    assert len(solver.calls) == 1
    return solver.calls[0], result


def _arrays(fields):
    return {name: np.array(f.dat.data_ro) for name, f in fields.items()
            if hasattr(f, "dat")}


def _expected_final_state(path):
    with firedrake.CheckpointFile(str(path), "r") as checkpoint:
        mesh = checkpoint.load_mesh()
        return {
            "velocity": np.array(checkpoint.load_function(mesh, "velocity", idx=_FINAL_IDX).dat.data_ro),
            "thickness": np.array(checkpoint.load_function(mesh, "thickness", idx=_FINAL_IDX).dat.data_ro),
            "surface": np.array(checkpoint.load_function(mesh, "surface", idx=_FINAL_IDX).dat.data_ro),
            "thickness_idx0": np.array(checkpoint.load_function(mesh, "thickness", idx=0).dat.data_ro),
        }


def test_compact_initialization_matches_full_history_final_state(tmp_path, monkeypatch):
    full = tmp_path / "full_history.h5"
    compact = tmp_path / "compact.h5"
    _write_full_history(full, _HISTORY)
    _load_build_module().build_compact_initialization(str(full), str(compact), _FINAL_IDX)

    expected = _expected_final_state(full)
    assert not np.array_equal(expected["thickness"], expected["thickness_idx0"])

    full_inputs, full_result = _run_initialize(full, compact=False, monkeypatch=monkeypatch)
    compact_inputs, compact_result = _run_initialize(compact, compact=True, monkeypatch=monkeypatch)

    full_arrays = _arrays(full_inputs)
    compact_arrays = _arrays(compact_inputs)
    checkpoint_fields = (
        "velocity", "thickness", "surface", "beta", "fluidity", "grounded", "floating",
    )
    assert set(full_arrays) == set(compact_arrays)
    for name in checkpoint_fields + ("uThresh",):
        assert np.array_equal(full_arrays[name], compact_arrays[name]), name

    # Both paths must have selected the final spin-up state, not any earlier index.
    assert np.array_equal(compact_arrays["velocity"], expected["velocity"])
    assert np.array_equal(compact_arrays["thickness"], expected["thickness"])
    assert np.array_equal(compact_arrays["surface"], expected["surface"])

    # Every initialized field returned by initializeRun (h, h0, s, s0, u, bed,
    # zF, grounded, floating, A0, beta0; smb is stubbed) is identical.
    assert len(full_result) == len(compact_result) == 12
    for position, (a, b) in enumerate(zip(full_result, compact_result)):
        if a is None and b is None:
            continue
        assert np.array_equal(np.array(a.dat.data_ro), np.array(b.dat.data_ro)), position


def test_compact_path_does_not_need_full_history_file(tmp_path, monkeypatch):
    full = tmp_path / "full_history.h5"
    compact = tmp_path / "compact.h5"
    _write_full_history(full, _HISTORY)
    _load_build_module().build_compact_initialization(str(full), str(compact), _FINAL_IDX)
    expected = _expected_final_state(full)
    os.remove(full)

    inputs, _ = _run_initialize(compact, compact=True, monkeypatch=monkeypatch)
    assert np.array_equal(np.array(inputs["thickness"].dat.data_ro), expected["thickness"])


def test_compact_file_size_is_independent_of_history_length(tmp_path):
    build = _load_build_module().build_compact_initialization
    short_full = tmp_path / "short_full.h5"
    long_full = tmp_path / "long_full.h5"
    _write_full_history(short_full, (_FINAL_IDX,))
    _write_full_history(long_full, tuple(range(0, 10)) + (_FINAL_IDX,))

    short_compact = tmp_path / "short_compact.h5"
    long_compact = tmp_path / "long_compact.h5"
    build(str(short_full), str(short_compact), _FINAL_IDX)
    build(str(long_full), str(long_compact), _FINAL_IDX)

    assert os.path.getsize(long_compact) == os.path.getsize(short_compact)
    assert os.path.getsize(long_compact) < os.path.getsize(long_full) / 10
