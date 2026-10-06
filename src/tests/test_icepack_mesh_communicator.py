# ==============================================================================
# @des: Regression tests for the Icepack Mode-3 world_size>1 hang fix.
#
# Background: Mode 3 partitions ensemble groups across independent
# sub-communicators (topology.spatial_comm), one group per rank/rank-group.
# `_icepack_native.py::_build_shared_context` correctly set
# `mesh_kwargs["comm"] = topology.spatial_comm` before calling
# `initializeMesh`, but `initializeMesh` (applications/icepack_model/
# examples/idealized_pig/_icepack_model.py) never forwarded that comm down
# to `setupMesh`/`argusToFiredrakeMesh`
# (applications/icepack_model/examples/idealized_pig/modelfunc/), whose
# `firedrake.Mesh(gmshFile)` call therefore always used Firedrake's own
# default communicator (COMM_WORLD). Since each independent ensemble
# group's rank reaches that call at a different, unsynchronized wall-clock
# time (each doing independent per-member work in between), this implicit
# COMM_WORLD collective mismatched across ranks and hung forever (observed
# stuck inside MPI_Comm_dup, which Firedrake's mesh construction performs
# internally). Confirmed via direct instrumentation: one rank alone
# entered the implicit-COMM_WORLD firedrake.Mesh() call while its "peer"
# was elsewhere in independent work.
#
# Fix: thread an explicit `comm` parameter through
# `argusToFiredrakeMesh(meshFile, savegmsh=False, comm=None)` ->
# `setupMesh(meshFile, ..., comm=None)` -> `initializeMesh` (now passes
# `comm=icesee_kwargs.get("comm")`), coalescing to `firedrake.COMM_WORLD`
# only when no comm is supplied at all (preserves old behavior for any
# caller that doesn't pass one).
#
# These tests mock only `firedrake.Mesh` itself (the expensive,
# hang-prone call) so they stay fast and deterministic; everything before
# it -- parsing the real production Argus mesh file and writing the real
# gmsh file -- runs for real (~0.1s), so the tests exercise the actual
# comm-forwarding code path, not a stand-in for it.
# ==============================================================================
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_IDEALIZED_PIG_DIR = _REPO_ROOT / "applications" / "icepack_model" / "examples" / "idealized_pig"
for _extra_path in (str(_REPO_ROOT), str(_REPO_ROOT.parent), str(_IDEALIZED_PIG_DIR)):
    if _extra_path not in sys.path:
        sys.path.insert(0, _extra_path)

_MESH_FILE = str(_IDEALIZED_PIG_DIR / "data" / "PigFull2017GeomFull.exp")

_ARGV_BACKUP = sys.argv[:]
_CWD_BACKUP = os.getcwd()
sys.argv = [sys.argv[0], "--data_path", tempfile.mkdtemp(prefix="icesee_test_data_path_")]
os.chdir(_IDEALIZED_PIG_DIR)
try:
    from applications.icepack_model.examples.idealized_pig import _icepack_model as model
    import applications.icepack_model.examples.idealized_pig.modelfunc.argusToFiredrakeMesh as argusToFiredrakeMesh
    import applications.icepack_model.examples.idealized_pig.modelfunc.setupMesh as setupMesh
finally:
    sys.argv = _ARGV_BACKUP
    os.chdir(_CWD_BACKUP)

import firedrake


class _CommCaptured(Exception):
    def __init__(self, comm):
        self.comm = comm


def _fake_mesh_capturing_comm(*args, **kwargs):
    raise _CommCaptured(kwargs.get("comm"))


def test_argus_to_firedrake_mesh_forwards_explicit_comm(monkeypatch):
    """argusToFiredrakeMesh must pass its own `comm` argument straight
    through to firedrake.Mesh(), not silently rely on Mesh()'s own
    COMM_WORLD default -- the exact bug that hung Mode-3 for world_size>1."""

    sentinel_comm = object()
    monkeypatch.setattr(argusToFiredrakeMesh.firedrake, "Mesh", _fake_mesh_capturing_comm)

    with pytest.raises(_CommCaptured) as excinfo:
        argusToFiredrakeMesh.argusToFiredrakeMesh(_MESH_FILE, comm=sentinel_comm)

    assert excinfo.value.comm is sentinel_comm


def test_argus_to_firedrake_mesh_falls_back_to_comm_world_when_comm_absent(monkeypatch):
    """No comm supplied at all (comm=None, the default) must still reach
    firedrake.Mesh() with a real communicator (COMM_WORLD), matching
    Firedrake's own convention and preserving pre-fix behavior for any
    caller that never passes a comm."""

    monkeypatch.setattr(argusToFiredrakeMesh.firedrake, "Mesh", _fake_mesh_capturing_comm)

    with pytest.raises(_CommCaptured) as excinfo:
        argusToFiredrakeMesh.argusToFiredrakeMesh(_MESH_FILE)

    assert excinfo.value.comm is firedrake.COMM_WORLD


def test_setup_mesh_forwards_comm_to_argus_to_firedrake_mesh(monkeypatch):
    """setupMesh must forward its own `comm` argument down to
    argusToFiredrakeMesh -- previously silently dropped here, defaulting
    the call three levels down to COMM_WORLD regardless of the caller's
    actual intended scope (e.g. topology.spatial_comm)."""

    sentinel_comm = object()
    monkeypatch.setattr(argusToFiredrakeMesh.firedrake, "Mesh", _fake_mesh_capturing_comm)

    with pytest.raises(_CommCaptured) as excinfo:
        setupMesh.setupMesh(_MESH_FILE, comm=sentinel_comm)

    assert excinfo.value.comm is sentinel_comm


def test_initialize_mesh_forwards_icesee_kwargs_comm(monkeypatch):
    """_icepack_model.initializeMesh (the function _icepack_native.py's
    _build_shared_context actually calls for Mode 3) must forward
    icesee_kwargs["comm"] all the way down to firedrake.Mesh(). Also
    stubs getMeshFromCheckPoint so this test never touches a real
    checkpoint file."""

    sentinel_comm = object()
    monkeypatch.setattr(argusToFiredrakeMesh.firedrake, "Mesh", _fake_mesh_capturing_comm)
    monkeypatch.setattr(model.mf, "getMeshFromCheckPoint", lambda *a, **k: None)

    icesee_kwargs = {
        "initFile": "unused.h5",
        "meshFile": _MESH_FILE,
        "comm": sentinel_comm,
    }

    with pytest.raises(_CommCaptured) as excinfo:
        model.initializeMesh(**icesee_kwargs)

    assert excinfo.value.comm is sentinel_comm
