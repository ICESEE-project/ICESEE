# ==============================================================================
# @des: Regression tests for
# applications/icepack_model/examples/idealized_pig/tools/
# build_compact_initialization.py and _icepack_model.py::initializeRun's
# corresponding compact_initialization read-side kwarg.
#
# Mode-3 scalability investigation finding this addresses: the production
# ~39GB initFile stores a thickness/velocity/surface TIME SERIES (idx
# 0..20000), but initializeRun only ever reads exactly one index
# (idx=20000). CheckpointFile.save_function(f, idx=N) measured to
# allocate on-disk storage proportional to N even when only that one
# index is ever written (34.67 GB for one field written at idx=20000, vs
# 22.2 MB for the same field in normal/non-timestepping mode) -- so the
# compact file must be written WITHOUT idx=, and initializeRun needs an
# explicit, default-False opt-in (compact_initialization) to read it that
# way. Numerical equivalence against the full file was verified directly
# (all 12 returned fields from initializeRun, including the derived
# diagnostic-solve velocity and SMB: max_abs_diff == 0.0 exactly).
#
# The real ~39GB source file only exists on machines with the production
# idealized_pig dataset -- every test here that needs it skips cleanly
# when it (or Firedrake) is unavailable, exactly like this repo's other
# real-data-gated tests (see test_icepack_idealized_pig_native.py's own
# pattern).
# ==============================================================================
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_IDEALIZED_PIG_DIR = _REPO_ROOT / "applications" / "icepack_model" / "examples" / "idealized_pig"
_TOOLS_DIR = _IDEALIZED_PIG_DIR / "tools"
_SOURCE_H5 = _IDEALIZED_PIG_DIR / "data" / "extended_beta1000yrs.h5"


def _firedrake_available():
    try:
        import firedrake  # noqa: F401
    except ImportError:
        return False
    return True


def _load_build_module():
    """Import build_compact_initialization.py without requiring it to sit
    on sys.path permanently (mirrors this repo's other standalone-script
    import patterns)."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "build_compact_initialization", _TOOLS_DIR / "build_compact_initialization.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_field_name_lists_match_initializeRun_exactly():
    # Pure-Python, no Firedrake needed: these two tuples are not free
    # choices -- they must exactly match the field names
    # _icepack_model.py::initializeRun actually reads, or the compact
    # file silently omits something production needs.
    module = _load_build_module()
    assert module.TIMESTEPPED_FIELDS == ("velocity", "thickness", "surface")
    assert module.STATIC_FIELDS == ("bed", "grounded", "floating", "fluidity", "extended_beta")

    init_run_source = (_IDEALIZED_PIG_DIR / "_icepack_model.py").read_text()
    for name in module.TIMESTEPPED_FIELDS + module.STATIC_FIELDS:
        assert f'"{name}"' in init_run_source, (
            f"field {name!r} not found in _icepack_model.py -- "
            "build_compact_initialization.py's field lists have drifted "
            "from what initializeRun actually reads"
        )


def test_rejects_source_equal_to_dest():
    module = _load_build_module()
    with pytest.raises(ValueError):
        module.build_compact_initialization("/tmp/same.h5", "/tmp/same.h5", 20000, comm=None)


def test_compact_initialization_kwarg_defaults_to_full_file_behavior_unchanged():
    # Static, no-Firedrake-needed guard: the default path (kwarg absent)
    # must still read velocity/thickness/surface at idx=20000 from the
    # full file -- i.e. the new kwarg must be opt-in, not a silent
    # behavior change for every existing idealized_pig configuration.
    source = (_IDEALIZED_PIG_DIR / "_icepack_model.py").read_text()
    assert 'icesee_kwargs.get("compact_initialization", False)' in source
    assert 'checkpoint.load_function(mesh, "velocity", idx=20000)' in source


@pytest.mark.skipif(not _firedrake_available(), reason="firedrake is not importable in this environment")
@pytest.mark.skipif(
    not _SOURCE_H5.exists(),
    reason=(
        "the real ~39GB idealized_pig production checkpoint "
        f"({_SOURCE_H5}) is not present in this environment"
    ),
)
def test_compact_file_is_dramatically_smaller_and_equivalent(tmp_path):
    """Real, end-to-end check (only runs where the production dataset
    exists): build a compact file from the real source, then verify (a)
    it is far smaller and (b) initializeRun returns numerically identical
    fields whether reading the full file at idx=20000 or the compact file
    via compact_initialization=True."""
    import numpy as np
    import firedrake
    import icepack

    module = _load_build_module()
    dest = tmp_path / "compact.h5"
    timings = module.build_compact_initialization(str(_SOURCE_H5), str(dest), 20000)

    source_bytes = os.path.getsize(_SOURCE_H5)
    dest_bytes = os.path.getsize(dest)
    assert dest_bytes < source_bytes / 100, (
        f"compact file ({dest_bytes} bytes) is not dramatically smaller "
        f"than the source ({source_bytes} bytes) -- regression guard for "
        "the idx=N-allocates-proportional-to-N finding"
    )

    sys.argv = [sys.argv[0], "--data_path", tempfile.mkdtemp(prefix="icesee_test_data_path_")]
    saved_cwd = os.getcwd()
    os.chdir(_IDEALIZED_PIG_DIR)
    # idealized_pig's own "modelfunc" package does a bare (non-relative)
    # internal import that only resolves when its own directory is on
    # sys.path (see icesee_da_distributed.py's _jit_import_mode3_runner
    # fix for the same root cause) -- required here too, or importing
    # _icepack_model.py below fails with ModuleNotFoundError: modelfunc.
    if str(_IDEALIZED_PIG_DIR) not in sys.path:
        sys.path.insert(0, str(_IDEALIZED_PIG_DIR))
    try:
        # Deliberately NOT "from ICESEE.config._utility_imports import
        # icesee_kwargs" here: that module builds a process-wide SINGLETON
        # dict once, the first time anything imports it in this pytest
        # session -- if an earlier-collected test file (any application,
        # not just icepack) imported it first, this test would silently
        # get THAT application's config instead of idealized_pig's own
        # (reproduced directly: KeyError: 'meshFile' when running the full
        # suite, because an unrelated app's params.yaml has no such key).
        # Loading idealized_pig's own params.yaml directly via
        # config_loader (the same two functions _utility_imports.py
        # itself uses to build icesee_kwargs, config/_utility_imports.py
        # lines 483-485) sidesteps the whole import-order dependency.
        sys.path.insert(0, str(_REPO_ROOT / "config"))
        from config_loader import load_yaml_to_dict, get_section

        _parameters = load_yaml_to_dict(str(_IDEALIZED_PIG_DIR / "params.yaml"))
        base_kwargs = {}
        base_kwargs.update(get_section(_parameters, "physical-parameters"))
        base_kwargs.update(get_section(_parameters, "modeling-parameters"))
        base_kwargs.update(get_section(_parameters, "enkf-parameters"))

        from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_model import (
            initializeMesh, initializeRun, schoofFriction, regViscosity,
        )

        def run(compact, init_file):
            kw = dict(base_kwargs)
            kw.update({
                "comm": firedrake.COMM_WORLD,
                "initFile": init_file,
                "compact_initialization": compact,
            })
            mesh, _meshOpts, Q, V = initializeMesh(**kw)
            forward_model = icepack.models.IceStream(friction=schoofFriction, viscosity=regViscosity)
            opts = {"dirichlet_ids": [1], "diagnostic_solver_parameters": {"max_iterations": 150, "tolerance": 1e-6}}
            forward_solver = icepack.solvers.FlowSolver(forward_model, **opts)
            return initializeRun(kw, forward_solver, mesh, Q, V)

        res_full = run(False, str(_SOURCE_H5))
        res_compact = run(True, str(dest))
    finally:
        os.chdir(saved_cwd)

    for a, b in zip(res_full, res_compact):
        if not hasattr(a, "dat"):
            continue
        av = a.dat.data_ro
        bv = b.dat.data_ro
        assert av.shape == bv.shape
        assert np.array_equal(av, bv), (
            "compact-file initialization diverged from full-file "
            "initialization -- max abs diff "
            f"{np.abs(av - bv).max()}"
        )
