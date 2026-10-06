# ==============================================================================
# @des: Lightweight, fully-mocked unit tests for idealized_pig's mode-3
#       DA-cycle runner (``mode3_runner.py``). Mirrors
#       ``test_icepack_idealized_pig_native.py``'s import dance and
#       ``test_lorenz_mode3_runner.py``'s orchestration-focused test intent,
#       but stays a pure unit test (no real Firedrake mesh, no real ~37GB
#       idealized_pig data) by monkeypatching every heavy dependency
#       ``run_icepack_execution_mode_3`` calls. This validates the runner's
#       own control flow (setup-phase call order, nd broadcast, checkpoint
#       directory/timestep scheme, observation scheduling, final timing
#       report) rather than the real Firedrake/icepack physics, which is
#       covered separately by ``test_icepack_idealized_pig_native.py`` and
#       real end-to-end smoke runs.
# @date: 2026-08-25
# @author: Brian Kyanjo
# ==============================================================================
from __future__ import annotations

import os
import sys
import tempfile
import types
from pathlib import Path

import h5py
import numpy as np
import pytest

# Same standalone-script import dance as test_icepack_idealized_pig_native.py
# (including the --data_path override that keeps config/_utility_imports.py's
# auto-clean guard off idealized_pig's own real _modelrun_datasets -- see
# that file's comment for why this is required, not just style):
# _icepack_model.py (imported transitively via mode3_runner.py) assumes it is
# imported from inside examples/idealized_pig with that directory on
# sys.path, and config._utility_imports parses sys.argv / loads params.yaml
# relative to the cwd / deletes+recreates data_path, all at import time.
_REPO_ROOT = Path(__file__).resolve().parents[2]
_IDEALIZED_PIG_DIR = _REPO_ROOT / "applications" / "icepack_model" / "examples" / "idealized_pig"
for _extra_path in (str(_REPO_ROOT), str(_REPO_ROOT.parent), str(_IDEALIZED_PIG_DIR)):
    if _extra_path not in sys.path:
        sys.path.insert(0, _extra_path)

_SAFE_DATA_PATH = tempfile.mkdtemp(prefix="icesee_icepack_mode3_runner_test_")
_ARGV_BACKUP = sys.argv[:]
_CWD_BACKUP = os.getcwd()
sys.argv = [sys.argv[0], "--data_path", _SAFE_DATA_PATH]
os.chdir(_IDEALIZED_PIG_DIR)
try:
    from applications.icepack_model.examples.idealized_pig import mode3_runner as runner
finally:
    sys.argv = _ARGV_BACKUP
    os.chdir(_CWD_BACKUP)

from mpi4py import MPI

# ``mode3_runner.py`` itself imports these two modules with the
# ``ICESEE.``-prefixed style (unlike some sibling test files), so importing
# them here without the prefix would bind a second, distinct module object
# with its own empty ``_REGISTRY`` -- match the runner's own style exactly so
# ``get_execution_mode_3`` sees the same registry it registered into.
from ICESEE.src.parallelization.distributed_mode3_registry import get_execution_mode_3
from ICESEE.src.parallelization.distributed_native_cycle import NativeCycleResult


def test_icepack_registers_its_mode3_runner():
    registration = get_execution_mode_3("icepack")
    assert registration is not None
    assert registration.runner is runner.run_icepack_execution_mode_3


class _FakeDat:
    def __init__(self, size):
        self.data = np.zeros(size)


class _FakeH0:
    def __init__(self, size):
        self.dat = _FakeDat(size)


def test_build_reference_setup_context_computes_global_nd_on_comm_self(monkeypatch):
    captured_comm = {}

    def fake_initialize_model(**kwargs):
        captured_comm["comm"] = kwargs["comm"]
        h0 = _FakeH0(4)
        return (
            "h", h0, "s", "s0", "u", "bed", "zF", "grounded", "floating",
            "A0", "beta0", "smb", "basal_melt_field", "Q", "V", "solver",
        )

    monkeypatch.setattr(runner, "initialize_model", fake_initialize_model)
    icesee_kwargs = {"total_state_param_vars": 3}

    nd = runner._build_reference_setup_context(icesee_kwargs)

    assert nd == 12
    assert icesee_kwargs["nd"] == 12
    assert icesee_kwargs["h0"].dat.data.size == 4
    assert icesee_kwargs["solver"] == "solver"
    assert captured_comm["comm"] is MPI.COMM_SELF


def test_run_icepack_execution_mode_3_orchestration(tmp_path, monkeypatch):
    """End-to-end orchestration check with every heavy dependency faked.

    Validates: setup phase runs once and broadcasts a correct ``nd``; the
    initial checkpoint plus one checkpoint per timestep are written to the
    documented directory scheme; the analysis cycle is only asked for a
    batch on the scheduled observation timestep; and the final timing
    report / ``save_all_data`` are each invoked exactly once (root-only,
    but real ``MPI.COMM_WORLD`` under plain pytest has size 1 so rank 0
    always runs this test).
    """

    calls = {"true_wrong": 0, "synth_obs": 0}

    def fake_initialize_model(**kwargs):
        h0 = _FakeH0(4)
        return (
            "h", h0, "s", "s0", "u", "bed", "zF", "grounded", "floating",
            "A0", "beta0", "smb", "basal_melt_field", "Q", "V", "solver",
        )

    def fake_generate_true_wrong_state(**kwargs):
        calls["true_wrong"] += 1
        return kwargs

    def fake_generate_synthetic_observations(**kwargs):
        calls["synth_obs"] += 1
        nd = kwargs["nd"]
        with h5py.File(kwargs["synthetic_obs_file"], "w") as f:
            f.create_dataset("hu_obs", data=np.ones((nd, 1)))
            f.create_dataset("R", data=np.eye(nd))
        return kwargs

    class _FakeModelModule:
        pass

    class _FakeSupportedModels:
        def __init__(self, model=None, verbose=None):
            self.model = model

        def call_model(self):
            return _FakeModelModule()

    fake_topology = types.SimpleNamespace(
        spatial_ranks=1, spatial_comm=None, world_size=1, ensemble_groups=1
    )

    class _FakePool:
        member_ids = [0, 1]

        def snapshot_owned(self):
            return "initial-ensemble-snapshot"

        def route_observation_ids(self, global_observation_ids, adapter, *, topology, icesee_kwargs):
            # Single spatial rank (fake_topology.spatial_ranks=1) owns
            # every requested observation -- identity routing.
            ids = np.asarray(global_observation_ids)
            return np.arange(ids.size, dtype=np.int64), ids

    class _FakeUtils:
        def __init__(self, icesee_kwargs):
            self.icesee_kwargs = icesee_kwargs

        def JObs_indices(self, nd):
            return np.arange(nd)

        def generate_observation_schedule(self, **kw):
            # One observation, scheduled at timestep k=1 (nt=2 -> k in {0,1}).
            return (np.array([1.0]), np.array([1]), 1)

    cycle_calls = []

    def fake_run_native_global_analysis_cycle(
        pool, adapter, k, batches, *, number_of_batches, topology, icesee_kwargs, error_mode
    ):
        cycle_calls.append((k, number_of_batches))
        return NativeCycleResult(
            local_forecast="forecast",
            local_analysis=f"analysis-{k}",
            transform=np.eye(2),
            observation_rows=len(batches[0].observation_ids) if batches else 0,
        )

    checkpoint_calls = []

    def fake_save_distributed_checkpoint(root, timestep, ensemble, topology, *, run_id, metadata=None):
        checkpoint_calls.append((root, timestep, ensemble))

    save_all_data_calls = []
    timing_calls = []

    monkeypatch.setattr(runner, "initialize_model", fake_initialize_model)
    monkeypatch.setattr(runner, "generate_true_wrong_state", fake_generate_true_wrong_state)
    monkeypatch.setattr(runner, "generate_synthetic_observations", fake_generate_synthetic_observations)
    monkeypatch.setattr(runner, "SupportedModels", _FakeSupportedModels)
    monkeypatch.setattr(runner, "create_distributed_topology", lambda world, spatial_ranks: fake_topology)
    monkeypatch.setattr(runner, "initialize_native_member_pool", lambda adapter, topo, kw: _FakePool())
    monkeypatch.setattr(runner, "UtilsFunctions", _FakeUtils)
    monkeypatch.setattr(runner, "run_native_global_analysis_cycle", fake_run_native_global_analysis_cycle)
    monkeypatch.setattr(runner, "save_distributed_checkpoint", fake_save_distributed_checkpoint)
    monkeypatch.setattr(runner, "save_all_data", lambda *a, **kw: save_all_data_calls.append(kw))
    monkeypatch.setattr(runner, "emit_performance_report", lambda *a, **kw: timing_calls.append(kw))

    icesee_kwargs = dict(
        model_name="icepack",
        Nens=2,
        nt=2,
        total_state_param_vars=1,
        data_path=str(tmp_path),
        execution_mode=3,
        obs_max_time=1,
        execution_flag=0,
        t=np.linspace(0, 1, 3),
        verbose=False,
    )

    result = runner.run_icepack_execution_mode_3(**icesee_kwargs)

    assert calls["true_wrong"] == 1
    assert calls["synth_obs"] == 1
    assert result["nd"] == 4

    # One initial checkpoint + one per timestep (nt=2).
    assert len(checkpoint_calls) == 1 + 2
    initial_root, initial_timestep, initial_ensemble = checkpoint_calls[0]
    assert initial_timestep == 0
    assert initial_ensemble == "initial-ensemble-snapshot"
    assert os.path.basename(initial_root) == "initial_condition"

    step_roots = [c[0] for c in checkpoint_calls[1:]]
    step_timesteps = [c[1] for c in checkpoint_calls[1:]]
    assert all(os.path.basename(r) == "steps" for r in step_roots)
    assert step_timesteps == [0, 1]
    assert [c[2] for c in checkpoint_calls[1:]] == ["analysis-0", "analysis-1"]

    # Analysis batch only requested on the scheduled observation timestep (k=1).
    assert cycle_calls == [(0, 0), (1, 1)]

    assert len(save_all_data_calls) == 1
    # Exactly one performance report per run. The fused native cycle has no
    # separate analysis timer, so the analysis step is reported as not
    # measured (None) rather than as a fake 0 when an analysis ran.
    assert len(timing_calls) == 1
    assert timing_calls[0]["phases"]["analysis_step"] is None
    assert timing_calls[0]["counts"]["analysis_step"] == 1
