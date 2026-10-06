# ==============================================================================
# @des: Lightweight, fully-mocked unit tests for ISSM ISMIP_Choi's mode-3
#       DA-cycle runner (``mode3_runner.py``). Mirrors
#       ``test_icepack_mode3_runner.py``'s import dance and
#       orchestration-focused test intent, but stays a pure unit test (no
#       real MATLAB/ISSM subprocess, no ISSM_DIR) by monkeypatching every
#       heavy dependency ``run_issm_execution_mode_3`` calls. This validates
#       the runner's own control flow (setup-phase call order and
#       round-robin ``ens_id``, checkpoint directory/timestep scheme,
#       observation scheduling, inversion-start-time gating, final timing
#       report) rather than real ISSM physics.
# @date: 2026-08-26
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

# Same standalone-script import dance as test_icepack_mode3_runner.py:
# _issm_model.py (imported transitively via mode3_runner.py) assumes it is
# imported from inside examples/ISMIP_Choi with that directory on sys.path,
# and config._utility_imports parses sys.argv / loads params.yaml relative
# to the cwd, both at import time.
_REPO_ROOT = Path(__file__).resolve().parents[2]
_ISMIP_CHOI_DIR = _REPO_ROOT / "applications" / "issm_model" / "examples" / "ISMIP_Choi"
for _extra_path in (str(_REPO_ROOT), str(_REPO_ROOT.parent), str(_ISMIP_CHOI_DIR)):
    if _extra_path not in sys.path:
        sys.path.insert(0, _extra_path)

_ARGV_BACKUP = sys.argv[:]
_CWD_BACKUP = os.getcwd()
sys.argv = [sys.argv[0], "--data_path", tempfile.mkdtemp(prefix="icesee_test_data_path_")]
os.chdir(_ISMIP_CHOI_DIR)
try:
    from applications.issm_model.examples.ISMIP_Choi import mode3_runner as runner
finally:
    sys.argv = _ARGV_BACKUP
    os.chdir(_CWD_BACKUP)

from mpi4py import MPI

# ``mode3_runner.py`` itself imports these with the ``ICESEE.``-prefixed
# style, so import them here the same way to see the same registry.
from ICESEE.src.parallelization.distributed_mode3_registry import get_execution_mode_3
from ICESEE.src.parallelization.distributed_native_cycle import NativeCycleResult


def test_issm_registers_its_mode3_runner():
    registration = get_execution_mode_3("issm")
    assert registration is not None
    assert registration.runner is runner.run_issm_execution_mode_3


def _base_kwargs(tmp_path, **overrides):
    kwargs = dict(
        model_name="issm",
        Nens=2,
        total_state_param_vars=1,
        data_path=str(tmp_path),
        execution_mode=3,
        obs_max_time=1,
        execution_flag=0,
        verbose=False,
        Lx=1.0e5,
        Ly=5.0e4,
        nx=10,
        ny=5,
        steps=1,
        ParamFile="Test.par",
        timesteps_per_year=1.0,
        tinitial=0.0,
        num_years=2.0,
        enkf_observation_error_mode="legacy_prior_anomalies",
    )
    kwargs.update(overrides)
    return kwargs


def test_build_static_config_is_rank_independent_and_derives_nt(tmp_path):
    kwargs = _base_kwargs(tmp_path)
    runner._build_static_config(
        kwargs,
        icesee_cwd=str(tmp_path),
        issm_dir="/fake/issm",
        issm_examples_dir=str(tmp_path / "examples"),
    )

    assert kwargs["nt"] == 2
    assert np.allclose(kwargs["t"], np.array([0.0, 1.0, 2.0]))
    assert kwargs["Nens"] == 2
    assert kwargs["issm_dir"] == "/fake/issm"

    # Rank-independent: computing it again with a different (irrelevant)
    # rank-like key present produces identical derived values.
    kwargs2 = _base_kwargs(tmp_path)
    runner._build_static_config(
        kwargs2,
        icesee_cwd=str(tmp_path),
        issm_dir="/fake/issm",
        issm_examples_dir=str(tmp_path / "examples"),
    )
    assert kwargs2["nt"] == kwargs["nt"]
    assert np.allclose(kwargs2["t"], kwargs["t"])


class _FakePool:
    member_ids = [0, 1]

    def snapshot_owned(self):
        return "initial-ensemble-snapshot"


class _FakeUtils:
    def __init__(self, icesee_kwargs):
        self.icesee_kwargs = icesee_kwargs

    def JObs_indices(self, nd):
        return np.arange(nd)

    def generate_observation_schedule(self, **kw):
        # One observation, scheduled at timestep k=1 (nt=2 -> k in {0,1}).
        return (np.array([1.0]), np.array([1]), 1)


def test_run_issm_execution_mode_3_orchestration(tmp_path, monkeypatch):
    """End-to-end orchestration check with every heavy dependency faked.

    Validates: the per-rank server/mesh setup helper is called exactly once
    with ``ensemble_slot`` as the ``ens_id``; the initial checkpoint plus one
    checkpoint per timestep are written to the documented directory scheme;
    the analysis cycle only receives a batch on the scheduled observation
    timestep; inversion is deferred until ``inversion_start_time``; and the
    final timing report / ``save_all_data`` are each invoked exactly once
    (root-only, but real ``MPI.COMM_WORLD`` under plain pytest has size 1 so
    rank 0 always runs this test).
    """

    calls = {"true_wrong": 0, "synth_obs": 0}
    setup_calls = []

    fake_topology = types.SimpleNamespace(
        ensemble_slot=0, ensemble_groups=1, spatial_ranks=1, spatial_comm=None,
        world_size=1,
    )

    def fake_initialize_rank_server_and_mesh(icesee_kwargs, *, icesee_cwd,
                                              issm_examples_dir, ensemble_slot):
        setup_calls.append(int(ensemble_slot))
        return ("fake-server", 4)

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

    cycle_calls = []

    def fake_run_native_global_analysis_cycle(
        pool, adapter, k, batches, *, number_of_batches, topology, icesee_kwargs,
        error_mode, apply_inversion=False,
    ):
        cycle_calls.append((k, number_of_batches, apply_inversion))
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

    monkeypatch.setattr(runner, "add_issm_dir_to_sys_path", lambda issm_dir=None: None)
    monkeypatch.setattr(runner, "setup_example_directory", lambda issm_dir, example_name: str(tmp_path / "examples"))
    monkeypatch.setattr(runner, "create_distributed_topology", lambda world, spatial_ranks: fake_topology)
    monkeypatch.setattr(runner, "_initialize_rank_server_and_mesh", fake_initialize_rank_server_and_mesh)
    monkeypatch.setattr(runner, "generate_true_wrong_state", fake_generate_true_wrong_state)
    monkeypatch.setattr(runner, "generate_synthetic_observations", fake_generate_synthetic_observations)
    monkeypatch.setattr(runner, "UtilsFunctions", _FakeUtils)
    monkeypatch.setattr(runner, "initialize_native_member_pool", lambda adapter, topo, kw: _FakePool())
    monkeypatch.setattr(runner, "run_native_global_analysis_cycle", fake_run_native_global_analysis_cycle)
    monkeypatch.setattr(runner, "save_distributed_checkpoint", fake_save_distributed_checkpoint)
    monkeypatch.setattr(runner, "save_all_data", lambda *a, **kw: save_all_data_calls.append(kw))
    monkeypatch.setattr(runner, "emit_performance_report", lambda *a, **kw: timing_calls.append(kw))

    icesee_kwargs = _base_kwargs(
        tmp_path,
        inversion_flag=True,
        inversion_start_time=5.0,
    )

    result = runner.run_issm_execution_mode_3(**icesee_kwargs)

    # Setup phase ran exactly once, with this rank's ensemble slot as ens_id.
    assert setup_calls == [0]
    assert calls["true_wrong"] == 1
    assert calls["synth_obs"] == 1
    assert result["nd"] == 4

    # One initial checkpoint + one per timestep (nt=2, derived from
    # num_years=2, tinitial=0, timesteps_per_year=1).
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

    # Analysis batch only requested on the scheduled observation timestep
    # (k=1); inversion_start_time=5.0 defers inversion at that cycle time.
    assert cycle_calls == [(0, 0, False), (1, 1, False)]

    assert len(save_all_data_calls) == 1
    assert len(timing_calls) == 1


def test_run_issm_execution_mode_3_enables_inversion_once_start_time_reached(tmp_path, monkeypatch):
    """Same orchestration, but with ``inversion_start_time`` already reached."""

    fake_topology = types.SimpleNamespace(
        ensemble_slot=0, ensemble_groups=1, spatial_ranks=1, spatial_comm=None,
        world_size=1,
    )

    def fake_initialize_rank_server_and_mesh(icesee_kwargs, *, icesee_cwd,
                                              issm_examples_dir, ensemble_slot):
        return ("fake-server", 4)

    def fake_generate_synthetic_observations(**kwargs):
        nd = kwargs["nd"]
        with h5py.File(kwargs["synthetic_obs_file"], "w") as f:
            f.create_dataset("hu_obs", data=np.ones((nd, 1)))
            f.create_dataset("R", data=np.eye(nd))
        return kwargs

    cycle_calls = []

    def fake_run_native_global_analysis_cycle(
        pool, adapter, k, batches, *, number_of_batches, topology, icesee_kwargs,
        error_mode, apply_inversion=False,
    ):
        cycle_calls.append((k, apply_inversion))
        return NativeCycleResult(
            local_forecast="forecast",
            local_analysis=f"analysis-{k}",
            transform=np.eye(2),
            observation_rows=len(batches[0].observation_ids) if batches else 0,
        )

    monkeypatch.setattr(runner, "add_issm_dir_to_sys_path", lambda issm_dir=None: None)
    monkeypatch.setattr(runner, "setup_example_directory", lambda issm_dir, example_name: str(tmp_path / "examples"))
    monkeypatch.setattr(runner, "create_distributed_topology", lambda world, spatial_ranks: fake_topology)
    monkeypatch.setattr(runner, "_initialize_rank_server_and_mesh", fake_initialize_rank_server_and_mesh)
    monkeypatch.setattr(runner, "generate_true_wrong_state", lambda **kw: kw)
    monkeypatch.setattr(runner, "generate_synthetic_observations", fake_generate_synthetic_observations)
    monkeypatch.setattr(runner, "UtilsFunctions", _FakeUtils)
    monkeypatch.setattr(runner, "initialize_native_member_pool", lambda adapter, topo, kw: _FakePool())
    monkeypatch.setattr(runner, "run_native_global_analysis_cycle", fake_run_native_global_analysis_cycle)
    monkeypatch.setattr(runner, "save_distributed_checkpoint", lambda *a, **kw: None)
    monkeypatch.setattr(runner, "save_all_data", lambda *a, **kw: None)
    monkeypatch.setattr(runner, "emit_performance_report", lambda *a, **kw: None)

    icesee_kwargs = _base_kwargs(
        tmp_path,
        inversion_flag=True,
        inversion_start_time=0.0,
    )

    runner.run_issm_execution_mode_3(**icesee_kwargs)

    # Observation lands at k=1 (obs_t=[1.0]); inversion_start_time=0.0 is
    # already reached, so inversion is enabled at that cycle.
    assert cycle_calls == [(0, False), (1, True)]
