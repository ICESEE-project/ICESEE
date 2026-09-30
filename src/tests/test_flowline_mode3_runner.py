# ==============================================================================
# @des: Lightweight, fully-mocked unit tests for flowline_1d's mode-3
#       DA-cycle runner (``mode3_runner.py``). Mirrors
#       ``test_icepack_mode3_runner.py``'s/``test_issm_mode3_runner.py``'s
#       import dance and orchestration-focused test intent. Unlike those two,
#       ``_build_static_config`` is NOT mocked here -- it is pure,
#       rank-independent float/int arithmetic plus one fast, deterministic
#       ``scipy.optimize.root``/JAX solve with no randomness, so exercising
#       it for real is both fast and a stronger check than faking it. This
#       validates the runner's own control flow (setup-phase call order,
#       checkpoint directory/timestep scheme, observation scheduling, final
#       timing report) rather than real long-horizon flowline physics.
# @date: 2026-08-26
# @author: Brian Kyanjo
# ==============================================================================
from __future__ import annotations

import importlib
import os
import sys
import tempfile
import types
from pathlib import Path

import h5py
import numpy as np
import pytest

# Same standalone-script import dance as test_icepack_mode3_runner.py /
# test_issm_mode3_runner.py: _flowline_model.py (imported transitively via
# mode3_runner.py) and config._utility_imports both assume they are imported
# from inside examples/flowline_1d with that directory on sys.path / as cwd.
_REPO_ROOT = Path(__file__).resolve().parents[2]
_FLOWLINE_1D_DIR = _REPO_ROOT / "applications" / "flowline_model" / "examples" / "flowline_1d"
for _extra_path in (str(_REPO_ROOT), str(_REPO_ROOT.parent), str(_FLOWLINE_1D_DIR)):
    if _extra_path not in sys.path:
        sys.path.insert(0, _extra_path)

# config._utility_imports deletes/recreates data_path at import time, so point
# it at a throwaway directory instead of flowline_1d's own _modelrun_datasets.
_SAFE_DATA_PATH = tempfile.mkdtemp(prefix="icesee_flowline_mode3_runner_test_")
_ARGV_BACKUP = sys.argv[:]
_CWD_BACKUP = os.getcwd()
sys.argv = [sys.argv[0], "--data_path", _SAFE_DATA_PATH]
os.chdir(_FLOWLINE_1D_DIR)
try:
    from applications.flowline_model.examples.flowline_1d import mode3_runner as runner
    # config._utility_imports executes its whole CLI/YAML pipeline as
    # module-level code and Python only runs that once per process --
    # a plain `import`/`from ... import` here would silently return
    # whatever app's config another test module (collected earlier in
    # the same pytest session, e.g. test_enkf_serial_process_noise.py
    # under Lorenz96) already cached, instead of flowline_1d's own.
    # Force a fresh run under this file's own argv/cwd via reload.
    import ICESEE.config._utility_imports as _utility_imports_module
    importlib.reload(_utility_imports_module)
    _real_icesee_kwargs = _utility_imports_module.icesee_kwargs
finally:
    sys.argv = _ARGV_BACKUP
    os.chdir(_CWD_BACKUP)

from mpi4py import MPI

# ``mode3_runner.py`` itself imports these with the ``ICESEE.``-prefixed
# style, so import them here the same way to see the same registry.
from ICESEE.src.parallelization.distributed_mode3_registry import get_execution_mode_3
from ICESEE.src.parallelization.distributed_native_cycle import NativeCycleResult


def test_flowline_registers_its_mode3_runner():
    registration = get_execution_mode_3("flowline")
    assert registration is not None
    assert registration.runner is runner.run_flowline_execution_mode_3


def _base_kwargs(tmp_path, **overrides):
    # Real params.yaml-derived physical/grid parameters (hscale, A, n, C,
    # N1, N2, sigGZ, ...) come from the module-level config; copy it fresh
    # so tests never mutate the shared global.
    kwargs = dict(_real_icesee_kwargs)
    kwargs.update(
        Nens=2,
        num_years=2.0,
        data_path=str(tmp_path),
        execution_mode=3,
        obs_max_time=1,
        execution_flag=0,
        verbose=False,
        enkf_observation_error_mode="legacy_prior_anomalies",
    )
    kwargs.update(overrides)
    return kwargs


def test_build_static_config_is_rank_independent_and_derives_nd(tmp_path):
    kwargs = _base_kwargs(tmp_path)
    runner._build_static_config(kwargs)

    assert kwargs["nd"] == 2 * kwargs["NX"] + 1
    assert kwargs["NX"] == kwargs["N1"] + kwargs["N2"]
    assert kwargs["nt"] == int(float(kwargs["num_years"]))
    assert kwargs["var_nd"] == {
        "h": kwargs["NX"],
        "u": kwargs["NX"],
        "xg": 1,
    }

    # Rank-independent (no MPI dependency at all): computing it again from
    # a fresh copy produces bit-identical derived values.
    kwargs2 = _base_kwargs(tmp_path)
    runner._build_static_config(kwargs2)
    assert kwargs2["nd"] == kwargs["nd"]
    assert np.array_equal(kwargs2["t"], kwargs["t"])


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


def test_run_flowline_execution_mode_3_orchestration(tmp_path, monkeypatch):
    """End-to-end orchestration check with every heavy dependency faked.

    Validates: the initial checkpoint plus one checkpoint per timestep are
    written to the documented directory scheme; the analysis cycle only
    receives a batch on the scheduled observation timestep and never
    requests inversion (flowline_1d has none); and the final timing report /
    ``save_all_data`` are each invoked exactly once (root-only, but real
    ``MPI.COMM_WORLD`` under plain pytest has size 1 so rank 0 always runs
    this test).
    """

    calls = {"true_wrong": 0, "synth_obs": 0}

    fake_topology = types.SimpleNamespace(
        spatial_ranks=1, spatial_comm=None, world_size=1, ensemble_groups=1
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

    monkeypatch.setattr(runner, "create_distributed_topology", lambda world, spatial_ranks: fake_topology)
    monkeypatch.setattr(runner, "initialize_native_member_pool", lambda adapter, topo, kw: _FakePool())
    monkeypatch.setattr(runner, "generate_true_wrong_state", fake_generate_true_wrong_state)
    monkeypatch.setattr(runner, "generate_synthetic_observations", fake_generate_synthetic_observations)
    monkeypatch.setattr(runner, "UtilsFunctions", _FakeUtils)
    monkeypatch.setattr(runner, "run_native_global_analysis_cycle", fake_run_native_global_analysis_cycle)
    monkeypatch.setattr(runner, "save_distributed_checkpoint", fake_save_distributed_checkpoint)
    monkeypatch.setattr(runner, "save_all_data", lambda *a, **kw: save_all_data_calls.append(kw))
    monkeypatch.setattr(runner, "emit_performance_report", lambda *a, **kw: timing_calls.append(kw))

    icesee_kwargs = _base_kwargs(tmp_path)

    result = runner.run_flowline_execution_mode_3(**icesee_kwargs)

    assert calls["true_wrong"] == 1
    assert calls["synth_obs"] == 1
    nd = result["nd"]
    assert nd == 2 * result["NX"] + 1

    nt = result["nt"]
    # One initial checkpoint + one per timestep.
    assert len(checkpoint_calls) == 1 + nt
    initial_root, initial_timestep, initial_ensemble = checkpoint_calls[0]
    assert initial_timestep == 0
    assert initial_ensemble == "initial-ensemble-snapshot"
    assert os.path.basename(initial_root) == "initial_condition"

    step_roots = [c[0] for c in checkpoint_calls[1:]]
    step_timesteps = [c[1] for c in checkpoint_calls[1:]]
    assert all(os.path.basename(r) == "steps" for r in step_roots)
    assert step_timesteps == list(range(nt))
    assert [c[2] for c in checkpoint_calls[1:]] == [f"analysis-{k}" for k in range(nt)]

    # Analysis batch only requested on the scheduled observation timestep
    # (k=1); flowline_1d has no inversion mechanism at all.
    assert cycle_calls == [(k, 1 if k == 1 else 0) for k in range(nt)]

    assert len(save_all_data_calls) == 1
    assert len(timing_calls) == 1
