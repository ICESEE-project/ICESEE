# ==============================================================================
# @des: Execution-mode-3 restart discovery, resume planning, bounded
#       checkpoint retention, and crash consistency
#       (src/parallelization/distributed_checkpoint.py,
#       src/parallelization/mode3_checkpointing.py,
#       restore_native_member_pool in distributed_native_runtime.py).
#       Single-process fakes; the real multi-rank path is covered by
#       test_mode3_checkpoint_recovery_mpi.py.
# ==============================================================================
from __future__ import annotations

import errno
import json
from pathlib import Path
import re
import types

import h5py
import numpy as np
import pytest

from ICESEE.src.parallelization import distributed_checkpoint as checkpoint_module
from ICESEE.src.parallelization.distributed_adapter import (
    DistributedBlockStateLayout,
    DistributedStateBlockLayout,
)
from ICESEE.src.parallelization.distributed_checkpoint import (
    CheckpointValidationError,
    discover_latest_checkpoint,
    load_distributed_checkpoint,
    remove_abandoned_checkpoint_artifacts,
    save_distributed_checkpoint,
    validate_distributed_checkpoint,
)
from ICESEE.src.parallelization.distributed_fields import DistributedFieldRegistry
from ICESEE.src.parallelization.distributed_native_runtime import (
    NativeDistributedMember,
    restore_native_member_pool,
)
from ICESEE.src.parallelization.distributed_runtime import LocalMemberEnsemble
from ICESEE.src.parallelization.mode3_checkpointing import (
    CheckpointPolicy,
    Mode3CheckpointManager,
    estimate_checkpoint_storage,
)

RUN_ID = "icepack-idealized-pig-mode3"


class _Comm:
    def Get_rank(self):
        return 0

    def Get_size(self):
        return 1

    def allgather(self, value):
        return [value]

    def gather(self, value, root=0):
        return [value]

    def bcast(self, value, root=0):
        return value

    def Barrier(self):
        return None


def _topology(ensemble_groups=1, ensemble_slot=0):
    comm = _Comm()
    return types.SimpleNamespace(
        world=comm,
        spatial_comm=comm,
        ensemble_comm=comm,
        world_rank=0,
        world_size=1,
        ensemble_groups=ensemble_groups,
        ensemble_slot=ensemble_slot,
        spatial_ranks=1,
        is_ensemble_root=True,
    )


def _layout():
    # Two opaque blocks; the generic code never interprets their names.
    return DistributedBlockStateLayout(
        (
            DistributedStateBlockLayout("a", 0, 6, 0, 6),
            DistributedStateBlockLayout("b", 6, 4, 0, 4),
        ),
        layout_id="opaque-state-v1",
    )


def _ensemble(step, members=3):
    return LocalMemberEnsemble(
        _layout(),
        {m: np.arange(10, dtype=float) + 100.0 * m + 0.5 * step for m in range(members)},
    )


def _write_legacy(root: Path, step: int, *, did_analysis: bool, members=3) -> Path:
    """A checkpoint exactly as the pre-retention runner wrote it.

    The writer of that version produced the same files; its manifest lacked
    only the per-shard ``nbytes`` and the restart metadata added since.
    """

    saved = save_distributed_checkpoint(
        root,
        step,
        _ensemble(step, members),
        _topology(),
        run_id=RUN_ID,
        metadata={
            "Nens": members,
            "observation_rows": 7 if did_analysis else 0,
            "error_mode": "legacy_prior_anomalies",
            "did_analysis": did_analysis,
        },
    )
    manifest_path = saved.path / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for shard in manifest["shards"]:
        shard.pop("nbytes", None)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return saved.path


def _legacy_history(tmp_path, *, upto, analysis_steps, members=3):
    steps = tmp_path / "_mode3_state_history" / "steps"
    for step in range(upto + 1):
        _write_legacy(steps, step, did_analysis=step in analysis_steps, members=members)
    return steps


def _manager(tmp_path, *, nt, analysis_steps, policy=None, members=3, t=None):
    return Mode3CheckpointManager(
        tmp_path / "_mode3_state_history",
        run_id=RUN_ID,
        topology=_topology(),
        number_of_members=members,
        nt=nt,
        analysis_steps=analysis_steps,
        policy=policy or CheckpointPolicy(),
        time_grid=t,
    )


def _committed_steps(root: Path):
    return sorted(
        int(p.name.split("_")[1]) for p in root.iterdir() if re.fullmatch(r"checkpoint_\d{8}", p.name)
    )


def _leave_interrupted_transaction(root: Path, step: int):
    """What a killed run (or ENOSPC on a non-root rank) leaves behind."""

    staging = root / f".checkpoint_{step:08d}.{RUN_ID}.staging"
    staging.mkdir(parents=True)
    (staging / "shard_00000000.h5").write_bytes(b"partial")
    (staging / ".shard_00000001.h5.tmp").write_bytes(b"\x89HDF")
    return staging


# --- 1-4: discovery -----------------------------------------------------------

def test_discovery_selects_highest_valid_committed_checkpoint(tmp_path):
    steps = _legacy_history(tmp_path, upto=4, analysis_steps=set())
    found = discover_latest_checkpoint(steps, expected_run_id=RUN_ID, number_of_members=3)
    assert found.timestep == 4
    assert found.path == steps / "checkpoint_00000004"
    assert found.rejected == ()


def test_discovery_ignores_incomplete_staging_transaction(tmp_path):
    # The student's failure: 525 committed, 526 died in staging (ENOSPC).
    steps = tmp_path / "steps"
    _write_legacy(steps, 525, did_analysis=False)
    _leave_interrupted_transaction(steps, 526)
    found = discover_latest_checkpoint(steps, expected_run_id=RUN_ID, number_of_members=3)
    assert found.timestep == 525
    assert found.rejected == ()  # staging is never even a candidate


def test_committed_directory_with_only_tmp_shard_is_rejected(tmp_path):
    steps = tmp_path / "steps"
    _write_legacy(steps, 3, did_analysis=False)
    bad = _write_legacy(steps, 4, did_analysis=False)
    (bad / "shard_00000000.h5").rename(bad / ".shard_00000000.h5.tmp")
    found = discover_latest_checkpoint(steps, expected_run_id=RUN_ID, number_of_members=3)
    assert found.timestep == 3
    assert found.rejected[0][0] == "checkpoint_00000004"
    assert "missing" in found.rejected[0][1]


@pytest.mark.parametrize(
    "damage",
    ["truncated_shard", "manifest_incomplete", "manifest_garbage", "missing_member", "no_manifest"],
)
def test_corrupt_latest_checkpoint_falls_back_to_previous_valid(tmp_path, damage):
    steps = tmp_path / "steps"
    _write_legacy(steps, 7, did_analysis=False)
    bad = _write_legacy(steps, 8, did_analysis=False)
    manifest_path = bad / "manifest.json"
    if damage == "truncated_shard":
        data = (bad / "shard_00000000.h5").read_bytes()
        (bad / "shard_00000000.h5").write_bytes(data[: len(data) // 2])
    elif damage == "manifest_incomplete":
        manifest = json.loads(manifest_path.read_text())
        manifest["complete"] = False
        manifest_path.write_text(json.dumps(manifest))
    elif damage == "manifest_garbage":
        manifest_path.write_text('{"format": ')
    elif damage == "missing_member":
        manifest = json.loads(manifest_path.read_text())
        manifest["shards"][0]["member_ids"] = [0, 1]
        manifest_path.write_text(json.dumps(manifest))
    elif damage == "no_manifest":
        manifest_path.unlink()
    found = discover_latest_checkpoint(steps, expected_run_id=RUN_ID, number_of_members=3)
    assert found.timestep == 7
    assert [name for name, _ in found.rejected] == ["checkpoint_00000008"]


def test_size_recorded_by_current_writer_detects_truncation_without_opening(tmp_path):
    saved = save_distributed_checkpoint(
        tmp_path, 2, _ensemble(2), _topology(), run_id=RUN_ID, metadata={"Nens": 3}
    )
    shard = saved.path / "shard_00000000.h5"
    shard.write_bytes(shard.read_bytes()[:-16])
    with pytest.raises(CheckpointValidationError, match="bytes"):
        validate_distributed_checkpoint(saved.path, deep=False)


def test_discovery_rejects_wrong_run_and_ensemble_size(tmp_path):
    steps = tmp_path / "steps"
    _write_legacy(steps, 1, did_analysis=False)
    assert discover_latest_checkpoint(steps, expected_run_id="other").path is None
    assert discover_latest_checkpoint(steps, number_of_members=30).path is None


# --- 5-8: resume semantics ----------------------------------------------------

def test_resume_continues_after_completed_step_with_correct_time_and_analyses(tmp_path):
    # Student-like schedule: nt=600, dt=0.05, analyses every 20 steps from 200.
    nt = 600
    t = np.linspace(0, 30, nt + 1)
    analysis = list(range(200, 600, 20))
    steps = tmp_path / "_mode3_state_history" / "steps"
    for step in (523, 524, 525):
        _write_legacy(steps, step, did_analysis=step in analysis)
    _leave_interrupted_transaction(steps, 526)

    manager = _manager(tmp_path, nt=nt, analysis_steps=analysis, t=t)
    plan = manager.prepare_resume()

    assert plan.resumed and plan.legacy_checkpoint
    assert plan.completed_step == 525
    assert plan.start_step == 526
    # Analyses at 200..520 are done; 540, 560, 580 remain.
    assert plan.analyses_completed == 17
    assert [s for s in analysis if s >= plan.start_step] == [540, 560, 580]
    # State after cycle k is the state at t[k + 1].
    assert manager.model_time_after(525) == pytest.approx(26.3)
    assert plan.removed_artifacts == (f".checkpoint_00000526.{RUN_ID}.staging",)
    assert not any(p.name.endswith(".staging") for p in steps.iterdir())


def test_resume_immediately_after_analysis_does_not_repeat_it(tmp_path):
    analysis = [3, 6, 9]
    _legacy_history(tmp_path, upto=6, analysis_steps=set(analysis))
    plan = _manager(tmp_path, nt=10, analysis_steps=analysis).prepare_resume()
    assert plan.completed_step == 6
    assert plan.start_step == 7
    # Analyses 3 and 6 are inside the restored state; only 9 remains.
    assert plan.analyses_completed == 2
    assert analysis[plan.analyses_completed:] == [9]


def test_resume_immediately_before_analysis_runs_it_once(tmp_path):
    analysis = [3, 6, 9]
    _legacy_history(tmp_path, upto=5, analysis_steps=set(analysis))
    plan = _manager(tmp_path, nt=10, analysis_steps=analysis).prepare_resume()
    # The analysis at step 6 has not happened yet; the resumed run's first
    # cycle performs it.
    assert (plan.completed_step, plan.start_step, plan.analyses_completed) == (5, 6, 1)
    assert analysis[plan.analyses_completed] == plan.start_step


def test_runner_loop_from_plan_replays_exactly_the_remaining_schedule(tmp_path):
    """The mode-3 loop condition, driven from a resume plan, visits each
    remaining analysis once and none of the completed ones."""

    nt, analysis = 12, [2, 5, 8, 11]
    _legacy_history(tmp_path, upto=5, analysis_steps=set(analysis))
    plan = _manager(tmp_path, nt=nt, analysis_steps=analysis).prepare_resume()
    km, done = plan.analyses_completed, []
    for k in range(plan.start_step, nt):
        if km < len(analysis) and k == analysis[km]:
            done.append(k)
            km += 1
    assert done == [8, 11]
    assert km == len(analysis)


def test_resume_rejects_configuration_that_changes_the_schedule(tmp_path):
    _legacy_history(tmp_path, upto=6, analysis_steps={3, 6})
    with pytest.raises(ValueError, match="did_analysis"):
        _manager(tmp_path, nt=10, analysis_steps=[4, 8]).prepare_resume()
    with pytest.raises(ValueError, match="members"):
        _manager(tmp_path, nt=10, analysis_steps=[3, 6], members=4).prepare_resume()


def test_resume_without_any_valid_checkpoint_refuses_to_start_over(tmp_path):
    steps = tmp_path / "_mode3_state_history" / "steps"
    _leave_interrupted_transaction(steps, 0)
    with pytest.raises(FileNotFoundError, match="no valid committed checkpoint"):
        _manager(tmp_path, nt=10, analysis_steps=[]).prepare_resume()


def test_resume_from_initial_condition_when_no_step_was_committed(tmp_path):
    manager = _manager(tmp_path, nt=10, analysis_steps=[4])
    manager.write_initial(_ensemble(-1))
    plan = manager.prepare_resume()
    assert (plan.completed_step, plan.start_step, plan.analyses_completed) == (-1, 0, 0)
    assert plan.checkpoint_path.parent.name == "initial_condition"


def test_invalid_newer_committed_directory_is_quarantined_not_deleted(tmp_path):
    steps = _legacy_history(tmp_path, upto=5, analysis_steps=set())
    bad = _write_legacy(steps, 6, did_analysis=False)
    (bad / "shard_00000000.h5").write_bytes(b"not hdf5")
    manager = _manager(tmp_path, nt=10, analysis_steps=[])
    plan = manager.prepare_resume()
    assert plan.completed_step == 5
    assert plan.quarantined == (".checkpoint_00000006.invalid",)
    assert (steps / ".checkpoint_00000006.invalid" / "shard_00000000.h5").exists()
    # Step 6 can now be written by the resumed run.
    assert manager.commit_step(6, _ensemble(6), did_analysis=False, analyses_completed=0)


def test_restore_reproduces_saved_state_without_initializing_members(tmp_path):
    class _Field:
        def __init__(self, name, offset, size):
            self.name, self.global_size = name, size
            self.owned_start, self.owned_stop = 0, size
            self.synchronization_key = (name, id(self))
            self._values = np.zeros(size)

        def read_owned(self):
            return self._values.copy()

        def write_owned(self, values):
            self._values[:] = values

        def synchronize_ghosts(self):
            return None

    class _Adapter:
        def __init__(self):
            self.allocated = []

        def initialize_native_member(self, member_id, **_):
            raise AssertionError("a resume must not regenerate the initial ensemble")

        def allocate_native_member(self, member_id, **_):
            self.allocated.append(member_id)
            registry = DistributedFieldRegistry(
                [_Field("a", 0, 6), _Field("b", 6, 4)], layout_id="opaque-state-v1"
            )
            return NativeDistributedMember(member_id, registry, None)

        forecast_native_member = observe_native_member = finalize_native_analysis = (
            lambda *a, **k: None
        )

    steps = tmp_path / "steps"
    saved = save_distributed_checkpoint(
        steps, 4, _ensemble(4), _topology(), run_id=RUN_ID, metadata={"Nens": 3}
    )
    adapter = _Adapter()
    pool, checkpoint = restore_native_member_pool(
        adapter, _topology(), {"Nens": 3}, saved.path, expected_run_id=RUN_ID
    )
    assert adapter.allocated == [0, 1, 2]
    assert checkpoint.timestep == 4
    for member_id, values in _ensemble(4).members.items():
        np.testing.assert_array_equal(pool.snapshot_owned().members[member_id], values)


# --- 9-11: retention and crash consistency ---------------------------------------

def _run(manager, steps_range, analysis):
    km = 0
    for k in steps_range:
        if k in analysis:
            km += 1
        manager.commit_step(k, _ensemble(k), did_analysis=k in analysis, analyses_completed=km)


def test_default_retention_keeps_only_the_newest_two(tmp_path):
    manager = _manager(tmp_path, nt=10, analysis_steps=[3, 6])
    _run(manager, range(10), {3, 6})
    assert _committed_steps(manager.steps_root) == [8, 9]
    assert not any(p.name.startswith(".") for p in manager.steps_root.iterdir())
    assert discover_latest_checkpoint(manager.steps_root).timestep == 9


def test_retention_can_pin_analysis_steps_and_keep_last_one(tmp_path):
    policy = CheckpointPolicy(keep_last=1, keep_analysis=True)
    manager = _manager(tmp_path, nt=10, analysis_steps=[3, 6], policy=policy)
    _run(manager, range(10), {3, 6})
    assert _committed_steps(manager.steps_root) == [3, 6, 9]


def test_keep_last_zero_keeps_full_history(tmp_path):
    manager = _manager(tmp_path, nt=6, analysis_steps=[], policy=CheckpointPolicy(keep_last=0))
    _run(manager, range(6), set())
    assert _committed_steps(manager.steps_root) == list(range(6))


def test_checkpoint_every_writes_interval_and_final_step(tmp_path):
    policy = CheckpointPolicy(every=4, keep_last=0)
    manager = _manager(tmp_path, nt=10, analysis_steps=[], policy=policy)
    _run(manager, range(10), set())
    assert _committed_steps(manager.steps_root) == [3, 7, 9]


def test_retention_never_deletes_pre_retention_checkpoints(tmp_path):
    # Resuming the student's run must not delete her existing history.
    steps = _legacy_history(tmp_path, upto=5, analysis_steps=set())
    manager = _manager(tmp_path, nt=12, analysis_steps=[])
    plan = manager.prepare_resume()
    _run(manager, range(plan.start_step, 12), set())
    assert _committed_steps(steps) == [0, 1, 2, 3, 4, 5, 10, 11]


def test_failed_checkpoint_write_keeps_previous_valid_checkpoints(tmp_path, monkeypatch):
    manager = _manager(tmp_path, nt=10, analysis_steps=[])
    _run(manager, range(5), set())
    assert _committed_steps(manager.steps_root) == [3, 4]

    real_file = h5py.File

    def full_disk(path, mode="r", *args, **kwargs):
        if mode == "w":
            raise OSError(errno.ENOSPC, "No space left on device")
        return real_file(path, mode, *args, **kwargs)

    monkeypatch.setattr(checkpoint_module.h5py, "File", full_disk)
    with pytest.raises(RuntimeError, match="No space left on device"):
        manager.commit_step(5, _ensemble(5), did_analysis=False, analyses_completed=0)
    monkeypatch.undo()

    assert _committed_steps(manager.steps_root) == [3, 4]
    found = discover_latest_checkpoint(manager.steps_root, expected_run_id=RUN_ID)
    assert found.timestep == 4
    loaded, _ = load_distributed_checkpoint(
        found.path, _layout(), _topology(), expected_run_id=RUN_ID, number_of_members=3
    )
    np.testing.assert_array_equal(loaded.members[1], _ensemble(4).members[1])


def test_killed_writer_and_killed_pruner_leave_previous_checkpoint_recoverable(tmp_path):
    manager = _manager(tmp_path, nt=10, analysis_steps=[])
    _run(manager, range(4), set())
    _leave_interrupted_transaction(manager.steps_root, 4)
    # A prune interrupted after its un-publishing rename.
    (manager.steps_root / "checkpoint_00000002").rename(
        manager.steps_root / ".checkpoint_00000002.deleting"
    )
    assert discover_latest_checkpoint(manager.steps_root).timestep == 3

    removed = remove_abandoned_checkpoint_artifacts(manager.steps_root)
    assert sorted(p.name for p in removed) == [
        ".checkpoint_00000002.deleting",
        f".checkpoint_00000004.{RUN_ID}.staging",
    ]
    assert _committed_steps(manager.steps_root) == [3]


def test_storage_estimate_is_bounded_for_rolling_policy():
    gb = 65_000_000
    rolling = estimate_checkpoint_storage(
        CheckpointPolicy(), nt=600, start_step=0, analysis_steps=range(200, 600, 20),
        bytes_per_checkpoint=gb, files_per_checkpoint=5,
    )
    full = estimate_checkpoint_storage(
        CheckpointPolicy(keep_last=0), nt=600, start_step=0,
        analysis_steps=range(200, 600, 20), bytes_per_checkpoint=gb, files_per_checkpoint=5,
    )
    assert rolling["checkpoints_retained"] == 3  # newest two + initial
    assert full["checkpoints_retained"] == 601
    assert full["retained_bytes"] > 100 * rolling["retained_bytes"]


def test_storage_preflight_warns_when_policy_exceeds_free_space(tmp_path, capsys):
    manager = _manager(tmp_path, nt=600, analysis_steps=[], policy=CheckpointPolicy(keep_last=0))
    manager.storage_preflight(
        start_step=0, bytes_per_checkpoint=10**15, files_per_checkpoint=5, include_initial=True
    )
    assert "WARNING" in capsys.readouterr().out


# --- 11: current writer stays readable by the original reader contract ------------

def test_new_manifest_keeps_format_and_fields_of_original_contract(tmp_path):
    manager = _manager(tmp_path, nt=4, analysis_steps=[1])
    checkpoint = manager.commit_step(1, _ensemble(1), did_analysis=True, analyses_completed=1)
    manifest = json.loads((checkpoint.path / "manifest.json").read_text())
    assert manifest["format"] == "icesee-distributed-block-checkpoint-v2"
    assert manifest["complete"] is True
    assert manifest["metadata"]["did_analysis"] is True
    assert manifest["metadata"]["completed_step"] == 1
    assert manifest["metadata"]["checkpoint_role"] == "restart"


def test_layout_fingerprint_mismatch_is_refused(tmp_path):
    saved = save_distributed_checkpoint(
        tmp_path, 0, _ensemble(0), _topology(), run_id=RUN_ID,
        metadata={"Nens": 3}, layout_fingerprint="mesh-A",
    )
    kwargs = dict(expected_run_id=RUN_ID, number_of_members=3)
    load_distributed_checkpoint(
        saved.path, _layout(), _topology(), expected_layout_fingerprint="mesh-A", **kwargs
    )
    with pytest.raises(ValueError, match="fingerprint"):
        load_distributed_checkpoint(
            saved.path, _layout(), _topology(), expected_layout_fingerprint="mesh-B", **kwargs
        )


def test_checkpoint_without_fingerprint_needs_undecomposed_numbering(tmp_path):
    steps = tmp_path / "steps"
    path = _write_legacy(steps, 3, did_analysis=False)
    kwargs = dict(expected_run_id=RUN_ID, number_of_members=3, expected_layout_fingerprint="mesh-A")
    load_distributed_checkpoint(path, _layout(), _topology(), **kwargs)
    decomposed = _topology()
    decomposed.spatial_ranks = 2
    with pytest.raises(ValueError, match="cannot be verified"):
        load_distributed_checkpoint(path, _layout(), decomposed, **kwargs)


# --- 13: generic code stays model-agnostic -------------------------------------

@pytest.mark.parametrize(
    "module",
    [
        "src/parallelization/distributed_checkpoint.py",
        "src/parallelization/mode3_checkpointing.py",
        "src/parallelization/distributed_native_runtime.py",
    ],
)
def test_generic_checkpoint_code_has_no_model_field_names(module):
    source = (Path(__file__).resolve().parents[2] / module).read_text()
    for token in ("thickness", "velocity", "basal_melt", "surface", "icepack", "firedrake",
                  '"h"', '"u"', '"v"', '"s"', "'h'", "'u'", "'v'", "'s'"):
        assert token not in source.lower(), f"{module} mentions {token}"
