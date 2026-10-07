"""Restart, retention, and storage policy for execution-mode-3 checkpoints.

Model-agnostic: this module sees only step indices, the observation
schedule, member counts, and opaque packed state handled by
``distributed_checkpoint``.  It never knows a model's field names.

Directory layout under ``<data_path>/_mode3_state_history/``::

    initial_condition/checkpoint_00000000/   pre-forecast ensemble (once)
    steps/checkpoint_<k>/                    ensemble after cycle k

Checkpoints written by this module carry ``metadata["checkpoint_role"]``:

``restart``
    rolling restart state; pruned to the newest ``checkpoint_keep_last``.
``analysis``
    an analysis-step checkpoint pinned by ``checkpoint_keep_analysis``.
``history``
    written with ``checkpoint_keep_last: 0`` (keep everything).

Only ``restart`` checkpoints are ever deleted automatically.  Checkpoints
without a role (written before retention existed) are never deleted.

Resume semantics: a checkpoint for step ``k`` holds the ensemble after the
complete cycle ``k`` (forecast plus, when scheduled, its analysis).  A
resumed run restores it and continues with cycle ``k + 1``; analyses at
steps ``<= k`` are never repeated.  A cycle interrupted before its
checkpoint commits is recomputed in full from step ``k``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import os
from pathlib import Path
import shutil
from typing import Any, Iterable, Mapping, Sequence

from .distributed_checkpoint import (
    CheckpointValidationError,
    DistributedCheckpoint,
    committed_checkpoint_paths,
    delete_committed_checkpoint,
    discover_latest_checkpoint,
    quarantine_checkpoint,
    remove_abandoned_checkpoint_artifacts,
    save_distributed_checkpoint,
    validate_distributed_checkpoint,
)

HISTORY_DIRNAME = "_mode3_state_history"
INITIAL_DIRNAME = "initial_condition"
STEPS_DIRNAME = "steps"

ROLE_INITIAL = "initial"
ROLE_RESTART = "restart"
ROLE_ANALYSIS = "analysis"
ROLE_HISTORY = "history"
PRUNABLE_ROLES = frozenset({ROLE_RESTART})


@dataclass(frozen=True)
class CheckpointPolicy:
    """When step checkpoints are written and how many are retained.

    ``every``: write after every ``every``-th completed cycle (the final
    cycle is always written).  ``keep_last``: retain the newest
    ``keep_last`` step checkpoints; ``0`` keeps every checkpoint written.
    ``keep_analysis``: additionally write and permanently keep the
    checkpoint of every analysis step.
    """

    every: int = 1
    keep_last: int = 2
    keep_analysis: bool = False

    def __post_init__(self) -> None:
        if int(self.every) < 1:
            raise ValueError("checkpoint_every must be >= 1")
        if int(self.keep_last) < 0:
            raise ValueError("checkpoint_keep_last must be >= 0 (0 keeps all)")

    @classmethod
    def from_kwargs(cls, icesee_kwargs: Mapping[str, Any]) -> "CheckpointPolicy":
        return cls(
            every=int(icesee_kwargs.get("checkpoint_every", 1)),
            keep_last=int(icesee_kwargs.get("checkpoint_keep_last", 2)),
            keep_analysis=bool(icesee_kwargs.get("checkpoint_keep_analysis", False)),
        )

    @property
    def keeps_all(self) -> bool:
        return int(self.keep_last) == 0

    def writes_step(self, step: int, nt: int, did_analysis: bool) -> bool:
        step = int(step)
        return (
            (step + 1) % int(self.every) == 0
            or step == int(nt) - 1
            or (bool(did_analysis) and bool(self.keep_analysis))
        )

    def role_for(self, did_analysis: bool) -> str:
        if self.keeps_all:
            return ROLE_HISTORY
        if did_analysis and self.keep_analysis:
            return ROLE_ANALYSIS
        return ROLE_RESTART

    def describe(self) -> str:
        retained = "all" if self.keeps_all else f"newest {int(self.keep_last)}"
        pinned = " + every analysis step" if self.keep_analysis else ""
        return f"every {int(self.every)} step(s), keep {retained}{pinned}"


def schedule_digest(nt: int, analysis_steps: Iterable[int]) -> str:
    """Identity of a run's step count and analysis schedule."""

    payload = json.dumps(
        {"nt": int(nt), "analysis_steps": sorted(int(s) for s in analysis_steps)}
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


@dataclass(frozen=True)
class ResumePlan:
    """Where a run (re)starts.

    ``completed_step`` is the last cycle whose result is restored (-1 when
    starting from the pre-forecast ensemble), ``start_step`` the first cycle
    to run, ``analyses_completed`` the analysis events already applied.
    """

    resumed: bool
    checkpoint_path: Path | None
    completed_step: int
    start_step: int
    analyses_completed: int
    legacy_checkpoint: bool = False
    rejected: tuple[tuple[str, str], ...] = ()
    quarantined: tuple[str, ...] = ()
    removed_artifacts: tuple[str, ...] = ()
    notes: tuple[str, ...] = field(default_factory=tuple)


def fresh_start_plan() -> ResumePlan:
    return ResumePlan(False, None, -1, 0, 0)


def _manifest_metadata(path: Path) -> Mapping[str, Any] | None:
    try:
        with (path / "manifest.json").open("r", encoding="utf-8") as stream:
            manifest = json.load(stream)
        return manifest
    except (OSError, ValueError):
        return None


def check_schedule_consistency(
    steps_root: Path,
    *,
    nt: int,
    analysis_steps: Sequence[int],
    number_of_members: int,
    upto_step: int,
) -> int:
    """Verify that committed checkpoints agree with this run's schedule.

    Every committed manifest up to ``upto_step`` is compared with the
    current configuration: member count, ``did_analysis`` against the
    analysis schedule, and (for checkpoints written by this module) the
    step count and schedule digest.  This is what makes a resume of an
    older checkpoint safe: the restored state is only continued if the
    configuration demonstrably reproduces the run that wrote it.  Returns
    the number of manifests checked; raises ValueError on a mismatch.
    """

    analysis = {int(step) for step in analysis_steps}
    digest = schedule_digest(nt, analysis)
    problems: list[str] = []
    checked = 0
    for step, path in committed_checkpoint_paths(steps_root):
        if step > upto_step:
            continue
        manifest = _manifest_metadata(path)
        if manifest is None:
            continue
        checked += 1
        metadata = manifest.get("metadata") or {}
        if step >= int(nt):
            problems.append(f"step {step} is beyond nt={nt}")
        if int(manifest.get("number_of_members", -1)) != int(number_of_members):
            problems.append(
                f"step {step} has {manifest.get('number_of_members')} members, "
                f"configured Nens={number_of_members}"
            )
        if "did_analysis" in metadata and bool(metadata["did_analysis"]) != (step in analysis):
            problems.append(
                f"step {step} did_analysis={bool(metadata['did_analysis'])} but the "
                f"current schedule says {step in analysis}"
            )
        if "nt" in metadata and int(metadata["nt"]) != int(nt):
            problems.append(f"step {step} was written for nt={metadata['nt']}, now nt={nt}")
        if "schedule_digest" in metadata and metadata["schedule_digest"] != digest:
            problems.append(f"step {step} was written under a different analysis schedule")
        if len(problems) >= 5:
            break
    if problems:
        raise ValueError(
            "existing checkpoints do not match this run's configuration; resume "
            "with exactly the original options (Nens, num_years, dt, observation "
            "window): " + "; ".join(problems)
        )
    return checked


def estimate_checkpoint_storage(
    policy: CheckpointPolicy,
    *,
    nt: int,
    start_step: int,
    analysis_steps: Iterable[int],
    bytes_per_checkpoint: int,
    files_per_checkpoint: int,
    include_initial: bool = True,
) -> dict[str, int]:
    """Projected checkpoint footprint of the steps this run will execute."""

    analysis = {int(step) for step in analysis_steps}
    steps = range(int(start_step), int(nt))
    written = sum(1 for k in steps if policy.writes_step(k, nt, k in analysis))
    if policy.keeps_all:
        retained = written
    else:
        pinned = (
            sum(1 for k in steps if k in analysis) if policy.keep_analysis else 0
        )
        retained = min(written, int(policy.keep_last) + pinned)
    retained += 1 if include_initial else 0
    # A new checkpoint is fully committed before an old one is removed.
    peak = retained + (0 if policy.keeps_all else 1)
    return {
        "checkpoints_written": written,
        "checkpoints_retained": retained,
        "retained_bytes": retained * int(bytes_per_checkpoint),
        "peak_bytes": peak * int(bytes_per_checkpoint),
        "retained_files": retained * int(files_per_checkpoint),
    }


def _gb(value: float) -> str:
    return f"{value / 1e9:.2f} GB"


class Mode3CheckpointManager:
    """Collective restart/retention driver for one mode-3 run.

    Filesystem scans, validation, and deletion happen on world rank 0 only;
    their outcome (or error) is broadcast so every rank takes the same
    branch.  Checkpoint writes themselves remain the collective
    ``save_distributed_checkpoint`` transaction.
    """

    def __init__(
        self,
        history_root: str | os.PathLike[str],
        *,
        run_id: str,
        topology: Any,
        number_of_members: int,
        nt: int,
        analysis_steps: Iterable[int],
        policy: CheckpointPolicy,
        time_grid: Sequence[float] | None = None,
        layout_fingerprint: str | None = None,
    ) -> None:
        self.history_root = Path(history_root)
        self.initial_root = self.history_root / INITIAL_DIRNAME
        self.steps_root = self.history_root / STEPS_DIRNAME
        self.run_id = str(run_id)
        self.topology = topology
        self.number_of_members = int(number_of_members)
        self.nt = int(nt)
        self.analysis_steps = tuple(sorted({int(s) for s in analysis_steps if int(s) < int(nt)}))
        self._analysis_set = frozenset(self.analysis_steps)
        self.policy = policy
        self.time_grid = None if time_grid is None else [float(v) for v in time_grid]
        self.layout_fingerprint = layout_fingerprint
        self.digest = schedule_digest(self.nt, self.analysis_steps)
        self.pruned: list[int] = []
        self.prune_warnings: list[str] = []

    # -- collective helpers -------------------------------------------------
    @property
    def _is_root(self) -> bool:
        return int(self.topology.world_rank) == 0

    def _on_root(self, function):
        result = None
        error = None
        if self._is_root:
            try:
                result = function()
            except Exception as exc:  # broadcast instead of deadlocking peers
                error = (type(exc).__name__, str(exc))
        result, error = self.topology.world.bcast((result, error), root=0)
        if error is not None:
            name, message = error
            exception = {"ValueError": ValueError, "FileNotFoundError": FileNotFoundError}.get(
                name, RuntimeError
            )
            raise exception(message)
        return result

    # -- planning -----------------------------------------------------------
    def analyses_through(self, step: int) -> int:
        return sum(1 for s in self.analysis_steps if s <= int(step))

    def model_time_after(self, step: int) -> float | None:
        """Physical time of the state after cycle ``step`` (``t[step + 1]``)."""

        if self.time_grid is None or not 0 <= int(step) + 1 < len(self.time_grid):
            return None
        return self.time_grid[int(step) + 1]

    def prepare_resume(self) -> ResumePlan:
        """Find the newest valid committed checkpoint and plan continuation.

        Removes abandoned staging/deleting debris, validates candidates from
        newest to oldest, checks the configuration against existing
        manifests, and un-publishes (renames, never deletes) any invalid
        committed directory newer than the selected checkpoint so its step
        can be rewritten.  Raises FileNotFoundError when nothing valid
        exists: an explicit resume never silently starts over.
        """

        return self._on_root(self._prepare_resume_root)

    def _prepare_resume_root(self) -> ResumePlan:
        removed = []
        for root in (self.steps_root, self.initial_root):
            removed += [p.name for p in remove_abandoned_checkpoint_artifacts(root)]
        # The ensemble size is checked after selection (below) so that a
        # mismatch is reported as such instead of as "nothing to resume".
        discovery = discover_latest_checkpoint(self.steps_root, expected_run_id=self.run_id)
        notes = []
        if discovery.path is not None:
            completed = int(discovery.timestep)
            check_schedule_consistency(
                self.steps_root,
                nt=self.nt,
                analysis_steps=self.analysis_steps,
                number_of_members=self.number_of_members,
                upto_step=completed,
            )
            metadata = (discovery.manifest or {}).get("metadata") or {}
            legacy = "completed_step" not in metadata
            done = self.analyses_through(completed)
            if "analysis_events_completed" in metadata and int(
                metadata["analysis_events_completed"]
            ) != done:
                raise ValueError(
                    f"checkpoint {discovery.path.name} records "
                    f"{metadata['analysis_events_completed']} completed analyses; "
                    f"the current schedule implies {done}"
                )
            path = discovery.path
        else:
            initial = discover_latest_checkpoint(self.initial_root, expected_run_id=self.run_id)
            if initial.path is None:
                reasons = "; ".join(f"{n}: {r}" for n, r in discovery.rejected + initial.rejected)
                raise FileNotFoundError(
                    f"resume requested but no valid committed checkpoint exists under "
                    f"{self.history_root}" + (f" (rejected: {reasons})" if reasons else "")
                )
            completed, done, path = -1, 0, initial.path
            legacy = "completed_step" not in ((initial.manifest or {}).get("metadata") or {})
            notes.append("no step checkpoint; resuming from the pre-forecast ensemble")
            members = int((initial.manifest or {}).get("number_of_members", -1))
            if members != self.number_of_members:
                raise ValueError(
                    f"initial checkpoint holds {members} members; configured "
                    f"Nens={self.number_of_members}"
                )
        quarantined = []
        for step, candidate in committed_checkpoint_paths(self.steps_root):
            if step > completed:
                # Discovery already rejected it.  Structurally damaged
                # directories are moved aside; an intact checkpoint of a
                # different run or ensemble size is left alone and stops the
                # resume, because continuing would overwrite its step.
                try:
                    validate_distributed_checkpoint(candidate, deep=True)
                except CheckpointValidationError:
                    quarantined.append(quarantine_checkpoint(candidate).name)
                    continue
                raise ValueError(
                    f"{candidate.name} is an intact checkpoint of a different run_id "
                    "or ensemble size; refusing to resume into the same directory"
                )
        return ResumePlan(
            resumed=True,
            checkpoint_path=path,
            completed_step=completed,
            start_step=completed + 1,
            analyses_completed=done,
            legacy_checkpoint=legacy,
            rejected=tuple(discovery.rejected),
            quarantined=tuple(quarantined),
            removed_artifacts=tuple(removed),
            notes=tuple(notes),
        )

    def report_plan(self, plan: ResumePlan) -> None:
        if not self._is_root:
            return
        if not plan.resumed:
            return
        lines = [
            f"checkpoint        = {plan.checkpoint_path}",
            f"completed step    = {plan.completed_step}"
            + (
                f" (t = {self.model_time_after(plan.completed_step):g})"
                if self.model_time_after(plan.completed_step) is not None
                else ""
            ),
            f"next step         = {plan.start_step} of nt = {self.nt}",
            f"analyses done     = {plan.analyses_completed} of {len(self.analysis_steps)}",
            "checkpoint format = "
            + ("pre-retention (no restart metadata)" if plan.legacy_checkpoint else "current"),
        ]
        lines += [f"rejected          = {name}: {reason}" for name, reason in plan.rejected]
        lines += [f"quarantined       = {name}" for name in plan.quarantined]
        lines += [f"removed debris    = {name}" for name in plan.removed_artifacts]
        lines += [f"note              = {note}" for note in plan.notes]
        print("[ICESEE] Mode-3 resume:\n  " + "\n  ".join(lines), flush=True)

    # -- writing ------------------------------------------------------------
    def _metadata(self, step: int, *, role: str, did_analysis: bool, analyses: int,
                  extra: Mapping[str, Any] | None) -> dict[str, Any]:
        metadata = {
            "Nens": self.number_of_members,
            "did_analysis": bool(did_analysis),
            "checkpoint_role": role,
            "completed_step": int(step),
            "next_step": int(step) + 1,
            "nt": self.nt,
            "analysis_events_completed": int(analyses),
            "schedule_digest": self.digest,
        }
        time = self.model_time_after(step)
        if time is not None:
            metadata["model_time"] = time
        metadata.update(dict(extra or {}))
        return metadata

    def write_initial(self, snapshot: Any, *, extra: Mapping[str, Any] | None = None):
        metadata = self._metadata(
            -1, role=ROLE_INITIAL, did_analysis=False, analyses=0, extra=extra
        )
        metadata["completed_step"] = -1
        metadata["next_step"] = 0
        return save_distributed_checkpoint(
            self.initial_root,
            0,
            snapshot,
            self.topology,
            run_id=self.run_id,
            metadata=metadata,
            layout_fingerprint=self.layout_fingerprint,
        )

    def should_write(self, step: int, did_analysis: bool) -> bool:
        return self.policy.writes_step(step, self.nt, did_analysis)

    def commit_step(
        self,
        step: int,
        snapshot: Any,
        *,
        did_analysis: bool,
        analyses_completed: int,
        extra: Mapping[str, Any] | None = None,
    ) -> DistributedCheckpoint | None:
        """Write step ``step`` if the policy asks for it, then apply retention.

        Retention runs only after ``save_distributed_checkpoint`` returned,
        i.e. after the new checkpoint is committed on every rank; a failed
        write raises before anything is deleted.
        """

        if not self.should_write(step, did_analysis):
            return None
        checkpoint = save_distributed_checkpoint(
            self.steps_root,
            step,
            snapshot,
            self.topology,
            run_id=self.run_id,
            metadata=self._metadata(
                step,
                role=self.policy.role_for(did_analysis),
                did_analysis=did_analysis,
                analyses=analyses_completed,
                extra=extra,
            ),
            layout_fingerprint=self.layout_fingerprint,
        )
        if self._is_root:
            self._prune_root(newest=int(step))
        return checkpoint

    def _prune_root(self, *, newest: int) -> None:
        """Delete surplus rolling restart checkpoints (rank 0, non-fatal)."""

        if self.policy.keeps_all:
            return
        committed = committed_checkpoint_paths(self.steps_root)
        keep = {step for step, _ in committed[: int(self.policy.keep_last)]}
        keep.add(int(newest))
        for step, path in committed:
            if step in keep:
                continue
            manifest = _manifest_metadata(path)
            role = ((manifest or {}).get("metadata") or {}).get("checkpoint_role")
            if role not in PRUNABLE_ROLES:
                continue
            try:
                delete_committed_checkpoint(path)
                self.pruned.append(step)
            except OSError as error:
                warning = f"could not remove {path.name}: {error}"
                self.prune_warnings.append(warning)
                print(f"[ICESEE] WARNING: checkpoint retention {warning}", flush=True)

    # -- storage ------------------------------------------------------------
    def storage_preflight(
        self,
        *,
        start_step: int,
        bytes_per_checkpoint: int,
        files_per_checkpoint: int,
        include_initial: bool,
        warn_fraction: float = 0.5,
    ) -> dict[str, int] | None:
        """Print the projected checkpoint footprint; warn if it is excessive.

        Advisory only: it never stops a run.
        """

        if not self._is_root:
            return None
        estimate = estimate_checkpoint_storage(
            self.policy,
            nt=self.nt,
            start_step=start_step,
            analysis_steps=self.analysis_steps,
            bytes_per_checkpoint=bytes_per_checkpoint,
            files_per_checkpoint=files_per_checkpoint,
            include_initial=include_initial,
        )
        probe = self.history_root
        while not probe.exists() and probe != probe.parent:
            probe = probe.parent
        free = shutil.disk_usage(probe).free
        print(
            "[ICESEE] Mode-3 checkpoints: "
            f"{self.policy.describe()}; {_gb(bytes_per_checkpoint)} each; "
            f"this run writes {estimate['checkpoints_written']}, retains "
            f"<= {estimate['checkpoints_retained']} (~{_gb(estimate['retained_bytes'])}, "
            f"~{estimate['retained_files']} files), peak ~{_gb(estimate['peak_bytes'])}; "
            f"free space {_gb(free)}",
            flush=True,
        )
        if estimate["peak_bytes"] > free:
            print(
                "[ICESEE] WARNING: the checkpoint retention policy is projected to "
                "exceed the free space on this filesystem. Reduce "
                "checkpoint_keep_last, disable checkpoint_keep_analysis, or raise "
                "checkpoint_every.",
                flush=True,
            )
        elif estimate["peak_bytes"] > warn_fraction * free:
            print(
                "[ICESEE] WARNING: checkpoints are projected to use more than "
                f"{int(100 * warn_fraction)}% of the free space on this filesystem.",
                flush=True,
            )
        return estimate


def _inspect(history_root: str) -> int:
    """Read-only summary of a mode-3 checkpoint directory (no MPI, no model)."""

    root = Path(history_root)
    steps_root = root / STEPS_DIRNAME
    committed = committed_checkpoint_paths(steps_root)
    debris = sorted(
        p.name for p in (steps_root.iterdir() if steps_root.is_dir() else [])
        if p.name.startswith(".")
    )
    size = sum(f.stat().st_size for f in root.rglob("*") if f.is_file()) if root.is_dir() else 0
    print(f"history directory : {root}")
    print(f"committed steps   : {len(committed)}"
          + (f" ({committed[-1][0]} .. {committed[0][0]})" if committed else ""))
    print(f"uncommitted debris: {', '.join(debris) if debris else 'none'}")
    print(f"total size        : {_gb(size)}")
    found = discover_latest_checkpoint(steps_root)
    for name, reason in found.rejected:
        print(f"rejected          : {name}: {reason}")
    if found.path is None:
        initial = discover_latest_checkpoint(root / INITIAL_DIRNAME)
        if initial.path is None:
            print("latest valid      : none -- nothing to resume")
            return 1
        print(f"latest valid      : {initial.path} (pre-forecast ensemble; resume at step 0)")
        return 0
    metadata = (found.manifest or {}).get("metadata") or {}
    print(f"latest valid      : {found.path}")
    print(f"members / run_id  : {found.manifest['number_of_members']} / {found.manifest['run_id']}")
    print(f"analysis at {found.timestep:<6}: {metadata.get('did_analysis')}")
    print(f"next step to run  : {found.timestep + 1} (nothing left if this equals nt)")
    return 0


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("usage: python -m ICESEE.src.parallelization.mode3_checkpointing "
              "<data_path>/_mode3_state_history")
        raise SystemExit(2)
    raise SystemExit(_inspect(sys.argv[1]))
