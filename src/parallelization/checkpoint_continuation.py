"""Bounded Slurm continuation with durable, reconcilable submission records.

This supervisor does no numerical work. Invoke it periodically, or keep its CLI
running, on an approved submission host. It never changes model/filter settings.
"""
from __future__ import annotations

from contextlib import contextmanager
from datetime import date
import getpass
import hashlib
import json
import os
from pathlib import Path
import subprocess
import uuid

from .run_control import atomic_json, atomic_text


def source_fingerprint(root):
    """Include tracked modifications and untracked source; exclude run products."""
    root = Path(root).resolve()
    files = subprocess.check_output(["git", "-C", str(root), "ls-files", "-z",
                                     "--cached", "--others", "--exclude-standard"]).split(b"\0")
    digest = hashlib.sha256()
    for name in sorted(set(n for n in files if n)):
        relative = Path(os.fsdecode(name))
        if any(part in {"pace_runs", "_modelrun_datasets", "__pycache__", ".pytest_cache"}
               for part in relative.parts):
            continue
        path = root / relative
        digest.update(name + b"\0")
        if path.is_symlink():
            digest.update(b"symlink:" + os.fsencode(os.readlink(path)))
        elif path.is_file():
            digest.update(str(path.stat().st_mode & 0o777).encode() + b"\0")
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024**2), b""):
                    digest.update(chunk)
        elif not path.exists():
            digest.update(b"deleted")
        else:
            raise ValueError(f"source pin cannot cover nested repository/directory: {relative}")
    return digest.hexdigest()


def latest_progress(data_dir, *, nt, members, dt, expected_run_id=None, time_origin=0.0):
    """Inspect committed manifests only; native restore validates HDF5 contents."""
    steps = Path(data_dir) / "_mode3_state_history" / "steps"
    paths = sorted(steps.glob("checkpoint_[0-9]" + "[0-9]" * 7))
    for path in reversed(paths):
        manifest = path / "manifest.json"
        if not manifest.exists():
            continue
        value = json.loads(manifest.read_text())
        if value.get("complete") is not True:
            raise ValueError(f"uncommitted manifest at {path}")
        meta = value.get("metadata", {})
        next_step = int(meta.get("next_step", -1))
        if (not 0 < next_step <= nt or int(value.get("timestep", -2)) != next_step - 1
                or path.name != f"checkpoint_{next_step-1:08d}"
                or int(meta.get("nt", -1)) != nt
                or int(value.get("number_of_members", -1)) != members
                or abs(float(meta.get("model_time", -1)) - (time_origin + next_step * dt)) > 1e-8
                or not value.get("run_id")
                or (expected_run_id is not None and value["run_id"] != expected_run_id)):
            raise ValueError(f"incompatible continuation checkpoint: {path}")
        shards = value.get("shards", [])
        if not shards:
            raise ValueError(f"checkpoint has no shards: {path}")
        for shard in shards:
            name = shard.get("file", shard.get("filename", ""))
            if not name or Path(name).name != name:
                raise ValueError("invalid checkpoint shard path")
            target = path / name
            if not target.is_file() or target.stat().st_size == 0:
                raise ValueError(f"missing/empty checkpoint shard: {target}")
            if "nbytes" in shard and target.stat().st_size != int(shard["nbytes"]):
                raise ValueError(f"checkpoint shard size disagrees with manifest: {target}")
        return dict(next_step=next_step, run_id=value["run_id"],
                    manifest=str(manifest), sha256=hashlib.sha256(manifest.read_bytes()).hexdigest())
    return dict(next_step=0, run_id=expected_run_id)


class Slurm:
    def find(self, name, since):
        queued = subprocess.check_output(
            ["squeue", "--me", "--noheader", "--format=%i|%100j|%T"], text=True)
        accounted = subprocess.check_output(["sacct", "--user", getpass.getuser(),
            "--starttime", since, "--name", name, "--noheader", "--parsable2",
            "--format=JobIDRaw,JobName%100,State,ExitCode"], text=True)
        found = {}
        for line in accounted.splitlines():
            row = [part.strip() for part in line.split("|")]
            if len(row) >= 4 and row[0].isdigit() and row[1] == name:
                found[row[0]] = dict(id=row[0], state=row[2].split()[0].rstrip("+"), exit=row[3])
        for line in queued.splitlines():
            row = [part.strip() for part in line.split("|")]
            if len(row) == 3 and row[1] == name:
                found[row[0]] = dict(id=row[0], state=row[2], exit=None)
        if len(found) > 1:
            raise RuntimeError("multiple jobs match one submission nonce; manual review required")
        return next(iter(found.values()), None)

    def submit(self, command):
        result = subprocess.check_output(command, text=True).strip().split(";")[0]
        if not result.isdigit():
            raise RuntimeError("unrecognized sbatch response; reconcile before retrying")
        return result


@contextmanager
def case_lock(directory):
    lock = Path(directory) / ".supervisor-lock"
    lock.mkdir()  # Atomic across submission hosts; stale locks require inspection.
    try:
        (lock / "owner.json").write_text(json.dumps(dict(pid=os.getpid(), host=os.uname().nodename)))
        yield
    finally:
        (lock / "owner.json").unlink(missing_ok=True)
        lock.rmdir()


def tick(directory, *, scheduler=None, fingerprint=source_fingerprint):
    """Reconcile once and submit at most one successor. Never resubmit a nonce."""
    import yaml
    directory = Path(directory).resolve()
    scheduler = scheduler or Slurm()
    with case_lock(directory):
        spec = json.loads((directory / "case.json").read_text())
        state_path = directory / "state.json"
        state = json.loads(state_path.read_text()) if state_path.exists() else dict(submissions=0, phase="ready", next_step=0)
        if state["phase"] in {"complete", "stopped"}:
            return state
        if fingerprint(spec["source"]) != spec["source_sha256"]:
            raise RuntimeError("source changed since this case was pinned; no job submitted")
        parameters = Path(spec["parameters"])
        if hashlib.sha256(parameters.read_bytes()).hexdigest() != spec["parameters_sha256"]:
            raise RuntimeError("parameters changed since this case was pinned; no job submitted")
        progress = latest_progress(spec["data"], nt=spec["nt"], members=spec["members"],
                                   dt=spec["dt"], expected_run_id=state.get("run_id"),
                                   time_origin=spec.get("time_origin",0.0))
        if state["phase"] in {"submitting", "waiting"}:
            job = scheduler.find(state["job_name"], state["since"])
            if job is None:
                # An accepted sbatch can be invisible while accounting catches up.
                return dict(state, attention="submission outcome unknown; reconciliation required")
            if state.get("job_id") not in (None, job["id"]):
                raise RuntimeError("scheduler identity disagrees with the durable submission record")
            state.update(job_id=job["id"], phase="waiting", scheduler_state=job["state"])
            if job["state"] in {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED", "REQUEUED", "RESIZING"}:
                atomic_json(state_path, state)
                return state
            allowed = job["state"] in {"TIMEOUT", "NODE_FAIL", "PREEMPTED"}
            if job["state"] == "COMPLETED" and job["exit"] == "0:0":
                allowed = True
            if not allowed or progress["next_step"] <= state["start_step"]:
                state.update(phase="stopped", reason="failure or no committed progress; inspect logs")
                atomic_json(state_path, state)
                return state
        if progress["next_step"] == spec["nt"]:
            state.update(phase="complete", **progress)
            atomic_json(state_path, state)
            return state
        if state["submissions"] >= spec["max_submissions"]:
            state.update(phase="stopped", reason="submission budget exhausted", **progress)
            atomic_json(state_path, state)
            return state
        cfg = yaml.safe_load(parameters.read_text())
        options = cfg["enkf-parameters"]
        options.update(data_path=spec["data"], resume_from_checkpoint=progress["next_step"] > 0,
                       run_segment_cycles=spec["segment_cycles"], run_walltime_seconds=spec["run_seconds"])
        if progress["next_step"] > 0:
            options.update(initialize_ensemble=0, generate_true_state=0,
                           generate_synthetic_obs=0, generate_nurged_state=0)
        number = state["submissions"] + 1
        run = directory / f"segment-{number:04d}"
        run.mkdir()
        config = run / "params.yaml"
        atomic_text(config, yaml.safe_dump(cfg, sort_keys=False))
        job_name = "icesee-cont-" + uuid.uuid4().hex[:20]
        batch_options = [value + ",ICESEE_LOG_PATH=" + str(run / "logs")
                         if value.startswith("--export=") else value
                         for value in spec["sbatch_args"]]
        command = ["sbatch", "--parsable", "--job-name=" + job_name,
                   "--output=" + str(run / "batch-%j.out"),
                   *batch_options, spec["launcher"], "-F", str(config)]
        # Persist intent before the external side effect. An ambiguous failure
        # leaves this nonce in 'submitting', never creates a replacement job.
        state.update(phase="submitting", submissions=number, job_name=job_name,
                     since=date.today().isoformat(), start_step=progress["next_step"],
                     run_id=progress.get("run_id"), next_step=progress["next_step"], command=command, job_id=None)
        atomic_json(run / "submission.json", state)
        atomic_json(state_path, state)
        job_id = scheduler.submit(command)
        state.update(phase="waiting", job_id=job_id)
        atomic_json(run / "submission.json", state)
        atomic_json(state_path, state)
        return state
