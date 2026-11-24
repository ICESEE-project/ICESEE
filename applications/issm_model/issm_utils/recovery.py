"""Transactional ISSM forecasts, including the MATLAB model cache and RNG."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import uuid

import numpy as np

from ICESEE.src.parallelization.forecast_recovery import (
    ForecastRecoveryPolicy, RecoverableForecastError, advance_with_recovery,
)
from ICESEE.src.parallelization.run_control import atomic_json, atomic_text

NUMERICAL_REJECTION = "ICESEE:ISSM:ForecastRejected"
HELPERS = Path(__file__).with_name("matlab2python")


def matlab_string(value):
    return "'" + str(value).replace("'", "''") + "'"


def checked_command(server, command, options):
    """Only an explicit numerical rejection is eligible for timestep retry."""
    root = Path(options["issm_examples_dir"]).resolve()
    for name in ("issm_recovery_command.m", "issm_recovery_rng.m"):
        target = root / name
        if not target.exists() or target.read_bytes() != (HELPERS / name).read_bytes():
            atomic_text(target, (HELPERS / name).read_text())
    timeout = float(options.get("issm_command_timeout", 600))
    if not np.isfinite(timeout) or timeout <= 0:
        raise ValueError("issm_command_timeout must be finite and positive")
    errors = root / "recovery_errors"
    errors.mkdir(exist_ok=True)
    error_file = errors / (uuid.uuid4().hex + ".json")
    wrapped = f"issm_recovery_command({matlab_string(command)}, {matlab_string(error_file)});"
    if server.send_command(wrapped, timeout=timeout):
        return
    if error_file.is_file():
        detail = json.loads(error_file.read_text())
        if detail.get("identifier") == NUMERICAL_REJECTION:
            raise RecoverableForecastError(detail.get("message", "ISSM rejected forecast"))
        raise RuntimeError("non-retryable MATLAB error: " + str(detail))
    raise RuntimeError("MATLAB command failed or timed out without a numerical rejection; no retry")


def _files(directory):
    if not directory.exists():
        return []
    result = []
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            raise ValueError("ISSM recovery refuses symlinks in member model state")
        if path.is_file():
            result.append(path)
    return result


class MatlabMemberTransaction:
    """Capture model files and MATLAB RNG before a potentially rejected step."""

    def __init__(self, options, member_id):
        self.options = options
        self.member_id = int(member_id)
        if self.member_id < 0:
            raise ValueError("ISSM member ID must be nonnegative")
        self.workspace = Path(options["issm_examples_dir"]).resolve()
        self.models = self.workspace / "Models" / f"ens_id_{self.member_id}"
        if self.models.is_symlink():
            raise ValueError("ISSM member model directory cannot be a symlink")
        self.data = (Path(options["icesee_path"]) / options["data_path"]).resolve()
        self.server = options["server"]
        self.marker = self.data / "_issm_recovery_accepted" / f"member_{self.member_id}.json"
        self.limit = float(options.get("issm_recovery_snapshot_mb", 1024)) * 1024**2
        if not np.isfinite(self.limit) or self.limit <= 0:
            raise ValueError("issm_recovery_snapshot_mb must be finite and positive")

    def output_files(self):
        return [self.data/name for name in (
            f"ensemble_output_{self.member_id}.h5", f"ensemble_out_scalar_{self.member_id}.h5")
            if (self.data/name).exists() or (self.data/name).is_symlink()]

    def model_files(self):
        return [self.models/name for name in ("enkf_state.mat", "initialize_ensemble.mat", "inverse_state.mat")
                if (self.models/name).exists() or (self.models/name).is_symlink()]

    def capture(self, destination):
        destination = Path(destination)
        source_files = self.model_files() + self.output_files()
        if any(path.is_symlink() for path in source_files):
            raise ValueError("ISSM recovery refuses symlinked member files")
        if sum(path.stat().st_size for path in source_files) > self.limit:
            raise MemoryError("ISSM member snapshot exceeds its configured byte budget")
        destination.mkdir(parents=True, exist_ok=False)
        record = dict(member_id=self.member_id, model_present=self.models.exists(), files={})
        if self.models.exists():
            (destination / "models").mkdir()
            for path in self.model_files():
                shutil.copyfile(path,destination/"models"/path.name)
        (destination / "outputs").mkdir()
        for path in self.output_files():
            shutil.copyfile(path, destination / "outputs" / path.name)
        if self.marker.exists():
            shutil.copyfile(self.marker, destination / "accepted.json")
        checked_command(self.server,
            f"issm_recovery_rng('save', {matlab_string(destination / 'rng.mat')});", self.options)
        for path in _files(destination):
            with path.open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            record["files"][str(path.relative_to(destination))] = dict(sha256=digest, size=path.stat().st_size)
            with path.open("rb") as stream:
                os.fsync(stream.fileno())
        atomic_json(destination / "snapshot.json", record)
        return destination

    def restore(self, snapshot):
        snapshot = Path(snapshot)
        record = json.loads((snapshot / "snapshot.json").read_text())
        if record["member_id"] != self.member_id:
            raise ValueError("ISSM backend snapshot belongs to another member")
        for name, signature in record["files"].items():
            source = snapshot / name
            if Path(name).is_absolute() or ".." in Path(name).parts or source.is_symlink():
                raise ValueError("unsafe ISSM backend snapshot path")
            if source.stat().st_size != signature["size"]:
                raise ValueError("ISSM backend snapshot size mismatch")
            with source.open("rb") as stream:
                if hashlib.file_digest(stream, "sha256").hexdigest() != signature["sha256"]:
                    raise ValueError("ISSM backend snapshot checksum mismatch")
        for path in self.model_files():
            path.unlink()
        if record["model_present"]:
            self.models.mkdir(parents=True,exist_ok=True)
            for path in (snapshot/"models").iterdir():
                shutil.copyfile(path,self.models/path.name)
        for path in self.output_files():
            path.unlink()
        for path in (snapshot / "outputs").iterdir():
            shutil.copyfile(path, self.data / path.name)
        if (snapshot / "accepted.json").exists():
            self.marker.parent.mkdir(exist_ok=True)
            shutil.copyfile(snapshot / "accepted.json", self.marker)
        else:
            self.marker.unlink(missing_ok=True)
        checked_command(self.server,
            f"issm_recovery_rng('restore', {matlab_string(snapshot / 'rng.mat')});", self.options)


def advance_issm_forecast(ensemble, *, forecast, options):
    """Run a whole-member physical forecast; noise stays in ICESEE's caller."""
    from ICESEE.src.utils.tools import icesee_get_index
    policy = ForecastRecoveryPolicy.from_options(options)
    if policy is None:
        raise ValueError("advance_issm_forecast requires an explicit forecast_recovery policy")
    comm = options.get("comm")
    if comm is not None and comm.Get_size() != 1:
        raise ValueError("ISSM rollback requires one Python rank per complete member")
    initial = np.asarray(ensemble).copy()
    working = [initial.copy()]
    result = [{}]
    _, indices, _ = icesee_get_index(initial, **options)
    step = int(options["k"])
    grid = np.asarray(options["t"], dtype=float)
    start, stop = float(grid[step]), float(grid[step + 1])
    transaction = MatlabMemberTransaction(options, options["ens_id"])
    scratch_root = transaction.data / "_recovery_work"
    scratch_root.mkdir(parents=True, exist_ok=True)
    old_cwd = Path.cwd()
    with tempfile.TemporaryDirectory(prefix=f"member-{transaction.member_id}-", dir=scratch_root) as temp:
        snapshot_path = Path(temp) / "snapshot"
        def capture():
            transaction.capture(snapshot_path)
            return initial.copy()
        def restore(values):
            transaction.restore(snapshot_path)
            working[0] = values.copy()
        def advance(time, dt):
            kwargs = dict(options, dt=dt, tinitial=time, tfinal=time + dt,
                          _issm_recovery_active=True)
            result[0] = forecast(working[0].copy(), **kwargs)
            if not isinstance(result[0], dict):
                raise RuntimeError("ISSM forecast returned no accepted state")
            for name in options["vec_inputs"]:
                if name in result[0]:
                    working[0][indices[name]] = np.asarray(result[0][name]).ravel()
                elif name in {"Thickness", "Surface", "Vx", "Vy"}:
                    raise RuntimeError("ISSM forecast omitted physical field: " + name)
        def validate():
            if not np.isfinite(working[0]).all():
                raise RecoverableForecastError("ISSM returned a non-finite forecast")
        try:
            event = advance_with_recovery(start_time=start, interval=stop-start,
                policy=policy, capture=capture, restore=restore, advance=advance,
                validate=validate, comm=comm, stable_dt=options.get("issm_forecast_max_dt"))
            print("[ICESEE] issm_forecast_recovery " + json.dumps(dict(event,
                member=transaction.member_id, timestep=step)), flush=True)
        finally:
            os.chdir(old_cwd)
    publish_scalar_diagnostics(options, step)
    return {name:working[0][indices[name]].copy() for name in options["vec_inputs"]}


def publish_scalar_diagnostics(options, step):
    """The last accepted member publishes scalar means, without MPI in trials."""
    import time
    import h5py
    names = options.get("scalar_inputs", [])
    if not names:
        return
    data = (Path(options["icesee_path"]) / options["data_path"]).resolve()
    markers = data / "_issm_recovery_accepted"
    atomic_json(markers / f"member_{int(options['ens_id'])}.json", dict(step=int(step)))
    lock = markers / "aggregate.lock"
    deadline = time.monotonic() + 30
    while True:
        try:
            lock.mkdir()
            break
        except FileExistsError:
            if time.monotonic() > deadline:
                raise RuntimeError("ISSM scalar diagnostic writer remained locked")
            time.sleep(.02)
    try:
        for member in range(int(options["Nens"])):
            marker = markers / f"member_{member}.json"
            if not marker.exists() or json.loads(marker.read_text())["step"] != step:
                return
        with h5py.File(data / "ensemble_scalar_output.h5", "a") as target:
            for name in names:
                values = []
                for member in range(int(options["Nens"])):
                    with h5py.File(data / f"ensemble_out_scalar_{member}.h5") as source:
                        value = np.asarray(source[name][:]).ravel()[step]
                        if not np.isfinite(value):
                            raise RuntimeError("non-finite accepted ISSM scalar diagnostic")
                        values.append(value)
                if name not in target:
                    target.create_dataset(name,shape=(int(options["nt"]),),dtype="f8",fillvalue=np.nan)
                target[name][step] = np.mean(values)
    finally:
        lock.rmdir()


def run_checked_forecast_command(command, options):
    """Opt-in transport used before the legacy shutdown-on-error handlers."""
    workspace = Path(options["issm_examples_dir"])
    payload = workspace / f"recovery-options-{uuid.uuid4().hex}.mat"
    from scipy.io import savemat
    execution = workspace / "recovery_execution" / payload.stem
    execution.mkdir(parents=True)
    savemat(payload, dict(recovery_active=True, recovery_execution=str(execution)))
    finished = False
    try:
        # Append one optional argument to the existing MATLAB run_model call.
        command = command.rstrip().rstrip(";").rstrip()
        if not command.endswith(")"):
            raise ValueError("unexpected ISSM MATLAB forecast command")
        checked_command(options["server"], command[:-1] + ", " + matlab_string(payload) + ");", options)
        finished = True
    except RecoverableForecastError:
        finished = True
        raise
    finally:
        payload.unlink(missing_ok=True)
        if finished:
            shutil.rmtree(execution)


def recovery_identity(options):
    """Keep operational stop/resume controls out of the experiment identity."""
    excluded = {"server", "comm", "icesee_comm", "model_module", "k", "km", "rank",
        "sub_rank", "ens_id", "icesee_path", "issm_examples_dir", "issm_dir", "data_path",
        "example_name", "run_stop_time", "run_segment_cycles", "run_walltime_seconds",
        "resume_from_checkpoint", "initialize_ensemble", "generate_true_state",
        "generate_nurged_state", "generate_synthetic_obs", "force_fresh_start",
        "verbose", "cli_overrides", "tinitial", "tfinal", "run_id"}
    excluded.update({"icesee_rank", "icesee_size", "cluster_name", "restart_enabled",
                     "true_nurged_file", "synthetic_obs_file"})
    payload = {}
    def encode(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value).__name__)
    for key, value in options.items():
        if key in excluded or key.startswith("_"):
            continue
        try:
            encoded = json.dumps(value, sort_keys=True, default=encode)
        except (TypeError, ValueError):
            continue
        payload[key] = json.loads(encoded)
    source = Path(options.get("icesee_path", "."))
    payload["source_files"] = {name:hashlib.sha256((source/name).read_bytes()).hexdigest()
        for name in ("run_model.m", "_issm_model.py", "_issm_enkf.py") if (source/name).is_file()}
    return "issm-backend-v1-" + hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()[:24]


def commit_backend(pool, topology, options, timestep):
    """Durable backend state precedes the vector checkpoint that references it."""
    world = topology.world
    token = world.bcast(uuid.uuid4().hex if topology.world_rank == 0 else None, root=0)
    data = (Path(options["icesee_path"]) / options["data_path"]).resolve()
    root = data / "_issm_backend" / f"step_{timestep+1:08d}_{token}"
    problem = None
    owned = []
    try:
        root.mkdir(parents=True,exist_ok=True)
        for member in pool.member_ids:
            MatlabMemberTransaction(options,member).capture(root/f"member_{member}")
            owned.append(int(member))
    except Exception as exc:
        problem = f"rank {topology.world_rank}: {type(exc).__name__}: {exc}"
    problems = world.allgather(problem)
    if any(problems):
        raise RuntimeError("ISSM backend checkpoint failed: " + "; ".join(p for p in problems if p))
    members = [m for part in world.allgather(owned) for m in part]
    problem = None
    descriptor = None
    if topology.world_rank == 0:
        try:
            if sorted(members) != list(range(int(options["Nens"]))):
                raise ValueError("ISSM backend checkpoint member coverage mismatch")
            scalar = data / "ensemble_scalar_output.h5"
            scalar_digest = None
            if scalar.exists():
                shutil.copyfile(scalar,root/scalar.name)
                scalar_digest = hashlib.sha256((root/scalar.name).read_bytes()).hexdigest()
            manifests = {str(m):hashlib.sha256((root/f"member_{m}/snapshot.json").read_bytes()).hexdigest() for m in members}
            atomic_json(root/"manifest.json",dict(format="icesee-issm-backend-v1",complete=True,members=manifests,
                python_ranks=int(world.Get_size()),timestep=int(timestep),scalar_sha256=scalar_digest))
            descriptor = dict(path=str(root.relative_to(data)),
                sha256=hashlib.sha256((root/"manifest.json").read_bytes()).hexdigest())
        except Exception as exc:
            problem = f"{type(exc).__name__}: {exc}"
    problem,descriptor = world.bcast((problem,descriptor),root=0)
    if problem:
        raise RuntimeError("ISSM backend manifest failed: " + problem)
    return descriptor


def prune_backend(manager, topology, options):
    """Prune only committed backend bundles no retained vector references."""
    if topology.world_rank != 0 or manager.policy.keeps_all:
        return
    try:
        _prune_backend_root(manager,options)
    except Exception as exc:
        print(f"[ICESEE] Backend retention kept extra snapshots: {type(exc).__name__}: {exc}",flush=True)


def _prune_backend_root(manager, options):
    data = (Path(options["icesee_path"]) / options["data_path"]).resolve()
    referenced = set()
    for root in (manager.initial_root,manager.steps_root):
        for path in root.glob("checkpoint_????????/manifest.json"):
            metadata=json.loads(path.read_text()).get("metadata",{})
            descriptor=metadata.get("issm_backend")
            if descriptor:referenced.add(descriptor["path"])
    for path in (data/"_issm_backend").glob("step_????????_*"):
        manifest=path/"manifest.json"
        if str(path.relative_to(data)) in referenced or not manifest.exists():
            continue
        value=json.loads(manifest.read_text())
        if value.get("format")=="icesee-issm-backend-v1" and value.get("complete"):
            shutil.rmtree(path)


def restore_backend_member(member_id, checkpoint, topology, options):
    descriptor = checkpoint.metadata.get("issm_backend")
    if not descriptor:
        raise ValueError("ISSM resume requires a checkpoint with complete MATLAB backend state")
    relative = Path(descriptor["path"])
    if relative.is_absolute() or ".." in relative.parts or relative.parts[0] != "_issm_backend":
        raise ValueError("unsafe ISSM backend checkpoint reference")
    data = (Path(options["icesee_path"]) / options["data_path"]).resolve()
    root = data / relative
    raw = (root/"manifest.json").read_bytes()
    if hashlib.sha256(raw).hexdigest() != descriptor["sha256"]:
        raise ValueError("ISSM backend manifest checksum mismatch")
    manifest = json.loads(raw)
    if not manifest["complete"] or manifest["python_ranks"] != topology.world.Get_size():
        raise ValueError("ISSM backend resume currently requires the original Python rank count")
    snapshot = root / f"member_{int(member_id)}"
    if hashlib.sha256((snapshot/"snapshot.json").read_bytes()).hexdigest() != manifest["members"][str(member_id)]:
        raise ValueError("ISSM member snapshot manifest checksum mismatch")
    MatlabMemberTransaction(options,member_id).restore(snapshot)
    if topology.world_rank == 0 and manifest["scalar_sha256"]:
        source = root / "ensemble_scalar_output.h5"
        if hashlib.sha256(source.read_bytes()).hexdigest() != manifest["scalar_sha256"]:
            raise ValueError("ISSM scalar checkpoint checksum mismatch")
        shutil.copyfile(source,data/source.name)
