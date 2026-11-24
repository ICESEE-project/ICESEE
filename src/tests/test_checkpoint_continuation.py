import hashlib
import json
from pathlib import Path
import pytest
import yaml

from ICESEE.src.parallelization.checkpoint_continuation import (
    tick, latest_progress, case_lock, source_fingerprint, Slurm,
)


class Scheduler:
    def __init__(self):
        self.commands = []
        self.job = None
        self.ambiguous = False
    def find(self, name, since):
        return self.job
    def submit(self, command):
        self.commands.append(command)
        if self.ambiguous:
            raise TimeoutError("connection lost after submission")
        self.job = dict(id="1234", state="RUNNING", exit=None)
        return "1234"


@pytest.fixture
def case(tmp_path):
    params = tmp_path / "parameters.yaml"
    params.write_text(yaml.safe_dump({"modeling-parameters":{"num_years":82,"timesteps_per_year":.05},
        "enkf-parameters":{"execution_mode":3,"Nens":32,"sig_Q":[1,2],"seed":11}}))
    spec = dict(source="/source", source_sha256="pin", parameters=str(params),
        parameters_sha256=hashlib.sha256(params.read_bytes()).hexdigest(),
        data=str(tmp_path / "data"), nt=1640, members=32, dt=.05,
        segment_cycles=20, run_seconds=3000, max_submissions=3,
        launcher="/source/run.sbatch", sbatch_args=["--ntasks=8"])
    (tmp_path / "case.json").write_text(json.dumps(spec))
    return tmp_path


def checkpoint(case, step, *, run_id="science-policy-v1"):
    path = case / "data" / "_mode3_state_history" / "steps" / f"checkpoint_{step:08d}"
    path.mkdir(parents=True)
    (path / "shard_0.h5").write_bytes(b"placeholder; real HDF5 is verified by native restore")
    manifest = dict(complete=True,timestep=step,number_of_members=32,run_id=run_id,
        metadata=dict(nt=1640,next_step=step+1,model_time=(step+1)*.05),shards=[{"file":"shard_0.h5"}])
    (path / "manifest.json").write_text(json.dumps(manifest))
    return path


def run(case, scheduler, fingerprint=lambda source:"pin"):
    return tick(case, scheduler=scheduler, fingerprint=fingerprint)


def test_no_duplicate_while_running_and_successor_preserves_science(case):
    scheduler = Scheduler()
    first = run(case, scheduler)
    assert first["phase"] == "waiting"
    assert run(case, scheduler)["submissions"] == 1
    assert len(scheduler.commands) == 1
    checkpoint(case, 19)
    scheduler.job = dict(id="1234",state="COMPLETED",exit="0:0")
    second = run(case, scheduler)
    assert second["submissions"] == 2
    cfg = yaml.safe_load((case / "segment-0002" / "params.yaml").read_text())
    options = cfg["enkf-parameters"]
    assert options["resume_from_checkpoint"] and options["initialize_ensemble"] == 0
    assert options["sig_Q"] == [1,2] and options["seed"] == 11
    assert cfg["modeling-parameters"]["num_years"] == 82
    assert "run_stop_time" not in options


def test_ambiguous_submit_is_reconciled_never_retried(case):
    scheduler = Scheduler()
    scheduler.ambiguous = True
    with pytest.raises(TimeoutError):
        run(case, scheduler)
    result = run(case, scheduler)
    assert result["phase"] == "submitting" and "attention" in result
    assert len(scheduler.commands) == 1
    scheduler.job = dict(id="1234",state="RUNNING",exit=None)
    assert run(case, scheduler)["job_id"] == "1234"
    assert len(scheduler.commands) == 1


@pytest.mark.parametrize("state,code,retry", [("TIMEOUT","0:15",True),
    ("NODE_FAIL","1:0",True),("PREEMPTED","0:15",True),
    ("FAILED","1:0",False),("OUT_OF_MEMORY","0:9",False),
    ("CANCELLED","0:15",False),("COMPLETED","1:0",False)])
def test_retry_classification_does_not_loop_on_scientific_or_input_failures(case,state,code,retry):
    scheduler = Scheduler()
    run(case, scheduler)
    checkpoint(case, 19)
    scheduler.job = dict(id="1234",state=state,exit=code)
    outcome = run(case, scheduler)
    assert len(scheduler.commands) == (2 if retry else 1)
    assert outcome["phase"] == ("waiting" if retry else "stopped")


def test_no_progress_stops_even_after_transient_scheduler_failure(case):
    scheduler = Scheduler()
    run(case, scheduler)
    scheduler.job = dict(id="1234",state="TIMEOUT",exit="0:15")
    assert run(case, scheduler)["phase"] == "stopped"
    assert len(scheduler.commands) == 1


def test_budgets_and_completion(case):
    scheduler = Scheduler()
    spec = json.loads((case / "case.json").read_text())
    spec["max_submissions"] = 1
    (case / "case.json").write_text(json.dumps(spec))
    run(case, scheduler)
    checkpoint(case, 19)
    scheduler.job = dict(id="1234",state="COMPLETED",exit="0:0")
    assert run(case, scheduler)["reason"] == "submission budget exhausted"


def test_complete_does_not_submit_another_job(case):
    scheduler = Scheduler()
    run(case, scheduler)
    checkpoint(case, 1639)
    scheduler.job = dict(id="1234",state="COMPLETED",exit="0:0")
    assert run(case, scheduler)["phase"] == "complete"
    assert len(scheduler.commands) == 1


def test_changed_source_parameters_or_run_identity_stop(case):
    scheduler = Scheduler()
    with pytest.raises(RuntimeError,match="source changed"):
        run(case,scheduler,fingerprint=lambda source:"different")
    run(case,scheduler)
    checkpoint(case,19)
    scheduler.job = dict(id="1234",state="COMPLETED",exit="0:0")
    run(case,scheduler)
    checkpoint(case,39,run_id="changed-noise")
    with pytest.raises(ValueError,match="incompatible"):
        run(case,scheduler)


def test_changed_parameters_prevent_submission(case):
    scheduler = Scheduler()
    (case / "parameters.yaml").write_text("changed: true\n")
    with pytest.raises(RuntimeError,match="parameters changed"):
        run(case,scheduler)
    assert not scheduler.commands


def test_missing_shard_is_never_selected(case):
    path = checkpoint(case, 19)
    (path / "shard_0.h5").unlink()
    with pytest.raises(ValueError,match="missing/empty"):
        latest_progress(case / "data",nt=1640,members=32,dt=.05)


def test_checkpoint_time_can_start_at_nonzero_model_time(case):
    path = checkpoint(case,19)
    value = json.loads((path/'manifest.json').read_text())
    value['metadata']['model_time'] += 10
    (path/'manifest.json').write_text(json.dumps(value))
    assert latest_progress(case/'data',nt=1640,members=32,dt=.05,time_origin=10)['next_step']==20
    with pytest.raises(ValueError,match='incompatible'):
        latest_progress(case/'data',nt=1640,members=32,dt=.05)


def test_truncated_shard_is_never_selected(case):
    path = checkpoint(case,19)
    manifest = json.loads((path / "manifest.json").read_text())
    manifest["shards"][0]["nbytes"] = 9999
    (path / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match="size disagrees"):
        latest_progress(case / "data",nt=1640,members=32,dt=.05)


def test_exclusive_case_lock(case):
    with case_lock(case):
        with pytest.raises(FileExistsError):
            with case_lock(case):
                pytest.fail("second supervisor acquired a live lock")
    assert not (case / ".supervisor-lock").exists()


def test_source_pin_includes_dirty_code_and_ignores_run_products(tmp_path):
    import subprocess
    subprocess.run(["git","init","-q",str(tmp_path)],check=True)
    code = tmp_path / "model.py"
    code.write_text("x=1\n")
    subprocess.run(["git","-C",str(tmp_path),"add","model.py"],check=True)
    original = source_fingerprint(tmp_path)
    output = tmp_path / "pace_runs" / "run1"
    output.mkdir(parents=True)
    (output / "log.txt").write_text("progress")
    assert source_fingerprint(tmp_path) == original
    code.write_text("x=2\n")
    assert source_fingerprint(tmp_path) != original


def test_slurm_reconciles_long_names_without_padding(monkeypatch):
    import subprocess
    name = "icesee-cont-01234567890123456789"
    monkeypatch.setattr(subprocess,"check_output",lambda args,**kw:
        f"1234|{name:100s}|RUNNING\n" if args[0]=="squeue" else "")
    assert Slurm().find(name,"2026-11-24")["id"] == "1234"


def test_cli_prepare_resolves_auto_spatial_count_without_submitting(tmp_path,monkeypatch):
    import importlib.util
    from types import SimpleNamespace
    root = Path(__file__).resolve().parents[2]
    path = root / "applications/icepack_model/examples/idealized_pig/tools/continue_case.py"
    spec = importlib.util.spec_from_file_location("continuation_cli_test",path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module,"source_fingerprint",lambda source:"frozen")
    parameters = tmp_path / "input.yaml"
    parameters.write_text(yaml.safe_dump({"modeling-parameters":{"num_years":5,"timesteps_per_year":.05},
        "enkf-parameters":{"execution_mode":3,"Nens":24,"model_nprocs":0}}))
    args = SimpleNamespace(parameters=parameters,case=tmp_path/"case",data=tmp_path/"data",
        dependencies=tmp_path,nodes=1,ranks=8,memory="48G",walltime="02:00:00",
        run_seconds=5400,segment_cycles=20,max_submissions=3)
    result = module.prepare(args)
    assert result["status"] == "prepared; no job submitted"
    frozen = json.loads((args.case / "case.json").read_text())
    assert frozen["nt"] == 100 and frozen["members"] == 24
    assert not (args.case / "state.json").exists()
