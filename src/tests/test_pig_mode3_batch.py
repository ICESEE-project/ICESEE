"""Exercise the PIG batch entry point without starting MPI or Slurm."""
import json
import os
from pathlib import Path
import subprocess
import shutil
import sys
import yaml

import pytest

SCRIPT = (Path(__file__).resolve().parents[2] /
          "applications/icepack_model/examples/idealized_pig/run_idealized_1.sbatch")


def _run(tmp_path, *, counts=None, existing=False, step_exit=0, finalizer_exit=0, discover=False, local_params=False, args=(), config_updates=None):
    lab = tmp_path / "lab space"
    source = lab / "ICESEE-Spack-source" / "ICESEE"
    pig = source / "applications/icepack_model/examples/idealized_pig"
    (pig / "pace").mkdir(parents=True)
    (pig / "pace/resolve_environment.sh").write_text(
        'ICESEE_ROOT="$ICESEE_SOURCE_ROOT"\n'
        'SPACK_REPO="$ICESEE_SPACK_ROOT"\n'
        'PIG_DIR="$ICESEE_ROOT/applications/icepack_model/examples/idealized_pig"\n'
        'cd "$PIG_DIR"\n')
    exporter = source / "src/parallelization/distributed_history_export.py"
    exporter.parent.mkdir(parents=True)
    exporter.write_text(
        'import json, os, sys\n'
        'with open(os.environ["FINALIZER_CAPTURE"], "w") as f: json.dump(sys.argv[1:], f)\n'
        'sys.exit(int(os.environ["FINALIZER_EXIT"]))\n')
    cfg = {"science": "tracked", "modeling-parameters": {"num_years": 5, "timesteps_per_year": 0.05},
           "enkf-parameters": {"Nens": 24, "model_nprocs": 2, "execution_mode": 3,
                               "checkpoint_every": 1, "checkpoint_keep_last": 0,
                               "mode3_export_ensemble_history": True}}
    for section, values in (config_updates or {}).items():
        cfg[section].update(values)
    (pig / "params.yaml").write_text(yaml.safe_dump(cfg))
    (source / "config").mkdir()
    shutil.copy2(SCRIPT.parents[4] / "config/config_loader.py", source / "config/config_loader.py")
    deps = lab / "ICESEE-Spack"
    (deps / "scripts").mkdir(parents=True)
    (deps / "scripts/activate.sh").write_text(":\n")
    (deps / ".spack-env/icesee").mkdir(parents=True)
    if local_params:
        config = lab / "ICESEE-run-configs/idealized_pig.yaml"
        config.parent.mkdir()
        config.write_text(yaml.safe_dump(dict(cfg, science="retained")))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "srun").write_text(
        '#!/usr/bin/env python3\nimport json, os, sys\n'
        'with open(os.environ["CAPTURE"], "w") as f: json.dump(sys.argv[1:], f)\n'
        'sys.exit(int(os.environ["STEP_EXIT"]))\n')
    (bin_dir / "git").write_text('#!/bin/sh\nprintf "test-source\\n"\n')
    for path in bin_dir.iterdir():
        path.chmod(0o755)
    (bin_dir / "python").symlink_to(sys.executable)
    data = pig / "_modelrun_datasets"
    if existing:
        data.mkdir()
        (data / "keep.txt").write_text("previous results")
        if cfg['enkf-parameters'].get('resume_from_checkpoint'):
            (data/'_mode3_state_history').mkdir()
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("SLURM_", "ICESEE_"))}
    env.update(ICESEE_SOURCE_ROOT=str(source), ICESEE_SPACK_ROOT=str(deps),
               SLURM_JOB_NUM_NODES="2", SLURM_NTASKS="16", SLURM_NTASKS_PER_NODE="8",
               SLURM_CPUS_PER_TASK="1", SLURM_JOB_ID="12345",
               PATH=str(bin_dir) + os.pathsep + env["PATH"],
               CAPTURE=str(tmp_path / "argv.json"), STEP_EXIT=str(step_exit),
               FINALIZER_EXIT=str(finalizer_exit),
               FINALIZER_CAPTURE=str(tmp_path / "finalizer-argv.json"))
    if discover:
        env.pop("ICESEE_SOURCE_ROOT")
        env.pop("ICESEE_SPACK_ROOT")
        env["SLURM_SUBMIT_DIR"] = str(pig)
    for key, value in (counts or {}).items():
        if value is None:
            env.pop(key, None)
        else:
            env[key] = value
    # GNU time is specific to PACE; mock only its wrapper on other hosts.
    launcher = tmp_path / "launcher.sbatch"
    launcher.write_text(SCRIPT.read_text().replace('/usr/bin/time -v', 'env'))
    result = subprocess.run(["bash", str(launcher), *args], env=env, text=True,
                            capture_output=True, timeout=10)
    return result, pig, tmp_path / "argv.json"


def test_validated_allocation_has_16_ranks_and_32_members(tmp_path):
    result, pig, capture = _run(tmp_path, config_updates={"enkf-parameters": {"Nens": 32}})
    assert result.returncode == 0, result.stderr
    args = json.loads(capture.read_text())
    for expected in ("--nodes=2", "--ntasks=16", "--ntasks-per-node=8",
                     "--cpus-per-task=1", "--distribution=block:block",
                     "--execution_mode=3",

                     "--Nens=32", "--model_nprocs=2",
                     "--data_path=" + str(pig / "_modelrun_datasets")):
        assert expected in args
    assert not any("mode3_experimental" in arg for arg in args)
    assert not any("overlap" in arg or "overcommit" in arg for arg in args)
    assert (pig / "pace_runs/mode3-12345/exit-status").read_text().strip() == "0"
    assert not (pig / "_modelrun_datasets.launch-lock").exists()


@pytest.mark.parametrize('nodes,ranks,per_node', [('1','1','1'), ('2','40','20')])
def test_missing_tasks_per_node_is_derived_from_allocation(tmp_path, nodes, ranks, per_node):
    result, _, capture = _run(tmp_path, counts={
        'SLURM_JOB_NUM_NODES':nodes, 'SLURM_NTASKS':ranks,
        'SLURM_NTASKS_PER_NODE':None,
    }, config_updates={'enkf-parameters':{'Nens':40, 'model_nprocs':1}})
    assert result.returncode == 0, result.stderr
    args = json.loads(capture.read_text())
    assert '--ntasks='+ranks in args
    assert '--ntasks-per-node='+per_node in args


def test_missing_tasks_per_node_rejects_nonuniform_allocation(tmp_path):
    result, _, capture = _run(tmp_path, counts={
        'SLURM_JOB_NUM_NODES':'2', 'SLURM_NTASKS':'3',
        'SLURM_NTASKS_PER_NODE':None,
    }, config_updates={'enkf-parameters':{'model_nprocs':1}})
    assert result.returncode == 2 and not capture.exists()
    assert 'uniform tasks/node' in result.stderr


def test_case_specific_log_path_keeps_the_original_allocation(tmp_path):
    logs = tmp_path/'paired case/logs'
    result, _, capture = _run(tmp_path, counts={"ICESEE_LOG_PATH": str(logs)})
    assert result.returncode == 0, result.stderr
    assert (logs/'params.yaml').is_file() and (logs/'exit-status').read_text().strip() == '0'
    args = json.loads(capture.read_text())
    assert '--nodes=2' in args and '--ntasks=16' in args


def test_log_path_inside_recreated_data_is_refused_before_launch(tmp_path):
    logs = tmp_path/'lab space/ICESEE-Spack-source/ICESEE/applications/icepack_model/examples/idealized_pig/_modelrun_datasets/logs'
    result, _, capture = _run(tmp_path, counts={"ICESEE_LOG_PATH": str(logs)})
    assert result.returncode != 0 and not capture.exists()
    assert 'Logs must be outside' in result.stderr


def test_reference_24_member_eight_rank_allocation_uses_mode_selector_only(tmp_path):
    result, _, capture = _run(tmp_path, counts={
        "SLURM_JOB_NUM_NODES": "1", "SLURM_NTASKS": "8",
        "SLURM_NTASKS_PER_NODE": "8", "ICESEE_NENS": "24",
    })
    assert result.returncode == 0, result.stderr
    args = json.loads(capture.read_text())
    for expected in ("--nodes=1", "--ntasks=8", "--Nens=24",
                     "--model_nprocs=2", "--execution_mode=3"):
        assert expected in args
    assert not any("mode3_experimental" in arg for arg in args)


@pytest.mark.parametrize("counts", [
    {"SLURM_NTASKS": "40"},
    {"SLURM_NTASKS": "14", "SLURM_NTASKS_PER_NODE": "7"},
])
def test_inconsistent_or_cross_node_groups_fail_before_launch(tmp_path, counts):
    result, _, capture = _run(tmp_path, counts=counts)
    assert result.returncode == 2
    assert not capture.exists()


def test_existing_scientific_results_are_preserved(tmp_path):
    result, pig, capture = _run(tmp_path, existing=True)
    assert result.returncode == 2
    assert not capture.exists()
    assert (pig / "_modelrun_datasets/keep.txt").read_text() == "previous results"


def test_explicit_native_resume_preserves_existing_data(tmp_path):
    result, pig, capture = _run(tmp_path, existing=True,
        config_updates={'enkf-parameters':{'resume_from_checkpoint':True}})
    assert result.returncode == 0, result.stderr
    assert capture.exists() and (pig/'_modelrun_datasets/keep.txt').read_text()=='previous results'


def test_resume_without_history_is_refused(tmp_path):
    result, _, capture = _run(tmp_path,
        config_updates={'enkf-parameters':{'resume_from_checkpoint':True}})
    assert result.returncode == 2 and not capture.exists()


def test_disabled_history_does_not_export_after_failure(tmp_path):
    result, _, capture = _run(tmp_path, step_exit=17,
        config_updates={'enkf-parameters':{'mode3_export_ensemble_history':False}})
    assert result.returncode == 17 and capture.exists()
    assert not (tmp_path/'finalizer-argv.json').exists()


def test_srun_failure_propagates_and_releases_launch_lock(tmp_path):
    result, pig, _ = _run(tmp_path, step_exit=17)
    assert result.returncode == 17
    assert (pig / "pace_runs/mode3-12345/exit-status").read_text().strip() == "17"
    assert not (pig / "_modelrun_datasets.launch-lock").exists()


def test_direct_submission_discovers_source_dependencies_and_saved_parameters(tmp_path):
    result, pig, capture = _run(tmp_path, discover=True, local_params=True, counts={
        "SLURM_JOB_NUM_NODES": "1", "SLURM_NTASKS": "8",
    })
    assert result.returncode == 0, result.stderr
    args = json.loads(capture.read_text())
    config = pig.parents[5] / "ICESEE-run-configs/idealized_pig.yaml"
    assert "--force-params=" + str(pig / "pace_runs/mode3-12345/params.yaml") in args
    assert "--Nens=24" in args
    assert "--ntasks=8" in args
    assert yaml.safe_load((pig / "pace_runs/mode3-12345/params.yaml").read_text())["science"] == "tracked"
    assert str(pig.parents[5] / "ICESEE-Spack") in result.stdout


def test_explicit_parameter_file_takes_precedence(tmp_path):
    config = tmp_path / "explicit.yaml"
    config.write_text(yaml.safe_dump({"science": "explicit",
        "modeling-parameters": {"num_years": 5, "timesteps_per_year": 0.05},
        "enkf-parameters": {"Nens": 24, "model_nprocs": 2, "execution_mode": 3}}))
    result, pig, capture = _run(tmp_path, local_params=True, args=["-F", str(config)])
    assert result.returncode == 0, result.stderr
    assert "--force-params=" + str(pig / "pace_runs/mode3-12345/params.yaml") in json.loads(capture.read_text())
    assert yaml.safe_load((pig / "pace_runs/mode3-12345/params.yaml").read_text())["science"] == "explicit"


def test_missing_explicit_parameters_fail_before_launch(tmp_path):
    result, _, capture = _run(tmp_path, args=["-F", str(tmp_path / "absent.yaml")])
    assert result.returncode == 2
    assert not capture.exists()


@pytest.mark.parametrize('dt,columns', [('0.1', 51), ('0.025', 201)])
def test_chosen_dt_reaches_driver_and_history_column_report(tmp_path, dt, columns):
    result, _, capture = _run(tmp_path, config_updates={'modeling-parameters': {'timesteps_per_year': dt}})
    assert result.returncode == 0, result.stderr
    assert '--dt=' + str(float(dt)) in json.loads(capture.read_text())
    assert f'history_columns={columns}' in result.stdout


@pytest.mark.parametrize('dt', ['0', '-0.1', 'nan', 'inf', 'invalid'])
def test_invalid_dt_is_rejected_before_launch(tmp_path, dt):
    result, _, capture = _run(tmp_path, config_updates={'modeling-parameters': {'timesteps_per_year': dt}})
    assert result.returncode == 2
    assert not capture.exists()


@pytest.mark.parametrize("finalizer_exit", [0, 4])
def test_failed_step_recovers_history_without_hiding_simulation_failure(tmp_path, finalizer_exit):
    result, pig, _ = _run(tmp_path, step_exit=17, finalizer_exit=finalizer_exit)
    assert result.returncode == 17
    log = pig / "pace_runs/mode3-12345"
    assert (log / "exit-status").read_text().strip() == "17"
    assert (log / "finalizer-exit-status").read_text().strip() == str(finalizer_exit)
    args = json.loads((tmp_path / "finalizer-argv.json").read_text())
    assert args == [str(pig / "_modelrun_datasets"), "--allow-partial", "--simulation-exit-code", "17"]
    assert (log / "finalizer.log").is_file()


def test_successful_step_does_not_run_failure_recovery(tmp_path):
    result, pig, _ = _run(tmp_path)
    assert result.returncode == 0
    assert not (tmp_path / "finalizer-argv.json").exists()
    assert not (pig / "pace_runs/mode3-12345/finalizer-exit-status").exists()


def test_old_environment_overrides_do_not_replace_selected_yaml(tmp_path):
    result, pig, capture = _run(tmp_path, local_params=True, counts={
        "ICESEE_PARAMS_FILE": "/not/the/config.yaml", "ICESEE_NENS": "99", "ICESEE_DT": "0.9"})
    assert result.returncode == 0, result.stderr
    args = json.loads(capture.read_text())
    assert "--Nens=24" in args and "--dt=0.05" in args
    assert "Ignoring ICESEE_PARAMS_FILE" in result.stderr
    assert yaml.safe_load((pig / "pace_runs/mode3-12345/params.yaml").read_text())["science"] == "tracked"


@pytest.mark.parametrize("args", [["-F"], ["--unknown"], ["--Nens=99"]])
def test_bad_arguments_fail_before_launch(tmp_path, args):
    result, _, capture = _run(tmp_path, args=args)
    assert result.returncode == 2 and not capture.exists()


@pytest.mark.parametrize("mode", [1, 2])
def test_yaml_selects_execution_mode_without_a_hidden_mode3_override(tmp_path, mode):
    result, pig, capture = _run(tmp_path, config_updates={"enkf-parameters": {"execution_mode": mode}})
    assert result.returncode == 0, result.stderr
    assert f"--execution_mode={mode}" in json.loads(capture.read_text())
    assert (pig / f"pace_runs/mode{mode}-12345/params.yaml").is_file()


def test_yaml_observation_and_profile_changes_are_preserved(tmp_path):
    result, pig, capture = _run(tmp_path, config_updates={
        "enkf-parameters": {"obs_start_time": 0.25, "obs_max_time": 0.75, "freq_obs": 0.25},
        "modeling-parameters": {"save_steps": [0.25, 0.75], "num_years": 1}})
    assert result.returncode == 0, result.stderr
    cfg = yaml.safe_load((pig / "pace_runs/mode3-12345/params.yaml").read_text())
    assert cfg["enkf-parameters"]["obs_start_time"] == 0.25
    assert cfg["enkf-parameters"]["obs_max_time"] == 0.75
    assert cfg["modeling-parameters"]["save_steps"] == [0.25, 0.75]
    assert cfg["modeling-parameters"]["nt"] == 20
    assert not any(a.startswith(("--obs_", "--freq_obs", "--save_steps")) for a in json.loads(capture.read_text()))


def test_forced_yaml_inheritance_is_resolved_in_the_saved_snapshot(tmp_path):
    base = tmp_path / "base.yaml"
    base.write_text(yaml.safe_dump({"modeling-parameters": {"num_years": 5, "timesteps_per_year": 0.05},
        "enkf-parameters": {"Nens": 24, "model_nprocs": 2, "execution_mode": 3}}))
    selected = tmp_path / "variant.yaml"
    selected.write_text("extends: base.yaml\nenkf-parameters:\n  Nens: 32\n")
    result, pig, capture = _run(tmp_path, args=["-F", str(selected)])
    assert result.returncode == 0, result.stderr
    cfg = yaml.safe_load((pig / "pace_runs/mode3-12345/params.yaml").read_text())
    assert "extends" not in cfg and cfg["enkf-parameters"]["Nens"] == 32
    assert "--Nens=32" in json.loads(capture.read_text())
