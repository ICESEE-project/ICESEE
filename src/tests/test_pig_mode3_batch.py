"""Exercise the student batch entry point without starting MPI or Slurm."""
import json
import os
from pathlib import Path
import subprocess

import pytest

SCRIPT = (Path(__file__).resolve().parents[2] /
          "applications/icepack_model/examples/idealized_pig/run_idealized_1.sbatch")


def _run(tmp_path, *, counts=None, existing=False, step_exit=0):
    source = tmp_path / "scientific source" / "ICESEE"
    pig = source / "applications/icepack_model/examples/idealized_pig"
    (pig / "pace").mkdir(parents=True)
    (pig / "pace/resolve_environment.sh").write_text(
        'ICESEE_ROOT="$ICESEE_SOURCE_ROOT"\n'
        'PIG_DIR="$ICESEE_ROOT/applications/icepack_model/examples/idealized_pig"\n'
        'cd "$PIG_DIR"\n')
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "srun").write_text(
        '#!/usr/bin/env python3\nimport json, os, sys\n'
        'with open(os.environ["CAPTURE"], "w") as f: json.dump(sys.argv[1:], f)\n'
        'sys.exit(int(os.environ["STEP_EXIT"]))\n')
    (bin_dir / "git").write_text('#!/bin/sh\nprintf "test-source\\n"\n')
    for path in bin_dir.iterdir():
        path.chmod(0o755)
    data = pig / "_modelrun_datasets"
    if existing:
        data.mkdir()
        (data / "keep.txt").write_text("previous results")
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("SLURM_", "ICESEE_"))}
    env.update(ICESEE_SOURCE_ROOT=str(source), ICESEE_SPACK_ROOT=str(tmp_path / "deps"),
               SLURM_JOB_NUM_NODES="2", SLURM_NTASKS="16", SLURM_NTASKS_PER_NODE="8",
               SLURM_CPUS_PER_TASK="1", SLURM_JOB_ID="12345",
               PATH=str(bin_dir) + os.pathsep + env["PATH"],
               CAPTURE=str(tmp_path / "argv.json"), STEP_EXIT=str(step_exit))
    env.update(counts or {})
    result = subprocess.run(["bash", str(SCRIPT)], env=env, text=True,
                            capture_output=True, timeout=10)
    return result, pig, tmp_path / "argv.json"


def test_validated_allocation_has_16_ranks_and_32_members(tmp_path):
    result, pig, capture = _run(tmp_path)
    assert result.returncode == 0, result.stderr
    args = json.loads(capture.read_text())
    for expected in ("--nodes=2", "--ntasks=16", "--ntasks-per-node=8",
                     "--cpus-per-task=1", "--distribution=block:block",
                     "--execution_mode=3", "--mode3_experimental=True",
                     "--Nens=32", "--model_nprocs=2", "--save_steps=[2,3,4]",
                     "--data_path=" + str(pig / "_modelrun_datasets")):
        assert expected in args
    assert not any("overlap" in arg or "overcommit" in arg for arg in args)
    assert (pig / "pace_runs/mode3-12345/exit-status").read_text().strip() == "0"
    assert not (pig / "_modelrun_datasets.launch-lock").exists()


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


def test_srun_failure_propagates_and_releases_launch_lock(tmp_path):
    result, pig, _ = _run(tmp_path, step_exit=17)
    assert result.returncode == 17
    assert (pig / "pace_runs/mode3-12345/exit-status").read_text().strip() == "17"
    assert not (pig / "_modelrun_datasets.launch-lock").exists()
