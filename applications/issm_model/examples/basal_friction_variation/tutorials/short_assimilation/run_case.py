#!/usr/bin/env python3
"""Preflight and delegate a compact ISSM teaching case to the existing driver.

This wrapper never implements a forecast, analysis, or inverse solver. Check-only
mode loads YAML without importing the side-effecting application configuration.
"""
from __future__ import annotations

import argparse
import datetime as dt
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


CASE_DIR = Path(__file__).resolve().parent
EXAMPLE_DIR = CASE_DIR.parent.parent
VARIANTS = ("fixed", "enkf", "hybrid")


def repository_root() -> Path:
    for parent in CASE_DIR.parents:
        if (parent / "pyproject.toml").is_file() and (parent / "src").is_dir():
            return parent
    raise RuntimeError("Run this installed tutorial from an ICESEE source checkout")


def load_profile(repo: Path, variant: str) -> dict:
    # config_loader has no output-cleanup/application initialization side effects.
    spec = importlib.util.spec_from_file_location(
        "_icesee_tutorial_yaml", repo / "config" / "config_loader.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.load_yaml_to_dict(str(CASE_DIR / f"params-{variant}.yaml"))


def inspect_case(repo: Path, variant: str, run_id: str, ranks: int, mode: int,
                 launcher: str) -> tuple[dict, list[str]]:
    config = load_profile(repo, variant)
    modeling = config["modeling-parameters"]
    enkf = config["enkf-parameters"]
    data_path = f"_modelrun_tutorials/{run_id}/{variant}"
    example_name = f"{EXAMPLE_DIR.name}_tutorial_{run_id}_{variant}"
    issues = []
    required = [
        EXAMPLE_DIR / "run_da_issm.py",
        EXAMPLE_DIR / "Domain.exp",
        EXAMPLE_DIR / modeling["ParamFile"],
        EXAMPLE_DIR / modeling["reference_data_dir"] / modeling["reference_data"],
    ]
    if EXAMPLE_DIR.name == "ISMIP_Choi":
        required.append(EXAMPLE_DIR / "bed_kriging_results.h5")
    for path in required:
        if not path.is_file():
            issues.append(f"Missing required scientific/source input: {path}")
    if (EXAMPLE_DIR / data_path).exists():
        issues.append(f"Refusing existing assimilation output: {EXAMPLE_DIR / data_path}")
    issm_root = os.environ.get("ISSM_DIR")
    if not issm_root or not Path(issm_root).is_dir():
        issues.append("ISSM_DIR must select the approved ISSM installation")
    elif (Path(issm_root) / "examples" / example_name).exists():
        issues.append("Refusing an existing ISSM example workspace for this run ID")
    matlab = shutil.which("matlab")
    if matlab is None:
        issues.append("MATLAB executable is not on PATH after environment activation")
    if ranks > 1 and shutil.which(launcher) is None:
        issues.append(f"MPI launcher is not on PATH: {launcher}")
    if mode == 0 and ranks != 1:
        issues.append("Mode0 must be launched as one Python process in this tutorial")
    if not 1 <= ranks <= int(enkf["Nens"]):
        issues.append("Choose between one and four Python ranks for the four-member case")
    command = [
        sys.executable, str(EXAMPLE_DIR / "run_da_issm.py"),
        "-F", str(CASE_DIR / f"params-{variant}.yaml"),
        f"--execution_mode={mode}", "--Nens=4", "--model_nprocs=1",
        f"--example_name={example_name}", f"--data_path={data_path}",
    ]
    if ranks > 1:
        command = [launcher, "-n", str(ranks), *command]
    info = {
        "application": EXAMPLE_DIR.name, "variant": variant, "run_id": run_id,
        "configuration_before_application_setup": config,
        "command": command, "cwd": str(EXAMPLE_DIR),
        "output_path": str(EXAMPLE_DIR / data_path),
        "issm_example_name": example_name,
        "python": sys.executable, "python_ranks": ranks, "model_nprocs": 1,
        "matlab_executable": matlab, "ISSM_DIR": issm_root,
        "required_files": [str(p) for p in required],
        "expected_forecast_cycles": 10,
        "expected_state_columns": 11,
        "expected_observation_times": [1.0, 2.0],
        "expected_observation_state_columns": [5, 10],
        "backend_executed": False,
        "verification_status": "not performed by this wrapper",
        "preflight_issues": issues,
    }
    return info, issues


def verify_output_artifacts(info: dict) -> dict:
    """Check structural completion separately from the child process status.

    These checks cannot validate physical time, geometry or scientific skill.
    Import h5py only after execution; preflight stays free of backend imports.
    """
    output = Path(info["output_path"])
    issues = []
    ensemble_path = output / "icesee_ensemble_data.h5"
    observation_path = output / "synthetic_obs.h5"
    for path in (ensemble_path, observation_path):
        if not path.is_file():
            issues.append(f"Missing expected output: {path.name}")
    shapes = {}
    if ensemble_path.is_file():
        try:
            import h5py
            with h5py.File(ensemble_path, "r") as handle:
                for name in ("ensemble", "ensemble_mean"):
                    if name not in handle:
                        issues.append(f"Missing ensemble dataset: {name}")
                    else:
                        shapes[name] = list(handle[name].shape)
                ensemble = shapes.get("ensemble", [])
                mean = shapes.get("ensemble_mean", [])
                expected_columns = int(info["expected_state_columns"])
                if (len(ensemble) != 3 or ensemble[0] < 1 or ensemble[1] != 4
                        or ensemble[2] != expected_columns):
                    issues.append("Ensemble shape does not match four members and expected stored columns")
                if (len(mean) != 2 or not ensemble or mean != [ensemble[0], expected_columns]):
                    issues.append("Ensemble-mean shape does not match the stored ensemble")
        except Exception as error:
            issues.append(f"Cannot inspect ensemble HDF5: {type(error).__name__}: {error}")
    return {
        "status": "failed" if issues else "passed structural completion checks",
        "issues": issues, "ensemble_shapes": shapes,
        "limitation": "Does not validate physical time, initial state, finite values, geometry, inversion convergence or scientific skill",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", choices=(*VARIANTS, "all"), default="fixed")
    parser.add_argument("--mode", type=int, choices=(0, 1, 2), default=1)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--launcher", default="mpiexec")
    parser.add_argument("--run-id", default=dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--record-root", type=Path)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.run_id):
        parser.error("run-id may contain only letters, numbers, underscores and hyphens")
    repo = repository_root()
    variants = VARIANTS if args.variant == "all" else (args.variant,)
    cases = [inspect_case(repo, v, args.run_id, args.ranks, args.mode, args.launcher)
             for v in variants]
    print(json.dumps([info for info, _ in cases], indent=2), flush=True)
    if any(issues for _, issues in cases):
        return 2
    if args.check_only:
        print("Preflight paths/commands are available; no model, license or MPI test was executed.")
        return 0
    record_root = (args.record_root or CASE_DIR / "run_records").resolve() / args.run_id
    if record_root.exists():
        raise FileExistsError(f"Refusing existing record directory: {record_root}")
    record_root.mkdir(parents=True)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        (str(repo.parent), str(repo), env.get("PYTHONPATH", ""))
    )
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[key] = "1"
    for info, _ in cases:
        case_records = record_root / info["variant"]
        case_records.mkdir()
        (case_records / "requested-run.json").write_text(json.dumps(info, indent=2) + "\n")
        with (case_records / "stdout.log").open("w") as stdout, (case_records / "stderr.log").open("w") as stderr:
            result = subprocess.run(info["command"], cwd=EXAMPLE_DIR, env=env,
                                    stdout=stdout, stderr=stderr, check=False)
        info["backend_executed"] = True
        info["simulation_exit_code"] = result.returncode
        info["output_verification"] = verify_output_artifacts(info)
        tutorial_exit_code = result.returncode or (1 if info["output_verification"]["issues"] else 0)
        info["tutorial_exit_code"] = tutorial_exit_code
        (case_records / "execution-status.json").write_text(json.dumps(info, indent=2) + "\n")
        print(f"{info['variant']}: process exit {result.returncode}, "
              f"{info['output_verification']['status']}; inspect {case_records}", flush=True)
        if tutorial_exit_code:
            return tutorial_exit_code
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
