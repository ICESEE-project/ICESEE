"""Prepare and supervise bounded PIG Mode-3 continuation jobs."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

PIG = Path(__file__).resolve().parents[1]
SOURCE = PIG.parents[3]
sys.path[:0] = [str(SOURCE.parent), str(SOURCE / "config")]

from ICESEE.src.parallelization.checkpoint_continuation import source_fingerprint, tick
from ICESEE.src.parallelization.run_control import atomic_json, atomic_text


def prepare(args):
    import yaml
    from config_loader import load_yaml_to_dict
    cfg = load_yaml_to_dict(str(args.parameters.resolve()))
    model, options = cfg["modeling-parameters"], cfg["enkf-parameters"]
    if int(options.get("execution_mode", -1)) != 3:
        raise ValueError("this continuation launcher requires native Mode 3")
    if int(options.get("checkpoint_every", 1)) != 1:
        raise ValueError("continuation requires checkpoint_every=1")
    if options.get("run_stop_time") is not None or options.get("reference_only"):
        raise ValueError("remove diagnostic run_stop_time/reference_only from the base configuration")
    dt, horizon = float(model["timesteps_per_year"]), float(model["num_years"])
    if not all(math.isfinite(v) and v > 0 for v in (dt, horizon)):
        raise ValueError("duration and dt must be finite and positive")
    nt = round(horizon / dt)
    if not math.isclose(nt * dt, horizon, rel_tol=0, abs_tol=1e-8):
        raise ValueError("duration must be an integral number of dt intervals")
    if min(args.segment_cycles, args.max_submissions, args.run_seconds, args.ranks, args.nodes) < 1:
        raise ValueError("segment, allocation and submission budgets must be positive")
    spatial = int(options.get("model_nprocs", 0)) or 1
    if spatial < 1 or args.ranks % args.nodes or (args.ranks // args.nodes) % spatial:
        raise ValueError("each node must contain complete model groups")
    parts = args.walltime.split(":")
    if len(parts) != 3 or any(not p.isdigit() for p in parts):
        raise ValueError("walltime must be HH:MM:SS")
    hours, minutes, seconds = map(int, parts)
    if minutes >= 60 or seconds >= 60 or hours * 3600 + minutes * 60 + seconds <= args.run_seconds:
        raise ValueError("allocation walltime must exceed run_seconds; reserve a full cycle plus output time")
    source = SOURCE.resolve()
    case = args.case.resolve()
    data = args.data.resolve()
    if source == case or source in case.parents and "pace_runs" not in case.relative_to(source).parts:
        raise ValueError("put case records outside source or beneath pace_runs")
    if case.exists():
        raise ValueError("case directory already exists; use tick/watch to continue it")
    case.mkdir(parents=True)
    parameters = case / "parameters.yaml"
    atomic_text(parameters, yaml.safe_dump(cfg, sort_keys=False))
    spec = dict(source=str(source), source_sha256=source_fingerprint(source),
        parameters=str(parameters), parameters_sha256=hashlib.sha256(parameters.read_bytes()).hexdigest(),
        data=str(data), nt=nt, members=int(options["Nens"]), dt=dt,
        segment_cycles=args.segment_cycles, run_seconds=args.run_seconds,
        max_submissions=args.max_submissions, launcher=str(PIG / "run_idealized_1.sbatch"),
        sbatch_args=["--nodes=" + str(args.nodes), "--ntasks=" + str(args.ranks),
            "--ntasks-per-node=" + str(args.ranks // args.nodes), "--mem=" + args.memory,
            "--time=" + args.walltime,
            "--chdir=" + str(PIG),
            "--export=ALL,ICESEE_SOURCE_ROOT=" + str(source) + ",ICESEE_SPACK_ROOT=" + str(args.dependencies.resolve())
            + ",ICESEE_DATA_PATH=" + str(data)])
    if any("," in str(v) or "\n" in str(v) for v in (source, data, args.dependencies.resolve())):
        raise ValueError("source/dependency/data paths cannot contain commas or newlines")
    atomic_json(case / "case.json", spec)
    return dict(case=str(case), status="prepared; no job submitted", max_submissions=args.max_submissions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    create = commands.add_parser("prepare", help="freeze a case without submitting")
    create.add_argument("case", type=Path)
    create.add_argument("-F", "--parameters", type=Path, default=PIG / "params.yaml")
    create.add_argument("--data", type=Path, required=True)
    create.add_argument("--dependencies", type=Path, required=True)
    create.add_argument("--segment-cycles", type=int, required=True)
    create.add_argument("--max-submissions", type=int, required=True)
    create.add_argument("--run-seconds", type=int, required=True)
    create.add_argument("--walltime", required=True)
    create.add_argument("--nodes", type=int, required=True)
    create.add_argument("--ranks", type=int, required=True)
    create.add_argument("--memory", default="48G")
    for action in ("tick", "watch"):
        sub = commands.add_parser(action)
        sub.add_argument("case", type=Path)
    args = parser.parse_args()
    if args.action == "prepare":
        print(json.dumps(prepare(args), indent=2))
        return
    while True:
        state = tick(args.case)
        print(json.dumps(state, sort_keys=True), flush=True)
        if state["phase"] == "stopped" or "attention" in state:
            raise SystemExit(2)
        if args.action != "watch" or state["phase"] in {"complete", "stopped"} or "attention" in state:
            return
        time.sleep(60)


if __name__ == "__main__":
    main()
