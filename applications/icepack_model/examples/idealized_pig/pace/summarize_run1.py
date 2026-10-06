# ==============================================================================
# @des: Post-run instrumentation summary for PACE RUN 1 (Idealized PIG
# Mode-3 calibration, 2026-09-28). Parses the run's own stdout log (which
# already contains ICESEE's built-in timing table -- src/utils/tools.py's
# display_timing_verbose -- plus the member-store stats line and the
# per-step forecast timing this reconciliation added to mode3_runner.py)
# and combines it with filesystem measurements (store file count/size,
# checkpoint size) and, where available, peak RSS, into one JSON summary.
#
# Does not re-run or re-instrument anything -- purely a log/filesystem
# parser, run once after RUN 1 finishes.
# ==============================================================================
from __future__ import annotations

import argparse
import ast
import json
import os
import re
from pathlib import Path


_ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")
# display_timing_verbose (src/utils/tools.py) prints one line per metric as
# "<ansi>[0m<box-char><ansi> NAME   VALUE <ansi><box-char><ansi>" -- name and
# value are separated by whitespace only (no internal column separator), and
# the whole line is bounded by a single U+2551 box-drawing char on each side.
_TIME_ROW_RE_BOX = re.compile(
    r"║\s*([A-Za-z/\(\)\sΣ]+?)\s{2,}(\d{2}:\d{2}:\d{2}:\d{2}\.\d{3})\s*║"
)
_STORE_STATS_RE = re.compile(r"\[ICESEE\] mode3 member-store stats \(world_rank=(\d+)\): (\{.*\})")
_STEP_RE = re.compile(
    r"\[ICESEE\] mode3 forecast step (\d+): wall=([\d.]+)s did_analysis=(True|False)"
)


def _dhms_to_seconds(s: str) -> float:
    days, hours, minutes, rest = s.split(":")
    seconds, ms = rest.split(".")
    return (
        int(days) * 86400 + int(hours) * 3600 + int(minutes) * 60
        + int(seconds) + int(ms) / 1000.0
    )


def parse_log(log_path: Path) -> dict:
    raw = log_path.read_text(errors="replace")
    text = _ANSI_RE.sub("", raw)

    timing = {}
    for name, value in _TIME_ROW_RE_BOX.findall(text):
        name = name.strip()
        if name:
            timing[name] = _dhms_to_seconds(value)

    store_stats_by_rank = {}
    for rank, blob in _STORE_STATS_RE.findall(text):
        try:
            store_stats_by_rank[int(rank)] = ast.literal_eval(blob)
        except (ValueError, SyntaxError):
            store_stats_by_rank[int(rank)] = {"_unparsed": blob}

    steps = []
    for k, wall, did_analysis in _STEP_RE.findall(text):
        steps.append({"step": int(k), "wall_seconds": float(wall), "did_analysis": did_analysis == "True"})

    first_step = steps[0]["wall_seconds"] if steps else None
    rest_steps = [s["wall_seconds"] for s in steps[1:]] if len(steps) > 1 else []
    rest_mean = sum(rest_steps) / len(rest_steps) if rest_steps else None

    return {
        "timing_table_seconds": timing,
        "store_stats_by_rank": store_stats_by_rank,
        "per_step": steps,
        "first_step_wall_seconds": first_step,
        "steady_state_mean_step_wall_seconds": rest_mean,
        "steady_state_step_count": len(rest_steps),
    }


def _dir_stats(path: Path) -> dict:
    if not path.is_dir():
        return {"exists": False}
    total_bytes = 0
    file_count = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            fp = Path(root) / name
            try:
                total_bytes += fp.stat().st_size
                file_count += 1
            except OSError:
                pass
    return {"exists": True, "file_count": file_count, "total_bytes": total_bytes}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--summary-dir", default=None, help="where run1_summary.json is written (default: --output-dir)")
    parser.add_argument("--wall-seconds", type=float, default=None, help="wall time of the DA run itself")
    parser.add_argument("--job-id", default="unknown")
    parser.add_argument("--nens", type=int, required=True)
    parser.add_argument("--model-nprocs", type=int, required=True)
    parser.add_argument("--world-size", type=int, required=True)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    log_parse = parse_log(Path(args.log_file))
    # ICESEE writes the end-of-run performance summary to
    # <data_path>/performance.json from the same data as the printed report;
    # its "legacy" block is the historical timing table (null = not measured).
    performance_json = out_dir / "performance.json"
    if performance_json.is_file():
        performance = json.loads(performance_json.read_text())
        log_parse["timing_table_seconds"] = performance.get("legacy", {})
        log_parse["performance_json"] = str(performance_json)

    summary = {
        "job_id": args.job_id,
        "config": {
            "Nens": args.nens,
            "model_nprocs": args.model_nprocs,
            "world_size": args.world_size,
            "ensemble_groups": args.world_size // args.model_nprocs,
            "rounds": -(-args.nens // (args.world_size // args.model_nprocs)),  # ceil division
        },
        "da_wall_seconds": args.wall_seconds,
        "data_path": str(out_dir),
        "timing": log_parse,
        "member_store_dir": _dir_stats(out_dir / "_mode3_member_store"),
        "checkpoint_dir": _dir_stats(out_dir / "_mode3_state_history"),
    }

    out_path = Path(args.summary_dir or out_dir) / "run1_summary.json"
    with out_path.open("w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"[summarize_run1] wrote {out_path}")
    print(f"[summarize_run1] first forecast step: {log_parse['first_step_wall_seconds']}s")
    print(f"[summarize_run1] steady-state mean ({log_parse['steady_state_step_count']} steps): {log_parse['steady_state_mean_step_wall_seconds']}s")
    print(f"[summarize_run1] timing table: {json.dumps(log_parse['timing_table_seconds'], indent=2)}")


if __name__ == "__main__":
    main()
