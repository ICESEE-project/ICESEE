"""Compare completed PIG Mode-3 checkpoints and physical flowline profiles.

Report differences without imposing a scientific acceptance tolerance.
"""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def inventory(data):
    paths = sorted((data / "_mode3_state_history/steps").glob("checkpoint_*/manifest.json"))
    completed = [(p, json.loads(p.read_text())) for p in paths]
    completed = [(p, m) for p, m in completed if m.get("complete")]
    if not completed:
        raise ValueError(f"No committed checkpoint under {data}")
    path, manifest = max(completed, key=lambda item: item[1]["timestep"])
    if manifest["metadata"]["next_step"] != manifest["metadata"]["nt"]:
        raise ValueError(f"Run is not complete: {data}")
    if manifest["timestep"] != manifest["metadata"]["nt"]-1:
        raise ValueError(f"Final step does not match run length: {data}")
    entries = {}
    coverage = defaultdict(list)
    sizes = {block["name"]: block["global_size"] for block in manifest["blocks"]}
    for shard in manifest["shards"]:
        shard_path = (path.parent / shard["file"]).resolve()
        if shard_path.parent != path.parent.resolve():
            raise ValueError("Shard path leaves checkpoint directory")
        for member in shard["member_ids"]:
            for block in shard["blocks"]:
                key = (member, block["name"], block["owned_start"], block["owned_stop"])
                if key in entries:
                    raise ValueError(f"Duplicate ownership: {key}")
                entries[key] = (shard_path, shard.get("layout_fingerprint"))
                coverage[(member, block["name"])].append((block["owned_start"], block["owned_stop"]))
    expected = {(member, field) for member in range(manifest["number_of_members"]) for field in sizes}
    if coverage.keys() != expected:
        raise ValueError("Missing or unexpected member/field ownership")
    for (member, field), intervals in coverage.items():
        end = 0
        for start, stop in sorted(intervals):
            if start != end or stop < start:
                raise ValueError(f"Gap/overlap in ownership: {member}/{field}")
            end = stop
        if end != sizes[field]:
            raise ValueError(f"Incomplete global field: {member}/{field}")
    return manifest, entries


def metrics(a, b):
    delta = b - a
    return dict(values=a.size, max_abs=float(np.max(np.abs(delta))),
                sum_squared=float(np.sum(delta * delta)), reference_squared=float(np.sum(a * a)))


def finish(stats):
    return {"values": stats["values"], "max_abs": stats["max_abs"],
            "rmse": float(np.sqrt(stats["sum_squared"] / stats["values"])),
            "relative_l2": float(np.sqrt(stats["sum_squared"] / stats["reference_squared"]))
            if stats["reference_squared"] else (0.0 if not stats["sum_squared"] else None)}


def compare_fields(reference, candidate):
    rm, ri = inventory(reference)
    cm, ci = inventory(candidate)
    for key in ("number_of_members", "timestep", "global_size", "blocks", "layout_id", "run_id"):
        if rm[key] != cm[key]:
            raise ValueError(f"Incompatible checkpoint {key}")
    for key in ("next_step", "nt", "model_time", "analysis_events_completed"):
        if rm["metadata"][key] != cm["metadata"][key]:
            raise ValueError(f"Incompatible run metadata {key}")
    if ri.keys() != ci.keys():
        raise ValueError("Owned intervals/member scheduling differ; physical mapping is required")
    stats = defaultdict(lambda: dict(values=0, max_abs=0., sum_squared=0., reference_squared=0.))
    for key in sorted(ri):
        member, field, start, stop = key
        rp, rf = ri[key]; cp, cf = ci[key]
        if not rf or rf != cf:
            raise ValueError(f"Different mesh numbering for {key}; do not compare raw rows")
        with h5py.File(rp, "r") as r, h5py.File(cp, "r") as c:
            rd, cd = r[f"member_{member}/{field}"], c[f"member_{member}/{field}"]
            if rd.shape != cd.shape or rd.shape != (stop-start,):
                raise ValueError(f"Invalid field shape: {key}")
            for pos in range(0, stop-start, 8192):
                a, b = rd[pos:pos+8192], cd[pos:pos+8192]
                if not (np.isfinite(a).all() and np.isfinite(b).all()):
                    raise ValueError(f"Nonfinite checkpoint field: {key}")
                chunk = metrics(a, b); target = stats[field]
                target["max_abs"] = max(target["max_abs"], chunk.pop("max_abs"))
                for name, value in chunk.items():
                    target[name] += value
    return {k: finish(v) for k, v in stats.items()}, rm["metadata"]


def read_profile(path):
    with h5py.File(path, "r") as f:
        result = {name: f[name][:] for name in ("h_profiles", "s_profiles", "distances", "bed_values")}
        for name in ("written_steps", "profile_times"):
            result[name] = f[name][:] if name in f else np.asarray(f.attrs[name])
    for name in ("distances", "bed_values", "profile_times", "written_steps"):
        if not np.isfinite(result[name]).all():
            raise ValueError(f"Nonfinite profile metadata: {path}/{name}")
    written = result["written_steps"] >= 0
    for name in ("h_profiles", "s_profiles"):
        if result[name].ndim != 2 or result[name].shape[1] != len(written):
            raise ValueError(f"Invalid profile shape: {path}/{name}")
        if not np.isfinite(result[name][:, written]).all():
            raise ValueError(f"Nonfinite written profile: {path}/{name}")
    return result


def compare_profiles(reference, candidate, output, members=None, labels=("Reference", "Candidate")):
    rp = {p.name: p for p in reference.glob("hs*.h5")}
    cp = {p.name: p for p in candidate.glob("hs*.h5")}
    if rp.keys() != cp.keys() or "hs_profiles_true.h5" not in rp:
        raise ValueError("Truth/member profile file sets differ")
    if members is not None:
        expected = {"hs_profiles_true.h5"} | {f"hs_ensemble_profiles{m}.h5" for m in range(members)}
        if rp.keys() != expected:
            raise ValueError("Profile files do not cover the checkpoint member IDs")
    reports = {}; ensembles = [[], []]; truths = []
    for name in sorted(rp):
        a, b = read_profile(rp[name]), read_profile(cp[name])
        for key in ("distances", "bed_values", "written_steps", "profile_times"):
            np.testing.assert_array_equal(a[key], b[key], err_msg=f"Profile geometry/time mismatch: {name}/{key}")
        valid = np.asarray(a["written_steps"]) >= 0
        if not valid.any():
            raise ValueError(f"No written profile columns: {name}")
        reports[name] = {}
        for key in ("h_profiles", "s_profiles"):
            if a[key].shape != b[key].shape or a[key].shape[1] != len(valid):
                raise ValueError(f"Profile shape mismatch: {name}/{key}")
            reports[name][key] = finish(metrics(a[key][:, valid], b[key][:, valid]))
            reports[name][key]["by_time"] = [
                dict(time=float(a["profile_times"][i]), step=int(a["written_steps"][i]),
                     **finish(metrics(a[key][:, i], b[key][:, i])))
                for i in np.flatnonzero(valid)
            ]
        if name == "hs_profiles_true.h5":
            truths = [a, b]
        else:
            ensembles[0].append(a); ensembles[1].append(b)
    times = truths[0]["profile_times"]; valid = truths[0]["written_steps"] >= 0
    x = np.ravel(truths[0]["distances"]) / 1000.
    for field, label in (("h_profiles", "Thickness (m)"), ("s_profiles", "Surface elevation (m)")):
        stacks = [np.stack([p[field] for p in members]) for members in ensembles]
        fig, axes = plt.subplots(2, int(valid.sum()), figsize=(4*int(valid.sum()), 6), squeeze=False)
        for col, index in enumerate(np.flatnonzero(valid)):
            top, bottom = axes[:, col]
            for k, (color, title) in enumerate(zip(("#255d9c", "#ce6828"), labels)):
                mean = stacks[k][:, :, index].mean(axis=0)
                sd = stacks[k][:, :, index].std(axis=0, ddof=1) if len(stacks[k]) > 1 else np.zeros_like(mean)
                top.plot(x, mean, color=color, label=title+" mean")
                top.fill_between(x, mean-sd, mean+sd, color=color, alpha=.15)
            top.plot(x, truths[0][field][:, index], "k--", linewidth=1, label="Truth")
            top.set(title=f"Profile time {times[index]:g} yr", ylabel=label)
            bottom.plot(x, stacks[1][:, :, index].mean(axis=0)-stacks[0][:, :, index].mean(axis=0), color="#ce6828", label="Mean difference")
            bottom.plot(x, truths[1][field][:, index]-truths[0][field][:, index], "k--", label="Truth difference")
            bottom.set(xlabel="Flowline distance (km)", ylabel=f"{labels[1]} − {labels[0]} (m)")
            top.grid(alpha=.2); bottom.grid(alpha=.2)
        axes[0, 0].legend(fontsize=8); axes[1, 0].legend(fontsize=8)
        fig.suptitle("Matched PIG runs: ensemble mean ± one ensemble standard deviation")
        fig.tight_layout(); fig.savefig(output / f"{field}_comparison.png", dpi=180); plt.close(fig)
    return reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path); parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference-label", default="Reference")
    parser.add_argument("--candidate-label", default="Candidate")
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    fields, metadata = compare_fields(args.reference, args.candidate)
    profiles = compare_profiles(args.reference, args.candidate, args.output, metadata.get("Nens"),
                                (args.reference_label, args.candidate_label))
    report = dict(reference=str(args.reference.resolve()), candidate=str(args.candidate.resolve()),
                  metadata=metadata, fields=fields, profiles=profiles,
                  interpretation="Measured differences; no scientific acceptance tolerance imposed.")
    (args.output / "comparison.json").write_text(json.dumps(report, indent=2)+"\n")
    with (args.output / "field_metrics.csv").open("w", newline="") as f:
        writer=csv.DictWriter(f, fieldnames=["field", "values", "max_abs", "rmse", "relative_l2"])
        writer.writeheader(); writer.writerows(dict(field=k, **v) for k, v in fields.items())
    print(json.dumps(fields, indent=2))


if __name__ == "__main__":
    main()
