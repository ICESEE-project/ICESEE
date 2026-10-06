# ==============================================================================
# @des: Reproducibility manifest for a PACE Mode-3 Idealized PIG run
# (2026-09-28, PACE calibration prep). Writes one JSON file capturing the
# environment/config context needed to interpret a run's timing/output,
# without dumping the full environment (no secrets, no raw env dump).
#
# Because the original ICESEE repo's .git object store is corrupted (see
# this session's git-corruption note) and this run's working tree carries
# uncommitted reconciliation work on top of a clean clone, "git revision"
# alone is not sufficient provenance. This script therefore also records:
#   - the clean clone's own git revision/branch (informational only -- it
#     predates the reconciliation and will not match the working tree);
#   - an explicit "working tree modified" flag plus a content-hash digest
#     of every tracked-relevant + untracked Idealized PIG/Mode-3 file, so
#     two runs can be compared for "same code" without needing a commit.
#
# Usage: python pace_run_manifest.py --out <path/manifest.json> [--icesee-kwargs-json <path>]
# ==============================================================================
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]

# The specific files this reconciliation touched or that materially affect
# a Mode-3 Idealized PIG run's behavior -- hashed individually (not the
# whole repo, which is large and includes unrelated applications) so two
# manifests can be diffed to see exactly which of these changed.
_TRACKED_RELEVANT_FILES = [
    "applications/icepack_model/examples/idealized_pig/_icepack_model.py",
    "applications/icepack_model/examples/idealized_pig/_icepack_enkf.py",
    "applications/icepack_model/examples/idealized_pig/_icepack_native.py",
    "applications/icepack_model/examples/idealized_pig/mode3_runner.py",
    "applications/icepack_model/examples/idealized_pig/run_da_icepack.py",
    "applications/icepack_model/examples/idealized_pig/params.yaml",
    "applications/icepack_model/icepack_utils/_distributed_fields.py",
    "config/_utility_imports.py",
    "src/parallelization/distributed_native_runtime.py",
    "src/parallelization/distributed_streaming_runtime.py",
    "src/parallelization/distributed_member_store.py",
    "src/parallelization/distributed_member_store_hdf5.py",
    "src/parallelization/distributed_checkpoint.py",
    "src/parallelization/mode3_projections.py",
]


def _run(cmd: list[str]) -> str | None:
    try:
        out = subprocess.run(
            cmd, cwd=str(_REPO_ROOT), capture_output=True, text=True, timeout=15,
        )
        if out.returncode != 0:
            return None
        return out.stdout.strip()
    except Exception:
        return None


def _git_provenance() -> dict:
    branch = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"])
    revision = _run(["git", "rev-parse", "HEAD"])
    status = _run(["git", "status", "--porcelain"])
    return {
        "branch": branch,
        "revision": revision,
        "revision_note": (
            "clean-clone HEAD only -- predates this session's uncommitted "
            "reconciliation work; use working_tree_digest below to identify "
            "the actual code that ran"
        ),
        "working_tree_modified": bool(status),
        "modified_file_count": len(status.splitlines()) if status else 0,
    }


def _file_digest(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()[:16]


def _working_tree_digest() -> dict:
    """SHA-256 of each relevant file's CURRENT on-disk content, regardless
    of git-tracked/untracked/modified status. Two manifests with identical
    per-file digests here ran byte-identical Idealized PIG/Mode-3 code,
    independent of commit state."""

    digests = {}
    for rel in _TRACKED_RELEVANT_FILES:
        digests[rel] = _file_digest(_REPO_ROOT / rel)
    combined = hashlib.sha256(
        json.dumps(digests, sort_keys=True).encode("utf-8")
    ).hexdigest()[:16]
    return {"per_file_sha256_16": digests, "combined_digest": combined}


def _python_env() -> dict:
    info = {"python_version": sys.version.replace("\n", " ")}
    for mod_name, attr_candidates in (
        ("mpi4py", ["__version__"]),
        ("h5py", ["__version__"]),
        ("firedrake", ["__version__"]),
        ("petsc4py", ["__version__"]),
        ("numpy", ["__version__"]),
    ):
        try:
            mod = __import__(mod_name)
            version = None
            for attr in attr_candidates:
                version = getattr(mod, attr, None)
                if version:
                    break
            info[f"{mod_name}_version"] = str(version) if version else "importable, no __version__"
        except ImportError:
            info[f"{mod_name}_version"] = None
    # HDF5 library version + MPI-parallel-I/O capability (h5py exposes both).
    try:
        import h5py
        info["hdf5_version"] = ".".join(str(v) for v in h5py.h5.get_libversion())
        info["h5py_has_mpi"] = bool(getattr(h5py.get_config(), "mpi", False))
    except Exception:
        info["hdf5_version"] = None
        info["h5py_has_mpi"] = None
    try:
        from mpi4py import MPI
        info["mpi_vendor"] = MPI.get_vendor()
        info["mpi_world_size"] = MPI.COMM_WORLD.Get_size()
    except Exception:
        info["mpi_vendor"] = None
        info["mpi_world_size"] = None
    return info


def build_manifest(icesee_kwargs: dict | None = None) -> dict:
    manifest = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "git": _git_provenance(),
        "working_tree_digest": _working_tree_digest(),
        "python_env": _python_env(),
        "slurm": {
            key: os.environ.get(key)
            for key in (
                "SLURM_JOB_ID", "SLURM_JOB_NAME", "SLURM_JOB_NUM_NODES",
                "SLURM_NTASKS", "SLURM_NTASKS_PER_NODE", "SLURM_CPUS_PER_TASK",
                "SLURM_SUBMIT_DIR", "SLURM_JOB_PARTITION", "SLURM_JOB_ACCOUNT",
            )
            if os.environ.get(key) is not None
        },
    }
    if icesee_kwargs:
        run_config_keys = (
            "Nens", "model_nprocs", "num_years", "timesteps_per_year",
            "num_state_vars", "num_param_vars", "joint_estimation",
            "vec_inputs", "random_field_method", "execution_mode",
            "use_member_streaming", "use_store_streaming_analysis",
            "member_store_backend", "member_store_root", "data_path",
            "compact_initialization", "initFile", "freq_obs",
            "obs_start_time", "obs_max_time", "run_id",
        )
        manifest["run_config"] = {
            k: icesee_kwargs[k] for k in run_config_keys if k in icesee_kwargs
        }
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="Path to write manifest.json")
    parser.add_argument(
        "--icesee-kwargs-json", default=None,
        help="Optional path to a JSON dump of the run's icesee_kwargs "
             "(only the run_config_keys subset is captured -- no secrets, "
             "no full environment)",
    )
    args = parser.parse_args()

    icesee_kwargs = None
    if args.icesee_kwargs_json and os.path.exists(args.icesee_kwargs_json):
        with open(args.icesee_kwargs_json) as f:
            icesee_kwargs = json.load(f)

    manifest = build_manifest(icesee_kwargs)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(manifest, f, indent=2, default=str)
    print(f"[pace_run_manifest] wrote {out_path}")


if __name__ == "__main__":
    main()
