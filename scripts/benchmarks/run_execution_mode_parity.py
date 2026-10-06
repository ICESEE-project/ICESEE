#!/usr/bin/env python3
"""Run one ICESEE case in modes 1 and 2, then compare every ensemble value.

The two runs are generated from one YAML file.  Only ``execution_mode`` and
``data_path`` differ.  Full history and float64 storage are forced because a
rolling history or float32 storage cannot satisfy the strict parity test.
The default gate permits only round-off-level differences caused by distinct
MPI reduction orders.  ``--bitwise`` is retained as a platform diagnostic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import uuid
from pathlib import Path

import yaml

from compare_execution_modes import compare


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_identity(repository: Path) -> dict[str, str | None]:
    """Record the exact source tree used by a parity run.

    The worktree hash is important here: parity development is commonly run
    before committing a fix, so a commit SHA alone cannot distinguish old and
    current implementations.
    """
    identity: dict[str, str | None] = {"git_head": None, "git_diff_sha256": None}
    try:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repository, check=True,
            capture_output=True,
        ).stdout
        diff = subprocess.run(
            ["git", "diff", "--binary", "HEAD"], cwd=repository, check=True,
            capture_output=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return identity
    identity["git_head"] = head.decode("ascii").strip()
    identity["git_diff_sha256"] = hashlib.sha256(diff).hexdigest()
    return identity


def _write_config(
    source: Path,
    destination: Path,
    mode: int,
    data_path: Path,
    run_cwd: Path,
    num_years: float | None = None,
):
    with source.open("r", encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    enkf = config.setdefault("enkf-parameters", {})
    enkf.update(
        execution_mode=mode,
        # All paired applications are launched from ``run_cwd``.  A relative
        # value works for both Python adapters and MATLAB-coupled ISSM, whose
        # helper scripts intentionally resolve data_path against the
        # application directory rather than accepting an arbitrary absolute
        # host path.
        data_path=os.path.relpath(data_path.resolve(), run_cwd.resolve()),
        ensemble_history_mode="full",
        ensemble_storage_dtype="float64",
        ensemble_finalize_mode="vds",
        restart_enabled=False,
        force_fresh_start=True,
        # Emit per-member input/model-output/post-noise traces.  This is
        # intentionally benchmark-only and is invaluable when the first
        # divergent stage must be identified without changing production I/O.
        execution_parity_trace=True,
    )
    if num_years is not None:
        modeling = config.setdefault("modeling-parameters", {})
        modeling["num_years"] = float(num_years)
        # Observation generation must not request columns beyond the shortened
        # model trajectory.  Preserve the source schedule otherwise.
        if "obs_max_time" in enkf:
            enkf["obs_max_time"] = min(
                float(enkf["obs_max_time"]), float(num_years)
            )
        # A shortened parity case must not map future bed surveys onto its
        # final observation column.  Keep only snapshots that physically
        # occur inside the shortened trajectory.  The same filtered schedule
        # is written to both mode configurations.
        snapshots = enkf.get("bed_obs_snapshot")
        if isinstance(snapshots, (list, tuple)):
            enkf["bed_obs_snapshot"] = [
                float(value)
                for value in snapshots
                if 0.0 <= float(value) <= float(num_years) + 1.0e-12
            ]
    with destination.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(config, stream, sort_keys=False)


def _run(command: list[str], cwd: Path):
    print("+", " ".join(command), flush=True)
    environment = os.environ.copy()
    environment["PYTHONUNBUFFERED"] = "1"
    subprocess.run(command, cwd=cwd, check=True, env=environment)


def _state_outputs(path: Path, mode: int) -> list[Path]:
    if mode == 1:
        expected = path / "icesee_ensemble_data.h5"
        return [expected] if expected.is_file() else []
    shards = sorted(path.glob("icesee_enkf_ens_*.h5"))
    if shards:
        return shards
    finalized = path / "icesee_ensemble_data.h5"
    return [finalized] if finalized.is_file() else []


def _require_mode_output(path: Path, mode: int, launched_at: float):
    """Catch ICESEE runners that report an internal failure but exit zero."""
    outputs = _state_outputs(path, mode)
    if not outputs:
        raise RuntimeError(
            f"Mode {mode} did not produce a recognizable ensemble in {path}"
        )
    stale = [item for item in outputs if item.stat().st_mtime < launched_at]
    if stale:
        names = ", ".join(str(item) for item in stale[:3])
        raise RuntimeError(
            f"Mode {mode} output predates this launch ({names}). Refusing to "
            "compare stale state files."
        )
    return outputs


def _write_provenance(path: Path, *, run_id: str, mode: int,
                      launched_at: float, outputs: list[Path],
                      source_config: Path, generated_config: Path,
                      application: Path, repository: Path):
    """Bind each mode directory to this exact benchmark invocation."""
    payload = {
        "schema": 2,
        "run_id": run_id,
        "execution_mode": mode,
        "launched_at_unix": launched_at,
        "completed_at_unix": time.time(),
        "outputs": [item.name for item in outputs],
        "output_mtime_ns": {
            item.name: item.stat().st_mtime_ns for item in outputs
        },
        "source_config": str(source_config),
        "source_config_sha256": _sha256_file(source_config),
        "generated_config": str(generated_config),
        "generated_config_sha256": _sha256_file(generated_config),
        "application": str(application),
        "application_sha256": _sha256_file(application),
        **_git_identity(repository),
    }
    with (path / "execution_parity_provenance.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)


def _write_launch_manifest(path: Path, *, run_id: str, source_config: Path,
                           application: Path, repository: Path):
    """Fingerprint a parity invocation before either expensive run starts."""
    payload = {
        "schema": 1,
        "run_id": run_id,
        "status": "running",
        "launched_at_unix": time.time(),
        "source_config": str(source_config),
        "source_config_sha256": _sha256_file(source_config),
        "application": str(application),
        "application_sha256": _sha256_file(application),
        **_git_identity(repository),
    }
    with (path / "execution_parity_launch.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)


def _validate_mpi_launcher(mpiexec: str, python: str, ranks: int, cwd: Path):
    """Reject launchers that create independent singleton MPI processes."""
    probe = (
        "from mpi4py import MPI; import sys; "
        f"sys.exit(0 if MPI.COMM_WORLD.Get_size() == {ranks} else 86)"
    )
    try:
        _run([mpiexec, "-np", str(ranks), python, "-c", probe], cwd)
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"MPI preflight failed for {mpiexec!r}. The launcher must belong "
            "to the same MPI implementation as mpi4py and must create one "
            f"communicator of size {ranks}. Pass the matching executable with "
            "--mpiexec (for example /opt/homebrew/bin/mpirun)."
        ) from error


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--application", required=True, type=Path,
                        help="Application entry point, for example run_da_lorenz96.py")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--ranks", type=int, default=2)
    parser.add_argument("--nens", type=int)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--row-chunk", type=int, default=8192)
    parser.add_argument(
        "--application-arg",
        action="append",
        default=[],
        help=(
            "Additional argument forwarded to both application launches. "
            "Repeat as needed; for an argument beginning with '--', use "
            "--application-arg=--name=value."
        ),
    )
    parser.add_argument("--atol", type=float, default=1.0e-10)
    parser.add_argument("--rtol", type=float, default=1.0e-8)
    parser.add_argument(
        "--bitwise", action="store_true",
        help="Require exact equality (diagnostic only; not MPI-portable).",
    )
    parser.add_argument(
        "--num-years",
        type=float,
        help=(
            "Optional shortened model horizon for fast parity diagnosis. "
            "For the Lorenz96 dt=0.01 case, 0.2 exercises 20 forecasts and "
            "at least one analysis."
        ),
    )
    args = parser.parse_args()

    if args.ranks < 1 or args.row_chunk < 1:
        parser.error("--ranks and --row-chunk must be positive")

    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error(
            f"--output must be new or empty; refusing to overwrite {output}"
        )
    output.mkdir(parents=True, exist_ok=True)
    run_id = str(uuid.uuid4())
    repository = Path(__file__).resolve().parents[2]
    source_config = args.config.resolve()
    partial = output / "mode1_partial"
    full = output / "mode2_full"
    config1 = output / "mode1.yaml"
    config2 = output / "mode2.yaml"
    if args.num_years is not None and args.num_years <= 0:
        parser.error("--num-years must be positive")
    application = args.application.resolve()
    application_cwd = application.parent
    _write_launch_manifest(
        output,
        run_id=run_id,
        source_config=source_config,
        application=application,
        repository=repository,
    )
    _write_config(
        source_config, config1, 1, partial, application_cwd,
        args.num_years,
    )
    _write_config(
        source_config, config2, 2, full, application_cwd,
        args.num_years,
    )

    print(f"[parity] run id: {run_id}")
    print(f"[parity] fresh output root: {output}")
    _validate_mpi_launcher(args.mpiexec, args.python, args.ranks, application_cwd)
    common = [
        args.mpiexec,
        "-np",
        str(args.ranks),
        args.python,
        str(application),
        *args.application_arg,
    ]
    nens = [] if args.nens is None else [f"--Nens={args.nens}"]
    mode1_started = time.time()
    _run(common + nens + ["-F", str(config1)], application_cwd)
    mode1_outputs = _require_mode_output(partial, 1, mode1_started)
    _write_provenance(
        partial, run_id=run_id, mode=1, launched_at=mode1_started,
        outputs=mode1_outputs, source_config=source_config,
        generated_config=config1, application=application,
        repository=repository,
    )
    mode2_started = time.time()
    _run(common + nens + ["-F", str(config2)], application_cwd)
    mode2_outputs = _require_mode_output(full, 2, mode2_started)
    _write_provenance(
        full, run_id=run_id, mode=2, launched_at=mode2_started,
        outputs=mode2_outputs, source_config=source_config,
        generated_config=config2, application=application,
        repository=repository,
    )

    rtol = 0.0 if args.bitwise else args.rtol
    atol = 0.0 if args.bitwise else args.atol
    passed = compare(
        partial,
        full,
        row_chunk=args.row_chunk,
        rtol=rtol,
        atol=atol,
        require_provenance=True,
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
