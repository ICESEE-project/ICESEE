# ==============================================================================
# @des: Level-3 (lightweight real-model integration) test: runs the actual
#       ICESEE mode-2 driver on Lorenz-96 across the P/Nens configurations
#       from Stage 3M, including P > Nens (multi-rank-per-member), and
#       validates genuine completion -- not merely "mpirun returned".
#
#       A prior run was observed to report success (mpirun exit code 0)
#       while every rank had actually crashed with an unhandled KeyError:
#       the exception was swallowed by icesee_da_full_parallel.py's error
#       handler with no re-raise and no abort. This test's assertions are
#       deliberately layered (exit code, absence of error markers in
#       output, exact expected shape, finite values, presence of the
#       ownership diagnostic) so that a similar false-success cannot slip
#       through on any single check.
# ==============================================================================
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pytest
import yaml

from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

REPO_ROOT = Path(__file__).resolve().parents[2]
LORENZ_DIR = REPO_ROOT / "applications" / "lorenz_model" / "examples" / "lorenz96"

_LAUNCHER = find_compatible_mpi_launcher()

pytestmark = pytest.mark.skipif(
    _LAUNCHER is None, reason="no MPI launcher compatible with this mpi4py build"
)

# (P, Nens, label) -- the exact matrix from Stage 3M, including every case
# that previously hung or false-passed.
MATRIX = [
    (1, 1, "P=1,Nens=1"),
    (2, 4, "P=2,Nens=4"),
    (3, 4, "P=3,Nens=4 (Nens>=P, non-divisible rounds)"),
    (4, 4, "P=4,Nens=4"),
    (8, 4, "P=8,Nens=4 (Nens<P)"),
    (10, 4, "P=10,Nens=4 (Nens<P, non-divisible)"),
    (16, 4, "P=16,Nens=4 (oversubscribed)"),
]

EXPECTED_NT = 1000  # num_years / dt from lorenz96/params.yaml (10 / 0.01)

_ERROR_MARKERS = ("Fatal error", "Traceback (most recent call last)", "KeyError", "IndexError")


def _make_config(nens, data_path, tmp_path):
    with (LORENZ_DIR / "params.yaml").open() as f:
        config = yaml.safe_load(f)
    enkf = config["enkf-parameters"]
    enkf["execution_mode"] = 2
    enkf["Nens"] = nens
    enkf["data_path"] = str(data_path)
    enkf["restart_enabled"] = False
    enkf["force_fresh_start"] = True
    dest = tmp_path / f"lorenz_nens{nens}.yaml"
    with dest.open("w") as f:
        yaml.safe_dump(config, f, sort_keys=False)
    return dest


def _run_case(p, nens, tmp_path, timeout=120):
    data_path = tmp_path / f"data_p{p}_nens{nens}"
    config = _make_config(nens, data_path, tmp_path)
    # --oversubscribe is an Open MPI-ism; only pass it when supported, so
    # this test does not assume a specific MPI implementation beyond what
    # find_compatible_mpi_launcher already verified works.
    cmd = [_LAUNCHER]
    if _supports_oversubscribe():
        cmd.append("--oversubscribe")
    cmd += [
        "-n", str(p), sys.executable, "-m",
        "ICESEE.applications.lorenz_model.examples.lorenz96.run_da_lorenz96",
        "-F", str(config),
    ]
    result = subprocess.run(
        cmd, cwd=LORENZ_DIR, capture_output=True, text=True, timeout=timeout,
    )
    return result, data_path


def _launcher_flavor():
    try:
        out = subprocess.run([_LAUNCHER, "--version"], capture_output=True, text=True, timeout=10)
        return out.stdout
    except Exception:
        return ""


def _supports_oversubscribe():
    return "Open MPI" in _launcher_flavor() or "open-mpi" in _launcher_flavor().lower()


@pytest.mark.parametrize("p,nens,label", MATRIX, ids=[m[2] for m in MATRIX])
def test_lorenz96_mode2_completes_genuinely(p, nens, label, tmp_path):
    try:
        result, data_path = _run_case(p, nens, tmp_path)
    except subprocess.TimeoutExpired as exc:
        pytest.fail(
            f"{label}: hung past the timeout instead of completing or "
            f"failing cleanly. Partial stdout tail: "
            f"{(exc.stdout or '')[-2000:]!r}"
        )

    combined = result.stdout + result.stderr

    # 1. Process exit code.
    assert result.returncode == 0, (
        f"{label}: non-zero exit ({result.returncode}).\n"
        f"stdout tail:\n{result.stdout[-3000:]}\nstderr tail:\n{result.stderr[-2000:]}"
    )

    # 2. No error/traceback markers, even if the exit code were somehow 0
    #    (the exact false-success shape previously observed).
    for marker in _ERROR_MARKERS:
        assert marker not in combined, (
            f"{label}: found error marker {marker!r} in output despite "
            f"exit code 0 -- this is the false-success pattern this test "
            f"exists to catch.\noutput tail:\n{combined[-3000:]}"
        )

    # 3. Ownership diagnostic present and reports the replicated invariant
    #    (Lorenz-96 is configured state_distribution: replicated -- see
    #    applications/lorenz_model/examples/lorenz96/params.yaml).
    ownership_lines = [ln for ln in result.stdout.splitlines() if "[ICESEE][ownership]" in ln]
    if p > nens:
        # Only the Nens < P branch prints this diagnostic today.
        assert ownership_lines, f"{label}: expected an ownership diagnostic line, found none"
        line = ownership_lines[0]
        assert "distribution=replicated" in line, f"{label}: {line}"
        # The defining replicated invariant: global_size == local_size
        # regardless of the model communicator's size.
        import re
        local = int(re.search(r"local_size=(\d+)", line).group(1))
        glob = int(re.search(r"global_size=(\d+)", line).group(1))
        assert glob == local, (
            f"{label}: replicated ownership must have global_size == "
            f"local_size regardless of model_comm_size; got {line}"
        )

    # 4. Output file exists with the expected shape and finite values.
    mean_file = data_path / "icesee_ensemble_data.h5"
    assert mean_file.exists(), f"{label}: missing output file {mean_file}"
    with h5py.File(mean_file, "r") as f:
        mean = f["ensemble_mean"][:]
    assert mean.shape == (3, EXPECTED_NT + 1), (
        f"{label}: expected ensemble_mean shape (3, {EXPECTED_NT + 1}) "
        f"(a full {EXPECTED_NT}-step run), got {mean.shape} -- a smaller "
        f"shape indicates the run silently truncated rather than "
        f"completing every timestep."
    )
    assert np.all(np.isfinite(mean)), f"{label}: non-finite values in ensemble_mean"


def test_explicit_ranks_per_model_1_with_spare_ranks_completes(tmp_path):
    """P=6, Nens=4, ranks_per_model=1 explicit -> 6 singleton model groups,
    only 4 used, 2 spare ranks.

    This specific combination (an *explicit* ranks_per_model=1 request
    together with world_size > Nens) was unreachable before Stage 4A's
    ResourcePlan generalization -- the legacy auto policy (ranks_per_model
    omitted) never produces a spare rank when it resolves to 1. Once it
    became reachable, it exposed two real, generic (non-Icepack-specific)
    bugs in the shared Mode-2 runtime: `_mpi_ensemble_intialization.py`'s
    ranks_per_model == 1 branch called icesee_get_index (which touches
    `comm.Get_size()`) before its `if color is not None:` guard, so a
    spare rank's MPI.COMM_NULL subcomm raised MPI_ERR_COMM; and
    `_mpi_forecast_functions.py`'s forecast step did the same for its own
    outer icesee_get_index call. Both are fixed generically (not with an
    Icepack-specific workaround) in this same continuation.
    """
    with (LORENZ_DIR / "params.yaml").open() as f:
        config = yaml.safe_load(f)
    enkf = config["enkf-parameters"]
    enkf["execution_mode"] = 2
    enkf["Nens"] = 4
    enkf["ranks_per_model"] = 1
    data_path = tmp_path / "data_spare_ranks"
    enkf["data_path"] = str(data_path)
    enkf["restart_enabled"] = False
    enkf["force_fresh_start"] = True
    config_path = tmp_path / "lorenz_spare_ranks.yaml"
    with config_path.open("w") as f:
        yaml.safe_dump(config, f, sort_keys=False)

    cmd = [_LAUNCHER]
    if _supports_oversubscribe():
        cmd.append("--oversubscribe")
    cmd += [
        "-n", "6", sys.executable, "-m",
        "ICESEE.applications.lorenz_model.examples.lorenz96.run_da_lorenz96",
        "-F", str(config_path),
    ]
    result = subprocess.run(
        cmd, cwd=LORENZ_DIR, capture_output=True, text=True, timeout=120,
    )
    combined = result.stdout + result.stderr
    assert result.returncode == 0, (
        f"non-zero exit ({result.returncode}).\n"
        f"stdout tail:\n{result.stdout[-3000:]}\nstderr tail:\n{result.stderr[-2000:]}"
    )
    for marker in _ERROR_MARKERS:
        assert marker not in combined, (
            f"found error marker {marker!r} despite exit code 0.\n"
            f"output tail:\n{combined[-3000:]}"
        )
    assert "spare=2" in result.stdout, (
        f"expected the topology line to report spare=2 (6 world ranks, "
        f"4 singleton model groups used, 2 spare); stdout tail:\n"
        f"{result.stdout[-2000:]}"
    )

    mean_file = data_path / "icesee_ensemble_data.h5"
    assert mean_file.exists(), f"missing output file {mean_file}"
    with h5py.File(mean_file, "r") as f:
        mean = f["ensemble_mean"][:]
    assert mean.shape == (3, EXPECTED_NT + 1)
    assert np.all(np.isfinite(mean))
