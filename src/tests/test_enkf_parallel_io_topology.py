# ==============================================================================
# @des: Real-MPI targeted tests for Stage 4B's HDF5/collective fix
# (the unconditional per-timestep priming collective added to
# _mpi_forecast_functions.py), using EnKF_fully_parallel_IO directly with
# synthetic deterministic data via _enkf_io_topology_worker.py.
#
# The real Lorenz96 P/Nens matrix (test_lorenz96_mode2_p_nens_matrix.py)
# already covers every topology combination at ranks_per_model=1,
# including spare ranks (P=10,Nens=4) and a partial-round case (P=3,
# Nens=4) end-to-end through the real DA cycle. No current application is
# registered supports_multi_rank_per_model=True (model_capabilities.py),
# so ranks_per_model > 1 cannot be exercised through a real application --
# this file covers exactly that gap using a synthetic worker that talks to
# EnKF_fully_parallel_IO directly, bypassing the model layer, per Stage
# 4B's own guidance ("use a synthetic topology/I/O worker rather than
# bypassing capability validation").
# ==============================================================================
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

_REPO_ROOT = Path(__file__).resolve().parents[2]
_WORKER = Path(__file__).resolve().parent / "parallel_mpi" / "_enkf_io_topology_worker.py"

_RESULT_RE = re.compile(
    r"RESULT rank=(?P<rank>\d+) is_spare=(?P<is_spare>True|False) "
    r"n_errors_seen_by_this_rank=(?P<n_errors>\d+) total_errors=(?P<total_errors>\d+)"
)


def _run_worker(tmp_path, world_size, nens, ranks_per_model, n_timesteps=5, batch_size=2):
    mpirun = find_compatible_mpi_launcher()
    if mpirun is None:
        pytest.skip("no compatible mpirun available in this environment")

    env = dict(os.environ)
    env["PYTHONPATH"] = f"{_REPO_ROOT.parent}:{_REPO_ROOT}"
    out_dir = tmp_path / "enkf_io"
    out_dir.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        [
            mpirun, "--oversubscribe", "-n", str(world_size), sys.executable, str(_WORKER),
            str(nens), str(ranks_per_model), str(out_dir), str(n_timesteps), str(batch_size),
        ],
        capture_output=True,
        text=True,
        timeout=90,
        env=env,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    results = {}
    for line in result.stdout.splitlines():
        m = _RESULT_RE.match(line)
        if m:
            results[int(m.group("rank"))] = m.groupdict()
    done_ranks = {
        int(line.split("rank=")[1])
        for line in result.stdout.splitlines()
        if line.startswith("DONE")
    }
    return results, done_ranks, out_dir, result.stdout


@pytest.mark.parametrize(
    "world_size,nens,ranks_per_model",
    [
        (8, 4, 2),    # multi-rank-per-model, no spares
        (10, 4, 2),   # multi-rank-per-model + spare ranks (the P=10/Nens=4 topology)
        (4, 5, 2),    # multi-rank-per-model + a partial final round
    ],
)
def test_enkf_io_completes_with_zero_data_errors(tmp_path, world_size, nens, ranks_per_model):
    results, done_ranks, out_dir, stdout = _run_worker(tmp_path, world_size, nens, ranks_per_model)

    assert set(results) == set(range(world_size)), stdout
    assert done_ranks == set(range(world_size)), stdout

    for rank, info in results.items():
        assert int(info["total_errors"]) == 0, (
            f"rank {rank} observed {info['total_errors']} data errors "
            f"(cross-member/timestep contamination or lost writes):\n{stdout}"
        )


def test_enkf_io_repeated_reads_writes_across_batch_boundary(tmp_path):
    # batch_size=2 with 5 timesteps forces at least two batch-window
    # transitions; a hidden reopen mismatch would surface as either a hang
    # (caught by the subprocess timeout) or a data error (caught above).
    results, done_ranks, out_dir, stdout = _run_worker(
        tmp_path, world_size=8, nens=4, ranks_per_model=2, n_timesteps=6, batch_size=2
    )
    assert done_ranks == set(range(8)), stdout
    for rank, info in results.items():
        assert int(info["total_errors"]) == 0, stdout


def test_enkf_io_spare_ranks_present_and_finalize_cleanly(tmp_path):
    results, done_ranks, out_dir, stdout = _run_worker(tmp_path, world_size=10, nens=4, ranks_per_model=2)
    spare = {r for r, info in results.items() if info["is_spare"] == "True"}
    assert spare == {8, 9}, stdout
    assert done_ranks == set(range(10)), stdout


def test_enkf_io_file_contents_have_expected_shape_and_no_missing_members(tmp_path):
    results, done_ranks, out_dir, stdout = _run_worker(
        tmp_path, world_size=8, nens=4, ranks_per_model=2, n_timesteps=3, batch_size=2
    )
    assert done_ranks == set(range(8)), stdout

    shard_files = sorted(out_dir.glob("topology_probe_ens_*.h5"))
    assert len(shard_files) >= 1, f"no shard files found in {out_dir}: {stdout}"

    nd = 6
    nens = 4
    for shard_path in shard_files:
        with h5py.File(shard_path, "r") as f:
            assert "states" in f, f"{shard_path} missing 'states' dataset"
            data = f["states"][:]
            assert data.shape == (nd, nens), (
                f"{shard_path}: expected shape {(nd, nens)}, got {data.shape}"
            )
            assert np.all(np.isfinite(data))
            # Every member column must be internally consistent (all nd
            # entries share the same synthetic fill value: 1000*ens_id+t)
            # and every one of the nens members must be present (no
            # all-zero "never written" column, which would indicate a
            # missing member).
            for ens_id in range(nens):
                column = data[:, ens_id]
                assert np.allclose(column, column[0]), (
                    f"{shard_path} member {ens_id}: cross-member contamination "
                    f"within one column: {column}"
                )
                assert not np.allclose(column, 0.0), (
                    f"{shard_path} member {ens_id}: looks unwritten (all zero)"
                )
