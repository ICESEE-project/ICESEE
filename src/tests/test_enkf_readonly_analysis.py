# ==============================================================================
# @des: Regression tests for mode-2 analysis reads through read-only views.
#
# On PACE's NFS project storage, MPI-IO row-slab reads through read-write
# HDF5 handles stalled indefinitely (the first Idealized PIG mode-2 analysis
# made no progress for over an hour), while the same reads through read-only
# handles took 0.01 s.  Analysis reads therefore release the read-write
# window and use a read-only view; forecast writes keep MPI-IO locking.
# The NFS stall itself cannot be reproduced on local storage; these tests
# pin the handle lifecycle and numerical results around it.
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
from mpi4py import MPI

from ICESEE.src.parallelization.EnKF_parallel_io import EnKF_fully_parallel_IO
from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

_REPO_ROOT = Path(__file__).resolve().parents[2]
_WORKER = (
    Path(__file__).resolve().parent / "parallel_mpi" / "_enkf_io_readonly_analysis_worker.py"
)
_RESULT_RE = re.compile(r"RESULT rank=(\d+) errors=(\d+) total_errors=(\d+)")


def _open_state_intents(prefix):
    intents = []
    for fid in h5py.h5f.get_obj_ids(types=h5py.h5f.OBJ_FILE):
        name = fid.name.decode() if isinstance(fid.name, bytes) else str(fid.name)
        if prefix in name:
            intents.append((Path(name).name, fid.get_intent() == h5py.h5f.ACC_RDONLY))
    return intents


def _single_rank_io(tmp_path):
    comm = MPI.COMM_SELF
    return EnKF_fully_parallel_IO(
        "ro_ens", 4, 2, 5, comm, comm,
        {"ensemble_history_mode": "full", "force_fresh_start": True,
         "collective_threshold": 10_000},
        serial_file_creation=True, base_path=str(tmp_path), batch_size=2,
    )


def test_read_only_view_releases_window_and_write_access_returns_on_demand(tmp_path):
    io = _single_rank_io(tmp_path)
    io._ensure_batch_range(0, 1)
    io.write_forecast(1, np.arange(4.0), 0)
    assert all(not ro for _, ro in _open_state_intents("ro_ens"))

    with io._read_only_state(1) as dset:
        # No read-write handle or dataset object survives the transition.
        assert io.files == [] and io.datasets == [] and io.current_batch_start == -1
        assert _open_state_intents("ro_ens") == [("ro_ens_0001.h5", True)]
        np.testing.assert_array_equal(dset[:, 0], np.arange(4.0))
        with pytest.raises((OSError, ValueError)):
            dset[0, 0] = 1.0

    assert _open_state_intents("ro_ens") == []

    # The next forecast reopens the read-write window covering t and t + 1.
    np.testing.assert_array_equal(io.read_forecast(1, 0), np.arange(4.0))
    io.write_forecast(2, np.full(4, 7.0), 0)
    assert ("ro_ens_0002.h5", False) in _open_state_intents("ro_ens")
    np.testing.assert_array_equal(io.read_forecast(2, 0), np.full(4, 7.0))
    io.close()


def test_release_flushes_pending_writes_before_read_only_reopen(tmp_path):
    io = _single_rank_io(tmp_path)
    io._ensure_batch_range(2, 3)
    io.write_forecast(3, np.full(4, 3.5), 1)

    with io._read_only_state(3) as dset:
        np.testing.assert_array_equal(dset[:, 1], np.full(4, 3.5))
    io.close()


def test_read_only_view_closes_on_error(tmp_path):
    io = _single_rank_io(tmp_path)
    io._ensure_batch(0)

    with pytest.raises(RuntimeError, match="synthetic"):
        with io._read_only_state(0):
            raise RuntimeError("synthetic")

    assert _open_state_intents("ro_ens") == []
    io.close()


@pytest.mark.parametrize("world_size", [1, 3, 4, 8])
def test_real_mpi_analysis_reads_writes_windows_and_resume(tmp_path, world_size):
    if not h5py.get_config().mpi:
        pytest.skip("parallel HDF5 is unavailable in this interpreter")
    mpirun = find_compatible_mpi_launcher()
    if mpirun is None:
        pytest.skip("no compatible mpirun available in this environment")
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{_REPO_ROOT.parent}:{_REPO_ROOT}:{env.get('PYTHONPATH', '')}"
    out_dir = tmp_path / "enkf_io"
    out_dir.mkdir()

    result = subprocess.run(
        [mpirun, "--oversubscribe", "-n", str(world_size), sys.executable,
         str(_WORKER), str(out_dir)],
        capture_output=True, text=True, timeout=180, env=env,
    )

    output = result.stdout + result.stderr
    results = _RESULT_RE.findall(result.stdout)
    assert result.returncode == 0, output
    assert len(results) == world_size, output
    assert all(total == "0" for _, _, total in results), output


def test_serial_row_shard_assembly_preserves_values_and_canonical_until_publish(tmp_path):
    io = _single_rank_io(tmp_path)
    io._ensure_batch(1)
    for member in range(2):
        io.write_forecast(1, np.full(4, -1.0), member)
    io._release_batch()
    path, handle, values = io._create_rank_analysis_file(1)
    expected = np.arange(8.0).reshape(4, 2)
    values[:] = expected
    handle.close()
    temp = io._assemble_rank_analysis_files(1, path)
    with h5py.File(temp, 'r') as f:
        np.testing.assert_array_equal(f['states'][:], expected)
    with h5py.File(io._state_file_name(1), 'r') as f:
        np.testing.assert_array_equal(f['states'][:], np.full((4, 2), -1.0))
    io._publish_analysis_file(1, temp, {})
    io._write_published_analysis_mean(1)
    with h5py.File(io._state_file_name(1), 'r') as f:
        np.testing.assert_array_equal(f['states'][:], expected)
    with h5py.File(f'{io.base_path}/{io.file_prefix}_mean.h5', 'r') as f:
        np.testing.assert_array_equal(f['mean'][:, 1], expected.mean(axis=1))
    io.close()


def test_invalid_serial_row_shard_does_not_replace_canonical(tmp_path):
    io = _single_rank_io(tmp_path)
    io._ensure_batch(1)
    io._release_batch()
    canonical = Path(io._state_file_name(1))
    original = canonical.read_bytes()
    path, handle, _ = io._create_rank_analysis_file(1)
    handle.close()
    class GapComm:
        def gather(self, value, root):
            return [(1, 4, path)]
        def bcast(self, value, root):
            return value
    comm = io.mpi_comm
    io.mpi_comm = GapComm()
    try:
        with pytest.raises(RuntimeError, match='gap'):
            io._assemble_rank_analysis_files(1, path)
        assert canonical.read_bytes() == original
    finally:
        io.mpi_comm = comm
        io.close()
