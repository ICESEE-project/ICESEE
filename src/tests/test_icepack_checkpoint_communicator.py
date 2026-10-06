# ==============================================================================
# @des: Regression tests for the Icepack checkpoint-communicator fix in
# `initializeRun` (applications/icepack_model/examples/idealized_pig/_icepack_model.py).
#
# Background: historical test tree passed comm=icesee_kwargs["comm"] to
# firedrake.CheckpointFile; main previously omitted it, silently falling
# back to CheckpointFile's own default (COMM_WORLD, baked into
# firedrake.checkpointing.CheckpointFile.__init__). A real, tiny 2-rank
# Firedrake probe (two independent single-rank "spatial groups" each
# writing and reading their own checkpoint file) confirmed this is a real
# correctness bug, not just theoretical: with no comm= (both ranks
# collectively opening CheckpointFile on COMM_WORLD despite each naming a
# different file), the read did not hang but returned silently WRONG
# values -- worse than a deadlock. With comm=<this rank's own group comm>,
# both ranks read back exactly the values they wrote, fast, no
# cross-contamination.
#
# These tests cover two things:
#   1. (real Firedrake, mocked CheckpointFile) initializeRun must still
#      forward icesee_kwargs["comm"] to CheckpointFile -- proven by
#      substituting a fake CheckpointFile that records the comm it was
#      given and aborts immediately (before touching any real physics),
#      so this test stays cheap regardless of initializeRun's own solve/
#      SMB-file work further down.
#   2. (real Firedrake, real CheckpointFile, real MPI) the actual
#      independent-spatial-group scenario mode 3 depends on: explicit
#      comm= gives each group its own correct data; the previous
#      default-comm behavior is demonstrated to corrupt data across
#      groups, run once here (single real subprocess) as the concrete,
#      reproducible proof behind the fix above -- not merely asserted.
# ==============================================================================
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

_REPO_ROOT = Path(__file__).resolve().parents[2]
_IDEALIZED_PIG_DIR = _REPO_ROOT / "applications" / "icepack_model" / "examples" / "idealized_pig"
for _extra_path in (str(_REPO_ROOT), str(_REPO_ROOT.parent), str(_IDEALIZED_PIG_DIR)):
    if _extra_path not in sys.path:
        sys.path.insert(0, _extra_path)

_SAFE_DATA_PATH = tempfile.mkdtemp(prefix="icesee_icepack_checkpoint_comm_test_")
_ARGV_BACKUP = sys.argv[:]
_CWD_BACKUP = os.getcwd()
sys.argv = [sys.argv[0], "--data_path", _SAFE_DATA_PATH]
os.chdir(_IDEALIZED_PIG_DIR)
try:
    from applications.icepack_model.examples.idealized_pig import _icepack_model as model
finally:
    sys.argv = _ARGV_BACKUP
    os.chdir(_CWD_BACKUP)

import firedrake


class _CommCaptured(Exception):
    def __init__(self, comm):
        self.comm = comm


def test_initialize_run_forwards_icesee_kwargs_comm_to_checkpoint_file(monkeypatch):
    """initializeRun must pass icesee_kwargs["comm"] to CheckpointFile, not
    silently rely on CheckpointFile's own COMM_WORLD default. Aborts inside
    the fake CheckpointFile's __init__ -- before initializeRun's real
    diagnostic_solve/readSMB calls -- so this test never touches real
    physics or an SMB file."""

    sentinel_comm = object()

    class _FakeCheckpointFile:
        def __init__(self, filename, mode, comm=None):
            raise _CommCaptured(comm)

    monkeypatch.setattr(model.firedrake, "CheckpointFile", _FakeCheckpointFile)

    icesee_kwargs = {"initFile": "unused.h5", "comm": sentinel_comm}

    with pytest.raises(_CommCaptured) as excinfo:
        model.initializeRun(icesee_kwargs, forward_solver=None, mesh=None, Q=None, V=None)

    assert excinfo.value.comm is sentinel_comm


def test_initialize_run_falls_back_to_comm_world_when_comm_key_absent(monkeypatch):
    """Matches Firedrake's own `comm or COMM_WORLD` convention: an
    icesee_kwargs without a 'comm' key (never expected in practice -- every
    real caller sets it -- but a safe, explicit fallback) must still reach
    CheckpointFile with a real communicator, not None."""

    class _FakeCheckpointFile:
        def __init__(self, filename, mode, comm=None):
            raise _CommCaptured(comm)

    monkeypatch.setattr(model.firedrake, "CheckpointFile", _FakeCheckpointFile)

    icesee_kwargs = {"initFile": "unused.h5"}  # no "comm" key at all

    with pytest.raises(_CommCaptured) as excinfo:
        model.initializeRun(icesee_kwargs, forward_solver=None, mesh=None, Q=None, V=None)

    assert excinfo.value.comm is firedrake.COMM_WORLD


@pytest.mark.parametrize("scoped", [True, False])
def test_independent_spatial_groups_checkpoint_io_real_mpi(tmp_path, scoped):
    """Real, tiny, 2-rank Firedrake/MPI reproduction of the exact scenario
    Icepack mode 3 depends on: two independent single-rank "spatial
    groups" each writing and reading their own checkpoint file.

    scoped=True (comm=<this rank's own group comm>, what the fix does):
    each rank must read back exactly the values it wrote.

    scoped=False (no comm=, CheckpointFile's own COMM_WORLD default, the
    pre-fix behavior): demonstrates the actual failure mode found during
    investigation -- not a hang, but silently WRONG values, read back
    faster than a real hang would even allow, i.e. this assertion is
    expected to catch corrupted data, proving the bug this fix addresses
    was real and reproducible, not theoretical.
    """
    script = tmp_path / "probe.py"
    script.write_text(_PROBE_SCRIPT)
    mode = "scoped" if scoped else "default"

    mpirun = find_compatible_mpi_launcher()
    if mpirun is None:
        pytest.skip("no compatible mpirun available in this environment")

    result = subprocess.run(
        [mpirun, "-n", "2", sys.executable, str(script), mode, str(tmp_path)],
        capture_output=True,
        text=True,
        timeout=60,
    )

    if scoped:
        # Correct behavior: must succeed, and each rank reads back exactly
        # what it wrote.
        assert result.returncode == 0, result.stdout + result.stderr
        rank0_line = next(l for l in result.stdout.splitlines() if l.startswith("[rank 0]"))
        rank1_line = next(l for l in result.stdout.splitlines() if l.startswith("[rank 1]"))
        rank0_sum = float(rank0_line.rsplit("sum=", 1)[1])
        rank1_sum = float(rank1_line.rsplit("sum=", 1)[1])
        assert rank0_sum == pytest.approx(10.0)
        assert rank1_sum == pytest.approx(510.0)
    else:
        # Pre-fix behavior (mismatched-collective CheckpointFile.__init__
        # calls, one per independent group, all defaulting to COMM_WORLD):
        # undefined behavior at the MPI/HDF5 layer. Investigation observed
        # this manifest as silently WRONG values rather than a hang or a
        # clean error, but accept either wrong-values or a nonzero exit as
        # confirmation this path is broken -- the one outcome that must
        # NOT happen is a clean success with both ranks' correct values,
        # which would mean the bug this fix addresses wasn't real.
        if result.returncode == 0:
            rank0_line = next(l for l in result.stdout.splitlines() if l.startswith("[rank 0]"))
            rank1_line = next(l for l in result.stdout.splitlines() if l.startswith("[rank 1]"))
            rank0_sum = float(rank0_line.rsplit("sum=", 1)[1])
            rank1_sum = float(rank1_line.rsplit("sum=", 1)[1])
            assert not (rank0_sum == pytest.approx(10.0) and rank1_sum == pytest.approx(510.0)), (
                "expected the pre-fix default-comm path to corrupt data or fail, "
                "but it returned exactly correct values"
            )


_PROBE_SCRIPT = '''
import sys
import numpy as np
from mpi4py import MPI
import firedrake
from firedrake import UnitIntervalMesh, FunctionSpace, Function

mode = sys.argv[1]
out_dir = sys.argv[2]
world = MPI.COMM_WORLD
rank = world.Get_rank()
group_comm = world.Split(color=rank, key=0)

filename = f"{out_dir}/checkpoint_probe_rank{rank}.h5"
mesh_name = f"probe_mesh_{rank}"

mesh = UnitIntervalMesh(4, comm=group_comm, name=mesh_name)
V = FunctionSpace(mesh, "CG", 1)
f = Function(V, name="probe")
f.dat.data[:] = np.arange(f.dat.data.size, dtype=float) + 100.0 * rank

with firedrake.CheckpointFile(filename, "w", comm=group_comm) as chk:
    chk.save_mesh(mesh)
    chk.save_function(f)
world.Barrier()

if mode == "scoped":
    with firedrake.CheckpointFile(filename, "r", comm=group_comm) as chk:
        mesh2 = chk.load_mesh(mesh_name)
        f2 = chk.load_function(mesh2, "probe")
else:
    with firedrake.CheckpointFile(filename, "r") as chk:
        mesh2 = chk.load_mesh(mesh_name)
        f2 = chk.load_function(mesh2, "probe")

print(f"[rank {rank}] sum={float(f2.dat.data.sum())}", flush=True)
'''
