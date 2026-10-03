# ==============================================================================
# @des: Regression tests for the Idealized PIG Mode-3 PACE hang (2026-10-03).
#
# A 40-rank, model_nprocs=2 job placed 14/13/13 ranks on three nodes, so one
# 2-rank model group (world ranks 26, 27) straddled two nodes. Firedrake's
# collective CheckpointFile open on such a group's communicator (parallel
# HDF5 over OMPIO, NFS) can deadlock inside MPI_File_open (inconsistent
# shared-file-pointer component selection across nodes), and every other
# rank then waits forever at the next world collective. (The same job was
# also exposed to the PyOP2 compiled-kernel cache deadlock covered by
# test_pyop2_compile_cache.py.)
#
# Icepack declares requires_node_local_model_groups, and startup refuses
# such a layout on every rank with one clear error instead of hanging. Models
# that do not declare it only get a warning; node-aligned layouts are
# unaffected.
# ==============================================================================
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from ICESEE.src.parallelization.parallel_mpi.model_capabilities import (
    get_model_capabilities,
    node_spanning_model_groups,
    validate_node_local_model_groups,
)
from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

_REPO_ROOT = Path(__file__).resolve().parents[2]
_WORKER = Path(__file__).resolve().parent / "parallel_mpi" / "_node_span_probe_worker.py"


def _hosts(*counts):
    return [f"node{i}" for i, count in enumerate(counts) for _ in range(count)]


def test_failing_pace_layout_has_a_node_spanning_group():
    # srun -n 40 in a 3-node allocation placed ranks 0-13 / 14-26 / 27-39:
    # group 13 = ranks 26 and 27 straddles the second and third node.
    spanning = node_spanning_model_groups(_hosts(14, 13, 13), lambda r: r // 2)
    assert spanning == {13: ["node1", "node2"]}


def test_node_aligned_layout_has_no_spanning_group():
    assert node_spanning_model_groups(_hosts(20, 20), lambda r: r // 2) == {}


def test_spare_ranks_are_ignored():
    # Ranks 4 and 5 are spares (group None) on another node.
    group = lambda r: r // 2 if r < 4 else None  # noqa: E731
    assert node_spanning_model_groups(_hosts(4, 2), group) == {}


class _FakeWorld:
    def __init__(self, hosts, rank=0):
        self.hosts = hosts
        self.rank = rank
        self.allgathers = 0

    def allgather(self, _value):
        self.allgathers += 1
        return list(self.hosts)

    def Get_rank(self):
        return self.rank


def test_icepack_declares_node_local_model_groups():
    assert get_model_capabilities("icepack").requires_node_local_model_groups is True


def test_icepack_node_spanning_layout_is_refused():
    world = _FakeWorld(_hosts(14, 13, 13))
    with pytest.raises(RuntimeError, match="requires node-local model groups"):
        validate_node_local_model_groups("icepack", world, lambda r: r // 2, 2)


def test_other_models_only_warn(capsys):
    world = _FakeWorld(_hosts(3, 3))
    spanning = validate_node_local_model_groups("lorenz", world, lambda r: r // 2, 2)
    assert spanning == {1: ["node0", "node1"]}
    assert "WARNING" in capsys.readouterr().out


def test_single_rank_groups_skip_the_check():
    world = _FakeWorld(_hosts(1, 1))
    assert validate_node_local_model_groups("icepack", world, lambda r: r, 1) == {}
    assert world.allgathers == 0


def _run_worker(tmp_path, ranks_per_node, requires):
    mpirun = find_compatible_mpi_launcher()
    if mpirun is None:
        pytest.skip("no compatible mpirun available in this environment")
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{_REPO_ROOT.parent}:{_REPO_ROOT}"
    # 6 ranks, Nens=3, 2 ranks/model -> groups {0,1}, {2,3}, {4,5}.
    result = subprocess.run(
        [mpirun, "--oversubscribe", "-n", "6", sys.executable, str(_WORKER),
         "3", "2", str(ranks_per_node), str(int(requires)), str(tmp_path)],
        capture_output=True, text=True, timeout=90, env=env,
    )
    return result


def test_real_mpi_node_spanning_icepack_like_layout_fails_fast_on_every_rank(tmp_path):
    result = _run_worker(tmp_path, ranks_per_node=3, requires=True)
    refused = [line for line in result.stdout.splitlines() if line.startswith("REFUSED")]
    assert len(refused) == 6, result.stdout + result.stderr
    assert all("group 1 on node0+node1" in line for line in refused)


def test_real_mpi_node_aligned_layout_starts(tmp_path):
    result = _run_worker(tmp_path, ranks_per_node=2, requires=True)
    started = [line for line in result.stdout.splitlines() if line.startswith("STARTED")]
    assert result.returncode == 0 and len(started) == 6, result.stdout + result.stderr
    assert "WARNING" not in result.stdout
