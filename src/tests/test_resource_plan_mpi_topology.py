# ==============================================================================
# @des: Real-MPI subprocess tests for Stage 4's hierarchical resource
# topology, exercising the actual ParallelManager.icesee_mpi_init() ->
# resource_plan.plan_resources() -> MPI.Comm.Split() chain under real
# MPI.COMM_WORLD (not mocked), via _topology_probe_worker.py. No
# Firedrake/ISSM dependency -- a synthetic, permissive test model.
#
# Complements test_resource_plan.py's pure-Python planner tests (which
# verify the *plan* is correct without launching MPI) by verifying the
# *real communicators* built from that plan actually have the sizes,
# membership, and COMM_NULL/spare behavior the plan promises, and that
# every rank -- including spare ones -- terminates cleanly.
# ==============================================================================
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

_REPO_ROOT = Path(__file__).resolve().parents[2]
_WORKER = Path(__file__).resolve().parent / "parallel_mpi" / "_topology_probe_worker.py"

_RESULT_RE = re.compile(
    r"RESULT rank=(?P<rank>\d+) world_size=(?P<world_size>\d+) "
    r"is_spare=(?P<is_spare>True|False) group_id=(?P<group_id>\S+) "
    r"group_size=(?P<group_size>\d+) rank_in_group=(?P<rank_in_group>\S+) "
    r"group_collective_sum=(?P<group_collective_sum>\S+) "
    r"ranks_per_model=(?P<ranks_per_model>\d+) num_groups=(?P<num_groups>\d+) "
    r"rounds=(?P<rounds>\d+) spare_ranks=(?P<spare_ranks>\d+) "
    r"round0_member=(?P<round0_member>\S+)"
)


def _run_topology_probe(tmp_path, world_size, nens, ranks_per_model):
    mpirun = find_compatible_mpi_launcher()
    if mpirun is None:
        pytest.skip("no compatible mpirun available in this environment")

    env = dict(os.environ)
    env["PYTHONPATH"] = f"{_REPO_ROOT.parent}:{_REPO_ROOT}"

    result = subprocess.run(
        [
            mpirun, "--oversubscribe", "-n", str(world_size),
            sys.executable, str(_WORKER), str(nens), str(ranks_per_model), str(tmp_path),
        ],
        capture_output=True,
        text=True,
        timeout=60,
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
    return results, done_ranks


@pytest.mark.parametrize(
    "world_size,nens,ranks_per_model,exp_groups,exp_rounds,exp_spare",
    [
        (4, 4, 1, 4, 1, 0),
        (8, 4, 2, 4, 1, 0),
        (10, 4, 2, 4, 1, 2),
        (4, 5, 2, 2, 3, 0),
    ],
)
def test_real_mpi_topology_matches_plan(
    tmp_path, world_size, nens, ranks_per_model, exp_groups, exp_rounds, exp_spare
):
    results, done_ranks = _run_topology_probe(tmp_path, world_size, nens, ranks_per_model)

    # Every launched rank reported a result and terminated cleanly.
    assert set(results) == set(range(world_size))
    assert done_ranks == set(range(world_size))

    for rank, info in results.items():
        assert int(info["world_size"]) == world_size
        assert int(info["ranks_per_model"]) == ranks_per_model
        assert int(info["num_groups"]) == exp_groups
        assert int(info["rounds"]) == exp_rounds
        assert int(info["spare_ranks"]) == exp_spare

    spare_ranks = {r for r, info in results.items() if info["is_spare"] == "True"}
    active_ranks = set(results) - spare_ranks
    assert len(spare_ranks) == exp_spare
    assert len(active_ranks) == exp_groups * ranks_per_model

    # Spare ranks got no group and no group-local collective result (their
    # subcomm is COMM_NULL -- confirmed by never entering the collective).
    for rank in spare_ranks:
        assert results[rank]["group_id"] == "spare"
        assert results[rank]["group_collective_sum"] == "None"
        assert results[rank]["round0_member"] == "None"

    # Every active rank's group-local collective sum equals exactly its
    # own group's size -- proof the group communicator contains exactly
    # (and only) its own ranks_per_model members, no more, no less.
    for rank in active_ranks:
        info = results[rank]
        assert int(info["group_size"]) == ranks_per_model
        assert int(info["group_collective_sum"]) == ranks_per_model

    # Block-contiguous membership: ranks [g*R, (g+1)*R) belong to group g.
    for rank in active_ranks:
        expected_group = rank // ranks_per_model
        assert int(results[rank]["group_id"]) == expected_group


def test_spare_ranks_do_not_enter_group_collectives_and_all_ranks_terminate(tmp_path):
    # P=10, Nens=4, R=2: the exact known Stage-4 P=10/Nens=4 topology --
    # verifies the topology layer itself (not the HDF5 layer, which is
    # explicitly out of scope for this stage) completes cleanly for every
    # one of the 10 ranks, including both spares.
    results, done_ranks = _run_topology_probe(tmp_path, world_size=10, nens=4, ranks_per_model=2)
    assert done_ranks == set(range(10))
    spare_ranks = {r for r, info in results.items() if info["is_spare"] == "True"}
    assert spare_ranks == {8, 9}
    for rank in spare_ranks:
        assert results[rank]["group_collective_sum"] == "None"
