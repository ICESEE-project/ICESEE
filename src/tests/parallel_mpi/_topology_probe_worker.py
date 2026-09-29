# ==============================================================================
# @des: Lightweight real-MPI worker for Stage 4 topology tests. Exercises
# ParallelManager.icesee_mpi_init/icesee_mpi_ens_distribution directly
# (real MPI.Comm.Split under real MPI.COMM_WORLD) with a synthetic,
# Firedrake/ISSM-free icesee_kwargs, and prints one machine-parseable line
# per rank so the parent test (test_resource_plan_mpi_topology.py) can
# verify communicator membership/size/spare status without needing to
# inspect MPI objects across process boundaries.
#
# Usage: mpirun -n <P> python _topology_probe_worker.py <nens> <ranks_per_model>
# ==============================================================================
import sys

from mpi4py import MPI

from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import (
    ParallelManager,
)
from ICESEE.src.parallelization.parallel_mpi.model_capabilities import (
    register_model_capabilities,
)

nens = int(sys.argv[1])
ranks_per_model = sys.argv[2]
ranks_per_model = None if ranks_per_model == "none" else int(ranks_per_model)

# Permissive synthetic model: this worker exists purely to exercise real
# communicator construction, not to validate any particular application's
# capability policy (that is model_capabilities.py's own, separately
# tested, concern).
register_model_capabilities("mpi-topology-test-model", supports_multi_rank_per_model=True)

world = MPI.COMM_WORLD
world_rank = world.Get_rank()
world_size = world.Get_size()

icesee_kwargs = {
    "Nens": nens,
    "ranks_per_model": ranks_per_model,
    "default_run": True,
    "sequential_run": False,
    "even_distribution": False,
    "execution_mode": 2,
    "data_path": sys.argv[3] if len(sys.argv) > 3 else "/tmp",
    "model_name": "mpi-topology-test-model",
    "verbose": False,
}

sub_rank, sub_size, subcomm, ens_id = ParallelManager().icesee_mpi_init(icesee_kwargs)

is_spare = subcomm is None or subcomm == MPI.COMM_NULL
if is_spare:
    group_id = "spare"
    group_size = 0
    rank_in_group = "spare"
else:
    group_id = icesee_kwargs["resource_plan"].group_id(world_rank)
    group_size = subcomm.Get_size()
    rank_in_group = subcomm.Get_rank()

# A group-local collective: every rank in a real (non-spare) group
# participates; the result must be exactly this group's own size, proving
# no rank outside the group leaked into it (and no rank inside it is
# missing).
if not is_spare:
    group_collective_sum = subcomm.allreduce(1, op=MPI.SUM)
else:
    group_collective_sum = None

plan = icesee_kwargs["resource_plan"]

print(
    f"RESULT rank={world_rank} world_size={world_size} is_spare={is_spare} "
    f"group_id={group_id} group_size={group_size} rank_in_group={rank_in_group} "
    f"group_collective_sum={group_collective_sum} "
    f"ranks_per_model={plan.ranks_per_model} num_groups={plan.num_model_groups} "
    f"rounds={plan.num_rounds} spare_ranks={plan.spare_ranks} "
    f"round0_member={ens_id}",
    flush=True,
)

# Every rank -- spare or not -- must reach this world-wide Barrier cleanly,
# proving spare ranks do not leave world-scoped control flow desynced from
# non-spare ranks even though they took the early-return path above.
world.Barrier()
print(f"DONE rank={world_rank}", flush=True)
