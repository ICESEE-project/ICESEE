# ==============================================================================
# @des: Real-MPI worker for test_node_local_model_groups.py. Fakes each
# rank's host name (rank // ranks_per_node) so a node-spanning model-group
# layout can be reproduced on one machine, then runs the real
# ParallelManager.icesee_mpi_init startup path and prints one line per rank.
#
# Usage: mpirun -n <P> python _node_span_probe_worker.py \
#            <nens> <ranks_per_model> <ranks_per_node> <requires_node_local 0|1> <data_path>
# ==============================================================================
import socket
import sys

from mpi4py import MPI

nens, ranks_per_model, ranks_per_node, requires = (int(a) for a in sys.argv[1:5])
data_path = sys.argv[5]

world = MPI.COMM_WORLD
world_rank = world.Get_rank()
socket.gethostname = lambda: f"node{world_rank // ranks_per_node}"

from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import (  # noqa: E402
    ParallelManager,
)
from ICESEE.src.parallelization.parallel_mpi.model_capabilities import (  # noqa: E402
    register_model_capabilities,
)

register_model_capabilities(
    "node-span-test-model",
    supports_multi_rank_per_model=True,
    requires_node_local_model_groups=bool(requires),
)

icesee_kwargs = {
    "Nens": nens,
    "ranks_per_model": ranks_per_model,
    "default_run": True,
    "sequential_run": False,
    "even_distribution": False,
    "execution_mode": 2,
    "data_path": data_path,
    "model_name": "node-span-test-model",
    "verbose": False,
}

try:
    ParallelManager().icesee_mpi_init(icesee_kwargs)
except RuntimeError as exc:
    print(f"REFUSED rank={world_rank} {exc}", flush=True)
else:
    # A collective after startup: completes only if every rank got here.
    world.Barrier()
    print(f"STARTED rank={world_rank}", flush=True)
