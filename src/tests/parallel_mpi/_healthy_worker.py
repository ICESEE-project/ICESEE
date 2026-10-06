# ==============================================================================
# @des: Negative-control worker for test_mpi_failure_handling.py -- a
#       trivial healthy 2-rank job (no exception, one real Barrier) that
#       must exit 0 quickly. Exists so the failure test proves the harness
#       distinguishes FAIL from PASS, not merely that it always reports
#       nonzero.
# ==============================================================================
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
comm.Barrier()
print(f"[rank {rank}] completed normally", flush=True)
