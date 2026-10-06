# ==============================================================================
# @des: Worker process for test_mpi_failure_handling.py. Not a pytest test
#       itself -- launched as a real multi-rank job via mpirun/mpiexec by
#       that test, to exercise the actual MPI failure pattern
#       (icesee_da_full_parallel.py's exception handler) in isolation from
#       any particular model driver.
#
#       Rank 1 raises deliberately. Every rank follows the same pattern the
#       driver's exception handler uses: print full diagnostics (never
#       suppress the original exception), then comm.Abort() -- never a
#       further collective that assumes every rank is still reachable.
#       Healthy ranks (here, rank 0) are simulated reaching a later
#       Barrier(), matching a real driver's per-timestep synchronization;
#       this proves Abort() prevents them from hanging there forever
#       waiting for the rank that never arrives.
# ==============================================================================
import sys
import traceback

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

try:
    if rank == 1:
        raise RuntimeError("intentional test failure on rank 1")
    # Healthy ranks proceed toward whatever the next real synchronization
    # point would be. Without Abort() on the failing rank, this Barrier()
    # would hang forever once rank 1 never arrives.
    comm.Barrier()
    print(f"[rank {rank}] reached the barrier unexpectedly (rank 1 should "
          "have aborted the job before this)", flush=True)
except Exception as e:
    tb_str = "".join(traceback.format_exception(*sys.exc_info()))
    print(f"[rank {rank}] Fatal error: {e}\n{tb_str}", flush=True)
    comm.Abort(1)
