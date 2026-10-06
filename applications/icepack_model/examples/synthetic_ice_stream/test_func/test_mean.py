# test_mean.py
import numpy as np, h5py, os
from mpi4py import MPI

class Dummy:
    pass

def build_world(comm, nd=61615, nens=54, nt=1):
    """Partition rows block-wise across ranks."""
    size = comm.Get_size()
    rank = comm.Get_rank()
    rows_per = nd // size
    extra = nd % size
    start = rank * rows_per + min(rank, extra)
    stop  = start + rows_per + (1 if rank < extra else 0)
    return start, stop

def create_synthetic(comm, path, nd, nens):
    rank = comm.Get_rank()
    if rank == 0:
        # dataset per timestep/batch, shape (nd, nens)
        if os.path.exists(path): os.remove(path)
        with h5py.File(path, "w") as f:
            d = f.create_dataset("X0", shape=(nd, nens), dtype="f8", chunks=(min(nd,4096), 1))
            # Deterministic pattern: row i, ens j
            i = np.arange(nd)[:,None]
            j = np.arange(nens)[None,:]
            d[...] = (i * 0.01) + (j * 1.0)  # easy to sanity-check
    comm.Barrier()

def test_parallel_mean():
    from pathlib import Path
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    nd, nens, nt = 10007, 17, 3
    base = "tmp_mean_test"
    Path(base).mkdir(exist_ok=True)
    src_path = f"{base}/synthetic.h5"

    # Create synthetic per-timestep dataset list: here one dataset for simplicity.
    create_synthetic(comm, src_path, nd, nens)

    # Build a dummy "self"
    self = Dummy()
    self.mpi_comm = comm
    self.nd = nd
    self.nens = nens
    self.nt = nt
    self.base_path = base
    self.file_prefix = "forecast"
    self.current_batch_start = 0

    nd0, nd1 = build_world(comm, nd, nens, nt)
    self.nd_start_world, self.nd_end_world = nd0, nd1

    # Open source file and mimic self.datasets[batch_idx]
    fsrc = h5py.File(src_path, "r", driver="mpio", comm=comm)
    self.datasets = [fsrc["X0"]]  # batch 0 → timestep k=0

    # Import the function from your module or paste it above.
    from test_func import compute_forecast_mean_v3, compute_forecast_mean_chunked_v2  # replace import as needed
    # compute_forecast_mean_v3(self, k=0, flag="initial")
    compute_forecast_mean_chunked_v2(self, k=0, flag="initial")

    # Rank 0: verify
    comm.Barrier()
    ok = True
    if rank == 0:
        with h5py.File(f"{self.base_path}/{self.file_prefix}_mean.h5", "r") as f:
            col = f["mean"][:, 0]                      # parallel result
        with h5py.File(src_path, "r") as f:
            X = f["X0"][...]                           # (nd, nens)
        ref = X.mean(axis=1, dtype=np.float64)         # numpy reference
        ok = np.allclose(col, ref, rtol=0, atol=0)     # bitwise identity expected
        print("allclose:", ok, "| max abs diff:", np.max(np.abs(col - ref)))
    ok = comm.bcast(ok, root=0)
    assert ok, "Parallel mean mismatch vs numpy reference."

if __name__ == "__main__":
    test_parallel_mean()
