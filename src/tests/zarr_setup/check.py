import numpy as np
from mpi4py import MPI
import logging

class ForecastProcessor:
    def __init__(self, num_rows, num_cols, num_batches):
        self.mpi_comm = MPI.COMM_WORLD
        self.rank = self.mpi_comm.Get_rank()
        self.size = self.mpi_comm.Get_size()
        self.logger = logging.getLogger(__name__)
        
        # Create synthetic dataset: list of 2D arrays (num_batches, num_rows, num_cols)
        self.datasets = [np.arange(i * num_rows * num_cols, (i + 1) * num_rows * num_cols)
                        .reshape(num_rows, num_cols).astype('f8') for i in range(num_batches)]
        
        # Distribute rows across ranks
        rows_per_rank = num_rows // self.size
        remainder = num_rows % self.size
        self.nd_start_world = self.rank * rows_per_rank + min(self.rank, remainder)
        if self.rank < remainder:
            rows_per_rank += 1
        self.nd_end_world = self.nd_start_world + rows_per_rank
        self.current_batch_start = 0

    def _ensure_batch(self, t):
        """Ensure the batch for time t is loaded."""
        batch_idx = t - self.current_batch_start
        if batch_idx < 0 or batch_idx >= len(self.datasets):
            raise ValueError(f"Invalid batch index {batch_idx}")
        
    def compute_forecast_mean(self, t):
        self._ensure_batch(t)
        comm = self.mpi_comm
        rank = comm.Get_rank()
        size = comm.Get_size()

        batch_idx = t - self.current_batch_start
        start = MPI.Wtime()

        # Each rank holds a disjoint block of rows in [nd_start_world:nd_end_world)
        local_data = self.datasets[batch_idx][self.nd_start_world:self.nd_end_world, :]
        # Per-row mean across columns for the local rows
        local_mean = np.mean(local_data, axis=1).astype('f8', copy=False)  # shape: (nd_local,)

        # Gather the number of rows each rank is processing
        local_row_count = np.array([local_mean.shape[0]], dtype='i8')
        global_row_counts = np.zeros(size, dtype='i8') if rank == 0 else None
        comm.Gather(local_row_count, global_row_counts, root=0)

        # Compute displacements for Gatherv
        if rank == 0:
            displacements = np.zeros(size, dtype='i8')
            displacements[1:] = np.cumsum(global_row_counts[:-1])
            total_rows = np.sum(global_row_counts)
            global_mean = np.zeros(total_rows, dtype='f8')
        else:
            displacements = None
            global_mean = None

        # Use Gatherv to collect local means into global_mean on root
        comm.Gatherv(local_mean, [global_mean, global_row_counts, displacements, MPI.DOUBLE], root=0)

        # Broadcast the global mean to all ranks
        if rank == 0:
            result = global_mean
        else:
            result = np.zeros(np.sum(comm.bcast(global_row_counts if rank == 0 else None, root=0)), dtype='f8')
        comm.Bcast(result, root=0)

        # Log timing
        end = MPI.Wtime()
        if rank == 0:
            self.logger.debug(f"compute_forecast_mean took {end - start:.3f} seconds")

        return result

def test_compute_forecast_mean():
    # Configure logging
    logging.basicConfig(level=logging.DEBUG)
    
    # Initialize MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    # Test parameters
    num_rows = 10
    num_cols = 5
    num_batches = 2
    t = 0  # Test for batch index 0

    # Create processor instance
    processor = ForecastProcessor(num_rows, num_cols, num_batches)
    
    # Compute result
    result = processor.compute_forecast_mean(t)
    
    # Expected result: mean across columns for each row in datasets[0]
    if rank == 0:
        expected = np.mean(processor.datasets[0], axis=1)
        np.testing.assert_array_almost_equal(result, expected, decimal=6,
                                            err_msg="Computed mean does not match expected mean")
        print(f"Rank {rank}: Test passed! Result shape: {result.shape}")
    else:
        print(f"Rank {rank}: Result shape: {result.shape}")

if __name__ == "__main__":
    test_compute_forecast_mean()