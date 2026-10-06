from mpi4py import MPI
import zarr
from zarr.storage import LocalStore
import numpy as np
import click
from pathlib import Path
import shutil
import logging
from contextlib import contextmanager
from time import perf_counter
import uuid
try:
    from numcodecs import Blosc
except ImportError:
    raise ImportError("numcodecs is required for Blosc compression. Install it with 'pip install numcodecs'.")

# Set up logging
logging.basicConfig(level=logging.INFO, format=f'Rank {MPI.COMM_WORLD.Get_rank()}: %(message)s')
logger = logging.getLogger(__name__)

# Log Zarr version for debugging
logger.info(f"Zarr version: {zarr.__version__}")

# Configuration
nd = 100000  # State dimension
nens = 100  # Number of ensemble members
nt = 500 # Total number of time steps (now number of observations to assimilate)
_shape = (nd, nens)
chunks = (nd // 10, 1)  # Adjusted chunk size for better distribution
stores = {
    'eta': 'eta.zarr',
    'ha': 'ha.zarr',
    'dprime': 'dprime.zarr',
    'ha_prime': 'ha_prime.zarr'
}
compressor = Blosc(cname='zstd', clevel=3, shuffle=2)

# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

@contextmanager
def zarr_array(store_path, mode='r+', shape=None, chunks=None, dtype='f8'):
    """Context manager for Zarr array to ensure proper handling."""
    try:
        array = zarr.open(
            store=store_path, 
            mode=mode, 
            shape=shape, 
            chunks=chunks, 
            dtype=dtype, 
            compressor=compressor, 
            fill_value=np.nan,
            zarr_format=2  # Explicitly use Zarr format 2
        )
        yield array
    except Exception as e:
        logger.error(f"Failed to open Zarr array at {store_path}: {e}")
        raise
    finally:
        if hasattr(array.store, 'close'):
            array.store.close()

def safe_delete_store(store_path):
    """Safely delete a Zarr store if it exists and is valid."""
    store_path = Path(store_path)
    if store_path.exists():
        try:
            if store_path.is_dir() and (store_path / '.zarr').exists():
                shutil.rmtree(store_path)
                logger.info(f"Deleted existing Zarr store at {store_path}")
            else:
                logger.warning(f"Path {store_path} is not a Zarr store, skipping deletion")
        except Exception as e:
            logger.error(f"Failed to delete Zarr store at {store_path}: {e}")
            raise

def create_zarr_stores():
    """Create Zarr stores on rank 0 and broadcast success."""
    success = False
    if rank == 0:
        try:
            for name, store_path in stores.items():
                safe_delete_store(store_path)
                with zarr_array(store_path, mode='w', shape=_shape, chunks=chunks) as array:
                    logger.info(f"Created Zarr store {name} at {store_path}")
            success = True
        except Exception as e:
            logger.error(f"Failed to create Zarr stores: {e}")
            raise
    success = comm.bcast(success, root=0)
    if not success:
        raise RuntimeError("Zarr store creation failed")
    comm.Barrier()

def get_rank_chunks(total_rows, rank, size):
    """Distribute rows across MPI ranks."""
    rows_per_rank = total_rows // size
    remainder = total_rows % size
    start = rank * rows_per_rank + min(rank, remainder)
    end = start + rows_per_rank + (1 if rank < remainder else 0)
    return start, end

def main():
    zarr_arrays = {}  # Initialize empty dictionary to avoid UnboundLocalError
    try:
        # Create Zarr stores
        create_zarr_stores()

        # Open Zarr arrays for all ranks
        for name, store_path in stores.items():
            zarr_arrays[name] = zarr.open(
                store=store_path, 
                mode='r+', 
                shape=_shape, 
                chunks=chunks, 
                dtype='f8', 
                compressor=compressor, 
                fill_value=np.nan,
                zarr_format=2  # Explicitly use Zarr format 2
            )
        
        # Distribute work across ranks
        start_row, end_row = get_rank_chunks(nd, rank, size)
        logger.info(f"Rank {rank} assigned rows {start_row} to {end_row}")

        # Initialize the state ensemble with random data
        np.random.seed(42)  # For reproducibility across ranks
        local_shape = (end_row - start_row, nens)
        local_X = np.random.rand(*local_shape)
        zarr_arrays['eta'][start_row:end_row, :] = local_X
        comm.Barrier()

        # Emulate EnKF analysis
        for i in range(nt):
            # Pick a random observation index (same on all ranks due to seed)
            np.random.seed(i)  # Different seed for each i
            obs_index = np.random.randint(0, nd)
            r_k = 1.0

            # Get h_j, bar_h, ha, var_h
            local_bar_h = 0.0
            local_var_h = 0.0
            local_ha = np.zeros(nens)
            owns_obs = (start_row <= obs_index < end_row)
            if owns_obs:
                local_obs_row = obs_index - start_row
                h_j = zarr_arrays['eta'][obs_index, :]
                local_bar_h = np.mean(h_j)
                local_ha = h_j - local_bar_h
                local_var_h = np.var(h_j, ddof=1)

            # Allreduce to share with all ranks
            bar_h = comm.allreduce(local_bar_h, op=MPI.SUM)
            var_h = comm.allreduce(local_var_h, op=MPI.SUM)
            ha = comm.allreduce(local_ha, op=MPI.SUM)

            # Generate synthetic observation y_k = bar_h + noise
            y_k = bar_h + np.random.randn() * np.sqrt(r_k)

            # Compute update
            local_X = zarr_arrays['eta'][start_row:end_row, :]
            local_mean = np.mean(local_X, axis=1)
            local_A = local_X - local_mean[:, None]
            local_cov = np.dot(local_A, ha) / (nens - 1)
            k = local_cov / (var_h + r_k)
            alpha = 1 / (1 + np.sqrt(r_k / (var_h + r_k)))
            tilde_k = k * alpha
            innovation = y_k - bar_h
            local_mean += k * innovation
            local_A -= tilde_k[:, None] * ha
            local_X = local_mean[:, None] + local_A

            # Write back updated state
            zarr_arrays['eta'][start_row:end_row, :] = local_X

            # For other arrays, keep original random write (using the shard logic)
            shard_start = i * (nd // 10)
            shard_end = (i + 1) * (nd // 10)
            if shard_end > start_row and shard_start < end_row:
                local_start = max(shard_start, start_row) - shard_start
                local_end = min(shard_end, end_row) - shard_start
                local_shape = (local_end - local_start, nens)

                # Generate random data for the local chunk
                ha = np.random.rand(*local_shape)
                dprime = np.random.rand(*local_shape)
                ha_prime = np.random.rand(*local_shape)

                # Write to Zarr arrays
                try:
                    zarr_arrays['ha'][shard_start + local_start:shard_start + local_end, :] = ha
                    zarr_arrays['dprime'][shard_start + local_start:shard_start + local_end, :] = dprime
                    zarr_arrays['ha_prime'][shard_start + local_start:shard_start + local_end, :] = ha_prime
                    logger.info(f"Rank {rank} wrote shard {i} rows {shard_start + local_start}:{shard_start + local_end} for auxiliary arrays")
                except Exception as e:
                    logger.error(f"Rank {rank} failed to write shard {i}: {e}")
                    raise

            logger.info(f"Rank {rank} completed EnKF update for observation {i}")

        # Synchronize after writes
        comm.Barrier()
        logger.info(f"Rank {rank} completed writing")

    except Exception as e:
        logger.error(f"Rank {rank} encountered error: {e}")
        raise
    finally:
        # Close Zarr stores
        for name, array in zarr_arrays.items():
            if hasattr(array.store, 'close'):
                array.store.close()

if __name__ == '__main__':
    start_time = MPI.Wtime()
    main()
    end_time = MPI.Wtime()
    elapsed_time = end_time - start_time
    wall_time = comm.reduce(elapsed_time, op=MPI.MAX, root=0)
    execution_time = comm.reduce(elapsed_time, op=MPI.SUM, root=0)
    if rank == 0:
        logger.info(f"Total execution time: {execution_time:.2f} seconds, Wall time: {wall_time:.2f} seconds")
        logger.info("Test completed successfully.")