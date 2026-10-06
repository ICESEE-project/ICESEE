from mpi4py import MPI
import zarr
import os
import numpy as np

# Persistent arrays
eta_local = "analysis_data/eta_local.zarr"
ha_local = "analysis_data/ha_local.zarr"
dprime_local = "analysis_data/dprime_local.zarr"
ha_prime_local = "analysis_data/ha_prime_local.zarr"
zarr_files = [eta_local, ha_local, dprime_local, ha_prime_local]
zarr_arrays = {}

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Parameters
nd = 100000  # State dimension
nens = 100  # Number of ensemble members
nt = 1      # Total number of time steps
tobserve = list(range(0, nt, 10))  # Observation time steps (every 10 steps)
nt_obs = len(tobserve)

# Divide nd among ranks
nd_local_base = nd // size
remainder = nd % size

# Assign extra row to first `remainder` ranks
if rank < remainder:
    nd_local = nd_local_base + 1
    nd_start = rank * (nd_local_base + 1)
else:
    nd_local = nd_local_base
    nd_start = remainder * (nd_local_base + 1) + (rank - remainder) * nd_local_base

nd_end = nd_start + nd_local

# Let one rank create the Zarr arrays
if rank == 0:
    try:
        # Create directory if it doesn't exist, remove if it exists
        path = "analysis_data"
        if os.path.exists(path):
            os.system(f"rm -rf {path}")
        os.makedirs(path, exist_ok=True)
        print(f"Rank {rank}: Created 'analysis_data' directory.")

        # Create Zarr arrays with the specified shape and chunks
        zarr.create_array(store=eta_local, shape=(nd, nens), chunks=(nd_local, 1), dtype='f8')
        zarr.create_array(store=ha_local, shape=(nd, nens), chunks=(nd_local, 1), dtype='f8')
        zarr.create_array(store=dprime_local, shape=(nd, nens), chunks=(nd_local, 1), dtype='f8')
        zarr.create_array(store=ha_prime_local, shape=(nd, nens), chunks=(nd_local, 1), dtype='f8')
        print(f"Rank {rank}: Successfully created Zarr arrays.")
    except Exception as e:
        print(f"Rank {rank}: Failed to create Zarr arrays: {e}")

comm.Barrier()  # Ensure all ranks wait until the arrays are created
print(f"Rank {rank}: Passed creation barrier, opening Zarr arrays.")

# Open the Zarr arrays in read/write mode
try:
    zarr_arrays['eta'] = zarr.open_array(eta_local, mode='r+')
    zarr_arrays['ha'] = zarr.open_array(ha_local, mode='r+')
    zarr_arrays['dprime'] = zarr.open_array(dprime_local, mode='r+')
    zarr_arrays['ha_prime'] = zarr.open_array(ha_prime_local, mode='r+')
    print(f"Rank {rank}: Successfully opened Zarr arrays.")
except Exception as e:
    print(f"Rank {rank}: Failed to open Zarr arrays: {e}")

# Write data to Zarr arrays
try:
    # Generate random data for the local chunk (replace with actual data if needed)
    # local_shape = (nd_local, nens)
    local_shape = (nd_local,1)
    eta_chunk = np.random.rand(*local_shape)
    ha_chunk = np.random.rand(*local_shape)
    dprime_chunk = np.random.rand(*local_shape)
    ha_prime_chunk = np.random.rand(*local_shape)

    # Write to Zarr arrays
    try:
        # zarr_arrays['eta'][nd_start:nd_end, :] = eta_chunk
        # zarr_arrays['ha'][nd_start:nd_end, :] = ha_chunk
        # zarr_arrays['dprime'][nd_start:nd_end, :] = dprime_chunk
        # zarr_arrays['ha_prime'][nd_start:nd_end, :] = ha_prime_chunk
        ens_idx = 1
        zarr_arrays['eta'][nd_start:nd_end, ens_idx] = eta_chunk
        zarr_arrays['ha'][nd_start:nd_end, ens_idx] = ha_chunk
        zarr_arrays['dprime'][nd_start:nd_end, ens_idx] = dprime_chunk
        zarr_arrays['ha_prime'][nd_start:nd_end, ens_idx] = ha_prime_chunk
        print(f"Rank {rank}: Successfully wrote chunk (rows {nd_start}:{nd_end}) to Zarr arrays.")
    except Exception as e:
        print(f"Rank {rank}: Failed to write chunk to Zarr arrays: {e}")
    comm.Barrier()  # Ensure all ranks synchronize after writing
    print(f"Rank {rank}: Passed writing barrier.")
except Exception as e:
    print(f"Rank {rank}: Encountered an error during data generation or writing: {e}")