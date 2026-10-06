import numpy as np
import h5py
from mpi4py import MPI
import sys
import psutil
import time
import os

# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Dataset dimensions (adjust to match your setup)
nd = 70000  # State dimension (e.g., total_state_param_vars * hdim)
Nens = 100      # Number of ensembles
nt = 100        # Number of timesteps
dtype = np.float32

# Calculate local chunk sizes for state dimension
local_nd = nd // size
remainder = nd % size
if rank < remainder:
    local_nd += 1
local_offset = rank * (nd // size) + min(rank, remainder)

# Calculate local ensemble counts
ens_per_rank = Nens // size
ens_remainder = Nens % size
if rank < ens_remainder:
    ens_per_rank += 1
local_ens_offset = rank * (Nens // size) + min(rank, ens_remainder)
local_ens_count = ens_per_rank

# Generate dummy data
local_ensemble = np.random.rand(local_nd, local_ens_count).astype(dtype) if local_nd > 0 and local_ens_count > 0 else None
if rank == 0:
    ensemble_mean = np.random.rand(nd).astype(dtype)
else:
    ensemble_mean = np.zeros(nd, dtype=dtype)
ensemble_mean = comm.bcast(ensemble_mean, root=0)

# Output file
output_file = "test_large_ensemble.h5"  # Single shared file
# output_file = f"test_large_ensemble_{rank}.h5"  # Uncomment for per-rank files

# Redirect output to rank-specific log file
sys.stdout = open(f'log_rank_{rank}.txt', 'w')

# Diagnostics
process = psutil.Process()
start = time.time()
mem_info = process.memory_info()
print(f"[DEBUG] Rank: {rank} memory usage before: {mem_info.rss / 1024**2:.2f} MB")

# Write to HDF5
with h5py.File(output_file, 'w', driver='mpio', comm=comm) as f:
    f.atomic = True  # Enable collective metadata
    # Create datasets
    create_start = time.time()
    dset = f.create_dataset('ensemble', (nd, Nens, nt+1), dtype=dtype)
    ens_mean = f.create_dataset('ensemble_mean', (nd, nt+1), dtype=dtype)
    print(f"[DEBUG] Rank: {rank} created datasets in {time.time() - create_start} seconds")

    # Write ensemble data
    write_start = time.time()
    if local_ensemble is not None:
        dset[local_offset:local_offset + local_nd, local_ens_offset:local_ens_offset + local_ens_count, 0] = local_ensemble
    print(f"[DEBUG] Rank: {rank} wrote ensemble chunk ({local_nd}, {local_ens_count}) at offset ({local_offset}, {local_ens_offset}) in {time.time() - write_start} seconds")

    # Write ensemble_mean collectively
    write_start = time.time()
    local_ens_mean = ensemble_mean if rank == 0 else np.zeros(nd, dtype=dtype)
    ens_mean[:, 0] = local_ens_mean
    print(f"[DEBUG] Rank: {rank} wrote ensemble_mean in {time.time() - write_start} seconds")

# Final diagnostics
mem_info = process.memory_info()
print(f"[DEBUG] Rank: {rank} memory usage after: {mem_info.rss / 1024**2:.2f} MB")
print(f"[DEBUG] Rank: {rank} total time: {time.time() - start} seconds")
sys.stdout.close()
