from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Fixed test dimensions and seed
n_state = 12
n_obs = 4
Nens = 3
np.random.seed(42)  # Ensure consistency across ranks

# ---- Step 1: Create and split ensemble_vec and H ----
# Serial full arrays (used only on rank 0 for validation)
if rank == 0:
    full_ensemble_vec = np.random.rand(n_state, Nens)
    full_H = np.random.rand(n_obs, n_state)
else:
    full_ensemble_vec = None
    full_H = None

# Compute counts and displacements
counts = [n_state // size + (1 if i < n_state % size else 0) for i in range(size)]
displs = np.cumsum([0] + counts[:-1])
local_n_state = counts[rank]
start_row = displs[rank]

# Allocate local buffers
local_ensemble = np.empty((local_n_state, Nens))
local_H = np.empty((n_obs, local_n_state))

# Scatter chunks of ensemble_vec and H to each rank
for j in range(Nens):
    column = np.empty(local_n_state)
    if rank == 0:
        column_data = [full_ensemble_vec[displs[i]:displs[i]+counts[i], j] for i in range(size)]
    else:
        column_data = None
    comm.Scatterv(column_data, column, root=0)
    local_ensemble[:, j] = column

for i in range(n_obs):
    row = np.empty(local_n_state)
    if rank == 0:
        row_data = [full_H[i, displs[r]:displs[r]+counts[r]] for r in range(size)]
    else:
        row_data = None
    comm.Scatterv(row_data, row, root=0)
    local_H[i, :] = row

# ---- Step 2: Compute local sum and global mean of ensemble_vec ----
local_sum = np.sum(local_ensemble, axis=1, keepdims=True)  # shape (local_n_state, 1)
global_mean = np.empty((local_n_state, 1))
comm.Allreduce(local_sum, global_mean, op=MPI.SUM)
global_mean /= Nens

# ---- Step 3: Compute local perturbations ----
local_perturb = local_ensemble - global_mean  # shape (local_n_state, Nens)

# ---- Step 4: Local matrix multiplication ----
local_Eta = np.dot(local_H, local_perturb)  # shape (n_obs, Nens)

# ---- Step 5: Global reduction ----
Eta = np.empty_like(local_Eta)
comm.Allreduce(local_Eta, Eta, op=MPI.SUM)

# ---- Step 6: Serial reference (on rank 0 only) ----
if rank == 0:
    ensemble_mean_serial = np.mean(full_ensemble_vec, axis=1, keepdims=True)
    perturb_serial = full_ensemble_vec - ensemble_mean_serial
    Eta_serial = np.dot(full_H, perturb_serial)

    # Compare
    print("Eta shape:", Eta.shape)
    diff = Eta - Eta_serial
    print("Diff:\n", np.round(diff, 10))
    print("Max abs diff:", np.max(np.abs(diff)))





# from mpi4py import MPI
# import numpy as np

# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()

# # Dimensions (can be very large in practice)
# n_state = 300000    # rows of ensemble_vec, cols of H
# n_obs = 1000       # rows of H, rows of Eta
# Nens = 1000         # ensemble size

# # Distribute n_state rows across processes
# counts = [n_state // size + (1 if i < n_state % size else 0) for i in range(size)]
# displs = np.cumsum([0] + counts[:-1])
# local_n_state = counts[rank]
# start_row = displs[rank]

# # Each rank generates its local ensemble_vec and H block
# local_ensemble = np.random.rand(local_n_state, Nens)
# local_H = np.random.rand(n_obs, local_n_state)

# # Step 1: Compute local mean
# local_mean = np.mean(local_ensemble, axis=1, keepdims=True)  # shape (local_n_state, 1)

# # Step 2: Allreduce to get global mean row-wise
# # Need to sum local means and divide by Nens globally
# local_sum = np.sum(local_ensemble, axis=1, keepdims=True)  # (local_n_state, 1)

# # Gather global sum (aligned on full n_state)
# global_sum = None
# if rank == 0:
#     global_sum = np.empty((n_state, 1))
# comm.Gatherv(local_sum, (global_sum, counts, displs, MPI.DOUBLE), root=0)

# if rank == 0:
#     global_mean = global_sum / Nens
# else:
#     global_mean = np.empty((n_state, 1))

# comm.Bcast(global_mean, root=0)

# # Extract local piece of global mean
# local_mean = global_mean[start_row:start_row+local_n_state]

# # Step 3: Compute local perturbations
# local_perturb = local_ensemble - local_mean  # shape: (local_n_state, Nens)

# # Step 4: Local matrix multiplication (H_i × X_i)
# local_Eta = np.dot(local_H, local_perturb)  # shape: (n_obs, Nens)

# # Step 5: Global sum of all contributions
# Eta = np.empty_like(local_Eta)
# comm.Allreduce(local_Eta, Eta, op=MPI.SUM)

# # Final result: Eta of shape (n_obs, Nens)
# if rank == 0:
#     print("Eta shape:", Eta.shape)