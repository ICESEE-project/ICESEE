import h5py
from mpi4py import MPI
import numpy as np
import zarr
import numcodecs
import os
import gc

# from forcast_class import EnKFIO  # Import the class above
from forecast_analysis_i_o_class import EnKFIO
# from EnKF_analysis_io import EnKFIO
# from forecast_analysis_i_o_class_v0 import EnKFIO
# from assimilation_io_zarr import EnKFIO_zarr

def forecast(state):
    # create random values between 4 to 10
    # np.random.seed(42)
    state = np.random.rand(*state.shape) * np.sin(10) + 45
    #add noise with zero mean
    state += np.random.normal(0, 1, state.shape)
    return state  # Placeholder: modify state as needed

def analyze(state, H, d, ens_idx, enkf_io, params=None):
    return None

# Example parameters
nd = 1700
nens = 40
nt = 200
dt = 1.0
comm = MPI.COMM_WORLD
serial_file_creation = True

params = {
    "joint_estimation": False,
    "total_state_param_vars": 2,
    'inflation_factor': 1.05,
    'vec_inputs': ['h', 'v'],
    'num_state_vars': 1,
    'freq_obs': 2,
    'obs_start_time': 4,
    'obs_max_time': 8,
    'sig_obs' : [0.2, 0.2],
    "default_run": True

}

kwargs = {'nd': nd, 'nens': nens, 'nt': nt, 
          'comm': comm, 'comm_world': comm, 'params': params,
         'serial_file_creation': serial_file_creation,
         't': np.linspace(0, nt, int(nt/dt)+1),
         'dim_list': [nd for _ in range(comm.Get_size())],
         'vec_inputs': ['h', 'v']}
kwargs.update(params)

start_time = MPI.Wtime()

if nens >= comm.Get_size():
    subcomm_size = min(comm.Get_size(), nens)
    color = comm.Get_rank() % subcomm_size
    key = comm.Get_rank() // subcomm_size
    rounds = (comm.Get_size() + subcomm_size - 1) // subcomm_size
else:
    color = comm.Get_rank() % nens if comm.Get_rank() < nens else MPI.UNDEFINED
    key = comm.Get_rank() // nens
    rounds = 1

subcomm = comm.Split(color, key)
enkf_io = EnKFIO('enkf', nd, nens, nt, subcomm, comm, params, serial_file_creation, batch_size=100)
# enkf_io = EnKFIO_zarr('enkf', nd, nens, nt, tobserve, subcomm, comm, params, serial_file_creation, base_path="enkf_data", batch_size=100)

# Generate dummy statevec_true
statevec_true_chunk_size = (min(1000, nd//10), 1)
# intialize zarr array on root only to avoid race conditions
store_path = "output/statevec_true.zarr"
rank = comm.Get_rank()
if rank == 0:
    statevec_true = zarr.create_array(store_path, shape=(nd,nt), chunks=statevec_true_chunk_size, dtype='f8', overwrite=True)

comm.Barrier()
statevec_true = zarr.open(store_path, mode='r+')


size = comm.Get_size()

# # Distribute time steps (columns) across processes
# nt_per_process = nt // size 
# start_t = rank * nt_per_process
# end_t = start_t + nt_per_process if rank < size - 1 else nt

# np.random.seed(1234 + rank)

# # Generate and write random data for this process's portion
# for t in range(start_t, end_t):
#     # Generate random data for one time step (column)
#     data = np.random.rand(nd, 1)
#     # Write to Zarr array (Zarr handles chunk alignment)
#     statevec_true[:, t] = data[:, 0]

# # Synchronize to ensure all writes are complete
# comm.Barrier()

# # Optional: Verify the array (e.g., on rank 0)
# if rank == 0:
#     print(f"Zarr array created with shape {statevec_true.shape}, chunks {statevec_true.chunks}")

# Calculate total number of chunks
chunks_per_col = nd // (nd // 10)  # Number of chunks per column (10 in this case)
total_chunks = chunks_per_col * nt  # Total chunks in the array

# Distribute chunks across processes
chunks_per_process = total_chunks // size
remainder = total_chunks % size
start_chunk = rank * chunks_per_process + min(rank, remainder)
if rank < remainder:
    num_chunks = chunks_per_process + 1
else:
    num_chunks = chunks_per_process

# Calculate chunk indices for this process
chunk_indices = list(range(start_chunk, start_chunk + num_chunks))

# Set unique random seed for each process
np.random.seed(1234 + rank)

# Write data to assigned chunks
for chunk_idx in chunk_indices:
    # Convert global chunk index to (row_chunk, col_idx)
    row_chunk = chunk_idx // nt
    col_idx = chunk_idx % nt

    # Calculate row range for this chunk
    row_start = row_chunk * (nd // 10)
    row_end = row_start + (nd // 10)

    # Generate random data for this chunk
    data = np.random.rand(nd // 10, 1)

    # Write to the specific chunk
    statevec_true[row_start:row_end, col_idx] = data[:, 0]

# Synchronize to ensure all writes are complete
comm.Barrier()

# # Optional: Verify the array (e.g., on rank 0)
# if rank == 0:
#     print(f"Zarr array created with shape {statevec_true.shape}, chunks {statevec_true.chunks}")

# generate synthetic observations to get d = H \in (mxnd) @synthetic_obs \in (ndxkm)
synthetic_obs_zarr_path="output/synthetic_observations.zarr"
error_R_zarr_path="output/error_R.zarr"
kwargs.update({'synthetic_obs_zarr_path': synthetic_obs_zarr_path, 'error_R_zarr_path': error_R_zarr_path})
tobserve, m_obs = enkf_io._create_synthetic_observations(**kwargs)
kwargs.update({'tobserve': tobserve, 'm_obs': m_obs})

# generate the H file
rank = comm.Get_rank()
if rank == 0:
    print("Generating H matrix and saving to Zarr...")
    H_matrix_zarr_path = "output/H_matrix.zarr"
    kwargs.update({'H_matrix_zarr_path': H_matrix_zarr_path})
    enkf_io.H_matrix(**kwargs)
comm.Barrier()

time_mean = 0.0
time_mean_chunked = 0.0
analysis_time = 0.0
km = 0
for k in range(nt):
    for round_idx in range(rounds):
        ens_id = color + round_idx * subcomm_size
        if ens_id < nens:
            state = enkf_io.read_forecast(k, ens_id)
            print(f"Rank {subcomm.Get_rank()}, time {k}, ens_id {ens_id}, state shape: {state.shape}")
            state = forecast(state)
            enkf_io.write_forecast(k + 1 if k < nt - 1 else k, state, ens_id)

#     # compute the forecast mean
#     start_mean_time = MPI.Wtime()
    # enkf_io.compute_forecast_mean(k)
#     time_mean += MPI.Wtime() - start_mean_time
    comm.Barrier()
    # use by chunked via ensemble dimension
    start_mean_chunked_time = MPI.Wtime()
    enkf_io.compute_forecast_mean_chunked(k)
    # enkf_io.compute_forecast_mean_chunked_v2(k)
    # enkf_io.compute_forecast_mean_chunked_gather(k)
    time_mean_chunked += MPI.Wtime() - start_mean_chunked_time

    comm.Barrier()
    # result = enkf_io.compare_forecast_means(k, ens_chunk_size=1)
    # if rank == 0:
    #     print(f"Comparison result: {result}")

    if km < m_obs and k+1 == tobserve[km]:
    #     # compute Eta for @ ens_idx
    #     # kwargs.update({'Eta_matrix_zarr_path': "output/Eta_matrix.zarr"})
    #     # ens_idx = 1
    #     # enkf_io.Eta_matrix(km, ens_idx, **kwargs)

    #     # compute d for @ ens_idx
    #     kwargs.update({'d_matrix_zarr_path': "output/d_matrix.zarr"})
        # X5 = enkf_io.compute_X5(km, **kwargs)
    #     # comm.Barrier()
        analysis_time_start = MPI.Wtime()
        enkf_io.compute_analysis_update(km,**kwargs)
        analysis_time += MPI.Wtime() - analysis_time_start
        # comm.Barrier()

        # compute the analysis mean
        # comm.Barrier()
        # enkf_io.compute_forecast_mean_chunked(km)
        km += 1

enkf_io.close()

end_time = MPI.Wtime()
walltime = MPI.COMM_WORLD.reduce(end_time - start_time, op=MPI.MAX, root=0)
exec_time = MPI.COMM_WORLD.reduce(end_time - start_time, op=MPI.SUM, root=0)
# forecast_mean_time = MPI.COMM_WORLD.reduce(time_mean, op=MPI.MAX, root=0)
forecast_mean_chunked_time = MPI.COMM_WORLD.reduce(time_mean_chunked, op=MPI.MAX, root=0)
# forecast_mean_chunked_time = time_mean_chunked
analysis_time = MPI.COMM_WORLD.reduce(analysis_time, op=MPI.MAX, root=0)
if MPI.COMM_WORLD.Get_rank() == 0:
    print(f"Total execution time: {exec_time:.2f} seconds, Wall time: {walltime:.2f} seconds")
    # print(f"Forecast mean time: {forecast_mean_time:.2f} seconds")
    print(f"Forecast mean chunked time: {forecast_mean_chunked_time:.2f} seconds")
    print(f"Analysis time: {analysis_time:.2f} seconds")
    print(f"Test completed successfully on {MPI.COMM_WORLD.Get_size()} ranks.")
