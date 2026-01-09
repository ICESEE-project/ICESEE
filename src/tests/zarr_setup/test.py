import h5py
from mpi4py import MPI
import numpy as np
import zarr
import numcodecs
import os
import gc

from forcast_class import EnKFIO  # Import the class above

def forecast(state):
    return state  # Placeholder: modify state as needed

def analyze(state, H, d, ens_idx, enkf_io, params=None):
    # Step 1: Compute important matrices
    start = MPI.Wtime()
    ens_perturbations_local = state - np.mean(state, axis=0)  # (nd_local, 1)
    Eta_local = np.dot(H, ens_perturbations_local)  # (m, 1)
    D_local = d + Eta_local  # (m, 1)
    HA_local = np.dot(H, state)  # (m, 1)
    Dprime_local = D_local - HA_local  # (m, 1)
    HAprime_local = HA_local - np.mean(HA_local, axis=0)  # (m, 1)

    # Write matrices to file
    # enkf_io.write_matrix(t, 'Eta', Eta_local, ens_idx)
    # enkf_io.write_matrix(t, 'HA', HA_local, ens_idx)
    # enkf_io.write_matrix(t, 'Dprime', Dprime_local, ens_idx)
    # enkf_io.write_matrix(t, 'HAprime', HAprime_local, ens_idx)

    # enkf_io.mpi_comm.Barrier()

    # # Step 2: Rank 0 computes X5
    # X5 = None
    # if enkf_io.mpi_comm.Get_rank() == 0:
    #     m_obs = d.shape[0]
    #     nrmin = min(m_obs, enkf_io.nens)

    #     # Gather full matrices
    #     Eta = enkf_io.gather_matrix(t, 'Eta')  # (m, nens)
    #     HA = enkf_io.gather_matrix(t, 'HA')  # (m, nens)
    #     Dprime = enkf_io.gather_matrix(t, 'Dprime')  # (m, nens)
    #     HAprime = enkf_io.gather_matrix(t, 'HAprime')  # (m, nens)

    #     # Compute HA' + eta
    #     HAprime_eta = HAprime + Eta

    #     # Compute SVD
    #     U, sig, _ = np.linalg.svd(HAprime_eta, full_matrices=False)
    #     sig = sig**2

    #     # Compute significant eigenvalues
    #     sigsum = np.sum(sig[:nrmin])
    #     sigsum1 = 0.0
    #     nrsigma = 0
    #     for i in range(nrmin):
    #         if sigsum1 / sigsum < 0.999:
    #             nrsigma += 1
    #             sigsum1 += sig[i]
    #             sig[i] = 1.0 / sig[i]
    #         else:
    #             sig[i:nrmin] = 0.0
    #             break

    #     # Compute X1 = sig * U^T
    #     X1 = np.empty((nrmin, m_obs))
    #     for j in range(m_obs):
    #         for i in range(nrmin):
    #             X1[i, j] = sig[i] * U[j, i]

    #     # Compute X2 = X1 * Dprime
    #     X2 = np.dot(X1, Dprime)

    #     # Compute X3 = U * X2
    #     X3 = np.dot(U, X2)

    #     # Compute X4 = HAprime.T * X3
    #     X4 = np.dot(HAprime.T, X3)

    #     # Compute X5 = X4 + I
    #     X5 = X4 + np.eye(enkf_io.nens)

    #     # Verify X5 column sums
    #     if not np.allclose(np.sum(X5, axis=0), 1.0):
    #         print(f"Sum of each X5 column is not 1.0: {np.sum(X5, axis=0)}")

    #     del Eta, HA, Dprime, HAprime, HAprime_eta, U, sig, X1, X2, X3, X4
    #     gc.collect()

    # # Step 3: Broadcast X5 to all ranks
    # X5 = enkf_io.mpi_comm.bcast(X5, root=0)

    # # Step 4: Each rank reads its portion of the state and computes analysis
    # ens_vec_local = state  # (nd_local, 1)
    # ens_vec_global = np.zeros((enkf_io.nd, enkf_io.nens), dtype='f8')
    # counts = enkf_io.mpi_comm.allgather(enkf_io.nd_local)
    # displacements = enkf_io.mpi_comm.allgather(enkf_io.nd_start)
    # enkf_io.mpi_comm.Allgatherv(ens_vec_local, [ens_vec_global, counts, displacements, MPI.DOUBLE])
    # analysis_vec = np.dot(ens_vec_global, X5)  # (nd, nens)

    # # Step 5: Each rank writes its portion of the analysis
    # analysis_vec_local = analysis_vec[enkf_io.nd_start:enkf_io.nd_end, ens_idx]
    # enkf_io.write_analysis(t, analysis_vec_local, ens_idx)

    # # Apply inflation if provided
    # if params and 'inflation_factor' in params:
    #     mean_params = np.mean(analysis_vec_local, axis=1)
    #     perturbations = analysis_vec_local - mean_params.reshape(-1, 1)
    #     inflated_perturbations = perturbations * params['inflation_factor']
    #     analysis_vec_local = mean_params.reshape(-1, 1) + inflated_perturbations
    #     enkf_io.write_analysis(t, analysis_vec_local, ens_idx)

    # # Check for negative thickness if 'h' is in vec_inputs
    # if params and 'vec_inputs' in params:
    #     for i, var in enumerate(params.get('vec_inputs', [])):
    #         if var == 'h':
    #             start = i * (enkf_io.nd_local // params.get('num_state_vars', 1))
    #             end = start + (enkf_io.nd_local // params.get('num_state_vars', 1))
    #             analysis_vec_local[start:end, :] = np.maximum(analysis_vec_local[start:end, :], 1e-2)
    #             enkf_io.write_analysis(t, analysis_vec_local, ens_idx)

    # del X5, ens_vec_global, analysis_vec
    # gc.collect()

    # return analysis_vec_local
    return None

class Model:
    def __init__(self):
        self.params = {
            "number_obs_instants": 10,
            "joint_estimation": True,
            "total_state_param_vars": 5,
            "num_state_vars": 3
        }

    def H_matrix(self, n_model, zarr_path="H_matrix.zarr"):
        """Observation operator matrix, saved to a Zarr file.

        Args:
            n_model (int): Size of the model state.
            zarr_path (str): Path to save the Zarr file.

        Returns:
            np.ndarray: The H matrix.
        """
        n = n_model

        # Initialize the H matrix
        H = np.zeros((self.params["number_obs_instants"] * 2 + 1, n))

        # Calculate distance between measurements
        di = int((n - 2) / (2 * self.params["number_obs_instants"]))

        # Fill the H matrix
        for i in range(1, self.params["number_obs_instants"] + 1):
            H[i - 1, i * di - 1] = 1
            H[self.params["number_obs_instants"] + i - 1, int((n - 2) / 2) + i * di - 1] = 1

        H[self.params["number_obs_instants"] * 2, n - 2] = 1  # Final element

        # Check if we have parameter estimation
        if self.params.get('joint_estimation', False):
            ndim = n // self.params["total_state_param_vars"]
            state_variables_size = ndim * self.params["num_state_vars"]
            num_params_size = n - state_variables_size
            H_param = np.zeros(num_params_size)
            H[:, state_variables_size:] = H_param

        # Ensure the output directory exists
        output_dir = os.path.dirname(zarr_path)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir, exist_ok=True)
            except Exception as e:
                print(f"Error creating output directory {output_dir}: {e}")
                raise

        # Save the H matrix to a Zarr file
        try:
            # Define chunk size for efficient storage (adjust based on matrix size)
            chunk_size = (min(1000, H.shape[0]), min(1000, H.shape[1]))
            
            # Create a Zarr array with compression, using Zarr format 2
            zarr_array = zarr.open(
                zarr_path,
                mode='w',
                shape=H.shape,
                chunks=chunk_size,
                dtype=H.dtype,
                zarr_format=2,  # Explicitly use Zarr format 2
                compressor=numcodecs.Blosc(cname='zstd', clevel=5, shuffle=numcodecs.Blosc.SHUFFLE)
            )
            
            # Write the H matrix to the Zarr array
            zarr_array[:] = H

            # clean up H matrix from memory
            del H; gc.collect()
            
        except Exception as e:
            print(f"Error saving H matrix to Zarr file: {e}")
            raise

        return None

    def Eta_matrix(self, state):
        # pseudo code for Eta matrix computation (nd >> m, nens)
        # read H (m,nd) from zarr
        # ens_pertubations = state - np.mean(state, axis=1).reshape(-1, 1) # (nd, nens)
        # Eta = np.dot(H, ens_pertubations) # (m, nens)
        pass

# if __name__ == "__main__":

# Example parameters
nd = 36231
nens = 10
nt = 50
tobserve = [0, 5, 10, 15, 19, 24, 30, 50, 60, 130, 450]
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
serial_file_creation = True

start_time = MPI.Wtime()

if nens >= comm.Get_size():
    subcomm_size = min(comm.Get_size(), nens)
    color = comm.Get_rank() % subcomm_size
    key = comm.Get_rank() // subcomm_size
    rounds = (nens + subcomm_size - 1) // subcomm_size
else:
    color = comm.Get_rank() % nens if comm.Get_rank() < nens else MPI.UNDEFINED
    key = comm.Get_rank() // nens
    rounds = 1

subcomm = comm.Split(color, key)
enkf_io = EnKFIO('enkf', nd, nens, nt, tobserve, subcomm, comm, serial_file_creation, batch_size=100)


model = Model()

# Dummy observation operator and data
# m_obs = 100  # Example observation dimension
m_obs = len(tobserve) 
# H = np.random.rand(m_obs, nd)  # Example observation operator
d = np.random.rand(m_obs, 1)  # Example observation vector
# d = np.random.rand(m_obs,)
params = {
    'inflation_factor': 1.05,
    'vec_inputs': ['h'],
    'num_state_vars': 1
}

for t in range(nt):
    for round_idx in range(rounds):
        ens_id = color + round_idx * subcomm_size
        if ens_id < nens:
            state = enkf_io.read_forecast(t, ens_id)
            # print(f"Rank {subcomm.Get_rank()}, time {t}, ens_id {ens_id}, state shape: {state.shape}")
            state = forecast(state)
            enkf_io.write_forecast(t + 1 if t < nt - 1 else t, state, ens_id)
            
    if t in tobserve:

        # generate the H file
        if rank == 0:
            print("Generating H matrix and saving to Zarr...")
            model.H_matrix(n_model=nd, zarr_path="output/H_matrix.zarr")
        comm.Barrier()

        # if nens >= size
        if nens >= comm.Get_size():
            for round_idx in range(rounds):
                ens_id = color + round_idx * subcomm_size
                if ens_id < nens:
                    state = enkf_io.read_analysis(t, ens_id)
                    print(f"Rank {subcomm.Get_rank()}, Analysis read at time {t}, ens_id {ens_id}, state shape: {state.shape}")
                    # state = analyze(state, H, d, ens_id, enkf_io, params)
                    # enkf_io.write_analysis(t, state, ens_id)
        # state = enkf_io.read_analysis(t, ens_id)
        # print(f"Rank {subcomm.Get_rank()}, Analysis read at time {t}, state shape: {state.shape}")

        # state = analyze(state, H, d, ens_id, enkf_io, params)
        # enkf_io.write_analysis(t, state, ens_id)

enkf_io.close()

end_time = MPI.Wtime()
walltime = MPI.COMM_WORLD.reduce(end_time - start_time, op=MPI.MAX, root=0)
exec_time = MPI.COMM_WORLD.reduce(end_time - start_time, op=MPI.SUM, root=0)
if MPI.COMM_WORLD.Get_rank() == 0:
    print(f"Total execution time: {exec_time:.2f} seconds, Wall time: {walltime:.2f} seconds")
    print(f"Test completed successfully on {MPI.COMM_WORLD.Get_size()} ranks.")
