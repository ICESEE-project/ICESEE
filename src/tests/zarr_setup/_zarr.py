import numpy as np
import zarr
import numcodecs
import os
import gc

from mpi4py import MPI

from forcast_class import EnKFIO

# Assuming `self.params` is defined with necessary keys
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

    def Eta_matrix(self):
        pass
        


# Main execution
if __name__ == "__main__":
    

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Example parameters
    nd = 36231
    nens = 100
    nt = 500
    tobserve = [0, 5, 10, 15, 19, 24, 30, 50, 60, 130, 450]
    comm = MPI.COMM_WORLD
    serial_file_creation = True

    start_time = MPI.Wtime()
    model = Model()

    # if rank == 0:
    #     print("Generating H matrix and saving to Zarr...")
    #     H = model.H_matrix(n_model=nd, zarr_path="output/H_matrix.zarr")

    #     # To read the Zarr file later
    #     # H_from_zarr = zarr.open("output/H_matrix.zarr", mode='r')[:]
    #     # print(f"Matrix equality check: {np.array_equal(H, H_from_zarr)}")  # Should print True
    # comm.Barrier()


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
    
    enkf_io = EnKFIO('enkf', nd, nens, nt, tobserve, subcomm, comm, serial_file_creation, batch_size=100)

    # # Dummy observation operator and data
    # # m_obs = 100  # Example observation dimension
    # m_obs = len(tobserve) 
    # H = np.random.rand(m_obs, nd)  # Example observation operator
    # d = np.random.rand(m_obs, 1)  # Example observation vector
    # # d = np.random.rand(m_obs,)
    # params = {
    #     'inflation_factor': 1.05,
    #     'vec_inputs': ['h'],
    #     'num_state_vars': 1
    # }

    # for t in range(nt):
    #     for round_idx in range(rounds):
    #         ens_id = color + round_idx * subcomm_size
    #         if ens_id < nens:
    #             state = enkf_io.read_forecast(t, ens_id)
    #             print(f"Rank {subcomm.Get_rank()}, time {t}, ens_id {ens_id}, state shape: {state.shape}")
    #             state = forecast(state)
    #             enkf_io.write_forecast(t + 1 if t < nt - 1 else t, state, ens_id)
                
    #     if t in tobserve:
    #         state = enkf_io.read_analysis(t, ens_id)
    #         print(f"Rank {subcomm.Get_rank()}, Analysis read at time {t}, state shape: {state.shape}")
    #         # state = analyze(state, H, d, ens_id, enkf_io, params)
    #         # enkf_io.write_analysis(t, state, ens_id)

    # enkf_io.close()

    # end_time = MPI.Wtime()
    # walltime = MPI.COMM_WORLD.reduce(end_time - start_time, op=MPI.MAX, root=0)
    # exec_time = MPI.COMM_WORLD.reduce(end_time - start_time, op=MPI.SUM, root=0)
    # if MPI.COMM_WORLD.Get_rank() == 0:
    #     print(f"Total execution time: {exec_time:.2f} seconds, Wall time: {walltime:.2f} seconds")
    #     print(f"Test completed successfully on {MPI.COMM_WORLD.Get_size()} ranks.")