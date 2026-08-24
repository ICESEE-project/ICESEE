import numpy as np
import zarr
import numcodecs
import os
import gc

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
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    try:
        model = Model()
        n_model = 100
        if rank == 0:
            print("Generating H matrix and saving to Zarr...")
            H = model.H_matrix(n_model=n_model, zarr_path="output/H_matrix.zarr")

            # To read the Zarr file later
            # H_from_zarr = zarr.open("output/H_matrix.zarr", mode='r')[:]
            # print(f"Matrix equality check: {np.array_equal(H, H_from_zarr)}")  # Should print True
        comm.Barrier()

    except Exception as e:
        print(f"Error in main execution: {e}")