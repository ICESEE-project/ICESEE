import numpy as np
import h5py
import os
from mpi4py import MPI
import gc
import uuid

# Mock forecast function to simulate model_module.forecast_step_single
def forecast_step_single(ensemble, **kwargs):
    # Simulate state update (e.g., increment by a small value)
    return {'state': ensemble + 0.1}

# Mock icesee_get_index function
def icesee_get_index(ensemble_vec, **kwargs):
    hdim = ensemble_vec.shape[0] // kwargs['num_vars']
    state_block_size = hdim * kwargs['num_vars']
    vecs = {'state': ensemble_vec[:state_block_size]}
    indx_map = {'state': slice(0, state_block_size)}
    dim_per_proc = state_block_size
    return vecs, indx_map, dim_per_proc

# Mock generate_enkf_field function
def generate_enkf_field(**kwargs):
    hdim = kwargs['hdim']
    return np.random.normal(0, 1, hdim)

# Create a sample HDF5 file for testing
def create_test_hdf5(filename, Nens, hdim, num_vars, k):
    data_shape = (hdim * num_vars, Nens, k)
    ensemble_data = np.random.rand(*data_shape)  # Random data for testing
    with h5py.File(filename, 'w') as f:
        f.create_dataset('ensemble', data=ensemble_data)

def main():
    # MPI setup
    comm_world = MPI.COMM_WORLD
    rank_world = comm_world.Get_rank()
    size_world = comm_world.Get_size()

    # Parameters
    params = {
        'total_state_param_vars': 2,  # Number of state variables
        'num_state_vars': 2,
        'sig_Q': [0.1, 0.1],  # Standard deviations for noise
    }
    Nens = 4  # Number of ensemble members
    rounds = 2
    subcomm_size = 2
    k = 0  # Time step index
    alpha = 0.5  # Noise mixing parameter
    dt = 0.1  # Time step
    Lx, Ly = 10.0, 10.0  # Domain size
    len_scale = 1.0  # Correlation length scale
    hdim = 10  # Horizontal dimension
    state_block_size = hdim * params['num_state_vars']

    # Create subcommunicator
    color = rank_world // subcomm_size
    subcomm = comm_world.Split(color, rank_world)
    sub_rank = subcomm.Get_rank()

    # Create test HDF5 file on rank 0
    input_file = f"test_ensemble_data_{uuid.uuid4()}.h5"
    if rank_world == 0:
        create_test_hdf5(input_file, Nens, hdim, params['total_state_param_vars'], k + 1)

    comm_world.Barrier()  # Ensure file is created before reading

    # Initialize results
    ens_list = []
    time_forecast_noise_generation = 0.0

    # Model kwargs
    model_kwargs = {
        'dt': dt,
        'obs_index': [0],
        'joint_estimation': False,
        'localization_flag': False,
        'num_vars': params['total_state_param_vars'],
    }

    # Main loop over rounds
    for round_id in range(rounds):
        ensemble_id = color + round_id * subcomm_size
        model_kwargs.update({'ens_id': ensemble_id, 'comm': subcomm})

        if ensemble_id < Nens:
            subcomm.Barrier()
            # Read from HDF5 file
            with h5py.File(input_file, 'r') as f:
                ensemble_vec = f['ensemble'][:, ensemble_id, k].copy()  # Copy to avoid read-only issues

            # Forecast step
            updated_state = forecast_step_single(ensemble=ensemble_vec, **model_kwargs)

            # Update ensemble vector
            vecs, indx_map, dim_per_proc = icesee_get_index(ensemble_vec, **model_kwargs)
            for key, value in updated_state.items():
                ensemble_vec[indx_map[key]] = value

            # Noise generation
            _time_forecast_noise_generation = MPI.Wtime()
            noise_all = []
            q0 = []
            for ii, sig in enumerate(params['sig_Q']):
                if ii < params['num_state_vars']:
                    model_kwargs.update({'ii_sig': ii, 'hdim': hdim})
                    W = generate_enkf_field(**model_kwargs)
                    noise_ = alpha * W  # Simplified: assume previous noise is zero
                    q0.append(noise_)
                    Z = np.sqrt(dt) * sig * noise_
                    noise_all.append(Z)
            noise_ = np.concatenate(noise_all, axis=0)
            ensemble_vec[:state_block_size] = ensemble_vec[:state_block_size] + noise_[:state_block_size]
            noise = np.concatenate(q0, axis=0)
            model_kwargs.update({'noise': noise})

            # Clean up memory
            del noise_all, q0, noise_, W
            time_forecast_noise_generation += MPI.Wtime() - _time_forecast_noise_generation

            # Gather results
            local_shape = ensemble_vec.size
            local_shapes = subcomm.gather(local_shape, root=0)

            if sub_rank == 0:
                total_size = sum(local_shapes)
                counts = local_shapes
                displs = [sum(counts[:i]) for i in range(subcomm.Get_size())]
                gathered_ensemble = np.empty(total_size, dtype=ensemble_vec.dtype)
            else:
                gathered_ensemble = None
                counts = None
                displs = None

            subcomm.Gatherv(ensemble_vec, [gathered_ensemble, counts, displs, MPI.DOUBLE], root=0)

            if sub_rank == 0:
                ens_list.append(gathered_ensemble)
            else:
                ens_list.append(None)

    # Clean up
    del ens_list
    gc.collect()

    # Finalize ensemble shape
    if rank_world == 0:
        gathered_ensemble_global = [ens for ens in ens_list if ens is not None]
        if gathered_ensemble_global:
            ensemble_vec = np.column_stack(gathered_ensemble_global)
            shape_ens = np.array(ensemble_vec.shape, dtype=np.int32)
        else:
            shape_ens = np.array([0, 0], dtype=np.int32)
    else:
        shape_ens = np.empty(2, dtype=np.int32)

    shape_ens = comm_world.bcast(shape_ens, root=0)

    # Clean up HDF5 file
    if rank_world == 0 and os.path.exists(input_file):
        os.remove(input_file)

    return shape_ens, time_forecast_noise_generation

if __name__ == '__main__':
    shape_ens, time_noise = main()
    if MPI.COMM_WORLD.Get_rank() == 0:
        print(f"Ensemble shape: {shape_ens}, Noise generation time: {time_noise:.4f} seconds")