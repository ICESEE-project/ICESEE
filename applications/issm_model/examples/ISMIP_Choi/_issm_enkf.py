# ==============================================================================
# @des: This file contains run functions for ISSM data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2025-03-25
# @author: Brian Kyanjo
# ==============================================================================

import os
import numpy as np
import h5py


# --- import utility functions ---
from ICESEE.applications.issm_model.examples.ISMIP_Choi._issm_model import *
from ICESEE.config._utility_imports import icesee_get_index
from ICESEE.applications.issm_model.issm_utils.matlab2python.mat2py_utils import setup_ensemble_intial_data

# --- Forecast step ---
def forecast_step_single(ensemble=None, **kwargs):
    """ensemble: packs the state variables and parameters of a single ensemble member
    Returns: ensemble: updated ensemble member
    """
    #  -- control time stepping   
    time = kwargs.get('t')
    k    = kwargs.get('k')
    
    kwargs.update({'tinitial': time[k], 'tfinal': time[k+1]})

    #  call the run_model fun to push the state forward in time
    return run_model(ensemble, **kwargs)


# --- generate true state ---
def generate_true_state(**kwargs):
    """des: generate the true state of the model
    Returns: true_state: the true state of the model
    """
    params = kwargs.get('params')
    time   = kwargs.get('t')
    server = kwargs.get('server')
    
    issm_examples_dir   = kwargs.get('issm_examples_dir')
    icesee_path         = kwargs.get('icesee_path')
    data_path           = kwargs.get('data_path')
    comm                = kwargs.get('comm')
    vec_inputs          = kwargs.get('vec_inputs')

    #  --- change directory to the issm directory ---
    os.chdir(issm_examples_dir)

    # --- filename for data saving
    fname = 'true_state.mat'
    kwargs.update({'fname': fname})
    ens_id = kwargs.get('ens_id')

    # try:
    if True:
        # --- fetch treu state vector
        statevec_true = kwargs.get('statevec_true')

        # -- call the icesee_get_index function to get the index of the state vector
        vecs, indx_map, dim_per_proc = icesee_get_index(statevec_true, **kwargs)

        # -- fetch data from inital state
        try: 
        # if True:
            output_filename = f'{icesee_path}/{data_path}/ensemble_init_{ens_id}.h5'
            # print(f"[DEBUG-0] Attempting to open file: {output_filename}")
            if not os.path.exists(output_filename):
                print(f"[ERROR] File does not exist: {output_filename}")
                return None
            with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
                # -- fetch state variables
                for key in vec_inputs:
                    statevec_true[indx_map[key],0] = f[key][0]
        except Exception as e:
            print(f"[Generate-True-State: read init file] Error reading the file: {e}")
                            
        # -- mimic the time integration loop to save vec on every time step
        for k in range(kwargs.get('nt')):
            kwargs.update({'k': k})
            time = kwargs.get('t')
            kwargs.update({'tinitial': time[k], 'tfinal': time[k+1]})
            # --- write the state back to h5 file for ISSM model
            input_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
            with h5py.File(input_filename, 'w', driver='mpio', comm=comm) as f:
                for key in vec_inputs:
                    f.create_dataset(key, data=statevec_true[indx_map[key],k])

            # -- call the run_model function to push the state forward in time
            ISSM_model(**kwargs)
           
            # try:
            if True:
                output_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
                with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
                    for key in vec_inputs:
                        statevec_true[indx_map[key],k+1] = f[key][0]
                    
            # except Exception as e:
            #     print(f"[Generate-True-State: read output file] Error reading the file: {e}")
                # return None

        # Precompute filename
        # output_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'

        # # Open HDF5 file once in read/write mode
        # with h5py.File(output_filename, 'a', driver='mpio', comm=comm) as f:
        #     # Time integration loop
        #     for k in range(kwargs.get('nt')):
        #         # Update kwargs with current time step parameters
        #         kwargs.update({'k': k, 'tinitial': time[k], 'tfinal': time[k + 1]})

        #         print(f"[Generate-True-State:] time step {k} , tinitial: {kwargs.get('tinitial')}, tfinal: {kwargs.get('tfinal')}")
                
        #         # Write current state to HDF5 file
        #         try:
        #             for key in vec_inputs:
        #                 # Create or overwrite dataset for the current state
        #                 if key in f:
        #                     del f[key]  # Delete existing dataset if it exists
        #                 f.create_dataset(key, data=statevec_true[indx_map[key], k])
        #         except Exception as e:
        #             print(f"[Generate-True-State:] Error writing to file at step {k}: {e}")
        #             raise
                
        #         # Run ISSM model to advance state
        #         ISSM_model(**kwargs)
                
        #         # Read updated state
        #         try:
        #             for key in vec_inputs:
        #                 statevec_true[indx_map[key], k + 1] = f[key][0]
        #         except Exception as e:
        #             print(f"[Generate-True-State:] Error reading from file at step {k}: {e}")
        #             raise
    
        updated_state = {}
        for key in vec_inputs:
            updated_state[key] = statevec_true[indx_map[key],:]

        #  --- change directory back to the original directory ---
        os.chdir(icesee_path)
        
        return updated_state


def generate_nurged_state(**kwargs):
    """generate the nurged state of the model"""
    params = kwargs.get('params')
    time   = kwargs.get('t')
    server = kwargs.get('server')
    issm_examples_dir   = kwargs.get('issm_examples_dir')
    icesee_path         = kwargs.get('icesee_path')
    data_path           = kwargs.get('data_path')
    comm                = kwargs.get('comm')
    vec_inputs          = kwargs.get('vec_inputs')       

    #  --- change directory to the issm directory ---
    os.chdir(issm_examples_dir)

    # get the rank of the current process
    rank = comm.Get_rank()

    # --- filename for data saving
    fname = 'nurged_state.mat'
    kwargs.update({'fname': fname})
    ens_id = kwargs.get('ens_id')

    try:
    # if True:
        # --- fetch treu state vector
        statevec_nurged = kwargs.get('statevec_nurged')

        # -- call the icesee_get_index function to get the index of the state vector
        vecs, indx_map, dim_per_proc = icesee_get_index(statevec_nurged, **kwargs)

        # --- create a bump -1== to 0

        # -- fetch data from inital state
        try: 
            output_filename = f'{icesee_path}/{data_path}/ensemble_init_{ens_id}.h5'
            # print(f"[DEBUG] Attempting to open file: {output_filename}")
            if not os.path.exists(output_filename):
                print(f"[ERROR] File does not exist: {output_filename}")
                return None
            with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
                # -- fetch state variables
                for key in vec_inputs:
                    statevec_nurged[indx_map[key],0] = f[key][0]

        except Exception as e:
            print(f"[DEBUG] Error reading the file: {e}")
                            
        # -- mimic the time integration loop to save vec on every time step
        for k in range(kwargs.get('nt')):
            kwargs.update({'k': k})
            time = kwargs.get('t')
            kwargs.update({'tinitial': time[k], 'tfinal': time[k+1]})

            # --- write the state back to h5 file for ISSM model
            input_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
            with h5py.File(input_filename, 'w', driver='mpio', comm=comm) as f:
                for key in vec_inputs:
                    f.create_dataset(key, data=statevec_nurged[indx_map[key],k])

            # -- call the run_model function to push the state forward in time
            ISSM_model(**kwargs)

            try:
                output_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
                with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
                    for key in vec_inputs:
                        statevec_nurged[indx_map[key],k+1] = f[key][0]

            except Exception as e:
                print(f"[DEBUG] Error reading the file: {e}")
                # return None

        # Precompute filename
        # output_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'

        # # Open HDF5 file once in read/write mode
        # with h5py.File(output_filename, 'a', driver='mpio', comm=comm) as f:
        #     # Time integration loop
        #     for k in range(kwargs.get('nt')):
        #         # Update kwargs with current time step parameters
        #         kwargs.update({'k': k, 'tinitial': time[k], 'tfinal': time[k + 1]})
                
        #         # Write current state to HDF5 file
        #         try:
        #             for key in vec_inputs:
        #                 # Create or overwrite dataset for the current state
        #                 if key in f:
        #                     del f[key]  # Delete existing dataset if it exists
        #                 f.create_dataset(key, data=statevec_nurged[indx_map[key], k])
        #         except Exception as e:
        #             print(f"[nurged state:] Error writing to file at step {k}: {e}")
        #             raise
                
        #         # Run ISSM model to advance state
        #         ISSM_model(**kwargs)
                
        #         # Read updated state
        #         try:
        #             for key in vec_inputs:
        #                 statevec_nurged[indx_map[key], k + 1] = f[key][0]
        #         except Exception as e:
        #             print(f"[nurged state:] Error reading from file at step {k}: {e}")
        #             raise

        #  --- change directory back to the original directory ---
        os.chdir(icesee_path)
        
        # return updated_state
        return statevec_nurged
    
    except Exception as e:
        print(f"[DEBUG] Error sending command: {e}")
        # Ensure directory is changed back even on error
        os.chdir(icesee_path)
        return None
        
#  --- initialize ensemble members ---
def initialize_ensemble(ens, **kwargs):
    """des: initialize the ensemble members
    Returns: ensemble: the ensemble members
    """
    import h5py
    import os, sys

    server              = kwargs.get('server')
    issm_examples_dir   = kwargs.get('issm_examples_dir')
    icesee_path         = kwargs.get('icesee_path')
    data_path           = kwargs.get('data_path')
    comm                = kwargs.get('comm')
    vec_inputs          = kwargs.get('vec_inputs')

    #  --- change directory to the issm directory ---
    os.chdir(issm_examples_dir)
    ens_id = kwargs.get('ens_id')

    #  -- control time stepping
    kwargs.update({'k':0}) 
    kwargs.update({'tinitial': 0, 'tfinal': 0.2})


    # --- filename for data saving
    fname = 'initialize_ensemble.mat'
    kwargs.update({'fname': fname})

    try:
        # -- call the run_model function to initialize the ensemble members
        ISSM_model(**kwargs)
    except Exception as e:
        print(f"[Initialize ensemble]] Error initializing ensemble: {e}")
        server.kill_matlab_processes()

    # if nprocs <= Nens then make fname available to all processes
    Nens = kwargs.get('Nens')
    size_world = kwargs.get('size_world', 1)
    if size_world <= Nens:
        data_dir = f'{issm_examples_dir}/Models/ens_id_0'
        setup_ensemble_intial_data(Nens, data_dir, fname)

   
    try:
        #  -- Read data from the ISSM side to be accessed by ICESEE on the python side
        output_filename = f'{icesee_path}/{data_path}/ensemble_out_{ens_id}.h5'
        updated_state = {}
        with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
            for key in vec_inputs:
                updated_state[key] = f[key][0]
    except Exception as e:
        print(f"[Initialize ensemble] Error reading the file: {e}")
        server.kill_matlab_processes()

    os.chdir(icesee_path)

    return updated_state
        