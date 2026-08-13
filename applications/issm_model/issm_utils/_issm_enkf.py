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
from ICESEE.applications.issm_model.examples.ISMIP._issm_model import *
from ICESEE.config._utility_imports import icesee_get_index
from ICESEE.applications.issm_model.issm_utils._coordinates import (
    register_issm_coordinate_provider,
)

register_issm_coordinate_provider()

# --- Forecast step ---
def forecast_step_single(ensemble=None, **icesee_kwargs):
    """ensemble: packs the state variables and parameters of a single ensemble member
    Returns: ensemble: updated ensemble member
    """
    #  -- control time stepping
    time = icesee_kwargs.get('t')
    k    = icesee_kwargs.get('k')

    icesee_kwargs.update({'tinitial': time[k], 'tfinal': time[k+1]})

    #  call the run_model fun to push the state forward in time
    return run_model(ensemble, **icesee_kwargs)


# --- generate true state ---
def generate_true_state(**icesee_kwargs):
    """des: generate the true state of the model
    Returns: true_state: the true state of the model
    """
    time   = icesee_kwargs.get('t')
    server = icesee_kwargs.get('server')

    issm_examples_dir   = icesee_kwargs.get('issm_examples_dir')
    icesee_path         = icesee_kwargs.get('icesee_path')
    data_path           = icesee_kwargs.get('data_path')
    comm                = icesee_kwargs.get('comm')
    vec_inputs          = icesee_kwargs.get('vec_inputs')

    #  --- change directory to the issm directory ---
    os.chdir(issm_examples_dir)

    # --- filename for data saving
    fname = 'true_state.mat'
    icesee_kwargs.update({'fname': fname})
    ens_id = icesee_kwargs.get('ens_id')

    # try:
    if True:
        # --- fetch treu state vector
        statevec_true = icesee_kwargs.get('statevec_true')

        # -- call the icesee_get_index function to get the index of the state vector
        vecs, indx_map, dim_per_proc = icesee_get_index(statevec_true, **icesee_kwargs)

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
        for k in range(icesee_kwargs.get('nt')):
            icesee_kwargs.update({'k': k})
            time = icesee_kwargs.get('t')
            icesee_kwargs.update({'tinitial': time[k], 'tfinal': time[k+1]})
            # --- write the state back to h5 file for ISSM model
            input_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
            with h5py.File(input_filename, 'w', driver='mpio', comm=comm) as f:
                for key in vec_inputs:
                    f.create_dataset(key, data=statevec_true[indx_map[key],k])

            # -- call the run_model function to push the state forward in time
            ISSM_model(**icesee_kwargs)

            # try:
            if True:
                output_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
                with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
                    for key in vec_inputs:
                        statevec_true[indx_map[key],k+1] = f[key][0]

            # except Exception as e:
            #     print(f"[Generate-True-State: read output file] Error reading the file: {e}")
                # return None

        updated_state = {}
        for key in vec_inputs:
            updated_state[key] = statevec_true[indx_map[key],:]

        #  --- change directory back to the original directory ---
        os.chdir(icesee_path)

        return updated_state


def generate_nurged_state(**icesee_kwargs):
    """generate the nurged state of the model"""
    time   = icesee_kwargs.get('t')
    server = icesee_kwargs.get('server')
    issm_examples_dir   = icesee_kwargs.get('issm_examples_dir')
    icesee_path         = icesee_kwargs.get('icesee_path')
    data_path           = icesee_kwargs.get('data_path')
    comm                = icesee_kwargs.get('comm')
    vec_inputs          = icesee_kwargs.get('vec_inputs')

    #  --- change directory to the issm directory ---
    os.chdir(issm_examples_dir)

    # get the rank of the current process
    rank = comm.Get_rank()

    # --- filename for data saving
    fname = 'nurged_state.mat'
    icesee_kwargs.update({'fname': fname})
    ens_id = icesee_kwargs.get('ens_id')

    try:
    # if True:
        # --- fetch treu state vector
        statevec_nurged = icesee_kwargs.get('statevec_nurged')

        # -- call the icesee_get_index function to get the index of the state vector
        vecs, indx_map, dim_per_proc = icesee_get_index(statevec_nurged, **icesee_kwargs)

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
        for k in range(icesee_kwargs.get('nt')):
            icesee_kwargs.update({'k': k})
            time = icesee_kwargs.get('t')
            icesee_kwargs.update({'tinitial': time[k], 'tfinal': time[k+1]})

            # --- write the state back to h5 file for ISSM model
            input_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
            with h5py.File(input_filename, 'w', driver='mpio', comm=comm) as f:
                for key in vec_inputs:
                    f.create_dataset(key, data=statevec_nurged[indx_map[key],k])

            # -- call the run_model function to push the state forward in time
            ISSM_model(**icesee_kwargs)

            try:
                output_filename = f'{icesee_path}/{data_path}/ensemble_output_{ens_id}.h5'
                with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
                    for key in vec_inputs:
                        statevec_nurged[indx_map[key],k+1] = f[key][0]

            except Exception as e:
                print(f"[DEBUG] Error reading the file: {e}")
                # return None

        updated_state = {'Vx': statevec_nurged[indx_map["Vx"],:],
                        'Vy': statevec_nurged[indx_map["Vy"],:],
                        'Vz': statevec_nurged[indx_map["Vz"],:],
                        'Pressure': statevec_nurged[indx_map["Pressure"],:]}

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
def initialize_ensemble(ens, **icesee_kwargs):
    """des: initialize the ensemble members
    Returns: ensemble: the ensemble members
    """
    import h5py
    import os, sys

    server              = icesee_kwargs.get('server')
    issm_examples_dir   = icesee_kwargs.get('issm_examples_dir')
    icesee_path         = icesee_kwargs.get('icesee_path')
    data_path           = icesee_kwargs.get('data_path')
    comm                = icesee_kwargs.get('comm')
    vec_inputs          = icesee_kwargs.get('vec_inputs')

    #  --- change directory to the issm directory ---
    os.chdir(issm_examples_dir)
    ens_id = icesee_kwargs.get('ens_id')

    # get the rank of the current process
    rank = comm.Get_rank()

    # try:
    if True:
        output_filename = f'{icesee_path}/{data_path}/ensemble_init_{ens_id}.h5'
        updated_state = {}
        with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
            for key in vec_inputs:
                updated_state[key] = f[key][0]

        os.chdir(icesee_path)

        return updated_state
