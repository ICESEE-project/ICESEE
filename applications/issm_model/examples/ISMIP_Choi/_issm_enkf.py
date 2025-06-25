# ==============================================================================
# @des: This file contains run functions for ISSM data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2025-03-25
# @author: Brian Kyanjo
# ==============================================================================

import os
import numpy as np
import h5py
import gstools as gs

# --- import utility functions ---
from ICESEE.applications.issm_model.examples.ISMIP_Choi._issm_model import *
from ICESEE.config._utility_imports import icesee_get_index
from ICESEE.applications.issm_model.issm_utils.matlab2python.mat2py_utils import setup_ensemble_intial_data
from ICESEE.src.run_model_da.run_models_da import generate_enkf_field, generate_pseudo_random_field_1d


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

    # Do the true state run on the matlab side and only read the output on the python side once matlab is done with the simulation
    # --- call the issm model to generate the true state
    try:
        # -- call the run_model function to generate the true state
        kwargs.update({'k': 0})  # Set the initial time step
        ISSM_model(**kwargs)
    except Exception as e:
        print(f"[ICESEE Generate-True-State] Error generating true state: {e}")
        server.kill_matlab_processes()
        return None
    
    # On completion now fetch the true state from the Matlab output file to the ICESEE side (.h5 file)
    # -- fetch the true state vector
    statevec_true = kwargs.get('statevec_true')

    # -- call the icesee_get_index function to get the index of the state vector
    vecs, indx_map, dim_per_proc = icesee_get_index(statevec_true, **kwargs)

    # get the data extracted from the matlab output file
    input_filename = f'{icesee_path}/{data_path}/ensemble_true_state_{ens_id}.h5'
    with h5py.File(input_filename, 'r', driver='mpio', comm=comm) as f:
        # -- fetch state variables
        for k in range(1, kwargs.get('nt') + 1):
            key = f'Thickness_{k}'
            statevec_true[indx_map['Thickness'], k-1] = f[key][0]
            statevec_true[indx_map['bed'], k-1] = f['bed'][0]
            statevec_true[indx_map['coefficient'], k-1] = f['coefficient'][0]

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
    params = kwargs.get('params', {})

    nd = params.get('nd', 0)

    # --- fetch treu state vector
    statevec_nurged = kwargs.get('statevec_nurged')

    # -- call the icesee_get_index function to get the index of the state vector
    vecs, indx_map, dim_per_proc = icesee_get_index(statevec_nurged, **kwargs)

    # -- friction
    sill_friction = kwargs.get('sill_friction')
    range_friction = kwargs.get('range_friction')
    mean_friction  = kwargs.get('mean_friction')
    nugget_friction = kwargs.get('nugget_friction')
    friction_model = gs.Gaussian(dim=1, var=sill_friction, len_scale=range_friction, nugget=nugget_friction)
    friction_srf = gs.SRF(friction_model, mean=mean_friction, seed=kwargs.get('seed', 42))
    fdim = nd//params.get('total_state_param_vars', 1)
    x = np.linspace(0, range_friction*2, fdim)
    friction_field = friction_srf.structured([x])

    # --bed
    sill_bed = kwargs.get('sill_bed')
    range_bed = kwargs.get('range_bed')
    nugget_bed = kwargs.get('nugget_bed')
    bed_model = gs.Exponential(dim=1, var=sill_bed-nugget_bed, len_scale=range_bed, nugget=nugget_bed)
    bed_srf = gs.SRF(bed_model, mean=0, seed=kwargs.get('seed', 42))
    x = np.linspace(0, range_bed*2, fdim)
    bed_field = bed_srf.structured([x])

    # write the wrong states to a .h5 file to be read by the ISSM model before nurging
    friction_bed_filename = f'{icesee_path}/{data_path}/friction_bed_{ens_id}.h5'
    with h5py.File(friction_bed_filename, 'w', driver='mpio', comm=comm) as f:
        # -- write the friction field
        f.create_dataset('coefficient', data=friction_field)
        # -- write the bed field
        f.create_dataset('bed', data=bed_field)

    # -- call the run_model function to generate the nurged state
    try:
        kwargs.update({'k': 0})  # Set the initial time step
        ISSM_model(**kwargs)
    except Exception as e:
        print(f"[ICESEE Generate-Nurged-State] Error generating nurged state: {e}")
        server.kill_matlab_processes()

    # -- fetch the nurged state vector
    nurged_filename = f'{icesee_path}/{data_path}/ensemble_nurged_state_{ens_id}.h5'
    with h5py.File(nurged_filename, 'r', driver='mpio', comm=comm) as f:
        # -- fetch state variables
        for k in range(1, kwargs.get('nt') + 1):
            key = f'Thickness_{k}'
            statevec_nurged[indx_map['Thickness'], k-1] = f[key][0]
            statevec_nurged[indx_map['bed'], k-1] = f['bed'][0]
            statevec_nurged[indx_map['coefficient'], k-1] = f['coefficient'][0]

    #  --- change directory back to the original directory ---
    os.chdir(icesee_path)
    
    # return updated_state
    return statevec_nurged

        
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
        print(f"[ICESEE Initialize ensemble]] Error initializing ensemble: {e}")
        server.kill_matlab_processes()

    # if nprocs <= Nens then make fname available to all processes
    Nens = kwargs.get('Nens')
    size_world = kwargs.get('size_world', 1)
    if size_world <= Nens:
        data_dir = f'{issm_examples_dir}/Models/ens_id_0'
        setup_ensemble_intial_data(Nens, data_dir, fname)

    #  -- Read data from the ISSM side to be accessed by ICESEE on the python side
    output_filename = f'{icesee_path}/{data_path}/ensemble_out_{ens_id}.h5'
    updated_state = {}
    with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
        for key in vec_inputs:
            updated_state[key] = f[key][0]

    os.chdir(icesee_path)

    return updated_state
        