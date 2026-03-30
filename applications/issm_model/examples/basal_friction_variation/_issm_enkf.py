# ==============================================================================
# @des: This file contains run functions for ISSM data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2025-03-25
# @author: Brian Kyanjo
# ==============================================================================

import os, re
import numpy as np
import h5py
# import netCDF4
import gstools as gs
from mpi4py import MPI

# --- import utility functions ---
from _issm_model import *
from ICESEE.config._utility_imports import icesee_get_index
# from ICESEE.applications.issm_model.issm_utils.matlab2python.mat2py_utils import setup_ensemble_intial_data, MatlabServer

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

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

# --- inverse model after or before analysis for friction or velocity ---
def inverse_step_single(ensemble=None, **kwargs):
    """ensemble: packs the state variables and parameters of a single ensemble member
    Returns: ensemble: updated ensemble member
    """
    #  -- control time stepping   
    time = kwargs.get('t')
    km    = kwargs.get('km')
    # k   = kwargs.get('ik')
    # km   = kwargs.get('km')  # km is the time step for inverse model (can be before or after analysis)
    
    # kwargs.update({'tinitial': time[k], 'tfinal': time[k+1]})

    #  call the run_model fun to push the state forward in time
    return run_model_inverse(ensemble, **kwargs)

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
    vecs, indx_map, dim_per_proc = icesee_get_index(**kwargs)

    # get the data extracted from the matlab output file
    input_filename = f'{icesee_path}/{data_path}/ensemble_true_state_{ens_id}.h5'
    with h5py.File(input_filename, 'r') as f:
        # -- fetch state variables
        for k in range(1, kwargs.get('nt') + 1):
            for key in kwargs.get('vec_inputs'):
                key_name = f'{key}_{k}'
                statevec_true[indx_map[key], k-1] = f[key_name][0]
            
            # for key in kwargs.get('scalar_inputs', []):
            #     if key in f:
            #         statevec_true[indx_map[key], k-1] = f[key][:].reshape(-1, order='F')[0]

    updated_state = {}
    for key in vec_inputs:
        updated_state[key] = statevec_true[indx_map[key],:]
    
    # scalar_inputs = kwargs.get('scalar_inputs', [])
    # file_path = f'{icesee_path}/{data_path}/ensemble_true_state_scalar_{ens_id}.h5'
    # times, scalars = read_scalar_timeseries(file_path, scalar_inputs)

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
    vecs, indx_map, dim_per_proc = icesee_get_index(**kwargs)

    Lx = kwargs.get('Lx',params.get('Lx', 640e3))
    Ly = kwargs.get('Ly',params.get('Ly', 80e3))
    nx = kwargs.get('nx', params.get('nx', 32))
    ny = kwargs.get('ny', params.get('ny', 4))
    fdim = nd//params.get('total_state_param_vars', 1)
    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    seed_base = kwargs.get('seed_base', 42)

    #  # -- friction
    sill_friction = kwargs.get('sill_friction')
    range_friction = kwargs.get('range_friction')
    mean_friction  = kwargs.get('mean_friction')
    nugget_friction = kwargs.get('nugget_friction')
    # xx = np.linspace(0, range_friction, fdim)
    # # var_fric = max(sill_friction - nugget_friction, 0.0)
    # friction_model = gs.Gaussian(dim=1, var=sill_friction, len_scale=range_friction, nugget=nugget_friction)
    # friction_srf = gs.SRF(friction_model, seed=42)
    # # friction_field = np.asarray(friction_srf.structured([x, y])).reshape(-1)[:fdim]
    # friction_field = np.asarray(friction_srf.structured([xx])).reshape(-1)

    file_path = f'{icesee_path}/{data_path}/mesh_idxy_{0}.h5'
    with h5py.File(file_path, 'r') as f:
        x_param = f['/fric_x'][:]   # shape (fdim,)
        y_param = f['/fric_y'][:]   # shape (fdim,)

    # scale coords by correlation length so len_scale ~ 1
    x_scaled = x_param / range_friction
    y_scaled = y_param / range_friction

    model = gs.Gaussian(
        dim=2,
        var=sill_friction,
        len_scale=range_friction,
        nugget=nugget_friction,
    )

    srf = gs.SRF(model, seed=seed_base)

    # unstructured evaluation at real node positions
    # friction_field = np.asarray(srf((x_scaled, y_scaled)))  # (fdim,)
    friction_field = np.asarray(srf((x_param, y_param)))  # (fdim,)

    # # --bed
    # sill_bed = kwargs.get('sill_bed')
    # range_bed = kwargs.get('range_bed')
    # nugget_bed = kwargs.get('nugget_bed')

    # bed_kriging_file = f'{icesee_path}/bed_kriging_results.h5'
    # with h5py.File(bed_kriging_file, 'r') as f:
    #     bed_field = f['bed_ens'][...]

    # # bed_field = np.mean(bed_field, axis=0)
    # bed_field = bed_field[0, :]

    # write the wrong states to a .h5 file to be read by the ISSM model before nurging
    friction_bed_filename = f'{icesee_path}/{data_path}/friction_bed_{ens_id}.h5'
    # with h5py.File(friction_bed_filename, 'w', driver='mpio', comm=comm) as f:
    with h5py.File(friction_bed_filename, 'w') as f:
        # -- write the friction field
        f.create_dataset('FrictionCoefficient', data=friction_field)
        # -- write the bed field
        # f.create_dataset('Bed', data=bed_field)

    # -- call the run_model function to generate the nurged state
    try:
        kwargs.update({'k': 0})  # Set the initial time step
        ISSM_model(**kwargs)
    except Exception as e:
        print(f"[ICESEE Generate-Nurged-State] Error generating nurged state: {e}")
        server.kill_matlab_processes()

    # -- fetch the nurged state vector
    nurged_filename = f'{icesee_path}/{data_path}/ensemble_nurged_state_{ens_id}.h5'
    # with h5py.File(nurged_filename, 'r', driver='mpio', comm=comm) as f:
    with h5py.File(nurged_filename, 'r') as f:
        # -- fetch state variables
        for k in range(1, kwargs.get('nt') + 1):
            for key in kwargs.get('vec_inputs'):
                key_name = f'{key}_{k}'
                statevec_nurged[indx_map[key], k-1] = f[key_name][0]

    updated_state = {}
    for key in vec_inputs:
        updated_state[key] = statevec_nurged[indx_map[key], :]    

    #  --- change directory back to the original directory ---
    os.chdir(icesee_path)
    
    return updated_state

        
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
    params              = kwargs.get('params', {})
    nd                  = params.get('nd', 0)

    #  --- change directory to the issm directory ---
    os.chdir(issm_examples_dir)
    # ens_id = kwargs.get('ens_id')
    ens_id =  ens
    kwargs.update({'ens_id': ens_id})

    #  -- control time stepping
    kwargs.update({'k':0}) 
    dt = kwargs.get('dt')
    kwargs.update({'tinitial': 0, 'tfinal': dt})


    # --- filename for data saving
    fname = 'initialize_ensemble.mat'
    kwargs.update({'fname': fname})

    #*-----------------------
    Lx = kwargs.get('Lx',params.get('Lx', 640e3))
    Ly = kwargs.get('Ly',params.get('Ly', 80e3))
    fdim = nd//params.get('total_state_param_vars', 1)
    nx = kwargs.get('nx', params.get('nx', 32))
    ny = kwargs.get('ny', params.get('ny', 4))
    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    seed_base = kwargs.get('seed_base', 42)

    #  # -- friction
    sill_friction = kwargs.get('sill_friction')
    range_friction = kwargs.get('range_friction')
    mean_friction  = kwargs.get('mean_friction')
    nugget_friction = kwargs.get('nugget_friction')

    file_path = f'{icesee_path}/{data_path}/mesh_idxy_{0}.h5'
    with h5py.File(file_path, 'r') as f:
        x_param = f['/fric_x'][:]   # shape (fdim,)
        y_param = f['/fric_y'][:]   # shape (fdim,)

    # scale coords by correlation length so len_scale ~ 1
    x_scaled = x_param / range_friction
    y_scaled = y_param / range_friction

    model = gs.Gaussian(
        dim=2,
        var=sill_friction,
        len_scale=range_friction,
        nugget=nugget_friction,
    )

    srf = gs.SRF(model, seed=seed_base+ ens)

    # unstructured evaluation at real node positions
    friction_field = np.asarray(srf((x_param, y_param)))  # (fdim,)


    # bed_kriging_file = f'{icesee_path}/bed_kriging_results.h5'
    # with h5py.File(bed_kriging_file, 'r') as f:
    #     bed_field = f['bed_ens'][ens, :]

    # write the wrong states to a .h5 file to be read by the ISSM model before nurging
    friction_bed_filename = f'{icesee_path}/{data_path}/friction_bed_{ens_id}.h5'
    # with h5py.File(friction_bed_filename, 'w', driver='mpio', comm=comm) as f:
    with h5py.File(friction_bed_filename, 'w') as f:
        # -- write the friction field
        f.create_dataset('FrictionCoefficient', data=friction_field)
        # -- write the bed field
        # f.create_dataset('Bed', data=bed_field)
    #*-----------------------

    enkf_scalar_file = f'{icesee_path}/{data_path}/ensemble_out_scalar_{ens_id}.h5'
    if os.path.exists(enkf_scalar_file):
        os.remove(enkf_scalar_file)
    with h5py.File(enkf_scalar_file, 'w') as f:
        nt = kwargs.get('nt')
        for key in kwargs.get('scalar_inputs', []):
            f.create_dataset(key, shape=(nt,), dtype=np.float64)
            f[key][:] = np.nan

    try:
        # -- call the run_model function to initialize the ensemble members
        ISSM_model(**kwargs)
    except Exception as e:
        print(f"[ICESEE Initialize ensemble]] Error initializing ensemble: {e}")
        server.kill_matlab_processes()

    
    #  -- Read data from the ISSM side to be accessed by ICESEE on the python side
    output_filename = f'{icesee_path}/{data_path}/ensemble_out_{ens_id}.h5'
    updated_state = {}
    # with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
    with h5py.File(output_filename, 'r') as f:
        for key in kwargs.get('vec_inputs'):
            if key in kwargs.get('joint_estimated_params') and not kwargs.get('joint_estimation'):
                print(f"[ICESEE run_model Warning] Skipping joint estimated parameter '{key}' in output file since joint_estimation is set to False.\n"); 
                continue # skip the joint estimated parameters if we are not doing joint estimation
            if key in f:
                updated_state[key] = f[key][:].reshape(-1, order='F')
            else:
                print(f"[ICESEE initialize ensemble Warning] Key '{key}' not found in output file: {output_filename}")

    # --- compute the mean of the scalar inputs across all ensemble members and save in a separate file
    # first check if all Nens scalar output files are available before trying to read and compute the mean
    expected_scalar_files = [f'{icesee_path}/{data_path}/ensemble_out_scalar_{i}.h5' for i in range(kwargs.get('Nens', 1))]
    actual_scalar_files = [f for f in expected_scalar_files if os.path.exists(f)]
    if len(actual_scalar_files) == kwargs.get('Nens', 1):
        comm = MPI.COMM_WORLD
        comm.Barrier()
        k = 0  # initial time step index
        if comm.Get_rank() == 0:
            nt = kwargs.get('nt')
            nens = kwargs.get('Nens', 1)
            scalar_inputs = kwargs.get('scalar_inputs', [])
            scalar_means_file = f'{icesee_path}/{data_path}/ensemble_scalar_output.h5'

            if os.path.exists(scalar_means_file):
                os.remove(scalar_means_file)

            with h5py.File(scalar_means_file, 'w') as f_out:
                for key in kwargs.get('scalar_inputs', []):
                    f_out.create_dataset(key, shape=(nt,), dtype=np.float64)
                    
                    vals = []
                    for i in range(nens):
                        enkf_scalar_file = f'{icesee_path}/{data_path}/ensemble_out_scalar_{i}.h5'

                        with h5py.File(enkf_scalar_file, 'r') as f_ens:
                            if key in f_ens:
                                arr = np.asarray(f_ens[key][:]).reshape(-1, order='F')

                                if arr.size > k:
                                    vals.append(arr[k])
                                else:
                                    vals.append(np.nan)
                            else:
                                print(f"[ICESEE Warning] Key '{key}' not found in scalar output file: {enkf_scalar_file}")
                                vals.append(np.nan)

                    f_out[key][k] = np.nanmean(vals)

        comm.Barrier()
       
    os.chdir(icesee_path)
    
    return updated_state

