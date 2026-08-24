# ==============================================================================
# @des: This file contains run functions for ISSM data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2025-03-25
# @author: Brian Kyanjo
# ==============================================================================

import os
import numpy as np
import h5py
# import netCDF4
import gstools as gs

# --- import utility functions ---
from ICESEE.applications.issm_model.examples.ISMIP_Choi._issm_model import *
from ICESEE.config._utility_imports import icesee_get_index
from ICESEE.applications.issm_model.issm_utils._coordinates import (
    register_issm_coordinate_provider,
)
# from ICESEE.applications.issm_model.issm_utils.matlab2python.mat2py_utils import setup_ensemble_intial_data, MatlabServer


register_issm_coordinate_provider()

# --- import kriging functions ---
def _load_bed_from_kriging_file(bed_kriging_file, fdim, ens=None, use_mean=False):
    with h5py.File(bed_kriging_file, "r") as f:
        bed_ens = f["bed_ens"][:]

        # expected Python shape: (npts, Ne), but be safe
        if bed_ens.shape[0] != fdim and bed_ens.shape[1] == fdim:
            bed_ens = bed_ens.T

        if bed_ens.shape[0] != fdim:
            raise ValueError(
                f"bed_ens node mismatch: got {bed_ens.shape}, expected first dim {fdim}"
            )

        if use_mean:
            bed_field = np.mean(bed_ens, axis=1)
        else:
            if ens is None:
                ens = 0
            ens = int(ens) % bed_ens.shape[1]
            bed_field = bed_ens[:, ens]

    return np.asarray(bed_field, dtype=float).ravel()

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

# --- inverse model after or before analysis for friction or velocity ---
def inverse_step_single(ensemble=None, **icesee_kwargs):
    """ensemble: packs the state variables and parameters of a single ensemble member
    Returns: ensemble: updated ensemble member
    """
    #  -- control time stepping
    time = icesee_kwargs.get('t')
    km    = icesee_kwargs.get('km')
    # k   = icesee_kwargs.get('ik')
    # km   = icesee_kwargs.get('km')  # km is the time step for inverse model (can be before or after analysis)

    # icesee_kwargs.update({'tinitial': time[k], 'tfinal': time[k+1]})

    #  call the run_model fun to push the state forward in time
    return run_model_inverse(ensemble, **icesee_kwargs)

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
    initial_state_only = bool(icesee_kwargs.get('initial_state_only', False))
    fname = 'initial_true_state.mat' if initial_state_only else 'true_state.mat'
    icesee_kwargs.update({'fname': fname})
    ens_id = icesee_kwargs.get('ens_id')

    # Do the true state run on the matlab side and only read the output on the python side once matlab is done with the simulation
    # --- call the issm model to generate the true state
    try:
        # -- call the run_model function to generate the true state
        icesee_kwargs.update({'k': 0})  # Set the initial time step
        ISSM_model(**icesee_kwargs)
    except Exception as e:
        print(f"[ICESEE Generate-True-State] Error generating true state: {e}")
        server.kill_matlab_processes()
        return None

    # On completion now fetch the true state from the Matlab output file to the ICESEE side (.h5 file)
    # -- fetch the true state vector
    statevec_true = icesee_kwargs.get('statevec_true')

    # -- call the icesee_get_index function to get the index of the state vector
    vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)

    # get the data extracted from the matlab output file
    input_filename = f'{icesee_path}/{data_path}/ensemble_true_state_{ens_id}.h5'
    with h5py.File(input_filename, 'r') as f:
        # -- fetch state variables
        output_count = 1 if initial_state_only else icesee_kwargs.get('nt')
        for k in range(1, output_count + 1):
            key_Thickness=f'Thickness_{k}'
            # key_base = f'Base_{k}'
            key_surface = f'Surface_{k}'
            key_u  = f'Vx_{k}'
            key_v  = f'Vy_{k}'
            key_bed = f'bed_{k}'
            key_coefficient = f'coefficient_{k}'
            statevec_true[indx_map['Thickness'], k-1] = f[key_Thickness][0]
            # statevec_true[indx_map['Base'], k-1] = f[key_base][0]
            statevec_true[indx_map['Surface'], k-1] = f[key_surface][0]
            statevec_true[indx_map['Vx'], k-1] = f[key_u][0]
            statevec_true[indx_map['Vy'], k-1] = f[key_v][0]
            statevec_true[indx_map['bed'], k-1] = f[key_bed][0]
            statevec_true[indx_map['coefficient'], k-1] = f[key_coefficient][0]

    # In fully-parallel mode ``statevec_true`` is the on-disk HDF5 target.
    # Re-slicing every variable here would materialize the complete trajectory
    # in memory after it has already been written successfully.
    if isinstance(statevec_true, h5py.Dataset):
        statevec_true.flush()
        os.chdir(icesee_path)
        return None

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
    initial_state_only = bool(icesee_kwargs.get('initial_state_only', False))
    fname = 'initial_nurged_state.mat' if initial_state_only else 'nurged_state.mat'
    icesee_kwargs.update({'fname': fname})
    ens_id = icesee_kwargs.get('ens_id')

    nd = icesee_kwargs.get('nd', 0)

    # --- fetch treu state vector
    statevec_nurged = icesee_kwargs.get('statevec_nurged')

    # -- call the icesee_get_index function to get the index of the state vector
    vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)

    Lx = icesee_kwargs.get('Lx', 640e3)
    Ly = icesee_kwargs.get('Ly', 80e3)
    nx = icesee_kwargs.get('nx', 32)
    ny = icesee_kwargs.get('ny', 4)
    fdim = nd//icesee_kwargs.get('total_state_param_vars', 1)
    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    seed_base = icesee_kwargs.get('seed_base', 42)

    #  # -- friction
    sill_friction = icesee_kwargs.get('sill_friction')
    range_friction = icesee_kwargs.get('range_friction')
    mean_friction  = icesee_kwargs.get('mean_friction')
    nugget_friction = icesee_kwargs.get('nugget_friction')
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

    # --bed
    sill_bed = icesee_kwargs.get('sill_bed')
    range_bed = icesee_kwargs.get('range_bed')
    nugget_bed = icesee_kwargs.get('nugget_bed')

    bed_kriging_file = f'{icesee_path}/bed_kriging_results.h5'
    bed_field = _load_bed_from_kriging_file(
        bed_kriging_file,
        fdim,
        use_mean=True
    )

    # write the wrong states to a .h5 file to be read by the ISSM model before nurging
    friction_bed_filename = f'{icesee_path}/{data_path}/friction_bed_{ens_id}.h5'
    # with h5py.File(friction_bed_filename, 'w', driver='mpio', comm=comm) as f:
    with h5py.File(friction_bed_filename, 'w') as f:
        # -- write the friction field
        f.create_dataset('coefficient', data=friction_field)
        # -- write the bed field
        f.create_dataset('bed', data=bed_field)

    # -- call the run_model function to generate the nurged state
    try:
        icesee_kwargs.update({'k': 0})  # Set the initial time step
        ISSM_model(**icesee_kwargs)
    except Exception as e:
        print(f"[ICESEE Generate-Nurged-State] Error generating nurged state: {e}")
        server.kill_matlab_processes()

    # -- fetch the nurged state vector
    nurged_filename = f'{icesee_path}/{data_path}/ensemble_nurged_state_{ens_id}.h5'
    # with h5py.File(nurged_filename, 'r', driver='mpio', comm=comm) as f:
    with h5py.File(nurged_filename, 'r') as f:
        # -- fetch state variables
        output_count = 1 if initial_state_only else icesee_kwargs.get('nt')
        for k in range(1, output_count + 1):
            # key_thickness=f'Thickness_{k}'
            key_Thickness=f'Thickness_{k}'
            # key_base = f'Base_{k}'
            key_surface = f'Surface_{k}'
            key_u = f'Vx_{k}'
            key_v = f'Vy_{k}'
            key_bed = f'bed_{k}'
            key_coefficient = f'coefficient_{k}'
            statevec_nurged[indx_map['Thickness'], k-1] = f[key_Thickness][0]
            # statevec_nurged[indx_map['Base'], k-1] = f[key_base][0]
            statevec_nurged[indx_map['Surface'], k-1] = f[key_surface][0]
            statevec_nurged[indx_map['Vx'], k-1] = f[key_u][0]
            statevec_nurged[indx_map['Vy'], k-1] = f[key_v][0]
            statevec_nurged[indx_map['bed'], k-1] = f[key_bed][0]
            statevec_nurged[indx_map['coefficient'], k-1] = f[key_coefficient][0]
            # statevec_nurged[indx_map['bed'], k-1] = f['bed'][0]
            # statevec_nurged[indx_map['coefficient'], k-1] = f['coefficient'][0]

    # See ``generate_true_state``: signal successful in-place completion
    # without reconstructing a second copy of the full trajectory in RAM.
    if isinstance(statevec_nurged, h5py.Dataset):
        statevec_nurged.flush()
        os.chdir(icesee_path)
        return None

    updated_state = {}
    for key in vec_inputs:
        updated_state[key] = statevec_nurged[indx_map[key], :]

    #  --- change directory back to the original directory ---
    os.chdir(icesee_path)

    return updated_state


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
    nd                  = icesee_kwargs.get('nd', 0)

    #  --- change directory to the issm directory ---
    os.chdir(issm_examples_dir)
    # ens_id = icesee_kwargs.get('ens_id')
    ens_id =  ens
    icesee_kwargs.update({'ens_id': ens_id})

    # -- Optional pre-assimilation dynamic spin-up.  Historical profiles
    # retain the former one-model-step initialization because the default
    # duration equals the configured model timestep.  Reviewer smoke profiles
    # can request a longer, more finely substepped transient without consuming
    # any observation time: the returned end state is still DA timestep zero.
    icesee_kwargs.update({'k': 0})
    model_dt = float(icesee_kwargs.get('dt', 0.2))
    spinup_dt = float(icesee_kwargs.get('ensemble_spinup_dt', model_dt))
    spinup_years = float(
        icesee_kwargs.get('ensemble_spinup_years', spinup_dt)
    )
    if spinup_dt <= 0.0:
        raise ValueError("ensemble_spinup_dt must be positive")
    if spinup_years < spinup_dt:
        raise ValueError(
            "ensemble_spinup_years must be at least ensemble_spinup_dt"
        )
    icesee_kwargs.update({
        'dt': spinup_dt,
        'tinitial': 0.0,
        'tfinal': spinup_years,
    })


    # --- filename for data saving
    fname = 'initialize_ensemble.mat'
    icesee_kwargs.update({'fname': fname})

    #*-----------------------
    Lx = icesee_kwargs.get('Lx', 640e3)
    Ly = icesee_kwargs.get('Ly', 80e3)
    fdim = nd//icesee_kwargs.get('total_state_param_vars', 1)
    nx = icesee_kwargs.get('nx', 32)
    ny = icesee_kwargs.get('ny', 4)
    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    seed_base = icesee_kwargs.get('seed_base', 42)

    #  # -- friction
    sill_friction = icesee_kwargs.get('sill_friction')
    range_friction = icesee_kwargs.get('range_friction')
    mean_friction  = icesee_kwargs.get('mean_friction')
    nugget_friction = icesee_kwargs.get('nugget_friction')

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


    bed_kriging_file = f'{icesee_path}/bed_kriging_results.h5'
    bed_field = _load_bed_from_kriging_file(
        bed_kriging_file,
        fdim,
        ens=ens,
        use_mean=False
    )

    # write the wrong states to a .h5 file to be read by the ISSM model before nurging
    friction_bed_filename = f'{icesee_path}/{data_path}/friction_bed_{ens_id}.h5'
    # with h5py.File(friction_bed_filename, 'w', driver='mpio', comm=comm) as f:
    with h5py.File(friction_bed_filename, 'w') as f:
        # -- write the friction field
        f.create_dataset('coefficient', data=friction_field)
        # -- write the bed field
        f.create_dataset('bed', data=bed_field)
    #*-----------------------


    try:
        # -- call the run_model function to initialize the ensemble members
        ISSM_model(**icesee_kwargs)
    except Exception as e:
        print(f"[ICESEE Initialize ensemble]] Error initializing ensemble: {e}")
        server.kill_matlab_processes()


    #  -- Read data from the ISSM side to be accessed by ICESEE on the python side
    output_filename = f'{icesee_path}/{data_path}/ensemble_out_{ens_id}.h5'
    updated_state = {}
    # with h5py.File(output_filename, 'r', driver='mpio', comm=comm) as f:
    with h5py.File(output_filename, 'r') as f:
        # for key in vec_inputs:
        #     updated_state[key] = f[key][0]
        updated_state['Thickness'] = f['Thickness'][:].reshape(-1, order='F')
        # updated_state['Base'] = f['Base'][0]
        updated_state['Surface'] = f['Surface'][:].reshape(-1, order='F')
        updated_state['Vx'] = f['Vx'][:].reshape(-1, order='F')
        updated_state['Vy'] = f['Vy'][:].reshape(-1, order='F')
        if icesee_kwargs.get('joint_estimation', False):
            updated_state['bed'] = f['bed'][:].reshape(-1, order='F')
            updated_state['coefficient'] = f['coefficient'][:].reshape(-1, order='F')

    os.chdir(icesee_path)

    return updated_state
