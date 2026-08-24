# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-03-24
# @description: ISSM Model with Data Assimilation using a Python Wrapper.
#
# =============================================================================

# --- Imports ---
import sys
import os
import shutil
import socket
import numpy as np
import scipy.io as sio
from pathlib import Path

# --- Set up paths ---
os.chdir(Path(__file__).resolve().parent)

# --- ICESEE imports ---
# from ICESEE.config._utility_imports import *
from ICESEE.config._utility_imports import UtilsFunctions, icesee_kwargs
from ICESEE.src.utils.icesee_context import matlab_icesee_kwargs
from ICESEE.src.run_model_da.run_models_da import icesee_model_data_assimilation
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager

#  model-specific imports
from _issm_model import *
from ICESEE.applications.issm_model.issm_utils.matlab2python.mat2py_utils import add_issm_dir_to_sys_path, MatlabServer, setup_example_directory
from ICESEE.applications.issm_model.issm_utils.matlab2python.server_utils import run_icesee_with_server, setup_server_shutdown

# --- Initialize MPI ---
icesee_rank, icesee_size, icesee_comm, ens_id = ParallelManager().icesee_mpi_init(icesee_kwargs)

# print(f"[DEBUG] MPI rank: {icesee_rank}, size: {icesee_size} ens_id: {ens_id}")

# --- get current working directory ---
icesee_cwd = os.getcwd()

# --- change directory to issm model directory: make sure ISSM_DIR is set in the environment
issm_dir = os.environ.get('ISSM_DIR')  # make sure ISSM_DIR is set in the environment
add_issm_dir_to_sys_path(issm_dir)     # add the issm directory to the system path

# --- make the examples directory available ---
issm_examples_dir = setup_example_directory(issm_dir, icesee_kwargs.get('example_name'))

# --- fetch the modeling parameters ---
icesee_kwargs.update({
               'Lx': int(float(icesee_kwargs.get('Lx'))), 'Ly': int(float(icesee_kwargs.get('Ly'))),
                'nx': int(float(icesee_kwargs.get('nx'))), 'ny': int(float(icesee_kwargs.get('ny'))),
                'ParamFile': icesee_kwargs.get('ParamFile'),
                'cluster_name': socket.gethostname().replace('-', ''),
                'steps': int(float(icesee_kwargs.get('steps'))),
                'dt': float(icesee_kwargs.get('timesteps_per_year')),
                'tinitial': float(icesee_kwargs.get('tinitial')),
                'tfinal': float(icesee_kwargs.get('num_years')),
                't': np.linspace(icesee_kwargs.get('tinitial'), icesee_kwargs.get('num_years'), int((icesee_kwargs.get('num_years') - icesee_kwargs.get('tinitial'))/icesee_kwargs.get('timesteps_per_year'))+1),
                'nt': int((icesee_kwargs.get('num_years') - icesee_kwargs.get('tinitial'))/icesee_kwargs.get('timesteps_per_year')),
                'icesee_path': icesee_cwd,
                'data_path': icesee_kwargs.get('data_path'),
                'issm_dir': issm_dir,
                'issm_examples_dir': issm_examples_dir,
                'rank': icesee_rank,
                'nprocs': icesee_size,
                'ens_id': ens_id,
                'hpcmode': icesee_kwargs.get('hpcmode', False),
                'devmode': icesee_kwargs.get('devmode', False),
                'use_reference_data': icesee_kwargs.get('use_reference_data', False),
                'reference_data_dir': icesee_kwargs.get('reference_data_dir', 'data'),
                'reference_data' : icesee_kwargs.get('reference_data'),
                'sill_friction': icesee_kwargs.get('sill_friction', 90000),
                'range_friction': icesee_kwargs.get('range_friction', 5000),
                'mean_friction': icesee_kwargs.get('mean_friction', 2500),
                'nugget_friction': icesee_kwargs.get('nugget_friction', 0),
                'sill_bed': icesee_kwargs.get('sill_bed', 4000),
                'range_bed': icesee_kwargs.get('range_bed', 50000),
                'nugget_bed': icesee_kwargs.get('nugget_bed', 200),
                'deepwater_melting_rate': float(icesee_kwargs.get('deepwater_melting_rate', 200)),
                'smb': float(icesee_kwargs.get('smb', 0.0)),
                'vel_idx': int(float(icesee_kwargs.get('vel_idx', 2))),
                'inversion_flag': icesee_kwargs.get('inversion_flag', False),
                'friction_idx': int(float(icesee_kwargs.get('friction_idx', 5))),
                'min_friction': float(icesee_kwargs.get('min_friction', 2000)),
                'max_friction': float(icesee_kwargs.get('max_friction', 4000)),
                'Nens': int(float(icesee_kwargs.get('Nens'))),
                'bed_relaxation_factor': float(icesee_kwargs.get('bed_relaxation_factor', 0.05)),
                'initial_bed_bias': float(icesee_kwargs.get('initial_bed_bias', 0.0015)),
                'abs_vel_weight': float(icesee_kwargs.get('abs_vel_weight', 1.0)),
                'rel_vel_weight': float(icesee_kwargs.get('rel_vel_weight', 1.0)),
                'tikhonov_regularization_weight': float(icesee_kwargs.get('tikhonov_regularization_weight', 1e-13)),
                'b_nurge': float(icesee_kwargs.get('b_nurge', 0)),
                's_nurge': float(icesee_kwargs.get('s_nurge', 0)),
                'vec_inputs': icesee_kwargs.get('vec_inputs'), # State vector inputs (ice surface and velocities)
                'scalar_inputs': icesee_kwargs.get('scalar_inputs', []),
})

# observation schedule
obs_t, obs_idx, num_observations = UtilsFunctions(icesee_kwargs).generate_observation_schedule(**icesee_kwargs)
icesee_kwargs["obs_index"] = obs_idx
icesee_kwargs["number_obs_instants"] = num_observations

# --- save model metadata and update the ICESEE runtime context ---
icesee_kwargs_file = f'icesee_kwargs_{ens_id}.mat'
sio.savemat(icesee_kwargs_file, matlab_icesee_kwargs(icesee_kwargs))

# copy the issm_env.m from icesee_cwd  file to the examples directory
shutil.copy(os.path.join(icesee_cwd,'..','..','issm_utils','matlab2python', 'issm_env.m'), issm_examples_dir)
shutil.copy(os.path.join(icesee_cwd,'..','..','issm_utils','matlab2python', 'matlab_server.m'), issm_examples_dir)
shutil.copy(os.path.join(icesee_cwd, icesee_kwargs_file), issm_examples_dir)
shutil.copy(os.path.join(icesee_cwd, f'Domain.exp'), issm_examples_dir)
shutil.copy(os.path.join(icesee_cwd, icesee_kwargs.get('ParamFile')), issm_examples_dir)

# --- change directory to the examples directory ---
os.chdir(issm_examples_dir)

# --- intialize the matlab server ---
server = MatlabServer(color=ens_id,
                      Nens = icesee_kwargs['Nens'],
                      comm = icesee_comm,
                      verbose=icesee_kwargs.get('verbose'))

# Set up global shutdown handler
# setup_server_shutdown(server, icesee_comm, verbose=False)

# --- load the model parameters ---
icesee_kwargs.update({'server': server, 'Nens': icesee_kwargs.get('Nens'), 'icesee_comm': icesee_comm,
                        'icesee_path': icesee_cwd, 'ens_id': ens_id,
                        'data_path': icesee_kwargs.get('data_path'),
                        'model_nprocs': icesee_kwargs.get('model_nprocs'),})

# --- initialize the model ---
variable_size = initialize_model(**icesee_kwargs)

icesee_kwargs.update({'nd': variable_size*icesee_kwargs.get('total_state_param_vars')})
var_nd = {var: (1 if var in icesee_kwargs['scalar_inputs'] else variable_size) for var in icesee_kwargs['vec_inputs']}

# --- change directory back to the original directory ---
os.chdir(icesee_cwd)

# --- run the model ---
icesee_kwargs.update({'var_nd': var_nd,
               'nd': icesee_kwargs.get('nd'),
               'server': server})


try:
    icesee_model_data_assimilation(**icesee_kwargs)
    server.shutdown()
except Exception as e:
    print(f"[run_da_issm] Error running the model: {e}")
    server.kill_matlab_processes()
    exit()
