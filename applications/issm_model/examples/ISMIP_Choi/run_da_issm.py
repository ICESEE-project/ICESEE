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

# --- ICESEE imports ---
# from ICESEE.config._utility_imports import *
from ICESEE.config._utility_imports import params, kwargs, modeling_params, enkf_params, physical_params,UtilsFunctions
from ICESEE.src.run_model_da.run_models_da import icesee_model_data_assimilation
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager

#  model-specific imports
from ICESEE.applications.issm_model.examples.ISMIP_Choi._issm_model import initialize_model
from ICESEE.applications.issm_model.issm_utils.matlab2python.mat2py_utils import add_issm_dir_to_sys_path, MatlabServer
from ICESEE.applications.issm_model.issm_utils.matlab2python.server_utils import run_icesee_with_server, setup_server_shutdown

# --- Initialize MPI ---
icesee_rank, icesee_size, icesee_comm, ens_id = ParallelManager().icesee_mpi_init(params)

# print(f"[DEBUG] MPI rank: {icesee_rank}, size: {icesee_size} ens_id: {ens_id}")

# --- get current working directory ---
icesee_cwd = os.getcwd()

# --- change directory to issm model directory: make sure ISSM_DIR is set in the environment
issm_dir = os.environ.get('ISSM_DIR')  # make sure ISSM_DIR is set in the environment
add_issm_dir_to_sys_path(issm_dir)     # add the issm directory to the system path 

# --- make the examples directory available ---
issm_examples_dir = os.path.join(issm_dir, 'examples',kwargs.get('example_name'))

# --- fetch the modeling parameters ---
model_kwargs = {
               'Lx': int(float(physical_params.get('Lx'))), 'Ly': int(float(physical_params.get('Ly'))),
                'nx': int(float(physical_params.get('nx'))), 'ny': int(float(physical_params.get('ny'))),
                'ParamFile': modeling_params.get('ParamFile'),
                'cluster_name': socket.gethostname().replace('-', ''),
                'steps': int(float(modeling_params.get('steps'))),
                'dt': float(modeling_params.get('timesteps_per_year')),
                'tinitial': float(modeling_params.get('tinitial')),
                'tfinal': float(modeling_params.get('num_years')),
                't': np.linspace(modeling_params.get('tinitial'), modeling_params.get('num_years'), int((modeling_params.get('num_years') - modeling_params.get('tinitial'))/modeling_params.get('timesteps_per_year'))+1),
                'nt': int((modeling_params.get('num_years') - modeling_params.get('tinitial'))/modeling_params.get('timesteps_per_year')),
                'icesee_path': icesee_cwd,
                'data_path': kwargs.get('data_path'),
                'issm_dir': issm_dir,
                'issm_examples_dir': issm_examples_dir,
                'rank': icesee_rank,
                'nprocs': icesee_size,
                'ens_id': ens_id,
                'hpcmode': modeling_params.get('hpcmode', False),
                'devmode': modeling_params.get('devmode', False),
                'use_reference_data': modeling_params.get('use_reference_data', False),
                'reference_data_dir': modeling_params.get('reference_data_dir', 'data'),
                'reference_data' : modeling_params.get('reference_data'),
                'sill_friction': enkf_params.get('sill_friction', 90000),
                'range_friction': enkf_params.get('range_friction', 5000),
                'mean_friction': enkf_params.get('mean_friction', 2500),
                'nugget_friction': enkf_params.get('nugget_friction', 0),
                'sill_bed': enkf_params.get('sill_bed', 4000),
                'range_bed': enkf_params.get('range_bed', 50000),
                'nugget_bed': enkf_params.get('nugget_bed', 200),      
}

# observation schedule
obs_t, obs_idx, num_observations = UtilsFunctions(params).generate_observation_schedule(**model_kwargs)
model_kwargs["obs_index"] = obs_idx
params["number_obs_instants"] = num_observations

# --- save model kwargs to file and update Icesee kwargs ---
sio.savemat(f'model_kwargs_{ens_id}.mat', model_kwargs)
kwargs.update(model_kwargs)

# copy the issm_env.m from icesee_cwd  file to the examples directory             
shutil.copy(os.path.join(icesee_cwd,'..','..','issm_utils','matlab2python', 'issm_env.m'), issm_examples_dir)
shutil.copy(os.path.join(icesee_cwd,'..','..','issm_utils','matlab2python', 'matlab_server.m'), issm_examples_dir)
shutil.copy(os.path.join(icesee_cwd, f'model_kwargs_{ens_id}.mat'), issm_examples_dir)

# --- change directory to the examples directory ---
os.chdir(issm_examples_dir)

# --- intialize the matlab server ---
server = MatlabServer(color=ens_id,
                      Nens = params['Nens'],
                      comm = icesee_comm,
                      verbose=params.get('verbose')) 

# Set up global shutdown handler
# setup_server_shutdown(server, icesee_comm, verbose=False)

# --- load the model parameters ---
kwargs.update({'server': server, 'Nens': params.get('Nens'), 'icesee_comm': icesee_comm,
                        'icesee_path': icesee_cwd, 'ens_id': ens_id,
                        'data_path': kwargs.get('data_path'),
                        'model_nprocs': params.get('model_nprocs'),})

# --- initialize the model ---
variable_size = initialize_model(**kwargs)

params.update({'nd': variable_size*params.get('total_state_param_vars')})

# --- change directory back to the original directory ---
os.chdir(icesee_cwd)

# --- run the model ---
kwargs.update({'params': params, 
               'server': server})

try:
    icesee_model_data_assimilation(**kwargs)
    server.shutdown()
except Exception as e:
    print(f"[run_da_issm] Error running the model: {e}")
    server.kill_matlab_processes()
    exit()
