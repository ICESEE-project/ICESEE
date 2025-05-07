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
from ICESEE.applications.issm_model.examples.ISMIP._issm_model import initialize_model
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
                'extrusion_layers': int(float(modeling_params.get('extrusion_layers'))),
                'extrusion_exponent': int(float(modeling_params.get('extrusion_exponent'))),
                'steps': int(float(modeling_params.get('steps'))),
                'flow_model': modeling_params.get('flow_model'),
                'sliding_vx': float(modeling_params.get('sliding_vx')),
                'sliding_vy': float(modeling_params.get('sliding_vy')),
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
}

# observation schedule
obs_t, obs_idx, num_observations = UtilsFunctions(params).generate_observation_schedule(**model_kwargs)
model_kwargs["obs_index"] = obs_idx
params["number_obs_instants"] = num_observations

# --- save model kwargs to file and update Icesee kwargs ---
sio.savemat('model_kwargs.mat', model_kwargs)
kwargs.update(model_kwargs)

# copy the issm_env.m from icesee_cwd  file to the examples directory             
shutil.copy(os.path.join(icesee_cwd,'..','..','issm_utils','matlab2python', 'issm_env.m'), issm_examples_dir)
shutil.copy(os.path.join(icesee_cwd,'..','..','issm_utils','matlab2python', 'matlab_server.m'), issm_examples_dir)
shutil.copy(os.path.join(icesee_cwd, 'model_kwargs.mat'), issm_examples_dir)

# --- change directory to the examples directory ---
os.chdir(issm_examples_dir)

# --- intialize the matlab server ---
server = MatlabServer(color=ens_id, verbose=1)
server.launch() # start the server

# Set up global shutdown handler
setup_server_shutdown(server, icesee_comm, verbose=False)

# --- intialize ISSM model ---
modeling_params.update({'server': server, 'Nens': params.get('Nens'),
                        'icesee_path': icesee_cwd, 'ens_id': ens_id,
                        'data_path': kwargs.get('data_path')})
# if icesee_rank == 0:
#     variable_size = initialize_model(physical_params, modeling_params, icesee_comm)
# else:
#     variable_size = 0.0

# # wait for rank 0 to write to file before proceeding
# icesee_comm.Barrier()
# variable_size = icesee_comm.bcast(variable_size, root=0)

if icesee_rank == 0:
    variable_size = initialize_model(physical_params, modeling_params, icesee_comm)
else:
    variable_size = 0.0

# wait for rank 0 to write to file before proceeding
icesee_comm.Barrier()
variable_size = icesee_comm.bcast(variable_size, root=0)

params.update({'nd': variable_size*params.get('total_state_param_vars')})

# --- change directory back to the original directory ---
os.chdir(icesee_cwd)

# --- run the model ---
kwargs.update({'params': params, 
               'server': server})


if False:
    try:
        result = run_icesee_with_server(
            icesee_model_data_assimilation(
            enkf_params["model_name"],
            enkf_params["filter_type"],
            **kwargs), server, True,icesee_comm,verbose=True
        )
    except Exception as e:
        print(f"[DEBUG] Error running the model: {e}")
        result = None
    finally:
        try:
            server.shutdown()
            server.reset_terminal()
        except Exception as e:
            print(f"[DEBUG] Error shutting down server: {e}")
        sys.exit(1)
else:
    # result = run_icesee_with_server(
    #     icesee_model_data_assimilation(
    #     enkf_params["model_name"],
    #     enkf_params["filter_type"],
    #     **kwargs), server, False,icesee_comm,verbose=False
    # )
    # try:
        icesee_model_data_assimilation(
            enkf_params["model_name"],
            enkf_params["filter_type"],
            **kwargs)
        server.shutdown()
    # except Exception as e:
    #     print(f"[run_da_issm] Error running the model: {e}")
    #     result = None
    # server.kill_matlab_processes()
#     print("Checking stdout:", sys.stdout, file=sys.stderr)  # Use stderr to avoid stdout issues
# sys.stdout.flush()