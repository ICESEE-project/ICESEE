# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-05-26
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
from ICESEE.config._utility_imports import UtilsFunctions, icesee_kwargs
from ICESEE.src.utils.icesee_context import matlab_icesee_kwargs
from ICESEE.src.run_model_da.run_models_da import icesee_model_data_assimilation
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager

#  model-specific imports
from ICESEE.applications.issm_model.examples.ISMIP._issm_model import initialize_model
from ICESEE.applications.issm_model.issm_utils.matlab2python.mat2py_utils import add_issm_dir_to_sys_path, MatlabServer
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
issm_examples_dir = os.path.join(issm_dir, 'examples',icesee_kwargs.get('example_name'))

# --- fetch the modeling parameters ---
icesee_kwargs.update({
               'Lx': int(float(icesee_kwargs.get('Lx'))), 'Ly': int(float(icesee_kwargs.get('Ly'))),
                'nx': int(float(icesee_kwargs.get('nx'))), 'ny': int(float(icesee_kwargs.get('ny'))),
                'ParamFile': icesee_kwargs.get('ParamFile'),
                'cluster_name': socket.gethostname().replace('-', ''),
                'extrusion_layers': int(float(icesee_kwargs.get('extrusion_layers'))),
                'extrusion_exponent': int(float(icesee_kwargs.get('extrusion_exponent'))),
                'steps': int(float(icesee_kwargs.get('steps'))),
                'flow_model': icesee_kwargs.get('flow_model'),
                'sliding_vx': float(icesee_kwargs.get('sliding_vx')),
                'sliding_vy': float(icesee_kwargs.get('sliding_vy')),
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

# --- change directory to the examples directory ---
os.chdir(issm_examples_dir)

# --- intialize the matlab server ---
server = MatlabServer(color=ens_id,
                      Nens = icesee_kwargs['Nens'],
                      comm = icesee_comm,
                       verbose=icesee_kwargs.get('verbose'))

# Set up global shutdown handler
setup_server_shutdown(server, icesee_comm, verbose=False)

# --- load the model parameters ---
icesee_kwargs.update({'server': server, 'Nens': icesee_kwargs.get('Nens'), 'icesee_comm': icesee_comm,
                        'icesee_path': icesee_cwd, 'ens_id': ens_id,
                        'data_path': icesee_kwargs.get('data_path'),
                        'model_nprocs': icesee_kwargs.get('model_nprocs'),})

# --- initialize the model ---
variable_size = initialize_model(**icesee_kwargs)

icesee_kwargs.update({'nd': variable_size*icesee_kwargs.get('total_state_param_vars')})

# --- change directory back to the original directory ---
os.chdir(icesee_cwd)

# --- run the model ---
icesee_kwargs.update({'server': server})


if False:
    try:
        result = run_icesee_with_server(
            icesee_model_data_assimilation(
            icesee_kwargs["model_name"],
            icesee_kwargs["filter_type"],
            **icesee_kwargs), server, True,icesee_comm,verbose=True
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
    #     icesee_kwargs["model_name"],
    #     icesee_kwargs["filter_type"],
    #     **icesee_kwargs), server, False,icesee_comm,verbose=False
    # )
    try:
        icesee_model_data_assimilation(**icesee_kwargs)
        server.shutdown()
    except Exception as e:
        print(f"[run_da_issm] Error running the model: {e}")
        server.kill_matlab_processes()
        exit()
#     print("Checking stdout:", sys.stdout, file=sys.stderr)  # Use stderr to avoid stdout issues
# sys.stdout.flush()
