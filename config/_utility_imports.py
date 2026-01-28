# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-01-13
# @description: ICESEE model applications utility imports
# =============================================================================

# --- Imports ---
import os
import sys
import h5py
import numpy as np
import warnings
from scipy.stats import norm, multivariate_normal
from tqdm import tqdm
import yaml
import argparse
from argparse import ArgumentParser

# Suppress warnings
warnings.filterwarnings('ignore')

from ICESEE.config.extract_flags import extract_flags, generate_flags_markdown

def get_project_root():
    '''Automatically determines the root of the project.'''
    current_dir = os.path.dirname(os.path.abspath(__file__))  # Get absolute path of the current script
    
    # Traverse upwards until we reach the root of the project (assuming 'src' folder exists at root)
    while not os.path.exists(os.path.join(current_dir, 'src')):
        current_dir = os.path.dirname(current_dir)  # Move one level up
    
    return current_dir

# Get the root of the project
project_root = get_project_root()

# Construct the path to 'src' from the project ro
utils_dir = os.path.join(project_root, 'src', 'utils')
run_model_da_dir = os.path.join(project_root, 'src', 'run_model_da')
config_loader_dir = os.path.join(project_root, 'config')
applications_dir = os.path.join(project_root, 'applications')
parallelization_dir = os.path.join(project_root, 'src', 'parallelization')

# Insert the models directory at the beginning of sys.path
# sys.path.insert(0, models_dir)
sys.path.insert(0, utils_dir)
sys.path.insert(0, run_model_da_dir)
sys.path.insert(0, config_loader_dir)
sys.path.insert(0, parallelization_dir)

# import the necessary modules
from tools import save_arrays_to_h5, extract_datasets_from_h5, save_all_data
from tools import icesee_get_index
from utils import UtilsFunctions
from config_loader import load_yaml_to_dict, get_section

# Check if running in Jupyter notebook (for visualization)
flag_jupyter = False
if 'ipykernel' in sys.modules:
    print('[ICESEE] Running in Jupyter - disabling command line arguments')
    # leave entire routine
    flag_jupyter = True

# =============================================================================
# --- Command Line Arguments ---
if not flag_jupyter:
    # Mapping for execution mode
    execution_modes_str = {
        'default_run': 0,
        'sequential_run': 1,
        'even_distribution': 2
    }
    execution_modes_int = {v: k for k, v in execution_modes_str.items()}  # Reverse mapping

    # CL args.
    parser = ArgumentParser(description='ICESEE: Ice Sheet Parameter and State Estimation model')
    parser.add_argument('--Nens', type=int, required=False, default=1, help='ensemble members')
    parser.add_argument('--verbose', action='store_true', help='verbose output')
    parser.add_argument('--default_run', action='store_true', help='default run')
    parser.add_argument('--sequential_run', action='store_true', help='sequential run')
    parser.add_argument('--even_distribution', action='store_true', help='even distribution')
    parser.add_argument('--data_path', type=str, required=False, default= '_modelrun_datasets', help='folder to save data for single or multiple runs')
    parser.add_argument('execution_mode', type=int, choices=[0, 1, 2], nargs='?', help='Execution mode: 0=default_run, 1=sequential_run, 2=even_distribution')
    parser.add_argument('--model_nprocs', type=int, required = False, default=0, help='number of processors for the coupled model')
    parser.add_argument('-F', '--force-params', type=str, required=False, default='params.yaml', help='Path to YAML parameter file (default: params.yaml)')

    args = parser.parse_args()

    # check if default run arugment is provided
    run_flag = False
    if (args.default_run or args.sequential_run or args.even_distribution):
        run_flag = True

    # Determine execution mode
    selected_mode = 'default_run'  # Default mode

    if args.execution_mode is not None:
        selected_mode = execution_modes_int[args.execution_mode]  # Convert int to string
    else:
        for mode in execution_modes_str.keys():
            if getattr(args, mode):
                selected_mode = mode
                break

    # Set flags explicitly
    args.default_run = (selected_mode == 'default_run')
    args.sequential_run = (selected_mode == 'sequential_run')
    args.even_distribution = (selected_mode == 'even_distribution')

    # Explicit use of parameters
    Nens = int(args.Nens)
    data_path = args.data_path
    model_nprocs = int(args.model_nprocs)
    _verbose = args.verbose
    parameters_file = args.force_params  # Use provided YAML file or default 'params.yaml'

    # Create params dictionary
    params = {
        'Nens': int(args.Nens),
        'default_run': args.default_run,
        'sequential_run': args.sequential_run,
        'even_distribution': args.even_distribution,
        'data_path': args.data_path,
        'model_nprocs': int(args.model_nprocs),
        'verbose': args.verbose,
    }

    # print(f'Execution mode selected: {selected_mode}')
    # print(f'Params: {params}')

    # Log which file is being loaded if verbose
    # if _verbose:
    #     print(f'[ICESEE] Loading parameters from {parameters_file}')

    # Verify if the specified parameters file exists
    if not os.path.exists(parameters_file):
        raise FileNotFoundError(f"Parameter file '{parameters_file}' not found. Please ensure the file exists.")

    # Load parameters from the specified YAML file
    parameters = load_yaml_to_dict(parameters_file)

    physical_params = get_section(parameters, 'physical-parameters')
    modeling_params = get_section(parameters, 'modeling-parameters')
    enkf_params     = get_section(parameters, 'enkf-parameters')

    # --- Ensemble Parameters ---
    params.update({
        'nt': int(float(modeling_params['num_years']) * float(modeling_params['timesteps_per_year'])), # number of time steps
        'dt': 1.0 / float(modeling_params['timesteps_per_year']),
        'num_state_vars': int(float(enkf_params.get('num_state_vars', 1))),
        'num_param_vars': int(float(enkf_params.get('num_param_vars', 0))),
        'number_obs_instants': int(int(float(enkf_params.get('obs_max_time', 1))) / float(enkf_params.get('freq_obs', 1))),
        'inflation_factor': float(enkf_params.get('inflation_factor', 1.0)),
        'freq_obs': float(enkf_params.get('freq_obs', 1)),
        'obs_max_time': int(float(enkf_params.get('obs_max_time', 1))),
        'obs_start_time': float(enkf_params.get('obs_start_time', 1)),
        'localization_flag': bool(enkf_params.get('localization_flag', False)),
        'parallel_flag': enkf_params.get('parallel_flag', 'serial'),
        'n_modeltasks': int(enkf_params.get('n_modeltasks', 1)),
        'execution_flag': int(enkf_params.get('execution_flag', 0)),
        'model_name': enkf_params.get('model_name', 'model'),
        'use_random_fields': bool(enkf_params.get('use_random_fields', False)),
        'execution_mode'   : int(enkf_params.get('execution_mode', 1)),  # 0 -> serial, 1 -> partial parallel_run, 2 -> fully parallel run
        'serial_file_creation': bool(enkf_params.get('serial_file_creation', True)),
        'chunk_size': int(enkf_params.get('chunk_size', 5000)),
        'joint_estimated_params': enkf_params.get('joint_estimated_params', []),
        'coupled_model_datasets_dir': enkf_params.get('coupled_model_datasets', 'data'),
        'vec_inputs': enkf_params['vec_inputs'],
        'collective_threshold': int(enkf_params.get('collective_threshold', 16)), # threshold for switching to collective I/O
    })

    params.update({'batch_size': min(int(enkf_params.get('batch_size', 50)), params['nt'])})  # number of time steps to process in each batch
    
    # --- incase CL args not provided ---
    if Nens == 1:
        params['Nens'] = int(float(enkf_params.get('Nens', 1)))

    if data_path == '_modelrun_datasets':
        params['data_path'] = enkf_params.get('data_path', '_modelrun_datasets')

    if model_nprocs == 0:
        params['model_nprocs'] = enkf_params.get('model_nprocs', 0) 
    
    if run_flag:
        execution_flag = params.get('execution_flag')

        if execution_flag == 1:
            params.update({'sequential_run': True, 'default_run': False})
        elif execution_flag == 2:
            params.update({'even_distribution': True, 'default_run': False})
        else:
            params['default_run'] = True

    #either way update the execution flag
    if params['sequential_run']:
        params['execution_flag'] = 1
    elif params['even_distribution']:
        params['execution_flag'] = 2
    else:
        params['execution_flag'] = 0

    # set run modes
    execution_mode = {
        'serial': 1 if params.get('execution_mode', 0) == 0  else 0,
        'partial': 1 if params.get('execution_mode', 0) == 1  else 0,
        'full': 1 if params.get('execution_mode', 0) == 2  else 0,
    }
    # if none of the above modes is set to True set partial to True
    if not any(execution_mode.values()):
        execution_mode['partial'] = True

    params.update({'mode': execution_mode})
    
    # update for time t
    params['t'] = np.linspace(0, int(float(modeling_params['num_years'])), params['nt'] + 1)

    # get verbose flag
    if args.verbose:
        _verbose = True
    else:
        _verbose  = modeling_params.get('verbose', False)

    # model kwargs
    kwargs = {
        't': params['t'],
        'nt': params['nt'],
        'dt': params['dt'],
        'obs_index': (np.linspace(int(params['freq_obs']/params['dt']), \
                            int(params['obs_max_time']/params['dt']), int(params['number_obs_instants']))).astype(int),
        'joint_estimation': bool(enkf_params.get('joint_estimation', False)),
        'parameter_estimation': bool(enkf_params.get('parameter_estimation', False)),
        'state_estimation': bool(enkf_params.get('state_estimation', False)),
        'joint_estimated_params': enkf_params.get('joint_estimated_params', []),
        'global_analysis': bool(enkf_params.get('global_analysis', True)),
        'local_analysis': bool(enkf_params.get('local_analysis', False)),
        'observed_params':enkf_params.get('observed_params', []),
        'verbose':_verbose,
        'param_ens_spread': enkf_params.get('param_ens_spread', []),
        'data_path': params['data_path'],
        'example_name': modeling_params.get('example_name', params.get('model_name')),
        'length_scale': enkf_params.get('length_scale', []),
        'Q_rho': enkf_params.get('Q_rho', 1.0),
        'generate_synthetic_obs': enkf_params.get('generate_synthetic_obs', True),
        'generate_true_state': enkf_params.get('generate_true_state', True),
        'generate_nurged_state': enkf_params.get('generate_nurged_state', True),
        'use_ensemble_pertubations': enkf_params.get('use_ensemble_pertubations', True),
        'sequential_ensemble_initialization': enkf_params.get('sequential_ensemble_initialization', False),
        'observations_available': enkf_params.get('observations_available', False),
        'obs_data_path': enkf_params.get('obs_data_path', params.get('coupled_model_datasets_dir', 'data') + '/observations_data.h5'),
        'create_ensemble_dataset': enkf_params.get('create_ensemble_dataset', True),
        'restart_enabled': enkf_params.get('restart_enabled', True),
        'force_fresh_start': enkf_params.get('force_fresh_start', False),
        'checkpoint_every': int(enkf_params.get('checkpoint_every', 1)),
        'base_seed': int(enkf_params.get('base_seed', 42)),
        'k_start_override': enkf_params.get('k_start_override', None),
        'ICESEE_PERFORMANCE_TEST': bool(enkf_params.get('ICESEE_PERFORMANCE_TEST', False)), # this is an environment variable
        'h5_file_compression': enkf_params.get('h5_file_compression', None), # e.g., 'gzip' or 'lzf' or 'szip' or None
        'h5_file_compression_level': int(enkf_params.get('h5_file_compression_level', 4)), # 0-9 for gzip, 1-9 for szip, ignored for lzf and None
        'h5_file_chunk_size': int(enkf_params.get('h5_file_chunk_size', 1000)),
        'bed_obs_snapshot':enkf_params.get('bed_obs_snapshot', []),# list of time snapshots to observe bed variables
        'bed_obs_stride':enkf_params.get('bed_obs_stride',None ), # spatial stride in km for bed observations
        'bed_obs_spacing':enkf_params.get('bed_obs_spacing', None), # observation spacing every n grid points {int}
        'bed_obs_indices':enkf_params.get('bed_obs_indices', None), # specific indices to observe {list} (bed subvector indices)
        'bed_obs_mask':enkf_params.get('bed_obs_mask', None), # boolean mask array for bed observations {np.array}
        'initialize_ensemble':enkf_params.get('initialize_ensemble', True),
        'initial_spread_factor': enkf_params.get('initial_spread_factor', 1.0),
        'observed_vars': enkf_params.get('observed_vars', []),
        'vel_idx': int(float(enkf_params.get('vel_idx', 2))),
        'inversion_flag': enkf_params.get('inversion_flag', False),
        'friction_idx': int(float(enkf_params.get('friction_idx', 5))),
        'bed_relaxation_factor': float(enkf_params.get('bed_relaxation_factor', 0.05)), # relaxation factor for bed elevation updates (-1 < factor <= 1) (when bed is not observed)
        'initial_bed_bias': float(enkf_params.get('initial_bed_bias', 0.0015)), # initial bias for bed elevation (in model units)
        'abs_vel_weight': float(enkf_params.get('abs_vel_weight', 1.0)), # weight for absolute velocity in inversion
        'rel_vel_weight': float(enkf_params.get('rel_vel_weight', 1.0)), # weight for relative velocity in inversion
        'tikhonov_regularization_weight': float(enkf_params.get('tikhonov_regularization_weight', 1e-13)), # Tikhonov regularization weight for inversion
        'var_nd': enkf_params.get('var_nd', None), # variable state dimension for each state variable in vec_inputs. Used when state variables have different dimensions
        'scalar_inputs': enkf_params.get('scalar_inputs', []), # list of scalar input variables
    }


    # # update kwargs dictonary with params
    kwargs.update({'physical_params': physical_params})
    kwargs.update({'modeling_params': modeling_params})
    kwargs.update({'enkf_params': enkf_params})

    # -- update the kwargs with physical, modeling and enkf parameters
    kwargs.update(physical_params)
    kwargs.update(modeling_params)
    kwargs.update(enkf_params)

    joint_estimated_params = len(kwargs.get('joint_estimated_params', []))
    if kwargs['joint_estimation']:
        params['total_state_param_vars'] = params['num_state_vars'] + joint_estimated_params
    else:
        params['total_state_param_vars'] = params['num_state_vars']

    # add joint estimation flag to params
    params['joint_estimation'] = kwargs['joint_estimation']

    # unpack standard deviations
    params.update({
        'sig_model': enkf_params.get('sig_model', np.array([0.01])*params['total_state_param_vars']),
        'sig_obs': enkf_params.get('sig_obs', np.array([0.01])*params['total_state_param_vars']),
        'sig_Q': enkf_params.get('sig_Q', np.array([0.01])*params['total_state_param_vars']),
        })
    
    if kwargs['joint_estimation']:
       kwargs.update({'vec_inputs': params['vec_inputs']})
    else:
        params.update({
            'sig_obs': np.array(params['sig_obs'][:params['num_state_vars']]),
            'sig_Q': np.array(params['sig_Q'][:params['num_state_vars']]),
            'sig_model': np.array(params['sig_model'][:params['num_state_vars']]),
        })
        params['vec_inputs'] = params['vec_inputs'][:params['num_state_vars']]
        kwargs.update({'vec_inputs': params['vec_inputs']})

    # --- Observations Parameters ---
    if kwargs.get('observations_available', False):
        # load observation data
        if not os.path.exists(kwargs.get('obs_data_path', 'observations_data.h5')):
            raise FileNotFoundError(f"Observation data file '{kwargs.get('obs_data_path', 'observations_data.h5')}' not found. Please ensure the file exists.")
        # Tell the user to load: kwargs['obs_index'] and params['number_obs_instants']
        print("[ICESEE] Please load 'obs_index' and 'number_obs_instants' from the observation data file into the model dictionary.")
        # obs_data = extract_datasets_from_h5(kwargs.get('obs_data_path', 'observations_data.h5'))
        # kwargs.update({'obs_data': obs_data})
    else:
        # generate observation schedule for synthetic observations
        obs_t, obs_idx, num_observations = UtilsFunctions(params).generate_observation_schedule(**kwargs)
        kwargs['obs_index'] = obs_idx
        params['number_obs_instants'] = num_observations
        kwargs['m_obs'] = num_observations

    kwargs['parallel_flag']       = enkf_params.get('parallel_flag', 'serial')
    kwargs['commandlinerun']      = enkf_params.get('commandlinerun', False)

    #  check available parameters in the obseve_params list that need to be observed 
    params_vec = []
    for i, vars in enumerate(kwargs['vec_inputs']):
        if i >= params['num_state_vars']:
            params_vec.append(vars)

    kwargs['params_vec'] = params_vec
    kwargs.update({'params': params})
    kwargs.update(params)

    import re

    # if re.match(r'\AMPI_model\Z', kwargs.get('parallel_flag'), re.IGNORECASE):
    #     # --- Initialize MPI ---
    #     from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager

    #     icesee_rank, icesee_size, icesee_comm, ens_id = ParallelManager().icesee_mpi_init(params)

    #     # check if _modelrun_datasets exists in path if not create one
    #     _modelrun_datasets = kwargs.get('data_path',None)
    #     if icesee_rank == 0 and not os.path.exists(_modelrun_datasets):
    #         os.makedirs(_modelrun_datasets, exist_ok=True)

    #     #  synchronize the processes
    #     icesee_comm.Barrier()

    # else:
    if not re.match(r'\AMPI_model\Z', kwargs.get('parallel_flag'), re.IGNORECASE):
        icesee_rank = 0
        icesee_size = 1
        icesee_comm = None

        # check if _modelrun_datasets exists in path if not create one
        _modelrun_datasets = kwargs.get('data_path',None)
        if not os.path.exists(_modelrun_datasets):
            os.makedirs(_modelrun_datasets, exist_ok=True)


