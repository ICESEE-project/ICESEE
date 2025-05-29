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
warnings.filterwarnings("ignore")

def get_project_root():
    """Automatically determines the root of the project."""
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
    print("Running in Jupyter - disabling command line arguments")
    # leave entire routine
    flag_jupyter = True

# =============================================================================
# --- Command Line Arguments ---
if not flag_jupyter:
    # Mapping for execution mode
    execution_modes_str = {
        "default_run": 0,
        "sequential_run": 1,
        "even_distribution": 2
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
    parser.add_argument('--model_nprocs', type=int, required = False, default=0, help='number of processors for the model')

    args = parser.parse_args()

    # check if default run arugment is provided
    run_flag = False
    if (args.default_run or args.sequential_run or args.even_distribution):
        run_flag = True

    # Determine execution mode
    selected_mode = "default_run"  # Default mode

    if args.execution_mode is not None:
        selected_mode = execution_modes_int[args.execution_mode]  # Convert int to string
    else:
        for mode in execution_modes_str.keys():
            if getattr(args, mode):
                selected_mode = mode
                break

    # Set flags explicitly
    args.default_run = (selected_mode == "default_run")
    args.sequential_run = (selected_mode == "sequential_run")
    args.even_distribution = (selected_mode == "even_distribution")

    # Explicit use of parameters
    Nens = int(args.Nens)
    data_path = args.data_path
    model_nprocs = int(args.model_nprocs)
    _verbose = args.verbose

    # Create params dictionary
    params = {
        "Nens": int(args.Nens),
        "default_run": args.default_run,
        "sequential_run": args.sequential_run,
        "even_distribution": args.even_distribution,
        "data_path": args.data_path,
        "model_nprocs": int(args.model_nprocs),
        "verbose": args.verbose,
    }

    # print(f"Execution mode selected: {selected_mode}")
    # print(f"Params: {params}")

    # Load parameters from a YAML file
    parameters_file = "params.yaml"
    parameters = load_yaml_to_dict(parameters_file)

    physical_params = get_section(parameters, "physical-parameters")
    modeling_params = get_section(parameters, "modeling-parameters")
    enkf_params     = get_section(parameters, "enkf-parameters")

    # --- Ensemble Parameters ---
    params.update({
        "nt": int(float(modeling_params["num_years"])) * int(float(modeling_params["timesteps_per_year"])), # number of time steps
        "dt": 1.0 / float(modeling_params["timesteps_per_year"]),
        "num_state_vars": int(float(enkf_params.get("num_state_vars", 1))),
        "num_param_vars": int(float(enkf_params.get("num_param_vars", 0))),
        "number_obs_instants": int(int(float(enkf_params.get("obs_max_time", 1))) / float(enkf_params.get("freq_obs", 1))),
        "inflation_factor": float(enkf_params.get("inflation_factor", 1.0)),
        "freq_obs": float(enkf_params.get("freq_obs", 1)),
        "obs_max_time": int(float(enkf_params.get("obs_max_time", 1))),
        "obs_start_time": float(enkf_params.get("obs_start_time", 1)),
        "localization_flag": bool(enkf_params.get("localization_flag", False)),
        "parallel_flag": enkf_params.get("parallel_flag", "serial"),
        "n_modeltasks": int(enkf_params.get("n_modeltasks", 1)),
        "execution_flag": int(enkf_params.get("execution_flag", 0)),
        "model_name": enkf_params.get("model_name", "model"),
    })

    # --- incase CL args not provided ---
    if Nens == 1:
        params["Nens"] = int(float(enkf_params.get("Nens", 1)))

    if data_path == '_modelrun_datasets':
        params["data_path"] = enkf_params.get("data_path", "_modelrun_datasets")

    if model_nprocs == 0:
        params["model_nprocs"] = enkf_params.get("model_nprocs", 1)
    
    if run_flag:
        execution_flag = params.get("execution_flag")

        if execution_flag == 1:
            params.update({"sequential_run": True, "default_run": False})
        elif execution_flag == 2:
            params.update({"even_distribution": True, "default_run": False})
        else:
            params["default_run"] = True

    #either way update the execution flag
    if params["sequential_run"]:
        params["execution_flag"] = 1
    elif params["even_distribution"]:
        params["execution_flag"] = 2
    else:
        params["execution_flag"] = 0
    
    # update for time t
    params["t"] = np.linspace(0, int(float(modeling_params["num_years"])), params["nt"] + 1)

    # get verbose flag
    if args.verbose:
        _verbose = True
    else:
        _verbose  = modeling_params.get("verbose", False)

    # model kwargs
    kwargs = {
        "t": params["t"],
        "nt": params["nt"],
        "dt": params["dt"],
        "obs_index": (np.linspace(int(params["freq_obs"]/params["dt"]), \
                            int(params["obs_max_time"]/params["dt"]), int(params["number_obs_instants"]))).astype(int),
        "joint_estimation": bool(enkf_params.get("joint_estimation", False)),
        "parameter_estimation": bool(enkf_params.get("parameter_estimation", False)),
        "state_estimation": bool(enkf_params.get("state_estimation", False)),
        "vec_inputs": enkf_params['vec_inputs'],
        "global_analysis": bool(enkf_params.get("global_analysis", True)),
        "local_analysis": bool(enkf_params.get("local_analysis", False)),
        "observed_params":enkf_params.get("observed_params", []),
        "verbose":_verbose,
        "param_ens_spread": enkf_params.get("param_ens_spread", []),
        "data_path": params["data_path"],
        'example_name': modeling_params.get('example_name', params.get('model_name')),
        'length_scale': enkf_params.get("length_scale", []),
        'Q_rho': enkf_params.get("Q_rho", 1.0),
    }


    if kwargs["joint_estimation"]:
        params["total_state_param_vars"] = params["num_state_vars"] + params["num_param_vars"]
    else:
        params["total_state_param_vars"] = params["num_state_vars"]
    # params["total_state_param_vars"] = params["num_state_vars"] + params["num_param_vars"]

    # add joint estimation flag to params
    params["joint_estimation"] = kwargs["joint_estimation"]

    # unpack standard deviations
    params.update({
        "sig_model": enkf_params.get("sig_model", np.array([0.01])*params["total_state_param_vars"]),
        "sig_obs": enkf_params.get("sig_obs", np.array([0.01])*params["total_state_param_vars"]),
        "sig_Q": enkf_params.get("sig_Q", np.array([0.01])*params["total_state_param_vars"]),
        })

    # --- Observations Parameters ---
    obs_t, obs_idx, num_observations = UtilsFunctions(params).generate_observation_schedule(**kwargs)
    kwargs["obs_index"] = obs_idx
    params["number_obs_instants"] = num_observations
    kwargs["parallel_flag"]       = enkf_params.get("parallel_flag", "serial")
    kwargs["commandlinerun"]      = enkf_params.get("commandlinerun", False)

    #  check available parameters in the obseve_params list that need to be observed 
    params_vec = []
    for i, vars in enumerate(kwargs["vec_inputs"]):
        if i >= params["num_state_vars"]:
            params_vec.append(vars)

    kwargs["params_vec"] = params_vec

    # update kwargs dictonary with params
    kwargs.update({'params': params})
    kwargs.update({'physical_params': physical_params})
    kwargs.update({'modeling_params': modeling_params})
    kwargs.update({'enkf_params': enkf_params})

    # -- update the kwargs with physical, modeling and enkf parameters
    kwargs.update(physical_params)
    kwargs.update(modeling_params)
    kwargs.update(enkf_params)
    kwargs.update(params)
    
    import re

    # if re.match(r"\AMPI_model\Z", kwargs.get('parallel_flag'), re.IGNORECASE):
    #     # --- Initialize MPI ---
    #     from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager

    #     icesee_rank, icesee_size, icesee_comm, ens_id = ParallelManager().icesee_mpi_init(params)

    #     # check if _modelrun_datasets exists in path if not create one
    #     _modelrun_datasets = kwargs.get("data_path",None)
    #     if icesee_rank == 0 and not os.path.exists(_modelrun_datasets):
    #         os.makedirs(_modelrun_datasets, exist_ok=True)

    #     #  synchronize the processes
    #     icesee_comm.Barrier()

    # else:
    if not re.match(r"\AMPI_model\Z", kwargs.get('parallel_flag'), re.IGNORECASE):
        icesee_rank = 0
        icesee_size = 1
        icesee_comm = None

        # check if _modelrun_datasets exists in path if not create one
        _modelrun_datasets = kwargs.get("data_path",None)
        if not os.path.exists(_modelrun_datasets):
            os.makedirs(_modelrun_datasets, exist_ok=True)
