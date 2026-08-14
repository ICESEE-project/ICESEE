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
    # Mapping for the intra-mode ensemble distribution strategy.  This is
    # intentionally separate from ``execution_mode`` (0/1/2 in YAML).
    distribution_modes_str = {
        'default_run': 0,
        'sequential_run': 1,
        'even_distribution': 2
    }
    distribution_modes_int = {v: k for k, v in distribution_modes_str.items()}

    # CL args.
    parser = ArgumentParser(description='ICESEE: Ice Sheet Parameter and State Estimation model')
    # ``None`` means "not supplied on the command line".  Do not use a valid
    # runtime value (for example Nens=1) as a sentinel: after the configuration
    # sections are flattened into ``icesee_kwargs`` we must still know which
    # values the user explicitly requested on the command line.
    parser.add_argument('--Nens', type=int, required=False, default=None, help='ensemble members')
    parser.add_argument('--verbose', action='store_true', help='verbose output')
    parser.add_argument('--default_run', action='store_true', help='default run')
    parser.add_argument('--sequential_run', action='store_true', help='sequential run')
    parser.add_argument('--even_distribution', action='store_true', help='even distribution')
    parser.add_argument('--data_path', type=str, required=False, default=None, help='folder to save data for single or multiple runs')
    parser.add_argument('distribution_mode', type=int, choices=[0, 1, 2], nargs='?', help='Ensemble distribution: 0=default, 1=sequential, 2=even')
    parser.add_argument('--model_nprocs', type=int, required = False, default=None, help='number of processors for the coupled model')
    parser.add_argument('-F', '--force-params', type=str, required=False, default='params.yaml', help='Path to YAML parameter file (default: params.yaml)')

    args = parser.parse_args()

    # check if default run arugment is provided
    run_flag = False
    if (args.default_run or args.sequential_run or args.even_distribution):
        run_flag = True

    # Determine the distribution strategy within the selected execution mode.
    selected_mode = 'default_run'  # Default mode

    if args.distribution_mode is not None:
        selected_mode = distribution_modes_int[args.distribution_mode]
    else:
        for mode in distribution_modes_str.keys():
            if getattr(args, mode):
                selected_mode = mode
                break

    # Set flags explicitly
    args.default_run = (selected_mode == 'default_run')
    args.sequential_run = (selected_mode == 'sequential_run')
    args.even_distribution = (selected_mode == 'even_distribution')

    # Explicit use of parameters
    Nens = int(args.Nens) if args.Nens is not None else None
    data_path = args.data_path
    model_nprocs = int(args.model_nprocs) if args.model_nprocs is not None else None
    _verbose = args.verbose
    parameters_file = args.force_params  # Use provided YAML file or default 'params.yaml'

    # Canonical ICESEE runtime context. All command-line, YAML-derived, and
    # runtime values are packed into this single flat dictionary.
    icesee_kwargs = {
        'Nens': Nens,
        'default_run': args.default_run,
        'sequential_run': args.sequential_run,
        'even_distribution': args.even_distribution,
        'data_path': data_path,
        'model_nprocs': model_nprocs,
        'verbose': args.verbose,
    }

    # print(f'Execution mode selected: {selected_mode}')
    # print(f'Params: {icesee_kwargs}')

    # Log which file is being loaded if verbose
    # if _verbose:
    #     print(f'[ICESEE] Loading parameters from {parameters_file}')

    # Verify if the specified parameters file exists
    if not os.path.exists(parameters_file):
        raise FileNotFoundError(f"Parameter file '{parameters_file}' not found. Please ensure the file exists.")

    # Load parameters from the specified YAML file
    _parameters = load_yaml_to_dict(parameters_file)

    _physical_section = get_section(_parameters, 'physical-parameters')
    _modeling_section = get_section(_parameters, 'modeling-parameters')
    _enkf_section     = get_section(_parameters, 'enkf-parameters')

    # Command-line values have the highest precedence.  Keep the resolved
    # values separately because the raw YAML sections are flattened below and
    # may contain keys with the same names.
    Nens = int(args.Nens) if args.Nens is not None else int(float(_enkf_section.get('Nens', 1)))
    data_path = args.data_path if args.data_path is not None else _enkf_section.get('data_path', '_modelrun_datasets')
    model_nprocs = int(args.model_nprocs) if args.model_nprocs is not None else int(float(_enkf_section.get('model_nprocs', 0)))
    icesee_kwargs.update({
        'Nens': Nens,
        'data_path': data_path,
        'model_nprocs': model_nprocs,
    })

    # --- Ensemble Parameters ---
    icesee_kwargs.update({
        'nt': int(float(_modeling_section['num_years']) * float(_modeling_section['timesteps_per_year'])), # number of time steps
        'dt': 1.0 / float(_modeling_section['timesteps_per_year']),
        'num_state_vars': int(float(_enkf_section.get('num_state_vars', 1))),
        'num_param_vars': int(float(_enkf_section.get('num_param_vars', 0))),
        'number_obs_instants': int(int(float(_enkf_section.get('obs_max_time', 1))) / float(_enkf_section.get('freq_obs', 1))),
        'inflation_factor': float(_enkf_section.get('inflation_factor', 1.0)),
        'state_inflation_factor': float(_enkf_section.get('state_inflation_factor', float(_enkf_section.get('inflation_factor', 1.0)))),
        'param_inflation_factor': float(_enkf_section.get('param_inflation_factor', float(_enkf_section.get('inflation_factor', 1.0)))),
        'bed_inflation_factor': float(_enkf_section.get('bed_inflation_factor', int(float(_enkf_section.get('inflation_factor', 1.0))))),
        'freq_obs': float(_enkf_section.get('freq_obs', 1)),
        'obs_max_time': int(float(_enkf_section.get('obs_max_time', 1))),
        'obs_start_time': float(_enkf_section.get('obs_start_time', 1)),
        'localization_flag': bool(_enkf_section.get('localization_flag', False)),
        'n_modeltasks': int(_enkf_section.get('n_modeltasks', 1)),
        'execution_flag': int(_enkf_section.get('execution_flag', 0)),
        'model_name': _enkf_section.get('model_name', 'model'),
        'use_random_fields': bool(_enkf_section.get('use_random_fields', False)),
        'execution_mode'   : int(_enkf_section.get('execution_mode', 1)),  # 0 -> serial, 1 -> partial parallel_run, 2 -> fully parallel run
        'serial_file_creation': bool(_enkf_section.get('serial_file_creation', True)),
        'chunk_size': int(_enkf_section.get('chunk_size', 5000)),
        'joint_estimated_params': _enkf_section.get('joint_estimated_params', []),
        'coupled_model_datasets_dir': _enkf_section.get('coupled_model_datasets', 'data'),
        'vec_inputs': _enkf_section['vec_inputs'],
        'collective_threshold': int(_enkf_section.get('collective_threshold', 16)), # threshold for switching to collective I/O
    })

    # Mode 2 uses a two-shard sliding window by default: the current and next
    # ensemble states remain open without reserving the full run history.
    # ``batch_size`` is consumed by execution mode 2.  The execution mode is
    # commonly selected by the launcher rather than stored in the YAML, so it
    # is not reliable to derive this default from ``_enkf_section`` here.
    # Keep the file window small unless the user explicitly opts into a wider
    # one; otherwise a mode-2 run can pre-create a large fraction of its full
    # ensemble history before the first forecast is evaluated.
    default_batch_size = 2
    icesee_kwargs.update({'batch_size': min(int(_enkf_section.get('batch_size', default_batch_size)), icesee_kwargs['nt'])})
    icesee_kwargs.update({
        'synthetic_observation_storage': str(_enkf_section.get('synthetic_observation_storage', 'compact')),
        'ensemble_storage_dtype': str(_enkf_section.get('ensemble_storage_dtype', 'float64')),
        'storage_safety_factor': float(_enkf_section.get('storage_safety_factor', 1.10)),
        'fail_on_insufficient_storage': bool(_enkf_section.get('fail_on_insufficient_storage', False)),
    })

    # Spatial random-field backend.  Keep FFT as the backward-compatible
    # default; graph mode is selected explicitly and consumes the physical
    # node coordinates registered by each model application.
    random_field_method = str(
        _enkf_section.get(
            'random_field_method',
            _enkf_section.get('enkf_field_method', 'fft'),
        )
    ).strip().lower()
    if random_field_method not in {'fft', 'graph'}:
        raise ValueError(
            "enkf-parameters.random_field_method must be either 'fft' or 'graph'; "
            f"got {random_field_method!r}"
        )
    icesee_kwargs['random_field_method'] = random_field_method

    if run_flag:
        execution_flag = icesee_kwargs.get('execution_flag')

        if execution_flag == 1:
            icesee_kwargs.update({'sequential_run': True, 'default_run': False})
        elif execution_flag == 2:
            icesee_kwargs.update({'even_distribution': True, 'default_run': False})
        else:
            icesee_kwargs['default_run'] = True

    #either way update the execution flag
    if icesee_kwargs['sequential_run']:
        icesee_kwargs['execution_flag'] = 1
    elif icesee_kwargs['even_distribution']:
        icesee_kwargs['execution_flag'] = 2
    else:
        icesee_kwargs['execution_flag'] = 0

    # update for time t
    icesee_kwargs['t'] = np.linspace(0, int(float(_modeling_section['num_years'])), icesee_kwargs['nt'] + 1)

    # get verbose flag
    if args.verbose:
        _verbose = True
    else:
        _verbose  = _modeling_section.get('verbose', False)

    # Add model and analysis options to the canonical runtime context.
    icesee_kwargs.update({
        't': icesee_kwargs['t'],
        'nt': icesee_kwargs['nt'],
        'dt': icesee_kwargs['dt'],
        'obs_index': (np.linspace(int(icesee_kwargs['freq_obs']/icesee_kwargs['dt']), \
                            int(icesee_kwargs['obs_max_time']/icesee_kwargs['dt']), int(icesee_kwargs['number_obs_instants']))).astype(int),
        'joint_estimation': bool(_enkf_section.get('joint_estimation', False)),
        'parameter_estimation': bool(_enkf_section.get('parameter_estimation', False)),
        'state_estimation': bool(_enkf_section.get('state_estimation', False)),
        'joint_estimated_params': _enkf_section.get('joint_estimated_params', []),
        'global_analysis': bool(_enkf_section.get('global_analysis', True)),
        'local_analysis': bool(_enkf_section.get('local_analysis', False)),
        'enkf_observation_error_mode': str(_enkf_section.get(
            'enkf_observation_error_mode',
            'legacy_prior_anomalies' if _enkf_section.get('use_ensemble_pertubations', True) else 'generated_R'
        )),
        'observed_params':_enkf_section.get('observed_params', []),
        'verbose':_verbose,
        'param_ens_spread': _enkf_section.get('param_ens_spread', []),
        'data_path': icesee_kwargs['data_path'],
        'example_name': _modeling_section.get('example_name', icesee_kwargs.get('model_name')),
        'length_scale': _enkf_section.get('length_scale', []),
        'random_field_method': random_field_method,
        'Q_rho': _enkf_section.get('Q_rho', 1.0),
        'generate_synthetic_obs': _enkf_section.get('generate_synthetic_obs', True),
        'generate_true_state': _enkf_section.get('generate_true_state', True),
        'generate_nurged_state': _enkf_section.get('generate_nurged_state', True),
        'use_ensemble_pertubations': _enkf_section.get('use_ensemble_pertubations', True),
        'sequential_ensemble_initialization': _enkf_section.get('sequential_ensemble_initialization', False),
        'observations_available': _enkf_section.get('observations_available', False),
        'obs_data_path': _enkf_section.get('obs_data_path', icesee_kwargs.get('coupled_model_datasets_dir', 'data') + '/observations_data.h5'),
        'create_ensemble_dataset': _enkf_section.get('create_ensemble_dataset', True),
        # Full-parallel large-data controls.  VDS exposes one logical 3-D
        # ensemble without copying the timestep shards.
        'ensemble_finalize_mode': str(_enkf_section.get('ensemble_finalize_mode', 'vds')),
        # ``auto`` preserves full history when it fits and switches to a
        # restartable rolling window when it does not.
        'ensemble_history_mode': str(_enkf_section.get('ensemble_history_mode', 'auto')),
        'analysis_memory_budget_mb': float(_enkf_section.get('analysis_memory_budget_mb', 256.0)),
        'analysis_row_chunk_size': int(_enkf_section.get('analysis_row_chunk_size', 0)),
        'observation_row_chunk_size': int(_enkf_section.get('observation_row_chunk_size', 0)),
        'analysis_finalize_memory_budget_mb': float(_enkf_section.get('analysis_finalize_memory_budget_mb', 2048.0)),
        'analysis_finalize_member_chunk_size': int(_enkf_section.get('analysis_finalize_member_chunk_size', 0)),
        'local_analysis_memory_budget_mb': float(_enkf_section.get('local_analysis_memory_budget_mb', 2048.0)),
        'observation_operator_storage': str(_enkf_section.get('observation_operator_storage', 'indices')),
        'analysis_svd_energy': float(_enkf_section.get('analysis_svd_energy', 0.999)),
        'finalize_row_chunk_size': int(_enkf_section.get('finalize_row_chunk_size', 16384)),
        'restart_enabled': _enkf_section.get('restart_enabled', True),
        'force_fresh_start': _enkf_section.get('force_fresh_start', False),
        'checkpoint_every': int(_enkf_section.get('checkpoint_every', 1)),
        'base_seed': int(_enkf_section.get('base_seed', 42)),
        # Keep process-noise timing independent of execution mode.  The
        # default matches the established mode-1 schedule.
        'process_noise_schedule': str(_enkf_section.get('process_noise_schedule', 'observations')),
        'k_start_override': _enkf_section.get('k_start_override', None),
        'ICESEE_PERFORMANCE_TEST': bool(_enkf_section.get('ICESEE_PERFORMANCE_TEST', False)), # this is an environment variable
        'h5_file_compression': _enkf_section.get('h5_file_compression', None), # e.g., 'gzip' or 'lzf' or 'szip' or None
        'h5_file_compression_level': int(_enkf_section.get('h5_file_compression_level', 4)), # 0-9 for gzip, 1-9 for szip, ignored for lzf and None
        'h5_file_chunk_size': int(_enkf_section.get('h5_file_chunk_size', 1000)),
        'bed_obs_snapshot':_enkf_section.get('bed_obs_snapshot', []),# list of time snapshots to observe bed variables
        'bed_obs_stride':_enkf_section.get('bed_obs_stride',None ), # spatial stride in km for bed observations
        'bed_obs_track_half_width_m': float(_enkf_section.get('bed_obs_track_half_width_m', 1000.0)), # half-width of cross-flow radar-track sampling band
        'bed_obs_spacing':_enkf_section.get('bed_obs_spacing', None), # observation spacing every n grid points {int}
        'bed_obs_indices':_enkf_section.get('bed_obs_indices', None), # specific indices to observe {list} (bed subvector indices)
        'bed_obs_mask':_enkf_section.get('bed_obs_mask', None), # boolean mask array for bed observations {np.array}
        'bed_update_domain': str(_enkf_section.get('bed_update_domain', 'all')),
        'initialize_ensemble':_enkf_section.get('initialize_ensemble', True),
        'initial_spread_factor': _enkf_section.get('initial_spread_factor', 1.0),
        'observed_vars': _enkf_section.get('observed_vars', []),
        'vel_idx': int(float(_enkf_section.get('vel_idx', 2))),
        'inversion_flag': _enkf_section.get('inversion_flag', False),
        'friction_idx': int(float(_enkf_section.get('friction_idx', 5))),
        'bed_relaxation_factor': float(_enkf_section.get('bed_relaxation_factor', 0.05)), # relaxation factor for bed elevation updates (-1 < factor <= 1) (when bed is not observed)
        'initial_bed_bias': float(_enkf_section.get('initial_bed_bias', 0.0015)), # initial bias for bed elevation (in model units)
        'abs_vel_weight': float(_enkf_section.get('abs_vel_weight', 1.0)), # weight for absolute velocity in inversion
        'rel_vel_weight': float(_enkf_section.get('rel_vel_weight', 1.0)), # weight for relative velocity in inversion
        'tikhonov_regularization_weight': float(_enkf_section.get('tikhonov_regularization_weight', 1e-13)), # Tikhonov regularization weight for inversion
        'var_nd': _enkf_section.get('var_nd', None), # variable state dimension for each state variable in vec_inputs. Used when state variables have different dimensions
        'scalar_inputs': _enkf_section.get('scalar_inputs', []), # list of scalar input variables
        'generate_true_wrong_state_only': _enkf_section.get('generate_true_wrong_state_only', False), # flag to only generate true and wrong state without running the assimilation
        'initial_state_only': _enkf_section.get('initial_state_only', False),
        'generate_synthetic_obs_only': bool(_enkf_section.get('generate_synthetic_obs_only', False)), # flag to only generate synthetic observations without running the assimilation
        'localized_vars': _enkf_section.get('localized_vars', []), # list of variables to localize (only used if localization_flag is True)
        'frozen_analysis_vars': _enkf_section.get('frozen_analysis_vars', []), # state-vector blocks held fixed during the analysis update
        'localization_radius': _enkf_section.get('localization_radius', None), # localization radius (float or dict {var_name: radius})
        'node_coords': _enkf_section.get('node_coords', {}), # dict {var_name: (n_i,2) node coordinates}
        'obs_node_coords': _enkf_section.get('obs_node_coords', {}), # dict {var_name: (m_i,2) active obs-node coords, static/union across the run}
        'taper_type': _enkf_section.get('taper_type', 'gaspari_cohn'), # 'gaspari_cohn' (default) or 'gaussian'
        'partitioned_io_flag': _enkf_section.get('partitioned_io_flag', False), # when true: no rank ever holds the full (nd, Nens) ensemble
        'adaptive_radius': _enkf_section.get('adaptive_radius', 1), # adaptive radius flag: True, False
        # Feature-flagged parameter-inference plugin
        'inference_plugin_enabled': bool(_enkf_section.get('inference_plugin_enabled', False)),

        # SMB: observed parameters use the EnKF posterior; unobserved SMB in
        # vec_inputs automatically uses the continuity-based inference branch.
        'physics_smb_inference': bool(_enkf_section.get('physics_smb_inference', False)),
        'smb_history_length': int(_enkf_section.get('smb_history_length', 5)),
        'smb_divergence_neighbors': int(_enkf_section.get('smb_divergence_neighbors', 24)),
        'smb_graph_neighbors': int(_enkf_section.get('smb_graph_neighbors', 12)),
        'smb_spatial_regularization': float(_enkf_section.get('smb_spatial_regularization', 25.0)),
        'smb_temporal_regularization': float(_enkf_section.get('smb_temporal_regularization', 4.0)),
        'smb_blend_factor': float(_enkf_section.get('smb_blend_factor', 0.35)),
        'smb_inference_start_time': float(_enkf_section.get('smb_inference_start_time', 0.0)),
        'smb_spinup_hold_factor': float(_enkf_section.get('smb_spinup_hold_factor', 0.0)),
        'smb_blend_ramp_time': float(_enkf_section.get('smb_blend_ramp_time', 0.0)),
        'smb_projection_basis': str(_enkf_section.get('smb_projection_basis', 'none')),
        'smb_physical_bounds': _enkf_section.get('smb_physical_bounds', None),
        'mesh_coordinate_scale_to_m': float(_enkf_section.get('mesh_coordinate_scale_to_m', 1.0)),

        # Bed: the EnKF supplies a raw increment and the plugin applies the
        # declared increment prior/constraints when this feature is enabled.
        'physics_bed_inference': bool(_enkf_section.get('physics_bed_inference', False)),
        'bed_update_mode': str(_enkf_section.get('bed_update_mode', 'legacy')),
        'bed_inference_start_time': float(_enkf_section.get('bed_inference_start_time', 0.0)),
        'bed_spinup_hold_factor': float(_enkf_section.get('bed_spinup_hold_factor', 1.0)),
        'bed_blend_ramp_time': float(_enkf_section.get('bed_blend_ramp_time', 0.0)),
        'bed_update_blend_factor': float(_enkf_section.get('bed_update_blend_factor', 0.15)),
        'bed_spatial_regularization': float(_enkf_section.get('bed_spatial_regularization', 40.0)),
        'bed_graph_neighbors': int(_enkf_section.get('bed_graph_neighbors', 12)),
        'bed_max_update_per_cycle': _enkf_section.get('bed_max_update_per_cycle', None),
        'bed_projection_basis': str(_enkf_section.get('bed_projection_basis', 'none')),
        'bed_physical_bounds': _enkf_section.get('bed_physical_bounds', None),
        'bed_enforce_below_surface': bool(_enkf_section.get('bed_enforce_below_surface', True)),
        'bed_min_surface_separation': float(_enkf_section.get('bed_min_surface_separation', 1.0)),
        'bed_update_mask': _enkf_section.get('bed_update_mask', None),

    })


    # Flatten all YAML sections into the one runtime context.  The section
    # dictionaries remain local to this loader and are not propagated as a
    # second configuration API.
    icesee_kwargs.update(_physical_section)
    icesee_kwargs.update(_modeling_section)
    icesee_kwargs.update(_enkf_section)

    # Re-apply normalized and command-line-resolved values after the raw YAML
    # update.  Without this step, e.g. ``--Nens=40`` is silently replaced by
    # the YAML value while constructing the canonical dictionary.
    icesee_kwargs.update({
        'Nens': Nens,
        'data_path': data_path,
        'model_nprocs': model_nprocs,
        'random_field_method': random_field_method,
    })

    joint_estimated_params = len(icesee_kwargs.get('joint_estimated_params', []))
    if icesee_kwargs['joint_estimation']:
        icesee_kwargs['total_state_param_vars'] = icesee_kwargs['num_state_vars'] + icesee_kwargs['num_param_vars']
    else:
        icesee_kwargs['total_state_param_vars'] = icesee_kwargs['num_state_vars'] + icesee_kwargs['num_param_vars'] - joint_estimated_params

    # unpack standard deviations
    icesee_kwargs.update({
        'sig_model': _enkf_section.get('sig_model', np.array([0.01])*icesee_kwargs['total_state_param_vars']),
        'sig_obs': _enkf_section.get('sig_obs', np.array([0.01])*icesee_kwargs['total_state_param_vars']),
        'sig_Q': _enkf_section.get('sig_Q', np.array([0.01])*icesee_kwargs['total_state_param_vars']),
        })

    if not icesee_kwargs['joint_estimation']:
        icesee_kwargs.update({
            'sig_obs': np.array(icesee_kwargs['sig_obs'][:icesee_kwargs['num_state_vars']]),
            'sig_Q': np.array(icesee_kwargs['sig_Q'][:icesee_kwargs['num_state_vars']]),
            'sig_model': np.array(icesee_kwargs['sig_model'][:icesee_kwargs['num_state_vars']]),
        })
        icesee_kwargs['vec_inputs'] = icesee_kwargs['vec_inputs'][:icesee_kwargs['num_state_vars']]

    # --- Observations Parameters ---
    if icesee_kwargs.get('observations_available', False):
        # load observation data
        if not os.path.exists(icesee_kwargs.get('obs_data_path', 'observations_data.h5')):
            raise FileNotFoundError(f"Observation data file '{icesee_kwargs.get('obs_data_path', 'observations_data.h5')}' not found. Please ensure the file exists.")
        # Tell the user to load: icesee_kwargs['obs_index'] and icesee_kwargs['number_obs_instants']
        print("[ICESEE] Please load 'obs_index' and 'number_obs_instants' from the observation data file into the model dictionary.")
        # obs_data = extract_datasets_from_h5(icesee_kwargs.get('obs_data_path', 'observations_data.h5'))
        # icesee_kwargs.update({'obs_data': obs_data})
    else:
        # generate observation schedule for synthetic observations
        obs_t, obs_idx, num_observations = UtilsFunctions(icesee_kwargs).generate_observation_schedule(**icesee_kwargs)
        icesee_kwargs['obs_index'] = obs_idx
        icesee_kwargs['number_obs_instants'] = num_observations
        icesee_kwargs['m_obs'] = num_observations

    icesee_kwargs['commandlinerun']      = _enkf_section.get('commandlinerun', False)

    #  check available parameters in the obseve_params list that need to be observed
    params_vec = []
    for i, vars in enumerate(icesee_kwargs['vec_inputs']):
        if i >= icesee_kwargs['num_state_vars']:
            params_vec.append(vars)

    icesee_kwargs['params_vec'] = params_vec

    if int(icesee_kwargs['execution_mode']) == 0:
        icesee_rank = 0
        icesee_size = 1
        icesee_comm = None

        # check if _modelrun_datasets exists in path if not create one
        _modelrun_datasets = icesee_kwargs.get('data_path',None)
        if not os.path.exists(_modelrun_datasets):
            os.makedirs(_modelrun_datasets, exist_ok=True)
