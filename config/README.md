# ICESEE FLAGS README

## All Main Flags used in ICESEE 
| Name | Description | Type | Default | Required | Choices | Source |
|------|-------------|------|---------|----------|---------|--------|
| `--Nens` | ensemble members | int | 1 | No | None | CLI |
| `--data_path` | folder to save data for single or multiple runs | str | _modelrun_datasets | No | None | CLI |
| `--default_run` | default run | str | None | No | None | CLI |
| `--even_distribution` | even distribution | str | None | No | None | CLI |
| `--model_nprocs` | number of processors for the coupled model | int | 0 | No | None | CLI |
| `--sequential_run` | sequential run | str | None | No | None | CLI |
| `--verbose` | verbose output | str | None | No | None | CLI |
| `-F` | Path to YAML parameter file (default: params.yaml) | str | params.yaml | No | None | CLI |
| `Nens` | Parameter for nens in dictionary | Unknown | Computed | No | None | Dictionary |
| `Q_rho` | YAML configuration parameter for q rho | float | 1.0 | No | None | YAML |
| `batch_size` | YAML configuration parameter for batch size | int | 50 | No | None | YAML |
| `chunk_size` | YAML configuration parameter for chunk size | int | 5000 | No | None | YAML |
| `commandlinerun` | Parameter for commandlinerun in dictionary | Unknown | Computed | No | None | Dictionary |
| `coupled_model_datasets` | YAML configuration parameter for coupled model datasets | str | data | No | None | YAML |
| `coupled_model_datasets_dir` | YAML configuration parameter for coupled model datasets dir | str | data | No | None | YAML |
| `data_path` | Parameter for data path in dictionary | Unknown | Computed | No | None | Dictionary |
| `default_run` | Parameter for default run in dictionary | bool | True | No | None | Dictionary |
| `example_name` | YAML configuration parameter for example name | Unknown | None | No | None | YAML |
| `execution_flag` | Controls execution flag behavior in script logic | int | 0 | No | None | Dictionary |
| `execution_mode` | YAML configuration parameter for execution mode | int | 0 | No | None | YAML |
| `flag_jupyter` | leave entire routine | bool | True | No | None | Internal |
| `freq_obs` | YAML configuration parameter for freq obs | int | 1 | No | None | YAML |
| `generate_nurged_state` | YAML configuration parameter for generate nurged state | bool | True | No | None | YAML |
| `generate_synthetic_obs` | YAML configuration parameter for generate synthetic obs | bool | True | No | None | YAML |
| `generate_true_state` | YAML configuration parameter for generate true state | bool | True | No | None | YAML |
| `global_analysis` | YAML configuration parameter for global analysis | bool | True | No | None | YAML |
| `inflation_factor` | YAML configuration parameter for inflation factor | float | 1.0 | No | None | YAML |
| `joint_estimated_params` | Variable used for joint estimated params in script logic | Unknown | Unknown | No | None | Variable |
| `joint_estimation` | add joint estimation flag to params | Unknown | Unknown | No | None | Dictionary |
| `length_scale` | YAML configuration parameter for length scale | list | [] | No | None | YAML |
| `local_analysis` | YAML configuration parameter for local analysis | bool | False | No | None | YAML |
| `localization_flag` | Controls localization flag behavior in script logic | bool | False | No | None | YAML |
| `model_name` | YAML configuration parameter for model name | Unknown | None | No | None | YAML |
| `model_nprocs` | Parameter for model nprocs in dictionary | Unknown | Computed | No | None | Dictionary |
| `n_modeltasks` | YAML configuration parameter for n modeltasks | int | 1 | No | None | YAML |
| `num_param_vars` | YAML configuration parameter for num param vars | int | 0 | No | None | YAML |
| `num_state_vars` | YAML configuration parameter for num state vars | int | 1 | No | None | YAML |
| `number_obs_instants` | Parameter for number obs instants in dictionary | Unknown | Unknown | No | None | Dictionary |
| `obs_data_path` | YAML configuration parameter for obs data path | str | observations_data.h5 | No | None | YAML |
| `obs_index` | Parameter for obs index in dictionary | Unknown | Unknown | No | None | Dictionary |
| `obs_max_time` | YAML configuration parameter for obs max time | int | 1 | No | None | YAML |
| `obs_start_time` | YAML configuration parameter for obs start time | int | 1 | No | None | YAML |
| `observations_available` | --- Observations Parameters --- | bool | False | No | None | YAML |
| `observed_params` | YAML configuration parameter for observed params | list | [] | No | None | YAML |
| `parallel_flag` | Controls parallel flag behavior in script logic | Unknown | Computed | No | None | Dictionary |
| `param_ens_spread` | YAML configuration parameter for param ens spread | list | [] | No | None | YAML |
| `parameter_estimation` | YAML configuration parameter for parameter estimation | bool | False | No | None | YAML |
| `params_vec` | check available parameters in the obseve_params list that need to be observed | list | [] | No | None | Variable |
| `partial` | Parameter for partial in dictionary | bool | True | No | None | Dictionary |
| `run_flag` | Controls run flag behavior in script logic | bool | True | No | None | Internal |
| `sequential_ensemble_initialization` | YAML configuration parameter for sequential ensemble initialization | bool | False | No | None | YAML |
| `serial_file_creation` | YAML configuration parameter for serial file creation | bool | True | No | None | YAML |
| `sig_Q` | YAML configuration parameter for sig q | Unknown | None | No | None | YAML |
| `sig_model` | YAML configuration parameter for sig model | Unknown | None | No | None | YAML |
| `sig_obs` | YAML configuration parameter for sig obs | Unknown | None | No | None | YAML |
| `state_estimation` | YAML configuration parameter for state estimation | bool | False | No | None | YAML |
| `t` | update for time t | Unknown | Computed | No | None | Dictionary |
| `total_state_param_vars` | Parameter for total state param vars in dictionary | Unknown | Unknown | No | None | Dictionary |
| `use_ensemble_pertubations` | YAML configuration parameter for use ensemble pertubations | bool | True | No | None | YAML |
| `use_random_fields` | YAML configuration parameter for use random fields | bool | False | No | None | YAML |
| `vec_inputs` | Parameter for vec inputs in dictionary | Unknown | Unknown | No | None | Dictionary |
| `verbose` | YAML configuration parameter for verbose | bool | False | No | None | YAML |

## All Main Flags used in ICESEE 
| Name | Description | Type | Default | Required | Choices | Source |
|------|-------------|------|---------|----------|---------|--------|
| `--Nens` | ensemble members | int | 1 | No | None | CLI |
| `--data_path` | folder to save data for single or multiple runs | str | _modelrun_datasets | No | None | CLI |
| `--default_run` | default run | str | None | No | None | CLI |
| `--even_distribution` | even distribution | str | None | No | None | CLI |
| `--model_nprocs` | number of processors for the coupled model | int | 0 | No | None | CLI |
| `--sequential_run` | sequential run | str | None | No | None | CLI |
| `--verbose` | verbose output | str | None | No | None | CLI |
| `-F` | Path to YAML parameter file (default: params.yaml) | str | params.yaml | No | None | CLI |
| `Nens` | Parameter for nens in dictionary | Unknown | Computed | No | None | Dictionary |
| `Q_rho` | YAML configuration parameter for q rho | float | 1.0 | No | None | YAML |
| `batch_size` | YAML configuration parameter for batch size | int | 50 | No | None | YAML |
| `chunk_size` | YAML configuration parameter for chunk size | int | 5000 | No | None | YAML |
| `commandlinerun` | Parameter for commandlinerun in dictionary | Unknown | Computed | No | None | Dictionary |
| `coupled_model_datasets` | YAML configuration parameter for coupled model datasets | str | data | No | None | YAML |
| `coupled_model_datasets_dir` | YAML configuration parameter for coupled model datasets dir | str | data | No | None | YAML |
| `data_path` | Parameter for data path in dictionary | Unknown | Computed | No | None | Dictionary |
| `default_run` | Parameter for default run in dictionary | bool | True | No | None | Dictionary |
| `example_name` | YAML configuration parameter for example name | Unknown | None | No | None | YAML |
| `execution_flag` | Controls execution flag behavior in script logic | int | 0 | No | None | Dictionary |
| `execution_mode` | YAML configuration parameter for execution mode | int | 0 | No | None | YAML |
| `flag_jupyter` | leave entire routine | bool | True | No | None | Internal |
| `freq_obs` | YAML configuration parameter for freq obs | int | 1 | No | None | YAML |
| `generate_nurged_state` | YAML configuration parameter for generate nurged state | bool | True | No | None | YAML |
| `generate_synthetic_obs` | YAML configuration parameter for generate synthetic obs | bool | True | No | None | YAML |
| `generate_true_state` | YAML configuration parameter for generate true state | bool | True | No | None | YAML |
| `global_analysis` | YAML configuration parameter for global analysis | bool | True | No | None | YAML |
| `inflation_factor` | YAML configuration parameter for inflation factor | float | 1.0 | No | None | YAML |
| `joint_estimated_params` | Variable used for joint estimated params in script logic | Unknown | Unknown | No | None | Variable |
| `joint_estimation` | add joint estimation flag to params | Unknown | Unknown | No | None | Dictionary |
| `length_scale` | YAML configuration parameter for length scale | list | [] | No | None | YAML |
| `local_analysis` | YAML configuration parameter for local analysis | bool | False | No | None | YAML |
| `localization_flag` | Controls localization flag behavior in script logic | bool | False | No | None | YAML |
| `model_name` | YAML configuration parameter for model name | Unknown | None | No | None | YAML |
| `model_nprocs` | Parameter for model nprocs in dictionary | Unknown | Computed | No | None | Dictionary |
| `n_modeltasks` | YAML configuration parameter for n modeltasks | int | 1 | No | None | YAML |
| `num_param_vars` | YAML configuration parameter for num param vars | int | 0 | No | None | YAML |
| `num_state_vars` | YAML configuration parameter for num state vars | int | 1 | No | None | YAML |
| `number_obs_instants` | Parameter for number obs instants in dictionary | Unknown | Unknown | No | None | Dictionary |
| `obs_data_path` | YAML configuration parameter for obs data path | str | observations_data.h5 | No | None | YAML |
| `obs_index` | Parameter for obs index in dictionary | Unknown | Unknown | No | None | Dictionary |
| `obs_max_time` | YAML configuration parameter for obs max time | int | 1 | No | None | YAML |
| `obs_start_time` | YAML configuration parameter for obs start time | int | 1 | No | None | YAML |
| `observations_available` | --- Observations Parameters --- | bool | False | No | None | YAML |
| `observed_params` | YAML configuration parameter for observed params | list | [] | No | None | YAML |
| `parallel_flag` | Controls parallel flag behavior in script logic | Unknown | Computed | No | None | Dictionary |
| `param_ens_spread` | YAML configuration parameter for param ens spread | list | [] | No | None | YAML |
| `parameter_estimation` | YAML configuration parameter for parameter estimation | bool | False | No | None | YAML |
| `params_vec` | check available parameters in the obseve_params list that need to be observed | list | [] | No | None | Variable |
| `partial` | Parameter for partial in dictionary | bool | True | No | None | Dictionary |
| `run_flag` | Controls run flag behavior in script logic | bool | True | No | None | Internal |
| `sequential_ensemble_initialization` | YAML configuration parameter for sequential ensemble initialization | bool | False | No | None | YAML |
| `serial_file_creation` | YAML configuration parameter for serial file creation | bool | True | No | None | YAML |
| `sig_Q` | YAML configuration parameter for sig q | Unknown | None | No | None | YAML |
| `sig_model` | YAML configuration parameter for sig model | Unknown | None | No | None | YAML |
| `sig_obs` | YAML configuration parameter for sig obs | Unknown | None | No | None | YAML |
| `state_estimation` | YAML configuration parameter for state estimation | bool | False | No | None | YAML |
| `t` | update for time t | Unknown | Computed | No | None | Dictionary |
| `total_state_param_vars` | Parameter for total state param vars in dictionary | Unknown | Unknown | No | None | Dictionary |
| `use_ensemble_pertubations` | YAML configuration parameter for use ensemble pertubations | bool | True | No | None | YAML |
| `use_random_fields` | YAML configuration parameter for use random fields | bool | False | No | None | YAML |
| `vec_inputs` | Parameter for vec inputs in dictionary | Unknown | Unknown | No | None | Dictionary |
| `verbose` | YAML configuration parameter for verbose | bool | False | No | None | YAML |

## All Main Flags used in ICESEE 
| Name | Description | Type | Default | Required | Choices | Source |
|------|-------------|------|---------|----------|---------|--------|
| `--Nens` | ensemble members | int | 1 | No | None | CLI |
| `--data_path` | folder to save data for single or multiple runs | str | _modelrun_datasets | No | None | CLI |
| `--default_run` | default run | str | None | No | None | CLI |
| `--even_distribution` | even distribution | str | None | No | None | CLI |
| `--model_nprocs` | number of processors for the coupled model | int | 0 | No | None | CLI |
| `--sequential_run` | sequential run | str | None | No | None | CLI |
| `--verbose` | verbose output | str | None | No | None | CLI |
| `-F` | Path to YAML parameter file (default: params.yaml) | str | params.yaml | No | None | CLI |
| `Nens` | Parameter for nens in dictionary | Unknown | Computed | No | None | Dictionary |
| `Q_rho` | YAML configuration parameter for q rho | float | 1.0 | No | None | YAML |
| `batch_size` | YAML configuration parameter for batch size | int | 50 | No | None | YAML |
| `chunk_size` | YAML configuration parameter for chunk size | int | 5000 | No | None | YAML |
| `commandlinerun` | Parameter for commandlinerun in dictionary | Unknown | Computed | No | None | Dictionary |
| `coupled_model_datasets` | YAML configuration parameter for coupled model datasets | str | data | No | None | YAML |
| `coupled_model_datasets_dir` | YAML configuration parameter for coupled model datasets dir | str | data | No | None | YAML |
| `data_path` | Parameter for data path in dictionary | Unknown | Computed | No | None | Dictionary |
| `default_run` | Parameter for default run in dictionary | bool | True | No | None | Dictionary |
| `example_name` | YAML configuration parameter for example name | Unknown | None | No | None | YAML |
| `execution_flag` | Controls execution flag behavior in script logic | int | 0 | No | None | Dictionary |
| `execution_mode` | YAML configuration parameter for execution mode | int | 0 | No | None | YAML |
| `flag_jupyter` | leave entire routine | bool | True | No | None | Internal |
| `freq_obs` | YAML configuration parameter for freq obs | int | 1 | No | None | YAML |
| `generate_nurged_state` | YAML configuration parameter for generate nurged state | bool | True | No | None | YAML |
| `generate_synthetic_obs` | YAML configuration parameter for generate synthetic obs | bool | True | No | None | YAML |
| `generate_true_state` | YAML configuration parameter for generate true state | bool | True | No | None | YAML |
| `global_analysis` | YAML configuration parameter for global analysis | bool | True | No | None | YAML |
| `inflation_factor` | YAML configuration parameter for inflation factor | float | 1.0 | No | None | YAML |
| `joint_estimated_params` | Variable used for joint estimated params in script logic | Unknown | Unknown | No | None | Variable |
| `joint_estimation` | add joint estimation flag to params | Unknown | Unknown | No | None | Dictionary |
| `length_scale` | YAML configuration parameter for length scale | list | [] | No | None | YAML |
| `local_analysis` | YAML configuration parameter for local analysis | bool | False | No | None | YAML |
| `localization_flag` | Controls localization flag behavior in script logic | bool | False | No | None | YAML |
| `model_name` | YAML configuration parameter for model name | Unknown | None | No | None | YAML |
| `model_nprocs` | Parameter for model nprocs in dictionary | Unknown | Computed | No | None | Dictionary |
| `n_modeltasks` | YAML configuration parameter for n modeltasks | int | 1 | No | None | YAML |
| `num_param_vars` | YAML configuration parameter for num param vars | int | 0 | No | None | YAML |
| `num_state_vars` | YAML configuration parameter for num state vars | int | 1 | No | None | YAML |
| `number_obs_instants` | Parameter for number obs instants in dictionary | Unknown | Unknown | No | None | Dictionary |
| `obs_data_path` | YAML configuration parameter for obs data path | str | observations_data.h5 | No | None | YAML |
| `obs_index` | Parameter for obs index in dictionary | Unknown | Unknown | No | None | Dictionary |
| `obs_max_time` | YAML configuration parameter for obs max time | int | 1 | No | None | YAML |
| `obs_start_time` | YAML configuration parameter for obs start time | int | 1 | No | None | YAML |
| `observations_available` | --- Observations Parameters --- | bool | False | No | None | YAML |
| `observed_params` | YAML configuration parameter for observed params | list | [] | No | None | YAML |
| `parallel_flag` | Controls parallel flag behavior in script logic | Unknown | Computed | No | None | Dictionary |
| `param_ens_spread` | YAML configuration parameter for param ens spread | list | [] | No | None | YAML |
| `parameter_estimation` | YAML configuration parameter for parameter estimation | bool | False | No | None | YAML |
| `params_vec` | check available parameters in the obseve_params list that need to be observed | list | [] | No | None | Variable |
| `partial` | Parameter for partial in dictionary | bool | True | No | None | Dictionary |
| `run_flag` | Controls run flag behavior in script logic | bool | True | No | None | Internal |
| `sequential_ensemble_initialization` | YAML configuration parameter for sequential ensemble initialization | bool | False | No | None | YAML |
| `serial_file_creation` | YAML configuration parameter for serial file creation | bool | True | No | None | YAML |
| `sig_Q` | YAML configuration parameter for sig q | Unknown | None | No | None | YAML |
| `sig_model` | YAML configuration parameter for sig model | Unknown | None | No | None | YAML |
| `sig_obs` | YAML configuration parameter for sig obs | Unknown | None | No | None | YAML |
| `state_estimation` | YAML configuration parameter for state estimation | bool | False | No | None | YAML |
| `t` | update for time t | Unknown | Computed | No | None | Dictionary |
| `total_state_param_vars` | Parameter for total state param vars in dictionary | Unknown | Unknown | No | None | Dictionary |
| `use_ensemble_pertubations` | YAML configuration parameter for use ensemble pertubations | bool | True | No | None | YAML |
| `use_random_fields` | YAML configuration parameter for use random fields | bool | False | No | None | YAML |
| `vec_inputs` | Parameter for vec inputs in dictionary | Unknown | Unknown | No | None | Dictionary |
| `verbose` | YAML configuration parameter for verbose | bool | False | No | None | YAML |

## All Main Flags used in ICESEE 
| Name | Description | Type | Default | Required | Choices | Source |
|------|-------------|------|---------|----------|---------|--------|
| `--Nens` | ensemble members | int | 1 | No | None | CLI |
| `--data_path` | folder to save data for single or multiple runs | str | _modelrun_datasets | No | None | CLI |
| `--default_run` | default run | str | None | No | None | CLI |
| `--even_distribution` | even distribution | str | None | No | None | CLI |
| `--model_nprocs` | number of processors for the coupled model | int | 0 | No | None | CLI |
| `--sequential_run` | sequential run | str | None | No | None | CLI |
| `--verbose` | verbose output | str | None | No | None | CLI |
| `-F` | Path to YAML parameter file (default: params.yaml) | str | params.yaml | No | None | CLI |
| `Nens` | Parameter for nens in dictionary | Unknown | Computed | No | None | Dictionary |
| `Q_rho` | YAML configuration parameter for q rho | float | 1.0 | No | None | YAML |
| `batch_size` | YAML configuration parameter for batch size | int | 50 | No | None | YAML |
| `chunk_size` | YAML configuration parameter for chunk size | int | 5000 | No | None | YAML |
| `commandlinerun` | Parameter for commandlinerun in dictionary | Unknown | Computed | No | None | Dictionary |
| `coupled_model_datasets` | YAML configuration parameter for coupled model datasets | str | data | No | None | YAML |
| `coupled_model_datasets_dir` | YAML configuration parameter for coupled model datasets dir | str | data | No | None | YAML |
| `data_path` | Parameter for data path in dictionary | Unknown | Computed | No | None | Dictionary |
| `default_run` | Parameter for default run in dictionary | bool | True | No | None | Dictionary |
| `example_name` | YAML configuration parameter for example name | Unknown | None | No | None | YAML |
| `execution_flag` | Controls execution flag behavior in script logic | int | 0 | No | None | Dictionary |
| `execution_mode` | YAML configuration parameter for execution mode | int | 0 | No | None | YAML |
| `flag_jupyter` | leave entire routine | bool | True | No | None | Internal |
| `freq_obs` | YAML configuration parameter for freq obs | int | 1 | No | None | YAML |
| `generate_nurged_state` | YAML configuration parameter for generate nurged state | bool | True | No | None | YAML |
| `generate_synthetic_obs` | YAML configuration parameter for generate synthetic obs | bool | True | No | None | YAML |
| `generate_true_state` | YAML configuration parameter for generate true state | bool | True | No | None | YAML |
| `global_analysis` | YAML configuration parameter for global analysis | bool | True | No | None | YAML |
| `inflation_factor` | YAML configuration parameter for inflation factor | float | 1.0 | No | None | YAML |
| `joint_estimated_params` | Variable used for joint estimated params in script logic | Unknown | Unknown | No | None | Variable |
| `joint_estimation` | add joint estimation flag to params | Unknown | Unknown | No | None | Dictionary |
| `length_scale` | YAML configuration parameter for length scale | list | [] | No | None | YAML |
| `local_analysis` | YAML configuration parameter for local analysis | bool | False | No | None | YAML |
| `localization_flag` | Controls localization flag behavior in script logic | bool | False | No | None | YAML |
| `model_name` | YAML configuration parameter for model name | Unknown | None | No | None | YAML |
| `model_nprocs` | Parameter for model nprocs in dictionary | Unknown | Computed | No | None | Dictionary |
| `n_modeltasks` | YAML configuration parameter for n modeltasks | int | 1 | No | None | YAML |
| `num_param_vars` | YAML configuration parameter for num param vars | int | 0 | No | None | YAML |
| `num_state_vars` | YAML configuration parameter for num state vars | int | 1 | No | None | YAML |
| `number_obs_instants` | Parameter for number obs instants in dictionary | Unknown | Unknown | No | None | Dictionary |
| `obs_data_path` | YAML configuration parameter for obs data path | str | observations_data.h5 | No | None | YAML |
| `obs_index` | Parameter for obs index in dictionary | Unknown | Unknown | No | None | Dictionary |
| `obs_max_time` | YAML configuration parameter for obs max time | int | 1 | No | None | YAML |
| `obs_start_time` | YAML configuration parameter for obs start time | int | 1 | No | None | YAML |
| `observations_available` | --- Observations Parameters --- | bool | False | No | None | YAML |
| `observed_params` | YAML configuration parameter for observed params | list | [] | No | None | YAML |
| `parallel_flag` | Controls parallel flag behavior in script logic | Unknown | Computed | No | None | Dictionary |
| `param_ens_spread` | YAML configuration parameter for param ens spread | list | [] | No | None | YAML |
| `parameter_estimation` | YAML configuration parameter for parameter estimation | bool | False | No | None | YAML |
| `params_vec` | check available parameters in the obseve_params list that need to be observed | list | [] | No | None | Variable |
| `partial` | Parameter for partial in dictionary | bool | True | No | None | Dictionary |
| `run_flag` | Controls run flag behavior in script logic | bool | True | No | None | Internal |
| `sequential_ensemble_initialization` | YAML configuration parameter for sequential ensemble initialization | bool | False | No | None | YAML |
| `serial_file_creation` | YAML configuration parameter for serial file creation | bool | True | No | None | YAML |
| `sig_Q` | YAML configuration parameter for sig q | Unknown | None | No | None | YAML |
| `sig_model` | YAML configuration parameter for sig model | Unknown | None | No | None | YAML |
| `sig_obs` | YAML configuration parameter for sig obs | Unknown | None | No | None | YAML |
| `state_estimation` | YAML configuration parameter for state estimation | bool | False | No | None | YAML |
| `t` | update for time t | Unknown | Computed | No | None | Dictionary |
| `total_state_param_vars` | Parameter for total state param vars in dictionary | Unknown | Unknown | No | None | Dictionary |
| `use_ensemble_pertubations` | YAML configuration parameter for use ensemble pertubations | bool | True | No | None | YAML |
| `use_random_fields` | YAML configuration parameter for use random fields | bool | False | No | None | YAML |
| `vec_inputs` | Parameter for vec inputs in dictionary | Unknown | Unknown | No | None | Dictionary |
| `verbose` | YAML configuration parameter for verbose | bool | False | No | None | YAML |

## All Main Flags used in ICESEE 
| Name | Description | Type | Default | Required | Choices | Source |
|------|-------------|------|---------|----------|---------|--------|
| `--Nens` | ensemble members | int | 1 | No | None | CLI |
| `--data_path` | folder to save data for single or multiple runs | str | _modelrun_datasets | No | None | CLI |
| `--default_run` | default run | str | None | No | None | CLI |
| `--even_distribution` | even distribution | str | None | No | None | CLI |
| `--model_nprocs` | number of processors for the coupled model | int | 0 | No | None | CLI |
| `--sequential_run` | sequential run | str | None | No | None | CLI |
| `--verbose` | verbose output | str | None | No | None | CLI |
| `-F` | Path to YAML parameter file (default: params.yaml) | str | params.yaml | No | None | CLI |
| `Nens` | Parameter for nens in dictionary | Unknown | Computed | No | None | Dictionary |
| `Q_rho` | YAML configuration parameter for q rho | float | 1.0 | No | None | YAML |
| `batch_size` | YAML configuration parameter for batch size | int | 50 | No | None | YAML |
| `chunk_size` | YAML configuration parameter for chunk size | int | 5000 | No | None | YAML |
| `commandlinerun` | Parameter for commandlinerun in dictionary | Unknown | Computed | No | None | Dictionary |
| `coupled_model_datasets` | YAML configuration parameter for coupled model datasets | str | data | No | None | YAML |
| `coupled_model_datasets_dir` | YAML configuration parameter for coupled model datasets dir | str | data | No | None | YAML |
| `data_path` | Parameter for data path in dictionary | Unknown | Computed | No | None | Dictionary |
| `default_run` | Parameter for default run in dictionary | bool | True | No | None | Dictionary |
| `example_name` | YAML configuration parameter for example name | Unknown | None | No | None | YAML |
| `execution_flag` | Controls execution flag behavior in script logic | int | 0 | No | None | Dictionary |
| `execution_mode` | YAML configuration parameter for execution mode | int | 0 | No | None | YAML |
| `flag_jupyter` | leave entire routine | bool | True | No | None | Internal |
| `freq_obs` | YAML configuration parameter for freq obs | int | 1 | No | None | YAML |
| `generate_nurged_state` | YAML configuration parameter for generate nurged state | bool | True | No | None | YAML |
| `generate_synthetic_obs` | YAML configuration parameter for generate synthetic obs | bool | True | No | None | YAML |
| `generate_true_state` | YAML configuration parameter for generate true state | bool | True | No | None | YAML |
| `global_analysis` | YAML configuration parameter for global analysis | bool | True | No | None | YAML |
| `inflation_factor` | YAML configuration parameter for inflation factor | float | 1.0 | No | None | YAML |
| `joint_estimated_params` | Variable used for joint estimated params in script logic | Unknown | Unknown | No | None | Variable |
| `joint_estimation` | add joint estimation flag to params | Unknown | Unknown | No | None | Dictionary |
| `length_scale` | YAML configuration parameter for length scale | list | [] | No | None | YAML |
| `local_analysis` | YAML configuration parameter for local analysis | bool | False | No | None | YAML |
| `localization_flag` | Controls localization flag behavior in script logic | bool | False | No | None | YAML |
| `model_name` | YAML configuration parameter for model name | Unknown | None | No | None | YAML |
| `model_nprocs` | Parameter for model nprocs in dictionary | Unknown | Computed | No | None | Dictionary |
| `n_modeltasks` | YAML configuration parameter for n modeltasks | int | 1 | No | None | YAML |
| `num_param_vars` | YAML configuration parameter for num param vars | int | 0 | No | None | YAML |
| `num_state_vars` | YAML configuration parameter for num state vars | int | 1 | No | None | YAML |
| `number_obs_instants` | Parameter for number obs instants in dictionary | Unknown | Unknown | No | None | Dictionary |
| `obs_data_path` | YAML configuration parameter for obs data path | str | observations_data.h5 | No | None | YAML |
| `obs_index` | Parameter for obs index in dictionary | Unknown | Unknown | No | None | Dictionary |
| `obs_max_time` | YAML configuration parameter for obs max time | int | 1 | No | None | YAML |
| `obs_start_time` | YAML configuration parameter for obs start time | int | 1 | No | None | YAML |
| `observations_available` | --- Observations Parameters --- | bool | False | No | None | YAML |
| `observed_params` | YAML configuration parameter for observed params | list | [] | No | None | YAML |
| `parallel_flag` | Controls parallel flag behavior in script logic | Unknown | Computed | No | None | Dictionary |
| `param_ens_spread` | YAML configuration parameter for param ens spread | list | [] | No | None | YAML |
| `parameter_estimation` | YAML configuration parameter for parameter estimation | bool | False | No | None | YAML |
| `params_vec` | check available parameters in the obseve_params list that need to be observed | list | [] | No | None | Variable |
| `partial` | Parameter for partial in dictionary | bool | True | No | None | Dictionary |
| `run_flag` | Controls run flag behavior in script logic | bool | True | No | None | Internal |
| `sequential_ensemble_initialization` | YAML configuration parameter for sequential ensemble initialization | bool | False | No | None | YAML |
| `serial_file_creation` | YAML configuration parameter for serial file creation | bool | True | No | None | YAML |
| `sig_Q` | YAML configuration parameter for sig q | Unknown | None | No | None | YAML |
| `sig_model` | YAML configuration parameter for sig model | Unknown | None | No | None | YAML |
| `sig_obs` | YAML configuration parameter for sig obs | Unknown | None | No | None | YAML |
| `state_estimation` | YAML configuration parameter for state estimation | bool | False | No | None | YAML |
| `t` | update for time t | Unknown | Computed | No | None | Dictionary |
| `total_state_param_vars` | Parameter for total state param vars in dictionary | Unknown | Unknown | No | None | Dictionary |
| `use_ensemble_pertubations` | YAML configuration parameter for use ensemble pertubations | bool | True | No | None | YAML |
| `use_random_fields` | YAML configuration parameter for use random fields | bool | False | No | None | YAML |
| `vec_inputs` | Parameter for vec inputs in dictionary | Unknown | Unknown | No | None | Dictionary |
| `verbose` | YAML configuration parameter for verbose | bool | False | No | None | YAML |
