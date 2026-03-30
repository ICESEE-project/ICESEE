<!-- BEGIN: ICESEE-FLAGS -->
## All Main Flags Used in ICESEE

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
| `ICESEE_PERFORMANCE_TEST` | this is an environment variable | bool | False | No | None | YAML |
| `Nens` | Parameter for nens in dictionary | Unknown | Computed | No | None | Dictionary |
| `Q_rho` | YAML configuration parameter for q rho | float | 1.0 | No | None | YAML |
| `abs_vel_weight` | weight for absolute velocity in inversion | float | 1.0 | No | None | YAML |
| `base_seed` | YAML configuration parameter for base seed | int | 42 | No | None | YAML |
| `batch_size` | number of time steps to process in each batch | Unknown | Computed | No | None | Dictionary |
| `bed_obs_indices` | specific indices to observe {list} (bed subvector indices) | NoneType | None | No | None | YAML |
| `bed_obs_mask` | boolean mask array for bed observations {np.array} | NoneType | None | No | None | YAML |
| `bed_obs_snapshot` | list of time snapshots to observe bed variables | list | [] | No | None | YAML |
| `bed_obs_spacing` | observation spacing every n grid points {int} | NoneType | None | No | None | YAML |
| `bed_obs_stride` | spatial stride in km for bed observations | NoneType | None | No | None | YAML |
| `bed_relaxation_factor` | relaxation factor for bed elevation updates (-1 < factor <= 1) (when bed is not observed) | float | 0.05 | No | None | YAML |
| `checkpoint_every` | YAML configuration parameter for checkpoint every | int | 1 | No | None | YAML |
| `chunk_size` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `collective_threshold` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `commandlinerun` | Parameter for commandlinerun in dictionary | Unknown | Computed | No | None | Dictionary |
| `coupled_model_datasets` | YAML configuration parameter for coupled model datasets | str | data | No | None | YAML |
| `coupled_model_datasets_dir` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `create_ensemble_dataset` | YAML configuration parameter for create ensemble dataset | bool | True | No | None | YAML |
| `data_path` | Parameter for data path in dictionary | Unknown | Computed | No | None | Dictionary |
| `default_run` | Parameter for default run in dictionary | bool | True | No | None | Dictionary |
| `dt` | --- Ensemble Parameters --- | Unknown | Unknown | No | None | Dictionary |
| `enkf_params` | Parameter for enkf params in dictionary | Unknown | Computed(enkf_params) | No | None | Dictionary |
| `enkf_params_keys` | Parameter for enkf params keys in dictionary | dict | Unknown | No | None | Dictionary |
| `even_distribution` | Parameter for even distribution in dictionary | bool | True | No | None | Dictionary |
| `example_name` | YAML configuration parameter for example name | Unknown | Computed | No | None | YAML |
| `execution_flag` | Controls execution flag behavior in script logic | int | 0 | No | None | Dictionary |
| `execution_mode` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `flag_jupyter` | leave entire routine | bool | True | No | None | Internal |
| `force_fresh_start` | YAML configuration parameter for force fresh start | bool | False | No | None | YAML |
| `freq_obs` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `friction_idx` | YAML configuration parameter for friction idx | int | 5 | No | None | YAML |
| `generate_nurged_state` | YAML configuration parameter for generate nurged state | bool | True | No | None | YAML |
| `generate_synthetic_obs` | YAML configuration parameter for generate synthetic obs | bool | True | No | None | YAML |
| `generate_synthetic_obs_only` | flag to only generate synthetic observations without running the assimilation | bool | False | No | None | YAML |
| `generate_true_state` | YAML configuration parameter for generate true state | bool | True | No | None | YAML |
| `generate_true_wrong_state_only` | flag to only generate true and wrong state without running the assimilation | bool | False | No | None | YAML |
| `global_analysis` | YAML configuration parameter for global analysis | bool | True | No | None | YAML |
| `h5_file_chunk_size` | YAML configuration parameter for h5 file chunk size | int | 1000 | No | None | YAML |
| `h5_file_compression` | e.g., 'gzip' or 'lzf' or 'szip' or None | NoneType | None | No | None | YAML |
| `h5_file_compression_level` | 0-9 for gzip, 1-9 for szip, ignored for lzf and None | int | 4 | No | None | YAML |
| `inflation_factor` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `initial_bed_bias` | initial bias for bed elevation (in model units) | float | 0.0015 | No | None | YAML |
| `initial_spread_factor` | YAML configuration parameter for initial spread factor | float | 1.0 | No | None | YAML |
| `initialize_ensemble` | YAML configuration parameter for initialize ensemble | bool | True | No | None | YAML |
| `inversion_flag` | Controls inversion flag behavior in script logic | bool | False | No | None | YAML |
| `joint_estimated_params` | Variable used for joint estimated params in script logic | Unknown | Computed | No | None | Variable |
| `joint_estimation` | add joint estimation flag to params | Unknown | Unknown | No | None | Dictionary |
| `k_start_override` | YAML configuration parameter for k start override | NoneType | None | No | None | YAML |
| `length_scale` | YAML configuration parameter for length scale | list | [] | No | None | YAML |
| `local_analysis` | YAML configuration parameter for local analysis | bool | False | No | None | YAML |
| `localization_flag` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `m_obs` | Parameter for m obs in dictionary | Unknown | Computed(num_observations) | No | None | Dictionary |
| `mode` | Parameter for mode in dictionary | Unknown | Computed(execution_mode) | No | None | Dictionary |
| `model_name` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `model_nprocs` | Parameter for model nprocs in dictionary | Unknown | Computed | No | None | Dictionary |
| `modeling_params` | Parameter for modeling params in dictionary | Unknown | Computed(modeling_params) | No | None | Dictionary |
| `modeling_params_keys` | Parameter for modeling params keys in dictionary | dict | Unknown | No | None | Dictionary |
| `n_modeltasks` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `nt` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `num_param_vars` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `num_state_vars` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `number_obs_instants` | Parameter for number obs instants in dictionary | Unknown | Computed(num_observations) | No | None | Dictionary |
| `obs_data_path` | YAML configuration parameter for obs data path | str | observations_data.h5 | No | None | YAML |
| `obs_index` | Parameter for obs index in dictionary | Unknown | Computed(obs_idx) | No | None | Dictionary |
| `obs_max_time` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `obs_start_time` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `observations_available` | --- Observations Parameters --- | bool | False | No | None | YAML |
| `observed_params` | YAML configuration parameter for observed params | list | [] | No | None | YAML |
| `observed_vars` | YAML configuration parameter for observed vars | list | [] | No | None | YAML |
| `parallel_flag` | Controls parallel flag behavior in script logic | Unknown | Computed | No | None | Dictionary |
| `param_ens_spread` | YAML configuration parameter for param ens spread | list | [] | No | None | YAML |
| `parameter_estimation` | YAML configuration parameter for parameter estimation | bool | False | No | None | YAML |
| `params` | Parameter for params in dictionary | Unknown | Computed(params) | No | None | Dictionary |
| `params_vec` | check available parameters in the obseve_params list that need to be observed | list | [] | No | None | Variable |
| `partial` | Parameter for partial in dictionary | bool | True | No | None | Dictionary |
| `physical_params` | # update kwargs dictonary with params | Unknown | Computed(physical_params) | No | None | Dictionary |
| `physical_params_keys` | -- update the kwargs with physical, modeling and enkf parameters | dict | Unknown | No | None | Dictionary |
| `rel_vel_weight` | weight for relative velocity in inversion | float | 1.0 | No | None | YAML |
| `restart_enabled` | YAML configuration parameter for restart enabled | bool | True | No | None | YAML |
| `run_flag` | Controls run flag behavior in script logic | bool | True | No | None | Internal |
| `scalar_inputs` | list of scalar input variables | list | [] | No | None | YAML |
| `sequential_ensemble_initialization` | YAML configuration parameter for sequential ensemble initialization | bool | False | No | None | YAML |
| `sequential_run` | Parameter for sequential run in dictionary | bool | True | No | None | Dictionary |
| `serial_file_creation` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `sig_Q` | Parameter for sig q in dictionary | Unknown | Computed | No | None | Dictionary |
| `sig_model` | Parameter for sig model in dictionary | Unknown | Computed | No | None | Dictionary |
| `sig_obs` | Parameter for sig obs in dictionary | Unknown | Computed | No | None | Dictionary |
| `state_estimation` | YAML configuration parameter for state estimation | bool | False | No | None | YAML |
| `t` | update for time t | Unknown | Computed | No | None | Dictionary |
| `tikhonov_regularization_weight` | Tikhonov regularization weight for inversion | float | 1e-13 | No | None | YAML |
| `total_state_param_vars` | Parameter for total state param vars in dictionary | Unknown | Unknown | No | None | Dictionary |
| `use_ensemble_pertubations` | YAML configuration parameter for use ensemble pertubations | bool | True | No | None | YAML |
| `use_random_fields` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `var_nd` | variable state dimension for each state variable in vec_inputs. Used when state variables have different dimensions | NoneType | None | No | None | YAML |
| `vec_inputs` | Parameter for vec inputs in dictionary | Unknown | Unknown | No | None | Dictionary |
| `vel_idx` | YAML configuration parameter for vel idx | int | 2 | No | None | YAML |
| `verbose` | YAML configuration parameter for verbose | bool | False | No | None | YAML |
<!-- END: ICESEE-FLAGS -->
