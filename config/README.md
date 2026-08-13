<!-- BEGIN: ICESEE-FLAGS -->
## All Main Flags Used in ICESEE

| Name | Description | Type | Default | Required | Choices | Source |
|------|-------------|------|---------|----------|---------|--------|
| `--Nens` | ensemble members | int | None | No | None | CLI |
| `--data_path` | folder to save data for single or multiple runs | str | None | No | None | CLI |
| `--default_run` | default run | str | None | No | None | CLI |
| `--even_distribution` | even distribution | str | None | No | None | CLI |
| `--model_nprocs` | number of processors for the coupled model | int | None | No | None | CLI |
| `--sequential_run` | sequential run | str | None | No | None | CLI |
| `--verbose` | verbose output | str | None | No | None | CLI |
| `-F` | Path to YAML parameter file (default: params.yaml) | str | params.yaml | No | None | CLI |
| `ICESEE_PERFORMANCE_TEST` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `Nens` | Re-apply normalized and command-line-resolved values after the raw YAML update.  Without this step, e.g. ``--Nens=40`` is silently replaced by the YAML value while constructing the canonical dictionary. | Unknown | Computed(Nens) | No | None | Dictionary |
| `Q_rho` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `_enkf_section_keys` | Parameter for  enkf section keys in dictionary | dict | Unknown | No | None | Dictionary |
| `_modeling_section_keys` | Parameter for  modeling section keys in dictionary | dict | Unknown | No | None | Dictionary |
| `_physical_section_keys` | Flatten all YAML sections into the one runtime context.  The section dictionaries remain local to this loader and are not propagated as a second configuration API. | dict | Unknown | No | None | Dictionary |
| `abs_vel_weight` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `adaptive_radius` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `base_seed` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `batch_size` | number of time steps to process in each batch | Unknown | Computed | No | None | Dictionary |
| `bed_blend_ramp_time` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_enforce_below_surface` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_graph_neighbors` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_inference_start_time` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_inflation_factor` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `bed_max_update_per_cycle` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_min_surface_separation` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_obs_indices` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_obs_mask` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_obs_snapshot` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_obs_spacing` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_obs_stride` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_obs_track_half_width_m` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_physical_bounds` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_projection_basis` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_relaxation_factor` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_spatial_regularization` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_spinup_hold_factor` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_update_blend_factor` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_update_domain` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_update_mask` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `bed_update_mode` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `checkpoint_every` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `chunk_size` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `collective_threshold` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `commandlinerun` | Parameter for commandlinerun in dictionary | Unknown | Computed | No | None | Dictionary |
| `coupled_model_datasets` | YAML configuration parameter for coupled model datasets | str | data | No | None | YAML |
| `coupled_model_datasets_dir` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `create_ensemble_dataset` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `data_path` | Re-apply normalized and command-line-resolved values after the raw YAML update.  Without this step, e.g. ``--Nens=40`` is silently replaced by the YAML value while constructing the canonical dictionary. | Unknown | Computed(data_path) | No | None | Dictionary |
| `default_run` | Parameter for default run in dictionary | bool | True | No | None | Dictionary |
| `dt` | Add model and analysis options to the canonical runtime context. | Unknown | Unknown | No | None | Dictionary |
| `enkf_field_method` | YAML configuration parameter for enkf field method | str | fft | No | None | YAML |
| `enkf_observation_error_mode` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `even_distribution` | Parameter for even distribution in dictionary | bool | True | No | None | Dictionary |
| `example_name` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `execution_flag` | Controls execution flag behavior in script logic | int | 0 | No | None | Dictionary |
| `execution_mode` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `flag_jupyter` | leave entire routine | bool | True | No | None | Internal |
| `force_fresh_start` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `freq_obs` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `friction_idx` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `frozen_analysis_vars` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `generate_nurged_state` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `generate_synthetic_obs` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `generate_synthetic_obs_only` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `generate_true_state` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `generate_true_wrong_state_only` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `global_analysis` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `h5_file_chunk_size` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `h5_file_compression` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `h5_file_compression_level` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `inference_plugin_enabled` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `inflation_factor` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `initial_bed_bias` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `initial_spread_factor` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `initial_state_only` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `initialize_ensemble` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `inversion_flag` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `joint_estimated_params` | Variable used for joint estimated params in script logic | Unknown | Computed | No | None | Variable |
| `joint_estimation` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `k_start_override` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `length_scale` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `local_analysis` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `localization_flag` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `localization_radius` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `localized_vars` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `m_obs` | Parameter for m obs in dictionary | Unknown | Computed(num_observations) | No | None | Dictionary |
| `mesh_coordinate_scale_to_m` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `mode` | Parameter for mode in dictionary | Unknown | Computed(execution_mode) | No | None | Dictionary |
| `model_name` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `model_nprocs` | Re-apply normalized and command-line-resolved values after the raw YAML update.  Without this step, e.g. ``--Nens=40`` is silently replaced by the YAML value while constructing the canonical dictionary. | Unknown | Computed(model_nprocs) | No | None | Dictionary |
| `n_modeltasks` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `node_coords` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `nt` | Add model and analysis options to the canonical runtime context. | Unknown | Unknown | No | None | Dictionary |
| `num_param_vars` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `num_state_vars` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `number_obs_instants` | Parameter for number obs instants in dictionary | Unknown | Computed(num_observations) | No | None | Dictionary |
| `obs_data_path` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `obs_index` | Parameter for obs index in dictionary | Unknown | Computed(obs_idx) | No | None | Dictionary |
| `obs_max_time` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `obs_node_coords` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `obs_start_time` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `observations_available` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `observed_params` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `observed_vars` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `parallel_flag` | Controls parallel flag behavior in script logic | Unknown | Computed | No | None | Dictionary |
| `param_ens_spread` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `param_inflation_factor` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `parameter_estimation` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `params_vec` | check available parameters in the obseve_params list that need to be observed | list | [] | No | None | Variable |
| `partial` | Parameter for partial in dictionary | bool | True | No | None | Dictionary |
| `partitioned_io_flag` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `physics_bed_inference` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `physics_smb_inference` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `random_field_method` | Re-apply normalized and command-line-resolved values after the raw YAML update.  Without this step, e.g. ``--Nens=40`` is silently replaced by the YAML value while constructing the canonical dictionary. | Unknown | Computed(random_field_method) | No | None | Dictionary |
| `rel_vel_weight` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `restart_enabled` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `run_flag` | Controls run flag behavior in script logic | bool | True | No | None | Internal |
| `scalar_inputs` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `sequential_ensemble_initialization` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `sequential_run` | Parameter for sequential run in dictionary | bool | True | No | None | Dictionary |
| `serial_file_creation` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `sig_Q` | Parameter for sig q in dictionary | Unknown | Computed | No | None | Dictionary |
| `sig_model` | Parameter for sig model in dictionary | Unknown | Computed | No | None | Dictionary |
| `sig_obs` | Parameter for sig obs in dictionary | Unknown | Computed | No | None | Dictionary |
| `smb_blend_factor` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `smb_blend_ramp_time` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `smb_divergence_neighbors` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `smb_graph_neighbors` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `smb_history_length` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `smb_inference_start_time` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `smb_physical_bounds` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `smb_projection_basis` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `smb_spatial_regularization` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `smb_spinup_hold_factor` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `smb_temporal_regularization` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `state_estimation` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `state_inflation_factor` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `t` | Add model and analysis options to the canonical runtime context. | Unknown | Unknown | No | None | Dictionary |
| `taper_type` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `tikhonov_regularization_weight` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `total_state_param_vars` | Parameter for total state param vars in dictionary | Unknown | Unknown | No | None | Dictionary |
| `use_ensemble_pertubations` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `use_random_fields` | --- Ensemble Parameters --- | Unknown | Computed | No | None | Dictionary |
| `var_nd` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `vec_inputs` | Parameter for vec inputs in dictionary | Unknown | Unknown | No | None | Dictionary |
| `vel_idx` | Add model and analysis options to the canonical runtime context. | Unknown | Computed | No | None | Dictionary |
| `verbose` | Add model and analysis options to the canonical runtime context. | Unknown | Computed(_verbose) | No | None | Dictionary |
<!-- END: ICESEE-FLAGS -->
