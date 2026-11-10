# Configuration reference and precedence

ICESEE reads a YAML file with `physical-parameters`, `modeling-parameters`, and `enkf-parameters` sections. The first two contain application physics/time definitions; the third contains ensemble, analysis, storage, and execution settings. [`config/_utility_imports.py`](../../config/_utility_imports.py) flattens them into one mutable `icesee_kwargs` dictionary. There is no independent nested public runtime configuration object.

## Resolution order

1. `-F PATH` or `--force-params=PATH` selects the YAML file; the default is `params.yaml` in the application's working directory.
2. `load_yaml_to_dict` resolves an optional top-level `extends`. A relative base path is resolved relative to the file declaring it. Dictionaries merge recursively; lists/scalars replace their base value. Cyclic inheritance raises. Archive the complete inheritance chain, or a fully resolved copy.
3. The loader computes defaults and derived values, then flattens the physical, modeling, and ensemble sections in that order. Later sections replace duplicate keys. A key repeated across sections is therefore ambiguous to a reader even if its software outcome is defined.
4. Explicit named CLI values for `Nens`, `data_path`, `model_nprocs`, and `execution_mode` are reapplied after flattening. The loader normalizes selected values such as the random-field method.
5. Other `--key=value` or `--key value` options override existing context keys. A bare `--key` becomes YAML `true`. Unknown keys raise instead of becoming unused parameters. Values are parsed with `yaml.safe_load`, then coerced against the existing value's type where applicable. Quote lists, dictionaries, and paths containing spaces in the shell.
6. The application resolves its final model grid, dimensions and inputs. The dispatcher checks recorded generic overrides and raises if an application silently replaced a requested value. Derived quantities may need a different application option rather than a direct override.

The loader also prepares `data_path` before applying all generic overrides. Its restart-preservation decision specially reads `resume_from_checkpoint` from raw argv. Do not import it to inspect defaults; use the generated key index or parse the YAML with `config_loader`. See [run directory protection](operations.md#protecting-run-directories).

## Canonical selectors and state definition

| Key | Loader fallback or rule | Operational meaning |
|---|---|---|
| `execution_mode` | `1`; named CLI choices `0,1,2,3` | Chooses the runner. Mode 3 requires the selected application's registration and supported policy. |
| `Nens` | `1`; example YAML normally overrides | Number of scientific members. A one-member topology test does not supply a well-sampled EnKF. |
| `model_name` | `model` | Registered family identifier: `icepack`, `issm`, `lorenz`, `flowline`. |
| `vec_inputs` | Required lookup | Ordered state/parameter block names. Changing order changes row meaning. |
| `num_state_vars`, `num_param_vars` | `1`, `0` | Block counts used by indexing and initialization. Must match adapter output. |
| `var_nd`, `scalar_inputs` | `None`, `[]` | Optional variable-size/scalar metadata. Check compatibility of every selected analysis control. |
| `state_distribution` | Generic ownership helper defaults to `distributed` | Declare `replicated` for model groups that each hold the same complete member. |
| `ranks_per_model` | `None` | Modes 1/2 resource-plan control; explicit integer or supported `auto` policy. |
| `model_nprocs` | `0` in loader unless requested | Application-specific solve resource setting; native PIG uses it for spatial topology. ISSM uses internal solver resources. |
| `filter_type` | Application-defined | Method name. Availability of serial variants does not establish their support in every parallel/native policy. |

`dt`, `nt`, `num_years` and `timesteps_per_year` are application-sensitive. The loader initially derives `dt=1/timesteps_per_year` and `nt=num_years*timesteps_per_year` unless `nt` is explicit; raw modeling values and application code can replace them. Lorenz explicitly computes `nt=num_years/dt`. PIG resolves its own interval and observation schedule. Inspect the final time vector and logs, rather than inferring a run duration from an early default.

## Analysis and uncertainty

| Keys | Fallback or behavior | What must be recorded |
|---|---|---|
| `base_seed` | `42` | Reproducible stream namespace; a model's `seed` can separately affect its own initialization. |
| `sig_model`, `sig_obs`, `sig_Q` | Loader-generated arrays when omitted | Initial/model uncertainty, observation factor scale, process noise; units and block order must be explained. |
| `initial_spread_factor` | `1.0` | Generic initialization scaling; does not imply every native adapter perturbs the same fields. |
| `random_field_method` | `fft`, alias fallback `enkf_field_method` | `fft` or `graph`; graph requires physical coordinates. |
| `length_scale`, `Q_rho` | `[]`, `1.0` | Spatial/temporal noise behavior, interpreted by the selected generator. |
| `process_noise_schedule` | `observations` | `observations`, `every_step`, or `none`; the schedule refers to forecast cycles ending at state columns. |
| `enkf_observation_error_mode` | Derived from `use_ensemble_pertubations` | Set explicitly to avoid changing scientific semantics when a legacy boolean changes. |
| `global_analysis`, `local_analysis` | `true`, `false` | Identify actual branches in the selected runner. They are not universal capability declarations. |
| `localization_flag`, `localized_vars`, `localization_radius` | `false`, `[]`, `None` | Variable neighborhoods, coordinate units, radius policy and active group counts. |
| `inflation_factor` | `1.0` | Common anomaly inflation. Per-state/parameter/bed values may override it. |
| `analysis_relaxation_factor`, `analysis_relaxation_factors` | Analysis helper default `1.0`, `{}` | Multiply increments before inflation; explicit opt-in changes the posterior. |
| `analysis_increment_limits`, `frozen_analysis_vars` | Helper `{}`, loader `[]` | Per-variable increment caps; forecast values restored for frozen analysis blocks. |
| `analysis_svd_energy` | `0.999` | Streamed spectral retained-energy threshold; inspect direct/local paths before changing. |

See [numerics](numerics.md) for exact factor equations and post-analysis ordering. The [complete loader index](configuration-key-index.md) lists the source fallback expressions. Those expressions describe absent-key behavior; they are not a universal runtime schema. Raw YAML flattening and model code remain part of the resolution process.

## Storage controls

| Key | Loader fallback | Effect and boundary |
|---|---|---|
| `data_path` | `_modelrun_datasets` | Dedicated run output path, interpreted from application cwd unless absolute. A fresh run cleans it. |
| `ensemble_history_mode` | `auto` | Mode-2 full or rolling retention; auto can choose rolling when storage projection is unsafe. Explicit `full` is needed for full-member history studies. |
| `ensemble_storage_dtype` | `float64` | On-disk precision, not a promise that analysis uses that precision. Strict parity uses float64. |
| `batch_size` | `2` before runtime normalization | Active I/O window, not ensemble mini-batch size or a physical timestep. |
| `analysis_memory_budget_mb` | `256` | Main row-block arithmetic budget. It does not bound model memory or all finalization work. |
| `analysis_row_chunk_size`, `observation_row_chunk_size` | `0`, `0` | Zero derives chunks from budgets; a positive explicit value bypasses that estimate. |
| `analysis_finalize_memory_budget_mb`, `local_analysis_memory_budget_mb` | `2048`, `2048` | Separate limits for finalization/local-analysis staging. |
| `observation_operator_storage`, `synthetic_observation_storage` | `indices`, `compact` | Avoid a dense global `H` and full-state synthetic observation arrays where supported. |
| `ensemble_finalize_mode` | `vds` | Full-history metadata view over Mode-2 shards; source shards remain necessary. |
| `storage_safety_factor`, `fail_on_insufficient_storage` | `1.10`, `false` | Capacity planning margin and refusal policy; production may explicitly request a conservative margin and hard refusal. |
| `h5_file_compression`, `h5_file_compression_level` | `None`, `4` | Storage filter; accepted normalization includes gzip, lzf and szip. gzip level must be 0–9. Parallel/backend support must be verified. |
| `partitioned_io_flag` | `false` | A distinct partitioned path; some global bed-domain/inference combinations require it false. |

`full` history and restart retention solve different needs. Full snapshots permit retrospective member spread; restart needs only enough complete committed states to resume. Deleting historical members prevents recovering their spread even if the mean time series remains.

## Restart and native checkpoint controls

| Key | Fallback | Meaning |
|---|---|---|
| `restart_enabled` | `true` | Runtime checkpoint/reuse capability; it does not protect data from loader cleanup. |
| `resume_from_checkpoint` | `false` | Explicitly preserve output and request supported resume. Not implemented by every registered native policy. |
| `force_fresh_start` | `false` | Reject/reinitialize reuse in paths that consume it; use a new output directory for a new experiment. |
| `checkpoint_every` | `1` | Native policy writes each selected completed cycle; final cycle is always eligible. |
| `checkpoint_keep_last` | `2` | Retain newest restart checkpoints. Zero means retain every written checkpoint. |
| `checkpoint_keep_analysis` | `false` | Pin analysis checkpoints in addition to rolling restart retention. |
| `mode3_export_ensemble_history` | `false` | Optional aggregate export. Full export requires every cycle written and retained. |
| `mode3_export_row_chunk_size` | `8192` | Bounds aggregate row reads. |
| `k_start_override` | `None` | Explicit restart position where implemented; does not make arbitrary old files compatible. |

Checkpoint identity includes run/layout/schedule policy. For native PIG it also distinguishes the corrected time grid and canonical process noise. Preserve the same input values, observation plan, ensemble size, time discretization and ownership assumptions when resuming. A portability or decomposition change must be validated separately.

## Observation and inference settings

`observed_vars` and `observed_params` select observed blocks. `observations_available` selects external observation handling; when true, the application must load schedule metadata as well as values. `obs_data_path` alone does not populate all necessary time indices. `generate_true_state`, `generate_nurged_state`, and `generate_synthetic_obs` control synthetic/reference artifacts; their generated files must not be confused with independently measured observations.

`obs_start_time`, `obs_max_time`, and `freq_obs` define requested physical times. Inspect the mapped schedule and analysis events. Bed surveys add snapshot/stride/spacing/index/mask controls, including track half width and coordinate scaling. Spatial strides in kilometers and geometry coordinates in meters need an explicit conversion; retain `mesh_coordinate_scale_to_m` where used.

Hybrid and inference options include `inversion_flag`, model-specific inversion activation and weights, `inference_plugin_enabled`, `physics_smb_inference`, `physics_bed_inference`, and `bed_update_mode`. The loader provides many regularization, blending, timing, projection, and physical-bound options. Their defaults are discoverable in the key index; a scientific protocol should list the subset actually active and explain their assumptions. A registered mode may refuse an unsupported feature rather than apply it.

## Change-control rule

Treat changes to seeds, uncertainty amplitudes, observation-error factor, localization, inference priors, physical bounds, state ordering, solver tolerances, or the physical time grid as scientific changes. Treat row chunks, equivalent I/O layout, or communicator scheduling as implementation changes only after same-input numerical comparisons establish that they preserve intended behavior. Record both categories in the run manifest, and keep a small deterministic control case available for each.
