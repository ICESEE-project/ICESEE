"""Independent truth/background generation without ensemble allocation."""
import json
from pathlib import Path


def build_reference_context(options, initialize_model, comm):
    names = ('h', 'h0', 's', 's0', 'u', 'bed', 'zF', 'grounded', 'floating',
             'A0', 'beta0', 'smb', 'basal_melt_field', 'Q', 'V', 'solver')
    values = initialize_model(**dict(options, comm=comm))
    if len(values) != len(names):
        raise ValueError('reference initializer returned an unexpected field layout')
    options.update(zip(names, values))
    options['comm'] = comm
    options['nd'] = int(options['h0'].dat.data.size*options['total_state_param_vars'])
    return options['nd']


def validate_reference_only(options, world_size):
    if not options.get('reference_generation_only', False):
        return False
    if world_size != 1:
        raise ValueError('reference_generation_only requires one MPI rank')
    if any(options.get(key) for key in ('resume_from_checkpoint', 'diagnostic_fork_checkpoint',
                                       'reference_data_path', 'run_stop_time')):
        raise ValueError('reference generation requires a fresh full trajectory')
    truth = bool(options.get('generate_true_state', False))
    background = bool(options.get('generate_nurged_state', False))
    if truth == background:
        raise ValueError('select exactly one of truth or background generation')
    if options.get('generate_synthetic_obs', False) and not truth:
        raise ValueError('generate observations in the truth job only')
    return True


def finish_reference_generation(options, h5_module):
    path = Path(options['data_path'])
    name = 'true_state' if options.get('generate_true_state') else 'nurged_state'
    expected = (int(options['nd']), int(options['nt']) + 1)
    with h5_module.File(options['true_nurged_file'], 'r') as handle:
        if handle[name].shape != expected:
            raise ValueError('generated reference shape disagrees with configured time grid')
    manifest = dict(generation_complete=True, dataset=name, shape=expected,
                    dt=float(options['dt']), final_time=float(options['nt'])*float(options['dt']),
                    basal_melt_spatial_policy=options.get('basal_melt_spatial_policy', 'evolving_depth'),
                    scientific_experiment_id=options.get('scientific_experiment_id'))
    temporary = path/'reference_manifest.json.tmp'
    temporary.write_text(json.dumps(manifest, indent=2)+'\n')
    temporary.replace(path/'reference_manifest.json')
    return manifest
