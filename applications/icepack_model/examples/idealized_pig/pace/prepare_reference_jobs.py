"""Prepare separate fresh truth and background configurations; submit nothing."""
import argparse
import copy
from pathlib import Path
import sys
import yaml


def reference_config(base, output, kind):
    if kind not in ('truth', 'background'):
        raise ValueError('unknown reference kind')
    result = copy.deepcopy(base)
    options = result['enkf-parameters']
    for key in ('reference_data_path', 'reference_params_file', 'run_stop_time',
                'diagnostic_fork_checkpoint'):
        options.pop(key, None)
    options.update(reference_generation_only=True, execution_mode=3, model_nprocs=1,
                   resume_from_checkpoint=False, mode3_export_ensemble_history=False,
                   analysis_stage_diagnostics=False, analysis_enabled=False,
                   basal_melt_spatial_policy='fixed_initial_depth', data_path=str(output),
                   generate_true_state=kind == 'truth', generate_nurged_state=kind == 'background',
                   generate_synthetic_obs=kind == 'truth')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('base', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parents[5]/'config'))
    from config_loader import load_yaml_to_dict
    base = load_yaml_to_dict(str(args.base.resolve()))
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    for kind in ('truth', 'background'):
        case = root/kind
        case.mkdir()
        (case/'params.yaml').write_text(yaml.safe_dump(
            reference_config(base, case/'data', kind), sort_keys=False))
    print(root)


if __name__ == '__main__':
    main()
