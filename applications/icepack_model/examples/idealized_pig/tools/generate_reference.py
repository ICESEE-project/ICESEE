"""Generate one independent PIG reference, without assimilation dispatch."""
import argparse
import os
from pathlib import Path
import sys

os.environ['OMP_NUM_THREADS'] = '1'
PIG = Path(__file__).resolve().parents[1]
SOURCE = PIG.parents[3]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-F', '--force-params', type=Path, default=PIG/'params.yaml')
    args = parser.parse_args()
    params = args.force_params.resolve()
    os.chdir(PIG)
    sys.path[:0] = [str(PIG), str(SOURCE.parent), str(SOURCE/'config')]
    from config_loader import load_yaml_to_dict
    cfg = load_yaml_to_dict(str(params))
    configured = {key:value for section in cfg.values() if isinstance(section,dict)
                  for key,value in section.items()}
    data = Path(configured['data_path']).resolve()
    if data.exists() or data.is_symlink():
        raise ValueError(f'refusing an existing reference output: {data}')
    sys.argv = [sys.argv[0], '--force-params='+str(params), '--data_path='+str(data)]
    from ICESEE.config._utility_imports import icesee_kwargs
    import firedrake
    import h5py
    import numpy as np
    from mpi4py import MPI
    from ICESEE.src.utils.pyop2_compile_cache import install_consistent_compile_cache
    install_consistent_compile_cache()
    from ICESEE.applications.icepack_model.examples.idealized_pig import _icepack_enkf as model
    from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_model import initialize_model
    from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import configure_melt_experiment
    from ICESEE.applications.icepack_model.examples.idealized_pig._reference_generation import (
        validate_reference_only, finish_reference_generation, build_reference_context,
    )
    from ICESEE.src.EnKF._generate_true_wrong_state import generate_true_wrong_state
    from ICESEE.src.EnKF._generate_synthetic_observations import generate_synthetic_observations
    from ICESEE.src.utils.timestep_contract import prepare_timestep_contract

    options = dict(icesee_kwargs)
    configure_melt_experiment(options)
    if not validate_reference_only(options, MPI.COMM_WORLD.Get_size()):
        raise ValueError('select reference_generation_only in the reference configuration')
    dt = float(options['timesteps_per_year'])
    years = float(options['num_years'])
    if not np.isfinite([dt, years]).all() or dt <= 0 or years <= 0:
        raise ValueError('reference timestep and duration must be finite and positive')
    nt = round(years/dt)
    if nt < 1 or not np.isclose(nt*dt, years, rtol=0, atol=1e-8):
        raise ValueError('reference duration must lie on its configured time grid')
    options.update(dt=dt, nt=nt, t=np.arange(nt+1)*dt, model_module=model,
                   true_nurged_file=str(data/'true_nurged_states.h5'),
                   synthetic_obs_file=str(data/'synthetic_obs.h5'))
    options['observed_vars_params'] = list(options.get('observed_vars', []))+list(options.get('observed_params', []))
    options['all_observed'] = options['observed_vars_params']
    os.environ['ICESEE_RESULTS_DIR'] = str(data)
    prepare_timestep_contract(MPI.COMM_WORLD, str(data), time_grid=options['t'],
                              process_noise_normalization=options.get('process_noise_normalization', 'legacy_run_length'),
                              scientific_experiment_id=options.get('scientific_experiment_id'))
    print(f'Reference only: {years:g} years, dt={dt:g}, {nt} cycles; no ensemble allocation', flush=True)
    build_reference_context(options, initialize_model, firedrake.COMM_SELF)
    if os.environ.get('ICESEE_REFERENCE_SETUP_ONLY') == '1':
        print('PASS: reference setup complete; zero timesteps advanced', flush=True)
        return
    generate_true_wrong_state(**options)
    if options.get('generate_synthetic_obs', False):
        generate_synthetic_observations(**options)
    result = finish_reference_generation(options, h5py)
    print(result, flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception:
        import traceback
        import re
        report = traceback.format_exc()
        print(report, file=sys.stderr, flush=True)
        for name in set(re.findall(r'Compile errors in (/tmp/[^\s]+)', report)):
            path = Path(name)
            if path.is_file() and path.stat().st_uid == os.getuid():
                print('COMPILER STDERR:', path.read_text(errors='replace')[:262144], file=sys.stderr, flush=True)
        raise
