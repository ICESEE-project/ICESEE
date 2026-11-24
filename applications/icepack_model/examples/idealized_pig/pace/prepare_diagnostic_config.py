"""Prepare a separate short diagnostic configuration from a frozen experiment.

No job is submitted. Noise amplitudes, physical scales and forcing drift must
be supplied explicitly; they are sensitivity settings, not calibrated defaults.
"""
import argparse
import copy
from pathlib import Path
import math
import yaml


def prepare_config(base, *, reference_data, reference_params, data_path, stop_time,
                   h_noise_sd=None, melt_noise_sd=None, noise_time=None, noise_length=None,
                   melt_trend=None, max_thickness_projection=None, variant="all", no_analysis=False,
                   fork_checkpoint=None, spread_relaxation=None, spread_fields=None,
                   stage_diagnostics=False):
    if variant not in {"gain", "balance", "forcing", "noise", "all"}:
        raise ValueError("unknown diagnostic variant")
    balance = variant != "gain"
    forcing = variant in {"forcing", "noise", "all"}
    noise = variant in {"noise", "all"}
    values = [stop_time]
    if balance: values += [max_thickness_projection]
    if forcing: values += [melt_trend]
    if noise: values += [h_noise_sd,melt_noise_sd,noise_time,noise_length]
    if any(value is None or not math.isfinite(float(value)) for value in values):
        raise ValueError("declare finite diagnostic settings for the selected variant")
    if (stop_time <= 0 or (balance and max_thickness_projection < 0)
            or (noise and (min(h_noise_sd,melt_noise_sd,noise_time)<0 or noise_length<=0))):
        raise ValueError("invalid diagnostic amplitude, time, length or projection limit")
    result=copy.deepcopy(base)
    model=result['modeling-parameters'];kw=result['enkf-parameters']
    if str(kw.get('enkf_observation_error_mode','')).lower()!='stochastic_r':
        raise ValueError("use an explicitly stochastic_R base to isolate gain changes")
    dt=float(model['timesteps_per_year']);horizon=float(model['num_years'])
    if stop_time>horizon or not math.isclose(stop_time/dt,round(stop_time/dt),rel_tol=0,abs_tol=1e-8):
        raise ValueError("stop_time must lie on the original time grid within its horizon")
    required=['h','u','v','s','basal_melt_field']
    if kw['vec_inputs']!=required or int(kw['num_state_vars'])!=5:
        raise ValueError("this diagnostic builder requires the five-field PIG state")
    if spread_relaxation is not None:
        if not math.isfinite(float(spread_relaxation)) or not 0 <= float(spread_relaxation) <= 1:
            raise ValueError("spread_relaxation must be finite and in [0, 1]")
        selected = ['h', 'basal_melt_field'] if spread_fields is None else list(spread_fields)
        if not selected or not set(selected).issubset(required):
            raise ValueError("spread_fields must identify nonempty PIG field names")
        kw.update(analysis_spread_relaxation=0.,
                  analysis_spread_relaxation_factors={name:float(spread_relaxation) for name in selected})
    elif spread_fields is not None:
        raise ValueError("spread_fields requires spread_relaxation")
    if stage_diagnostics:
        kw['analysis_stage_diagnostics'] = True
    for name in ('basal_melt_prior_scale_bounds','basal_melt_bounds'):
        if name not in kw:
            raise ValueError(f"declare {name} in the base experiment")
    reference_data=Path(reference_data).resolve();data_path=Path(data_path).resolve()
    if data_path.exists():
        raise ValueError("diagnostic data_path must be fresh")
    if (data_path==reference_data or data_path in reference_data.parents
            or reference_data in data_path.parents):
        raise ValueError("reference and output directories must be separate and non-nested")
    for name in ('true_nurged_states.h5','synthetic_obs.h5'):
        if not (reference_data/name).is_file():raise ValueError(f"missing reference {name}")
    if not Path(reference_params).is_file():raise ValueError("missing frozen reference configuration")
    kw.update(execution_mode=3,filter_type='EnKF',enkf_observation_error_mode='stochastic_R',
        enkf_gain_formulation='prescribed_r',local_analysis=False,localization_flag=False,
        analysis_enabled=not no_analysis,reference_data_path=str(reference_data),
        reference_params_file=str(Path(reference_params).resolve()),data_path=str(data_path),
        run_stop_time=float(stop_time),checkpoint_every=1,checkpoint_keep_last=0,
        mode3_export_ensemble_history=True,resume_from_checkpoint=False,
        generate_true_state=False,generate_synthetic_obs=False,generate_nurged_state=False)
    if balance:
        kw.update(pig_analysis_state_policy='diagnostic_balance',pig_report_analysis=True,
                  pig_max_thickness_projection=float(max_thickness_projection))
    if forcing:
        kw.update(basal_melt_forecast='estimated_depth_law',basal_melt_amplitude_trend=float(melt_trend))
    if noise:
        kw.update(process_noise_normalization='physical_ar1',process_noise_start='stationary',
            process_noise_correlation_time=float(noise_time),process_noise_schedule='every_step',
            process_noise_vars=['h','basal_melt_field'],sig_Q=[float(h_noise_sd),0.,0.,0.,float(melt_noise_sd)],
            graph_physical_length=True,random_field_method='graph',length_scale=[float(noise_length)]*5,
            graph_max_neighbors=256)
    if fork_checkpoint is not None:
        parent = Path(fork_checkpoint).resolve()
        if not (parent / 'manifest.json').is_file():
            raise ValueError("fork_checkpoint must be a closed committed checkpoint directory")
        kw.update(diagnostic_fork_checkpoint=str(parent),mode3_export_ensemble_history=False)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('base',type=Path)
    parser.add_argument('--variant',choices=('gain','balance','forcing','noise'),required=True)
    parser.add_argument('--no-analysis',action='store_true')
    parser.add_argument('--fork-checkpoint',type=Path)
    parser.add_argument('--spread-relaxation',type=float)
    parser.add_argument('--spread-fields',nargs='+')
    parser.add_argument('--stage-diagnostics',action='store_true')
    for name in ('reference-data','reference-params','data-path','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--stop-time',type=float,required=True)
    for name in ('h-noise-sd','melt-noise-sd','noise-time','noise-length','melt-trend','max-thickness-projection'):
        parser.add_argument('--'+name,type=float)
    args=vars(parser.parse_args());base=yaml.safe_load(args.pop('base').read_text())
    output=args.pop('output');result=prepare_config(base,**args)
    with output.open('x') as stream:yaml.safe_dump(result,stream,sort_keys=False)
    print(f'Prepared {output}; forcing horizon preserved. No job submitted.')


if __name__=='__main__':main()
