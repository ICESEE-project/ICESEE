"""Reference safety and diagnostic configuration construction without model imports."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest
import yaml
from ICESEE.src.utils.reference_data import validate_reference_paths
from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import validate_reference_config

ROOT=Path(__file__).resolve().parents[2]


def test_reference_data_cannot_be_cleaned_as_output(tmp_path):
    ref=tmp_path/'reference';ref.mkdir()
    for output in (ref,tmp_path,ref/'new'):
        with pytest.raises(ValueError,match='separate'):validate_reference_paths(output,ref)
    alias=tmp_path/'alias';alias.symlink_to(ref,target_is_directory=True)
    with pytest.raises(ValueError):validate_reference_paths(alias,ref)
    validate_reference_paths(tmp_path/'new',ref)


def test_frozen_reference_configuration_is_required_and_checked(tmp_path):
    ref=tmp_path/'reference';ref.mkdir()
    for name in ('true_nurged_states.h5','synthetic_obs.h5'):(ref/name).touch()
    options=dict(reference_data_path=str(ref),data_path=str(tmp_path/'new'),num_years=82,
        timesteps_per_year=.05,sig_obs=[5,.5,.1,.1,.1],obs_start_time=15,obs_max_time=65,freq_obs=2)
    with pytest.raises(ValueError,match='requires'):validate_reference_config(options)
    original={key:value for key,value in options.items() if key not in ('reference_data_path','data_path')}
    config=tmp_path/'params.yaml';config.write_text(yaml.safe_dump({'modeling-parameters':original}))
    options['reference_params_file']=str(config)
    validate_reference_config(options)
    with pytest.raises(ValueError,match='num_years'):validate_reference_config(dict(options,num_years=19))
    with pytest.raises(ValueError,match='sig_obs'):validate_reference_config(dict(options,sig_obs=[1]*5))


def test_builder_preserves_physical_horizon_and_original_yaml(tmp_path):
    path=ROOT/'applications/icepack_model/examples/idealized_pig/pace/prepare_diagnostic_config.py'
    spec=importlib.util.spec_from_file_location('diagnostic_config',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    base={'modeling-parameters':dict(num_years=82,timesteps_per_year=.05,hThresh=30),
        'enkf-parameters':dict(Nens=32,num_state_vars=5,vec_inputs=['h','u','v','s','basal_melt_field'],
            basal_melt_bounds=[0,120],basal_melt_prior_scale_bounds=[0,2],inflation_factor=1.15,
            enkf_observation_error_mode='stochastic_R')}
    ref=tmp_path/'references';ref.mkdir()
    for name in ('true_nurged_states.h5','synthetic_obs.h5'):(ref/name).touch()
    config=tmp_path/'params.yaml';config.write_text(yaml.safe_dump(base))
    result=module.prepare_config(base,reference_data=ref,reference_params=config,data_path=tmp_path/'out',
        stop_time=19,h_noise_sd=.01,melt_noise_sd=.1,noise_time=2,noise_length=1000,
        melt_trend=0,max_thickness_projection=1)
    assert base['enkf-parameters'].get('run_stop_time') is None
    assert result['modeling-parameters']==base['modeling-parameters']
    kw=result['enkf-parameters']
    assert kw['run_stop_time']==19 and kw['Nens']==32 and kw['filter_type']=='EnKF'
    assert kw['enkf_gain_formulation']=='prescribed_r' and kw['process_noise_vars']==['h','basal_melt_field']
    assert kw['sig_Q']==[.01,0,0,0,.1]
    spread = module.prepare_config(base, reference_data=ref, reference_params=config,
        data_path=tmp_path/'spread', stop_time=15.05, variant='balance', max_thickness_projection=1,
        spread_relaxation=.5, spread_fields=['h','basal_melt_field'], stage_diagnostics=True)
    assert spread['enkf-parameters']['analysis_spread_relaxation_factors']=={'h':.5,'basal_melt_field':.5}
    assert spread['enkf-parameters']['analysis_stage_diagnostics'] is True
    assert spread['modeling-parameters']==base['modeling-parameters']
    assert spread['enkf-parameters']['inflation_factor']==1.15
    assert 'analysis_spread_relaxation_factors' not in base['enkf-parameters']
    with pytest.raises(ValueError, match='spread_relaxation'):
        module.prepare_config(base, reference_data=ref, reference_params=config,
            data_path=tmp_path/'bad', stop_time=15.05, variant='balance', max_thickness_projection=1,
            spread_relaxation=1.1)


def test_preanalysis_diagnostic_fork_has_explicit_lineage_and_boundary(tmp_path):
    import json
    from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import diagnostic_fork_info
    ref=tmp_path/'reference';ref.mkdir()
    for name in ('true_nurged_states.h5','synthetic_obs.h5'):(ref/name).touch()
    original=dict(num_years=82,timesteps_per_year=.05,sig_obs=[5,.5,.1,.1,.1],
        obs_start_time=15,obs_max_time=65,freq_obs=2,process_noise_schedule='observations')
    config=tmp_path/'params.yaml';config.write_text(yaml.safe_dump({'enkf-parameters':original}))
    parent=ref/'history'/'steps'/'checkpoint_00000298';parent.mkdir(parents=True)
    manifest=dict(complete=True,number_of_members=32,run_id='original',timestep=298,
        metadata=dict(nt=1640,next_step=299,model_time=14.95,analysis_events_completed=0))
    (parent/'manifest.json').write_text(json.dumps(manifest))
    options=dict(original,nt=1640,dt=.05,Nens=32,obs_index=[300,340,380],number_obs_instants=3,
        run_stop_time=15.05,checkpoint_every=1,diagnostic_fork_checkpoint=str(parent),
        reference_data_path=str(ref),reference_params_file=str(config),data_path=str(tmp_path/'new'),
        mode3_export_ensemble_history=False)
    fork=diagnostic_fork_info(options)
    assert fork['next_step']==299 and fork['parent_run_id']=='original' and len(fork['manifest_sha256'])==64
    with pytest.raises(ValueError,match='complete history'):
        diagnostic_fork_info(dict(options,mode3_export_ensemble_history=True))
    manifest['metadata']['analysis_events_completed']=1
    (parent/'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match='pre-first-analysis'):diagnostic_fork_info(options)


@pytest.mark.parametrize('reference_from_cli',[False,True])
def test_config_refuses_cleanup_before_touching_reference_files(tmp_path,reference_from_cli):
    import os,subprocess,sys
    example=tmp_path/'example';example.mkdir()
    source=ROOT/'applications/lorenz_model/examples/lorenz96/params.yaml'
    base=yaml.safe_load(source.read_text())
    data=example/'_modelrun_datasets';data.mkdir();marker=data/'retained.h5';marker.write_bytes(b'preserve')
    if not reference_from_cli:base['enkf-parameters']['reference_data_path']=str(data)
    (example/'params.yaml').write_text(yaml.safe_dump(base))
    command=[sys.executable,'-c','import ICESEE.config._utility_imports',f'--data_path={data}']
    if reference_from_cli:command.append(f'--reference_data_path={data}')
    result=subprocess.run(command,cwd=example,env=dict(os.environ),capture_output=True,text=True,timeout=30)
    assert result.returncode!=0
    assert 'separate' in result.stderr, result.stdout+result.stderr
    assert marker.read_bytes()==b'preserve'
