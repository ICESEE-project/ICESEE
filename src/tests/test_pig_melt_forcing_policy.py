"""Melt-prior reproducibility, bounds and experiment identities."""
import numpy as np
import pytest
from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import (
    initial_melt_scale, initial_melt_values, bounded_melt_values,
    estimated_melt_settings, melt_experiment_identity,
)


def settings(**overrides):
    return dict(basal_melt_forecast="estimated_field", vec_inputs=["h","u","v","s","basal_melt_field"],
                basal_melt_prior_scale_bounds=[0,2], basal_melt_bounds=[0,120],base_seed=42,**overrides)


def test_prior_is_member_keyed_and_independent_of_spatial_partition_or_order():
    options=settings()
    nominal=np.array([0.,5.,10.,20.])
    expected=initial_melt_values(nominal,options,7)
    partitioned=np.concatenate([initial_melt_values(part,options,7) for part in np.array_split(nominal,3)])
    np.testing.assert_array_equal(partitioned,expected)
    a={mid:initial_melt_scale(options,mid) for mid in range(32)}
    b={mid:initial_melt_scale(options,mid) for mid in reversed(range(32))}
    assert a==b and len(set(a.values()))==32
    assert all(0<=value<=2 for value in a.values())
    assert initial_melt_scale(dict(options,Nens=40),7)==a[7]


def test_melt_prior_never_mutates_nominal_forcing():
    nominal=np.array([0.,10.,20.]);before=nominal.copy()
    result=initial_melt_values(nominal,settings(),3)
    np.testing.assert_array_equal(nominal,before)
    assert not np.shares_memory(nominal,result)
    assert result[0]==0


def test_estimated_melt_bounds_reject_nonfinite_state():
    np.testing.assert_array_equal(bounded_melt_values([-1,20,121],settings()),[0,20,120])
    with pytest.raises(FloatingPointError):bounded_melt_values([np.nan],settings())


@pytest.mark.parametrize("key", ["basal_melt_prior_scale_bounds","basal_melt_bounds"])
@pytest.mark.parametrize("bad", [None,[0],[1,0],[-1,1],[0,np.inf]])
def test_estimation_requires_explicit_valid_scientific_settings(key,bad):
    options=settings();options[key]=bad
    with pytest.raises(ValueError,match=key):estimated_melt_settings(options)


def test_prescribed_policy_preserves_legacy_experiment_identity():
    assert initial_melt_scale({},3)==1
    assert melt_experiment_identity({}) is None
    assert melt_experiment_identity(settings())!=melt_experiment_identity(dict(settings(),basal_melt_bounds=[0,100]))
    assert melt_experiment_identity(settings())!=melt_experiment_identity(dict(settings(),enkf_observation_error_mode="stochastic_R"))


def test_spread_policy_changes_restart_identity_only_when_enabled():
    base = settings()
    old = melt_experiment_identity(base)
    assert melt_experiment_identity(dict(base, analysis_spread_relaxation=0.)) == old
    assert melt_experiment_identity(dict(base, analysis_stage_diagnostics=True)) == old
    half = melt_experiment_identity(dict(base, analysis_spread_relaxation=.5))
    assert half != old
    assert half != melt_experiment_identity(dict(base, analysis_spread_relaxation=.7))
    assert half == melt_experiment_identity(dict(base, analysis_spread_relaxation_factors={name:.5 for name in base['vec_inputs']}))


def test_diagnostics_refuse_unimplemented_runtime_paths():
    from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import configure_melt_experiment
    for override in ({'execution_mode':1}, {'execution_mode':3,'use_store_streaming_analysis':True}):
        with pytest.raises(NotImplementedError, match='persistent native'):
            configure_melt_experiment(dict(settings(), analysis_stage_diagnostics=True, **override))


def test_restart_identity_includes_error_weights_and_analysis_settings():
    from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import configure_melt_experiment
    original=configure_melt_experiment(settings(sig_obs=[5,.5,.1,.1,.1],inflation_factor=1.15))
    identity=original['scientific_experiment_id']
    changed=dict(original,sig_obs=[5,.5,.1,.2,.1])
    assert configure_melt_experiment(changed)['scientific_experiment_id']!=identity
    changed=dict(original,inflation_factor=1.1)
    assert configure_melt_experiment(changed)['scientific_experiment_id']!=identity


def test_balanced_analysis_is_explicit_and_restart_incompatible():
    from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import configure_melt_experiment
    old=configure_melt_experiment(settings())
    new=configure_melt_experiment(dict(settings(),pig_analysis_state_policy='diagnostic_balance',hThresh=30,execution_mode=3))
    assert new['frozen_analysis_vars']==['s','u','v']
    assert new['scientific_experiment_id']!=old['scientific_experiment_id']
    changed=configure_melt_experiment(dict(new,hThresh=31))
    assert changed['scientific_experiment_id']!=new['scientific_experiment_id']
    with pytest.raises(ValueError,match='execution_mode=3'):
        configure_melt_experiment(dict(new,execution_mode=2))


def test_depth_law_has_explicit_predictor_and_changes_restart_identity():
    from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import configure_melt_experiment, advance_melt_amplitude, melt_scenario
    options=dict(settings(),basal_melt_forecast='estimated_depth_law',execution_mode=3,
        pig_analysis_state_policy='diagnostic_balance',hThresh=30,basal_melt_amplitude_trend=.5)
    configured=configure_melt_experiment(options)
    assert advance_melt_amplitude(400.,20.,2.,configured)==21.
    assert advance_melt_amplitude(4000.,20.,2.,configured)==120.
    with pytest.raises(ValueError,match='unidentifiable'):advance_melt_amplitude(0.,0.,.05,configured)
    assert melt_scenario(5.95)=='control' and melt_scenario(6.)=='warm' and melt_scenario(15.)=='control'
    changed=configure_melt_experiment(dict(configured,basal_melt_amplitude_trend=.2))
    assert changed['scientific_experiment_id']!=configured['scientific_experiment_id']
    with pytest.raises(ValueError,match='explicit'):configure_melt_experiment(dict(configured,basal_melt_amplitude_trend=None))


def test_stop_time_preserves_full_forcing_horizon():
    from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import resolve_stop_step
    options=dict(nt=1640,dt=.05,num_years=82,run_stop_time=19.,checkpoint_every=1)
    assert resolve_stop_step(options)==380
    assert options['num_years']==82 and options['nt']==1640
    for invalid in (0,19.01,83,np.nan):
        with pytest.raises(ValueError):resolve_stop_step(dict(options,run_stop_time=invalid))
