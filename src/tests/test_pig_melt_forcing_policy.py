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


def test_restart_identity_includes_error_weights_and_analysis_settings():
    from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import configure_melt_experiment
    original=configure_melt_experiment(settings(sig_obs=[5,.5,.1,.1,.1],inflation_factor=1.15))
    identity=original['scientific_experiment_id']
    changed=dict(original,sig_obs=[5,.5,.1,.2,.1])
    assert configure_melt_experiment(changed)['scientific_experiment_id']!=identity
    changed=dict(original,inflation_factor=1.1)
    assert configure_melt_experiment(changed)['scientific_experiment_id']!=identity
