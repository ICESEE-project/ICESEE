"""Exercise the actual melt expression with evolving and frozen geometry."""
import ast
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import (
    bind_melt_reference, melt_spatial_policy, melt_experiment_identity, validate_reference_config,
)

PIG = Path(__file__).resolve().parents[2]/'applications/icepack_model/examples/idealized_pig'


class Field:
    def __init__(self, space):
        self.space = space
        self.values = np.zeros(4)
    def interpolate(self, expression):
        self.values[:] = np.asarray(expression)
        return self
    def function_space(self): return self.space
    def __array__(self, dtype=None): return np.asarray(self.values, dtype=dtype)
    def __sub__(self, other): return self.values-np.asarray(other)
    def __ge__(self, other): return self.values>=other
    def __le__(self, other): return self.values<=other
    def __rsub__(self, other): return other-self.values


def setup(policy='fixed_initial_depth'):
    node = next(n for n in ast.parse((PIG/'_icepack_model.py').read_text()).body
                if isinstance(n, ast.FunctionDef) and n.name=='BasalMeltRate')
    ns = dict(firedrake=SimpleNamespace(Function=Field, conditional=np.where),
              _melt_max_constant=lambda value:value, melt_spatial_policy=melt_spatial_policy,
              EXPERIMENT_TRUE='true', EXPERIMENT_WRONG='wrong')
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(PIG/'_icepack_model.py'),'exec'),ns)
    options = dict(basal_melt_spatial_policy=policy, dt=.05, num_years=82, bmr_increase_time=0,
                   s0=Field('native'), h0=Field('native'))
    options['h0'].values[:] = [600,475,425,550]
    current_h = Field('native');current_h.values[:] = [200,800,600,350]
    current_s = Field('native')
    return ns['BasalMeltRate'], options, current_s, current_h


@pytest.mark.parametrize('year',[0,6,15,40,60,64,66,69,70,76,82])
def test_fixed_pattern_scales_linearly_across_all_former_scenario_boundaries(year):
    function, options, s, h = setup()
    initial_h = options['h0'].values.copy()
    melt, maximum = function(options, round(year/.05), np.ones(4), 'native', s, h, scenario='warm')
    expected = 20+80*year/82
    assert maximum == pytest.approx(expected)
    np.testing.assert_allclose(melt.values, expected*np.array([1.,.5,0.,1.]))
    np.testing.assert_array_equal(options['h0'].values, initial_h)


def test_wrong_forcing_keeps_same_pattern_with_constant_twenty_amplitude():
    function, options, s, h = setup()
    for year in (0,40,64,69,82):
        melt, maximum = function(options, round(year/.05), np.ones(4), 'native', s, h,
                                 scenario='warm',experiment='wrong')
        assert maximum == 20
        np.testing.assert_allclose(melt.values,[20,10,0,20])


def test_current_grounding_mask_remains_authoritative():
    function, options, s, h = setup()
    floating = np.array([0.,1.,1.,1.])
    melt, maximum = function(options, 1400, floating, 'native', s, h, scenario='warm')
    np.testing.assert_allclose(melt.values,maximum*np.array([0,.5,0,1]))


def test_default_evolving_depth_policy_retains_current_geometry_and_warm_thresholds():
    function, options, s, h = setup('evolving_depth')
    options.pop('basal_melt_spatial_policy')
    melt, maximum = function(options,0,np.ones(4),'native',s,h,scenario='warm')
    np.testing.assert_allclose(melt.values,maximum*np.array([0,1,1,0]))
    assert melt_spatial_policy({})=='evolving_depth'


def test_fixed_pattern_refuses_missing_or_foreign_mesh_references():
    function, options, s, h = setup()
    with pytest.raises(ValueError,match='immutable'):
        function(dict(options,h0=None),0,np.ones(4),'native',s,h)
    with pytest.raises(ValueError,match='current model space'):
        function(dict(options,s0=Field('reference-comm-self')),0,np.ones(4),'native',s,h)


def test_native_binding_never_reuses_world_reference_geometry_or_mutates_options():
    _, options, s, h = setup()
    context = SimpleNamespace(h0=h,s0=s)
    bound = bind_melt_reference(options,context)
    assert bound['h0'] is h and bound['s0'] is s
    assert options['h0'] is not h
    legacy = {'basal_melt_spatial_policy':'evolving_depth'}
    assert bind_melt_reference(legacy,context) is legacy


def test_restart_identity_and_cached_reference_compatibility(tmp_path):
    legacy = dict(num_years=82,timesteps_per_year=.05,sig_obs=[.1],obs_start_time=15,obs_max_time=65,freq_obs=2)
    fixed = dict(legacy,basal_melt_spatial_policy='fixed_initial_depth')
    assert melt_experiment_identity(legacy) is None
    assert melt_experiment_identity(fixed) is not None
    config = tmp_path/'params.yaml';config.write_text(yaml.safe_dump({'modeling-parameters':legacy}))
    ref = tmp_path/'old';ref.mkdir()
    for name in ('true_nurged_states.h5','synthetic_obs.h5'):(ref/name).touch()
    options = dict(fixed,reference_params_file=str(config),reference_data_path=str(ref),data_path=str(tmp_path/'fresh'))
    with pytest.raises(ValueError,match='basal_melt_spatial_policy'):
        validate_reference_config(options)
    config.write_text(yaml.safe_dump({'modeling-parameters':fixed}))
    validate_reference_config(options)
    with pytest.raises(ValueError,match='basal_melt_spatial_policy'):
        validate_reference_config(dict(options,basal_melt_spatial_policy='evolving_depth'))
