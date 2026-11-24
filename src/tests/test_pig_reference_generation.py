"""Separate reference jobs preserve the experiment and avoid ensemble setup."""
import ast
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import h5py
import numpy as np
import pytest
from ICESEE.applications.icepack_model.examples.idealized_pig._reference_generation import (
    validate_reference_only, finish_reference_generation, build_reference_context,
)

ROOT = Path(__file__).resolve().parents[2]
PIG = ROOT/'applications/icepack_model/examples/idealized_pig'
spec = importlib.util.spec_from_file_location('reference_builder', PIG/'pace/prepare_reference_jobs.py')
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


@pytest.mark.parametrize('kind', ['truth', 'background'])
def test_independent_generation_uses_same_full_horizon_and_physics(tmp_path, kind):
    base = {'modeling-parameters': {'num_years':82, 'timesteps_per_year':.05},
            'enkf-parameters': {'Nens':32, 'reference_data_path':'old', 'run_stop_time':30,
                                'sig_obs':.1, 'solver_setting':'unchanged'}}
    cfg = builder.reference_config(base, tmp_path/'data', kind)
    kw = cfg['enkf-parameters']
    assert validate_reference_only(kw, 1)
    assert cfg['modeling-parameters'] == base['modeling-parameters']
    assert kw['sig_obs'] == .1 and kw['solver_setting'] == 'unchanged'
    assert base['enkf-parameters']['Nens'] == 32
    with pytest.raises(ValueError, match='one MPI rank'):
        validate_reference_only(kw, 2)
    assert kw['generate_synthetic_obs'] == (kind == 'truth')


def test_background_can_create_fresh_hdf5_and_publish_manifest(tmp_path):
    source = ROOT/'src/EnKF/_generate_true_wrong_state.py'
    node = next(n for n in ast.parse(source.read_text()).body
                if isinstance(n, ast.FunctionDef) and n.name=='generate_true_wrong_state')
    import os, time
    ns = dict(np=np,h5py=h5py,os=os,time=time,record_phase=lambda *a,**kw:None)
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(source), 'exec'), ns)
    def generate(**kw):
        kw['statevec_nurged'][:] = np.arange(12).reshape(4,3)
    options = dict(nd=4,nt=2,dt=.05,data_path=str(tmp_path),
                   true_nurged_file=str(tmp_path/'true_nurged_states.h5'),
                   joint_estimation=False,localization_flag=False,num_state_vars=2,
                   generate_true_state=False,generate_nurged_state=True,
                   basal_melt_spatial_policy='fixed_initial_depth',
                   model_module=SimpleNamespace(generate_nurged_state=generate))
    ns['generate_true_wrong_state'](**options)
    result = finish_reference_generation(options, h5py)
    assert result['dataset']=='nurged_state' and result['generation_complete']
    with h5py.File(options['true_nurged_file']) as f:
        assert 'true_state' not in f
        np.testing.assert_array_equal(f['nurged_state'][:],np.arange(12).reshape(4,3))
    assert json.loads((tmp_path/'reference_manifest.json').read_text())['shape']==[4,3]


def test_reference_only_returns_before_observation_read_and_pool_allocation():
    text=(PIG/'mode3_runner.py').read_text()
    assert text.index('return finish_reference_generation') < text.index('hu_obs = f["hu_obs"][:]')


def test_direct_context_preserves_initial_native_fields_and_communicator():
    h0 = SimpleNamespace(dat=SimpleNamespace(data=np.ones(3)))
    fields = [object(),h0,*[object() for _ in range(14)]]
    comm=object(); seen=[]
    def initialize(**kw):
        seen.append(kw)
        return tuple(fields)
    options={'total_state_param_vars':5,'dt':.05}
    assert build_reference_context(options,initialize,comm)==15
    assert seen[0]['comm'] is comm and seen[0]['dt']==.05
    assert options['h0'] is h0 and options['solver'] is fields[-1]
