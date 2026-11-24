"""Model-free integration checks using the actual PIG forecast bodies."""
import ast
import __future__
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest

from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import (
    melt_forecast_mode, initial_melt_values, bounded_melt_values,
    configure_melt_experiment, bind_melt_reference,
)

PIG=Path(__file__).resolve().parents[2]/'applications/icepack_model/examples/idealized_pig'


class Field:
    def __init__(self,shape):
        self.values=np.zeros(shape)
        self.dat=SimpleNamespace(data=self.values,data_ro=self.values)
    def copy(self,deepcopy=False):
        result=Field(self.values.shape);result.values[:]=self.values;return result
    def interpolate(self,value):
        self.values[:]=value;return self
    def __mul__(self,other):
        return self.values*np.asarray(other)


def field(values):
    result=Field(np.shape(values));result.values[:]=values;return result


def load_functions(filename,names,namespace):
    namespace["bind_melt_reference"] = bind_melt_reference
    tree=ast.parse((PIG/filename).read_text())
    nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert len(nodes)==len(names)
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(PIG/filename),'exec',flags=__future__.annotations.compiler_flag),namespace)


def options():
    return dict(basal_melt_forecast='estimated_field',basal_melt_prior_scale_bounds=[0,2],
                basal_melt_bounds=[0,120],vec_inputs=['h','u','v','s','basal_melt_field'],
                base_seed=42,dt=.05,joint_estimation=False)


def model_namespace():
    ns=dict(Function=Field,initial_melt_values=initial_melt_values,bounded_melt_values=bounded_melt_values)
    load_functions('_icepack_model.py',{'initial_estimated_melt','estimated_melt_forecast'},ns)
    return ns


def test_native_forecast_uses_the_analyzed_forcing_and_current_mask():
    ns=model_namespace();captured=[]
    ctx=SimpleNamespace(Q=(3,),forward_solver=object(),smb=None,bed=None,h0=None)
    def solver(solver,h,u,smb,melt,bed,dt,h0,kwargs):
        captured.append(melt.dat.data_ro.copy())
        return h,u,field([1,2,3]),None,None
    ns.update(_shared_context=lambda *a:ctx,_recompute_floating=lambda *a:np.array([1,0,1]),
              _icepack_run_kwargs=lambda *a:{},melt_forecast_mode=melt_forecast_mode,
              _basal_melt_for_step=lambda *a,**kw:field([20,20,20]),
              Icepack=solver,IcepackForecastFields=lambda **kw:SimpleNamespace(**kw),EXPERIMENT_WRONG='wrong')
    load_functions('_icepack_native.py',{'forecast_member'},ns)
    state=SimpleNamespace(thickness=field([100]*3),surface=field([10]*3),velocity=field([1]*3),basal_melt=field([5,90,130]))
    ns['forecast_member'](state,10,topology=object(),icesee_kwargs=options())
    np.testing.assert_array_equal(captured[-1],[5,0,120])
    np.testing.assert_array_equal(state.basal_melt.dat.data_ro,[5,90,130])
    prescribed=dict(options(),basal_melt_forecast='prescribed_wrong')
    ns['forecast_member'](state,10,topology=object(),icesee_kwargs=prescribed)
    np.testing.assert_array_equal(captured[-1],[20,20,20])


def test_generic_forecast_uses_the_same_estimated_melt_provider():
    ns=model_namespace();captured=[];n=3
    indices={name:np.arange(i*n,(i+1)*n) for i,name in enumerate(options()['vec_inputs'])}
    def solver(solver,h,u,smb,melt,bed,dt,h0,kwargs):
        captured.append(melt.dat.data_ro.copy());return h,u,s,None,None
    h=field([100]*3);s=field([10]*3)
    ns.update(np=np,melt_forecast_mode=melt_forecast_mode,icesee_get_index=lambda **kw:(None,indices,None),
              mf=SimpleNamespace(flotationHeight=lambda *a:0,flotationMask=lambda *a:(np.array([1,0,1]),None)),
              Icepack=solver,profile_columns=lambda *a:[],EXPERIMENT_WRONG='wrong')
    load_functions('_icepack_model.py',{'run_model'},ns)
    kwargs=dict(options(),Q=(3,),V=(3,2),k=10,nt=20,h=h,h0=h,s=s,bed=object(),save_steps=[],
                ens_id=0,t=np.arange(21)*.05,water_to_ice=1.,total_state_param_vars=5,num_state_vars=5)
    vector=np.concatenate(([100]*3,[1]*3,[2]*3,[10]*3,[5,90,130]))
    result=ns['run_model'](vector,**kwargs)
    np.testing.assert_array_equal(captured[-1],[5,0,120])
    np.testing.assert_array_equal(result['basal_melt_field'],[5,0,120])


def test_melt_configuration_preserves_fixed_initial_physics_and_builds_bounds():
    cfg=configure_melt_experiment(options())
    assert cfg['initial_perturbation_vars']==[]
    assert cfg['analysis_bounds']=={'basal_melt_field':[0.,120.]}
    identity=cfg['scientific_experiment_id']
    assert identity and configure_melt_experiment(cfg)['scientific_experiment_id']==identity
    cfg['basal_melt_prior_scale_bounds']=[.5,1.5]
    assert configure_melt_experiment(cfg)['scientific_experiment_id']!=identity
    with pytest.raises(ValueError,match='disagree'):
        configure_melt_experiment(dict(options(),analysis_bounds={'basal_melt_field':[0,100]}))
