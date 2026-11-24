"""Exercise physical finalization without a Firedrake installation."""
import ast
from pathlib import Path
from types import SimpleNamespace as NS
import numpy as np
from ICESEE.applications.icepack_model.examples.idealized_pig._melt_forcing import bind_melt_reference


def test_finalization_bounds_geometry_and_solves_before_checkpoint_state():
    path=Path(__file__).resolve().parents[2]/'applications/icepack_model/examples/idealized_pig/_icepack_native.py'
    node=next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='finalize_analysis_state')
    class Field:
        def __init__(self,a):self.a=np.array(a,dtype=float);self.dat=NS(data=self.a,data_ro=self.a)
        def copy(self,deepcopy=False):return Field(self.a.copy())
    calls=[];updates=[]
    def solve(**kwargs):calls.append(kwargs);return Field([[1,2],[3,4]])
    ctx=NS(Q=None,bed=None,s0=None,beta0=np.ones(2),A0=None,forward_solver=NS(diagnostic_solve=solve))
    state=NS(thickness=Field([-70,100]),surface=Field([-90,20]),velocity=Field([[9,9],[9,9]]),basal_melt=Field([120,20]),update_from_forecast=updates.append)
    ns={'__package__':'ICESEE.applications.icepack_model.examples.idealized_pig','np':np,
        '_shared_context':lambda *a:ctx,'icepack':NS(compute_surface=lambda **kw:Field(kw['thickness'].a*.1),interpolate=lambda a,Q:a),
        'mf':NS(flotationHeight=lambda *a:0,flotationMask=lambda *a:(np.array([1,0]),np.array([0,1])),reduceNearGLBeta=lambda *a:np.ones(2)),
        'melt_forecast_mode':lambda opts:'estimated_field',
        'estimated_melt_forecast':lambda melt,mask,Q,opts:Field(np.clip(melt.a,0,120)*mask),
        'IcepackForecastFields':lambda **kw:NS(**kw)}
    ns['bind_melt_reference'] = bind_melt_reference
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),'exec'),ns)
    options=dict(pig_analysis_state_policy='diagnostic_balance',obs_index=[5],number_obs_instants=1,hThresh=30,GLThresh=10,uThresh=1)
    ns['finalize_analysis_state'](state,None,3,topology=None,icesee_kwargs=options)
    assert not calls and not updates
    ns['finalize_analysis_state'](state,None,4,topology=None,icesee_kwargs=options)
    assert len(calls)==len(updates)==1
    np.testing.assert_array_equal(updates[0].thickness.a,[30,100])
    np.testing.assert_array_equal(updates[0].surface.a,[3,10])
    np.testing.assert_array_equal(updates[0].basal_melt.a,[120,0])
    assert set(calls[0])=={'velocity','thickness','surface','beta','fluidity','uThresh','floating','grounded'}
    ns['finalize_analysis_state'](state,None,4,topology=None,icesee_kwargs={})
    assert len(calls)==1


def test_excessive_projection_is_refused_before_solving():
    import pytest
    path=Path(__file__).resolve().parents[2]/'applications/icepack_model/examples/idealized_pig/_icepack_native.py'
    node=next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='finalize_analysis_state')
    state=NS(thickness=NS(dat=NS(data_ro=np.array([-70.,100.]))))
    ns={'__package__':'ICESEE.applications.icepack_model.examples.idealized_pig','np':np,'_shared_context':lambda *args:object()}
    ns['bind_melt_reference'] = bind_melt_reference
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),'exec'),ns)
    with pytest.raises(RuntimeError,match='projection'):
        ns['finalize_analysis_state'](state,None,4,topology=None,icesee_kwargs=dict(
            pig_analysis_state_policy='diagnostic_balance',obs_index=[5],hThresh=30,pig_max_thickness_projection=1))


def test_melt_only_noise_does_not_repeat_velocity_solve():
    path=Path(__file__).resolve().parents[2]/'applications/icepack_model/examples/idealized_pig/_icepack_native.py'
    cls=next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.ClassDef) and n.name=='_IdealizedPigAdapter')
    method=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='finalize_forecast_noise')
    calls=[];updates=[]
    state=NS(thickness=object(),velocity=object(),surface=object(),basal_melt=object(),
        registry=NS(synchronize_ghosts=lambda:calls.append('halo')),update_from_forecast=updates.append)
    ns={'__package__':'ICESEE.applications.icepack_model.examples.idealized_pig',
        '_shared_context':lambda *args:NS(Q=None),'_recompute_floating':lambda *args:None,
        'melt_forecast_mode':lambda opts:'estimated_field','estimated_melt_forecast':lambda *args:'bounded-melt',
        'IcepackForecastFields':lambda **kw:NS(**kw),
        'finalize_analysis_state':lambda *args,**kw:calls.append('solve')}
    ns['bind_melt_reference'] = bind_melt_reference
    exec(compile(ast.Module(body=[method],type_ignores=[]),str(path),'exec'),ns)
    options=dict(pig_analysis_state_policy='diagnostic_balance',num_state_vars=5,
        vec_inputs=['h','u','v','s','basal_melt_field'],process_noise_vars=['h','basal_melt_field'],sig_Q=[0,0,0,0,.1])
    ns['finalize_forecast_noise'](NS(_state=lambda member:state),NS(member_id=2),299,topology=None,icesee_kwargs=options)
    assert calls==['halo'] and updates[0].basal_melt=='bounded-melt'
    ns['finalize_forecast_noise'](NS(_state=lambda member:state),NS(member_id=2),299,topology=None,
        icesee_kwargs=dict(options,sig_Q=[.01,0,0,0,.1]))
    assert calls==['halo','solve','halo']


def test_scalar_amplitude_noise_is_applied_in_the_current_depth_basis():
    path=Path(__file__).resolve().parents[2]/'applications/icepack_model/examples/idealized_pig/_icepack_native.py'
    cls=next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.ClassDef) and n.name=='_IdealizedPigAdapter')
    method=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='apply_native_process_noise')
    ns={'np':np,'_shared_context':lambda *a:NS(Q=None),'_recompute_floating':lambda *a:None,
        'depth_law_melt':lambda *a:NS(dat=NS(data_ro=np.array([0.,.5,1.])))}
    ns['bind_melt_reference'] = bind_melt_reference
    exec(compile(ast.Module(body=[method],type_ignores=[]),str(path),'exec'),ns)
    layout=NS(iter_owned_global_intervals=lambda:iter([('h',slice(0,3),(0,3)),('basal_melt_field',slice(3,6),(3,6))]))
    written=[];member=NS(fields=NS(layout=layout),unpack_owned=written.append)
    state=NS(surface=None,thickness=None)
    increment=np.r_[np.zeros(3),np.repeat(.2,3)]
    ns['apply_native_process_noise'](NS(_state=lambda m:state),member,np.zeros(6),increment,0,
        topology=None,icesee_kwargs={'basal_melt_noise_basis':'depth_amplitude'})
    np.testing.assert_array_equal(written[0],[0,0,0,0,.1,.2])
    np.testing.assert_array_equal(increment,[0,0,0,.2,.2,.2])
