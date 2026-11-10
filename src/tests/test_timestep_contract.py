"""State-column timing, initial-state preservation and restart contract tests."""
import ast
import json
from pathlib import Path
import types
import h5py
import numpy as np
import pytest
from ICESEE.src.utils.timestep_contract import (
    TIME_GRID_CONTRACT, observation_cycle_steps, observation_is_due,
    completed_observations, prepare_timestep_contract,
)


@pytest.mark.parametrize("dt", [.1,.05,.025,.06])
def test_observations_at_first_interior_and_final_saved_states(dt):
    nt = round(1/dt)
    columns = np.array([1, nt//2, nt])
    applied = [k+1 for k in range(nt)
               if observation_is_due(k, completed_observations(columns,k), columns)]
    np.testing.assert_array_equal(applied, columns)
    np.testing.assert_array_equal(np.array(applied)*dt, columns*dt)
    np.testing.assert_array_equal(observation_cycle_steps(columns,nt),columns-1)


def test_resume_before_and_after_analysis_never_repeats_or_skips_observations():
    columns = np.array([1,10,15,20])
    for start in [0,1,9,10,14,15,19,20]:
        km = completed_observations(columns,start)
        applied=[]
        for k in range(start,20):
            if observation_is_due(k,km,columns):
                applied.append(k+1);km+=1
        assert applied == columns[columns>start].tolist()


@pytest.mark.parametrize("indices", [[0,1], [-1,2], [2,2], [3,1], [.5,2]])
def test_unsupported_or_invalid_schedule_is_not_applied_after_wrong_forecast(indices):
    with pytest.raises(ValueError):observation_cycle_steps(indices,3)


def test_empty_observations_and_out_of_range_endpoint():
    assert not observation_is_due(0,0,np.array([],dtype=int))
    with pytest.raises(ValueError):observation_cycle_steps(np.array([4]),3)


class Comm:
    def __init__(self, rank=0, payload=None):self.rank,self.payload=rank,payload
    def Get_rank(self):return self.rank
    def bcast(self,value,root=0):
        if self.rank==0:self.payload=value
        return self.payload


def test_new_contract_is_atomic_reusable_and_rejects_old_state_collectively(tmp_path):
    comm=Comm();prepare_timestep_contract(comm,tmp_path)
    marker=tmp_path/'time_grid_contract.json'
    assert json.loads(marker.read_text())['time_grid_contract']==TIME_GRID_CONTRACT
    before=marker.read_bytes();prepare_timestep_contract(comm,tmp_path)
    assert marker.read_bytes()==before
    old=tmp_path/'old';old.mkdir();state=old/'icesee_ensemble_data.h5';state.write_bytes(b'preserve')
    with pytest.raises(ValueError,match='existing state'):prepare_timestep_contract(comm,old)
    with pytest.raises(ValueError,match='existing state'):prepare_timestep_contract(Comm(1,comm.payload),old)
    assert state.read_bytes()==b'preserve'
    assert not (old/'time_grid_contract.json').exists()


def test_generic_pig_initialization_copies_physical_fields_without_forecasting(tmp_path):
    import os
    source=Path(__file__).resolve().parents[2]/'applications/icepack_model/examples/idealized_pig/_icepack_enkf.py'
    node=next(n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='initialize_ensemble')
    class Field:
        def __init__(self,a):self.dat=types.SimpleNamespace(data_ro=np.array(a,dtype=float))
        def copy(self,deepcopy=False):return Field(self.dat.data_ro)
    h=Field([30,40,50]);u=Field([[1,2],[3,4],[5,6]]);surface=Field([10,20,30]);melt=Field([0,10,20])
    def forbidden(*a,**k):raise AssertionError('initialization advanced the model')
    ns=dict(np=np,h5py=h5py,os=os,Icepack=forbidden,
            initial_flowline_profile=lambda kw:(h.dat.data_ro,surface.dat.data_ro,np.zeros((3,2)),np.arange(3),np.zeros(3)))
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),ns)
    result=ns['initialize_ensemble'](0,joint_estimation=False,h0=h,u=u,s0=surface,
        basal_melt_field=melt,save_steps=[.5,1.],data_path=str(tmp_path))
    for name,expected in [('h',h.dat.data_ro),('s',surface.dat.data_ro),('u',u.dat.data_ro[:,0]),('v',u.dat.data_ro[:,1])]:
        np.testing.assert_array_equal(result[name],expected)
        assert not np.shares_memory(result[name],expected)
    with h5py.File(tmp_path/'hs_ensemble_profiles0.h5','r') as f:
        np.testing.assert_array_equal(f['written_steps'][:],[0,-1,-1])
        assert np.isnan(f['h_profiles'][:,1:]).all()


def test_restart_cannot_reinterpret_states_on_a_changed_dt_grid(tmp_path):
    comm=Comm();prepare_timestep_contract(comm,tmp_path,time_grid=np.arange(4)*.05)
    before=(tmp_path/'time_grid_contract.json').read_bytes()
    with pytest.raises(ValueError,match='incompatible'):
        prepare_timestep_contract(comm,tmp_path,time_grid=np.arange(4)*.1)
    assert (tmp_path/'time_grid_contract.json').read_bytes()==before


def test_analysis_audit_records_completed_state_and_rejects_misalignment():
    from ICESEE.src.utils.timestep_contract import format_analysis_event
    kw=dict(t=np.arange(21)*.05,obs_index=np.array([10,15,20]))
    record=json.loads(format_analysis_event(19,2,kw).split(' ',1)[1])
    assert record==dict(cycle=19,state_column=20,model_time=1.,observation_column=2,observation_state_column=20)
    with pytest.raises(ValueError,match='disagree'):format_analysis_event(18,2,kw)


def test_process_noise_and_analysis_use_the_same_completed_state_column():
    source=Path(__file__).resolve().parents[1]/'parallelization/_mpi_forecast_functions.py'
    node=next(n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='process_noise_is_due')
    namespace=dict(np=np,observation_cycle_steps=observation_cycle_steps)
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),namespace)
    due=namespace['process_noise_is_due']
    settings=dict(obs_index=[10,15,20],number_obs_instants=3)
    assert [k for k in range(20) if due(settings,k)]==[9,14,19]
    assert due(dict(settings,k=9))
    assert due(dict(process_noise_schedule='every_step'),0)
    assert not due(dict(process_noise_schedule='none'),9)
    with pytest.raises(ValueError):due(dict(obs_index=[0]),0)
