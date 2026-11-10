"""Inspect closed PIG outputs for the state-column timing contract."""
import json
from pathlib import Path
import sys
import time
import h5py
import numpy as np

def main():
    run=Path(sys.argv[1]); years=float(sys.argv[2]);dt=float(sys.argv[3]);nens=int(sys.argv[4])
    steps=round(years/dt)
    obs_years=[.5,.75,1.] if years==1 else [2.,3.,4.]
    columns=[round(v/dt) for v in obs_years]
    report=dict(run=str(run),nt=steps,nens=nens,modes={},initial_velocity_abs_tolerance=1e-8)
    for mode in (list(map(int,sys.argv[5].split(','))) if len(sys.argv)>5 else [1,2,3]):
        case=run/f'mode{mode}';data=case/'data'
        assert (case/'exit-status').read_text().strip()=='0',f'mode {mode} failed'
        contract=json.loads((data/'time_grid_contract.json').read_text())
        assert contract['time_grid_contract']=='state-columns-v2'
        events=[]
        for line in (case/'stdout.log').read_text().splitlines():
            if line.startswith('[ICESEE_ANALYSIS] '):events.append(json.loads(line.split(' ',1)[1]))
        assert [e['state_column'] for e in events]==columns,(mode,events)
        assert [e['observation_state_column'] for e in events]==columns
        assert [e['observation_column'] for e in events]==[0,1,2]
        np.testing.assert_allclose([e['model_time'] for e in events],obs_years,rtol=0,atol=1e-12)
        stats=dict(events=events,fields={},profiles=[])
        with h5py.File(data/'icesee_ensemble_data.h5','r') as f:
            ensemble_dataset=f['ensemble']
            mean_dataset=f['ensemble_mean']
            nd,n,times=ensemble_dataset.shape
            assert n==nens and times==steps+1
            assert mean_dataset.shape==(nd,times)
            if 't' in f:np.testing.assert_allclose(f['t'][:],np.arange(times)*dt,rtol=0,atol=1e-12)
            if mode==3:
                assert f.attrs['complete'] and f.attrs.get('export_complete',True)
                np.testing.assert_array_equal(f['analysis_applied'][:].nonzero()[0],columns)
                blocks=json.loads(f.attrs['blocks'])
            else:
                names=['h','u','v','s','basal_melt_field'];assert nd%5==0
                blocks=[dict(name=name,global_offset=i*(nd//5),global_size=nd//5) for i,name in enumerate(names)]
            with h5py.File(data/'true_nurged_states.h5','r') as truth:
                stats['initial_physical_fields']={}
                for b in blocks:
                    name,first,size=b['name'],b['global_offset'],b['global_size']
                    sl=slice(first,first+size)
                    initial=ensemble_dataset[sl,:,0]
                    if name in ['h','u','v','s']:
                        target=truth['true_state'][sl,0]
                        difference=float(np.abs(initial-target[:,None]).max())
                        stats['initial_physical_fields'][name]=dict(max_abs_from_truth=difference,max_member_range=float(np.ptp(initial,axis=1).max()))
                        if mode==3:
                            if name in ['h','s']:np.testing.assert_array_equal(initial,np.broadcast_to(target[:,None],initial.shape))
                            else:assert difference<=1e-8,(mode,name,difference)
            for b in blocks:
                name,first,size=b['name'],b['global_offset'],b['global_size']
                field=dict(min=float('inf'),max=-float('inf'),final_min=None,final_max=None)
                for start in range(first,first+size,8192):
                    stop=min(start+8192,first+size)
                    values=ensemble_dataset[start:stop,:,:]
                    assert np.isfinite(values).all(),(mode,name,start)
                    mean=mean_dataset[start:stop,:]
                    assert np.isfinite(mean).all()
                    np.testing.assert_allclose(mean,values.mean(axis=1),rtol=1e-12,atol=1e-10)
                    field['min']=min(field['min'],float(values.min()))
                    field['max']=max(field['max'],float(values.max()))
                final=ensemble_dataset[first:first+size,:,-1]
                field.update(final_min=float(final.min()),final_max=float(final.max()))
                stats['fields'][name]=field
        paths=[data/'hs_profiles_true.h5',*[data/f'hs_ensemble_profiles{i}.h5' for i in range(nens)]]
        for p in paths:
            with h5py.File(p,'r') as f:
                assert f['h_profiles'].shape==f['s_profiles'].shape==(415,4)
                assert np.isfinite(f['h_profiles'][:]).all() and np.isfinite(f['s_profiles'][:]).all()
                np.testing.assert_array_equal(f['written_steps'][:],[0,*columns])
                np.testing.assert_allclose(f['profile_times'][:],[0,*obs_years],rtol=0,atol=1e-12)
            stats['profiles'].append(p.name)
        if mode==3:
            manifests=sorted((data/'_mode3_state_history/steps').glob('checkpoint_*/manifest.json'))
            assert len(manifests)==steps
            initial=json.loads((data/'_mode3_state_history/initial_condition/checkpoint_00000000/manifest.json').read_text())
            assert initial['metadata']['model_time']==0 and initial['metadata']['completed_step']==-1
            assert ':time-grid-v2:' in initial['run_id']
            actual=[]
            for p in manifests:
                m=json.loads(p.read_text());md=m['metadata'];k=md['completed_step']
                assert m['complete'] and md['next_step']==k+1
                assert abs(md['model_time']-(k+1)*dt)<1e-12
                if md['did_analysis']:actual.append(md['next_step'])
            assert actual==columns
            assert md['analysis_events_completed']==3 and md['next_step']==steps
        report['modes'][str(mode)]=stats
        print(f'Mode {mode}: initial fields, analyses {obs_years}, {times} saved states, profiles and finite fields checked',flush=True)
    report['timing_checks_passed']=True
    (run/('timing-verification.json' if len(sys.argv)<=5 else 'timing-selected-modes.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({mode:record['fields'] for mode,record in report['modes'].items()},indent=2),flush=True)


if __name__ == '__main__':
    main()
