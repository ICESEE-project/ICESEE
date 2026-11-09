import numpy as np
import pytest
from mpi4py import MPI
from ICESEE.src.parallelization.distributed_topology import create_distributed_topology
from ICESEE.src.parallelization.distributed_native_runtime import initialize_native_member_pool
from ICESEE.src.parallelization.distributed_whole_member_adapter import WholeMemberNativeAdapter,WholeMemberState
from ICESEE.src.parallelization.distributed_native_cycle import NativeObservationBatch,run_native_global_analysis_cycle
from ICESEE.src.parallelization._mpi_analysis_functions import apply_analysis_controls_local


@pytest.mark.parametrize('analyze',[False,True])
def test_cycle_applies_shared_controls_only_at_analysis(analyze):
    topology=create_distributed_topology(MPI.COMM_SELF,spatial_ranks=1)
    options=dict(Nens=2,nd=4,vec_inputs=['h','parameter'],num_state_vars=1,
                 state_inflation_factor=1.15,param_inflation_factor=1.4)
    adapter=WholeMemberNativeAdapter(
        initialize_member=lambda member_id,**kw:WholeMemberState(
            {'h':np.array([1.+member_id,3.-member_id]),
             'parameter':np.array([5.+member_id,7.-member_id])}),
        forecast_member=lambda state,timestep,**kw:None)
    pool=initialize_native_member_pool(adapter,topology,options)
    calls=[]
    def control(**kwargs):
        calls.append(kwargs['global_rows'].copy())
        return apply_analysis_controls_local(**kwargs)
    batches=[NativeObservationBatch(observation_ids=np.array([0,1]),
                                     values=np.array([1.2,2.4]))] if analyze else []
    result=run_native_global_analysis_cycle(
        pool,adapter,0,batches,number_of_batches=int(analyze),topology=topology,
        icesee_kwargs=options,error_mode='legacy_prior_anomalies',
        state_row_chunk_size=1,analysis_row_controls=control)
    forecast=np.column_stack([result.local_forecast.members[i] for i in range(2)])
    expected=forecast@result.transform
    if analyze:
        expected=apply_analysis_controls_local(analysis_vec=expected,forecast_vec=forecast,
                  global_rows=np.arange(4),icesee_kwargs=options)
        assert len(calls)==4
    else:
        assert not calls
    actual=np.column_stack([result.local_analysis.members[i] for i in range(2)])
    np.testing.assert_allclose(actual,expected,rtol=1e-14,atol=1e-14)
    np.testing.assert_allclose(pool.snapshot_owned().members[0],actual[:,0])
