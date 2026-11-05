from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from ICESEE.src.parallelization.distributed_adapter import DistributedBlockStateLayout,DistributedStateBlockLayout
from ICESEE.src.parallelization.distributed_runtime import LocalMemberEnsemble,transform_local_members
from ICESEE.src.parallelization.distributed_streaming_runtime import transform_members_via_store
from ICESEE.src.parallelization.distributed_member_store import MemoryInactiveMemberStore
from ICESEE.src.parallelization._mpi_analysis_functions import apply_analysis_controls_local


def test_segmented_native_and_store_controls_match_dense_reference():
    layout=DistributedBlockStateLayout((
        DistributedStateBlockLayout("h",0,5,1,4),
        DistributedStateBlockLayout("parameter",5,5,2,5)))
    forecast=np.arange(18,dtype=float).reshape(6,3)
    transform=np.array([[.8,.1,.2],[.1,.7,.2],[.1,.2,.6]])
    options=dict(nd=10,Nens=3,vec_inputs=["h","parameter"],num_state_vars=1,
                 state_inflation_factor=1.3,param_inflation_factor=1.1,
                 analysis_relaxation_factor=.8,bed_update_policy="always")
    global_rows=layout.local_to_global_rows(np.arange(6))
    expected=apply_analysis_controls_local(
        analysis_vec=forecast@transform,forecast_vec=forecast,
        global_rows=global_rows,icesee_kwargs=options)
    def controls(rows,prior,analysis):
        mapped=layout.local_to_global_rows(np.arange(rows.start,rows.stop))
        return apply_analysis_controls_local(analysis_vec=analysis,
                    forecast_vec=prior,global_rows=mapped,icesee_kwargs=options)
    local=LocalMemberEnsemble(layout,{i:forecast[:,i].copy() for i in range(3)})
    actual=transform_local_members(local,transform,
        SimpleNamespace(ensemble_comm=MPI.COMM_SELF),options,row_chunk_size=2,
        row_transform=controls)
    np.testing.assert_allclose(np.column_stack(list(actual.members.values())),expected,rtol=1e-14,atol=1e-14)
    store=MemoryInactiveMemberStore()
    for i in range(3):store.put(i,forecast[:,i].copy())
    transform_members_via_store(store,(0,1,2),6,transform,MPI.COMM_SELF,3,
                                 row_chunk_size=2,row_transform=controls)
    np.testing.assert_allclose(np.column_stack([store.get(i) for i in range(3)]),
                               expected,rtol=1e-14,atol=1e-14)
