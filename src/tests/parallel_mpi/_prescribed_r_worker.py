"""MPI row ownership and prescribed-R gain parity; no model solve or files."""
import json
import numpy as np
from mpi4py import MPI
from ICESEE.src.parallelization.distributed_analysis import StochasticAnalysisProducts, ensemble_transform_from_products

comm=MPI.COMM_WORLD
rng=np.random.default_rng(73)
x=rng.normal(size=(7,32));h=rng.normal(size=(91,7));d=rng.normal(size=91)
sigma=np.linspace(.2,1.3,91)
eta=rng.normal(size=(91,32))*sigma[:,None];eta-=eta.mean(axis=1,keepdims=True)
y=h@x
products=StochasticAnalysisProducts.zeros(32,gain_formulation='prescribed_r')
owned=np.arange(91)[comm.rank::comm.size][::-1]
for start in range(0,len(owned),7):
    rows=owned[start:start+7]
    products.add_chunk(y[rows],d[rows],error_mode='stochastic_R',
        observation_errors=eta[rows],observation_variances=sigma[rows]**2)
products=products.reduced(comm)
result=x@ensemble_transform_from_products(products)
p=np.cov(x);gain=np.linalg.solve(h@p@h.T+np.diag(sigma**2),h@p).T
expected=x+gain@(d[:,None]+eta-y)
np.testing.assert_allclose(result,expected,rtol=2e-11,atol=2e-11)
assert products.observation_rows==91
if comm.rank==0:print(json.dumps(dict(ranks=comm.size,members=32,observations=91,max_member_error=float(np.max(np.abs(result-expected))))),flush=True)

from types import SimpleNamespace
from ICESEE.src.parallelization.distributed_native_cycle import _collective_native_call
def fails_on_one_rank():
    if comm.rank==0:raise ValueError('projection limit')
try:
    _collective_native_call(SimpleNamespace(world=comm),fails_on_one_rank,'analysis finalization')
except RuntimeError as exc:
    assert 'projection limit' in str(exc)
else:
    raise AssertionError('every rank must receive the model-group failure')
if comm.rank==0:print('Collective failure propagation passed',flush=True)

# Exercise the native cycle on a genuine ensemble-by-space process grid.
from ICESEE.src.parallelization.distributed_topology import create_distributed_topology
from ICESEE.src.parallelization.distributed_fields import DistributedFieldRegistry
from ICESEE.src.parallelization.distributed_native_runtime import NativeDistributedMember,NativeDistributedMemberPool
from ICESEE.src.parallelization.distributed_native_cycle import NativeObservationBatch,run_native_global_analysis_cycle
from ICESEE.src.parallelization.distributed_runtime import members_for_ensemble_slot

topology=create_distributed_topology(comm,spatial_ranks=2 if comm.size%2==0 else 1)
size=6;first=size*topology.spatial_rank//topology.spatial_ranks
last=size*(topology.spatial_rank+1)//topology.spatial_ranks
rng=np.random.default_rng(91);initial=rng.normal(size=(size,32));d=rng.normal(size=size)
sigma=np.linspace(.2,.7,size);eta=rng.normal(size=(size,32))*sigma[:,None];eta-=eta.mean(axis=1,keepdims=True)
class Field:
    name='state';global_size=size;owned_start=first;owned_stop=last
    def __init__(self,values):self.values=values.copy()
    def read_owned(self):return self.values
    def write_owned(self,values):self.values[:]=values
    def synchronize_ghosts(self):pass
class Adapter:
    def forecast_native_member(self,member,timestep,**kwargs):member.unpack_owned(member.pack_owned()+.1)
    def observe_native_member(self,member,rows,**kwargs):return member.fields.observe_owned_rows(rows)
    def finalize_native_analysis(self,*args,**kwargs):pass
ids=members_for_ensemble_slot(32,topology.ensemble_groups,topology.ensemble_slot)
pool=NativeDistributedMemberPool({mid:NativeDistributedMember(mid,
    DistributedFieldRegistry([Field(initial[first:last,mid])],layout_id='prescribed-r-mpi')) for mid in ids})
batch=NativeObservationBatch(np.arange(first,last),d[first:last],
    {mid:eta[first:last,mid] for mid in ids},sigma[first:last]**2)
result=run_native_global_analysis_cycle(pool,Adapter(),0,[batch],number_of_batches=1,topology=topology,
    icesee_kwargs=dict(Nens=32,enkf_gain_formulation='prescribed_r',enkf_observation_error_mode='stochastic_R'),
    error_mode='stochastic_R',state_row_chunk_size=2)
forecast=initial+.1;p=np.cov(forecast);gain=np.linalg.solve(p+np.diag(sigma**2),p).T
expected=forecast+gain@(d[:,None]+eta-forecast)
error=max(float(np.max(np.abs(result.local_analysis.members[mid]-expected[first:last,mid]))) for mid in ids)
errors=comm.allgather(error)
assert max(errors)<2e-11
if comm.rank==0:print(json.dumps(dict(native_cycle_ranks=comm.size,spatial_ranks=topology.spatial_ranks,
    ensemble_groups=topology.ensemble_groups,max_member_error=max(errors))),flush=True)
