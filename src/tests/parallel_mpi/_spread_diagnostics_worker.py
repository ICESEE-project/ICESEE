"""Real ensemble-by-space RTPS and same-cycle diagnostic parity."""
import json
import sys

import numpy as np
from mpi4py import MPI

from ICESEE.src.parallelization.analysis_diagnostics import NativeAnalysisDiagnostics
from ICESEE.src.parallelization.distributed_topology import create_distributed_topology
from ICESEE.src.parallelization.distributed_fields import DistributedFieldRegistry
from ICESEE.src.parallelization.distributed_runtime import members_for_ensemble_slot
from ICESEE.src.parallelization.distributed_native_runtime import NativeDistributedMember, NativeDistributedMemberPool
from ICESEE.src.parallelization.distributed_native_cycle import NativeObservationBatch, run_native_global_analysis_cycle
from ICESEE.src.parallelization._mpi_analysis_functions import apply_analysis_controls_local

comm = MPI.COMM_WORLD
spatial = int(sys.argv[1])
size = 3 if spatial == 4 else 7
topology = create_distributed_topology(comm, spatial_ranks=spatial)
first = size*topology.spatial_rank//spatial
last = size*(topology.spatial_rank+1)//spatial
rng = np.random.default_rng(193)
initial = rng.normal(size=(size, 32))
observed = rng.normal(size=size)
sigma = np.linspace(.03, .07, size)
errors = rng.normal(size=(size, 32))*sigma[:, None]
errors -= errors.mean(axis=1, keepdims=True)
truth = np.linspace(1., 2., size)

class Field:
    name = 'state'
    global_size = size
    owned_start, owned_stop = first, last
    def __init__(self, values): self.values = values.copy()
    def read_owned(self): return self.values
    def write_owned(self, values): self.values[:] = values
    def synchronize_ghosts(self): pass

class Adapter:
    def forecast_native_member(self, member, timestep, **kwargs):
        member.unpack_owned(member.pack_owned()+.1)
    def observe_native_member(self, member, rows, **kwargs):
        return member.fields.observe_owned_rows(rows)
    def finalize_native_analysis(self, member, forecast_owned, timestep, **kwargs):
        member.unpack_owned(member.pack_owned()+.2)

ids = members_for_ensemble_slot(32, topology.ensemble_groups, topology.ensemble_slot)
pool = NativeDistributedMemberPool({mid:NativeDistributedMember(mid,
    DistributedFieldRegistry([Field(initial[first:last, mid])], layout_id='rtps-mpi')) for mid in ids})
batch = NativeObservationBatch(np.arange(first, last), observed[first:last],
    {mid: errors[first:last, mid] for mid in ids}, sigma[first:last]**2)
diagnostic = NativeAnalysisDiagnostics(pool.layout, 32, topology,
    truth[first:last] if topology.is_ensemble_root else None)
options = dict(Nens=32, nd=size, vec_inputs=['state'], num_state_vars=1,
               enkf_gain_formulation='prescribed_r', enkf_observation_error_mode='stochastic_R',
               analysis_spread_relaxation=.5, inflation_factor=1.)
result = run_native_global_analysis_cycle(pool, Adapter(), 0, [batch], number_of_batches=1,
    topology=topology, icesee_kwargs=options, error_mode='stochastic_R', state_row_chunk_size=2,
    analysis_row_controls=apply_analysis_controls_local, analysis_diagnostics=diagnostic)
forecast = initial+.1
covariance = np.cov(forecast)
gain = np.linalg.solve(covariance+np.diag(sigma**2), covariance).T
raw = forecast+gain@(observed[:, None]+errors-forecast)
sd = .5*raw.std(axis=1, ddof=1)+.5*forecast.std(axis=1, ddof=1)
controlled = raw.mean(axis=1, keepdims=True)+(raw-raw.mean(axis=1, keepdims=True))* (sd/raw.std(axis=1, ddof=1))[:, None]
balanced = controlled+.2
local_error = max((float(np.max(np.abs(result.local_analysis.members[mid]-balanced[first:last, mid])))
                   for mid in ids if last>first), default=0.)
assert max(comm.allgather(local_error)) < 3e-11
for stage, matrix in [('forecast', forecast), ('raw_analysis', raw),
                       ('controlled_analysis', controlled), ('balanced_analysis', balanced)]:
    stats = result.diagnostics['stages'][stage]['state']
    assert stats['rows'] == size
    np.testing.assert_allclose(stats['rms_member_sd'], np.sqrt(np.var(matrix, axis=1, ddof=1).mean()), atol=2e-12)
    np.testing.assert_allclose(stats['mean_rmse'], np.sqrt(np.mean((matrix.mean(axis=1)-truth)**2)), atol=2e-12)
if comm.rank == 0:
    print(json.dumps(dict(ranks=comm.size, spatial_ranks=spatial, members=32,
        zero_owned_rank=spatial==4, max_member_error=max(comm.allgather(local_error)))))
else:
    comm.allgather(local_error)
