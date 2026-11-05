"""Real-MPI ownership and deterministic restart check; no model solves."""
from dataclasses import dataclass
from types import SimpleNamespace
import tempfile
import numpy as np
from mpi4py import MPI
from ICESEE.src.parallelization.native_process_noise import NativeProcessNoiseAdapter
from ICESEE.src.parallelization.distributed_fields import DistributedFieldRegistry
from ICESEE.src.parallelization.distributed_native_runtime import NativeDistributedMember
from ICESEE.src.parallelization import _mpi_forecast_functions as canonical

world=MPI.COMM_WORLD
rank,size=world.rank,world.size
n=7
start=n*rank//size;stop=n*(rank+1)//size
perm=np.array([4,0,6,2,1,5,3])
rows=np.concatenate((perm[start:stop],n+perm[start:stop]))
kw=dict(nd=2*n,nt=6,dt=.05,num_state_vars=1,total_state_param_vars=2,
        joint_estimation=True,localization_flag=False,Lx=10.,Ly=2.,
        sig_Q=[.8,0.],base_seed=42,obs_index=np.array([1,3,5]),
        number_obs_instants=3)
canonical.generate_enkf_field=lambda **options:options["rng"].standard_normal(options["noise_dim"])
topology=SimpleNamespace(spatial_comm=world)
base=SimpleNamespace(forecast_native_member=lambda *args,**options:None)
@dataclass
class Field:
 name:str
 global_size:int
 owned_start:int
 owned_stop:int
 values:np.ndarray
 def read_owned(self):return self.values
 def write_owned(self,values):self.values[:]=values
 def synchronize_ghosts(self):pass

def member(values=None):
 if values is None:values=np.zeros(2*(stop-start))
 fields=[Field(name,n,start,stop,values[i*(stop-start):(i+1)*(stop-start)].copy())
         for i,name in enumerate(("h","parameter"))]
 return NativeDistributedMember(9,DistributedFieldRegistry(fields,layout_id="noise-worker"))

def policy():
 return NativeProcessNoiseAdapter(base,topology=topology,canonical_rows=rows,
             reference_coords=np.arange(n)[:,None],icesee_kwargs=kw)

errors=0
with tempfile.TemporaryDirectory() as scratch:
 wrapped=policy();physical=member()
 expected=np.zeros(2*n)
 options=dict(wrapped.options,data_path=scratch,_modelrun_datasets=scratch)
 for k in (1,3):
  options["k"]=k
  expected=canonical.add_member_process_noise(expected,9,options)
  wrapped.forecast_native_member(physical,k,topology=topology,icesee_kwargs=kw)
  errors+=not np.array_equal(physical.pack_owned(),expected[rows])
 saved=physical.pack_owned().copy()
 options["k"]=5
 expected=canonical.add_member_process_noise(expected,9,options)
 restarted=policy();restored=member(saved)
 restarted.forecast_native_member(restored,5,topology=topology,icesee_kwargs=kw)
 errors+=not np.array_equal(restored.pack_owned(),expected[rows])
 wrapped.close();restarted.close()
total=world.allreduce(int(errors))
print(f"NOISE_RESULT rank={rank} errors={errors} total_errors={total}",flush=True)
raise SystemExit(1 if total else 0)
