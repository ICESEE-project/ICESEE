"""MPI output-only regression: real PIG dimensions, no model solve or noise."""
import sys,time
from pathlib import Path
import h5py,numpy as np
from mpi4py import MPI
from ICESEE.src.parallelization.EnKF_parallel_io import EnKF_fully_parallel_IO
world=MPI.COMM_WORLD
root=Path(sys.argv[1]); root.mkdir(parents=True,exist_ok=True)
world.Barrier()
nd,nens=270615,32
io=EnKF_fully_parallel_IO('shape_probe',nd,nens,2,world,world,
 {'ensemble_history_mode':'full','force_fresh_start':True},
 serial_file_creation=True,base_path=str(root),batch_size=1)
start=time.monotonic()
path,handle,dset=io._create_rank_analysis_file(1)
first,last=io.nd_start_world,io.nd_end_world
for row in range(first,last,1000):
 stop=min(last,row+1000)
 dset[row-first:stop-first,:]=np.arange(row,stop)[:,None]+np.arange(nens)[None,:]*0.125
handle.close()
temp=io._assemble_rank_analysis_files(1,path)
io._publish_analysis_file(1,temp,{})
io._write_published_analysis_mean(1)
errors=0
with h5py.File(io._state_file_name(1),'r',driver='mpio',comm=world) as f:
 for row in range(first,last,1000):
  stop=min(last,row+1000)
  expected=np.arange(row,stop)[:,None]+np.arange(nens)[None,:]*0.125
  errors+=int(not np.array_equal(f['states'][row:stop,:],expected))
with h5py.File(root/'shape_probe_mean.h5','r',driver='mpio',comm=world) as f:
 expected=np.arange(first,last)+(nens-1)*0.125/2
 errors+=int(not np.array_equal(f['mean'][first:last,1],expected))
io.close()
total=world.allreduce(errors)
print(f'WRITE_RESULT rank={world.rank} errors={errors} total_errors={total} elapsed={time.monotonic()-start:.3f}s',flush=True)
assert total==0
