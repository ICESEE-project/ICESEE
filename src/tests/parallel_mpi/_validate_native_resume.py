"""Inspect completed native runs; does not alter checkpoints."""
from pathlib import Path
import json,sys
import h5py,numpy as np
first,second=map(Path,sys.argv[1:3])
a=first/"_mode3_state_history/steps/checkpoint_00000019"
b=second/"_mode3_state_history/steps/checkpoint_00000019"
for directory in (a,b):
 manifest=json.loads((directory/"manifest.json").read_text())
 assert manifest["complete"] and manifest["metadata"]["completed_step"]==19
 assert manifest["metadata"]["analysis_events_completed"]==2
maxdiff=0.;normdiff=0.;normbase=0.
for path in sorted(a.glob("shard_*.h5")):
 with h5py.File(path,"r") as f,h5py.File(b/path.name,"r") as g:
  def visit(name,obj):
   global maxdiff,normdiff,normbase
   if isinstance(obj,h5py.Dataset):
    x=obj[:];y=g[name][:]
    assert np.isfinite(x).all() and np.isfinite(y).all(),name
    diff=x-y
    if x.size:maxdiff=max(maxdiff,float(np.max(np.abs(diff))))
    normdiff+=float(np.sum(diff**2));normbase+=float(np.sum(x**2))
    np.testing.assert_allclose(y,x,rtol=1e-9,atol=1e-7,err_msg=name)
  f.visititems(visit)
print(f"RESUME_RESULT max_abs={maxdiff:.12g} relative_l2={np.sqrt(normdiff/normbase):.12g} finite=True")
