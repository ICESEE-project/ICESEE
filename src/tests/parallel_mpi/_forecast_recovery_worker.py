"""One-rank rejection must roll every spatial participant back together."""
from mpi4py import MPI
import numpy as np
from ICESEE.src.parallelization.forecast_recovery import (
    ForecastRecoveryPolicy, RecoverableForecastError, advance_with_recovery,
)
from ICESEE.src.parallelization.run_control import CycleStopController

comm = MPI.COMM_WORLD
owned = np.array([float(comm.rank)])
def advance(time,dt):
    owned[:] += dt
    if comm.rank == comm.size-1 and dt > .25:
        raise RecoverableForecastError("one rank rejected")
result = advance_with_recovery(start_time=14,interval=1,
    policy=ForecastRecoveryPolicy(3,.0625),capture=lambda:owned.copy(),
    restore=lambda saved:owned.__setitem__(slice(None),saved),advance=advance,
    validate=lambda:None,comm=comm,stable_dt=.75 if comm.rank==0 else .9)
assert result["substeps"] == 4
assert owned[0] == comm.rank+1
assert len(set(comm.allgather((result["attempts"],result["dt"])))) == 1
def broken_capture():
    if comm.rank == comm.size-1:
        raise OSError("snapshot unavailable")
    return owned.copy()
try:
    advance_with_recovery(start_time=15,interval=1,
        policy=ForecastRecoveryPolicy(2,.125),capture=broken_capture,
        restore=lambda saved:None,advance=lambda *a:None,validate=lambda:None,comm=comm)
except RuntimeError as exc:
    assert "forecast capture failed" in str(exc)
else:
    raise AssertionError("rank-local capture failure was ignored")
assert owned[0] == comm.rank+1
with CycleStopController({"run_segment_cycles":3},comm) as controller:
    if comm.rank == comm.size-1:
        controller.reason = "SIGUSR1"
    assert controller.after_commit(1)
    assert controller.reason == "SIGUSR1"
if comm.rank == 0:
    print(f"PASS: {comm.size} ranks; rollback, adaptive subdivision, CFL minimum, capture failure, collective stop",flush=True)
