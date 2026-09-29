# ==============================================================================
# @des: Synthetic, small-file real-MPI worker directly exercising
# EnKF_fully_parallel_IO's batch-window read/write collective lifecycle
# under Stage 4A resource-plan topologies that no current real application
# can exercise in Mode 2 (ranks_per_model > 1 -- every registered
# application is conservatively supports_multi_rank_per_model=False, see
# model_capabilities.py). Bypasses the model layer entirely: deterministic
# synthetic per-member vectors stand in for a real forecast, applying the
# exact same "unconditional priming collective before the round-loop
# guard" pattern _mpi_forecast_functions.py now uses, so this test
# validates the HDF5/collective fix generalizes to R>1 + spare ranks +
# partial rounds together, not just the R=1 cases the real Lorenz96
# matrix covers.
#
# Usage: mpirun -n <P> python _enkf_io_topology_worker.py <nens> <ranks_per_model> <out_dir> <n_timesteps> <batch_size>
# ==============================================================================
import sys

import numpy as np
from mpi4py import MPI

from ICESEE.src.parallelization.parallel_mpi.resource_plan import plan_resources
from ICESEE.src.parallelization.EnKF_parallel_io import EnKF_fully_parallel_IO

nens = int(sys.argv[1])
ranks_per_model = int(sys.argv[2])
out_dir = sys.argv[3]
n_timesteps = int(sys.argv[4])
batch_size = int(sys.argv[5])

world = MPI.COMM_WORLD
world_rank = world.Get_rank()
world_size = world.Get_size()

plan = plan_resources(world_size, nens, ranks_per_model)
is_spare = plan.is_spare(world_rank)
color = None if is_spare else plan.group_id(world_rank)
key = world_rank
subcomm = world.Split(MPI.UNDEFINED if is_spare else color, key)
sub_rank = None if is_spare else subcomm.Get_rank()
subcomm_size_min = plan.num_model_groups
rounds = plan.num_rounds

nd = 6  # small, deterministic synthetic state size

icesee_kwargs = {
    "ensemble_history_mode": "full",
    "restart_enabled": False,
    "force_fresh_start": True,
    "collective_threshold": 10_000,  # force independent I/O for this tiny test
}

enkf_io = EnKF_fully_parallel_IO(
    "topology_probe_ens", nd, nens, n_timesteps, subcomm, world, icesee_kwargs,
    serial_file_creation=True, base_path=out_dir, batch_size=batch_size,
)

# Deterministic synthetic per-(member, timestep) vector: every entry
# encodes member id and timestep, so a read-back mismatch, a
# cross-member/timestep contamination, or a lost batch transition is
# immediately visible and unambiguous.
def synthetic_vector(ens_id, t):
    # +1 offset so member 0 at t=0 is never legitimately all-zero --
    # zero is reserved for "never written" in this test's own file-content
    # validation (see test_enkf_parallel_io_topology.py).
    return np.full(nd, fill_value=1000.0 * ens_id + t + 1.0, dtype=np.float64)


# Initial condition (t=0): only active ranks with a round-0 member write;
# every rank (spare included) must still join the collective priming open,
# exactly mirroring _mpi_ensemble_intialization.py's already-correct
# pattern for this same file/timestep.
enkf_io._ensure_batch(0)
if color is not None:
    # Every round's member gets its t=0 initial condition written here,
    # matching _mpi_ensemble_intialization.py's real pattern (which loops
    # every round, not just round 0) -- a group scheduled in a later
    # round still needs its own initial condition on the shared shard 0.
    for round_id in range(rounds):
        ens_id0 = plan.member_for(round_id, color)
        if ens_id0 is not None and sub_rank == 0:
            enkf_io.write_forecast(0, synthetic_vector(ens_id0, 0), ens_id0)
world.Barrier()

errors = []
for t in range(n_timesteps):
    # The Stage 4B fix: every world rank -- spare or active, with or
    # without a member this round -- joins this one priming collective for
    # timestep t before any round-dependent branching.
    enkf_io._ensure_batch_range(t, min(t + 1, enkf_io.nt - 1))

    if color is not None:
        for round_id in range(rounds):
            ens_id = color + round_id * subcomm_size_min
            if ens_id < 0 or ens_id >= nens:
                continue
            subcomm.Barrier()
            state = enkf_io.read_forecast(t, ens_id)
            expected = synthetic_vector(ens_id, t)
            if not np.allclose(state, expected):
                errors.append(
                    f"rank={world_rank} t={t} ens_id={ens_id}: "
                    f"read {state.tolist()} expected {expected.tolist()}"
                )
            if sub_rank == 0:
                enkf_io.write_forecast(t + 1, synthetic_vector(ens_id, t + 1), ens_id)
            subcomm.Barrier()

world.Barrier()
enkf_io.close()
world.Barrier()

errors = world.allgather(errors)
flat_errors = [e for rank_errors in errors for e in rank_errors]

print(
    f"RESULT rank={world_rank} is_spare={is_spare} "
    f"n_errors_seen_by_this_rank={len(errors[world_rank])} "
    f"total_errors={len(flat_errors)}",
    flush=True,
)
if world_rank == 0 and flat_errors:
    for e in flat_errors[:20]:
        print(f"ERROR {e}", flush=True)

print(f"DONE rank={world_rank}", flush=True)
