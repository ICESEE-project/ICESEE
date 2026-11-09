# ==============================================================================
# @des: Real-MPI worker for the mode-2 analysis read path in
# EnKF_fully_parallel_IO.  On NFS, MPI-IO reads through read-write handles
# can stall indefinitely, so analysis reads must go through a read-only view
# of the forecast shard while forecast writes keep their read-write window.
#
# One run covers, on every rank:
#   * analysis reads: only read-only handles on the analysed shard while the
#     transform is computed, and a published analysis equal to a serial
#     reference built from the same functions;
#   * subsequent writes: the read-write window returns lazily and forecast
#     writes after the analysis read back exactly;
#   * batch-window transitions: batch_size=2, so forecasts and analyses
#     cross window boundaries;
#   * checkpoint/resume: a new I/O object resumes from the driver's JSON
#     checkpoint, sees the published analysis, and runs a second analysis.
#
# Usage: mpirun -n <P> python _enkf_io_readonly_analysis_worker.py <out_dir>
# ==============================================================================
import sys
import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(90, repeat=True)

import h5py
import numpy as np
from mpi4py import MPI

from ICESEE.src.parallelization.EnKF_parallel_io import (
    EnKF_fully_parallel_IO,
    ensemble_transform_from_products,
)
from ICESEE.src.parallelization._mpi_analysis_functions import (
    apply_analysis_controls_local,
)
from ICESEE.src.utils.tools import load_checkpoint, save_checkpoint

out_dir = sys.argv[1]
world = MPI.COMM_WORLD
rank = world.Get_rank()
size = world.Get_size()
subcomm = world.Split(rank, 0)

PREFIX = "readonly_probe_ens"
nd, nens, nt_model, batch_size = 12, 4, 7, 2
obs_rows = np.array([0, 2, 5, 7, 8, 11], dtype=np.int64)
obs_steps = [3, 6]          # analysed shards; the first is before the restart
restart_after = 4           # last shard written before the simulated restart

errors = []


def check(name, condition, detail=""):
    if not condition:
        errors.append(f"rank={rank} {name}: {detail}")


def base_kwargs(fresh):
    return {
        "ensemble_history_mode": "full",
        "restart_enabled": True,
        "force_fresh_start": fresh,
        "collective_threshold": 10_000,
        "nd": nd,
        "vec_inputs": ["h", "u"],
        "num_state_vars": 2,
        "enkf_observation_error_mode": "legacy_prior_anomalies",
        "compute_analysis_mean": True,
        "model_name": "synthetic",
    }


def forecast_value(previous, ens_id, t):
    """Deterministic member-dependent forecast increment."""
    rng = np.random.default_rng(1000 * ens_id + t)
    return previous + 0.5 + rng.standard_normal(nd)


def open_state_intents():
    """Intents of every open HDF5 file id on a state shard, by file name."""
    found = []
    for fid in h5py.h5f.get_obj_ids(types=h5py.h5f.OBJ_FILE):
        name = fid.name.decode() if isinstance(fid.name, bytes) else str(fid.name)
        if PREFIX in name and not name.endswith("_mean.h5"):
            readonly = fid.get_intent() == h5py.h5f.ACC_RDONLY
            found.append((name.rsplit("/", 1)[-1], readonly))
    return found


class ProbeIO(EnKF_fully_parallel_IO):
    """Records open handles and a serial reference at the transform step."""

    def compute_X5_modified(self, state_dset=None, **icesee_kwargs):
        k = icesee_kwargs["k"] + 1
        intents = open_state_intents()
        shard = f"{PREFIX}_{k:04d}.h5"
        check("analysis handles are read-only", all(ro for _, ro in intents), intents)
        check("only the analysed shard is open",
              sorted({name for name, _ in intents}) == [shard], intents)
        check("no stale read-write window",
              self.files == [] and self.datasets == [] and self.current_batch_start == -1,
              (len(self.files), len(self.datasets), self.current_batch_start))
        check("analysis dataset is a read-only view",
              state_dset is not None and state_dset.file.mode == "r",
              None if state_dset is None else state_dset.file.mode)
        forecast = np.asarray(state_dset[:, :], dtype=np.float64)
        observations = icesee_kwargs["_probe_obs"]
        y = forecast[obs_rows]
        yprime = y - y.mean(axis=1, keepdims=True)
        factor = 2.0 * yprime
        X5_ref = ensemble_transform_from_products(
            yprime.T @ factor, factor.T @ factor,
            factor.T @ (observations[:, None] - y), 0.999,
        )
        self.reference = apply_analysis_controls_local(
            analysis_vec=forecast @ X5_ref,
            forecast_vec=forecast,
            global_rows=np.arange(nd, dtype=np.int64),
            icesee_kwargs=dict(icesee_kwargs),
        )
        return super().compute_X5_modified(state_dset=state_dset, **icesee_kwargs)


def make_io(fresh):
    return ProbeIO(PREFIX, nd, nens, nt_model, subcomm, world, base_kwargs(fresh),
                   serial_file_creation=True, base_path=out_dir, batch_size=batch_size)


def my_members():
    return [e for e in range(nens) if e % size == rank]


def forecast_step(io, t, expected):
    """Every rank reads t (collective window), owners write t + 1."""
    if rank == 0:
        print(f"PROGRESS forecast {t}", flush=True)
    io._ensure_batch_range(t, min(t + 1, io.nt - 1))
    for ens_id in my_members():
        state = io.read_forecast(t, ens_id)
        if t in expected:
            check(f"read t={t} member={ens_id}",
                  np.allclose(state, expected[t][:, ens_id]),
                  (state, expected[t][:, ens_id]))
        nxt = forecast_value(state, ens_id, t + 1)
        io.write_forecast(t + 1, nxt, ens_id)
        check(f"write t={t + 1} member={ens_id} restores write access",
              any(name == f"{PREFIX}_{t + 1:04d}.h5" and not ro
                  for name, ro in open_state_intents()),
              open_state_intents())
        # read_forecast anticipates t+2 and may collectively reopen files.
        # This owner-only assertion must use the already-open output shard.
        io._require_batch(t + 1)
        back = np.asarray(io.datasets[t + 1 - io.current_batch_start][:, ens_id])
        check(f"write-back t={t + 1} member={ens_id}", np.allclose(back, nxt))
    world.Barrier()


def analyse(io, shard, km, observations):
    if rank == 0:
        print(f"PROGRESS analysis {shard}", flush=True)
    settings = base_kwargs(False)
    settings.update({
        "k": shard - 1, "km": km,
        "time_analysis_file_writing": 0.0,
        "time_analysis_ensemble_mean_generation": 0.0,
        "_probe_obs": observations,
    })
    io.compute_analysis_update(**settings)
    check("window released after analysis",
          io.files == [] and open_state_intents() == [],
          (len(io.files), open_state_intents()))
    with h5py.File(f"{out_dir}/{PREFIX}_{shard:04d}.h5", "r", driver="mpio", comm=world) as f:
        published = np.asarray(f["states"][:, :])
    check(f"published analysis {shard} matches serial reference",
          np.allclose(published, io.reference, rtol=1e-12, atol=1e-12),
          np.max(np.abs(published - io.reference)))
    with h5py.File(f"{out_dir}/{PREFIX}_mean.h5", "r", driver="mpio", comm=world) as f:
        mean = np.asarray(f["mean"][:, shard])
    check(f"analysis mean {shard}", np.allclose(mean, published.mean(axis=1)))
    world.Barrier()
    return published


# Truth-free synthetic observations, identical on every rank.
obs_all = {s: np.linspace(10.0, 20.0, obs_rows.size) + s for s in obs_steps}
if rank == 0:
    with h5py.File(f"{out_dir}/H_matrix.h5", "w") as f:
        f["obs_indices"] = obs_rows
    with h5py.File(f"{out_dir}/synthetic_obs.h5", "w") as f:
        hu = np.zeros((nd, len(obs_steps)))
        for km, s in enumerate(obs_steps):
            hu[obs_rows, km] = obs_all[s]
        f["hu_obs"] = hu
world.Barrier()

# ---- Phase 1: fresh run with one analysis -----------------------------------
io = make_io(fresh=True)
io._ensure_batch(0)
for ens_id in my_members():
    io.write_forecast(0, np.arange(nd, dtype=np.float64) * (1.0 + 0.1 * ens_id), ens_id)
world.Barrier()

expected = {}
for t in range(restart_after):
    forecast_step(io, t, expected)
    if t + 1 in obs_steps:
        expected[t + 1] = analyse(io, t + 1, obs_steps.index(t + 1), obs_all[t + 1])

# Driver-style checkpoint after the last completed step, then shut down.
last_done = restart_after - 1
if rank == 0:
    save_checkpoint(out_dir, last_done_k=last_done, km=1, nt=nt_model, nd=nd,
                    nens=nens, ensemble_state_t=restart_after)
io._ensure_batch_range(restart_after, min(restart_after + 1, io.nt - 1))
pre_restart = {ens_id: io.read_forecast(restart_after, ens_id) for ens_id in my_members()}
world.Barrier()
io.close()
world.Barrier()

# ---- Phase 2: resume from the checkpoint ------------------------------------
ckpt = world.bcast(load_checkpoint(out_dir) if rank == 0 else None, root=0)
check("checkpoint loaded", ckpt is not None and ckpt["last_done_k"] == last_done, ckpt)
k_start = int(ckpt["last_done_k"]) + 1
io = make_io(fresh=False)
io._ensure_batch(obs_steps[0])
for ens_id in my_members():
    check("resumed run sees published analysis",
          np.allclose(io.read_forecast(obs_steps[0], ens_id),
                      expected[obs_steps[0]][:, ens_id]))
# All ranks, including ranks without members, prepare the resumed read/write
# window before any member owner calls read_forecast.
io._ensure_batch_range(k_start, min(k_start + 1, io.nt - 1))
for ens_id, state in pre_restart.items():
    check("resumed run sees pre-restart forecast",
          np.allclose(io.read_forecast(k_start, ens_id), state))
world.Barrier()

for t in range(k_start, nt_model):
    forecast_step(io, t, expected)
    if t + 1 in obs_steps:
        expected[t + 1] = analyse(io, t + 1, obs_steps.index(t + 1), obs_all[t + 1])

io.close()
world.Barrier()

faulthandler.cancel_dump_traceback_later()
all_errors = [e for rank_errors in world.allgather(errors) for e in rank_errors]
print(f"RESULT rank={rank} errors={len(errors)} total_errors={len(all_errors)}", flush=True)
if rank == 0:
    for e in all_errors[:20]:
        print(f"ERROR {e}", flush=True)
