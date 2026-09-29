# ==============================================================================
# @des: Standalone worker (not a pytest file) that builds a real Icepack
# mesh/model at whatever (world_size, ranks_per_model) it is launched
# with, computes the GLOBAL physical mesh coordinates and a member's
# initial-perturbation increment (generate_initial_member_increment) --
# and, optionally, a forecast-time process-noise increment
# (add_member_process_noise's own field-generation call) -- then saves
# both, on world rank 0, to an .npz file for offline coordinate-aware
# comparison between an R=1 and an R=2 run of this same script.
#
# This is validation/test tooling only: it does not touch any production
# code path, and its coordinate-mapped comparison lives entirely outside
# the hot DA loop (Part 21's constraint).
#
# Usage: mpirun -n <P> python _icepack_coordinate_noise_worker.py \
#            <out.npz> <nx> <ny> <ranks_per_model> <base_seed> <ens_id> \
#            <random_field_method: fft|graph> <mode: init|process>
# ==============================================================================
import sys
import os
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
for _p in (str(_REPO_ROOT), str(_REPO_ROOT.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

os.environ["OMP_NUM_THREADS"] = "1"

# _utility_imports (imported below) parses sys.argv itself (both a named
# argparse parser and a generic --key=value override mechanism) -- this
# script's own positional CLI args are not meant for it. Capture them and
# truncate sys.argv to just argv[0] first, exactly as run_da_icepack.py's
# own argv is naturally just argv[0] when launched without CLI overrides.
_worker_argv = sys.argv[1:9]
sys.argv = sys.argv[:1]

_SYNTHETIC_ICE_STREAM_DIR = (
    _REPO_ROOT / "applications" / "icepack_model" / "examples" / "synthetic_ice_stream"
)
os.chdir(_SYNTHETIC_ICE_STREAM_DIR)

import numpy as np
from mpi4py import MPI

import firedrake  # noqa: F401 -- must be imported before ICESEE's Icepack pieces

from ICESEE.config._utility_imports import icesee_kwargs
from ICESEE.applications.icepack_model.examples.synthetic_ice_stream._icepack_model import (
    initialize_model,
)
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import (
    ParallelManager,
)
from ICESEE.src.parallelization.parallel_mpi.model_capabilities import (
    register_model_capabilities,
)
from ICESEE.src.utils.localization import get_mesh_coordinates
from ICESEE.src.run_model_da._error_generation import generate_initial_member_increment
from ICESEE.src.parallelization._mpi_forecast_functions import add_member_process_noise

register_model_capabilities(
    "icepack", supports_multi_rank_per_model=True, supports_distributed_state=True,
    notes="TEST-ONLY override for the coordinate/noise comparator worker.",
)

(out_path, nx, ny, ranks_per_model, base_seed, ens_id, method, mode) = _worker_argv
nx, ny, ranks_per_model, base_seed, ens_id = (
    int(nx), int(ny), int(ranks_per_model), int(base_seed), int(ens_id),
)

icesee_kwargs.update({
    "default_run": True,
    "Nens": max(ens_id + 1, 2),
    "nx": nx, "ny": ny,
    "execution_mode": 2,
    "ranks_per_model": ranks_per_model,
    "random_field_method": method,
    "base_seed": base_seed,
})

rank, size, comm, _ = ParallelManager().icesee_mpi_init(icesee_kwargs)
if rank is None:
    # Spare rank: nothing to compute or save.
    MPI.COMM_WORLD.Barrier()
    sys.exit(0)

icesee_kwargs.update({
    "comm": comm,
    "nt": int(float(icesee_kwargs["num_years"])) * int(float(icesee_kwargs["timesteps_per_year"])),
    "dt": 1.0 / float(icesee_kwargs["timesteps_per_year"]),
})
icesee_kwargs = initialize_model(**icesee_kwargs)
icesee_kwargs["nd"] = icesee_kwargs["h0"].dat.data.size * icesee_kwargs["total_state_param_vars"]
icesee_kwargs.update({
    "seed": float(icesee_kwargs["seed"]),
    "subcomm": comm,
    "color": 0,
    "sub_rank": comm.Get_rank(),
})

# Global per-member state size and hdim, matching icesee_da_full_parallel.py's
# own resolve_state_ownership-based computation.
from ICESEE.src.utils.state_ownership import resolve_state_ownership
ownership = resolve_state_ownership(icesee_kwargs, comm)
global_shape = ownership.global_size
hdim = global_shape // icesee_kwargs["total_state_param_vars"]

coords = get_mesh_coordinates(icesee_kwargs)

if mode == "init":
    increment, raw = generate_initial_member_increment(
        hdim, icesee_kwargs, ens_id, global_shape
    )
elif mode == "process":
    icesee_kwargs.update({
        "k": 1, "dt": icesee_kwargs["dt"], "alpha": 0.5, "rho": 1.0,
        "sub_rank": comm.Get_rank(),
        "process_noise_schedule": "every_step",
        "data_path": os.path.dirname(out_path) or ".",
    })
    ensemble_vec = np.zeros(global_shape, dtype=np.float64)
    updated = add_member_process_noise(ensemble_vec, ens_id, icesee_kwargs)
    increment = updated
    raw = updated
else:
    raise ValueError(f"unknown mode {mode!r}")

if comm.Get_rank() == 0:
    np.savez(
        out_path,
        coords=np.asarray(coords),
        increment=np.asarray(increment),
        raw=np.asarray(raw),
        hdim=hdim,
        total_state_param_vars=icesee_kwargs["total_state_param_vars"],
        vec_inputs=np.asarray(icesee_kwargs.get("vec_inputs", []), dtype=object),
    )

MPI.COMM_WORLD.Barrier()
