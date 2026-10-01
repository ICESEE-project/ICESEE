# =============================================================================
# @author: Brian Kyanjo
# @date: 2024-11-06
# @description: Synthetic ice stream with data assimilation
# =============================================================================

# --- Imports ---
import sys
import os
import time
import numpy as np

# --- Configuration ---
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["PETSC_CONFIGURE_OPTIONS"] = "--download-mpich-device=ch3:sock"

# --- ICESEE configuration, before Firedrake: it consumes ICESEE's own
# command-line arguments, so PETSc (initialized by the Firedrake import)
# only sees PETSc options. ---
from ICESEE.config._utility_imports import icesee_kwargs

# --- firedrake imports ---
import firedrake
from firedrake.petsc import PETSc

from modelfunc import myerror
import modelfunc as mf
from modelfunc import firedrakeSmooth, flotationHeight, flotationMask

from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_model import initialize_model, initialState, initializeMesh
from ICESEE.src.run_model_da.run_models_da import icesee_model_data_assimilation
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager
from ICESEE.src.utils.performance import record_phase, register_package_versions
from ICESEE.src.utils.run_schedule import resolve_run_schedule
from mpi4py import MPI

# Record this application's model-stack versions in the end-of-run
# performance summary (every execution mode).
register_package_versions("petsc4py", "firedrake", "icepack")

# --- Register execution_mode 3 support (import-time side effect) ---
from ICESEE.applications.icepack_model.examples.idealized_pig import mode3_runner  # noqa: F401


# --- Initialize MPI ---

rank, size, comm, _ = ParallelManager().icesee_mpi_init(icesee_kwargs)

PETSc.Sys.Print("Fetching the model parameters ...")



# --- Time discretization: one authoritative model step (years) ---
# In this example `timesteps_per_year` is the step length in years, despite
# its name. `--dt`, when given on the command line, is that same step and
# takes precedence; both keys then hold the single resolved value.

num_years = float(icesee_kwargs["num_years"])
_cli = icesee_kwargs.get("cli_overrides", {})
if "dt" in _cli and "timesteps_per_year" in _cli and float(_cli["dt"]) != float(_cli["timesteps_per_year"]):
    raise ValueError(
        "--dt and --timesteps_per_year both set this example's model step and disagree: "
        f"{_cli['dt']} vs {_cli['timesteps_per_year']}"
    )
dt = float(_cli["dt"]) if "dt" in _cli else float(icesee_kwargs["timesteps_per_year"])
nt = int(round(num_years / dt))     # total number of time steps

# Physical-year lookup array (icesee/features convention, 2026-09-28
# reconciliation): t[k] gives the physical year of timestep k, used by
# run_model/generate_true_state/generate_nurged_state to match
# `save_steps` (given in physical years) to timestep indices when
# dumping flowline profiles. Does not feed BasalMeltRate (which uses
# `step`/`dt` directly) -- purely an I/O bookkeeping array.
t = np.linspace(0, int(num_years), nt + 1)

icesee_kwargs.update({"nt": nt, "dt": dt, "timesteps_per_year": dt, "t": t})

# Report the resolved schedule and initialization source, and stop now if
# the observation window gives no analysis, before any model work.
_init_file = icesee_kwargs["initFile"]
_compact = bool(icesee_kwargs.get("compact_initialization", False))
resolve_run_schedule(
    icesee_kwargs,
    rank=MPI.COMM_WORLD.Get_rank(),
    details={
        "initialization": (
            f"{_init_file} ({'compact final spin-up state' if _compact else 'full spin-up history, idx=20000'}"
            + (f", {os.path.getsize(_init_file) / 1e6:.1f} MB)" if os.path.exists(_init_file) else ", not found)")
        ),
        "truth count": 1 if icesee_kwargs.get("generate_true_state", True) else 0,
        "wrong reference": "enabled" if icesee_kwargs.get("generate_nurged_state", True) else "disabled",
    },
)


# --- Model initialization ---

PETSc.Sys.Print("Initializing icepack model ...")

icesee_kwargs.update({
    "comm": comm,
    "initFile": icesee_kwargs["initFile"],
    "paramsFile": icesee_kwargs["paramsFile"],
    "meshFile": icesee_kwargs["meshFile"],
    "SMBFile": icesee_kwargs["SMBFile"],
    "dt": dt,
    "num_years": icesee_kwargs["num_years"],
    "bmr_increase_time": int(icesee_kwargs["bmr_increase_time"]),
    "save_steps": icesee_kwargs["save_steps"],
    #"hThresh": icesee_kwargs["hThresh"]
})

_t_init = time.perf_counter()
h, h0, s, s0, u, bed, zF, grounded, floating, A0, beta0, smb, basal_melt_field, Q, V, forward_solver = initialize_model(**icesee_kwargs)
record_phase("model_initialization", time.perf_counter() - _t_init)



# ----- Update the parameters ----

icesee_kwargs["nd"] = h0.dat.data.size * icesee_kwargs["total_state_param_vars"] # get the size of the entire vector


icesee_kwargs.update({

    "smb": smb,
    "h": h,
    "h0": h0,
    "s": s,
    "s0": s0,
    "u": u,
    "basal_melt_field": basal_melt_field,
    "A0": A0,
    "beta0": beta0,
    "Q":Q,
    "V":V,
    "bed": bed,
    "seed":float(icesee_kwargs["seed"]),
    "zF": zF,
    "grounded": grounded,
    "floating": floating,
    "wrong_basal_melt_field": float(icesee_kwargs["wrong_basal_melt_field"]),
    "solver": forward_solver,
    "nd": icesee_kwargs["nd"],
    #"Lx":float(icesee_kwargs["Lx"]),
    #"Ly":float(icesee_kwargs["Ly"]),
    #"nx":float(icesee_kwargs["nx"]),
    #"ny":float(icesee_kwargs["ny"]),

})



# ----- Nudge the basal melt rate field -----

wrong_bmr = firedrake.Constant(icesee_kwargs["wrong_basal_melt_field"])

# firedrake.interpolate(expr, Q) (the free-function form) now returns a
# lazy symbolic Interpolate object in the installed Firedrake/UFL version,
# not an evaluated Function -- Function(Q).interpolate(expr) is the
# eager-evaluation equivalent (same math). Confirmed this matters here
# specifically: icesee_kwargs["basal_melt_field"] is itself produced by
# BasalMeltRate's own interpolate() call, so a lazy result fed into a
# second interpolate() call is exactly the "Non-contiguous argument
# numbers in interpolate" failure reproduced while diagnosing this.
bmr_nudged = firedrake.Function(icesee_kwargs["Q"]).interpolate(
    icesee_kwargs["basal_melt_field"] + wrong_bmr
)

icesee_kwargs.update({"bmr_nudged": bmr_nudged})



# --- Run Data Assimilation ---
PETSc.Sys.Print("Data assimilation with ICESEE ...")
icesee_model_data_assimilation(**icesee_kwargs)
