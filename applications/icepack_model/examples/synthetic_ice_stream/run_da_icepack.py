# =============================================================================
# @author: Brian Kyanjo
# @date: 2024-11-06
# @description: Synthetic ice stream with data assimilation
# =============================================================================

# --- Imports ---
import sys
import os
import numpy as np
from pathlib import Path

# --- Set up paths ---
os.chdir(Path(__file__).resolve().parent)

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

from ICESEE.applications.icepack_model.examples.synthetic_ice_stream._icepack_model import initialize_model
from ICESEE.src.run_model_da.run_models_da import icesee_model_data_assimilation
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager
from ICESEE.src.utils.performance import register_package_versions

# Record this application's model-stack versions in the end-of-run
# performance summary (every execution mode).
register_package_versions("petsc4py", "firedrake", "icepack")

# --- Initialize MPI ---
rank, size, comm, _ = ParallelManager().icesee_mpi_init(icesee_kwargs)

PETSc.Sys.Print("Fetching the model parameters ...")

# --- Ensemble Parameters ---
icesee_kwargs.update({
"nt": int(float(icesee_kwargs["num_years"])) * int(float(icesee_kwargs["timesteps_per_year"])),
"dt": 1.0 / float(icesee_kwargs["timesteps_per_year"])
})

# --- Model intialization ---
# A spare rank under the hierarchical resource plan (resource_plan.py,
# reached via ParallelManager().icesee_mpi_init above) belongs to no
# model group and is handed `rank=None`/`comm=MPI.COMM_NULL` -- it has no
# mesh to build and must never call any Firedrake mesh/function-space
# constructor (those are collectives over the model-group communicator;
# a spare rank calling one alone, or a mesh constructor being handed
# COMM_NULL, both fail/hang). Every other rank (rank is not None) builds
# its model exactly as before.
if rank is not None:
    PETSc.Sys.Print("Initializing icepack model ...")
    icesee_kwargs.update({'comm': comm})
    icesee_kwargs = initialize_model(**icesee_kwargs)   # icesee_kwargs now already has nx,ny,Lx,Ly,x,y,h,u,a,a_p,b,b_in,b_out,
                                           # h0,u0,solver_weertman,A,C,Q,V,mesh

    icesee_kwargs["nd"] = icesee_kwargs["h0"].dat.data.size * icesee_kwargs["total_state_param_vars"]

    # only genuinely new additions remain:
    icesee_kwargs.update({
        "da": float(icesee_kwargs["da"]),
        "dt": icesee_kwargs["dt"],
        "seed": float(icesee_kwargs["seed"]),
        "h_nurge_ic": float(icesee_kwargs["h_nurge_ic"]),
        "u_nurge_ic": float(icesee_kwargs["u_nurge_ic"]),
        "nurged_entries_percentage": float(icesee_kwargs["nurged_entries_percentage"]),
        "a_in_p": float(icesee_kwargs["a_in_p"]),
        "da_p": float(icesee_kwargs["da_p"]),
        "solver": icesee_kwargs["solver_weertman"],
        "nd": icesee_kwargs["nd"],
    })

    # --- nurged smb
    a_in = firedrake.Constant(icesee_kwargs["a_in_p"])
    da_p = firedrake.Constant(icesee_kwargs["da_p"])
    a_nuged = firedrake.Function(icesee_kwargs["Q"]).interpolate(a_in + da_p*icesee_kwargs["x"]/icesee_kwargs["Lx"])
    icesee_kwargs.update({"a_nuged":a_nuged})
else:
    # Placeholder-only, never dereferenced for real physics: `nd` (0) is
    # read unconditionally near the top of
    # icesee_model_data_assimilation_full_parallel (every world rank,
    # spare included, must have a value present -- the driver's own
    # early comm_world.bcast of the real global nd, added for Stage 4C,
    # supplies the value every rank actually uses downstream). Lx/Ly are
    # ALSO read unconditionally by every world rank early in that same
    # driver (Lx_dim = sqrt(Lx*Ly), for process-noise length-scale setup
    # -- reached before any color-is-not-None guard), so a spare rank
    # needs a real numeric value there too, not just a default that a
    # pre-existing raw-string config value would shadow via .get() --
    # match initialize_model's own int(float(...)) cast exactly so every
    # rank agrees on the same value.
    icesee_kwargs.update({
        'comm': comm,
        "nd": 0,
        "a_nuged": None,
        "Lx": int(float(icesee_kwargs["Lx"])),
        "Ly": int(float(icesee_kwargs["Ly"])),
    })

# --- Run Data Assimilation ---
PETSc.Sys.Print("Data assimilation with ICESEE ...")
icesee_model_data_assimilation(**icesee_kwargs)
