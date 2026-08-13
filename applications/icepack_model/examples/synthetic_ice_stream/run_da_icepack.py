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

# --- firedrake imports ---
import firedrake
from firedrake.petsc import PETSc

from ICESEE.config._utility_imports import icesee_kwargs
from ICESEE.applications.icepack_model.examples.synthetic_ice_stream._icepack_model import initialize_model
from ICESEE.src.run_model_da.run_models_da import icesee_model_data_assimilation
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager

# --- Initialize MPI ---
rank, size, comm, _ = ParallelManager().icesee_mpi_init(icesee_kwargs)

PETSc.Sys.Print("Fetching the model parameters ...")

# --- Ensemble Parameters ---
icesee_kwargs.update({
"nt": int(float(icesee_kwargs["num_years"])) * int(float(icesee_kwargs["timesteps_per_year"])),
"dt": 1.0 / float(icesee_kwargs["timesteps_per_year"])
})

# --- Model intialization ---
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

# --- Run Data Assimilation ---
PETSc.Sys.Print("Data assimilation with ICESEE ...")
icesee_model_data_assimilation(**icesee_kwargs)
