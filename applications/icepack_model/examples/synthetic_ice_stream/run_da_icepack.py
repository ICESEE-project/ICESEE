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

from ICESEE.config._utility_imports import *
from ICESEE.config._utility_imports import params, kwargs, modeling_params, enkf_params, physical_params
from ICESEE.applications.icepack_model.examples.synthetic_ice_stream._icepack_model import initialize_model
from ICESEE.src.run_model_da.run_models_da import icesee_model_data_assimilation
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager

# --- Initialize MPI ---
rank, size, comm, _ = ParallelManager().icesee_mpi_init(params)

PETSc.Sys.Print("Fetching the model parameters ...")

# --- Ensemble Parameters ---
params.update({
"nt": int(float(modeling_params["num_years"])) * int(float(modeling_params["timesteps_per_year"])),
"dt": 1.0 / float(modeling_params["timesteps_per_year"])
})

# --- Model intialization --- 
PETSc.Sys.Print("Initializing icepack model ...")
kwargs.update({'comm': comm})
kwargs = initialize_model(**kwargs)   # kwargs now already has nx,ny,Lx,Ly,x,y,h,u,a,a_p,b,b_in,b_out,
                                       # h0,u0,solver_weertman,A,C,Q,V,mesh

params["nd"] = kwargs["h0"].dat.data.size * params["total_state_param_vars"]

# only genuinely new additions remain:
kwargs.update({
    "da": float(modeling_params["da"]),
    "dt": params["dt"],
    "seed": float(enkf_params["seed"]),
    "h_nurge_ic": float(enkf_params["h_nurge_ic"]),
    "u_nurge_ic": float(enkf_params["u_nurge_ic"]),
    "nurged_entries_percentage": float(enkf_params["nurged_entries_percentage"]),
    "a_in_p": float(modeling_params["a_in_p"]),
    "da_p": float(modeling_params["da_p"]),
    "solver": kwargs["solver_weertman"],   
    "nd": params["nd"],
})

# --- nurged smb
a_in = firedrake.Constant(kwargs["a_in_p"])
da_p = firedrake.Constant(kwargs["da_p"])
a_nuged = firedrake.Function(kwargs["Q"]).interpolate(a_in + da_p*kwargs["x"]/kwargs["Lx"])
kwargs.update({"a_nuged":a_nuged})

# --- Run Data Assimilation ---
kwargs.update({'params': params}) # update the kwargs with the parameters

PETSc.Sys.Print("Data assimilation with ICESEE ...")
icesee_model_data_assimilation(**kwargs)


