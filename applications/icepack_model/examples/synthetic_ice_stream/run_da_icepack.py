# =============================================================================
# @author: Brian Kyanjo
# @date: 2024-11-06
# @description: Synthetic ice stream with data assimilation
# =============================================================================

# --- Imports ---
import sys
import os
import numpy as np

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
nx,ny,Lx,Ly,x,y,h,u,a,a_p,b,b_in,b_out,h0,u0,solver_weertman,A,C,Q,V = initialize_model(
    physical_params, modeling_params,comm
)

# update the parameters
params["nd"]=  h0.dat.data.size * params["total_state_param_vars"] # get the size of the entire vector
kwargs.update({"a":a, "h0":h0, "u0":u0, "C":C, "A":A,"Q":Q,"V":V, "da":float(modeling_params["da"]),
        "b":b, "dt":params["dt"], "seed":float(enkf_params["seed"]), "x":x, "y":y,
        "Lx":Lx, "Ly":Ly, "nx":nx, "ny":ny, "h_nurge_ic":float(enkf_params["h_nurge_ic"]), 
        "u_nurge_ic":float(enkf_params["u_nurge_ic"]),"nurged_entries_percentage":float(enkf_params["nurged_entries_percentage"]),
        "a_in_p":float(modeling_params["a_in_p"]), "da_p":float(modeling_params["da_p"]),
        "solver":solver_weertman,
        "a_p":a_p, "b_in":b_in, "b_out":b_out,
})

# --- nurged smb
a_in = firedrake.Constant(kwargs["a_in_p"])
da_p = firedrake.Constant(kwargs["da_p"])
a_nuged = firedrake.interpolate(a_in + da_p*kwargs["x"]/kwargs["Lx"], kwargs["Q"])
kwargs.update({"a_nuged":a_nuged})

# --- Run Data Assimilation ---
kwargs.update({'params': params}) # update the kwargs with the parameters

PETSc.Sys.Print("Data assimilation with ICESEE ...")
icesee_model_data_assimilation(
    enkf_params["model_name"],
    enkf_params["filter_type"],
    **kwargs  
)



