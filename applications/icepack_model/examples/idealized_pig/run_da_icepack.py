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

from modelfunc import myerror
import modelfunc as mf
from modelfunc import firedrakeSmooth, flotationHeight, flotationMask

from ICESEE.config._utility_imports import *
#from ICESEE.config._utility_imports import params, kwargs, paths, modeling_params, enkf_params
from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_model import initialize_model, initialState, initializeMesh
from ICESEE.src.run_model_da.run_models_da import icesee_model_data_assimilation
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager


# --- Initialize MPI ---

rank, size, comm, _ = ParallelManager().icesee_mpi_init(params)

PETSc.Sys.Print("Fetching the model parameters ...")



# --- Ensemble Parameters ---

num_years = float(modeling_params["num_years"])
dt = float(modeling_params["timesteps_per_year"])   # time step size
nt = int(round(num_years / dt))     # total number of time steps

params.update({"nt": nt, "dt": dt}) # update the parameter dictionary
kwargs.update({"nt": nt, "dt": dt}) # update kwargs to use in other icepack functions (e.g. BasalMeltRate)




# --- Model initialization --- 

PETSc.Sys.Print("Initializing icepack model ...")

kwargs.update({
    "comm": comm,
    "initFile": modeling_params["initFile"],
    "paramsFile": modeling_params["paramsFile"],
    "meshFile": modeling_params["meshFile"],
    "SMBFile": modeling_params["SMBFile"],
    "dt": modeling_params["timesteps_per_year"],
    "num_years": modeling_params["num_years"],
    "bmr_increase_time": int(modeling_params["bmr_increase_time"]),
    "save_steps": modeling_params["save_steps"],
    #"hThresh": modeling_params["hThresh"]
})

h, h0, s, s0, u, bed, zF, grounded, floating, A0, beta0, smb, basal_melt_field, Q, V, forward_solver = initialize_model(**kwargs)



# ----- Update the parameters ----

params["nd"] = h0.dat.data.size * params["total_state_param_vars"] # get the size of the entire vector


kwargs.update({

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
    "seed":float(enkf_params["seed"]),  
    "zF": zF,
    "grounded": grounded,
    "floating": floating,
    "wrong_basal_melt_field": float(enkf_params["wrong_basal_melt_field"]), 
    "solver": forward_solver,
    "nd": params["nd"],
    #"Lx":float(physical_params["Lx"]), 
    #"Ly":float(physical_params["Ly"]), 
    #"nx":float(physical_params["nx"]), 
    #"ny":float(physical_params["ny"]),
 
})



# ----- Nudge the basal melt rate field -----

wrong_bmr = firedrake.Constant(kwargs["wrong_basal_melt_field"])

bmr_nudged = firedrake.interpolate(kwargs["basal_melt_field"] + wrong_bmr, kwargs["Q"])

kwargs.update({"bmr_nudged": bmr_nudged})



# --- Run Data Assimilation ---
kwargs.update({'params': params}) # update the kwargs with the parameters

PETSc.Sys.Print("Data assimilation with ICESEE ...")
icesee_model_data_assimilation(**kwargs)


