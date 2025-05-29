# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-01-13
# @description: Lorenz96 model with data assimilation
# =============================================================================

# --- Imports ---
import sys
import os
import numpy as np

# --- ICESEE imports ---
from ICESEE.config._utility_imports import *
from ICESEE.config._utility_imports import params, kwargs, modeling_params, enkf_params, physical_params
from ICESEE.src.run_model_da.run_models_da import icesee_model_data_assimilation
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager

# --- Lorenz96 model imports ---
from ICESEE.applications.lorenz_model.examples.lorenz96._lorenz96_model import initialize_model

# --- Initialize MPI ---
rank, size, comm, _ = ParallelManager().icesee_mpi_init(params)

# --- Ensemble Parameters ---
params.update({
    "nt": int(float(modeling_params["num_years"])/float(modeling_params["dt"])),
    "dt": float(modeling_params["dt"]),
    "nd": int(float(enkf_params["num_state_vars"])),
})

# --- model parameters ---
kwargs.update({ "dt":params["dt"], "seed":float(enkf_params["seed"]),
          "t":np.linspace(0, int(float(modeling_params["num_years"])), params["nt"] + 1), 
          "u0True": np.array([1,1,1]), "u0b": np.array([2.0,3.0,4.0]), 
          "sigma":float(physical_params["sigma"]), "beta":eval(physical_params["beta"]),
          "rho":float(physical_params["rho"]),
})


# --- Run the model with data assimilation ---
kwargs.update({'params': params}) # update the kwargs with the parameters

# call ICESEE data assimilation function 
icesee_model_data_assimilation(**kwargs)
