# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-01-13
# @description: Lorenz96 model with data assimilation
# =============================================================================

# --- Imports ---
import sys
import os
import numpy as np
from pathlib import Path

# --- Set up paths ---
os.chdir(Path(__file__).resolve().parent)

# --- ICESEE imports ---
from ICESEE.config._utility_imports import icesee_kwargs
from ICESEE.src.run_model_da.run_models_da import icesee_model_data_assimilation
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager

# --- Lorenz96 model imports ---
from ICESEE.applications.lorenz_model.examples.lorenz96._lorenz96_model import initialize_model

# --- Initialize MPI ---
rank, size, comm, _ = ParallelManager().icesee_mpi_init(icesee_kwargs)

# --- Ensemble Parameters ---
icesee_kwargs.update({"nt": int(float(icesee_kwargs["num_years"])/float(icesee_kwargs["dt"])),"nd": int(float(icesee_kwargs["num_state_vars"]))})

# --- model parameters ---
icesee_kwargs.update({ "nt": icesee_kwargs["nt"],
                "nd": icesee_kwargs["nd"],
                "dt": float(icesee_kwargs["dt"]), "seed":float(icesee_kwargs["seed"]),
                "t":np.linspace(0, int(float(icesee_kwargs["num_years"])), icesee_kwargs["nt"] + 1),
                "u0True": np.array([1,1,1]), "u0b": np.array([2.0,3.0,4.0]),
                "sigma_96":float(icesee_kwargs["sigma_96"]), "beta_96":eval(icesee_kwargs["beta_96"]),
                "rho_96":float(icesee_kwargs["rho_96"]),
})


# --- Run the model with data assimilation ---
# call ICESEE data assimilation function
icesee_model_data_assimilation(**icesee_kwargs)
