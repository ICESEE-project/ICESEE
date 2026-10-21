# =============================================================================
# @author: Brian Kyanjo
# @date: 2026-01-23
# @description: Flowline 1D model with data assimilation
# =============================================================================

# --- Imports ---
import sys
import os
import numpy as np
from pathlib import Path

# --- Set up paths ---
os.chdir(Path(__file__).resolve().parent)

# --- JAX configuration ---
import jax
# Set the precision in JAX to use float64
jax.config.update("jax_enable_x64", True)

# --- ICESEE imports ---
from ICESEE.config._utility_imports import icesee_kwargs
from ICESEE.src.run_model_da.run_models_da import icesee_model_data_assimilation
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager

# --- Flowline 1D model imports ---
from ICESEE.applications.flowline_model.examples.flowline_1d._flowline_config import build_flowline_config

# --- Register execution_mode 3 support (import-time side effect) ---
from ICESEE.applications.flowline_model.examples.flowline_1d import mode3_runner  # noqa: F401

# --- Initialize MPI ---
rank, size, comm, _ = ParallelManager().icesee_mpi_init(icesee_kwargs)

# --- Model configuration (shared with mode3_runner.py) ---
build_flowline_config(icesee_kwargs)

# The MPI analysis of modes 1 and 2 assumes every state variable has the same
# number of entries; flowline_1d's grounding line is a single value next to
# the NX-point h and u fields, so those modes would assimilate incorrectly.
if int(icesee_kwargs.get("execution_mode", 3)) in (1, 2):
    raise SystemExit(
        "[ICESEE] flowline_1d supports execution_mode 0 (serial) and 3 "
        "(distributed); modes 1 and 2 need equal-sized state variables."
    )

# --- Run the model with data assimilation ---
# call ICESEE data assimilation function
icesee_model_data_assimilation(**icesee_kwargs)
