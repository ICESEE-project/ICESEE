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
from ICESEE.applications.flowline_model.examples.flowline_1d._flowline_model import initialize_model

# --- Initialize MPI ---
rank, size, comm, _ = ParallelManager().icesee_mpi_init(icesee_kwargs)

# --- Ensemble Parameters ---
icesee_kwargs.update({"nd": int(float(icesee_kwargs["num_state_vars"]))})

# --- model parameters ---
icesee_kwargs.update({ "nt": int(float(icesee_kwargs["num_years"])),
               "NT": int(float(icesee_kwargs["num_years"])),
                "nd": icesee_kwargs["nd"],
                "seed":float(icesee_kwargs["seed"]),
                "t":np.linspace(0, int(float(icesee_kwargs["num_years"])), icesee_kwargs["nt"] + 1),
                'hscale': float(icesee_kwargs['hscale']),
                'A': float(icesee_kwargs['A']),
                'n': int(icesee_kwargs['n']),
                'C': float(icesee_kwargs['C']),
                'rho_ice': float(icesee_kwargs['rho_ice']),
                'rho_water': float(icesee_kwargs['rho_water']),
                'g': float(icesee_kwargs['g']),
                'accum': float(icesee_kwargs['accum'])/float(icesee_kwargs['year']),
                'facemelt': float(icesee_kwargs['facemelt'])/float(icesee_kwargs['year']),
                'm': 1/int(icesee_kwargs['n']),
                'B': float(icesee_kwargs['A']) ** (-1 / int(icesee_kwargs['n'])),
                'ascale': 1.0 / float(icesee_kwargs['year']),
                'N1': int(icesee_kwargs['N1']),
                'N2': int(icesee_kwargs['N2']),
                'NX': int(icesee_kwargs['N1']) + int(icesee_kwargs['N2']),
                'TF': float(icesee_kwargs['year']),
                'sigGZ': float(icesee_kwargs['sigGZ']),
                'sigma1': np.linspace(float(icesee_kwargs['sigGZ']) / (int(icesee_kwargs['N1']) + 0.5), float(icesee_kwargs['sigGZ']), int(icesee_kwargs['N1'])),
                'sigma2': np.linspace(float(icesee_kwargs['sigGZ']), 1, int(icesee_kwargs['N2'] + 1)),
                'sillamp': float(icesee_kwargs['sillamp']),
                'sillsmooth': float(icesee_kwargs['sillsmooth']),
                'xsill': float(icesee_kwargs['xsill']),
                'tcurrent': int(icesee_kwargs['tcurrent']),
                'transient': int(icesee_kwargs['transient']),
                'uscale': float(icesee_kwargs['rho_ice']) * float(icesee_kwargs['g']) * float(icesee_kwargs['hscale']) * (1.0 / float(icesee_kwargs['year'])) / float(icesee_kwargs['C']),
                'scalar_inputs': icesee_kwargs.get('scalar_inputs', []),
})

# Add application-derived values to the ICESEE runtime context.
xscale = icesee_kwargs['uscale'] * icesee_kwargs['hscale'] / icesee_kwargs['ascale']
sigma = np.concatenate((icesee_kwargs['sigma1'], icesee_kwargs['sigma2'][1:icesee_kwargs['N2'] + 1]))
icesee_kwargs.update({'xscale': xscale,
                'tscale': (xscale / icesee_kwargs['uscale']),
                'eps': icesee_kwargs['B'] * ((icesee_kwargs['uscale'] / xscale) ** (1 / icesee_kwargs['n'])) / (2 * icesee_kwargs['rho_ice'] * icesee_kwargs['g'] * icesee_kwargs['hscale']),
                'lambda': 1 - (icesee_kwargs['rho_ice'] / icesee_kwargs['rho_water']),
                'dt': icesee_kwargs['TF'] / int(float(icesee_kwargs['num_years'])),
                'sigma': sigma,
                'grid': { 'sigma': sigma,
                          'sigma_elem':np.concatenate(([0], (sigma[:-1] + sigma[1:]) / 2)),
                          'dsigma': np.diff(sigma)
                        }
               })

# initialize the flowline model initial condition
huxg_out0 = initialize_model(**icesee_kwargs)
icesee_kwargs.update({'nd': huxg_out0.shape[0]})
var_nd = {var: (1 if var in icesee_kwargs['scalar_inputs'] else icesee_kwargs['NX']) for var in icesee_kwargs['vec_inputs']}
icesee_kwargs.update({'var_nd': var_nd})

# --- Run the model with data assimilation ---
# call ICESEE data assimilation function
icesee_model_data_assimilation(**icesee_kwargs)
