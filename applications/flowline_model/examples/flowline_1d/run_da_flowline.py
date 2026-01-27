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
from ICESEE.config._utility_imports import *
from ICESEE.config._utility_imports import params, kwargs, modeling_params, enkf_params, physical_params
from ICESEE.src.run_model_da.run_models_da import icesee_model_data_assimilation
from ICESEE.src.parallelization.parallel_mpi.icesee_mpi_parallel_manager import ParallelManager

# --- Flowline 1D model imports ---
from ICESEE.applications.flowline_model.examples.flowline_1d._flowline_model import initialize_model

# --- Initialize MPI ---
rank, size, comm, _ = ParallelManager().icesee_mpi_init(params)

# --- Ensemble Parameters ---
params.update({"nd": int(float(enkf_params["num_state_vars"]))})

# --- model parameters ---
kwargs.update({ "nt": int(float(modeling_params["num_years"])),
               "NT": int(float(modeling_params["num_years"])),
                "nd": params["nd"],
                "seed":float(enkf_params["seed"]),
                "t":np.linspace(0, int(float(modeling_params["num_years"])), params["nt"] + 1), 
                'hscale': float(physical_params['hscale']),
                'A': float(physical_params['A']),
                'n': int(physical_params['n']),
                'C': float(physical_params['C']),
                'rho_ice': float(physical_params['rho_ice']),
                'rho_water': float(physical_params['rho_water']),
                'g': float(physical_params['g']),
                'accum': float(physical_params['accum'])/float(modeling_params['year']),
                'facemelt': float(physical_params['facemelt'])/float(modeling_params['year']),
                'm': 1/int(physical_params['n']),
                'B': float(physical_params['A']) ** (-1 / int(physical_params['n'])),
                'ascale': 1.0 / float(modeling_params['year']),
                'N1': int(physical_params['N1']),
                'N2': int(physical_params['N2']),
                'NX': int(physical_params['N1']) + int(physical_params['N2']),
                'TF': float(modeling_params['year']),
                'sigGZ': float(physical_params['sigGZ']),
                'sigma1': np.linspace(float(physical_params['sigGZ']) / (int(physical_params['N1']) + 0.5), float(physical_params['sigGZ']), int(physical_params['N1'])),
                'sigma2': np.linspace(float(physical_params['sigGZ']), 1, int(physical_params['N2'] + 1)),
                'sillamp': float(physical_params['sillamp']),
                'sillsmooth': float(physical_params['sillsmooth']),
                'xsill': float(physical_params['xsill']),
                'tcurrent': int(physical_params['tcurrent']),
                'transient': int(physical_params['transient']),
                'uscale': float(physical_params['rho_ice']) * float(physical_params['g']) * float(physical_params['hscale']) * (1.0 / float(modeling_params['year'])) / float(physical_params['C']),
                'scalar_inputs': enkf_params.get('scalar_inputs', []),
})

# kwargs-kwargs update
xscale = kwargs['uscale'] * kwargs['hscale'] / kwargs['ascale']
sigma = np.concatenate((kwargs['sigma1'], kwargs['sigma2'][1:kwargs['N2'] + 1]))
kwargs.update({'xscale': xscale,
                'tscale': (xscale / kwargs['uscale']),
                'eps': kwargs['B'] * ((kwargs['uscale'] / xscale) ** (1 / kwargs['n'])) / (2 * kwargs['rho_ice'] * kwargs['g'] * kwargs['hscale']),
                'lambda': 1 - (kwargs['rho_ice'] / kwargs['rho_water']),
                'dt': kwargs['TF'] / int(float(modeling_params['num_years'])),
                'sigma': sigma,
                'grid': { 'sigma': sigma,
                          'sigma_elem':np.concatenate(([0], (sigma[:-1] + sigma[1:]) / 2)),
                          'dsigma': np.diff(sigma)
                        }
               })

# initialize the flowline model initial condition
huxg_out0 = initialize_model(**kwargs)
kwargs.update({'nd': huxg_out0.shape[0]})
var_nd = {var: (1 if var in kwargs['scalar_inputs'] else kwargs['NX']) for var in kwargs['vec_inputs']}
kwargs.update({'var_nd': var_nd})

# --- Run the model with data assimilation ---
kwargs.update({'params': params}) # update the kwargs with the parameters

# call ICESEE data assimilation function 
icesee_model_data_assimilation(**kwargs)
