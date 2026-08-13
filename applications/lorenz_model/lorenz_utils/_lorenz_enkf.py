# ==============================================================================
# @des: This file contains run functions for lorenz data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2025-01-13
# @author: Brian Kyanjo
# ==============================================================================

import sys
import os
import numpy as np

# --- import run_simulation function from the lorenz96 model ---
from ICESEE.applications.lorenz_model.examples.lorenz96._lorenz96_model import *
from ICESEE.config._utility_imports import icesee_get_index

# --- Forecast step for the Lorenz96 model ---
def forecast_step_single(ensemble=None, **icesee_kwargs):
    """inputs: run_simulation - function that runs the model
                ensemble - current state of the model
                dt - time step
                *args - additional arguments for the model
         outputs: uai - updated state of the model after one time step
    """

    #  call the run_model fun to push the state forward in time
    return run_model(ensemble, **icesee_kwargs)

# --- Background step for the Lorenz96 model ---
def background_step(k=None,statevec_bg=None, hdim=None, **icesee_kwargs):
    """inputs: k - current time step
                run_simulation - function that runs the model
                state - current state of the model
                dt - time step
                *args - additional arguments for the model
        outputs: state - updated state of the model after one time step
    """
    # Call the run_simulationfunction to push the state forward in time
    statevec_bg[:,k+1] = run_model(statevec_bg[:,k], **icesee_kwargs)

    return statevec_bg


# --- generate true state ---
def generate_true_state(**icesee_kwargs):
    """generate the true state of the model"""

    # Unpack the parameters
    statevec_true = icesee_kwargs["statevec_true"]

    nd = icesee_kwargs['nd']
    nt = icesee_kwargs['nt']
    dt = icesee_kwargs['dt']
    num_state_vars = icesee_kwargs['num_state_vars']
    u0True = icesee_kwargs.get('u0True', None)

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(statevec_true, **icesee_kwargs)


    # Set the initial condition
    statevec_true[:, 0] = u0True

    # Run the model forward in time
    for k in range(nt):
        state = run_model(statevec_true[:, k], **icesee_kwargs)
        statevec_true[indx_map['x'], k + 1] = state['x']
        statevec_true[indx_map['y'], k + 1] = state['y']
        statevec_true[indx_map['z'], k + 1] = state['z']


    updated_state = {'x' : statevec_true[indx_map['x'],:],
                     'y' : statevec_true[indx_map['y'],:],
                    'z' : statevec_true[indx_map['z'],:]}
    return updated_state

def generate_nurged_state(**icesee_kwargs):
    """generate the nurged state of the model"""

    # Unpack the parameters
    statevec_nurged = icesee_kwargs["statevec_nurged"]

    nd = icesee_kwargs['nd']
    nt = icesee_kwargs['nt']
    dt = icesee_kwargs['dt']
    num_state_vars = icesee_kwargs['num_state_vars']
    u0True = icesee_kwargs.get('u0True', None)

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(statevec_nurged, **icesee_kwargs)

    # Set the initial condition
    statevec_nurged[:, 0] = u0True

    # Run the model forward in time
    for k in range(nt):
        statevec_nurged[:, k + 1] = run_model(statevec_nurged[:, k], **icesee_kwargs)

    updated_state = {'x' : statevec_nurged[indx_map['x'],:],
                     'y' : statevec_nurged[indx_map['y'],:],
                    'z' : statevec_nurged[indx_map['z'],:]}
    return updated_state

# --- initialize the ensemble members ---
def initialize_ensemble(ens, **icesee_kwargs):
    """initialize the ensemble members"""
    # Unpack the parameters
    statevec_ens    = icesee_kwargs["statevec_ens"]

    nd, N = statevec_ens.shape
    hdim = nd // icesee_kwargs["num_state_vars"]

    u0b = icesee_kwargs.get('u0b', None)
    intialized_state = {'x' : u0b[1],
                        'y' : u0b[2],
                        'z' : u0b[3]}
    return intialized_state
