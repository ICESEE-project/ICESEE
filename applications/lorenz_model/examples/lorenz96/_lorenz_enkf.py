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
def forecast_step_single(ensemble=None, **kwargs):
    """inputs: run_simulation - function that runs the model
                ensemble - current state of the model
                dt - time step
                *args - additional arguments for the model
         outputs: uai - updated state of the model after one time step
    """

    #  call the run_model fun to push the state forward in time
    return run_model(ensemble, **kwargs)

# --- Background step for the Lorenz96 model ---
def background_step(k=None,statevec_bg=None, hdim=None, **kwargs):
    """inputs: k - current time step
                run_simulation - function that runs the model
                state - current state of the model
                dt - time step
                *args - additional arguments for the model
        outputs: state - updated state of the model after one time step
    """
    # Call the run_simulationfunction to push the state forward in time
    statevec_bg[:,k+1] = run_model(statevec_bg[:,k], **kwargs)

    return statevec_bg


# --- generate true state ---
def generate_true_state(**kwargs):
    """generate the true state of the model"""

    # Unpack the parameters
    params = kwargs["params"]
    statevec_true = kwargs["statevec_true"]

    nd = params['nd']
    # nt = params['nt']
    dt = params['dt']
    num_state_vars = params['num_state_vars']
    u0True = kwargs.get('u0True', None)

    nd, nt = statevec_true.shape

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(**kwargs)
    

    # Set the initial condition
    statevec_true[:, 0] = u0True

    # Run the model forward in time
    for k in range(nt-1):
        state = run_model(statevec_true[:, k], **kwargs)
        statevec_true[indx_map['x'], k + 1] = state['x']
        statevec_true[indx_map['y'], k + 1] = state['y']
        statevec_true[indx_map['z'], k + 1] = state['z']
         

    updated_state = {'x' : statevec_true[indx_map['x'],:],
                     'y' : statevec_true[indx_map['y'],:],
                    'z' : statevec_true[indx_map['z'],:]}
    return updated_state

def generate_nurged_state(**kwargs):
    """generate the nurged state of the model"""

    # Unpack the parameters
    params = kwargs["params"]
    statevec_nurged = kwargs["statevec_nurged"]

    # nd = params['nd']
    # nt = params['nt']
    nd,nt = statevec_nurged.shape
    dt = params['dt']
    num_state_vars = params['num_state_vars']
    u0b = kwargs.get('u0b', None)

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(**kwargs)

    # Set the initial condition
    statevec_nurged[:, 0] = u0b

    # Run the model forward in time
    for k in range(nt-1):
        state = run_model(statevec_nurged[:, k], **kwargs)
        statevec_nurged[indx_map['x'], k + 1] = state['x']
        statevec_nurged[indx_map['y'], k + 1] = state['y']
        statevec_nurged[indx_map['z'], k + 1] = state['z']

    updated_state = {'x' : statevec_nurged[indx_map['x'],:],
                     'y' : statevec_nurged[indx_map['y'],:],
                    'z' : statevec_nurged[indx_map['z'],:]}
    
    return updated_state
    
# --- initialize the ensemble members ---
def initialize_ensemble(ens, **kwargs):
    """initialize the ensemble members"""
    # Unpack the parameters
    params = kwargs["params"]
    statevec_ens    = kwargs["statevec_ens"]
    
    nd, N = statevec_ens.shape
    hdim = nd // params["num_state_vars"]

    u0b = kwargs.get('u0b', None)
    intialized_state = {'x' : u0b[0],
                        'y' : u0b[1],
                        'z' : u0b[2]}
    return intialized_state

# the user should able to override the default observation operator and its jacobian
def Obs_fun(virtual_obs=None):
    w = virtual_obs.copy()
    return w

def JObs_fun(nd=None):
    H_jac = np.eye(nd)
    return H_jac

def Cov_Obs_fun(sig_obs=None,nd=None, kwargs=None):
    R_cov = np.eye(nd) * (sig_obs ** 2)
    return R_cov