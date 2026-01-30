# ==============================================================================
# @des: This file contains run functions for flowline data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2026-01-26
# @author: Brian Kyanjo
# ==============================================================================

import sys
import os
import h5py
import numpy as np

# --- import run_simulation function from the flowline model ---
from ICESEE.applications.flowline_model.examples.flowline_1d._flowline_model import *
from ICESEE.config._utility_imports import icesee_get_index

# --- Forecast step for the flowline model ---
def forecast_step_single(ensemble=None, **kwargs):
    """inputs: run_simulation - function that runs the model
                ensemble - current state of the model
                dt - time step
                *args - additional arguments for the model
         outputs: uai - updated state of the model after one time step
    """

    #  call the run_model fun to push the state forward in time
    return run_model(ensemble, **kwargs)

# --- generate true state ---
def generate_true_state(**kwargs):
    """generate the true state of the model"""

    # Unpack the parameters
    params = kwargs["params"]
    statevec_true = kwargs["statevec_true"]

    nd, nt = statevec_true.shape

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(**kwargs)
    
    # Set the initial condition
    huxg_out0 = initialize_model(**kwargs)
    
    kwargs['facemelt'] = np.linspace(5, 85, kwargs["NT"]+1)/float(kwargs['year'])  # Varying terminus melt over time
    fm_dist = np.random.normal(0,20.0)
    fm_truth = kwargs["facemelt"]
    kwargs['transient'] = 1

    # print("[DEBUG] Initial huxg_out0:", huxg_out0, statevec_true[:, 0].shape )  # Debug print statement 
    statevec_true[:, 0] = huxg_out0

    NX = kwargs['NX']

    # Run the model forward in time
    for k in range(nt-1):
        kwargs["tcurrent"] = k+1
        huxg_out0 = run_model(statevec_true[:, k], **kwargs)
        for key, value in huxg_out0.items():
            statevec_true[indx_map[key], k + 1] = value
  
    updated_state = {}
    for key in kwargs['vec_inputs']:
        updated_state[key] = statevec_true[indx_map[key],:]
    return updated_state

def generate_nurged_state(**kwargs):
    """generate the nurged state of the model"""

    # Unpack the parameters
    params = kwargs["params"]
    statevec_nurged = kwargs["statevec_nurged"]

    nd,nt = statevec_nurged.shape
    # Set the initial condition
    huxg_out0 = initialize_model(**kwargs)

    statevec_nurged[:, 0] = huxg_out0
    kwargs['facemelt'] = np.linspace(5, 45, kwargs["NT"]+1)/float(kwargs['year'])  # Varying terminus melt over time
    fm_dist = np.random.normal(0,20.0)
    fm_truth = kwargs["facemelt"]
    kwargs['transient'] = 1
    NX = kwargs['NX']

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(**kwargs)

    # Run the model forward in time
    for k in range(nt-1):
        kwargs["tcurrent"] = k+1
        huxg_out0 = run_model(statevec_nurged[:, k], **kwargs)
        for key, value in huxg_out0.items():
            statevec_nurged[indx_map[key], k + 1] = value

    updated_state = {}
    for key in kwargs['vec_inputs']:
        updated_state[key] = statevec_nurged[indx_map[key],:]
    return updated_state
    
# --- initialize the ensemble members ---
def initialize_ensemble(ens, **kwargs):
    """initialize the ensemble members"""
    # Unpack the parameters
    params = kwargs["params"]
    statevec_ens    = kwargs["statevec_ens"]
    
    nd, N = statevec_ens.shape
    hdim = nd // params["num_state_vars"]

    vecs, indx_map, dim_per_proc = icesee_get_index(**kwargs)

    # Set the initial condition
    huxg_out0 = initialize_model(**kwargs)
    kwargs['facemelt'] = np.linspace(5, 85, kwargs["NT"]+1)/float(kwargs['year'])  # Varying terminus melt over time
    fm_dist = np.random.normal(0,20.0)
    fm_truth = kwargs["facemelt"]
    kwargs['transient'] = 1
    NX = kwargs['NX']

    initail_state = np.concatenate((huxg_out0[:-1], [kwargs['facemelt'][0]/ kwargs['uscale']]))
    initailized_state={}
    # print(f"[ICESEE] Rank: {ens} Initial state size: {initail_state.shape}, ensemble shape: {statevec_ens[:,ens].shape} huxg_out0 shape: {huxg_out0.shape}")
    for ii, key in enumerate(kwargs['vec_inputs']):
        initailized_state[key] = initail_state[indx_map[key]]

    return initailized_state

# the user should able to override the default observation operator and its jacobian
def H(m,n):

    H = np.zeros((m*2+1,n))

    # calculate the distance between measurements
    di = int((n-2)/(2*m))
    for i in range(m):
        H[i, i*di] = 1
        H[m+i, int((n-2)/2) + i*di] = 1

    H[2*m, n-2] = 1
    return H

def Obs_fun(virtual_obs=None):
    obs_file = '_modelrun_datasets/synthetic_obs.h5'
    with h5py.File(obs_file, 'r') as f:
       obs = f['hu_obs'][:]
       n,m = obs.shape
    z =  H(m,n) @ virtual_obs
    return z

def JObs_fun(nd=None):
    obs_file = '_modelrun_datasets/synthetic_obs.h5'
    with h5py.File(obs_file, 'r') as f:
       obs = f['hu_obs'][:]
       nd,m = obs.shape
    return H(m,nd)

def Cov_Obs_fun(sig_obs=None,nd=None, kwargs=None):
    R_cov = (sig_obs**2) * np.eye(2 * kwargs["m_obs"] + 1)
    return R_cov

def localization_function(**kwargs):
    grid = kwargs['grid']
    statevec_sig = np.concatenate((grid['sigma_elem'], grid['sigma'], np.array([1, 1])))
    taper = np.ones((statevec_sig.shape[0], statevec_sig.shape[0]))
    taper[-1, -3] = 2  
    taper[-3, -1] = 2  
    taper[-1, -1] = 10  
    taper[-2, -1] = 10  
    taper[-1, -2] = 10 

    return taper

def post_analysis_update(ensemble=None, **kwargs):
    """update the model state after the analysis step"""

    k = kwargs['k']

    kwargs['facemelt'][k+1:] = ensemble[-1] * kwargs['uscale'] * np.ones_like(kwargs['facemelt'][k+1:])
    return kwargs, ensemble
