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
from ICESEE.applications.flowline_model.examples.flowline_1d._flowline_config import (
    FACEMELT_TRUE_END, FACEMELT_WRONG_END, facemelt_schedule,
)

# --- Forecast step for the flowline model ---
def forecast_step_single(ensemble=None, **icesee_kwargs):
    """Advance one member from model step k to k+1.

    Members are forced with the wrong ("nurged") terminus melt; the
    observations of the true state correct them at analysis steps.
    """
    k = int(icesee_kwargs.get("k", 0))
    icesee_kwargs["facemelt"] = facemelt_schedule(FACEMELT_WRONG_END, **icesee_kwargs)
    icesee_kwargs["transient"] = 1
    icesee_kwargs["tcurrent"] = k + 1
    return run_model(np.asarray(ensemble, dtype=float), **icesee_kwargs)


def _run_trajectory(statevec, facemelt_end, **icesee_kwargs):
    """Fill ``statevec`` (nd x nt+1) from the initial condition under the
    given terminus-melt ramp."""
    nd, nt = statevec.shape
    vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)

    state = np.zeros((nd, nt))
    state[:, 0] = initialize_model(**icesee_kwargs)
    icesee_kwargs["facemelt"] = facemelt_schedule(facemelt_end, **icesee_kwargs)
    icesee_kwargs["transient"] = 1
    for k in range(nt - 1):
        icesee_kwargs["tcurrent"] = k + 1
        for key, value in run_model(state[:, k], **icesee_kwargs).items():
            state[indx_map[key], k + 1] = value
    statevec[:, :] = state
    return {key: state[indx_map[key], :] for key in icesee_kwargs["vec_inputs"]}


# --- generate true state ---
def generate_true_state(**icesee_kwargs):
    """generate the true state of the model"""
    return _run_trajectory(icesee_kwargs["statevec_true"], FACEMELT_TRUE_END, **icesee_kwargs)


def generate_nurged_state(**icesee_kwargs):
    """generate the nurged state of the model"""
    return _run_trajectory(icesee_kwargs["statevec_nurged"], FACEMELT_WRONG_END, **icesee_kwargs)


# --- initialize the ensemble members ---
def initialize_ensemble(ens, **icesee_kwargs):
    """Every member starts from the steady initial state; ICESEE adds each
    member's perturbation."""
    vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)
    huxg_out0 = initialize_model(**icesee_kwargs)
    return {key: huxg_out0[indx_map[key]] for key in icesee_kwargs["vec_inputs"]}


# --- observation operator: every state entry is observed (as in mode 3) ---
def Obs_fun(virtual_obs=None):
    return np.asarray(virtual_obs).copy()


def JObs_fun(nd=None):
    return np.eye(nd)


def Cov_Obs_fun(sig_obs=None, nd=None, icesee_kwargs=None):
    # the serial runner passes hdim * num_state_vars (99), not the real state
    # size: the scalar xg block is shorter than h and u
    nd = int((icesee_kwargs or {}).get("nd", nd))
    return np.eye(nd) * (sig_obs ** 2)
