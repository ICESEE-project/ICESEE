# ==============================================================================
# @des: This file contains run functions for icepack data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2024-11-4
# @author: Brian Kyanjo
# ==============================================================================

import numpy as np

# --- import run_simulation function from the available examples ---
from ICESEE.applications.icepack_model.examples.synthetic_ice_stream._icepack_model import *
from ICESEE.config._utility_imports import icesee_get_index


# --- Forecast step ---
def forecast_step_single(ensemble=None, **kwargs):
    """ensemble: packs the state variables:h,u,v of a single ensemble member
                 where h is thickness, u and v are the x and y components 
                 of the velocity field
    Returns: ensemble: updated ensemble member
    """
    #  call the run_model fun to push the state forward in time
    return run_model(ensemble, **kwargs)

# --- Background step ---
def background_step(**kwargs):
    """ computes the background state of the model
    Args:
        k: time step index
        background_vec: background state of the model
        hdim: dimension of the state variables
    Returns:
        background_vec: updated background state of the model
    """
    # unpack the **kwargs
    # a = kwargs.get('a', None)
    b = kwargs.get('b', None)
    dt = kwargs.get('dt', None)
    h0 = kwargs.get('h0', None)
    A = kwargs.get('A', None)
    C = kwargs.get('C', None)
    Q = kwargs.get('Q', None)
    V = kwargs.get('V', None)
    a_nuged = kwargs.get('a_nuged', None)
    solver = kwargs.get('solver', None)
    background_vec = kwargs.get('background_vec', None)

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(background_vec, **kwargs)

    # fetch the state variables
    hb = Function(Q)
    hb.dat.data[:]   = background_vec[indx_map["h"]]
    ub = Function(V)
    ub.dat.data[:,0] = background_vec[indx_map["u"]]
    ub.dat.data[:,1] = background_vec[indx_map["v"]]

    # call the ice stream model to update the state variables
    hb, ub = Icepack(solver, hb, ub,  a_nuged, b, dt, h0, fluidity = A, friction = C)

    # update the background state at the next time step
    updated_state = {'h': hb.dat.data_ro,
                    'u': ub.dat.data_ro[:,0],
                    'v': ub.dat.data_ro[:,1]}

    if kwargs["joint_estimation"]:
        updated_state['smb'] = a_nuged.dat.data_ro

    return updated_state

# --- generate true state ---
def generate_true_state(**kwargs):
    """generate the true state of the model"""

    # unpack the **kwargs
    a  = kwargs.get('a', None)
    b  = kwargs.get('b', None)
    dt = kwargs.get('dt', None)
    A  = kwargs.get('A', None)
    C  = kwargs.get('C', None)
    Q  = kwargs.get('Q', None)
    V  = kwargs.get('V', None)
    h0 = kwargs.get('h0', None)
    u0 = kwargs.get('u0', None)
    solver = kwargs.get('solver', None)
    statevec_true = kwargs["statevec_true"]
    params = kwargs["params"]

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(statevec_true, **kwargs)
    
    # --- fetch the state variables ---
    statevec_true[indx_map["h"],0] = h0.dat.data_ro
    statevec_true[indx_map["u"],0] = u0.dat.data_ro[:,0]
    statevec_true[indx_map["v"],0] = u0.dat.data_ro[:,1]

    # intialize the accumulation rate if joint estimation is enabled at the initial time step
    if kwargs["joint_estimation"]:
        statevec_true[indx_map["smb"],0] = a.dat.data_ro

    h = h0.copy(deepcopy=True)
    u = u0.copy(deepcopy=True)
    for k in range(params['nt']):
        # call the ice stream model to update the state variables
        h, u = Icepack(solver, h, u, a, b, dt, h0, fluidity = A, friction = C)

        statevec_true[indx_map["h"],k+1] = h.dat.data_ro
        statevec_true[indx_map["u"],k+1] = u.dat.data_ro[:,0]
        statevec_true[indx_map["v"],k+1] = u.dat.data_ro[:,1]

        # update the accumulation rate if joint estimation is enabled
        if kwargs["joint_estimation"]:
            statevec_true[indx_map["smb"],k+1] = a.dat.data_ro

    # update_state = {'h': statevec_true[indx_map["h"],:], 
    #                 'u': statevec_true[indx_map["u"],:], 
    #                 'v': statevec_true[indx_map["v"],:]}
    # # -- for joint estimation --
    # if kwargs["joint_estimation"]:
    #     update_state['smb'] = statevec_true[indx_map["smb"],:]
    # return update_state

# --- initialize the ensemble members ---
def initialize_ensemble(ens, **kwargs):
    
    """initialize the ensemble members"""

    # unpack the **kwargs
    h0 = kwargs.get('h0', None)
    u0 = kwargs.get('u0', None)
    params = kwargs["params"]
    # a  = kwargs.get('a', None)
    b  = kwargs.get('b', None)
    dt = kwargs.get('dt', None)
    A  = kwargs.get('A', None)
    C  = kwargs.get('C', None)
    Q  = kwargs.get('Q', None)
    V  = kwargs.get('V', None)
    a_nuged = kwargs.get('a_nuged', None)
    solver = kwargs.get('solver', None)
    h_nurge_ic      = kwargs.get('h_nurge_ic', None)
    u_nurge_ic      = kwargs.get('u_nurge_ic', None)
    nurged_entries_percentage  = kwargs.get('nurged_entries_percentage', None)
    statevec_ens    = kwargs["statevec_ens"]
    x = kwargs.get('x', None)
    Lx = kwargs.get('Lx', None)


    # initialize the ensemble members
    # hdim = vecs['h'].shape[0]
    hdim = h0.dat.data_ro.size
    # h_indx = int(np.ceil(nurged_entries_percentage*hdim+1))

    # # # create a bump -100 to 0
    # h_bump = np.linspace(-h_nurge_ic,0,h_indx)
    # h_with_bump = h_bump + h0.dat.data_ro[:h_indx]
    # h_perturbed = np.concatenate((h_with_bump, h0.dat.data_ro[h_indx:]))
    # statevec_ens[:hdim,ens] = h_perturbed 
    # h_perturbed = h0.dat.data_ro

    if u_nurge_ic != 0 or h_nurge_ic != 0:
        h_indx = int(np.ceil(nurged_entries_percentage*hdim+1))
   
        # u_indx = int(np.ceil(u_nurge_ic+1))
        u_indx = 1
        h_bump = np.linspace(-h_nurge_ic,0,h_indx)
        u_bump = np.linspace(-u_nurge_ic,0,h_indx)
        # h_bump = np.random.uniform(-h_nurge_ic,0,h_indx)
        # u_bump = np.random.uniform(-u_nurge_ic,0,h_indx)
        # print(f"hdim: {hdim}, h_indx: {h_indx}")
        # print(f"[Debug]: h_bump shape: {h_bump.shape} h0_index: {h0.dat.data_ro[:h_indx].shape}")
        h_with_bump = h_bump + h0.dat.data_ro[:h_indx]
        u_with_bump = u_bump + u0.dat.data_ro[:h_indx,0]
        v_with_bump = u_bump + u0.dat.data_ro[:h_indx,1]

        h_perturbed = np.concatenate((h_with_bump, h0.dat.data_ro[h_indx:]))
        u_perturbed = np.concatenate((u_with_bump, u0.dat.data_ro[h_indx:,0]))
        v_perturbed = np.concatenate((v_with_bump, u0.dat.data_ro[h_indx:,1]))

        h = Function(Q)
        u = Function(V)
        h.dat.data[:]   = h_perturbed
        u.dat.data[:,0] = u_perturbed
        u.dat.data[:,1] = v_perturbed
        h0 = h.copy(deepcopy=True)
        # call the solver
        h, u = Icepack(solver, h, u,  a_nuged, b, dt, h0, fluidity = A, friction = C)

        # update the nurged state with the solution
        h_perturbed = h.dat.data_ro
        u_perturbed = u0.dat.data_ro[:,0]
        v_perturbed = u0.dat.data_ro[:,1]
    else: 
        h_perturbed = h0.dat.data_ro + np.random.normal(0, 0.1, h0.dat.data_ro.size)
        u_perturbed = u0.dat.data_ro[:,0]
        v_perturbed = u0.dat.data_ro[:,1]

    initialized_state = {'h': h_perturbed, 
                         'u': u.dat.data_ro[:,0], 
                         'v': u.dat.data_ro[:,1]}
    
    # -- for joint estimation --
    if kwargs["joint_estimation"]:
        initialized_state['smb'] =  a_nuged.dat.data_ro
       
    return initialized_state

# --- generate the nurged state ---
def generate_nurged_state(**kwargs):
    """generate the nurged state of the model"""
    
    params = kwargs["params"]
    nt = params["nt"] - 1

    # unpack the **kwargs
    a = kwargs.get('a_p', None)
    t = kwargs.get('t', None)
    x = kwargs.get('x', None)
    Lx = kwargs.get('Lx', None)
    b = kwargs.get('b', None)
    dt = kwargs.get('dt', None)
    A = kwargs.get('A', None)
    C = kwargs.get('C', None)
    Q = kwargs.get('Q', None)
    V = kwargs.get('V', None)
    h0 = kwargs.get('h0', None)
    u0 = kwargs.get('u0', None)
    solver = kwargs.get('solver', None)
    a_in_p = kwargs.get('a_in_p', None)
    da_p = kwargs.get('da_p', None)
    da = kwargs.get('da', None)
    h_nurge_ic      = kwargs.get('h_nurge_ic', None)
    u_nurge_ic      = kwargs.get('u_nurge_ic', None)
    nurged_entries_percentage  = kwargs.get('nurged_entries_percentage', None)

    statevec_nurged = kwargs["statevec_nurged"]

     # --- define the state variables list ---
    vec_inputs = kwargs["vec_inputs"]

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(statevec_nurged, **kwargs)

    #  create a bump -100 to 0
    # h_indx = int(np.ceil(nurged_entries+1))
    # hdim = vecs['h'].shape[0]
    h_index_map = indx_map["h"]
    hdim = indx_map["h"].shape[0]

     # intialize the accumulation rate if joint estimation is enabled at the initial time step
    if kwargs["joint_estimation"]:
        tnur = np.linspace(.1, 2, nt)
        # aa   = a_in_p*(np.sin(tnur[0]) + 1)
        # daa  = da_p*(np.sin(tnur[0]) + 1)
        aa = a_in_p
        daa = da_p
        a_in = firedrake.Constant(aa)
        da_  = firedrake.Constant(daa)
        a    = firedrake.interpolate(a_in + da_ * x / Lx, Q)
        statevec_nurged[indx_map["smb"],0] = a.dat.data_ro

    # if velocity is nurged, then run to get a solution to be used as am initial guess for velocity.
    if u_nurge_ic != 0.0 or h_nurge_ic != 0.0:
        h_indx = int(np.ceil(nurged_entries_percentage*hdim+1))
   
        # u_indx = int(np.ceil(u_nurge_ic+1))
        u_indx = 1
        h_bump = np.linspace(-h_nurge_ic,0,h_indx)
        u_bump = np.linspace(-u_nurge_ic,0,h_indx)
        # h_bump = np.random.uniform(-h_nurge_ic,0,h_indx)
        # u_bump = np.random.uniform(-u_nurge_ic,0,h_indx)
        # print(f"hdim: {hdim}, h_indx: {h_indx}")
        # print(f"[Debug]: h_bump shape: {h_bump.shape} h0_index: {h0.dat.data_ro[:h_indx].shape}")
        h_with_bump = h_bump + h0.dat.data_ro[:h_indx]
        u_with_bump = u_bump + u0.dat.data_ro[:h_indx,0]
        v_with_bump = u_bump + u0.dat.data_ro[:h_indx,1]

        h_perturbed = np.concatenate((h_with_bump, h0.dat.data_ro[h_indx:]))
        u_perturbed = np.concatenate((u_with_bump, u0.dat.data_ro[h_indx:,0]))
        v_perturbed = np.concatenate((v_with_bump, u0.dat.data_ro[h_indx:,1]))

        h = Function(Q)
        u = Function(V)
        h.dat.data[:]   = h_perturbed
        u.dat.data[:,0] = u_perturbed
        u.dat.data[:,1] = v_perturbed
        h0 = h.copy(deepcopy=True)
        # call the solver
        h, u = Icepack(solver, h, u, a, b, dt, h0, fluidity = A, friction = C)

        # update the nurged state with the solution
        h_perturbed = h.dat.data_ro
        u_perturbed = u.dat.data_ro[:,0]
        v_perturbed = u.dat.data_ro[:,1]
    else: 
        h_perturbed = h0.dat.data_ro + np.random.normal(0, 0.1, h0.dat.data_ro.size)
        u_perturbed = u0.dat.data_ro[:,0]
        v_perturbed = u0.dat.data_ro[:,1]

    statevec_nurged[indx_map["h"],0]   = h_perturbed
    statevec_nurged[indx_map["u"],0]   = u_perturbed
    statevec_nurged[indx_map["v"],0]   = v_perturbed

    h = Function(Q)
    u = Function(V)
    h.dat.data[:] = h_perturbed
    u.dat.data[:,0] = u_perturbed
    u.dat.data[:,1] = v_perturbed
    h0 = h0.copy(deepcopy=True)
    

    for k in range(params['nt']):
        # aa   = a_in_p*(np.sin(tnur[k]) + 1)
        # daa  = da_p*(np.sin(tnur[k]) + 1)
        
        # call the ice stream model to update the state variables
        h, u = Icepack(solver, h, u, a, b, dt, h0, fluidity = A, friction = C)

        statevec_nurged[indx_map["h"],k+1] = h.dat.data_ro
        statevec_nurged[indx_map["u"],k+1] = u.dat.data_ro[:,0]
        statevec_nurged[indx_map["v"],k+1] = u.dat.data_ro[:,1]

        if kwargs["joint_estimation"]:
            # aa = a_in_p
            # daa = da_p
            a_in = firedrake.Constant(aa)
            da_  = firedrake.Constant(daa)
            a    = firedrake.interpolate(a_in + da_ * x / Lx, Q)
            statevec_nurged[indx_map["smb"],k+1] = a.dat.data_ro

    # return statevec_nurged