# ==============================================================================
# @des: This file contains run functions for icepack data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2024-11-4
# @author: Brian Kyanjo
# ==============================================================================

import numpy as np
import tqdm
import h5py

# --- import run_simulation function from the available examples ---
from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_model import *
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


# --- generate true state ---
def generate_true_state(**kwargs):
    """generate the true state of the model"""
    
    # # unpack the **kwargs
    smb  = kwargs.get('smb', None)
    basal_melt_field = kwargs.get('basal_melt_field', None)
    bed  = kwargs.get('bed', None)
    dt = kwargs.get('dt', None)
    nt = kwargs.get('nt', None)
    A0  = kwargs.get('A0', None)
    beta0  = kwargs.get('beta0', None)
    Q  = kwargs.get('Q', None)
    V  = kwargs.get('V', None)
    h0 = kwargs.get('h0', None)
    u0 = kwargs.get('u', None)
    s0 = kwargs.get('s0', None)
    solver = kwargs.get('solver', None)
    statevec_true = kwargs["statevec_true"]
    w2i     = float(kwargs.get('water_to_ice', 1.0))    # ratio of the density of water to ice
    save_steps = kwargs.get('save_steps', None)

    # # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(**kwargs)

    # ---- net accumulation used by prognostic step ------
    a_net = icepack.interpolate((smb - basal_melt_field) * w2i, Q)

    
    # # --- fetch the state variables ---
    statevec_true[indx_map["h"],0] = h0.dat.data_ro
    statevec_true[indx_map["u"],0] = u0.dat.data_ro[:,0]
    statevec_true[indx_map["v"],0] = u0.dat.data_ro[:,1]
    statevec_true[indx_map["s"],0] = s0.dat.data_ro


    # # intialize the accumulation rate if joint estimation is enabled at the initial time step
    if kwargs["joint_estimation"]:
        statevec_true[indx_map["basal_melt_field"],0] = basal_melt_field.dat.data_ro

    #h = h0.copy(deepcopy=True)
    #u = u0.copy(deepcopy=True) 
    h = h0.copy(deepcopy=True)
    u = u0.copy(deepcopy=True)
    s = s0.copy(deepcopy=True)

    # --- extract a profile of the flowline at the initial state ---
    h_profiles, s_profiles, profile_points, distances, bed_values = initial_flowline_profile(kwargs)

    # --- steps at which to extract flowline profiles during the simulation -- 
    flowline_profile_steps = {t/dt for t in kwargs["save_steps"]}

    hs_files = f"_modelrun_datasets/hs_profiles_true"
    
    with h5py.File(hs_files, "w") as F:
        dataset_h = F.create_dataset("h_profiles", (len(profile_points), len(save_steps)), dtype = "f8")
        dataset_s = F.create_dataset("s_profiles", (len(profile_points), len(save_steps)), dtype = "f8")
    
    for k in range(kwargs['nt']):
        
        # call the ice stream model to update the state variables
        h, u, s = Icepack(solver, h, u, a_net, bed, dt, h0, kwargs)
        #print("\ninside true state\n")
        #exit(0)

        statevec_true[indx_map["h"],k+1] = h.dat.data_ro
        statevec_true[indx_map["u"],k+1] = u.dat.data_ro[:,0]
        statevec_true[indx_map["v"],k+1] = u.dat.data_ro[:,1]
        statevec_true[indx_map["s"],k+1] = s.dat.data_ro

        # update the basal melt rate if joint estimation is enabled
        if kwargs["joint_estimation"]:
            statevec_true[indx_map["basal_melt_field"],k+1] = basal_melt_field.dat.data_ro
        
        if k in flowline_profile_steps:
            
            h_profiles, s_profiles = flowline_profile(kwargs, h_profiles, s_profiles, profile_points)
            
            
            with h5py.File(hs_files, "a") as F:
                F["h_profiles"] = h_profiles
                F["s_profiles"] = s_profiles

        updated_state = {} 
        
        for key in kwargs["vec_inputs"]:
            updated_state[key] = statevec_true[indx_map[key], :]

        return updated_state


# --- initialize the ensemble members ---
def initialize_ensemble(ens, **kwargs):
    
    """initialize the ensemble members"""

    # unpack the **kwargs
    smb  = kwargs.get('smb', None)
    basal_melt_field = kwargs.get('basal_melt_field', None)
    wrong_basal_melt_field = kwargs.get('"wrong_basal_melt_field"', None)
    bed  = kwargs.get('bed', None)
    dt = kwargs.get('dt', None)
    nt = kwargs.get('nt', None)
    A0  = kwargs.get('A0', None)
    beta0  = kwargs.get('beta0', None)
    Q  = kwargs.get('Q', None)
    V  = kwargs.get('V', None)
    h0 = kwargs.get('h0', None)
    u0 = kwargs.get('u', None)
    solver = kwargs.get('solver', None)
    


    # initialize the ensemble members
    basal_melt_field = basal_melt_field.dat.data_ro
    basal_melt_field_nudged = basal_melt_field + kwargs["wrong_basal_melt_field"]
    
    h = checkThickness(kwargs) 
    s = icepack.compute_surface(thickness=h, bed=bed)

    initialized_state = {'h': h.dat.data_ro,
                         'u': u0.dat.data_ro[:,0], 
                         'v': u0.dat.data_ro[:,1],
                         's': s.dat.data_ro}
    
    # -- for joint estimation --
    if kwargs["joint_estimation"]:
        initialized_state['basal_melt_field'] =  basal_melt_field_nudged
       
    return initialized_state

# --- generate the nurged state ---
def generate_nurged_state(**kwargs):
    """generate the nudged state of the model"""
    
    params = kwargs["params"]
    nt = params["nt"] - 1

    # unpack the **kwargs
    smb = kwargs.get('smb', None)
    #a = kwargs.get('a', None)
    #t = kwargs.get('t', None)
    #x = kwargs.get('x', None)
    #Lx = kwargs.get('Lx', None)
    bed = kwargs.get('bed', None)
    dt = kwargs.get('dt', None)
    A0 = kwargs.get('A0', None)
    beta0 = kwargs.get('beta0', None)
    Q = kwargs.get('Q', None)
    V = kwargs.get('V', None)
    h0 = kwargs.get('h0', None)
    u0 = kwargs.get('u0', None)
    solver = kwargs.get('solver', None)
    #a_in_p = kwargs.get('a_in_p', None)
    wrong_basal_melt_field = kwargs.get('wrong_basal_melt_field', None),
    basal_melt_field = kwargs.get('basal_melt_field'),
    #da_p = kwargs.get('da_p', None)
    #da = kwargs.get('da', None)

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

     # intialize the basal melt rate field if joint estimation is enabled at the initial time step
    if kwargs["joint_estimation"]:
        #tnur = np.linspace(.1, 2, nt)
        # aa   = a_in_p*(np.sin(tnur[0]) + 1)
        # daa  = da_p*(np.sin(tnur[0]) + 1)
        #aa = a_in_p
        #daa = da_p
        #a_in = firedrake.Constant(aa)
        #da_  = firedrake.Constant(daa)
        #a  = firedrake.interpolate(a_in + da_ * x / Lx, Q)
        #statevec_nurged[indx_map["smb"],0] = a.dat.data_ro
       basal_melt_field = basal_melt_field.dat.data_ro
       basal_melt_field_nudged = basal_melt_field + kwargs["wrong_basal_melt_field"]
       statevec_nurged[indx_map["basal_melt_field"],0] = basal_melt_field_nudged


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
        h, u = Icepack(kwargs, solver, h, u, a, b, dt, h0)

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
        h, u = Icepack(solver, h, u, a, b, dt, h0, fluidity = A0, friction = beta0)

        statevec_nurged[indx_map["h"],k+1] = h.dat.data_ro
        statevec_nurged[indx_map["u"],k+1] = u.dat.data_ro[:,0]
        statevec_nurged[indx_map["v"],k+1] = u.dat.data_ro[:,1]

        if kwargs["joint_estimation"]:
            # aa = a_in_p
            # daa = da_p
            #a_in = firedrake.Constant(aa)
            #da_  = firedrake.Constant(daa)
            #a    = firedrake.interpolate(a_in + da_ * x / Lx, Q)
            #statevec_nurged[indx_map["smb"],k+1] = a.dat.data_ro
            basal_melt_field = basal_melt_field.dat.data_ro
            basal_melt_field_nudged = basal_melt_field + kwargs["wrong_basal_melt_field"]
            statevec_nurged[indx_map["basal_melt_field"],0] = basal_melt_field_nudged

    # return statevec_nurged