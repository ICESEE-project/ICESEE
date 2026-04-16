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
    floating = kwargs.get('floating', None)
    grounded = kwargs.get('grounded', None)
    solver = kwargs.get('solver', None)
    statevec_true = kwargs["statevec_true"]
    save_steps = kwargs.get('save_steps', None)

    # # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(**kwargs)

    
    # # --- fetch the state variables ---
    statevec_true[indx_map["h"],0] = h0.dat.data_ro
    statevec_true[indx_map["u"],0] = u0.dat.data_ro[:,0]
    statevec_true[indx_map["v"],0] = u0.dat.data_ro[:,1]
    statevec_true[indx_map["s"],0] = s0.dat.data_ro


    # # add BMR field to EnKF state vector if joint estimation is enabled 
    if kwargs["joint_estimation"]:
        statevec_true[indx_map["basal_melt_field"],0] = basal_melt_field.dat.data_ro

    h = h0.copy(deepcopy=True)
    u = u0.copy(deepcopy=True)
    s = s0.copy(deepcopy=True)

    # --- extract a profile of the flowline at the initial state ---
    h_profiles, s_profiles, valid_points, distances, bed_values = initial_flowline_profile(kwargs)

    # --- step numbers at which to extract flowline profiles during the simulation -- 
    flowline_profile_steps = [t/dt for t in kwargs["save_steps"]]
    print(flowline_profile_steps)

    hs_files = f"_modelrun_datasets/hs_profiles_true"
    
    with h5py.File(hs_files, "w") as F:
        dataset_h = F.create_dataset("h_profiles", (len(valid_points), len(save_steps) + 1), dtype = "f8")
        dataset_s = F.create_dataset("s_profiles", (len(valid_points), len(save_steps) + 1), dtype = "f8")
        dataset_s[:,0] = s_profiles # save the initial surface elevation profile
        dataset_h[:,0] = h_profiles # save the initial thickness profile
        dataset_bed = F.create_dataset("bed_values", data = bed_values)
        dataset_distances = F.create_dataset("distances", data = distances)

    
   
    kk = 0
    
    # loop through each step (k = 0 to k = 100)
    for k in range(nt):

        #step = k * dt
        step = k
        

        ### Conditionals for depth-dependent basal melt rate function
        ### Select forcing scenario between 1935 - 2017 
        if step < (6/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 1935 - 1941
   
        elif (6 / dt) <= step < (15/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 1941 - 1950
    
        elif (15 / dt) <= step < (18/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 1950 - 1953
    
        elif (18/ dt) <= step < (20/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 1953 - 1955
        
        elif (20/ dt) <= step < (25/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 1955 - 1960
        
        elif (25/ dt) <= step < (27/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 1960 - 1962
        
        elif (27/ dt) <= step < (31/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 1962 - 1966
        
        elif (31/ dt) <= step < (40/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 1966 - 1975
        
        elif (40/ dt) <= step < (48/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 1975 - 1983
        
        elif (48/dt) <= step < (50/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 1983 - 1985
        
        elif (50/ dt) <= step < (59/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 1985 - 1994
        
        elif (59/ dt) <= step < (64/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 1994 - 2000
        
        elif (64/ dt) <= step < (69/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 2000 - 2005
        
        elif (69/ dt) <= step < (76/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 2005 - 2012
        
        elif (76/ dt) <= step:
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 2012 - 2017
        
        print(f"maximum_bmr=",{melt_max})
        #print("basal melt field=",{np.mean(basal_melt_field.dat.data_ro)})
        #print("floating", {np.mean(floating.dat.data_ro)})
        #print("thickness", {np.mean(h.dat.data_ro)})
        #print("sfc elevation", {np.mean(s.dat.data_ro)})
        # print(f"\ninside true state s_mean = {np.mean(s.dat.data_ro)}, h_mean = {np.mean(h.dat.data_ro)}, u_mean = {np.mean(u.dat.data_ro[:,0])}, v_mean = {np.mean(u.dat.data_ro[:,1])}\n")
        # exit(1)


        # call the ice stream model to update the state variables
        h, u, s = Icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, kwargs)
        

        statevec_true[indx_map["h"],k+1] = h.dat.data_ro
        statevec_true[indx_map["u"],k+1] = u.dat.data_ro[:,0]
        statevec_true[indx_map["v"],k+1] = u.dat.data_ro[:,1]
        statevec_true[indx_map["s"],k+1] = s.dat.data_ro



        # update the basal melt rate if joint estimation is enabled
        if kwargs["joint_estimation"]:
            statevec_true[indx_map["basal_melt_field"],k+1] = basal_melt_field.dat.data_ro
       

       
        if (kk <= len(flowline_profile_steps)-1):
            if (k == int(flowline_profile_steps[kk])):
                
                h_profiles, s_profiles = flowline_profile(h, s, valid_points)
                print(s_profiles.shape,h_profiles.shape,"\n")
                
                
                with h5py.File(hs_files, "a") as F:
                    F["h_profiles"][:,kk] = h_profiles
                    F["s_profiles"][:,kk] = s_profiles
                
                kk += 1

    updated_state = {}      
    for key in kwargs["vec_inputs"]:
        updated_state[key] = statevec_true[indx_map[key], :]

    return updated_state


# --- initialize the ensemble members ---
def initialize_ensemble(ens, **kwargs):
    
    """initialize the ensemble members"""

    ### unpack the **kwargs
    # smb  = kwargs.get('smb', None)
    # basal_melt_field = kwargs.get('basal_melt_field', None)
    # wrong_basal_melt_field = kwargs.get('"wrong_basal_melt_field"', None)
    # bed  = kwargs.get('bed', None)
    # dt = kwargs.get('dt', None)
    # nt = kwargs.get('nt', None)
    # A0  = kwargs.get('A0', None)
    # beta0  = kwargs.get('beta0', None)
    Q  = kwargs.get('Q', None)
    # V  = kwargs.get('V', None)
    # h0 = kwargs.get('h0', None)
    # u0 = kwargs.get('u', None)
    # solver = kwargs.get('solver', None)

    ### UNPACK THE **kwargs DICTIONARY
    smb  = kwargs.get('smb', None)
    basal_melt_field = kwargs.get('basal_melt_field', None)
    bed  = kwargs.get('bed', None)
    dt = kwargs.get('dt', None)
    #nt = kwargs.get('nt', None)
    #A0  = kwargs.get('A0', None)
    #beta0  = kwargs.get('beta0', None)
    #Q  = kwargs.get('Q', None)
    #V  = kwargs.get('V', None)
    h0 = kwargs.get('h0', None)
    u0 = kwargs.get('u', None)
    s0 = kwargs.get('s0', None)
    #floating = kwargs.get('floating', None)
    #grounded = kwargs.get('grounded', None)
    solver = kwargs.get('solver', None)
    save_steps = kwargs.get('save_steps', None)

    
    ### INITIALIZE THE ENSEMBLE MEMBERS
 

    h_perturbed = h0.dat.data_ro + (15 * np.random.normal(0, 1, h0.dat.data_ro.size))

    h_p = Function(Q)
    h_p.dat.data[:] = h_perturbed

    
    #s0 = icepack.compute_surface(thickness = h_p, bed = bed)


    # # -- update parameter if joint estimation is enabled
    if kwargs["joint_estimation"]:
        basal_melt_field_nudged = basal_melt_field.dat.data_ro + kwargs["wrong_basal_melt_field"]
        #basal_melt_field_nudged = basal_melt_field.dat.data_ro
        basal_melt_field = Function(Q)
        basal_melt_field.dat.data[:] = basal_melt_field_nudged

    h, u, s = Icepack(solver, h_p, u0, smb, basal_melt_field, bed, dt, h0, kwargs)

    initialized_state = {'h': h.dat.data_ro,
                         'u': u.dat.data_ro[:,0], 
                         'v': u.dat.data_ro[:,1],
                         's': s.dat.data_ro}
    
    # -- for joint estimation --
    if kwargs["joint_estimation"]:
        initialized_state['basal_melt_field'] =  basal_melt_field_nudged

    # --- create file to save ensemble flowline profiles ---
    h_profiles, s_profiles, valid_points, distances, bed_values = initial_flowline_profile(kwargs)
    hs_ensemble_files = f"_modelrun_datasets/hs_ensemble_profiles{ens}"

    with h5py.File(hs_ensemble_files, "w") as F:
        dataset_h_ensemble = F.create_dataset("h_profiles", (len(valid_points), len(save_steps) + 1), dtype = "f8")
        dataset_s_ensemble = F.create_dataset("s_profiles", (len(valid_points), len(save_steps) + 1), dtype = "f8")
        dataset_s_ensemble[:,0] = s_profiles # save the initial surface elevation profile
        dataset_h_ensemble[:,0] = h_profiles # save the initial thickness profile
        dataset_bed_ensemble = F.create_dataset("bed_values", data = bed_values)
        dataset_distances_ensemble = F.create_dataset("distances", data = distances)
        dataset_valid_points_ensemble = F.create_dataset("valid_points", data = valid_points)
    
    #print(f"\ndt = {dt}, h0_mean = {np.mean(h0.dat.data_ro)}, h_mean = {np.mean(initialized_state["h"])}, u_mean = {np.mean(u.dat.data_ro[:,0])}, v_mean = {np.mean(u.dat.data_ro[:,1])}\n")

    return initialized_state

# --- generate the nurged state ---
def generate_nurged_state(**kwargs):
    """generate the nudged state of the model"""
    
    #params = kwargs["params"]
    #nt = params["nt"] - 1

    # unpack the **kwargs
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
    #s0 = kwargs.get('s0', None)
    floating = kwargs.get('floating', None)
    grounded = kwargs.get('grounded', None)
    solver = kwargs.get('solver', None)
    save_steps = kwargs.get('save_steps', None)
     
    wrong_basal_melt_field = kwargs.get('wrong_basal_melt_field', None)
 

    statevec_nurged = kwargs["statevec_nurged"]

     # --- define the state variables list ---
    vec_inputs = kwargs["vec_inputs"]

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(statevec_nurged, **kwargs)

    #  create a bump -100 to 0
    # h_indx = int(np.ceil(nurged_entries+1))
    # hdim = vecs['h'].shape[0]
    #h_index_map = indx_map["h"]
    #hdim = indx_map["h"].shape[0]


    # if velocity is nurged, then run to get a solution to be used as am initial guess for velocity.
    # if u_nurge_ic != 0.0 or h_nurge_ic != 0.0:
    #     h_indx = int(np.ceil(nurged_entries_percentage*hdim+1))
   
    #     # u_indx = int(np.ceil(u_nurge_ic+1))
    #     u_indx = 1
    #     h_bump = np.linspace(-h_nurge_ic,0,h_indx)
    #     u_bump = np.linspace(-u_nurge_ic,0,h_indx)
    #     # h_bump = np.random.uniform(-h_nurge_ic,0,h_indx)
    #     # u_bump = np.random.uniform(-u_nurge_ic,0,h_indx)
    #     # print(f"hdim: {hdim}, h_indx: {h_indx}")
    #     # print(f"[Debug]: h_bump shape: {h_bump.shape} h0_index: {h0.dat.data_ro[:h_indx].shape}")
    #     h_with_bump = h_bump + h0.dat.data_ro[:h_indx]
    #     u_with_bump = u_bump + u0.dat.data_ro[:h_indx,0]
    #     v_with_bump = u_bump + u0.dat.data_ro[:h_indx,1]

    #     h_perturbed = np.concatenate((h_with_bump, h0.dat.data_ro[h_indx:]))
    #     u_perturbed = np.concatenate((u_with_bump, u0.dat.data_ro[h_indx:,0]))
    #     v_perturbed = np.concatenate((v_with_bump, u0.dat.data_ro[h_indx:,1]))

    #     h = Function(Q)
    #     u = Function(V)
    #     h.dat.data[:]   = h_perturbed
    #     u.dat.data[:,0] = u_perturbed
    #     u.dat.data[:,1] = v_perturbed
    #     h0 = h.copy(deepcopy=True)
    #     # call the solver
    #     h, u = Icepack(kwargs, solver, h, u, a, b, dt, h0)

    #     # update the nurged state with the solution
    #     h_perturbed = h.dat.data_ro
    #     u_perturbed = u.dat.data_ro[:,0]
    #     v_perturbed = u.dat.data_ro[:,1]
    # else: 
    #     h_perturbed = h0.dat.data_ro + np.random.normal(0, 0.1, h0.dat.data_ro.size)
    #     u_perturbed = u0.dat.data_ro[:,0]
    #     v_perturbed = u0.dat.data_ro[:,1]

    h_perturbed = h0.dat.data_ro + (15 * np.random.normal(0, 1, h0.dat.data_ro.size))

    h_p = Function(Q)
    h_p.dat.data[:] = h_perturbed

    
    s0 = icepack.compute_surface(thickness = h_p, bed = bed)

    s_perturbed = s0.dat.data_ro
    
    u_perturbed = u0.dat.data_ro[:,0]
    v_perturbed = u0.dat.data_ro[:,1]

    statevec_nurged[indx_map["h"],0]   = h_perturbed
    statevec_nurged[indx_map["u"],0]   = u_perturbed
    statevec_nurged[indx_map["v"],0]   = v_perturbed
    statevec_nurged[indx_map["s"],0]   = s_perturbed

    # # --  initialize the basal melt rate field if joint estimation is enabled at the initial time step
    # if kwargs["joint_estimation"]:
        
    #    basal_melt_field = basal_melt_field.dat.data_ro
    #    basal_melt_field_nudged = basal_melt_field + kwargs["wrong_basal_melt_field"]
    #    statevec_nurged[indx_map["basal_melt_field"],0] = basal_melt_field_nudged

    # # -- update parameter if joint estimation is enabled
    if kwargs["joint_estimation"]:
        basal_melt_field_nudged = basal_melt_field.dat.data_ro + kwargs["wrong_basal_melt_field"]
        basal_melt_field = Function(Q)
        basal_melt_field.dat.data[:] = basal_melt_field_nudged
        
        statevec_nurged[indx_map["basal_melt_field"],0] = basal_melt_field_nudged

    h = Function(Q)
    s = Function(Q)
    u = Function(V)
    h.dat.data[:] = h_perturbed
    u.dat.data[:,0] = u_perturbed
    u.dat.data[:,1] = v_perturbed
    s.dat.data[:] = s_perturbed
    

    # --- extract a profile of the flowline at the initial state ---
    h_profiles, s_profiles, valid_points, distances, bed_values = initial_flowline_profile(kwargs)

    # --- step numbers at which to extract flowline profiles during the simulation -- 
    flowline_profile_steps = [t/dt for t in kwargs["save_steps"]]

    hs_nudged_files = f"_modelrun_datasets/hs_profiles_wrong"

    with h5py.File(hs_nudged_files, "w") as F:
        dataset_h = F.create_dataset("h_nudged_profiles", (len(valid_points), len(save_steps) + 1), dtype = "f8")
        dataset_s = F.create_dataset("s_nudged_profiles", (len(valid_points), len(save_steps) + 1), dtype = "f8")
        dataset_s[:,0] = s_profiles # save the initial surface elevation profile
        dataset_h[:,0] = h_profiles # save the initial thickness profile
        dataset_bed = F.create_dataset("bed_values", data = bed_values)
        dataset_distances = F.create_dataset("distances", data = distances)

    # # -- index for saving flowline profiles
    kk = 0
    

    # loop through each step (k = 0 to k = 100)
    for k in range(nt):

        step = k

        ### Conditionals for depth-dependent basal melt rate function
        ### Select forcing scenario between 1935 - 2017 
        if step < (6/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 1935 - 1941
   
        elif (6 / dt) <= step < (15/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 1941 - 1950
    
        elif (15 / dt) <= step < (18/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 1950 - 1953
    
        elif (18/ dt) <= step < (20/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 1953 - 1955
        
        elif (20/ dt) <= step < (25/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 1955 - 1960
        
        elif (25/ dt) <= step < (27/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 1960 - 1962
        
        elif (27/ dt) <= step < (31/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 1962 - 1966
        
        elif (31/ dt) <= step < (40/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 1966 - 1975
        
        elif (40/ dt) <= step < (48/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 1975 - 1983
        
        elif (48/dt) <= step < (50/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 1983 - 1985
        
        elif (50/ dt) <= step < (59/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 1985 - 1994
        
        elif (59/ dt) <= step < (64/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 1994 - 2000
        
        elif (64/ dt) <= step < (69/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 2000 - 2005
        
        elif (69/ dt) <= step < (76/dt):
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='warm') # 2005 - 2012
        
        elif (76/ dt) <= step:
            basal_melt_field, melt_max = BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control') # 2012 - 2017
        
        # # -- DEBUGGING
        #print(f"maximum_bmr=",{melt_max})

        
        
        
        # call the ice stream model to update the state variables
        h, u, s = Icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, kwargs)

        statevec_nurged[indx_map["h"],k+1] = h.dat.data_ro
        statevec_nurged[indx_map["u"],k+1] = u.dat.data_ro[:,0]
        statevec_nurged[indx_map["v"],k+1] = u.dat.data_ro[:,1]
        statevec_nurged[indx_map["s"],k+1] = s.dat.data_ro

        # # -- update parameter if joint estimation is enabled
        if kwargs["joint_estimation"]:
            basal_melt_field_nudged = basal_melt_field.dat.data_ro + kwargs["wrong_basal_melt_field"]
            basal_melt_field = Function(Q)
            basal_melt_field.dat.data[:] = basal_melt_field_nudged
            
            statevec_nurged[indx_map["basal_melt_field"],k+1] = basal_melt_field_nudged
        
       # # -- saving flowline profile at certain steps
        if (kk <= len(flowline_profile_steps)-1):
            if (k == int(flowline_profile_steps[kk])):
                
                h_profiles, s_profiles = flowline_profile(h, s, valid_points)
                print(s_profiles.shape,h_profiles.shape,"\n")
                
                
                with h5py.File(hs_nudged_files, "a") as F:
                    F["h_nudged_profiles"][:,kk] = h_profiles
                    F["s_nudged_profiles"][:,kk] = s_profiles
                
                kk += 1

    updated_state = {}      
    for key in kwargs["vec_inputs"]:
        updated_state[key] = statevec_nurged[indx_map[key], :]

    return updated_state