# ==============================================================================
# @des: This file contains run functions for icepack data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2024-11-4
# @author: Brian Kyanjo
# ==============================================================================

import os

import numpy as np
import tqdm
import h5py

# --- import run_simulation function from the available examples ---
from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_model import *
from ICESEE.config._utility_imports import icesee_get_index
from ICESEE.applications.icepack_model.icepack_utils._coordinates import (
    register_icepack_coordinate_provider,
)

register_icepack_coordinate_provider()


# --- Forecast step ---
def forecast_step_single(ensemble=None, **icesee_kwargs):
    """ensemble: packs the state variables:h,u,v,s,basal_melt_rate of a single ensemble member
                 where h is thickness, u and v are the x and y components
                 of the velocity field, s is ice elevation, basal_melt_rate is ice shelf basal melt
    Returns: ensemble: updated ensemble member
    """
    #  call the run_model fun to push the state forward in time
    return run_model(ensemble, **icesee_kwargs)


# --- generate true state ---
def generate_true_state(**icesee_kwargs):
    """generate the true state of the model"""

    # # unpack the **icesee_kwargs
    smb  = icesee_kwargs.get('smb', None)
    basal_melt_field = icesee_kwargs.get('basal_melt_field', None)
    bed  = icesee_kwargs.get('bed', None)
    dt = icesee_kwargs.get('dt', None)
    nt = icesee_kwargs.get('nt', None)
    A0  = icesee_kwargs.get('A0', None)
    beta0  = icesee_kwargs.get('beta0', None)
    Q  = icesee_kwargs.get('Q', None)
    V  = icesee_kwargs.get('V', None)
    h0 = icesee_kwargs.get('h0', None)
    u0 = icesee_kwargs.get('u', None)
    s0 = icesee_kwargs.get('s0', None)
    floating = icesee_kwargs.get('floating', None)
    grounded = icesee_kwargs.get('grounded', None)
    solver = icesee_kwargs.get('solver', None)
    statevec_true = icesee_kwargs["statevec_true"]
    save_steps = icesee_kwargs.get('save_steps', None)

    # # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)


    # # --- fetch the state variables ---
    statevec_true[indx_map["h"],0] = h0.dat.data_ro
    statevec_true[indx_map["u"],0] = u0.dat.data_ro[:,0]
    statevec_true[indx_map["v"],0] = u0.dat.data_ro[:,1]
    statevec_true[indx_map["s"],0] = s0.dat.data_ro
    # basal_melt_field is always part of the packed state vector regardless of the state/parameter bookkeeping split (see _icepack_model.py's 2026-09-28 reconciliation note).
    statevec_true[indx_map["basal_melt_field"],0] = basal_melt_field.dat.data_ro

    h = h0.copy(deepcopy=True)
    u = u0.copy(deepcopy=True)
    s = s0.copy(deepcopy=True)
    basal_melt_field = basal_melt_field.copy(deepcopy=True)

    # --- extract a profile of the flowline at the initial state ---
    h_profiles, s_profiles, valid_points, distances, bed_values = initial_flowline_profile(icesee_kwargs)

    #### DEBUGGING
    #print(f"\n INITIAL h_profiles_mean = {np.mean(h_profiles)}, s_profiles_mean = {np.mean(s_profiles)} \n")
    #############

    # --- step numbers at which to extract flowline profiles during the simulation -- 
    flowline_profile_steps = [t/dt for t in icesee_kwargs["save_steps"]]
    print(f"\n steps where profiles are sampled = {flowline_profile_steps} \n")

    # Path-isolation fix (2026-09-28): was a hardcoded, non-run-specific relative path -- two jobs sharing this working directory would collide on the same file. 
    # Now scoped under this run's own data_path, matching mode3_runner.py's existing pattern.
    _modelrun_datasets = icesee_kwargs.get("data_path") or "_modelrun_datasets"
    os.makedirs(_modelrun_datasets, exist_ok=True)
    # Earlier runs wrote this file without the .h5 extension; readers fall
    # back to that name (see read_results.ipynb).
    hs_files = f"{_modelrun_datasets}/hs_profiles_true.h5"

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
        print(f"Time step = {step}")        

        ### Conditionals for depth-dependent basal melt rate function
        ### Select forcing scenario between 1935 - 2017
        if step < (6/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_TRUE) # 1935 - 1941

        elif (6 / dt) <= step < (15/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_TRUE) # 1941 - 1950

        elif (15 / dt) <= step < (18/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_TRUE) # 1950 - 1953

        elif (18/ dt) <= step < (20/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_TRUE) # 1953 - 1955

        elif (20/ dt) <= step < (25/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_TRUE) # 1955 - 1960

        elif (25/ dt) <= step < (27/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_TRUE) # 1960 - 1962

        elif (27/ dt) <= step < (31/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_TRUE) # 1962 - 1966

        elif (31/ dt) <= step < (40/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_TRUE) # 1966 - 1975

        elif (40/ dt) <= step < (48/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_TRUE) # 1975 - 1983

        elif (48/dt) <= step < (50/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_TRUE) # 1983 - 1985

        elif (50/ dt) <= step < (59/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_TRUE) # 1985 - 1994

        elif (59/ dt) <= step < (64/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_TRUE) # 1994 - 2000

        elif (64/ dt) <= step < (69/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_TRUE) # 2000 - 2005

        elif (69/ dt) <= step < (76/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_TRUE) # 2005 - 2012

        elif (76/ dt) <= step:
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_TRUE) # 2012 - 2017


        #### DEBUGGING 
        #print(f"\n Time step = {step}, h_mean = {np.mean(h.dat.data_ro)}, u_mean = {np.mean(u.dat.data_ro[:,0])}, v_mean = {np.mean(u.dat.data_ro[:,1])}, smb_mean = {np.mean(smb.dat.data_ro)}, bmr_mean = {np.mean(basal_melt_field.dat.data_ro)} \n")
        ###############
        


        # call the ice stream model to update the state variables
        h, u, s, floating, grounded = Icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, icesee_kwargs)

        statevec_true[indx_map["h"],k+1] = h.dat.data_ro
        statevec_true[indx_map["u"],k+1] = u.dat.data_ro[:,0]
        statevec_true[indx_map["v"],k+1] = u.dat.data_ro[:,1]
        statevec_true[indx_map["s"],k+1] = s.dat.data_ro

        # basal_melt_field is always part of the packed state vector.
        statevec_true[indx_map["basal_melt_field"],k+1] = basal_melt_field.dat.data_ro



        if (kk <= len(flowline_profile_steps)-1):
            if (k == int(flowline_profile_steps[kk])):

                h_profiles, s_profiles = flowline_profile(h, s, valid_points)
                print(s_profiles.shape,h_profiles.shape,"\n")


                # Column 0 holds the initial profile; the kk-th requested
                # save step goes to column kk + 1.
                with h5py.File(hs_files, "a") as F:
                    F["h_profiles"][:,kk+1] = h_profiles
                    F["s_profiles"][:,kk+1] = s_profiles

                kk += 1

    # Fully parallel execution supplies an HDF5 dataset as the write target.
    # The complete trajectory is already on disk; slicing every variable here
    # would materialize the entire (potentially tens-of-GB) trajectory in RAM.
    if isinstance(statevec_true, h5py.Dataset):
        statevec_true.flush()
        return None

    updated_state = {}
    for key in icesee_kwargs["vec_inputs"]:
        updated_state[key] = statevec_true[indx_map[key], :]

    return updated_state




# --- initialize the ensemble members ---
def initialize_ensemble(ens, **icesee_kwargs):

    """initialize the ensemble members"""

    ### UNPACK THE **icesee_kwargs DICTIONARY
    smb  = icesee_kwargs.get('smb', None)
    basal_melt_field = icesee_kwargs.get('basal_melt_field', None)
    bed  = icesee_kwargs.get('bed', None)
    dt = icesee_kwargs.get('dt', None)
    #nt = icesee_kwargs.get('nt', None)
    #A0  = icesee_kwargs.get('A0', None)
    #beta0  = icesee_kwargs.get('beta0', None)
    Q  = icesee_kwargs.get('Q', None)
    #V  = icesee_kwargs.get('V', None)
    h0 = icesee_kwargs.get('h0', None)
    u0 = icesee_kwargs.get('u', None)
    s0 = icesee_kwargs.get('s0', None)
    #floating = icesee_kwargs.get('floating', None)
    #grounded = icesee_kwargs.get('grounded', None)
    solver = icesee_kwargs.get('solver', None)
    save_steps = icesee_kwargs.get('save_steps', None)


    ### INITIALIZE THE ENSEMBLE MEMBERS



    # CORRECTED (2026-09-28, second pass): joint_estimation means "are we performing joint state-parameter estimation," not "does basal_melt_field exist." 
    # This nudge initializes the ensemble's WRONG prior belief about a jointly-estimated parameter -- genuinely joint-estimation-algorithm-specific, so the gate is kept. 
    # Under the current config (joint_estimation=0; true-vs-wrong divergence instead comes from BasalMeltRate's EXPERIMENT_TRUE/EXPERIMENT_WRONG selector -- see _icepack_model.py), 
    # this gate is False and the nudge does not fire, so the two mechanisms never operate simultaneously in practice.
    if icesee_kwargs["joint_estimation"]:
        basal_melt_field_nudged = basal_melt_field.dat.data_ro + icesee_kwargs["wrong_basal_melt_field"]
        basal_melt_field = Function(Q)
        basal_melt_field.dat.data[:] = basal_melt_field_nudged


    h, u, s, floating, grounded = Icepack(solver, h0, u0, smb, basal_melt_field, bed, dt, h0, icesee_kwargs)

    # basal_melt_field is always part of the packed state (field existence is driven by the application's configured state definition, vec_inputs
    # -- never by joint_estimation).
    initialized_state = {'h': h.dat.data_ro,
                         'u': u.dat.data_ro[:,0],
                         'v': u.dat.data_ro[:,1],
                         's': s.dat.data_ro,
                         'basal_melt_field': basal_melt_field.dat.data_ro}

    # --- create file to save ensemble flowline profiles ---
    h_profiles, s_profiles, valid_points, distances, bed_values = initial_flowline_profile(icesee_kwargs)
    
    # Path-isolation fix (2026-09-28): see generate_true_state's hs_files.
    _modelrun_datasets = icesee_kwargs.get("data_path") or "_modelrun_datasets"
    os.makedirs(_modelrun_datasets, exist_ok=True)
    hs_ensemble_files = f"{_modelrun_datasets}/hs_ensemble_profiles{ens}.h5"

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
def generate_nurged_state(**icesee_kwargs):

    """generate the nudged state of the model"""
    
    # unpack the **kwargs
    smb  = icesee_kwargs.get('smb', None)
    basal_melt_field = icesee_kwargs.get('basal_melt_field', None)
    bed  = icesee_kwargs.get('bed', None)
    dt = icesee_kwargs.get('dt', None)
    nt = icesee_kwargs.get('nt', None)
    A0  = icesee_kwargs.get('A0', None)
    beta0  = icesee_kwargs.get('beta0', None)
    Q  = icesee_kwargs.get('Q', None)
    V  = icesee_kwargs.get('V', None)
    h0 = icesee_kwargs.get('h0', None)
    u0 = icesee_kwargs.get('u', None)
    s0 = icesee_kwargs.get('s0', None)
    floating = icesee_kwargs.get('floating', None)
    grounded = icesee_kwargs.get('grounded', None)
    solver = icesee_kwargs.get('solver', None)
    statevec_nurged = icesee_kwargs["statevec_nurged"]
    save_steps = icesee_kwargs.get('save_steps', None)
     

    # --- define the state variables list ---
    vec_inputs = icesee_kwargs["vec_inputs"]

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)

    statevec_nurged[indx_map["h"],0]   = h0.dat.data_ro[:,0]
    statevec_nurged[indx_map["u"],0]   = u0.dat.data_ro[:,0]
    statevec_nurged[indx_map["v"],0]   = u0.dat.data_ro[:,1]
    statevec_nurged[indx_map["s"],0]   = s0.dat.data_ro[:,0]

    # CORRECTED (2026-09-28, second pass): joint_estimation gates the joint-estimation-algorithm-specific nudge (initializing the ensemble's 
    # WRONG prior belief about a jointly-estimated parameter), not whether basal_melt_field exists in the packed state -- that is decided by the 
    # application's configured state definition (vec_inputs) and is therefore unconditional below. Under the current config (joint_estimation=0), 
    # this gate is False and the nudge does not fire; the true-vs-wrong divergence instead comes from the forecast loop's BasalMeltRate(..., experiment=EXPERIMENT_WRONG) 
    # calls, so the two mechanisms never operate simultaneously in practice.
    if icesee_kwargs["joint_estimation"]:
        basal_melt_field_nudged = basal_melt_field.dat.data_ro + icesee_kwargs["wrong_basal_melt_field"]
        basal_melt_field = Function(Q)
        basal_melt_field.dat.data[:] = basal_melt_field_nudged

    statevec_nurged[indx_map["basal_melt_field"],0] = basal_melt_field.dat.data_ro

    h = h0.copy(deepcopy=True)
    u = u0.copy(deepcopy=True)
    s = s0.copy(deepcopy=True)
    basal_melt_field = basal_melt_field.copy(deepcopy=True)
         

    # --- extract a profile of the flowline at the initial state ---
    h_profiles, s_profiles, valid_points, distances, bed_values = initial_flowline_profile(icesee_kwargs)

    # --- step numbers at which to extract flowline profiles during the simulation --
    flowline_profile_steps = [t/dt for t in icesee_kwargs["save_steps"]]

    # Path-isolation fix (2026-09-28): see generate_true_state's hs_files.
    _modelrun_datasets = icesee_kwargs.get("data_path") or "_modelrun_datasets"
    os.makedirs(_modelrun_datasets, exist_ok=True)
    hs_nudged_files = f"{_modelrun_datasets}/hs_profiles_wrong.h5"

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
        print(f"Time step = {step}") 

        ### Conditionals for depth-dependent basal melt rate function
        ### Select forcing scenario between 1935 - 2017
        if step < (6/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_WRONG) # 1935 - 1941

        elif (6 / dt) <= step < (15/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_WRONG) # 1941 - 1950

        elif (15 / dt) <= step < (18/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_WRONG) # 1950 - 1953

        elif (18/ dt) <= step < (20/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_WRONG) # 1953 - 1955

        elif (20/ dt) <= step < (25/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_WRONG) # 1955 - 1960

        elif (25/ dt) <= step < (27/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_WRONG) # 1960 - 1962

        elif (27/ dt) <= step < (31/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_WRONG) # 1962 - 1966

        elif (31/ dt) <= step < (40/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_WRONG) # 1966 - 1975

        elif (40/ dt) <= step < (48/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_WRONG) # 1975 - 1983

        elif (48/dt) <= step < (50/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_WRONG) # 1983 - 1985

        elif (50/ dt) <= step < (59/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_WRONG) # 1985 - 1994

        elif (59/ dt) <= step < (64/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_WRONG) # 1994 - 2000

        elif (64/ dt) <= step < (69/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_WRONG) # 2000 - 2005

        elif (69/ dt) <= step < (76/dt):
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='warm', experiment=EXPERIMENT_WRONG) # 2005 - 2012

        elif (76/ dt) <= step:
            basal_melt_field, melt_max = BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_WRONG) # 2012 - 2017


        
        # call the ice stream model to update the state variables
        h, u, s, floating, grounded = Icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, icesee_kwargs)

        statevec_nurged[indx_map["h"],k+1] = h.dat.data_ro
        statevec_nurged[indx_map["u"],k+1] = u.dat.data_ro[:,0]
        statevec_nurged[indx_map["v"],k+1] = u.dat.data_ro[:,1]
        statevec_nurged[indx_map["s"],k+1] = s.dat.data_ro


        # CORRECTED (2026-09-28, second pass): joint_estimation gates only
        # the joint-estimation-algorithm-specific re-nudge of this step's
        # basal_melt_field; the packing into statevec_nurged is
        # unconditional (field existence comes from vec_inputs, not
        # joint_estimation). Under the current config this gate is False,
        # so BasalMeltRate's own experiment=EXPERIMENT_WRONG output (the
        # constant-melt "wrong" forcing) is packed unmodified.
        if icesee_kwargs["joint_estimation"]:
            basal_melt_field_nudged = basal_melt_field.dat.data_ro + icesee_kwargs["wrong_basal_melt_field"]
            basal_melt_field = Function(Q)
            basal_melt_field.dat.data[:] = basal_melt_field_nudged

        statevec_nurged[indx_map["basal_melt_field"],k+1] = basal_melt_field.dat.data_ro
        
        
       # # -- saving flowline profile at certain steps
        if (kk <= len(flowline_profile_steps)-1):
            if (k == int(flowline_profile_steps[kk])):

                h_profiles, s_profiles = flowline_profile(h, s, valid_points)
                print(s_profiles.shape,h_profiles.shape,"\n")


                # Column 0 holds the initial profile; the kk-th requested
                # save step goes to column kk + 1.
                with h5py.File(hs_nudged_files, "a") as F:
                    F["h_nudged_profiles"][:,kk+1] = h_profiles
                    F["s_nudged_profiles"][:,kk+1] = s_profiles

                kk += 1

    # See generate_true_state: mode 2 writes this trajectory in place.
    if isinstance(statevec_nurged, h5py.Dataset):
        statevec_nurged.flush()
        return None

    updated_state = {}
    for key in icesee_kwargs["vec_inputs"]:
        updated_state[key] = statevec_nurged[indx_map[key], :]

    return updated_state
