
# ==============================================================================
# @des: This file contains run functions for icepack data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2024-11-4
# @author: Brian Kyanjo
# ==============================================================================

# --- python imports ---
import sys
import os
import h5py

os.environ["OMP_NUM_THREADS"] = "1"

# --- import model functions --- 

import ICESEE.applications.icepack_model.examples.idealized_pig.modelfunc as mf


# ---- firedrake imports ----
import firedrake
from firedrake import Constant, interpolate, assemble, inner, dx, ds, conditional, Function
from firedrake.petsc import PETSc


# ---- icepack imports ----

import icepack
import icepack.models.friction
from icepack.constants import ice_density as rhoI, weertman_sliding_law as m, glen_flow_law as n
import icepack.models


# --- miscellaneous imports -- 

from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import math
import yaml
from mpl_toolkits.axes_grid1 import make_axes_locatable
import tqdm
import pandas as pd
import rasterio
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


# --- Utility imports ---
from ICESEE.config._utility_imports import icesee_get_index


# ---- model initial state ---

def initialState(h0, s0, u0, zb, grounded0, floating0, Q):
    
    
    h, s = h0.copy(deepcopy=True), s0.copy(deepcopy=True)
    hLast = h0.copy(deepcopy=True)
    u = u0.copy(deepcopy=True)
    zF = mf.flotationHeight(zb, Q)
    grounded = grounded0.copy(deepcopy=True)
    floating = floating0.copy(deepcopy=True)

    return h, hLast, s, u, zF, grounded, floating




# ---- initial mesh ---

def initializeMesh(**kwargs):
    
    initFile = kwargs["initFile"]
    meshFile = kwargs["meshFile"]
    meshI = mf.getMeshFromCheckPoint(initFile,kwargs)
    mesh, Q, V, meshOpts = \
        mf.setupMesh(meshFile, degree=1,
                     meshOversample=2,
                     newMesh=meshI)

    return mesh, meshOpts, Q, V




# ---- initializing the run ---- 

def initializeRun(kwargs, forward_solver, mesh, Q, V):


    with firedrake.CheckpointFile(kwargs["initFile"],'r') as checkpoint:
        velocity = checkpoint.load_function(mesh, "velocity", idx=20000) # idx index set to the LAST time step of the steady-state run (20,000 for 1000 year run)
        h0 = checkpoint.load_function(mesh, "thickness", idx=20000)
        s0 = checkpoint.load_function(mesh, "surface", idx=20000)
        bed = checkpoint.load_function(mesh, "bed")
        grounded0 = checkpoint.load_function(mesh, "grounded")
        floating0 = checkpoint.load_function(mesh, "floating")
        A0 = checkpoint.load_function(mesh, "fluidity")
        beta0 = checkpoint.load_function(mesh, "extended_beta")

    """ find initial velocity """
    uThresh = firedrake.Constant(kwargs['uThresh'])
    u0 = forward_solver.diagnostic_solve(velocity=velocity, thickness=h0,
                                                surface=s0,
                                                beta=beta0, fluidity=A0, 
                                                grounded=grounded0,
                                                floating=floating0,
                                                uThresh=uThresh)
    
    """ define initial state """
    h, hLast, s, u, zF, grounded, floating = initialState(h0, s0, u0, bed, grounded0, floating0, Q)

    """ define smb and melt """
    smb = readSMB(kwargs,Q)
    
    return h, h0, s, s0, u, bed, zF, grounded, floating, A0, beta0, smb




# ---- basal friction model --- 
def schoofFriction(velocity, grounded, beta, uThresh):
    
    
    C = grounded * beta**2
    mExp = (1./m + 1.)
    U = firedrake.sqrt(firedrake.inner(velocity, velocity))

    return C * ((uThresh**mExp + U**mExp)**(m/(m+1.)) - uThresh)



# --- rheology model --- 

def regViscosity(**kwargs):
  
    u = kwargs["velocity"]
    h = kwargs["thickness"]
    A = kwargs["fluidity"]

    return icepack.models.viscosity.viscosity_depth_averaged(velocity=u, thickness=h, fluidity= A)




# --- read in the SMB file --- 

def readSMB(kwargs, Q):

    
    if not os.path.exists(kwargs["SMBFile"]):
        myerror(f'readSMB: SMB file  ({SMBfile}) does not exist')

    SMB = mf.getModelVarFromTiff(kwargs["SMBFile"], Q)
    SMB = icepack.interpolate(
        firedrake.max_value(firedrake.min_value(SMB, 6), -6), Q)

    return SMB



# --- Basal melt rate field --- 

def BasalMeltRate(kwargs, step, floating, Q, s, h, scenario='control'):


    # Draft depth (negative below sea level): z = s - h
    draft = firedrake.Function(Q)
    draft.interpolate(s - h)


    beginning_bmr = 20
    final_bmr = 100 #from Dutriex, 2013

    ####### TIMING OF LINEAR INCREASE IN BMR 
    if (step * kwargs['dt']) < kwargs["bmr_increase_time"]:
        melt_max = beginning_bmr
    
    # The period over which BMR increases shortens by 'x' years
    else:
        melt_max = ((final_bmr - beginning_bmr)/(kwargs["num_years"] - kwargs["bmr_increase_time"])) * ((step - (kwargs["bmr_increase_time"]/kwargs["dt"])) * kwargs["dt"]) + beginning_bmr
    
    

    ##########################################
    
    # Define thresholds from Reed et al. (2024)
    if scenario == 'control':
        
        z_min = -450  # no melt above this
        z_max = -500  # max melt below this
    
    elif scenario == 'warm':

        z_min = -400
        z_max = -450
    
    else:
        raise ValueError("Invalid scenario: must be 'control' or 'warm'")


    # Build depth-dependent piecewise melt profile
    z = draft

    # Create melt expression 
    melt_expr = firedrake.conditional(
        z >= z_min, 
        0.0,
        firedrake.conditional(
            z <= z_max, 
            melt_max, 
            melt_max * (z_min - z) / (z_min - z_max)
        )
    )

    # Apply floating mask
    melt_expr *= floating

    # Interpolate into scalar function space
    basal_melt_rate_field = firedrake.interpolate(melt_expr, Q)

    return basal_melt_rate_field, melt_max



# --- Checking the ice thickness ---

def checkThickness(kwargs):
    Q = kwargs["Q"] 
    h = kwargs["h"]
    hthresh = kwargs["hThresh"]
    h = icepack.interpolate(firedrake.max_value(hthresh, h), Q)

    return h



# --- Compute the ice volume ---

def findVolume(h):

    vol = assemble(h * dx) / 1e3**3 # in km^3

    return vol



# --- Determine area of grounded ice --- 


def findGroundedArea(bed,surface,Q):


    zF = flotationHeight(bed,Q)
    floating,grounded = flotationMask(surface,zF,Q)

    grounded_area = assemble(grounded * dx) / 1e3**2

    return grounded_area




# --- model initialization ---

def initialize_model(**kwargs):

    comm = kwargs.get('comm')


    # ---- Mesh, spaces ----
    mesh, meshOpts, Q, V = initializeMesh(**kwargs)



    # ---- Model and solver ----
    
    forward_model = icepack.models.IceStream(friction=schoofFriction,viscosity=regViscosity)
    
    opts = {"dirichlet_ids": [1],"diagnostic_solver_parameters": {"max_iterations": 150, "tolerance": 1e-6},}
    # opts = {"dirichlet_ids": meshOpts["dirichlet_ids"],"diagnostic_solver_parameters": {"max_iterations": 150, "tolerance": 1e-6},}
    
    forward_solver = icepack.solvers.FlowSolver(forward_model, **opts)



    # ---- Initial fields ----
    h, h0, s, s0, u, bed, zF, grounded, floating, A0, beta0, smb = initializeRun(kwargs, forward_solver, mesh, Q, V)



    # ---- Initial basal melt  ----
    k = 0 # time step (ICESEE uses 'k' as the time-stepping index)
    basal_melt_field, melt_max0 = BasalMeltRate(kwargs, step=k, floating=floating, Q=Q, s=s0, h=h0, scenario="control")
    


    return (h, h0, s, s0, u, bed, zF, grounded, floating, A0, beta0, smb, basal_melt_field, Q, V, forward_solver)




# --- Visualizing the model with a flowline profile DURING THE SIMULATION  ---
def initial_flowline_profile(kwargs): 
    
    # -- arrays to save the ALL of the flowline profile outputs (including the initial state) throughout the simulation -- 
    h_profiles=[]
    s_profiles = []
    
    """ define flowline in EPSG:3031 to track the profile of PIG across a seafloor ridge """
    start = np.array([-1587750., -200000.]) 
    end = np.array([-1610250., -500000.]) 

    num_points = 1000
    t_values = np.linspace(0, 1, num_points)
    profile_points = np.outer(1 - t_values, start) + np.outer(t_values, end)

    """ try extracting values from the flowline and filter out points that fail """
    valid_points = []
    bed_values, surface_values, thickness_values = [], [], []

    for p in profile_points:
        try:
            bed_values.append(kwargs["bed"].at(tuple(p)))
            surface_values.append(kwargs["s"].at(tuple(p)))
            thickness_values.append(kwargs["h"].at(tuple(p)))
            valid_points.append(tuple(p))
        except firedrake.PointNotInDomainError:
            pass  # Ignore points that are outside the mesh

    # Convert valid points back to NumPy array
    valid_points = np.array(valid_points)
    distances = np.linspace(0, np.linalg.norm(end - start), len(valid_points))

    h_profiles.append(thickness_values)
    s_profiles.append(surface_values)
    
    return(h_profiles, s_profiles, valid_points, distances, bed_values)




# --- Visualizing the model with a flowline profile DURING THE SIMULATION  ---
def flowline_profile(h, s, valid_points): 
    
    # -- arrays to save the flowline profile outputs at a particular time step -- 
    surface_values, thickness_values = [], []
    
    # -- sample values from the flowline points -- 
    for p in valid_points:
        try:
            surface_values.append(s.at(tuple(p)))
            thickness_values.append(h.at(tuple(p)))
        except firedrake.PointNotInDomainError:
            pass  # Ignore points that are outside the mesh

    
    return np.array(thickness_values), np.array(surface_values)


#########################################################################################



# --- icepack model ---
def Icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, kwargs):
    """inputs: solver - icepack solver
                h - ice thickness
                u - ice velocity
                smb - ice accumulation field
                basal_melt_field - basal melt rate 
                b - ice bed
                dt - time step
                h0 - ice thickness inflow
                kwargs - additional arguments for the model
        outputs: h - updated ice thickness
                 u - updated ice velocity
                 s - updated ice surface elevation
    """
    w2i     = float(kwargs.get('water_to_ice', 1.0))  

     # ---- net accumulation used by prognostic step ------
    a = icepack.interpolate((smb - basal_melt_field) * w2i, kwargs["Q"])


    h = solver.prognostic_solve(
        dt = dt,
        thickness = h,
        velocity = u,
        accumulation = a,
        thickness_inflow = h0,
    )
    
    #print(f"\ndt = {dt}, h0_mean = {np.mean(h0.dat.data_ro)}, h_mean = {np.mean(h.dat.data_ro)}, u_mean = {np.mean(u.dat.data_ro[:,0])}, v_mean = {np.mean(u.dat.data_ro[:,1])}, Mean_net_accumlation_field:{np.mean(a.dat.data_ro)}\n")
    #exit(1)

    kwargs["h"] = h
    h = checkThickness(kwargs)

    s = icepack.compute_surface(thickness = h, bed = bed)

    # recompute flotation height 
    zF = mf.flotationHeight(kwargs["bed"], kwargs["Q"])
    floating, grounded = mf.flotationMask(s, zF, kwargs["Q"]) # for the floating mask, 1 = floating, 0 = grounded 

    # update basal friction to reduce near the grounding line
    betaScale = mf.reduceNearGLBeta(s, kwargs["s0"], zF, grounded, kwargs["Q"], kwargs["GLThresh"])
    beta = icepack.interpolate(kwargs["beta0"] * betaScale, kwargs["Q"])

    u = solver.diagnostic_solve(
        velocity = u,
        thickness = h,
        surface = s,
        beta = beta,
        fluidity = kwargs["A0"],
        uThresh = kwargs["uThresh"],
        floating = floating,
        grounded = grounded,
    )

    return h, u, s




# --- Run model for the icepack model ---
def run_model(ensemble, **kwargs):
    
    """des: icepack model function
        inputs: ensemble - current state of the model
                **kwargs - additional arguments for the model
        outputs: model run
    """

    # unpack the **kwargs
    k       = kwargs.get('k')                               # step number in the time loop
    smb     = kwargs.get('smb', None)                       # accumulation rate field
    basal_melt_field = kwargs.get('basal_melt_field', None) # basal melt rate field
    bed     = kwargs.get('bed', None)                       # bed topography
    dt      = kwargs.get('dt', None)                        # time step size 
    nt = kwargs.get('nt', None)                             # total number of time steps
    A0      = kwargs.get('A0', None)                        # fluidity parameter
    beta0   = kwargs.get('beta0', None)                     # basal friction coefficient parameter
    Q       = kwargs.get('Q', None)                         # scalar function space
    V       = kwargs.get('V', None)                         # vector function space
    h      = kwargs.get('h', None)                          # thickness
    h0      = kwargs.get('h0', None)                        # initial thickness
    u = kwargs.get('u', None)                               # velocity 
    s      = kwargs.get('s', None)                          # elevation 
    floating = kwargs.get('floating')
    grounded = kwargs.get('grounded', None)
    solver  = kwargs.get('solver', None)                    # flow solver 
    w2i     = float(kwargs.get('water_to_ice', 1.0))        # ratio of the density of water to ice
    save_steps = kwargs.get('save_steps', None)
    ens = kwargs.get('ens_id')
    t = kwargs.get('t')


    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(**kwargs)

    # calculate the time step using ICESEE's time-stepping index (k) and the step size (dt)
    step = k 

    # -- steps at which to save profiles of the ensemble --
    flowline_profile_steps = [t/dt for t in kwargs["save_steps"]]

    h_vec = ensemble[indx_map["h"]]
    u_vec = ensemble[indx_map["u"]]
    v_vec = ensemble[indx_map["v"]]
    s_vec = ensemble[indx_map["s"]]

    #print(f"h_vec_mean={np.mean(h_vec)}, u_vec_mean = {np.mean(u_vec)}, v_vec_mean = {np.mean(v_vec)}\n")

    if kwargs["joint_estimation"]:
        basal_melt_vec = ensemble[indx_map["basal_melt_field"]]

    h = Function(Q)
    u = Function(V)
    s = Function(Q)
    basal_melt_field = Function(Q)
    
    h.dat.data[:] = h_vec
    u.dat.data[:,0] = u_vec
    u.dat.data[:,1] = v_vec
    s.dat.data[:] = s_vec
    basal_melt_field.dat.data[:] = basal_melt_vec

    print(f"\ndt = {dt}, h0_mean = {np.mean(h0.dat.data_ro)}, h_mean = {np.mean(h.dat.data_ro)}, u_mean = {np.mean(u.dat.data_ro[:,0])}, v_mean = {np.mean(u.dat.data_ro[:,1])}\n")

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

    
    #print(f"\ndt = {dt}, h0_mean = {np.mean(h0.dat.data_ro)}, h_mean = {np.mean(h.dat.data_ro)}, u_mean = {np.mean(u.dat.data_ro[:,0])}, v_mean = {np.mean(u.dat.data_ro[:,1])}\n")
     
    h, u, s = Icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, kwargs)

    # ----- joint estimation --------

    # if kwargs["joint_estimation"]:

    #     bmr_vec = ensemble[indx_map["basal_melt_field"]]
    #     basal_melt = Function(Q)
    #     basal_melt.dat.data[:] = bmr_vec.copy()
   
    # else:
    #     basal_melt = kwargs.get('basal_melt_field', None)
        
    #     if basal_melt is None:
    #         raise ValueError("basal_melt_field missing in kwargs when joint_estimation=False")

    # return a list of the updated state variables
    updated_state = {'h': h.dat.data_ro,
                     'u': u.dat.data_ro[:,0],
                     'v': u.dat.data_ro[:,1],
                     's': s.dat.data_ro}
    
    if kwargs["joint_estimation"]:
        updated_state['basal_melt_field'] = basal_melt_field.dat.data_ro

    
    # -- extract the flowline profile at particular steps -- 
    
    
    #if step in int(flowline_profile_steps):
       
    hs_ensemble_files = f"_modelrun_datasets/hs_ensemble_profiles{ens}"
    with h5py.File(hs_ensemble_files, "r") as F:
        valid_points = F["valid_points"][:]

    h_profiles, s_profiles = flowline_profile(h, s, valid_points)
    print(s_profiles.shape, h_profiles.shape,"\n")

    print(f"t[{k}] = {t[k]}")

    # -- manually saving at the time steps but need to find an efficient way to determine this --
    if t[k] == 5:
        print(f"t[{k}] = {t[k]}")
        with h5py.File(hs_ensemble_files, "a") as F:
            F["h_profiles"][:,1] = h_profiles
            F["s_profiles"][:,1] = s_profiles
    
    if t[k] == 10:
        print(f"t[{k}] = {t[k]}")
        with h5py.File(hs_ensemble_files, "a") as F:
            F["h_profiles"][:,2] = h_profiles
            F["s_profiles"][:,2] = s_profiles
    
    if t[k] == 20:
        print(f"t[{k}] = {t[k]}")
        with h5py.File(hs_ensemble_files, "a") as F:
            F["h_profiles"][:,3] = h_profiles
            F["s_profiles"][:,3] = s_profiles



    return updated_state



