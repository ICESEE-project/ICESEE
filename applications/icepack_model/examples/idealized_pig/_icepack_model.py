
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

from ICESEE.src.utils.init_trace import traced_init

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

def initializeMesh(**icesee_kwargs):
    initFile = icesee_kwargs["initFile"]
    meshFile = icesee_kwargs["meshFile"]
    comm = icesee_kwargs.get("comm")
    if comm.Get_size() == 1:
        comm = firedrake.COMM_SELF
        icesee_kwargs["comm"] = comm
    with traced_init("checkpoint mesh load", comm=comm):
        meshI = mf.getMeshFromCheckPoint(initFile,icesee_kwargs)
    with traced_init("mesh construction (argus + function spaces)", comm=comm):
        mesh, Q, V, meshOpts = \
            mf.setupMesh(meshFile, degree=1,
                         meshOversample=2,
                         newMesh=meshI,
                         comm=comm)

    return mesh, meshOpts, Q, V




# ---- initializing the run ----

def initializeRun(icesee_kwargs, forward_solver, mesh, Q, V):

    # Checkpoint collectives must use the communicator that owns this mesh.
    checkpoint_comm = icesee_kwargs.get("comm") or firedrake.COMM_WORLD
    if checkpoint_comm.Get_size() == 1:
        checkpoint_comm = firedrake.COMM_SELF
    # Compact inputs store fields without the spin-up history's timestep index.
    compact_initialization = bool(icesee_kwargs.get("compact_initialization", False))
    with traced_init("checkpoint field load", comm=checkpoint_comm, compact=compact_initialization), \
            firedrake.CheckpointFile(icesee_kwargs["initFile"], 'r', comm=checkpoint_comm) as checkpoint:
        if compact_initialization:
            velocity = checkpoint.load_function(mesh, "velocity")
            h0 = checkpoint.load_function(mesh, "thickness")
            s0 = checkpoint.load_function(mesh, "surface")
        else:
            velocity = checkpoint.load_function(mesh, "velocity", idx=20000) # idx index set to the LAST time step of the steady-state run (20,000 for 1000 year run)
            h0 = checkpoint.load_function(mesh, "thickness", idx=20000)
            s0 = checkpoint.load_function(mesh, "surface", idx=20000)
        bed = checkpoint.load_function(mesh, "bed")
        grounded0 = checkpoint.load_function(mesh, "grounded")
        floating0 = checkpoint.load_function(mesh, "floating")
        A0 = checkpoint.load_function(mesh, "fluidity")
        beta0 = checkpoint.load_function(mesh, "extended_beta")

    """ find initial velocity """
    uThresh = firedrake.Constant(icesee_kwargs['uThresh'])
    with traced_init("initial diagnostic solve", comm=checkpoint_comm):
        u0 = forward_solver.diagnostic_solve(velocity=velocity, thickness=h0,
                                                    surface=s0,
                                                    beta=beta0, fluidity=A0,
                                                    grounded=grounded0,
                                                    floating=floating0,
                                                    uThresh=uThresh)

    """ define initial state """
    h, hLast, s, u, zF, grounded, floating = initialState(h0, s0, u0, bed, grounded0, floating0, Q)

    """ define smb and melt """
    smb = readSMB(icesee_kwargs,Q)

    return h, h0, s, s0, u, bed, zF, grounded, floating, A0, beta0, smb




# ---- basal friction model ---
def schoofFriction(velocity, grounded, beta, uThresh):


    C = grounded * beta**2
    mExp = (1./m + 1.)
    U = firedrake.sqrt(firedrake.inner(velocity, velocity))

    return C * ((uThresh**mExp + U**mExp)**(m/(m+1.)) - uThresh)



# --- rheology model ---

def regViscosity(**icesee_kwargs):

    u = icesee_kwargs["velocity"]
    h = icesee_kwargs["thickness"]
    A = icesee_kwargs["fluidity"]

    return icepack.models.viscosity.viscosity_depth_averaged(velocity=u, thickness=h, fluidity= A)




# --- read in the SMB file ---

def readSMB(icesee_kwargs, Q):


    if not os.path.exists(icesee_kwargs["SMBFile"]):
        myerror(f'readSMB: SMB file  ({SMBfile}) does not exist')

    SMB = mf.getModelVarFromTiff(icesee_kwargs["SMBFile"], Q)
    SMB = icepack.interpolate(
        firedrake.max_value(firedrake.min_value(SMB, 6), -6), Q)

    return SMB



# --- Basal melt rate field ---

# Canonical forcing selectors: truth ramps melt_max from 20 to 100;
# ensemble forecasts retain 20. Validate string selectors to avoid ambiguous
# boolean/string comparisons.
EXPERIMENT_TRUE = "true"    # reference/true trajectory: melt_max ramps 20 -> 100 over the run
EXPERIMENT_WRONG = "wrong"  # ensemble/forecast trajectory: melt_max held constant at 20

# melt_max changes every step of the true trajectory. Carrying it in one
# mutable Constant keeps the melt expression identical from step to step, so
# Firedrake compiles its interpolation kernel once instead of once per step
# (a changing Python float is a new literal in the form each time).
_MELT_MAX = None


def _melt_max_constant(value):
    global _MELT_MAX
    if _MELT_MAX is None:
        _MELT_MAX = firedrake.Constant(value)
    else:
        _MELT_MAX.assign(value)
    return _MELT_MAX

def BasalMeltRate(icesee_kwargs, step, floating, Q, s, h, scenario='control', experiment=EXPERIMENT_TRUE):

    if experiment not in (EXPERIMENT_TRUE, EXPERIMENT_WRONG):
        raise ValueError(
            f"BasalMeltRate: experiment must be {EXPERIMENT_TRUE!r} or "
            f"{EXPERIMENT_WRONG!r}, got {experiment!r}"
        )

    # Draft depth (negative below sea level): z = s - h
    draft = firedrake.Function(Q)
    draft.interpolate(s - h)


    beginning_bmr = 20
    final_bmr = 100 #from Dutriex, 2013

    ####### TIMING OF LINEAR INCREASE IN BMR (true trajectory only) #######
    if experiment == EXPERIMENT_WRONG:
        melt_max = beginning_bmr
    elif (step * icesee_kwargs['dt']) < icesee_kwargs["bmr_increase_time"]:
        melt_max = beginning_bmr

    # The period over which BMR increases shortens by 'x' years
    else:
        melt_max = ((final_bmr - beginning_bmr)/(icesee_kwargs["num_years"] - icesee_kwargs["bmr_increase_time"])) * ((step - (icesee_kwargs["bmr_increase_time"]/icesee_kwargs["dt"])) * icesee_kwargs["dt"]) + beginning_bmr



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
    melt_max_c = _melt_max_constant(melt_max)

    # Create melt expression
    melt_expr = firedrake.conditional(
        z >= z_min,
        0.0,
        firedrake.conditional(
            z <= z_max,
            melt_max_c,
            melt_max_c * (z_min - z) / (z_min - z_max)
        )
    )

    # Apply floating mask
    melt_expr *= floating

    # Interpolate into scalar function space.
    #
    # firedrake.interpolate(expr, Q) (the free-function form) now returns a
    # lazy symbolic Interpolate object in the installed Firedrake/UFL
    # version, not an evaluated Function -- Function(Q).interpolate(expr)
    # is the eager-evaluation equivalent (same math; this call runs on
    # every forecast step via BasalMeltRate, so it must stay eager).
    basal_melt_rate_field = firedrake.Function(Q).interpolate(melt_expr)

    return basal_melt_rate_field, melt_max



def profile_path(data_path, name):
    """Path of the flowline-profile file ``<data_path>/<name>.h5``. Runs made
    before the .h5 extension was added wrote ``<data_path>/<name>``; that
    file is returned instead when it is the only one present."""
    path = os.path.join(data_path, f"{name}.h5")
    legacy = os.path.join(data_path, name)
    if not os.path.exists(path) and os.path.exists(legacy):
        return legacy
    return path


# --- Checking the ice thickness ---

def checkThickness(icesee_kwargs):
    Q = icesee_kwargs["Q"]
    h = icesee_kwargs["h"]
    hthresh = icesee_kwargs["hThresh"]
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

def initialize_model(**icesee_kwargs):

    comm = icesee_kwargs.get('comm')
    if comm.Get_size() == 1:
        comm = firedrake.COMM_SELF
        icesee_kwargs["comm"] = comm


    # ---- Mesh, spaces ----
    mesh, meshOpts, Q, V = initializeMesh(**icesee_kwargs)



    # ---- Model and solver ----

    forward_model = icepack.models.IceStream(friction=schoofFriction,viscosity=regViscosity)

    opts = {"dirichlet_ids": [1],"diagnostic_solver_parameters": {"max_iterations": 150, "tolerance": 1e-6},}
    # opts = {"dirichlet_ids": meshOpts["dirichlet_ids"],"diagnostic_solver_parameters": {"max_iterations": 150, "tolerance": 1e-6},}

    with traced_init("Icepack solver construction", comm=comm):
        forward_solver = icepack.solvers.FlowSolver(forward_model, **opts)



    # ---- Initial fields ----
    h, h0, s, s0, u, bed, zF, grounded, floating, A0, beta0, smb = initializeRun(icesee_kwargs, forward_solver, mesh, Q, V)



    # ---- Initial basal melt  ----
    k = 0 # time step (ICESEE uses 'k' as the time-stepping index)
    basal_melt_field, melt_max0 = BasalMeltRate(icesee_kwargs, step=k, floating=floating, Q=Q, s=s0, h=h0, scenario="control", experiment=EXPERIMENT_TRUE)



    return (h, h0, s, s0, u, bed, zF, grounded, floating, A0, beta0, smb, basal_melt_field, Q, V, forward_solver)




# --- Visualizing the model with a flowline profile DURING THE SIMULATION  ---
def initial_flowline_profile(icesee_kwargs):
    from ._mode3_profiles import sample_fields
    start = np.array([-1587750., -200000.])
    end = np.array([-1610250., -500000.])
    points = start + np.linspace(0., 1., 1000)[:, None] * (end - start)
    values = sample_fields((icesee_kwargs["h"], icesee_kwargs["s"], icesee_kwargs["bed"]), points)
    valid = np.array([all(v[i] is not None for v in values) for i in range(len(points))])
    valid_points = points[valid]
    h, s, bed = [np.asarray([v[i] for i in np.flatnonzero(valid)], dtype=float).reshape(-1) for v in values]
    distances = np.linalg.norm(valid_points - start, axis=1)
    return [h], [s], valid_points, distances, bed


def flowline_profile(h, s, valid_points):
    from ._mode3_profiles import sample_fields
    values = sample_fields((h, s), valid_points)
    if any(value is None for field_values in values for value in field_values):
        raise ValueError("A saved flowline point is outside the current model mesh")
    return tuple(np.asarray(value, dtype=float).reshape(-1) for value in values)


#########################################################################################



# --- icepack model ---
def Icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, icesee_kwargs):
    """inputs: solver - icepack solver
                h - ice thickness
                u - ice velocity
                smb - ice accumulation field
                basal_melt_field - basal melt rate
                b - ice bed
                dt - time step
                h0 - ice thickness inflow
                icesee_kwargs - additional arguments for the model
        outputs: h - updated ice thickness
                 u - updated ice velocity
                 s - updated ice surface elevation
                 floating - recomputed floating mask (diagnostic, not DA
                   state -- a pure function of (s, bed); callers should
                   recompute it fresh each call rather than cache it)
                 grounded - recomputed grounded mask (diagnostic, complement of floating)
    """
    w2i     = float(icesee_kwargs.get('water_to_ice', 1.0))

     # ---- net accumulation used by prognostic step ------
    a = icepack.interpolate((smb - basal_melt_field) * w2i, icesee_kwargs["Q"])


    h = solver.prognostic_solve(
        dt = dt,
        thickness = h,
        velocity = u,
        accumulation = a,
        thickness_inflow = h0,
    )

    #exit(1)

    icesee_kwargs["h"] = h
    h = checkThickness(icesee_kwargs)

    s = icepack.compute_surface(thickness = h, bed = bed)

    # recompute flotation height
    zF = mf.flotationHeight(icesee_kwargs["bed"], icesee_kwargs["Q"])
    floating, grounded = mf.flotationMask(s, zF, icesee_kwargs["Q"]) # for the floating mask, 1 = floating, 0 = grounded

    # update basal friction to reduce near the grounding line
    betaScale = mf.reduceNearGLBeta(s, icesee_kwargs["s0"], zF, grounded, icesee_kwargs["Q"], icesee_kwargs["GLThresh"])
    beta = icepack.interpolate(icesee_kwargs["beta0"] * betaScale, icesee_kwargs["Q"])

    u = solver.diagnostic_solve(
        velocity = u,
        thickness = h,
        surface = s,
        beta = beta,
        fluidity = icesee_kwargs["A0"],
        uThresh = icesee_kwargs["uThresh"],
        floating = floating,
        grounded = grounded,
    )

    return h, u, s, floating, grounded




# --- Run model for the icepack model ---
def run_model(ensemble, **icesee_kwargs):

    """des: icepack model function
        inputs: ensemble - current state of the model
                **icesee_kwargs - additional arguments for the model
        outputs: model run
    """

    # unpack the **icesee_kwargs
    k       = icesee_kwargs.get('k')                               # step number in the time loop
    smb     = icesee_kwargs.get('smb', None)                       # accumulation rate field
    basal_melt_field = icesee_kwargs.get('basal_melt_field', None) # basal melt rate field
    bed     = icesee_kwargs.get('bed', None)                       # bed topography
    dt      = icesee_kwargs.get('dt', None)                        # time step size
    nt = icesee_kwargs.get('nt', None)                             # total number of time steps
    A0      = icesee_kwargs.get('A0', None)                        # fluidity parameter
    beta0   = icesee_kwargs.get('beta0', None)                     # basal friction coefficient parameter
    Q       = icesee_kwargs.get('Q', None)                         # scalar function space
    V       = icesee_kwargs.get('V', None)                         # vector function space
    h      = icesee_kwargs.get('h', None)                          # thickness
    h0      = icesee_kwargs.get('h0', None)                        # initial thickness
    u = icesee_kwargs.get('u', None)                               # velocity
    s      = icesee_kwargs.get('s', None)                          # elevation
    floating = icesee_kwargs.get('floating')
    grounded = icesee_kwargs.get('grounded', None)
    solver  = icesee_kwargs.get('solver', None)                    # flow solver
    w2i     = float(icesee_kwargs.get('water_to_ice', 1.0))        # ratio of the density of water to ice
    save_steps = icesee_kwargs.get('save_steps', None)
    ens = icesee_kwargs.get('ens_id')
    t = icesee_kwargs.get('t')


    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(**icesee_kwargs)

    # calculate the time step using ICESEE's time-stepping index (k) and the step size (dt)
    step = k

    print(step)

    # -- steps at which to save profiles of the ensemble --
    flowline_profile_steps = [t/dt for t in icesee_kwargs["save_steps"]]

    h_vec = ensemble[indx_map["h"]]
    u_vec = ensemble[indx_map["u"]]
    v_vec = ensemble[indx_map["v"]]
    s_vec = ensemble[indx_map["s"]]


    # The configured state vector includes basal_melt_field independently of joint estimation.
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



    h, u, s, floating, grounded = Icepack(solver, h, u, smb, basal_melt_field, bed, dt, h0, icesee_kwargs)

    # ----- joint estimation --------

    # if icesee_kwargs["joint_estimation"]:

    #     bmr_vec = ensemble[indx_map["basal_melt_field"]]
    #     basal_melt = Function(Q)
    #     basal_melt.dat.data[:] = bmr_vec.copy()

    # else:
    #     basal_melt = icesee_kwargs.get('basal_melt_field', None)

    #     if basal_melt is None:
    #         raise ValueError("basal_melt_field missing in icesee_kwargs when joint_estimation=False")

    # return a list of the updated state variables. basal_melt_field is
    # always part of the packed vector (see the unconditional read above).
    updated_state = {'h': h.dat.data_ro,
                     'u': u.dat.data_ro[:,0],
                     'v': u.dat.data_ro[:,1],
                     's': s.dat.data_ro,
                     'basal_melt_field': basal_melt_field.dat.data_ro}


    # Sample only configured physical-year save points. Keep the established
    # per-member .h5 schema and respect the run's data_path.
    columns = [
        index + 1 for index, year in enumerate(save_steps)
        if np.isclose(float(step + 1) * float(dt), float(year), rtol=0.0, atol=1e-9)
    ]
    if columns:
        _modelrun_datasets = icesee_kwargs.get("data_path") or "_modelrun_datasets"
        hs_ensemble_files = profile_path(_modelrun_datasets, f"hs_ensemble_profiles{ens}")
        with h5py.File(hs_ensemble_files, "r") as F:
            valid_points = F["valid_points"][:]
        h_profiles, s_profiles = flowline_profile(h, s, valid_points)
        with h5py.File(hs_ensemble_files, "a") as F:
            for column in columns:
                F["h_profiles"][:, column] = h_profiles
                F["s_profiles"][:, column] = s_profiles
                if "written_steps" in F:
                    F["written_steps"][column] = step + 1

    return updated_state



