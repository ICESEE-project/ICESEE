# ==============================================================================
# @des: This file contains run functions for icepack data assimilation.
#       - contains different options of the EnKF data assimilation schemes.
# @date: 2024-11-4
# @author: Brian Kyanjo
# ==============================================================================
#_icepack_model.py
# --- python imports ---
import sys
import os
import copy

os.environ["OMP_NUM_THREADS"] = "1"

# firedrake imports
import firedrake
from firedrake import *
from firedrake.petsc import PETSc

# icepack imports
import icepack
import icepack.models.friction
from icepack.constants import (
    ice_density as rho_I,
    water_density as rho_W,
    gravity as g,
    weertman_sliding_law as m
)

# --- Utility imports ---
from ICESEE.config._utility_imports import icesee_get_index

# --- model initialization ---
def initialize_model(**icesee_kwargs):
    """des: initialize the icepack model"""
    # --- get the communicator from icesee_kwargs
    comm = icesee_kwargs.get('comm')

    # get size and rank of the communicator
    size = comm.Get_size()
    rank = comm.Get_rank()

    # --- Geometry and Mesh ---
    PETSc.Sys.Print('Setting up mesh across %d processes' % size)
    Lx, Ly = int(float(icesee_kwargs["Lx"])), int(float(icesee_kwargs["Ly"]))
    nx, ny = int(float(icesee_kwargs["nx"])), int(float(icesee_kwargs["ny"]))
    PETSc.Sys.Print(f"Mesh dimensions: {Lx} x {Ly} with {nx} x {ny} elements")

    # --- make the comm object available to the mesh function
    mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly, quadrilateral=True, comm=comm)

    # -- get the degree of the finite element space
    degree = int(float(icesee_kwargs["degree"]))
    Q = firedrake.FunctionSpace(mesh, "CG",degree)
    V = firedrake.VectorFunctionSpace(mesh, "CG", degree)
    x,y = firedrake.SpatialCoordinate(mesh)

    # --- Bedrock and Surface Elevations ---
    b_in, b_out = (float(icesee_kwargs["b_in"])), (float(icesee_kwargs["b_out"]))
    s_in, s_out = (float(icesee_kwargs["s_in"])), (float(icesee_kwargs["s_out"]))

    b = firedrake.Function(Q).interpolate(b_in - (b_in - b_out) * x / Lx)
    s0 = firedrake.Function(Q).interpolate(s_in - (s_in - s_out) * x / Lx)
    h0 = firedrake.Function(Q).interpolate(s0 - b)

    # --- Driving Stress ---
    h_in = s_in - b_in
    ds_dx = (s_out - s_in) / Lx
    tau_D = -rho_I * g * h_in * ds_dx
    PETSc.Sys.Print(f"Driving stress = {1000*tau_D} kPa")

    # --- Initial Velocity ---
    u_in, u_out = float(icesee_kwargs["u_in"]), float(icesee_kwargs["u_out"])
    velocity_x = u_in + (u_out - u_in) * (x / Lx) ** 2
    u0 = firedrake.Function(V).interpolate(firedrake.as_vector((velocity_x, 0)))

    # --- Friction Coefficient ---
    PETSc.Sys.Print("Importing icepack ...")
    T = firedrake.Constant(float(icesee_kwargs["T"]))
    A = icepack.rate_factor(T)

    expr = (0.95 - 0.05 * x / Lx) * tau_D / u_in**(1 / m)
    C = firedrake.Function(Q).interpolate(expr)

    p_W = rho_W * g * firedrake.max_value(0, h0 - s0)
    p_I = rho_I * g * h0
    phi = 1 - p_W / p_I

    # --- Friction Law ---
    def weertman_friction_with_ramp(**icesee_kwargs):
        u = icesee_kwargs["velocity"]
        h = icesee_kwargs["thickness"]
        s = icesee_kwargs["surface"]
        C = icesee_kwargs["friction"]

        p_W = rho_W * g * firedrake.max_value(0, h - s)
        p_I = rho_I * g * h
        phi = 1 - p_W / p_I
        return icepack.models.friction.bed_friction(
            velocity=u,
            friction=C * phi,
        )

    # --- Ice Stream Model ---
    model_weertman = icepack.models.IceStream(friction=weertman_friction_with_ramp)

    opts = {"dirichlet_ids": [1], "side_wall_ids": [3, 4]}
    solver_weertman = icepack.solvers.FlowSolver(model_weertman, **opts)

    u0 = solver_weertman.diagnostic_solve(
        velocity=u0,
        thickness=h0,
        surface=s0,
        fluidity=A,
        friction=C,
    )

    expr = -1e3 * C * phi * sqrt(inner(u0, u0)) ** (1 / m - 1) * u0
    tau_b = firedrake.Function(V).interpolate(expr)

    # --- Accumulation ---
    a_in = firedrake.Constant(float(icesee_kwargs["a_in"]))
    da   = firedrake.Constant(float(icesee_kwargs["da"]))
    a    = firedrake.Function(Q).interpolate(a_in + da * x / Lx)

    # nurged accumulation
    a_in_p  = firedrake.Constant(float(icesee_kwargs["a_in_p"]))
    da_p    = firedrake.Constant(float(icesee_kwargs["da_p"]))
    a_p     = firedrake.Function(Q).interpolate(a_in_p + da_p * x / Lx)

    # --- Update h and u ---
    h = h0.copy(deepcopy=True)
    u = u0.copy(deepcopy=True)

    # print size h
    # print(f"Size of the function space: {h.dat.data.size} on rank {rank}")
    icesee_kwargs.update({"nx":nx, "ny":ny, "Lx":Lx, "Ly":Ly, "x":x, "y":y, "h":h, "u":u, "a":a, "a_p":a_p, "b":b, "b_in":b_in, "b_out":b_out, "h0":h0, "u0":u0, "solver_weertman":solver_weertman, "A":A, "C":C, "Q":Q, "V":V, "mesh":mesh})
    return icesee_kwargs

# --- icepack model ---
def Icepack(solver, h, u, a, b, dt, h0, **icesee_kwargs):
    """inputs: solver - icepack solver
                h - ice thickness
                u - ice velocity
                a - ice accumulation
                b - ice bed
                dt - time step
                h0 - ice thickness inflow
                *args - additional arguments for the model
        outputs: h - updated ice thickness
                 u - updated ice velocity
    """
    h = solver.prognostic_solve(
        dt = dt,
        thickness = h,
        velocity = u,
        accumulation = a,
        thickness_inflow = h0,
    )

    s = icepack.compute_surface(thickness = h, bed = b)

    u = solver.diagnostic_solve(
        velocity = u,
        thickness = h,
        surface = s,
        **icesee_kwargs
    )

    return h, u

# --- Run model for the icepack model ---
def run_model(ensemble, **icesee_kwargs):
    """des: icepack model function
        inputs: ensemble - current state of the model
                **icesee_kwargs - additional arguments for the model
        outputs: model run
    """

    # unpack the **icesee_kwargs
    # a = icesee_kwargs.get('a', None)
    b  = icesee_kwargs.get('b', None)
    dt = icesee_kwargs.get('dt', None)
    h0 = icesee_kwargs.get('h0', None)
    A  = icesee_kwargs.get('A', None)
    C  = icesee_kwargs.get('C', None)
    Q  = icesee_kwargs.get('Q', None)
    V  = icesee_kwargs.get('V', None)
    solver = icesee_kwargs.get('solver', None)

    # --- define the state variables list ---
    global vec_inputs

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(ensemble,**icesee_kwargs)

    # joint estimation
    if icesee_kwargs["joint_estimation"]:
        # Use analysis step to update the accumulation rate
        # - pack accumulation rate with the state variables to
        #   get ensemble = [h,u,v,a]
        # a_vec = ensemble[indx_map["smb"],ens]
        a_vec = ensemble[indx_map["smb"]]

        a = Function(Q)
        # print(f"a size: {a.dat.data.size} avec {a_vec.shape} ensemble shape: {ensemble.shape}")
        a.dat.data[:] = a_vec.copy()
    else:
        # don't update the accumulation rate (updates smb)
        a = icesee_kwargs.get('a', None)

    # create firedrake functions from the ensemble members
    h = Function(Q)
    # h.dat.data[:] = ensemble[indx_map["h"],ens]
    # h.dat.data[:] = ensemble[indx_map["h"]]
    h.dat.data[:] = copy.deepcopy(ensemble[indx_map["h"]])

     # create firedrake functions from the ensemble members
     # u = Function(V)
     #

    u = Function(V)
    # u.dat.data[:,0] = ensemble[indx_map["u"],ens]
    # u.dat.data[:,1] = ensemble[indx_map["v"],ens]
    u.dat.data[:,0] = ensemble[indx_map["u"]]
    u.dat.data[:,1] = ensemble[indx_map["v"]]

    # call the ice stream model to update the state variables
    h, u = Icepack(solver, h, u, a, b, dt, h0, fluidity = A, friction = C)

    # return a list of the updated state variables
    updated_state = {'h': copy.deepcopy(h.dat.data_ro),
                     'u': u.dat.data_ro[:,0],
                     'v': u.dat.data_ro[:,1]}

    if icesee_kwargs["joint_estimation"]:
        updated_state['smb'] = a_vec
    # else:
    #     updated_state['smb'] = a.dat.data_ro

    return updated_state

from ICESEE.applications.icepack_model.icepack_utils._coordinates import (
    register_icepack_coordinate_provider,
)

register_icepack_coordinate_provider()
