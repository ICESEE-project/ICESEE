import firedrake
from firedrake import *
import icepack
import icepack.models.friction
import tqdm
import numpy as np
from icepack.constants import weertman_sliding_law as m
import os
import h5py
from firedrake.petsc import PETSc

# import mpi4py
# Deactivate MPI auto-initialization and finalization 
# mpi4py.rc.initialize = False
# mpi4py.rc.finalize = False
# init_MPI = False
from mpi4py import MPI

Lx, Ly = 50e3, 12e3
nx, ny = 32,48

num_years = 3
timesteps_per_year = 2
b_in, b_out = 20, -40
s_in, s_out = 850, 50

# def forecast_step(t, ens, ensemble_local, a, b, h0, u0, solver_weertman, δt,A,C,Q,V):
def forecast_step(t, ens, ensemble_local, a, b, h0, u0, solver_weertman, δt,A,C,Q,V):

    # hdim = h.dat.data.size
    # hdim_loc = ensemble_local.shape[0]
    
    # extract h 
    h = firedrake.Function(Q)
    u = firedrake.Function(V)
    hdim_loc = h.dat.data.size
    # print (f"[DEBUG_in] Rank {rank} - hdim: {hdim_loc} loc ens size: {ensemble_local.shape} ens: {ens}")
   
    if t==0:
        hvec = h0.dat.data_ro.copy()
        uvec = u0.dat.data_ro[:,0].copy()
        vvec = u0.dat.data_ro[:,1].copy()
    else:
        hvec = ensemble_local[:hdim_loc,ens]
        uvec = ensemble_local[hdim_loc:2*hdim_loc,ens]
        vvec = ensemble_local[2*hdim_loc:3*hdim_loc,ens]

    h.dat.data[:] = hvec.copy()
    u.dat.data[:,0] = uvec.copy()
    u.dat.data[:,1] = vvec.copy()

    # print (f"[DEBUG] Rank {rank} - hdim_: {h_.dat.data.size}")
    
    h = solver_weertman.prognostic_solve(
        δt,
        thickness=h,
        velocity=u,
        accumulation=a,
        thickness_inflow=h0,
    )
    s = icepack.compute_surface(thickness=h, bed=b)

    u = solver_weertman.diagnostic_solve(
        velocity=u,
        thickness=h,
        surface=s,
        fluidity=A,
        friction=C,
    )

    # update the ensemble
    ensemble_local[:hdim_loc,ens] = h.dat.data_ro
    ensemble_local[hdim_loc:2*hdim_loc,ens] = u.dat.data_ro[:,0]
    ensemble_local[2*hdim_loc:3*hdim_loc,ens] = u.dat.data_ro[:,1]

    # add some noise to the ensemble
    # ensemble[:hdim,ens] += 1e-4 * np.random.randn(hdim)

    return  ensemble_local[:3*hdim_loc,ens]


def intialize(nx,ny,Lx,Ly,b_in,b_out,s_in,s_out,rank,comm,size):
    # --- to only be ran on rank 0 and made available to all ranks ---
    # os.environ["MPIEXEC"] = " " 
    # mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly)

    # comm = COMM_WORLD.Split(COMM_WORLD.rank % 2)
    # comm = MPI.COMM_WORLD
    # comm = comm.Split(rank % 2)

    # -----working code for mesh -----
    comm = comm.Split(rank % size)
    mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly, quadrilateral=True, comm=comm)
    # --------------------------------

    # PETSc.Sys.Print('setting up mesh across %d processes' % comm.Get_size())
    # if COMM_WORLD.rank % 2 == 0:
    # if rank % 2 == 0:
    #     # Even ranks create a quad mesh
    #     mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly, quadrilateral=True, comm=comm)
    # else:
    #     # Odd ranks create a triangular mesh
    #     mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly, comm=comm)

    Q = firedrake.FunctionSpace(mesh, "CG", 2)
    V = firedrake.VectorFunctionSpace(mesh, "CG", 2)
    x, y = firedrake.SpatialCoordinate(mesh)

    b = firedrake.interpolate(b_in - (b_in - b_out) * x / Lx, Q)

    s0 = firedrake.interpolate(s_in - (s_in - s_out) * x / Lx, Q)

    h0 = firedrake.interpolate(s0 - b, Q)

    from icepack.constants import (
        ice_density as ρ_I,
        water_density as ρ_W,
        gravity as g,
    )
    # comm.Barrier()
    h_in = s_in - b_in
    δs_δx = (s_out - s_in) / Lx
    τ_D = -ρ_I * g * h_in * δs_δx
    print(f"{1000 * τ_D} kPa")

    u_in, u_out = 20, 2400
    velocity_x = u_in + (u_out - u_in) * (x / Lx) ** 2
    u0 = firedrake.interpolate(firedrake.as_vector((velocity_x, 0)), V)

    T = firedrake.Constant(255.0)
    A = icepack.rate_factor(T)

    expr = (0.95 - 0.05 * x / Lx) * τ_D / u_in**(1 / m)
    C = firedrake.interpolate(expr, Q)

    p_W = ρ_W * g * firedrake.max_value(0, h0 - s0)
    p_I = ρ_I * g * h0
    ϕ = 1 - p_W / p_I


    def weertman_friction_with_ramp(**kwargs):
        u = kwargs["velocity"]
        h = kwargs["thickness"]
        s = kwargs["surface"]
        C = kwargs["friction"]

        p_W = ρ_W * g * firedrake.max_value(0, h - s)
        p_I = ρ_I * g * h
        ϕ = 1 - p_W / p_I
        return icepack.models.friction.bed_friction(
            velocity=u,
            friction=C * ϕ,
        )

    model_weertman = icepack.models.IceStream(friction=weertman_friction_with_ramp)
    opts = {"dirichlet_ids": [1], "side_wall_ids": [3, 4]}
    solver_weertman = icepack.solvers.FlowSolver(model_weertman, **opts)

    u0 = solver_weertman.diagnostic_solve(
        velocity=u0, thickness=h0, surface=s0, fluidity=A, friction=C
    )

    δt = 1.0 / timesteps_per_year
    num_timesteps = num_years * timesteps_per_year

    a_in = firedrake.Constant(1.7)
    δa = firedrake.Constant(-2.7)
    a = firedrake.interpolate(a_in + δa * x / Lx, Q)

    h = h0.copy(deepcopy=True)
    u = u0.copy(deepcopy=True)

    return h, u, a, b, h0,u0,solver_weertman, δt,A,C,Q,V,num_timesteps


def synthetic_ice_stream(rank,size,comm,Nens):

    # --- call the initialization function on rank 0 and write the function spaces to a file ---
    h, u, a, b, h0, u0, solver_weertman, δt,A,C,Q,V,num_timesteps = intialize(nx,ny,Lx,Ly,b_in,b_out,s_in,s_out,rank,comm,size)

    # --- Properly Distribute Tasks for All Cases ---
    # Case 1: More ensembles than processes → Distribute as evenly as possible
    mem_per_task = Nens // size  # Base number of tasks per process
    remainder = Nens % size       # Extra tasks to distribute

    if rank < remainder:
        # the first remainder gets mem_per_task+1 tasks each
        start = rank * (mem_per_task + 1)
        stop = start + (mem_per_task + 1)
        # stop = start + mem_per_task
    else:
        #  the remaining (size - remainder) get mem_per_task tasks each
        start = rank * mem_per_task + remainder
        stop = start + mem_per_task
        # stop = start + mem_per_task-1

    # --- Ensure All Ranks Participate (No Deadlocks) ---
    if start == stop:
        print(f"[Rank {rank}] No work assigned. Waiting at barrier.")
    else:
        print(f"[Rank {rank}] Processing ensembles {start} to {stop}")
        

    sizes_rank = []
    for size_loc in range (start,stop):
        hdim_loc = h.dat.data.size
        sizes_rank.append(hdim_loc)

    #  form local ensemble based on the coresponding rank
    # print (f"[DEBUG] Rank {rank} - sizes_rank: {sizes_rank}")
    hdim_rank = sizes_rank[-1]
    ensemble_local = np.zeros((3*hdim_rank, stop-start))
    
    # -- loop over the timesteps --
    for step in tqdm.trange(num_timesteps):

        # # if rank == 0:
        # print(f"[Rank {rank}] Processing ensembles {start} to {stop}, expected shape: {ensemble_local.shape}")
        # t = 0
        # global_dic = {}
        for ens in range(ensemble_local.shape[1]):
            # unique tag for each ens member corresponding to the rank
            # tag=int((rank * np.ceil(Nens/size)) + ens)
            ensemble_local[:3*hdim_rank,ens] = forecast_step(step , ens, ensemble_local, a, b, h0, u0,solver_weertman, δt,A,C,Q,V)
            # print (f"[DEBUG] Rank {rank} Processing ensembles {tag} of {Nens-1}, expected shape: {ensemble_local[:3*hdim_rank,ens].shape}")
            # add the ensemble to the global dictionary
            # global_dic[tag] = ensemble_local[:3*hdim_rank,ens]

        # Now we need to gather the ensemble
        gathered_ensemble = comm.allgather(ensemble_local)
        # Gather all local array shapes
        # all_shapes = comm.gather(ensemble_local.shape, root=0)
        # ============================================================
        # comm.Barrier() 
        if rank == 0:
            print(f"[Rank {rank}] Gathered shapes: {[arr.shape for arr in gathered_ensemble]}")
            # hstack the gathered ensemble
            # ensemble = np.hstack(gathered_ensemble)
            # print(f"[Rank {rank}] ensemble shape after gather: {ensemble.shape}")
            # print last 10 ensembles
            # print("Ensemble Matrix Data (Sample):", ensemble[-10:,:])
   

# main execution
if __name__ == "__main__":
    import time
    # import mpi4py.MPI as MPI
    from firedrake import *
    # import mpi4py
    # Deactivate MPI auto-initialization and finalization 
    # mpi4py.rc.initialize = False
    import warnings
    warnings.filterwarnings("ignore")

    Nens = 64         # Number of ensemble members
    n_modeltasks = 1 # Number of model tasks

    start_time = time.time()


    # Initialize parallel manager
    from icesee_mpi_parallel_manager import icesee_mpi_parallelization
    parallel_manager = icesee_mpi_parallelization(Nens, n_modeltasks, screen_output=0)

    # Get MPI information
    comm = parallel_manager.COMM_model
    rank = parallel_manager.rank_model
    size = parallel_manager.size_model
    task_id = parallel_manager.task_id
    # #  ----
    print(f"[Rank {rank}] Task ID: {task_id}")

    # --- Execute Assigned Tasks ---
    # for ens in range(start, stop):
    ensemble = synthetic_ice_stream(rank, size, comm,Nens)

    # --- Synchronization Before Finalizing ---
    comm.Barrier()  # Ensures all ranks reach this point before proceeding
    if rank == 0:
        print(ensemble)

    # --- Report Execution Time ---
    if rank == 0:
        print(f"\n[ICESEE] Time taken: {(time.time() - start_time) / 60:.2f} minutes")