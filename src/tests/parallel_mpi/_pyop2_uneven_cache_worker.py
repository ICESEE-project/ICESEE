# ==============================================================================
# @des: Real-MPI worker for test_pyop2_compile_cache.py. Reproduces the
# Idealized PIG model_nprocs=2 startup deadlock in miniature: every rank has a
# private PyOP2 cache, world rank 0 first compiles a kernel alone on
# COMM_SELF (as the mode-3 truth run does), then both ranks compile the same
# kernel on their shared 2-rank communicator. Rank 0 hits its cache, rank 1
# misses; unpatched PyOP2 deadlocks in make_so's broadcast.
#
# Usage: mpirun -n 2 python _pyop2_uneven_cache_worker.py <fix 0|1> <cache_root> <salt>
# ==============================================================================
import os
import sys

from mpi4py import MPI

fix, cache_root, salt = int(sys.argv[1]), sys.argv[2], float(sys.argv[3])
world = MPI.COMM_WORLD
private = os.path.join(cache_root, f"rank{world.rank}")
os.environ["PYOP2_CACHE_DIR"] = os.path.join(private, "pyop2")
os.environ["FIREDRAKE_TSFC_KERNEL_CACHE_DIR"] = os.path.join(private, "tsfc")
os.environ["XDG_CACHE_HOME"] = private

import firedrake as fd  # noqa: E402

if fix:
    from ICESEE.src.utils.pyop2_compile_cache import install_consistent_compile_cache

    assert install_consistent_compile_cache()


def integral(comm):
    mesh = fd.UnitSquareMesh(6, 6, comm=comm)
    x, y = fd.SpatialCoordinate(mesh)
    return fd.assemble((x * x + fd.Constant(salt) * y) * fd.dx)


if world.rank == 0:
    alone = integral(MPI.COMM_SELF)
world.Barrier()
together = integral(world)
if world.rank == 0:
    assert abs(alone - together) < 1e-12, (alone, together)
print(f"OK rank={world.rank} value={together:.12f}", flush=True)
