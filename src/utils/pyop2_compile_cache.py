# ==============================================================================
# @des: Communicator-consistent PyOP2 compiled-kernel (.so) disk cache.
#
#   PyOP2 2025.x ``pyop2.compilation.make_so`` is a ``parallel_cache(bcast=False)``
#   around a body that is collective over the node-local compilation
#   communicator (rank 0 compiles, then ``ccomm.bcast(soname)``). Each rank
#   looks the kernel up in the disk cache on its own, so when the ranks of one
#   model communicator see different cache contents, the ranks that miss
#   enter the broadcast while the ranks that hit do not, and the model group
#   deadlocks. Under ICESEE this happens whenever a model group has more than
#   one rank:
#
#   * rank-private caches -- the world root has already compiled the model's
#     kernels alone (truth/reference setup on COMM_SELF), so it hits where its
#     group partners miss: a deterministic deadlock (observed on PACE with
#     2 ranks, Nens=1, model_nprocs=2);
#   * caches shared by several model groups (node-local /tmp or the default
#     $VIRTUAL_ENV/.cache on NFS) -- another group can write the kernel between
#     two partners' lookups: an intermittent deadlock.
#
#   ``install_consistent_compile_cache()`` replaces ``make_so`` with a wrapper
#   that agrees on hit/miss across the compilation communicator (one
#   ``allreduce``) before anyone enters the collective body, and then uses the
#   same on-disk cache (same directory, keys and atomic write) as PyOP2. The
#   compiled code and the kernels loaded are unchanged; only the deadlock is
#   removed. It also records per-rank compile counts/time for the ICESEE
#   performance report ("jit_compile").
#
#   Set ICESEE_PYOP2_CONSISTENT_CACHE=0 to leave PyOP2 untouched. If the
#   installed PyOP2 does not have the expected structure the wrapper is not
#   installed and a warning is printed once.
# ==============================================================================

from __future__ import annotations

import os
import sys
import time

_INSTALLED = False


def _disabled() -> bool:
    return os.environ.get("ICESEE_PYOP2_CONSISTENT_CACHE", "1").strip().lower() in (
        "0", "false", "no", "off",
    )


def install_consistent_compile_cache() -> bool:
    """Install the communicator-consistent ``make_so`` (idempotent).

    Must be called after ``import firedrake`` (or ``pyop2``) and before the
    first kernel is compiled on a multi-rank communicator. Returns True when
    the wrapper is active.
    """

    global _INSTALLED
    if _INSTALLED:
        return True
    if _disabled():
        return False

    try:
        from mpi4py import MPI
        import pyop2.compilation as compilation
        from pyop2 import mpi as pyop2_mpi
        from pyop2.caching import _as_hexdigest
        from pyop2.configuration import configuration

        original = compilation.make_so
        # make_so = mpi.collective(parallel_cache(...)(EventDecorator(body))):
        # unwrap down to the parallel_cache wrapper and take what it wraps --
        # the collective compile itself, without PyOP2's per-rank lookup.
        cached = original
        while getattr(getattr(cached, "__code__", None), "co_filename", "").endswith(
            os.path.join("pyop2", "caching.py")
        ) is False:
            cached = cached.__wrapped__
        body = cached.__wrapped__
        hashkey = compilation._make_so_hashkey
        disk_access = compilation.CompilerDiskAccess
        compilation_comm = pyop2_mpi.compilation_comm
        internal_comm = pyop2_mpi.internal_comm
        qualname = body.__qualname__
    except (ImportError, AttributeError) as exc:
        print(
            f"[ICESEE] WARNING: PyOP2 make_so layout not recognised ({exc!r}); "
            "the compiled-kernel cache is left unchanged and multi-rank model "
            "groups may deadlock on uneven cache contents.",
            file=sys.stderr,
            flush=True,
        )
        return False

    from ICESEE.src.utils.performance import record_phase

    def make_so(compiler, code, extension, comm, filename=None):
        icomm = internal_comm(comm, compiler)
        ccomm = compilation_comm(icomm, compiler)
        key = (_as_hexdigest(*hashkey(compiler, code, extension, comm)), qualname)
        cache = disk_access(configuration["cache_dir"], extension=".so")
        try:
            path = cache[key]
        except KeyError:
            path = None
        # Every rank of the compilation communicator takes the same branch.
        if ccomm.allreduce(path is not None, op=MPI.LAND):
            return path
        start = time.perf_counter()
        soname = body(compiler, code, extension, comm, filename)
        cache[key] = soname  # atomic mkstemp + rename; duplicate writes are harmless
        record_phase("jit_compile", time.perf_counter() - start)
        return cache[key]

    make_so.__wrapped__ = body
    make_so.__doc__ = original.__doc__
    compilation.make_so = make_so
    _INSTALLED = True
    return True
