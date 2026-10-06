# ==============================================================================
# @des: Runtime capability detection for Mode 3's state-storage backend
# selection (Stage 4D.3 scalability design). Pure Python, no Firedrake/
# mpi4py import required at module level, so it is importable and unit-
# testable on any machine regardless of what scientific stack is
# installed.
#
# This module answers exactly one question -- "what can this process's
# environment actually do" -- and answers nothing about which backend to
# USE for a given run; that policy decision belongs to whatever selects a
# backend (not yet implemented -- this is deliberately the first, narrow,
# independently-testable piece of Stage 4D.3's backend-neutral interface,
# per the explicit instruction not to force a large untested rewrite).
#
# Every check here is a real runtime probe, not a platform name/OS guess:
# macOS vs Linux, or "this looks like PACE", is never asked. A future
# backend selector should key off these capabilities, not off hostname/
# platform string matching -- that is precisely the HPC-implementation-
# coupled-to-one-platform's-quirks failure mode Stage 4D.3 explicitly
# rules out.
# ==============================================================================
from __future__ import annotations

from dataclasses import dataclass
import os


@dataclass(frozen=True)
class StateBackendCapabilities:
    """Snapshot of what this process's environment actually supports.

    ``mpi_world_size`` is 1 both for a genuinely serial run and for an
    MPI run of exactly one rank -- callers that need to distinguish
    those should check ``mpi_available`` separately.
    """

    h5py_mpi_enabled: bool
    mpi_available: bool
    mpi_world_size: int
    node_local_scratch: str | None

    def supports_collective_hdf5(self) -> bool:
        """Whether an MPI-collective (parallel) HDF5 file open is usable
        right now: h5py must have been built against an MPI-enabled HDF5,
        AND mpi4py must report more than one rank actually participating
        (a 1-rank "MPI" run gains nothing from and should not pay the
        collective-open code path's extra synchronization cost)."""

        return self.h5py_mpi_enabled and self.mpi_available and self.mpi_world_size > 1


def detect_state_backend_capabilities(comm=None) -> StateBackendCapabilities:
    """Probe this process's actual environment. No platform/hostname
    guessing anywhere in this function -- every field is a real runtime
    check, so the identical code path runs on macOS, a Linux workstation,
    or PACE and simply reports different (correct) capabilities.

    Args:
        comm: an ``mpi4py.MPI.Comm``-like object (must expose
            ``Get_size()``), or ``None``. Only imported/used if passed
            explicitly, so calling this with no arguments never requires
            mpi4py to be importable.
    """

    try:
        import h5py

        h5py_mpi_enabled = bool(h5py.get_config().mpi)
    except Exception:
        h5py_mpi_enabled = False

    mpi_available = False
    mpi_world_size = 1
    if comm is not None:
        try:
            mpi_world_size = int(comm.Get_size())
            mpi_available = True
        except Exception:
            mpi_available = False
            mpi_world_size = 1
    else:
        try:
            from mpi4py import MPI

            mpi_available = True
            mpi_world_size = int(MPI.COMM_WORLD.Get_size())
        except Exception:
            mpi_available = False
            mpi_world_size = 1

    # Node-local scratch: never assumed. Only reported if an explicit,
    # already-configured environment variable names one AND it exists on
    # disk right now. ICESEE code must treat this as optional and fall
    # back to the configured shared/working directory when absent --
    # never invented from a guessed path convention (Stage 4D.3's own
    # "do not assume node-local scratch exists" instruction).
    node_local_scratch = None
    for _env_name in ("ICESEE_NODE_LOCAL_SCRATCH", "TMPDIR", "SLURM_TMPDIR"):
        candidate = os.environ.get(_env_name)
        if candidate and os.path.isdir(candidate):
            node_local_scratch = candidate
            break

    return StateBackendCapabilities(
        h5py_mpi_enabled=h5py_mpi_enabled,
        mpi_available=mpi_available,
        mpi_world_size=mpi_world_size,
        node_local_scratch=node_local_scratch,
    )
