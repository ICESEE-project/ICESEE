# ==============================================================================
# @des: Shared helper for tests that launch real multi-rank MPI jobs via
#       subprocess (Level-2 synthetic MPI tests and Level-3 real-model
#       integration tests -- see docs/testing.md for the tier definitions).
#       Not a test module itself; nothing here is collected by pytest.
# ==============================================================================
import shutil
import subprocess
import sys
from pathlib import Path


def _launcher_actually_works(candidate):
    """Functional compatibility check: does this launcher, paired with
    *this* Python's mpi4py, actually produce a real N-rank communicator?

    A version-string/vendor-name check is not reliable here: a
    PETSc/Firedrake build commonly bundles its own internal copy of Open
    MPI, so ``mpirun --version`` can genuinely say "Open MPI" while still
    being ABI-incompatible with the separately-built Open MPI that this
    Python's mpi4py links against (confirmed while building this helper:
    both "ranks" of a 2-rank launch identified as rank 0 with no real
    communicator between them, warned only via a "suspicious MPI execution
    environment" message on stderr -- the job silently ran as 2
    disconnected singletons instead of a coordinated 2-rank job). The only
    check that actually distinguishes a compatible launcher from this is
    running a real small job and checking every rank reports correctly.
    """
    probe = (
        "from mpi4py import MPI; c=MPI.COMM_WORLD; "
        "print(f'RANKPROBE {c.Get_rank()} {c.Get_size()}')"
    )
    try:
        result = subprocess.run(
            [candidate, "-n", "2", sys.executable, "-c", probe],
            capture_output=True, text=True, timeout=15,
        )
    except Exception:
        return False
    ranks_seen = set()
    for line in result.stdout.splitlines():
        if line.startswith("RANKPROBE"):
            _, rank_str, size_str = line.split()
            if size_str != "2":
                return False
            ranks_seen.add(rank_str)
    return ranks_seen == {"0", "1"}


def find_compatible_mpi_launcher():
    """Pick an mpirun/mpiexec whose MPI implementation actually matches the
    one this Python's mpi4py was built against.

    A sandbox/HPC environment can easily have more than one MPI
    installation on PATH (e.g. one bundled with a PETSc/Firedrake build,
    another from Homebrew or the system package manager). ``shutil.which``
    picks whichever comes first in PATH, with no guarantee it matches
    mpi4py's build. Any test that launches real MPI jobs via subprocess
    must use this instead of a bare ``shutil.which`` call, or it risks
    silently testing 2 disconnected singleton processes instead of a real
    multi-rank job -- passing for the wrong reason, which is exactly the
    kind of false result the MPI test suite exists to catch, not produce.

    Returns ``None`` if no working launcher is found (tests should skip,
    not fail, in that case -- this is an environment gap, not a code bug).
    """
    candidates = []
    for name in ("mpirun", "mpiexec"):
        found = shutil.which(name)
        if found and found not in candidates:
            candidates.append(found)
    for extra in ("/opt/homebrew/bin/mpirun", "/usr/local/bin/mpirun", "/usr/bin/mpirun"):
        if extra not in candidates and Path(extra).exists():
            candidates.append(extra)

    for candidate in candidates:
        if _launcher_actually_works(candidate):
            return candidate
    return None
