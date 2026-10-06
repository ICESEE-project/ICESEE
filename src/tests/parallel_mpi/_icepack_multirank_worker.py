# ==============================================================================
# @des: Real-MPI worker for test_icepack_multirank_analysis.py.
#
# Runs the real synthetic_ice_stream Icepack application end to end (mesh,
# true/nurged state, synthetic obs, ensemble init, forecast, EnKF analysis,
# save) with ranks_per_model > 1, after temporarily overriding Icepack's
# model_capabilities registration to supports_multi_rank_per_model=True for
# THIS PROCESS ONLY. The production registration in model_capabilities.py
# stays False (see its own notes) until R=1-vs-R=2 scientific equivalence
# and spare/round coverage are established for the full noisy ensemble
# path; this worker exists to give the already-validated forecast+first-
# analysis distributed-state fixes (Stage 4C) permanent, real regression
# coverage without flipping that production flag.
#
# Usage: mpirun -n <P> python _icepack_multirank_worker.py [run_da_icepack.py CLI overrides...]
# ==============================================================================
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
for _p in (str(_REPO_ROOT), str(_REPO_ROOT.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from ICESEE.src.parallelization.parallel_mpi.model_capabilities import (
    register_model_capabilities,
)

register_model_capabilities(
    "icepack",
    supports_multi_rank_per_model=True,
    supports_distributed_state=True,
    notes=(
        "TEST-ONLY override (see "
        "src/tests/parallel_mpi/_icepack_multirank_worker.py); the "
        "production registration in model_capabilities.py stays False."
    ),
)

_SCRIPT = (
    _REPO_ROOT
    / "applications"
    / "icepack_model"
    / "examples"
    / "synthetic_ice_stream"
    / "run_da_icepack.py"
)
sys.argv = [str(_SCRIPT)] + sys.argv[1:]
exec(
    compile(_SCRIPT.read_text(), str(_SCRIPT), "exec"),
    {"__name__": "__main__", "__file__": str(_SCRIPT)},
)
