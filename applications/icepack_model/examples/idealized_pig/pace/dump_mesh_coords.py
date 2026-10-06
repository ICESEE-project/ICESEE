# ==============================================================================
# @des: Companion to scripts/benchmarks/compare_pig_checkpoints.py
# (2026-09-28, PACE calibration prep, item E). Dumps this run's owned
# global-order physical mesh coordinates to a .npy file, reusing the EXACT
# same construction _icepack_native.py's _build_shared_context uses (so
# the coordinate ordering matches the checkpoint's own owned_start/
# owned_stop ranges for a P_model=1 shared context). Run once per
# invocation whose checkpoint you intend to compare across a SEPARATE
# invocation (different P_model, or a different backend) -- never
# imported by the production DA path.
#
# Usage: mpirun -n <P> python dump_mesh_coords.py --out coords.npy [run_da_icepack.py-style overrides]
# ==============================================================================
from __future__ import annotations

import sys
from pathlib import Path

_IDEALIZED_PIG_DIR = Path(__file__).resolve().parents[1]
_REPO_ROOT = _IDEALIZED_PIG_DIR.parents[3]
for _p in (str(_REPO_ROOT), str(_REPO_ROOT.parent), str(_IDEALIZED_PIG_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import argparse
import numpy as np
from mpi4py import MPI

_argv = sys.argv[1:]
_out = None
_rest = []
i = 0
while i < len(_argv):
    if _argv[i] == "--out":
        _out = _argv[i + 1]
        i += 2
    else:
        _rest.append(_argv[i])
        i += 1
sys.argv = [sys.argv[0]] + _rest

from ICESEE.config._utility_imports import icesee_kwargs
from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_native import (
    _shared_context,
)

world = MPI.COMM_WORLD
icesee_kwargs["comm"] = world
topology = type("Topo", (), {"spatial_comm": world})()

ctx = _shared_context(topology, icesee_kwargs)
coords = ctx.coords  # this rank's OWNED (x, y) pairs, shape (owned_size, 2)

all_coords = world.gather(coords, root=0)
if world.Get_rank() == 0:
    combined = np.concatenate(all_coords, axis=0)
    out_path = _out or "mesh_coords.npy"
    np.save(out_path, combined)
    print(f"[dump_mesh_coords] wrote {combined.shape} to {out_path}")
