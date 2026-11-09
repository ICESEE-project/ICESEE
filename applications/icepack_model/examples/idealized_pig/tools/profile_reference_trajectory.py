#!/usr/bin/env python3
# Time reference-trajectory setup, solves and I/O on one MPI rank.
# Usage: python tools/profile_reference_trajectory.py --output DIR [--steps N]
# Refuses an existing non-empty output directory.

import os

# Must precede NumPy: its BLAS reads the thread count when loaded.
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import hashlib
import json
import platform
import socket
import subprocess
import sys
import time
from pathlib import Path

_T_START = time.perf_counter()

PIG_DIR = Path(__file__).resolve().parents[1]


def _parse():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--dt", type=float, default=None,
                        help="model step in years (default: params.yaml timesteps_per_year)")
    parser.add_argument("--num_years", type=float, default=None,
                        help="forcing-ramp horizon (default: steps*dt, as run_da_icepack.py)")
    parser.add_argument("--experiment", choices=("true", "wrong"), default="true")
    parser.add_argument("--params", type=Path, default=PIG_DIR / "params.yaml")
    return parser.parse_args()


ARGS = _parse()

import numpy as np
import yaml
from mpi4py import MPI

if MPI.COMM_WORLD.Get_size() != 1:
    sys.exit("profile_reference_trajectory.py is serial; launch it without mpirun.")

out = ARGS.output.resolve()
if out.exists() and any(out.iterdir()):
    sys.exit(f"refusing to write into non-empty {out}")
out.mkdir(parents=True, exist_ok=True)

os.chdir(PIG_DIR)
sys.path.insert(0, str(PIG_DIR))

_t = time.perf_counter()
import firedrake
from firedrake.petsc import PETSc
import icepack
import h5py
from ICESEE.src.utils import performance
from ICESEE.src.utils.pyop2_compile_cache import install_consistent_compile_cache
install_consistent_compile_cache()
from ICESEE.applications.icepack_model.examples.idealized_pig import _icepack_model as pig
from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_native import (
    _basal_melt_for_step,
)
from ICESEE.applications.icepack_model.examples.idealized_pig._mode3_profiles import (
    profile_columns,
)
import ICESEE.applications.icepack_model.examples.idealized_pig.modelfunc as mf
IMPORT_S = time.perf_counter() - _t

# Recording
PHASES = {}
SOLVES = []
FLOW_SOLVER = None


def _jit_total():
    return performance._RECORDED_PHASES.get("jit_compile", [0.0, 0])[0]


class phase:
    """Accumulate wall time and JIT-compile time of one named phase."""

    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.t0, self.j0 = time.perf_counter(), _jit_total()
        return self

    def __exit__(self, *exc):
        entry = PHASES.setdefault(self.name, {"seconds": 0.0, "jit_seconds": 0.0, "calls": 0})
        entry["seconds"] += time.perf_counter() - self.t0
        entry["jit_seconds"] += _jit_total() - self.j0
        entry["calls"] += 1
        self.elapsed = time.perf_counter() - self.t0
        return False


_original_solve = firedrake.NonlinearVariationalSolver.solve


def _recording_solve(self, *args, **kwargs):
    snes = self.snes
    snes.setConvergenceHistory(length=10000, reset=True)
    record = {"step": CURRENT["step"]}
    rtol, atol, stol, max_it = snes.getTolerances()
    record["tolerances"] = {"rtol": rtol, "atol": atol, "stol": stol, "max_it": max_it}
    t0, j0 = time.perf_counter(), _jit_total()
    error = None
    try:
        return _original_solve(self, *args, **kwargs)
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}".splitlines()[0]
        raise
    finally:
        history = snes.getConvergenceHistory()[0]
        record.update(
            role=_role_of(self),
            seconds=time.perf_counter() - t0,
            jit_seconds=_jit_total() - j0,
            iterations=int(snes.getIterationNumber()),
            linear_iterations=int(snes.getLinearSolveIterations()),
            reason=int(snes.getConvergedReason()),
            residual_history=[float(v) for v in history],
            error=error,
        )
        SOLVES.append(record)


firedrake.NonlinearVariationalSolver.solve = _recording_solve
CURRENT = {"step": None}


def _role_of(nonlinear_solver):
    """Which Icepack solver owns this Firedrake solver (created lazily)."""
    for role, attr in (("diagnostic", "_diagnostic_solver"), ("prognostic", "_prognostic_solver")):
        if getattr(getattr(FLOW_SOLVER, attr, None), "_solver", None) is nonlinear_solver:
            return role
    return "other"


def _git(*args):
    try:
        return subprocess.run(["git", "-C", str(PIG_DIR), *args], capture_output=True,
                              text=True, timeout=30).stdout.strip()
    except Exception:
        return None


def _sha256(path, limit=64 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        digest.update(handle.read(limit))
    return digest.hexdigest()[:16]


def _field_stats(h, u, s):
    hv, uv, sv = h.dat.data_ro, u.dat.data_ro, s.dat.data_ro
    speed = np.hypot(uv[:, 0], uv[:, 1])
    return {
        "h_min": float(hv.min()), "h_max": float(hv.max()), "h_mean": float(hv.mean()),
        "s_max": float(sv.max()), "speed_max": float(speed.max()),
        "finite": bool(np.isfinite(hv).all() and np.isfinite(uv).all() and np.isfinite(sv).all()),
    }


# Configuration
with open(ARGS.params) as handle:
    params = yaml.safe_load(handle)
kw = {}
for section in params.values():
    kw.update(section)
dt = float(ARGS.dt if ARGS.dt is not None else kw["timesteps_per_year"])
nsteps = int(ARGS.steps)
num_years = float(ARGS.num_years) if ARGS.num_years is not None else nsteps * dt
kw.update(dt=dt, num_years=num_years, nt=nsteps, comm=MPI.COMM_SELF,
          bmr_increase_time=int(kw["bmr_increase_time"]))
experiment = pig.EXPERIMENT_TRUE if ARGS.experiment == "true" else pig.EXPERIMENT_WRONG

report = {
    "command": [sys.executable, *sys.argv],
    "host": socket.gethostname(),
    "platform": platform.platform(),
    "source": {"pig_dir": str(PIG_DIR), "revision": _git("rev-parse", "HEAD"),
               "dirty_files": (_git("status", "--short", "--", ".") or "").splitlines()},
    "versions": {
        "python": platform.python_version(), "numpy": np.__version__, "h5py": h5py.__version__,
        "firedrake": getattr(firedrake, "__version__", None), "icepack": getattr(icepack, "__version__", None),
        "petsc": ".".join(map(str, PETSc.Sys.getVersion())),
        "mpi": MPI.Get_library_version().splitlines()[0],
    },
    "environment": {k: os.environ.get(k) for k in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS",
        "PYOP2_CACHE_DIR", "FIREDRAKE_TSFC_KERNEL_CACHE_DIR", "TMPDIR")},
    "inputs": {name: {"path": kw[name], "bytes": os.path.getsize(kw[name]),
                      "sha256_first_64MiB": _sha256(kw[name])}
               for name in ("initFile", "meshFile", "SMBFile")},
    "settings": {"steps": nsteps, "dt": dt, "num_years": num_years, "experiment": experiment,
                 "compact_initialization": kw.get("compact_initialization"),
                 "save_steps": kw.get("save_steps")},
    "import_seconds": IMPORT_S,
}

# Setup
with phase("checkpoint_mesh_load"):
    mesh_from_checkpoint = mf.getMeshFromCheckPoint(kw["initFile"], kw)
with phase("argus_mesh_and_spaces"):
    mesh, Q, V, _opts = mf.setupMesh(kw["meshFile"], degree=1, meshOversample=2,
                                     newMesh=mesh_from_checkpoint, comm=MPI.COMM_SELF)
report["mesh"] = {"nodes": int(Q.dim()), "cells": int(mesh.num_cells())}
with phase("solver_construction"):
    model = icepack.models.IceStream(friction=pig.schoofFriction, viscosity=pig.regViscosity)
    solver_options = {"dirichlet_ids": [1],
                      "diagnostic_solver_parameters": {"max_iterations": 150, "tolerance": 1e-6}}
    flow_solver = icepack.solvers.FlowSolver(model, **solver_options)
FLOW_SOLVER = flow_solver
report["configured_solver_options"] = solver_options

CURRENT["step"] = "initial"
with phase("checkpoint_field_load_and_initial_diagnostic_solve"):
    h, h0, s, s0, u, bed, zF, grounded, floating, A0, beta0, smb = pig.initializeRun(
        kw, flow_solver, mesh, Q, V)
kw.update(Q=Q, V=V, h=h, h0=h0, s=s, s0=s0, u=u, bed=bed, A0=A0, beta0=beta0,
          floating=floating, grounded=grounded, smb=smb)
initial = _field_stats(h, u, s)
report["initial_state"] = initial

state_file = out / "reference_state.h5"
n = int(h.dat.data_ro.size)
with h5py.File(state_file, "w") as handle:
    handle.create_dataset("state", shape=(5 * n, nsteps + 1), dtype="f8", chunks=(n, 1))
points = None
save_steps = tuple(kw.get("save_steps") or ())

# Trajectory
STEPS = []
first_abnormal = None
failure = None
for k in range(nsteps):
    CURRENT["step"] = k
    row = {"step": k, "model_time_after": (k + 1) * dt}
    solves_before = len(SOLVES)
    try:
        with phase("basal_melt_update") as p:
            melt = _basal_melt_for_step(kw, k, floating, Q, s, h, experiment=experiment)
        row["melt_s"] = p.elapsed
        with phase("icepack_step_total") as p:
            h, u, s, floating, grounded = pig.Icepack(flow_solver, h, u, smb, melt, bed, dt, h0, kw)
        row["icepack_step_s"] = p.elapsed
    except Exception as exc:
        failure = {"step": k, "error": f"{type(exc).__name__}: {exc}".splitlines()[0]}
    step_solves = SOLVES[solves_before:]
    for role in ("prognostic", "diagnostic"):
        mine = [r for r in step_solves if r["role"] == role]
        row[f"{role}_s"] = sum(r["seconds"] for r in mine)
        row[f"{role}_iterations"] = [r["iterations"] for r in mine]
    if failure:
        STEPS.append(row)
        break
    row["geometry_update_s"] = row["icepack_step_s"] - row["prognostic_s"] - row["diagnostic_s"]
    with phase("hdf5_state_write") as p:
        with h5py.File(state_file, "a") as handle:
            data = handle["state"]
            data[0:n, k + 1] = h.dat.data_ro
            data[n:2 * n, k + 1] = u.dat.data_ro[:, 0]
            data[2 * n:3 * n, k + 1] = u.dat.data_ro[:, 1]
            data[3 * n:4 * n, k + 1] = s.dat.data_ro
            data[4 * n:5 * n, k + 1] = melt.dat.data_ro
    row["hdf5_write_s"] = p.elapsed
    if save_steps and profile_columns(k + 1, dt, save_steps):
        with phase("profile_sampling") as p:
            if points is None:
                points = pig.initial_flowline_profile(kw)[2]
            pig.flowline_profile(h, s, points)
        row["profile_s"] = p.elapsed
    row.update(_field_stats(h, u, s))
    STEPS.append(row)
    if first_abnormal is None:
        reasons = []
        if not row["finite"]:
            reasons.append("non-finite state")
        if row["h_max"] > 1.5 * initial["h_max"]:
            reasons.append(f"h_max {row['h_max']:.4g} > 1.5x initial {initial['h_max']:.4g}")
        if row["speed_max"] > 10.0 * initial["speed_max"]:
            reasons.append(f"speed_max {row['speed_max']:.4g} > 10x initial")
        max_it = max((r["tolerances"]["max_it"] for r in step_solves), default=0)
        if any(r["iterations"] >= 0.8 * max_it for r in step_solves if r["role"] == "diagnostic"):
            reasons.append(f"diagnostic solve used >= 80% of max_it={max_it}")
        if reasons:
            first_abnormal = {"step": k, "reasons": reasons}
    print(f"[profile] step {k:4d} t={row['model_time_after']:.3f} total={row['icepack_step_s']:.2f}s "
          f"prog={row['prognostic_s']:.2f}s{row['prognostic_iterations']} "
          f"diag={row['diagnostic_s']:.2f}s{row['diagnostic_iterations']} "
          f"h_max={row['h_max']:.1f} speed_max={row['speed_max']:.1f}", flush=True)

# Report
report.update(
    phases=PHASES,
    jit_compile_seconds_total=_jit_total(),
    steps=STEPS,
    solves=SOLVES,
    first_abnormal_step=first_abnormal,
    failure=failure,
    wall_seconds=time.perf_counter() - _T_START,
)
done = [r for r in STEPS if "geometry_update_s" in r]
if done:
    later = done[1:] or done
    report["per_step_summary_excluding_first"] = {
        key: {"mean": float(np.mean([r[key] for r in later])), "max": float(np.max([r[key] for r in later]))}
        for key in ("icepack_step_s", "prognostic_s", "diagnostic_s", "geometry_update_s",
                    "melt_s", "hdf5_write_s")
    }
    report["diagnostic_iterations"] = [it for r in done for it in r["diagnostic_iterations"]]
with open(out / "profile.json", "w") as handle:
    json.dump(report, handle, indent=1)
print(f"[profile] wrote {out / 'profile.json'} (wall {report['wall_seconds']:.1f}s, "
      f"first abnormal step: {first_abnormal}, failure: {failure})", flush=True)
sys.exit(1 if failure else 0)
