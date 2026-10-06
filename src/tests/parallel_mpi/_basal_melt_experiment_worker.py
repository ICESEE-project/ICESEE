# ==============================================================================
# @des: Worker for test_basal_melt_true_wrong_experiment.py. Prints
# JSON {"true": [melt_max_at_step0, melt_max_at_last_step],
#       "wrong": [melt_max_at_step0, melt_max_at_last_step]}
# computed via the ACTUAL BasalMeltRate() path (applications/icepack_model/
# examples/idealized_pig/_icepack_model.py) on a minimal UnitSquareMesh --
# never a duplicated formula. Single-rank; launched under mpirun -n 1 to
# match this repo's established real-Firedrake test convention (BasalMeltRate
# itself needs no multi-rank behavior, but Firedrake here is otherwise only
# ever exercised under mpirun in this test suite).
# ==============================================================================
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
_IDEALIZED_PIG_DIR = _REPO_ROOT / "applications" / "icepack_model" / "examples" / "idealized_pig"
for _extra_path in (str(_REPO_ROOT), str(_REPO_ROOT.parent), str(_IDEALIZED_PIG_DIR)):
    if _extra_path not in sys.path:
        sys.path.insert(0, _extra_path)

os.environ.setdefault("OMP_NUM_THREADS", "1")

import firedrake  # noqa: E402

from ICESEE.applications.icepack_model.examples.idealized_pig._icepack_model import (  # noqa: E402
    EXPERIMENT_TRUE,
    EXPERIMENT_WRONG,
    BasalMeltRate,
)

mesh = firedrake.UnitSquareMesh(2, 2)
Q = firedrake.FunctionSpace(mesh, "CG", 1)

floating = firedrake.Function(Q)
floating.dat.data[:] = 1.0  # entirely floating, so the melt mask never zeroes out melt_max's effect

s = firedrake.Function(Q)
s.dat.data[:] = -400.0  # draft = s - h; well inside the "control"/"warm" melt bands below

h = firedrake.Function(Q)
h.dat.data[:] = 100.0

# Mirrors params.yaml's icesee/features-convention Idealized PIG config
# (2026-09-28 reconciliation): num_years=82, timesteps_per_year=0.05,
# bmr_increase_time=0.
num_years = 82.0
dt = 0.05
bmr_increase_time = 0.0
nt = int(round(num_years / dt))

icesee_kwargs = {
    "dt": dt,
    "num_years": num_years,
    "bmr_increase_time": bmr_increase_time,
}

results = {}
for label, experiment in (("true", EXPERIMENT_TRUE), ("wrong", EXPERIMENT_WRONG)):
    _field0, melt_max_first = BasalMeltRate(
        icesee_kwargs, 0, floating, Q, s, h, scenario="control", experiment=experiment
    )
    _field1, melt_max_last = BasalMeltRate(
        icesee_kwargs, nt - 1, floating, Q, s, h, scenario="control", experiment=experiment
    )
    results[label] = [float(melt_max_first), float(melt_max_last)]

# Guards against the exact fragility this reconciliation removed: a caller
# passing a Python bool (or anything other than the two canonical strings)
# must fail loudly, never silently fall into the wrong-trajectory branch
# via a str/bool type mismatch.
try:
    BasalMeltRate(icesee_kwargs, 0, floating, Q, s, h, scenario="control", experiment=True)
    invalid_experiment_raises = False
except ValueError:
    invalid_experiment_raises = True
results["invalid_experiment_raises"] = invalid_experiment_raises

print(json.dumps(results))
