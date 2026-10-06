# ==============================================================================
# @des: Standalone (single-process, no MPI needed) worker for
# test_icepack_physical_nudge.py. Imports the REAL
# _physical_nudge_expr from the production synthetic_ice_stream
# _icepack_enkf.py module (not a reimplementation) and evaluates it on a
# small real Firedrake mesh, checking the physical-coordinate mask and
# taper properties Stage 4C's Option A fix is required to have:
#   - exactly zero for x/Lx > threshold_fraction (the mask boundary)
#   - a linear ramp from -amplitude at x=0 to 0 at x = Lx*threshold_fraction
#   - monotonically non-decreasing (less negative) as x increases inside
#     the nudged region
#   - amplitude=0 or threshold_fraction<=0 collapses to the zero field
#
# Run directly with plain python (no mpirun needed -- this is a
# single-rank R=1 mesh, only testing the pure taper formula, not
# decomposition behavior, which test_icepack_multirank_analysis.py's
# R1-vs-R2 tests already cover separately).
#
# Prints "TAPER_CHECK: PASS" and exits 0 on success, else prints the
# failing assertion and exits 1.
# ==============================================================================
import os
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
for _p in (str(_REPO_ROOT), str(_REPO_ROOT.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# _icepack_enkf.py imports ICESEE.config._utility_imports at module level,
# which parses sys.argv for --key=value CLI overrides; strip pytest's own
# argv down to just the program name so that import doesn't choke on
# unrecognized options.
sys.argv = [sys.argv[0]]

# ICESEE.config._utility_imports also looks for params.yaml relative to
# cwd at import time -- match run_da_icepack.py's own convention.
os.chdir(
    _REPO_ROOT
    / "applications"
    / "icepack_model"
    / "examples"
    / "synthetic_ice_stream"
)

import numpy as np
from firedrake import RectangleMesh, FunctionSpace, Function, SpatialCoordinate

from ICESEE.applications.icepack_model.examples.synthetic_ice_stream._icepack_enkf import (
    _physical_nudge_expr,
)

Lx, Ly = 50e3, 20e3
mesh = RectangleMesh(8, 4, Lx, Ly)
Q = FunctionSpace(mesh, "CG", 1)
x, y = SpatialCoordinate(mesh)

coords_fn = Function(FunctionSpace(mesh, "CG", 1)).interpolate(x)
xs = coords_fn.dat.data_ro.copy()

amplitude = 123.0
threshold_fraction = 0.3

expr = _physical_nudge_expr(x, Lx, amplitude, threshold_fraction)
field = Function(Q).interpolate(expr)
vals = field.dat.data_ro

failures = []

xi = xs / Lx
inside = xi <= threshold_fraction
outside = ~inside

# 1. Exactly zero outside the nudged region.
if not np.allclose(vals[outside], 0.0, atol=1e-9):
    failures.append(
        f"nonzero values outside threshold_fraction region: "
        f"max |v|={np.abs(vals[outside]).max():.6e}"
    )

# 2. Inside the region: value in [-amplitude, 0], exact linear formula.
expected_inside = -amplitude * (1.0 - xi[inside] / threshold_fraction)
if not np.allclose(vals[inside], expected_inside, atol=1e-6):
    failures.append(
        f"inside-region values do not match the linear taper formula: "
        f"max abs diff={np.abs(vals[inside] - expected_inside).max():.6e}"
    )
if vals[inside].min() < -amplitude - 1e-6:
    failures.append(f"taper exceeded -amplitude: min={vals[inside].min():.6e}")
if vals[inside].max() > 1e-6:
    failures.append(f"taper exceeded 0 at the inner edge: max={vals[inside].max():.6e}")

# 3. Monotonic: as x increases inside the region, the (negative) value
#    increases toward zero (i.e. magnitude decreases).
order = np.argsort(xs[inside])
sorted_vals = vals[inside][order]
diffs = np.diff(sorted_vals)
if np.any(diffs < -1e-9):
    failures.append(
        f"taper is not monotonically non-decreasing in x: "
        f"min diff={diffs.min():.6e}"
    )

# 4. amplitude=0 collapses to the zero field.
zero_amp_field = Function(Q).interpolate(
    _physical_nudge_expr(x, Lx, 0.0, threshold_fraction)
)
if not np.allclose(zero_amp_field.dat.data_ro, 0.0, atol=1e-12):
    failures.append("amplitude=0 did not collapse to the zero field")

# 5. threshold_fraction<=0 collapses to the zero field.
zero_thresh_field = Function(Q).interpolate(
    _physical_nudge_expr(x, Lx, amplitude, 0.0)
)
if not np.allclose(zero_thresh_field.dat.data_ro, 0.0, atol=1e-12):
    failures.append("threshold_fraction<=0 did not collapse to the zero field")

if failures:
    print("TAPER_CHECK: FAIL")
    for f in failures:
        print(" -", f)
    sys.exit(1)

print("TAPER_CHECK: PASS")
sys.exit(0)
