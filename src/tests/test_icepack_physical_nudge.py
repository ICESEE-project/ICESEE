# ==============================================================================
# @des: Real-MPI/Firedrake regression coverage for Stage 4C's Option A fix:
# physical x-coordinate-based nudging (_physical_nudge_expr in both
# applications/icepack_model/examples/synthetic_ice_stream/_icepack_enkf.py
# and applications/icepack_model/icepack_utils/_icepack_enkf.py), replacing
# the historical array-index-slice nudge that was decomposition-dependent
# (confirmed empirically to only have ~59% overlap with "the lowest-x
# nodes" even at R=1).
#
# Covers:
#   1. The taper/mask formula itself (via a lightweight standalone
#      Firedrake worker, no MPI/DA cycle needed).
#   2. Deterministic (no stochastic noise) R=1-vs-R=2 nudge equivalence,
#      to floating-point precision once DOF-ordering permutation is
#      accounted for.
#   3. Noisy (graph-method) initialization R=1-vs-R=2 equivalence, same
#      precision -- isolating that the stochastic contribution is also
#      decomposition-invariant, not just the deterministic nudge.
#   4. smb is untouched by the nudge (Stage 4C Part 5: only h/u/v are
#      affected).
#   5. A full small graph-method DA cycle stays R1-vs-R2 equivalent within
#      a generous tolerance that would still catch a real regression back
#      toward the old ~300-500% divergence, while tolerating the expected
#      small PETSc/Firedrake nonlinear-solver-path noise under different
#      decompositions (see model_capabilities.py's icepack notes).
# ==============================================================================
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

_REPO_ROOT = Path(__file__).resolve().parents[2]
_TAPER_WORKER = Path(__file__).resolve().parent / "parallel_mpi" / "_icepack_nudge_taper_worker.py"
_WORKER = Path(__file__).resolve().parent / "parallel_mpi" / "_icepack_multirank_worker.py"

_MPIRUN = find_compatible_mpi_launcher()


def _skip_if_unavailable():
    if _MPIRUN is None:
        pytest.skip("no compatible mpirun available in this environment")
    try:
        import firedrake  # noqa: F401
    except ImportError:
        pytest.skip("firedrake is not importable in this environment")


def _env():
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{_REPO_ROOT.parent}:{_REPO_ROOT}"
    env["PATH"] = f"/opt/homebrew/bin:{env.get('PATH', '')}"
    return env


def test_physical_nudge_taper_mask_and_shape():
    """_physical_nudge_expr: zero outside the mask, exact linear taper
    inside, monotonic, and correctly collapses to zero for
    amplitude=0 / threshold_fraction<=0."""
    _skip_if_unavailable()
    result = subprocess.run(
        [sys.executable, str(_TAPER_WORKER)],
        capture_output=True, text=True, timeout=120, env=_env(),
    )
    assert result.returncode == 0, (
        f"taper check failed:\nstdout:\n{result.stdout[-4000:]}\n"
        f"stderr:\n{result.stderr[-2000:]}"
    )
    assert "TAPER_CHECK: PASS" in result.stdout


def _run_worker(tmp_path, world_size, ranks_per_model, out_name, extra_args=()):
    out_dir = tmp_path / out_name
    out_dir.mkdir(parents=True, exist_ok=True)
    args = [
        _MPIRUN, "--oversubscribe", "-n", str(world_size),
        sys.executable, str(_WORKER),
        "--default_run",
        "--Nens=1",
        "--nx=4", "--ny=4",
        "--num_years=2", "--timesteps_per_year=4",
        "--freq_obs=1", "--obs_max_time=2", "--obs_start_time=1",
        "--execution_mode=2",
        f"--ranks_per_model={ranks_per_model}",
        f"--data_path={out_dir}",
        "--create_ensemble_dataset=False",
        *extra_args,
    ]
    result = subprocess.run(args, capture_output=True, text=True, timeout=180, env=_env())
    assert result.returncode == 0, (
        f"world_size={world_size},ranks_per_model={ranks_per_model}: "
        f"non-zero exit ({result.returncode}).\n"
        f"stdout tail:\n{result.stdout[-4000:]}\nstderr tail:\n{result.stderr[-2000:]}"
    )
    return out_dir


def _written_timesteps(mean_arr):
    return [t for t in range(mean_arr.shape[1]) if np.count_nonzero(mean_arr[:, t]) > 0]


def _load_mean(out_dir):
    with h5py.File(out_dir / "icesee_enkf_ens_mean.h5", "r") as f:
        return f["mean"][:]


def _per_variable_sorted_reldiff(m1, m2, t):
    """Coordinate/permutation-agnostic comparison: R=1 and R=2 assign
    physical mesh nodes to different flat array positions (a benign
    DOF-ordering artifact of Firedrake's own partitioner, confirmed
    separately in test_icepack_multirank_analysis.py's
    test_true_state_scientifically_identical_between_r1_and_r2), so a
    genuine apples-to-apples comparison must compare the SET of values
    per variable block, not raw flat-index position."""
    nd = m1.shape[0]
    q = nd // 4
    names = ["h", "u", "v", "smb"]
    out = {}
    for i, name in enumerate(names):
        a = np.sort(m1[i * q:(i + 1) * q, t])
        b = np.sort(m2[i * q:(i + 1) * q, t])
        d = np.abs(a - b)
        scale = max(np.abs(a).max(), 1e-8)
        out[name] = d.max() / scale
    return out


def test_deterministic_nudge_r1_r2_equivalence_to_machine_precision(tmp_path):
    """Stage 4C Part 9: with stochastic noise fully disabled (sig_Q=0),
    the physical nudge alone must produce scientifically identical
    R=1/R=2 fields (as a set, per variable) to floating-point precision,
    through the full deterministic forecast+analysis cycle -- confirming
    the nudge fix itself (not noise, not the nonlinear solver's
    decomposition-dependent path) is exactly decomposition-invariant."""
    _skip_if_unavailable()
    out_r1 = _run_worker(
        tmp_path, 1, 1, "det_r1",
        extra_args=["--sig_Q=[0,0,0,0]", "--random_field_method=graph"],
    )
    out_r2 = _run_worker(
        tmp_path, 2, 2, "det_r2",
        extra_args=["--sig_Q=[0,0,0,0]", "--random_field_method=graph"],
    )
    m1 = _load_mean(out_r1)
    m2 = _load_mean(out_r2)
    ts = _written_timesteps(m1)
    assert ts == _written_timesteps(m2)
    assert len(ts) >= 1
    for t in ts:
        rel = _per_variable_sorted_reldiff(m1, m2, t)
        for name, r in rel.items():
            assert r < 1e-8, (
                f"deterministic t={t} variable={name}: relative diff "
                f"{r:.3e} exceeds machine-precision tolerance -- the "
                "physical nudge is no longer exactly decomposition-"
                "invariant with noise disabled"
            )


def test_graph_noisy_init_r1_r2_equivalence_to_machine_precision(tmp_path):
    """Stage 4C Part 10: with the default (nonzero) sig_Q and
    random_field_method='graph', the COMBINED deterministic-nudge +
    stochastic-init initial ensemble mean must still be scientifically
    identical between R=1 and R=2 to floating-point precision -- the
    graph method's coordinate-keyed white noise (coordinate_keyed_white_
    noise) is decomposition-invariant by construction, and this confirms
    that holds end to end through real Icepack ensemble initialization."""
    _skip_if_unavailable()
    out_r1 = _run_worker(
        tmp_path, 1, 1, "noisy_r1",
        extra_args=["--random_field_method=graph"],
    )
    out_r2 = _run_worker(
        tmp_path, 2, 2, "noisy_r2",
        extra_args=["--random_field_method=graph"],
    )
    m1 = _load_mean(out_r1)
    m2 = _load_mean(out_r2)
    ts = _written_timesteps(m1)
    t0 = ts[0]
    rel = _per_variable_sorted_reldiff(m1, m2, t0)
    for name, r in rel.items():
        assert r < 1e-8, (
            f"noisy initial ensemble mean (t={t0}) variable={name}: "
            f"relative diff {r:.3e} exceeds machine-precision tolerance"
        )


def test_full_graph_da_cycle_r1_r2_stays_within_expected_tolerance(tmp_path):
    """Stage 4C Part 11: a full multi-timestep graph-method DA cycle
    (forecast + two analysis updates) stays R1-vs-R2 equivalent within a
    tolerance that tolerates expected PETSc/Firedrake nonlinear-solver-
    path noise under different decompositions but would still catch a
    real regression back toward the historical ~300-500% divergence this
    fix resolved."""
    _skip_if_unavailable()
    out_r1 = _run_worker(
        tmp_path, 1, 1, "full_r1",
        extra_args=["--random_field_method=graph"],
    )
    out_r2 = _run_worker(
        tmp_path, 2, 2, "full_r2",
        extra_args=["--random_field_method=graph"],
    )
    m1 = _load_mean(out_r1)
    m2 = _load_mean(out_r2)
    ts = _written_timesteps(m1)
    assert ts == _written_timesteps(m2)
    for t in ts:
        rel = _per_variable_sorted_reldiff(m1, m2, t)
        for name, r in rel.items():
            assert r < 0.10, (
                f"t={t} variable={name}: relative diff {r:.3e} exceeds "
                "the 10% regression-guard tolerance -- this is well "
                "above the ~0.1-3% expected solver-path noise observed "
                "during Stage 4C validation and likely indicates a real "
                "decomposition-invariance regression"
            )


def test_smb_unaffected_by_nudge_amplitude(tmp_path):
    """Stage 4C Part 5: the nudge only ever touches h/u/v -- smb must be
    bit-for-bit identical whether or not the nudge is active, since
    _physical_nudge_expr is never applied to the smb field."""
    _skip_if_unavailable()
    out_nudge_off = _run_worker(
        tmp_path, 1, 1, "smb_nudge_off",
        extra_args=[
            "--sig_Q=[0,0,0,0]", "--random_field_method=graph",
            "--h_nurge_ic=0", "--u_nurge_ic=0",
        ],
    )
    out_nudge_on = _run_worker(
        tmp_path, 1, 1, "smb_nudge_on",
        extra_args=["--sig_Q=[0,0,0,0]", "--random_field_method=graph"],
    )
    m_off = _load_mean(out_nudge_off)
    m_on = _load_mean(out_nudge_on)
    nd = m_off.shape[0]
    q = nd // 4
    t0 = _written_timesteps(m_off)[0]
    smb_off = m_off[3 * q:4 * q, t0]
    smb_on = m_on[3 * q:4 * q, t0]
    assert np.array_equal(smb_off, smb_on), (
        "smb changed when the nudge was toggled on -- the nudge must "
        "only ever affect h/u/v"
    )
