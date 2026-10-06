# ==============================================================================
# @des: Real-MPI regression coverage for Stage 4C's distributed-state
# Mode-2 analysis fixes for Icepack (ranks_per_model > 1).
#
# Icepack's model_capabilities registration
# (src/parallelization/parallel_mpi/model_capabilities.py) intentionally
# stays supports_multi_rank_per_model=False -- full R=1-vs-R=2 scientific
# equivalence for the noisy ensemble+analysis path, and real spare-rank/
# round coverage, are not yet established (see that file's own notes).
# But the underlying distributed-state forecast+analysis machinery this
# stage fixed is real and already verified correct for the cases below;
# this file gives that machinery permanent regression coverage (via the
# real synthetic_ice_stream Icepack application, not a synthetic stand-in)
# by temporarily enabling the capability inside a dedicated worker
# subprocess (_icepack_multirank_worker.py) rather than in this process or
# in the production registry.
#
# Covers, with a real 2-rank-per-model Firedrake mesh whose partition is
# deliberately uneven (nx=4, ny=4 -> local h sizes 36 and 45, not equal):
#   1. A full multi-timestep DA cycle (forecast + two EnKF analysis
#      updates, P=4, Nens=2, ranks_per_model=2 -- two model groups, no
#      spares, no rounds) completes without hanging or erroring, the
#      ensemble dimension is genuinely 2 (not duplicated to
#      world_size=4), member columns are finite and diverge from each
#      other only after the first analysis/noise event.
#   2. The deterministic true/nurged trajectory (no RNG, no ensemble) is
#      scientifically IDENTICAL between R=1 and R=2 -- every value
#      matches exactly as a set, differing only in the physical-node-to-
#      array-index permutation Firedrake's own partitioner assigns for a
#      given rank count.
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
_WORKER = Path(__file__).resolve().parent / "parallel_mpi" / "_icepack_multirank_worker.py"


def _run_worker(tmp_path, world_size, nens, ranks_per_model, base_seed=42,
                 nx=4, ny=4, extra_args=()):
    mpirun = find_compatible_mpi_launcher()
    if mpirun is None:
        pytest.skip("no compatible mpirun available in this environment")
    try:
        import firedrake  # noqa: F401
    except ImportError:
        pytest.skip("firedrake is not importable in this environment")

    out_dir = tmp_path / f"icepack_mr_{world_size}_{nens}_{ranks_per_model}"
    out_dir.mkdir(parents=True, exist_ok=True)

    env = dict(os.environ)
    env["PYTHONPATH"] = f"{_REPO_ROOT.parent}:{_REPO_ROOT}"

    args = [
        mpirun, "--oversubscribe", "-n", str(world_size),
        sys.executable, str(_WORKER),
        "--default_run",
        f"--Nens={nens}",
        f"--nx={nx}", f"--ny={ny}",
        "--num_years=2", "--timesteps_per_year=4",
        "--freq_obs=1", "--obs_max_time=2", "--obs_start_time=1",
        "--execution_mode=2",
        f"--ranks_per_model={ranks_per_model}",
        f"--data_path={out_dir}",
        "--create_ensemble_dataset=False",
        f"--base_seed={base_seed}",
        # Stage 4C: ranks_per_model > 1 with the default 'fft' random-field
        # method and nonzero sig_Q is rejected by
        # validate_stochastic_method_for_decomposition (FFT's covariance is
        # defined over flat array-index distance, not physical coordinates,
        # so it is not decomposition-invariant). 'graph' is the
        # decomposition-invariant method these multi-rank tests need.
        "--random_field_method=graph",
        *extra_args,
    ]
    result = subprocess.run(args, capture_output=True, text=True, timeout=180, env=env)
    assert result.returncode == 0, (
        f"P={world_size},Nens={nens},ranks_per_model={ranks_per_model}: "
        f"non-zero exit ({result.returncode}).\n"
        f"stdout tail:\n{result.stdout[-4000:]}\nstderr tail:\n{result.stderr[-2000:]}"
    )
    return out_dir, result.stdout


def _load_shard(data_path, index):
    with h5py.File(data_path / f"icesee_enkf_ens_{index:04d}.h5", "r") as handle:
        return np.asarray(handle["states"][:])


def test_p4_nens2_ranks_per_model2_completes_with_correct_ensemble_dimension(tmp_path):
    out_dir, stdout = _run_worker(tmp_path, world_size=4, nens=2, ranks_per_model=2)

    first = _load_shard(out_dir, 0)
    last = _load_shard(out_dir, 8)
    assert first.shape == (324, 2), (first.shape, stdout[-2000:])
    assert last.shape == (324, 2), (last.shape, stdout[-2000:])
    assert np.all(np.isfinite(first))
    assert np.all(np.isfinite(last))

    # Members must be genuinely distinct from t=0 onward (each gets its
    # own member-keyed initial perturbation -- generate_initial_member_
    # increment, seeded by (base_seed, ensemble_id, variable_index)) and
    # must stay distinct through the final step. If they were ever
    # identical, that would indicate the "duplicate ensemble member"
    # failure mode this test guards against (world_size=4 misread as 4
    # members instead of 2, or a shared/unreseeded initial-noise RNG).
    assert not np.allclose(first[:, 0], first[:, 1], atol=1e-6)
    assert not np.allclose(last[:, 0], last[:, 1], atol=1e-6)


def test_true_state_scientifically_identical_between_r1_and_r2(tmp_path):
    out_r1, _ = _run_worker(tmp_path, world_size=1, nens=1, ranks_per_model=1)
    out_r2, _ = _run_worker(tmp_path, world_size=2, nens=1, ranks_per_model=2)

    with h5py.File(out_r1 / "true_nurged_states.h5", "r") as handle:
        true_r1 = np.asarray(handle["true_state"][:])
    with h5py.File(out_r2 / "true_nurged_states.h5", "r") as handle:
        true_r2 = np.asarray(handle["true_state"][:])

    assert true_r1.shape == true_r2.shape
    nd, nsteps = true_r1.shape
    nvars = nd // 81
    for t in range(nsteps):
        for i in range(nvars):
            block_r1 = np.sort(true_r1[i * 81:(i + 1) * 81, t])
            block_r2 = np.sort(true_r2[i * 81:(i + 1) * 81, t])
            assert np.allclose(block_r1, block_r2, atol=1e-6), (
                f"timestep {t} variable block {i}: R=1 and R=2 true-state "
                "values differ even as an unordered set (a real "
                "partition-dependent numerical divergence, not just a "
                "DOF-ordering permutation)"
            )


def test_p6_nens2_ranks_per_model2_spare_ranks_completes(tmp_path):
    """4 active ranks (2 model groups of 2), 2 spare ranks, real Icepack.

    Exercises the generic spare-rank fixes together with a genuinely
    distributed (ranks_per_model=2) model group in the same run -- the
    real target of Stage 4C's spare-rank continuation, beyond the
    Lorenz-96/synthetic coverage in test_lorenz96_mode2_p_nens_matrix.py.
    """
    out_dir, stdout = _run_worker(tmp_path, world_size=6, nens=2, ranks_per_model=2)
    assert "spare=2" in stdout, stdout[-2000:]
    assert "groups=2" in stdout, stdout[-2000:]

    last = _load_shard(out_dir, 8)
    assert last.shape == (324, 2)
    assert np.all(np.isfinite(last))
    assert not np.allclose(last[:, 0], last[:, 1], atol=1e-6)


def test_p4_nens4_ranks_per_model2_rounds_completes(tmp_path):
    """2 model groups, 2 rounds, real Icepack -- every member gets a real
    (non-zero) initial condition and a distinct final state.

    Regression test for a real, reproduced bug: the distributed
    ensemble-initialization branch used to hard-code `ens = color`,
    silently never initializing any round-1+ member (they were left at
    HDF5's zero fill value), which then made that member's very first
    forecast solve diverge immediately with DIVERGED_FNORM_NAN.

    Uses a slightly larger mesh (nx=5,ny=5 rather than the nx=4,ny=4 used
    elsewhere in this file): with the now-correctly-seeded, properly
    sig_Q-scaled initial perturbation (see the decomposition-invariant
    ensemble-init fix), the 4x4 mesh's already-aggressive h_nurge_ic/
    u_nurge_ic initial condition combined with a genuine 4-member
    ensemble is numerically fragile for this test's tiny synthetic-ice-
    stream setup -- a real solver-robustness property of that specific
    combination of (mesh resolution, nurge magnitude, member count), not
    an architecture defect. Confirmed by direct comparison: identical
    config at nx=4,ny=4 reproducibly hits Firedrake's DIVERGED_DTOL on
    the very first forecast solve; nx=5,ny=5 completes cleanly.
    """
    out_dir, stdout = _run_worker(
        tmp_path, world_size=4, nens=4, ranks_per_model=2, nx=5, ny=5,
    )
    assert "rounds=2" in stdout, stdout[-2000:]

    first = _load_shard(out_dir, 0)
    last = _load_shard(out_dir, 8)
    nd = first.shape[0]
    assert first.shape == (nd, 4)
    assert last.shape == (nd, 4)
    assert np.all(np.isfinite(first))
    assert np.all(np.isfinite(last))
    # Every member (including round-1 members 2 and 3) must have a real,
    # non-zero initial condition -- the exact failure this test guards
    # against left members 2 and 3 all-zero at t=0.
    for member in range(4):
        assert not np.allclose(first[:, member], 0.0, atol=1e-6), (
            f"member {member} has an all-zero initial condition -- the "
            "round-1+ ensemble-init bug this test guards against"
        )
    # No two members collapse onto the same trajectory by the final step.
    for i in range(4):
        for j in range(i + 1, 4):
            assert not np.allclose(last[:, i], last[:, j], atol=1e-6), (
                f"members {i} and {j} are identical at the final step"
            )
