# ==============================================================================
# @des: Regression coverage for scripts/benchmarks/mode3_large_state_benchmark.py
# -- the synthetic Nx/Ne/P large-state scalability benchmark built to give
# ICESEE's distributed mode-3 architecture a serious, repeatable memory/scaling
# regression tool (motivated by the hypothetical 37 GB/member x 40-member =
# 1.48 TB stress case; see that script's own module docstring).
#
# These tests run the REAL benchmark script (which itself drives the REAL
# production distributed primitives -- distributed_native_runtime.py,
# distributed_native_cycle.py, distributed_analysis.py, distributed_checkpoint.py
# -- against a Firedrake-free synthetic model) at tiny, fast, CI-safe sizes.
# They exist to catch a future regression where some change accidentally
# introduces a full-Nx or full-Nx*Ne allocation on one rank -- the exact
# failure mode this whole benchmark was built to guard against.
# ==============================================================================
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from ICESEE.src.tests._mpi_launcher import find_compatible_mpi_launcher

_REPO_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT = _REPO_ROOT / "scripts" / "benchmarks" / "mode3_large_state_benchmark.py"

_MPIRUN = find_compatible_mpi_launcher()


def _skip_if_unavailable():
    if _MPIRUN is None:
        pytest.skip("no compatible mpirun available in this environment")


def _env():
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{_REPO_ROOT.parent}:{_REPO_ROOT}"
    env["PATH"] = f"/opt/homebrew/bin:{env.get('PATH', '')}"
    return env


def _run(world_size, extra_args, tmp_path, ckpt_name="ckpt"):
    ckpt = tmp_path / ckpt_name
    args = [
        _MPIRUN, "--oversubscribe", "-n", str(world_size),
        sys.executable, str(_SCRIPT),
        "--checkpoint-root", str(ckpt),
        *extra_args,
    ]
    result = subprocess.run(args, capture_output=True, text=True, timeout=120, env=_env())
    return result


def test_benchmark_refuses_dangerously_large_nx_without_opt_in():
    """The safety guard must refuse a --nx implying a full-global-array
    allocation beyond the default threshold, with no MPI needed."""
    result = subprocess.run(
        [sys.executable, str(_SCRIPT), "--nx", "999999999999", "--ne", "1"],
        capture_output=True, text=True, timeout=30, env=_env(),
    )
    assert result.returncode != 0
    assert "Refusing to run" in result.stderr


def test_benchmark_single_rank_owned_bytes_match_prediction(tmp_path):
    """P=1 (no spatial decomposition): owned bytes must equal the full Nx
    exactly (this rank owns everything) -- the trivial baseline case."""
    _skip_if_unavailable()
    result = _run(1, ["--nx", "10000", "--ne", "1", "--p-model", "1",
                      "--row-chunk-size", "512", "--nobs", "20", "--n-cycles", "1"], tmp_path)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-2000:]
    summary = json.loads(result.stdout)
    assert summary["peak_owned_array_bytes_this_rank"] == summary["predicted_owned_bytes_this_rank"]
    assert summary["peak_owned_array_bytes_this_rank"] == 10000 * 8


def test_benchmark_spatial_decomposition_bounds_owned_bytes_to_nx_over_p(tmp_path):
    """Genuine spatial decomposition (P_model=2): each rank's owned array
    must be bounded to (roughly) Nx/P_model -- never the full Nx -- and
    must exactly match this run's own prediction (the core "no hidden
    full-Nx allocation" regression guard)."""
    _skip_if_unavailable()
    result = _run(4, ["--nx", "40000", "--ne", "2", "--p-model", "2",
                      "--row-chunk-size", "512", "--nobs", "50", "--n-cycles", "2"], tmp_path)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-2000:]
    summary = json.loads(result.stdout)
    assert summary["peak_owned_array_bytes_this_rank"] == summary["predicted_owned_bytes_this_rank"]
    assert summary["peak_owned_array_elements_this_rank"] < 40000
    assert summary["peak_owned_array_elements_this_rank"] == 40000 // 2


def test_benchmark_checkpoint_files_written_and_restart_succeeds(tmp_path):
    """A full cycle including checkpoint write and this rank's own
    restart-read must complete without error and produce nonzero file/byte
    counts."""
    _skip_if_unavailable()
    result = _run(2, ["--nx", "5000", "--ne", "2", "--p-model", "1",
                      "--row-chunk-size", "256", "--nobs", "10", "--n-cycles", "1"], tmp_path)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-2000:]
    summary = json.loads(result.stdout)
    assert summary["checkpoint_files_total"] > 0
    assert summary["checkpoint_bytes_total"] > 0
    assert summary["t_restart_s"] is not None and summary["t_restart_s"] >= 0.0


def test_benchmark_analysis_workspace_grows_with_ensemble_size(tmp_path):
    """Empirical monotonicity check standing in for the full O(B*Ne + Ne^2)
    sweep documented in the module docstring: at fixed row_chunk_size, a
    larger Ne must not shrink the measured analysis-workspace peak."""
    _skip_if_unavailable()
    small = _run(2, ["--nx", "20000", "--ne", "2", "--p-model", "1",
                     "--row-chunk-size", "1024", "--nobs", "50", "--n-cycles", "1"],
                 tmp_path, ckpt_name="ckpt_small")
    assert small.returncode == 0, small.stdout[-3000:] + small.stderr[-2000:]
    large = _run(8, ["--nx", "20000", "--ne", "8", "--p-model", "1",
                     "--row-chunk-size", "1024", "--nobs", "50", "--n-cycles", "1"],
                 tmp_path, ckpt_name="ckpt_large")
    assert large.returncode == 0, large.stdout[-3000:] + large.stderr[-2000:]
    small_summary = json.loads(small.stdout)
    large_summary = json.loads(large.stdout)
    assert (
        large_summary["analysis_workspace_measured_traced_peak_bytes_max"]
        > small_summary["analysis_workspace_measured_traced_peak_bytes_max"]
    )
