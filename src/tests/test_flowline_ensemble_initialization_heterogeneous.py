# ==============================================================================
# @des: Real, application-level regression test proving flowline_1d's
# scalar 'xg' variable is now handled correctly by Mode 0's ensemble
# initialization (src/EnKF/_ensemble_initialization.py), through the
# actual flowline model/EnKF functions -- not a mock.
#
# flowline_1d's default configuration (execution_mode: 0, vec_inputs=
# ['h','u','xg'], scalar_inputs=['xg']) previously crashed the very first
# time ensemble_initialization ran: the old scalar_inputs branch computed
# hdim as a naive average (state_size // total_state_param_vars) instead
# of the real per-variable block sizes (h:NX, u:NX, xg:1), producing a
# noise vector far too short to broadcast into the ensemble array.
# ==============================================================================
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_FLOWLINE_1D_DIR = _REPO_ROOT / "applications" / "flowline_model" / "examples" / "flowline_1d"
for _extra_path in (str(_REPO_ROOT), str(_REPO_ROOT.parent), str(_FLOWLINE_1D_DIR)):
    if _extra_path not in sys.path:
        sys.path.insert(0, _extra_path)

# config/_utility_imports.py (pulled in transitively) parses sys.argv, loads
# params.yaml relative to cwd, and deletes/recreates data_path at import
# time -- see test_flowline_mode3_runner.py for the identical pattern.
_SAFE_DATA_PATH = tempfile.mkdtemp(prefix="icesee_flowline_ensemble_init_test_")
_ARGV_BACKUP = sys.argv[:]
_CWD_BACKUP = os.getcwd()
sys.argv = [sys.argv[0], "--data_path", _SAFE_DATA_PATH]
os.chdir(_FLOWLINE_1D_DIR)
try:
    from applications.flowline_model.examples.flowline_1d import _flowline_enkf as model_module
    # _build_static_config is pure, rank-independent physical/grid-parameter
    # derivation (see its own docstring in mode3_runner.py) -- not mode-3
    # specific despite living in that file; reused here rather than
    # duplicating flowline's fairly involved parameter recipe.
    from applications.flowline_model.examples.flowline_1d.mode3_runner import (
        _build_static_config,
    )
finally:
    sys.argv = _ARGV_BACKUP
    os.chdir(_CWD_BACKUP)

from ICESEE.src.EnKF._ensemble_initialization import ensemble_initialization


def _small_flowline_kwargs(tmp_path, nens=3, sig_q=0.02):
    tmp_path.mkdir(parents=True, exist_ok=True)
    kwargs = dict(
        model_module=model_module,
        Nens=nens,
        num_years=2.0,
        num_state_vars=3,
        num_param_vars=0,
        total_state_param_vars=3,
        vec_inputs=["h", "u", "xg"],
        scalar_inputs=["xg"],
        sig_Q=[sig_q, sig_q, sig_q],
        base_seed=42,
        data_path=str(tmp_path),
        default_run=True,
        even_distribution=False,
        joint_estimation=False,
        localization_flag=False,
        # Small grid for a fast real solve.
        N1=4,
        N2=2,
        hscale=1000.0,
        A=4e-12,
        n=3,
        C=3e6,
        rho_ice=900.0,
        rho_water=1000.0,
        g=9.81,
        accum=0.65,
        facemelt=5,
        sigGZ=0.97,
        transient=0,
        tcurrent=1,
        xsill=50e3,
        sillamp=500,
        sillsmooth=1e-5,
        year=31536000,
        seed=1,
    )
    _build_static_config(kwargs)
    return kwargs


def test_flowline_ensemble_initialization_does_not_crash_and_shapes_correctly(tmp_path):
    kwargs = _small_flowline_kwargs(tmp_path)
    NX = kwargs["NX"]
    nd = kwargs["nd"]
    assert nd == 2 * NX + 1  # h:NX + u:NX + xg:1

    result = ensemble_initialization(**kwargs)
    icesee_kwargs_out, ensemble_vec, *_ = result

    assert ensemble_vec.shape == (nd, kwargs["Nens"])
    assert np.all(np.isfinite(ensemble_vec))


def test_flowline_xg_block_is_a_single_perturbed_scalar_per_member(tmp_path):
    kwargs = _small_flowline_kwargs(tmp_path, nens=4)
    NX = kwargs["NX"]
    _, ensemble_vec, *_ = ensemble_initialization(**kwargs)

    xg_row = ensemble_vec[2 * NX, :]  # the single 'xg' entry, one per member
    # Members must not all collapse to the identical xg value (the exact
    # symptom of the old shared/unreseeded-RNG + wrong-block-size bug).
    assert len(set(np.round(xg_row, 12))) > 1
    assert np.all(np.isfinite(xg_row))


def test_flowline_h_and_u_blocks_have_nonzero_spread(tmp_path):
    kwargs = _small_flowline_kwargs(tmp_path, nens=4)
    NX = kwargs["NX"]
    _, ensemble_vec, *_ = ensemble_initialization(**kwargs)

    h_block = ensemble_vec[:NX, :]
    u_block = ensemble_vec[NX:2 * NX, :]
    assert np.all(h_block.std(axis=1) > 0.0)
    assert np.all(u_block.std(axis=1) > 0.0)


def test_flowline_zero_sig_q_leaves_xg_at_the_noise_free_value(tmp_path):
    kwargs_perturbed = _small_flowline_kwargs(tmp_path / "perturbed", nens=2, sig_q=0.05)
    kwargs_zero = _small_flowline_kwargs(tmp_path / "zero", nens=2, sig_q=0.0)

    _, ensemble_perturbed, *_ = ensemble_initialization(**kwargs_perturbed)
    _, ensemble_zero, *_ = ensemble_initialization(**kwargs_zero)

    NX = kwargs_zero["NX"]
    # With sig_Q=0 every member's xg (and h/u) increment is exactly zero,
    # so both ensemble members must be identical at the xg row.
    xg_zero = ensemble_zero[2 * NX, :]
    np.testing.assert_allclose(xg_zero, xg_zero[0])
    # The perturbed run must actually differ member to member.
    xg_perturbed = ensemble_perturbed[2 * NX, :]
    assert not np.allclose(xg_perturbed, xg_perturbed[0])
