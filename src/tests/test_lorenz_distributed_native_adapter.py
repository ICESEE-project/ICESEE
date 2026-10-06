"""Unit tests for the production execution_mode 3 Lorenz96 native adapter."""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

# ``distributed_native_adapter.py`` imports ``_lorenz96_model.py``, which
# imports ``config._utility_imports`` -- a module that, at import time,
# parses ``sys.argv`` with ``argparse`` and reads ``params.yaml`` relative to
# the current working directory. Satisfy both here so this test runs under
# plain ``pytest`` from the repo root, following the same pattern used by
# ``test_icepack_idealized_pig_native.py`` for other application modules with
# this same import-time coupling.
_REPO_ROOT = Path(__file__).resolve().parents[2]
_LORENZ96_DIR = _REPO_ROOT / "applications" / "lorenz_model" / "examples" / "lorenz96"
for _extra_path in (str(_REPO_ROOT), str(_REPO_ROOT.parent)):
    if _extra_path not in sys.path:
        sys.path.insert(0, _extra_path)

_ARGV_BACKUP = sys.argv[:]
_CWD_BACKUP = os.getcwd()
sys.argv = [sys.argv[0], "--data_path", tempfile.mkdtemp(prefix="icesee_test_data_path_")]
os.chdir(_LORENZ96_DIR)
try:
    from applications.lorenz_model.lorenz_utils import (
        distributed_native_adapter as native,
    )
finally:
    sys.argv = _ARGV_BACKUP
    os.chdir(_CWD_BACKUP)

from src.parallelization.distributed_native_runtime import (
    validate_native_distributed_adapter,
)


class _Topology:
    def __init__(self, spatial_ranks=1):
        self.spatial_ranks = spatial_ranks


def _kwargs(**overrides):
    kwargs = dict(
        native_initial_ensemble={
            0: np.array([1.0, 1.0, 1.0]),
            1: np.array([2.0, 2.0, 2.0]),
        },
        dt=0.01,
        nt=20,
        sig_Q=[0.01, 0.01, 0.01],
        vec_inputs=["x", "y", "z"],
        sigma_96=10.0,
        beta_96=8.0 / 3.0,
        rho_96=28.0,
        seed=1,
        length_scale=[2, 2, 2],
        nd=3,
        num_state_vars=3,
        default_run=True,
        base_seed=42,
    )
    kwargs.update(overrides)
    return kwargs


def test_adapter_satisfies_native_distributed_protocol():
    validate_native_distributed_adapter(native.LorenzNativeAdapter())


def test_ar1_parameters_matches_serial_formula():
    dt, nt = 0.01, 2000
    alpha, rho = native.ar1_parameters({"dt": dt, "nt": nt})

    tau = max(500, max(200, dt))
    expected_alpha = 1 - dt / tau
    n = nt
    expected_rho = np.sqrt(
        (1 / dt)
        * ((1 - expected_alpha) ** 2)
        * (
            1
            / (
                n
                - (2 * expected_alpha)
                - (n * expected_alpha**2)
                + (2 * expected_alpha ** (n + 1))
            )
        )
    )
    assert alpha == pytest.approx(expected_alpha)
    assert rho == pytest.approx(expected_rho)


def test_initialize_native_member_requires_spatial_ranks_one():
    adapter = native.LorenzNativeAdapter()
    with pytest.raises(ValueError, match="spatial_ranks == 1"):
        adapter.initialize_native_member(
            0, topology=_Topology(spatial_ranks=2), icesee_kwargs=_kwargs()
        )


def test_initialize_native_member_packs_initial_state():
    adapter = native.LorenzNativeAdapter()
    member = adapter.initialize_native_member(
        1, topology=_Topology(), icesee_kwargs=_kwargs()
    )
    np.testing.assert_allclose(member.pack_owned(), [2.0, 2.0, 2.0])
    np.testing.assert_allclose(member.model_context["noise_state"], [0.0, 0.0, 0.0])


def test_forecast_native_member_gives_distinct_independent_member_noise():
    adapter = native.LorenzNativeAdapter()
    kwargs = _kwargs()
    topo = _Topology()
    m0 = adapter.initialize_native_member(0, topology=topo, icesee_kwargs=kwargs)
    m1 = adapter.initialize_native_member(1, topology=topo, icesee_kwargs=kwargs)

    adapter.forecast_native_member(m0, 0, topology=topo, icesee_kwargs=kwargs)
    adapter.forecast_native_member(m1, 0, topology=topo, icesee_kwargs=kwargs)

    # Different initial states already diverge under the noise-free dynamics,
    # but the persisted per-member AR(1) noise state must also differ (it is
    # keyed by member_id, not shared) -- this is the parity mechanism itself.
    assert not np.allclose(
        m0.model_context["noise_state"], m1.model_context["noise_state"]
    )
    assert not np.allclose(m0.pack_owned(), [1.0, 1.0, 1.0])


def test_forecast_native_member_is_reproducible_across_two_runs():
    adapter = native.LorenzNativeAdapter()
    kwargs = _kwargs()
    topo = _Topology()

    m_a = adapter.initialize_native_member(0, topology=topo, icesee_kwargs=kwargs)
    m_b = adapter.initialize_native_member(0, topology=topo, icesee_kwargs=kwargs)
    adapter.forecast_native_member(m_a, 5, topology=topo, icesee_kwargs=kwargs)
    adapter.forecast_native_member(m_b, 5, topology=topo, icesee_kwargs=kwargs)

    np.testing.assert_allclose(m_a.pack_owned(), m_b.pack_owned())
    np.testing.assert_allclose(
        m_a.model_context["noise_state"], m_b.model_context["noise_state"]
    )


def test_zero_sig_q_gives_no_noise_increment():
    adapter = native.LorenzNativeAdapter()
    kwargs = _kwargs(sig_Q=[])
    member = adapter.initialize_native_member(
        0, topology=_Topology(), icesee_kwargs=kwargs
    )
    increment = adapter._process_noise_increment(member, 0, kwargs)
    np.testing.assert_allclose(increment, [0.0, 0.0, 0.0])


def test_observe_native_member_returns_requested_rows():
    adapter = native.LorenzNativeAdapter()
    kwargs = _kwargs()
    member = adapter.initialize_native_member(
        0, topology=_Topology(), icesee_kwargs=kwargs
    )
    values = adapter.observe_native_member(
        member, np.array([0, 2]), topology=_Topology(), icesee_kwargs=kwargs
    )
    np.testing.assert_allclose(values, [1.0, 1.0])


def test_finalize_native_analysis_is_a_no_op():
    adapter = native.LorenzNativeAdapter()
    kwargs = _kwargs()
    member = adapter.initialize_native_member(
        0, topology=_Topology(), icesee_kwargs=kwargs
    )
    result = adapter.finalize_native_analysis(
        member,
        np.array([1.0, 1.0, 1.0]),
        0,
        topology=_Topology(),
        icesee_kwargs=kwargs,
    )
    assert result is None
