"""The CI comparison rejects member, clock and observation disagreements."""
import copy
import importlib.util
from pathlib import Path

import numpy as np
import pytest


def _checker():
    path = Path(__file__).resolve().parents[2] / "scripts/ci/run_lorenz96_ci.py"
    spec = importlib.util.spec_from_file_location("lorenz_ci_parity", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.validate_parallel_parity


def _history():
    t = np.linspace(0, 10, 1001)
    members = np.broadcast_to(np.sin(t), (3, 4, t.size)).copy()
    return dict(t=t, ind_m=np.array([10, 30, 50]),
                ensemble_true_state=np.sin(t)[None, :],
                w=np.ones((3, 3)), ensemble=members,
                ensemble_vec_mean=members.mean(axis=1))


def test_identical_histories_match():
    history = _history()
    metrics = _checker()(copy.deepcopy(history), history)
    assert metrics["ensemble"]["max_abs"] == 0


def test_compensating_member_errors_are_rejected_even_when_means_match():
    expected = _history()
    actual = copy.deepcopy(expected)
    actual["ensemble"][0, 0, 30] += 0.01
    actual["ensemble"][0, 1, 30] -= 0.01
    np.testing.assert_array_equal(actual["ensemble_vec_mean"], expected["ensemble_vec_mean"])
    with pytest.raises(AssertionError, match="assimilation-window"):
        _checker()(actual, expected)


def test_free_forecast_uses_trajectory_scale_at_a_zero_crossing():
    expected = _history()
    expected["ensemble"][..., 900] = 0
    expected["ensemble_vec_mean"][:, 900] = 0
    actual = copy.deepcopy(expected)
    actual["ensemble"][..., 900] = 5e-7
    actual["ensemble_vec_mean"][:, 900] = 5e-7
    _checker()(actual, expected)


def test_large_late_forecast_error_is_rejected():
    expected = _history()
    actual = copy.deepcopy(expected)
    actual["ensemble"][..., 900] += 0.01
    actual["ensemble_vec_mean"][:, 900] += 0.01
    with pytest.raises(AssertionError, match="full-history norm"):
        _checker()(actual, expected)


@pytest.mark.parametrize("key", ["t", "ind_m", "w", "ensemble_true_state"])
def test_changed_clock_observations_or_truth_are_rejected(key):
    expected = _history()
    actual = copy.deepcopy(expected)
    actual[key] = actual[key].copy()
    actual[key].flat[0] += 1
    with pytest.raises(AssertionError):
        _checker()(actual, expected)
