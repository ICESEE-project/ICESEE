# ==============================================================================
# @des: Tests for the generic resolved run schedule and CLI-override checks
#       (src/utils/run_schedule.py): analysis events counted from the final
#       time grid, early failure when a synthetic-observation run would do no
#       analysis, and command-line values that an application replaced.
# ==============================================================================
from __future__ import annotations

import numpy as np
import pytest

from ICESEE.src.utils import performance
from ICESEE.src.utils.run_schedule import resolve_run_schedule, verify_cli_overrides_respected


@pytest.fixture(autouse=True)
def _clean_performance_registry():
    performance.clear_registry()
    yield
    performance.clear_registry()


def _kwargs(num_years, dt, obs_start, obs_max, freq_obs=1.0, **extra):
    nt = int(round(num_years / dt))
    kwargs = {
        "num_years": num_years, "dt": dt, "nt": nt,
        "t": np.linspace(0, num_years, nt + 1),
        "obs_start_time": obs_start, "obs_max_time": obs_max, "freq_obs": freq_obs,
    }
    kwargs.update(extra)
    return kwargs


def test_schedule_is_rebuilt_from_the_final_time_grid(capsys):
    kwargs = _kwargs(10, 0.5, obs_start=2, obs_max=10, obs_index=np.array([200]))
    events = resolve_run_schedule(kwargs, rank=0, details={"truth count": 1})
    # Observations at years 2..10 fall on steps 4, 6, ..., 20; step 20 == nt
    # is never reached by a forecast loop over k < nt.
    assert list(kwargs["obs_index"]) == [4, 6, 8, 10, 12, 14, 16, 18, 20]
    assert events == 8
    out = capsys.readouterr().out
    assert "nt              = 20" in out and "len(t)          = 21" in out
    assert "analysis events = 8" in out and "truth count     = 1" in out


def test_production_schedule_is_unchanged():
    kwargs = _kwargs(82, 0.05, obs_start=10, obs_max=40)
    assert resolve_run_schedule(kwargs, rank=1) == 31
    assert kwargs["obs_index"][0] == 200 and kwargs["obs_index"][-1] == 800


def test_zero_analysis_synthetic_run_fails_before_model_work():
    kwargs = _kwargs(10, 0.5, obs_start=10, obs_max=40)  # only t=10, i.e. step nt
    with pytest.raises(ValueError, match="no analysis event"):
        resolve_run_schedule(kwargs, rank=1)


def test_zero_analysis_is_allowed_when_no_synthetic_observations_are_generated():
    kwargs = _kwargs(10, 0.5, obs_start=10, obs_max=40, generate_synthetic_obs=False)
    assert resolve_run_schedule(kwargs, rank=1) == 0


def test_mismatched_time_grid_is_rejected():
    kwargs = _kwargs(10, 0.5, obs_start=2, obs_max=10)
    kwargs["nt"] = 200
    with pytest.raises(ValueError, match="time grid"):
        resolve_run_schedule(kwargs, rank=1)


def test_replaced_cli_value_is_reported():
    verify_cli_overrides_respected({"dt": 0.5, "Nens": 4, "cli_overrides": {"dt": 0.5, "Nens": 4}})
    verify_cli_overrides_respected({"seed": 5.0, "cli_overrides": {"seed": 5}})
    with pytest.raises(ValueError, match=r"--dt=0\.5 was replaced by 0\.05"):
        verify_cli_overrides_respected({"dt": 0.05, "cli_overrides": {"dt": 0.5}})


def test_final_state_snapshot_is_counted_but_not_assimilated(capsys):
    # 2 yr at dt=0.05 with observations every 0.5 yr from 0.5 to 2: snapshots
    # at steps 10, 20, 30, 40; step 40 is the final state (nt), so 3 analyses.
    kwargs = _kwargs(2, 0.05, obs_start=0.5, obs_max=2, freq_obs=0.5)
    assert resolve_run_schedule(kwargs, rank=0) == 3
    out = capsys.readouterr().out
    assert "obs snapshots   = 4 (steps [10, 20, 30, 40])" in out
    assert "analysis events = 3 (steps [10, 20, 30]); a snapshot at the final state" in out
    assert performance._RUN_METADATA["observation_snapshots"] == 4
