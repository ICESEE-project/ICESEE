import random
import numpy as np
import pytest

from ICESEE.src.parallelization.forecast_recovery import (
    ForecastRecoveryPolicy, ForecastRecoveryExhausted, RecoverableForecastError,
    advance_with_recovery,
)


def test_failed_trial_restores_state_and_random_streams_and_lands_on_boundary():
    state = np.array([2.0])
    draws, times, events = [], [], []
    np.random.seed(71)
    random.seed(72)
    expected_np, expected_py = np.random.random(2), [random.random(), random.random()]
    np.random.seed(71)
    random.seed(72)

    def advance(time, dt):
        draws.append((np.random.random(), random.random()))
        times.append((time, dt))
        state[:] += 10 * dt
        if dt > .025:
            raise RecoverableForecastError("nonlinear solve rejected")

    outcome = advance_with_recovery(start_time=14.95, interval=.05,
        policy=ForecastRecoveryPolicy(3, .005), capture=lambda:state.copy(),
        restore=lambda saved:state.__setitem__(slice(None), saved),
        advance=advance, validate=lambda:None, report=events.append)
    assert outcome["substeps"] == 2 and outcome["attempts"] == 2
    assert outcome["end_time"] == 15.0
    np.testing.assert_allclose(state, [2.5])
    np.testing.assert_allclose(draws[1:], list(zip(expected_np, expected_py)))
    assert draws[0] == draws[1]
    assert times == [(14.95, .05), (14.95, .025), (14.975, .025)]
    assert [e["status"] for e in events] == ["rejected", "accepted"]


@pytest.mark.parametrize("kind", ["retry", "fatal"])
def test_rejection_restores_original_even_after_partial_progress(kind):
    state = [10.]
    def advance(time, dt):
        state[0] = -900
        if kind == "retry":
            raise RecoverableForecastError("stopped")
        raise FileNotFoundError("compiler header missing")
    error = ForecastRecoveryExhausted if kind == "retry" else RuntimeError
    with pytest.raises(error):
        advance_with_recovery(start_time=0, interval=1,
            policy=ForecastRecoveryPolicy(2, .1), capture=lambda:state[:],
            restore=lambda saved:state.__setitem__(slice(None), saved),
            advance=advance, validate=lambda:None)
    assert state == [10.]


def test_rank_local_rejection_is_agreed_before_any_rank_advances():
    class Comm:
        checks = 0
        def allgather(self, value):
            if value is None:
                self.checks += 1
            # Capture, initial restore, then the first completed local substep.
            if value is None and self.checks == 3:
                return [value, ("retry", "remote rank rejected")]
            return [value, value]
    comm = Comm()
    result = advance_with_recovery(start_time=0, interval=1,
        policy=ForecastRecoveryPolicy(2, .1), capture=lambda:0, restore=lambda x:None,
        advance=lambda *args:None, validate=lambda:None, comm=comm)
    assert result["substeps"] == 2


def test_cfl_bound_and_total_work_limits():
    times = []
    outcome = advance_with_recovery(start_time=2, interval=1,
        policy=ForecastRecoveryPolicy(4, .05, 8), stable_dt=.3,
        capture=lambda:0, restore=lambda x:None,
        advance=lambda t,dt:times.append((t,dt)), validate=lambda:None)
    assert outcome["substeps"] == 4
    assert times[-1] == (2.75, .25)
    with pytest.raises(ForecastRecoveryExhausted, match="0 attempts"):
        advance_with_recovery(start_time=0, interval=1,
            policy=ForecastRecoveryPolicy(4, .05, 2), stable_dt=.1,
            capture=lambda:0, restore=lambda x:None,
            advance=lambda *a:None, validate=lambda:None)


@pytest.mark.parametrize("raw", [True, {}, {"max_halvings":99,"min_dt":.1},
    {"max_halvings":2,"min_dt":0}, {"max_halvings":2.5,"min_dt":.1},
    {"max_halvings":2,"min_dt":.1,"unknown":1}])
def test_policy_rejects_unbounded_or_ambiguous_settings(raw):
    with pytest.raises((ValueError, TypeError)):
        ForecastRecoveryPolicy.from_options({"forecast_recovery":raw})


def test_recovery_is_disabled_unless_declared():
    assert ForecastRecoveryPolicy.from_options({}) is None
    assert ForecastRecoveryPolicy.from_options({"forecast_recovery":False}) is None


def test_rounded_external_interval_accepts_equal_minimum_dt():
    result = advance_with_recovery(start_time=.4, interval=.6-.4,
        policy=ForecastRecoveryPolicy(2,.05,4),stable_dt=.05,
        capture=lambda:None,restore=lambda saved:None,
        advance=lambda time,dt:None,validate=lambda:None)
    assert result['substeps']==4
