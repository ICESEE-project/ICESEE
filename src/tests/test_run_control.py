import json
import signal
from types import SimpleNamespace
import pytest
from ICESEE.src.parallelization.run_control import CycleStopController, atomic_json


COMM = SimpleNamespace(Get_rank=lambda:0, allgather=lambda value:[value])


def test_cycle_budget_only_signals_stop_after_a_committed_boundary():
    previous = signal.getsignal(signal.SIGUSR1)
    with CycleStopController({"run_segment_cycles":2}, COMM) as controller:
        assert not controller.after_commit(1)
        assert controller.after_commit(2)
        assert controller.reason == "cycle_budget"
    assert signal.getsignal(signal.SIGUSR1) == previous


def test_walltime_includes_initialization_and_preserves_signal_handlers():
    with CycleStopController({"run_walltime_seconds":10}, COMM,
                             elapsed=9, clock=lambda:100) as controller:
        assert not controller.after_commit(1)
        controller._signal(signal.SIGUSR1, None)
        assert controller.after_commit(2)
        assert controller.reason == "SIGUSR1"
    with CycleStopController({"run_walltime_seconds":10}, COMM,
                             elapsed=12, clock=lambda:100) as controller:
        assert controller.after_commit(1)
        assert controller.reason == "walltime_budget"


@pytest.mark.parametrize("options", [{"run_segment_cycles":0}, {"run_segment_cycles":1.1},
    {"run_walltime_seconds":float("nan")}, {"run_segment_cycles":2,"checkpoint_every":3}])
def test_invalid_limits_fail_before_handler_installation(options):
    with pytest.raises(ValueError):
        CycleStopController(options, COMM)


def test_report_is_atomically_replaced(tmp_path):
    path = tmp_path / "segment.json"
    atomic_json(path, {"next_step":2})
    atomic_json(path, {"next_step":4})
    assert json.loads(path.read_text()) == {"next_step":4}
    assert sorted(p.name for p in tmp_path.iterdir()) == ["segment.json"]
