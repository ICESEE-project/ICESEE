from types import SimpleNamespace
import numpy as np
from ICESEE.src.parallelization.native_file_analysis import FileBackedNativeAnalysis
from ICESEE.src.parallelization.distributed_whole_member_adapter import build_whole_member_field_registry
from ICESEE.src.parallelization.distributed_runtime import LocalMemberEnsemble


def test_cycle_publishes_members_and_restores_actual_analysis(monkeypatch):
    layout = build_whole_member_field_registry({"h": np.zeros(2)}).layout
    forecast = LocalMemberEnsemble(layout, {3: np.array([2.,4.])})
    analyzed = LocalMemberEnsemble(layout, {3: np.array([5.,7.])})
    events = []
    transform = np.array([[.8,.2],[.2,.8]])
    io = SimpleNamespace(
        _ensure_batch=lambda step: events.append(("open",step)),
        write_forecast=lambda step,values,member: events.append(("write",step,member,values.copy())),
        compute_analysis_update=lambda **kw: events.append(("analysis",kw["k"])),
        last_analysis_transform=transform)
    policy = FileBackedNativeAnalysis(); policy.io = io
    monkeypatch.setattr(policy,"_restore",lambda pool,step: analyzed)
    pool = SimpleNamespace(forecast_with_adapter=lambda *args,**kw: forecast)
    topology = SimpleNamespace(world=SimpleNamespace(Barrier=lambda: events.append("barrier")))
    result = policy(pool,object(),4,[SimpleNamespace(observation_ids=[0,1])],
                    number_of_batches=1,topology=topology,
                    icesee_kwargs={"Nens":2,"k":4},error_mode="legacy_prior_anomalies")
    assert result.local_analysis is analyzed
    np.testing.assert_equal(result.transform,transform)
    assert result.observation_rows == 2
    assert events[0] == ("open",5)
    assert events[1][1:3] == (5,3)
    assert events[2:] == ["barrier",("analysis",4)]


def test_forecast_only_cycle_does_not_analyze_or_restore(monkeypatch):
    layout = build_whole_member_field_registry({"h": np.zeros(1)}).layout
    forecast = LocalMemberEnsemble(layout,{0:np.array([9.])})
    policy = FileBackedNativeAnalysis()
    policy.io = SimpleNamespace(_ensure_batch=lambda step:None, write_forecast=lambda *args:None)
    pool = SimpleNamespace(forecast_with_adapter=lambda *args,**kw:forecast)
    topology = SimpleNamespace(world=SimpleNamespace(Barrier=lambda:None))
    result = policy(pool,object(),2,[],number_of_batches=0,topology=topology,
                    icesee_kwargs={"Nens":1},error_mode="legacy_prior_anomalies")
    assert result.local_analysis is forecast
    assert result.observation_rows == 0
    np.testing.assert_equal(result.transform,np.eye(1))
