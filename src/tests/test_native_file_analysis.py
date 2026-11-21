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


def test_whole_member_initialization_uses_duration_independent_noise(monkeypatch):
    from ICESEE.src.parallelization import native_file_analysis as module
    world = SimpleNamespace(rank=0,size=1,Barrier=lambda:None)
    topology = SimpleNamespace(spatial_ranks=1,world=world,spatial_comm=object(),ensemble_slot=0)
    initial = SimpleNamespace(members={0:np.zeros(2)})
    pool = SimpleNamespace(snapshot_owned=lambda:initial)
    io = SimpleNamespace(H_matrix=lambda **kw:None, _ensure_batch=lambda step:None,
                         write_forecast=lambda *args:None)
    monkeypatch.setattr(module, "EnKF_fully_parallel_IO", lambda *args,**kw:io)
    monkeypatch.setattr(module, "apply_initial_spread_filebacked", lambda *args:None)
    results = []
    for nt in (100,1640):
        policy = FileBackedNativeAnalysis()
        monkeypatch.setattr(policy, "_restore", lambda *args:initial)
        options = dict(nd=2,Nens=1,nt=nt,dt=.05,data_path="unused",
                       process_noise_normalization="per_step")
        policy.initialize(pool,object(),topology=topology,icesee_kwargs=options)
        results.append((options["alpha"],options["rho"]))
    assert results[0] == results[1] == (.9999,1.)
