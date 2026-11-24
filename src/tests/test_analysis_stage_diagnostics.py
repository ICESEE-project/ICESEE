"""Same-cycle diagnostics preserve the update and distinguish finalizer changes."""
import json
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from src.tests.test_distributed_native_cycle import _Adapter, _member, _topology
from src.parallelization.distributed_native_runtime import NativeDistributedMemberPool
from src.parallelization.distributed_native_cycle import NativeObservationBatch, run_native_global_analysis_cycle
from src.parallelization.analysis_diagnostics import (
    NativeAnalysisDiagnostics, prepare_analysis_diagnostics, publish_analysis_diagnostics,
)
from src.parallelization.distributed_adapter import DistributedBlockStateLayout, DistributedStateBlockLayout
from src.utils.ensemble_spread import relax_analysis_spread


def run_cycle(diagnostics):
    pool = NativeDistributedMemberPool({0: _member(0), 1: _member(1), 2: _member(2)})
    class Balance(_Adapter):
        def finalize_native_analysis(self, member, forecast_owned, timestep, **kwargs):
            member.unpack_owned(member.pack_owned()+7)
    adapter = Balance(); topology = _topology()
    observer = NativeAnalysisDiagnostics(pool.layout, 3, topology, np.arange(5.)) if diagnostics else None
    def controls(analysis_vec, forecast_vec, global_rows, **kwargs):
        return relax_analysis_spread(analysis_vec, forecast_vec, global_rows, {"state": (0, 5)},
                                      {"analysis_spread_relaxation": .5})
    result = run_native_global_analysis_cycle(pool, adapter, 1,
        [NativeObservationBatch(np.array([0, 4]), np.array([2., 7.]))], number_of_batches=1,
        topology=topology, icesee_kwargs={"Nens": 3}, error_mode="legacy_prior_anomalies",
        state_row_chunk_size=2, analysis_row_controls=controls, analysis_diagnostics=observer)
    return result


def test_rejected_partial_stages_report_bound_violations_without_changing_members():
    pool=NativeDistributedMemberPool({0:_member(0),1:_member(1),2:_member(2)})
    observer=NativeAnalysisDiagnostics(pool.layout,3,_topology(),field_lower_bounds={'state':30.})
    members=np.tile([29.,30.,31.],(5,1));original=members.copy()
    observer.observe('controlled_analysis',slice(0,5),members)
    report=observer.finish()
    item=report['stages']['controlled_analysis']['state']
    assert item['below_lower_bound_entries']==5 and item['member_min']==29.
    assert report['stages']['balanced_analysis']['state']['rows']==0
    np.testing.assert_array_equal(members,original)


def test_boundary_metrics_are_same_time_and_do_not_change_members():
    observed, plain = run_cycle(True), run_cycle(False)
    for member in observed.local_analysis.members:
        np.testing.assert_array_equal(observed.local_analysis.members[member], plain.local_analysis.members[member])
    stages = observed.diagnostics["stages"]
    forecast = np.column_stack(list(observed.local_forecast.members.values()))
    raw = forecast @ observed.transform
    expected = relax_analysis_spread(raw, forecast, np.arange(5), {"state": (0, 5)}, {"analysis_spread_relaxation": .5})
    for name, values in [("forecast", forecast), ("raw_analysis", raw),
                         ("controlled_analysis", expected), ("balanced_analysis", expected+7)]:
        report = stages[name]["state"]
        assert report["rows"] == 5
        np.testing.assert_allclose(report["rms_member_sd"], np.sqrt(np.var(values, axis=1, ddof=1).mean()))
        np.testing.assert_allclose(report["mean_rmse"], np.sqrt(np.mean((values.mean(axis=1)-np.arange(5))**2)))


def test_reference_snapshot_respects_noncontiguous_owned_row_order(tmp_path):
    layout = DistributedBlockStateLayout((DistributedStateBlockLayout("h", 0, 4, 1, 3),
                                          DistributedStateBlockLayout("melt", 4, 4, 0, 1)))
    truth = np.arange(24.).reshape(8, 3)
    path = tmp_path/"truth.h5"
    with h5py.File(path, "w") as handle:handle["true_state"] = truth
    observer = prepare_analysis_diagnostics(layout, 3, _topology(), truth_file=path,
                                            reference_rows=np.array([3, 1, 6]), column=2)
    np.testing.assert_array_equal(observer.truth, truth[[3, 1, 6], 2])
    with pytest.raises(RuntimeError, match="reference rows"):
        prepare_analysis_diagnostics(layout, 3, _topology(), truth_file=path,
                                      reference_rows=np.array([3, 3, 6]), column=2)


def test_zero_owned_fields_and_no_truth_are_reported_without_nan():
    layout = DistributedBlockStateLayout((DistributedStateBlockLayout("empty", 0, 3, 0, 0),
                                          DistributedStateBlockLayout("h", 3, 1, 0, 1)))
    observer = NativeAnalysisDiagnostics(layout, 3, _topology())
    for stage in observer.stats:observer.observe(stage, slice(0, 1), np.array([[1., 2., 3.]]))
    result = observer.finish()
    assert result["stages"]["forecast"]["empty"]["rows"] == 0
    assert result["stages"]["forecast"]["h"]["mean_rmse"] is None
    json.dumps(result, allow_nan=False)


def test_deferred_observer_error_is_collectively_reported():
    layout = DistributedBlockStateLayout((DistributedStateBlockLayout("h", 0, 1, 0, 1),))
    observer = NativeAnalysisDiagnostics(layout, 3, _topology())
    observer.observe("forecast", slice(0, 1), np.zeros((1, 2)))
    with pytest.raises(RuntimeError, match="member block"):
        observer.finish()


def test_atomic_report_replay_and_foreign_run_refusal(tmp_path):
    topology = _topology();report = run_cycle(True).diagnostics
    publish_analysis_diagnostics(tmp_path, report, topology=topology, run_id="run", timestep=4, model_time=.5)
    path = tmp_path/"analysis_00000004.json"
    before = path.read_bytes()
    publish_analysis_diagnostics(tmp_path, report, topology=topology, run_id="run", timestep=4, model_time=.5)
    assert path.read_bytes() == before
    with pytest.raises(RuntimeError, match="different run"):
        publish_analysis_diagnostics(tmp_path, report, topology=topology, run_id="other", timestep=4, model_time=.5)
    assert path.read_bytes() == before
