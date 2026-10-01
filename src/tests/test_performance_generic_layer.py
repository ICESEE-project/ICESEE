# ==============================================================================
# @des: Tests for the generic end-of-run performance summary
#       (src/utils/performance.py): model/mode independence, the optional
#       registration interface, historical metric compatibility, and that the
#       text and JSON outputs come from the same summary. All inputs are
#       synthetic per-rank records; nothing depends on wall-clock speed.
# ==============================================================================
from __future__ import annotations

import ast
import json
import re
from pathlib import Path

import pytest

from ICESEE.src.utils import performance
from ICESEE.src.utils.performance import (
    LEGACY_LABELS,
    aggregate_rank_records,
    emit_performance_report,
    register_io_provider,
    register_metrics,
    register_run_metadata,
    render_performance_report,
    write_performance_json,
)

_REPO_ROOT = Path(__file__).resolve().parents[2]
_GENERIC_MODULE = _REPO_ROOT / "src" / "utils" / "performance.py"


@pytest.fixture(autouse=True)
def _clean_registry():
    performance.clear_registry()
    yield
    performance.clear_registry()


def _record(elapsed, phases, rss=1_000_000, host="node-a", counts=None, io=None):
    return {
        "elapsed_s": elapsed,
        "phases": phases,
        "counts": counts or {},
        "peak_rss_bytes": rss,
        "io": io or {},
        "host": host,
    }


# --- 1/2: no application imports, no model/mode/field knowledge ----------------
def test_generic_module_imports_no_application_package():
    tree = ast.parse(_GENERIC_MODULE.read_text())
    imported = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported += [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
    assert not [name for name in imported if "applications" in name or "_model" in name], imported


def test_generic_module_has_no_model_mode_or_field_knowledge():
    source = _GENERIC_MODULE.read_text().lower()
    code = "\n".join(
        line.split("#", 1)[0] for line in source.splitlines()
    )
    for word in ("icepack", "firedrake", "issm", "lorenz", "flowline", "basal_melt",
                 "pig", "member_store", "ranks_per_model", "resource_plan"):
        assert word not in code, word
    assert not re.search(r"execution_mode\s*==", code)
    assert not re.search(r"model(_name)?\s*==", code)


# --- 3/4: optional registration, and missing optional data --------------------
def test_optional_metrics_are_rendered_without_reporter_changes():
    register_run_metadata(some_backend_parameter=7, another_label="abc")
    register_metrics("Backend counters", {"widgets_processed": 12})
    register_io_provider("scratch_store", lambda: {"bytes_read": 4e9, "read_time_s": 2.0, "reads": 3})

    record = _record(5.0, {"forecast_step": 3.0}, io=performance._local_io_counters())
    summary = aggregate_rank_records([record, dict(record)])
    text = render_performance_report(summary)

    assert summary["run"]["some_backend_parameter"] == 7
    assert "some_backend_parameter" in text and "abc" in text
    assert summary["metrics"]["Backend counters"]["widgets_processed"] == 12
    assert "Backend counters" in text and "widgets_processed" in text
    io = summary["io"]["scratch_store"]
    assert io["bytes_read"] == 8e9
    assert io["read_bytes_per_s"] == pytest.approx(8e9 / 2.0)
    assert "scratch_store read" in text
    # The provider measures only reads: no write line, not "0 B written".
    assert "bytes_written" not in io and "scratch_store write" not in text


def test_missing_optional_information_does_not_break_reporting():
    summary = aggregate_rank_records([{"elapsed_s": 1.0}])
    text = render_performance_report(summary)
    assert summary["memory"] is None
    assert summary["io"] == {}
    assert summary["legacy"]["Forecast Step Time"] is None
    assert "[ICESEE] Performance Metrics (1 ranks)" in text

    unknown_rss = aggregate_rank_records([_record(1.0, {}, rss=None), _record(1.0, {})])
    assert unknown_rss["memory"] is None


def test_failing_io_provider_is_skipped():
    def broken():
        raise RuntimeError("backend unavailable")

    register_io_provider("broken", broken)
    register_io_provider("ok", lambda: {"bytes_written": 10})
    assert list(performance._local_io_counters()) == ["ok"]


# --- 5: same reporter for mode 0/1/2/3-shaped inputs, no branching ------------
@pytest.mark.parametrize(
    "mode,world_size,phases",
    [
        (0, 1, {"forecast_step": 4.0, "analysis_step": 1.0, "init_file_io": 0.1}),
        (1, 3, {"forecast_step": 4.0, "analysis_step": 1.0, "forecast_noise": 0.2}),
        (2, 4, {"forecast_step": 4.0, "setup_file_io": 0.3, "final_file_io": 0.2,
                "forecast_file_io": 0.5}),
        (3, 4, {"forecast_step": 4.0, "true_wrong_state": 2.0, "analysis_file_io": 0.4}),
    ],
)
def test_reporter_accepts_every_execution_mode_shape(mode, world_size, phases):
    register_run_metadata(execution_mode=mode, ensemble_size=4)
    records = [_record(10.0 + rank, phases, host=f"n{rank % 2}") for rank in range(world_size)]
    summary = aggregate_rank_records(records)
    text = render_performance_report(summary)
    assert summary["ranks"]["world_size"] == world_size
    assert summary["run"]["execution_mode"] == mode
    assert set(phases) <= set(summary["time"]["phases"])
    assert f"{'execution_mode':<34}{mode}" in text


# --- rank statistics, legacy metrics -----------------------------------------
def test_rank_statistics_and_wall_vs_computational_time():
    records = [
        _record(10.0, {"forecast_step": 2.0}, rss=1e9),
        _record(12.0, {"forecast_step": 6.0}, rss=3e9),
    ]
    summary = aggregate_rank_records(records)
    forecast = summary["time"]["phases"]["forecast_step"]
    assert (forecast["min"], forecast["mean"], forecast["max"], forecast["sum"]) == (2.0, 4.0, 6.0, 8.0)
    assert forecast["imbalance"] == pytest.approx(1.5)
    assert summary["legacy"]["Wall-Clock Time (max)"] == 12.0
    assert summary["legacy"]["Computational Time (Σ)"] == 22.0
    rss = summary["memory"]["peak_rss_bytes"]
    assert (rss["min"], rss["max"], rss["sum"]) == (1e9, 3e9, 4e9)


def test_operation_counts_give_time_per_operation():
    records = [
        _record(1.0, {"forecast_step": 8.0}, counts={"forecast_step": 4}),
        _record(1.0, {"forecast_step": 4.0}, counts={"forecast_step": 4}),
    ]
    forecast = aggregate_rank_records(records)["time"]["phases"]["forecast_step"]
    assert forecast["operations"] == 4
    assert forecast["mean_s_per_operation"] == pytest.approx(1.5)


# --- 6: every historical metric remains available ----------------------------
def test_all_historical_metrics_are_preserved_with_historical_definitions():
    phases = {
        "true_wrong_state": 1.0, "ensemble_init": 2.0, "forecast_step": 3.0,
        "analysis_step": 4.0, "init_file_io": 0.5, "forecast_file_io": 0.25,
        "analysis_file_io": 0.125, "forecast_noise": 0.1, "init_ensemble_mean": 0.01,
        "forecast_ensemble_mean": 0.02, "analysis_ensemble_mean": 0.03,
    }
    legacy = aggregate_rank_records([_record(20.0, phases)])["legacy"]
    assert tuple(legacy) == LEGACY_LABELS
    assert legacy["Assimilation Time"] == 2.0 + 3.0 + 4.0
    assert legacy["Total File I/O Time"] == 0.5 + 0.25 + 0.125
    assert legacy["Forecast Noise Time"] == 0.1

    # Additional "*_file_io" phases (mode 2 passes setup/final file I/O)
    # count toward total file I/O, matching that driver's own table.
    extra = aggregate_rank_records([_record(1.0, {"init_file_io": 1.0, "final_file_io": 2.0})])
    assert extra["legacy"]["Total File I/O Time"] == 3.0


def test_legacy_labels_match_the_historical_table():
    source = (_REPO_ROOT / "src" / "utils" / "tools.py").read_text()
    table = source.split("def display_timing_verbose", 1)[1].split("def get_grid_dimensions", 1)[0]
    historical = re.findall(r'\("([^"]+)", format_time', table)
    assert tuple(historical) == LEGACY_LABELS


# --- 7: text and JSON from the same summary object ---------------------------
def test_text_and_json_are_rendered_from_the_same_summary(tmp_path):
    register_run_metadata(ensemble_size=4)
    summary = aggregate_rank_records([
        _record(3.0, {"forecast_step": 1.25}),
        _record(4.0, {"forecast_step": 2.5}),
    ], metadata={"git_revision": "abc123def4567"})
    path = write_performance_json(summary, tmp_path / "performance.json")
    loaded = json.loads(path.read_text())
    assert loaded == json.loads(json.dumps(summary, default=str))
    text = render_performance_report(loaded)
    assert text == render_performance_report(summary)
    assert "abc123def456" in text and "ensemble_size" in text


def test_emit_writes_json_and_prints_summary(tmp_path):
    class SingleRankComm:
        def gather(self, value, root=0):
            return [value]

        def Get_rank(self):
            return 0

    printed = []
    summary = emit_performance_report(
        SingleRankComm(), elapsed_s=2.0, phases={"forecast_step": 1.0},
        counts={"forecast_step": 2}, output_dir=tmp_path, emit=printed.append,
    )
    assert len(printed) == 1 and "[ICESEE] Performance Metrics" in printed[0]
    written = json.loads((tmp_path / "performance.json").read_text())
    assert written["schema"] == summary["schema"] == "icesee.performance/1"
    assert written["time"]["phases"]["forecast_step"]["operations"] == 2
    assert "git_revision" in written["metadata"] and "versions" in written["metadata"]


def test_non_root_ranks_return_none_and_write_nothing(tmp_path):
    class NonRootComm:
        def gather(self, value, root=0):
            return None

        def Get_rank(self):
            return 1

    assert emit_performance_report(
        NonRootComm(), elapsed_s=1.0, phases={}, output_dir=tmp_path, emit=lambda _: None
    ) is None
    assert not (tmp_path / "performance.json").exists()


# --- measured vs not measured vs no events -----------------------------------
def test_not_measured_no_events_and_measured_zero_are_distinct():
    summary = aggregate_rank_records([
        _record(
            5.0,
            {"forecast_step": 2.0, "analysis_step": None, "forecast_noise": 0.0,
             "init_file_io": 0.0},
            counts={"forecast_step": 4, "analysis_step": 1},
        ),
    ])
    legacy = summary["legacy"]
    assert legacy["Analysis Step Time"] is None            # instrumented nowhere
    assert legacy["Ensemble Init Time"] is None            # not passed at all
    assert legacy["Forecast Noise Time"] == 0.0            # genuinely measured zero
    assert "analysis_step" not in summary["time"]["phases"]
    # Sums use only measured components instead of treating gaps as zero.
    assert legacy["Assimilation Time"] == 2.0
    assert legacy["Total File I/O Time"] == 0.0

    text = render_performance_report(summary)
    assert re.search(r"Analysis Step Time\s+not measured", text)
    assert re.search(r"Ensemble Init Time\s+not measured", text)
    assert re.search(r"Forecast Noise Time\s+00:00:00:00\.000", text)

    no_events = aggregate_rank_records([
        _record(5.0, {"forecast_step": 2.0, "analysis_step": 0.0},
                counts={"forecast_step": 4, "analysis_step": 0}),
    ])
    assert "Analysis Step Time" in no_events["legacy_no_events"]
    assert re.search(r"Analysis Step Time\s+no events", render_performance_report(no_events))

    nothing_measured = aggregate_rank_records([_record(1.0, {})])["legacy"]
    assert nothing_measured["Assimilation Time"] is None
    assert nothing_measured["Total File I/O Time"] is None


# --- the historical table is part of the single report -----------------------
def test_single_report_contains_the_historical_15_row_table_in_order():
    phases = {key: 1.5 for key, _ in performance.LEGACY_PHASES}
    records = [_record(90061.25, phases), _record(90061.25, phases)]
    text = render_performance_report(aggregate_rank_records(records))

    assert text.count("[ICESEE] Performance Metrics (2 ranks)") == 1
    assert "(DAY:HR:MIN:SEC.ms)" in text
    positions = [text.index(label) for label in LEGACY_LABELS]
    assert positions == sorted(positions)
    assert re.search(r"Wall-Clock Time \(max\)\s+01:01:01:01\.250", text)
    assert re.search(r"Forecast Step Time\s+00:00:00:01\.500", text)


def test_emit_prints_exactly_one_report():
    class SingleRankComm:
        def gather(self, value, root=0):
            return [value]

        def Get_rank(self):
            return 0

    printed = []
    summary = emit_performance_report(
        SingleRankComm(), elapsed_s=1.0,
        phases={"forecast_step": 0.5, "analysis_step": None},
        counts={"analysis_step": 1},
        emit=printed.append,
    )
    assert len(printed) == 1
    assert printed[0].count("[ICESEE] Performance Metrics") == 1
    assert summary["legacy"]["Analysis Step Time"] is None


_RUN_DRIVERS = (
    "src/run_model_da/icesee_da_serial.py",
    "src/run_model_da/icesee_da_partial_parallel.py",
    "src/run_model_da/icesee_da_full_parallel.py",
    "applications/lorenz_model/lorenz_utils/mode3_runner.py",
    "applications/flowline_model/examples/flowline_1d/mode3_runner.py",
    "applications/issm_model/examples/ISMIP_Choi/mode3_runner.py",
    "applications/icepack_model/examples/idealized_pig/mode3_runner.py",
)


@pytest.mark.parametrize("relative_path", _RUN_DRIVERS)
def test_every_run_driver_emits_one_report_through_the_generic_layer(relative_path):
    tree = ast.parse((_REPO_ROOT / relative_path).read_text())
    called = [
        node.func.id if isinstance(node.func, ast.Name) else getattr(node.func, "attr", None)
        for node in ast.walk(tree) if isinstance(node, ast.Call)
    ]
    assert called.count("emit_performance_report") == 1
    assert "display_timing_verbose" not in called
    assert "display_timing_default" not in called
    assert not (_REPO_ROOT / "src" / "utils" / "performance_report.py").exists()


def test_not_measured_values_survive_the_json_round_trip(tmp_path):
    summary = aggregate_rank_records([
        _record(3.0, {"forecast_step": 1.0, "analysis_step": None},
                counts={"analysis_step": 2}),
    ])
    loaded = json.loads(write_performance_json(summary, tmp_path / "p.json").read_text())
    assert loaded["legacy"]["Analysis Step Time"] is None
    assert render_performance_report(loaded) == render_performance_report(summary)


def test_generic_module_has_no_distributed_topology_vocabulary():
    code = "\n".join(
        line.split("#", 1)[0] for line in _GENERIC_MODULE.read_text().lower().splitlines()
    )
    for word in ("spatial_ranks", "ensemble_groups", "topology", "checkpoint",
                 "execution_mode", "petsc", "rounds"):
        assert word not in code, word


# --- phases recorded inside the pipeline --------------------------------------
def test_recorded_phases_merge_into_the_single_report_and_skips_are_no_events():
    class SingleRankComm:
        def gather(self, value, root=0):
            return [value]

        def Get_rank(self):
            return 0

    performance.record_phase("truth_generation", 2.0)
    performance.record_phase("truth_generation", 1.0)
    performance.record_phase("wrong_reference_generation", 0.0, operations=0)
    performance.record_phase("forecast_step", 99.0)  # the driver's own value wins
    printed = []
    summary = emit_performance_report(
        SingleRankComm(), elapsed_s=5.0, phases={"forecast_step": 1.5},
        emit=printed.append,
    )
    phases = summary["time"]["phases"]
    assert phases["truth_generation"]["max"] == 3.0
    assert phases["truth_generation"]["operations"] == 2
    assert phases["wrong_reference_generation"]["operations"] == 0
    assert phases["forecast_step"]["max"] == 1.5
    assert re.search(r"wrong_reference_generation\s+no events", printed[0])
    assert not re.search(r"wrong_reference_generation\s+0\.000", printed[0])


def test_phase_run_on_one_rank_keeps_its_operation_count():
    records = [
        _record(5.0, {"truth_generation": 4.0}, counts={"truth_generation": 1}),
        _record(5.0, {}),
        _record(5.0, {}),
    ]
    truth = aggregate_rank_records(records)["time"]["phases"]["truth_generation"]
    assert truth["operations"] == 1
    assert truth["mean_s_per_operation"] == pytest.approx(4.0)
    assert (truth["min"], truth["max"]) == (0.0, 4.0)


def test_subset_phase_is_marked_and_rendered_under_its_parent(tmp_path):
    records = [
        _record(10.0, {"forecast_step": 8.0, "forecast_step/with_analysis": 2.0},
                counts={"forecast_step": 40, "forecast_step/with_analysis": 3}),
    ]
    summary = aggregate_rank_records(records)
    child = summary["time"]["phases"]["forecast_step/with_analysis"]
    assert child["subset_of"] == "forecast_step"
    assert "subset_of" not in summary["time"]["phases"]["forecast_step"]
    text = render_performance_report(summary)
    lines = text.splitlines()
    parent_at = next(i for i, l in enumerate(lines) if l.lstrip().startswith("forecast_step "))
    assert lines[parent_at + 1].lstrip().startswith("of which with_analysis")
    assert text.count("with_analysis") == 1
    # A subset adds nothing to totals: Assimilation Time is the forecast alone.
    assert summary["legacy"]["Assimilation Time"] == 8.0
    loaded = json.loads(write_performance_json(summary, tmp_path / "p.json").read_text())
    assert render_performance_report(loaded) == text
