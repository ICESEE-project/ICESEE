# ==============================================================================
# @des: Generic end-of-run performance summary for every ICESEE run.
#
#       Extends (does not replace) the historical "[ICESEE] Performance
#       Metrics" table printed by src/utils/tools.py::display_timing_verbose.
#       Each execution driver passes its own per-rank phase timers (the same
#       local values it reduces for that table) to ``emit_performance_report``,
#       which reduces them over the given communicator (min/mean/max/sum),
#       adds peak RSS, registered I/O counters and registered run metadata,
#       prints a human-readable summary, and writes ``performance.json``. The
#       text and the JSON are rendered from one summary dict.
#
#       Model and execution-mode agnostic by construction: this module never
#       imports an application package and never branches on a model name,
#       execution mode, or field name. Anything backend- or model-specific
#       (resource plans, member stores, application counters) is *registered*
#       by the component that owns it via ``register_run_metadata``,
#       ``register_metrics`` or ``register_io_provider``; the reporter only
#       renders what it is given.
#
#       Timer semantics: phase timers may overlap and may not cover the whole
#       run, so no percentages or additive breakdown of wall time are
#       reported. "Wall-clock" is the maximum elapsed time over ranks;
#       "computational" is the sum of elapsed time over ranks.
#
#       Measurement semantics: a phase a driver does not instrument is either
#       omitted from ``phases`` or passed as None, and is reported as "not
#       measured", never as 0. A phase whose operation count is registered as
#       0 is reported as "no events". Only a timed phase shows a number.
# ==============================================================================
from __future__ import annotations

import datetime
import json
import os
import platform
import socket
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Any, Callable, Mapping, Optional

SCHEMA = "icesee.performance/1"

# The historical rows of display_timing_verbose, in their original order.
# ``key`` is the per-rank phase name drivers pass in; the legacy value of a
# phase row is its max over ranks, exactly as every driver computes it.
LEGACY_PHASES = (
    ("true_wrong_state", "True/Wrong State Time"),
    ("ensemble_init", "Ensemble Init Time"),
    ("forecast_step", "Forecast Step Time"),
    ("analysis_step", "Analysis Step Time"),
    ("init_file_io", "Init file I/O Time"),
    ("forecast_file_io", "Forecast File I/O Time"),
    ("analysis_file_io", "Analysis File I/O Time"),
    ("forecast_noise", "Forecast Noise Time"),
    ("init_ensemble_mean", "Init Ensemble Mean Computation"),
    ("forecast_ensemble_mean", "Forecast Ensemble Mean Computation"),
    ("analysis_ensemble_mean", "Analysis Ensemble Mean Computation"),
)
LEGACY_LABELS = (
    "Computational Time (Σ)",
    "Wall-Clock Time (max)",
    "True/Wrong State Time",
    "Ensemble Init Time",
    "Forecast Step Time",
    "Analysis Step Time",
    "Assimilation Time",
    "Init file I/O Time",
    "Forecast File I/O Time",
    "Analysis File I/O Time",
    "Total File I/O Time",
    "Forecast Noise Time",
    "Init Ensemble Mean Computation",
    "Forecast Ensemble Mean Computation",
    "Analysis Ensemble Mean Computation",
)

_IO_KEYS = ("bytes_read", "bytes_written", "read_time_s", "write_time_s", "reads", "writes")

_RUN_METADATA: "OrderedDict[str, Any]" = OrderedDict()
_METRICS: "OrderedDict[str, OrderedDict[str, Any]]" = OrderedDict()
_IO_PROVIDERS: "OrderedDict[str, Callable[[], Mapping[str, float]]]" = OrderedDict()
_EXTRA_VERSION_PACKAGES: "list[str]" = []
_RECORDED_PHASES: "OrderedDict[str, list]" = OrderedDict()


# ------------------------------------------------------------------------------
# Registration interface (optional; used by execution backends/applications)
# ------------------------------------------------------------------------------
def register_run_metadata(**items: Any) -> None:
    """Record run-describing values (e.g. ensemble size, ranks per model).
    Later registrations of the same key overwrite earlier ones."""
    _RUN_METADATA.update(items)


def register_metrics(section: str, values: Mapping[str, Any]) -> None:
    """Record a named group of additional numeric/text metrics, rendered
    verbatim under ``section``."""
    _METRICS.setdefault(section, OrderedDict()).update(values)


def register_io_provider(name: str, provider: Callable[[], Mapping[str, float]]) -> None:
    """Register a callable returning cumulative per-rank I/O counters, using
    any of the keys ``bytes_read``, ``bytes_written``, ``read_time_s``,
    ``write_time_s``, ``reads``, ``writes``. Evaluated once, at report time."""
    _IO_PROVIDERS[name] = provider


def register_package_versions(*names: str) -> None:
    """Add installed packages whose versions the reproducibility metadata
    should record, beyond ICESEE's own core dependencies."""
    for name in names:
        if name not in _EXTRA_VERSION_PACKAGES:
            _EXTRA_VERSION_PACKAGES.append(name)


def record_phase(name: str, seconds: float, operations: int = 1) -> None:
    """Accumulate this rank's time for a phase measured inside the pipeline
    (outside the driver that emits the report). ``operations=0`` records
    that the phase was skipped, which is reported as "no events"."""
    entry = _RECORDED_PHASES.setdefault(name, [0.0, 0])
    entry[0] += float(seconds)
    entry[1] += int(operations)


def clear_registry() -> None:
    _RUN_METADATA.clear()
    _METRICS.clear()
    _IO_PROVIDERS.clear()
    _EXTRA_VERSION_PACKAGES.clear()
    _RECORDED_PHASES.clear()


# ------------------------------------------------------------------------------
# Local (per-rank) measurements
# ------------------------------------------------------------------------------
def peak_rss_bytes() -> Optional[int]:
    """Peak resident set size of this process, or None where the standard
    library cannot report it (e.g. Windows)."""
    try:
        import resource
    except ImportError:
        return None
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # ru_maxrss is in bytes on macOS and in kilobytes on Linux/BSD.
    return int(peak) if sys.platform == "darwin" else int(peak) * 1024


def _local_io_counters() -> "OrderedDict[str, dict]":
    counters: "OrderedDict[str, dict]" = OrderedDict()
    for name, provider in _IO_PROVIDERS.items():
        try:
            raw = provider() or {}
        except Exception:  # a broken optional provider must not break the run
            continue
        # Keep only what the provider reports: a direction it does not
        # measure is "not measured", not zero bytes.
        counters[name] = {key: float(raw[key]) for key in _IO_KEYS if key in raw}
    return counters


# ------------------------------------------------------------------------------
# Aggregation
# ------------------------------------------------------------------------------
def _stats(values):
    values = [float(v) for v in values]
    total = sum(values)
    mean = total / len(values)
    return {
        "min": min(values),
        "mean": mean,
        "max": max(values),
        "sum": total,
        "imbalance": (max(values) / mean) if mean > 0 else None,
    }


def aggregate_rank_records(records, metadata: Optional[Mapping[str, Any]] = None) -> dict:
    """Build the summary dict from one record per rank.

    Each record: {"elapsed_s": float, "phases": {name: seconds},
    "counts": {name: int}, "peak_rss_bytes": int|None,
    "io": {provider: {io counters}}, "host": str}. Pure function (no MPI), so
    it is unit-testable with synthetic records.
    """
    records = list(records)
    world_size = len(records)

    phase_names = []
    for record in records:
        for name in record.get("phases", {}):
            if name not in phase_names:
                phase_names.append(name)
    phases = OrderedDict()
    for name in phase_names:
        values = [record.get("phases", {}).get(name) for record in records]
        counts = [record.get("counts", {}).get(name) for record in records]
        if all(value is None for value in values):
            continue  # not instrumented on any rank: reported as "not measured"
        entry = _stats(0.0 if value is None else value for value in values)
        # Counts come from the ranks that report one, so a phase run on a
        # subset of ranks (e.g. one root rank) keeps its operation count.
        reported = [(value, count) for value, count in zip(values, counts)
                    if isinstance(count, (int, float))]
        if reported:
            entry["operations"] = int(max(count for _, count in reported))
            if all(count > 0 for _, count in reported):
                per_op = [(value or 0.0) / count for value, count in reported]
                entry["mean_s_per_operation"] = sum(per_op) / len(per_op)
        # "parent/child" names a subset of the phase "parent" (its time is
        # already included there).
        if "/" in name:
            entry["subset_of"] = name.split("/", 1)[0]
        phases[name] = entry

    elapsed = _stats(record.get("elapsed_s", 0.0) for record in records)

    def _measured_sum(entries):
        entries = [entry for entry in entries if entry is not None]
        return sum(entry["max"] for entry in entries) if entries else None

    legacy = OrderedDict()
    legacy["Computational Time (Σ)"] = elapsed["sum"]
    legacy["Wall-Clock Time (max)"] = elapsed["max"]
    phase_max = {key: phases[key]["max"] if key in phases else None for key, _ in LEGACY_PHASES}
    for key, label in LEGACY_PHASES[:4]:
        legacy[label] = phase_max[key]
    legacy["Assimilation Time"] = _measured_sum(
        phases.get(key) for key in ("ensemble_init", "forecast_step", "analysis_step")
    )
    for key, label in LEGACY_PHASES[4:7]:
        legacy[label] = phase_max[key]
    # Every phase named "*_file_io" counts toward total file I/O, so a driver
    # with additional file-I/O phases reports them without changes here.
    legacy["Total File I/O Time"] = _measured_sum(
        entry for name, entry in phases.items() if name.endswith("_file_io")
    )
    for key, label in LEGACY_PHASES[7:]:
        legacy[label] = phase_max[key]
    # A phase with zero registered operations was not exercised in this run.
    no_events = [
        label for key, label in LEGACY_PHASES
        if phases.get(key, {}).get("operations") == 0
    ]

    rss = [record.get("peak_rss_bytes") for record in records]
    memory = None
    if all(isinstance(v, (int, float)) for v in rss):
        memory = {"peak_rss_bytes": _stats(rss)}

    io = OrderedDict()
    for record in records:
        for name, counters in record.get("io", {}).items():
            io.setdefault(name, [])
            io[name].append(counters)
    io_summary = OrderedDict()
    for name, per_rank in io.items():
        entry = OrderedDict()
        for direction, bytes_key, time_key, ops_key, rate_key in (
            ("read", "bytes_read", "read_time_s", "reads", "read_bytes_per_s"),
            ("write", "bytes_written", "write_time_s", "writes", "write_bytes_per_s"),
        ):
            reported = any(key in c for c in per_rank for key in (bytes_key, time_key, ops_key))
            if not reported:
                continue  # this provider does not measure this direction
            total_bytes = sum(c.get(bytes_key, 0.0) for c in per_rank)
            time_max = max(c.get(time_key, 0.0) for c in per_rank)
            entry[bytes_key] = total_bytes
            entry[ops_key] = sum(c.get(ops_key, 0.0) for c in per_rank)
            entry[time_key] = _stats(c.get(time_key, 0.0) for c in per_rank)
            # Aggregate bytes over the slowest rank's time: effective
            # throughput of the whole job for this category.
            entry[rate_key] = (total_bytes / time_max) if time_max > 0 else None
        io_summary[name] = entry

    hosts = sorted({record.get("host", "") for record in records if record.get("host")})

    return OrderedDict(
        schema=SCHEMA,
        metadata=dict(metadata or {}),
        run=OrderedDict(_RUN_METADATA),
        ranks=OrderedDict(world_size=world_size, hosts=hosts, host_count=len(hosts)),
        time=OrderedDict(elapsed_s=elapsed, phases=phases),
        legacy=legacy,
        legacy_no_events=no_events,
        memory=memory,
        io=io_summary,
        metrics=OrderedDict((k, OrderedDict(v)) for k, v in _METRICS.items()),
        notes=[
            "Wall-clock is the maximum elapsed time over ranks; computational is "
            "the sum of elapsed time over ranks.",
            "Phase timers may overlap and need not cover the whole run; they are "
            "not an additive breakdown of wall time.",
        ],
    )


# ------------------------------------------------------------------------------
# Reproducibility metadata (root rank only; cheap, never raises)
# ------------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parents[2]
_VERSION_PACKAGES = ("numpy", "scipy", "mpi4py", "h5py")
_SCHEDULER_ENV = (
    "SLURM_JOB_ID", "SLURM_JOB_NAME", "SLURM_CLUSTER_NAME", "SLURM_JOB_PARTITION",
    "SLURM_JOB_NUM_NODES", "SLURM_NTASKS", "SLURM_CPUS_PER_TASK", "SLURM_JOB_NODELIST",
)


def _git(*args):
    try:
        result = subprocess.run(
            ["git", "-C", str(_REPO_ROOT), *args],
            capture_output=True, text=True, timeout=5,
        )
    except Exception:
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def collect_reproducibility_metadata(extra: Optional[Mapping[str, Any]] = None) -> dict:
    from importlib import metadata as importlib_metadata

    versions = {"python": platform.python_version()}
    for package in (*_VERSION_PACKAGES, *_EXTRA_VERSION_PACKAGES):
        try:
            versions[package] = importlib_metadata.version(package)
        except Exception:
            pass
    if "h5py" in sys.modules:
        try:
            versions["hdf5"] = sys.modules["h5py"].version.hdf5_version
        except Exception:
            pass
    try:
        from mpi4py import MPI
        versions["mpi_library"] = MPI.Get_library_version().splitlines()[0].strip()
    except Exception:
        pass

    revision = _git("rev-parse", "HEAD")
    status = _git("status", "--porcelain", "--untracked-files=no")
    info = OrderedDict(
        timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        hostname=socket.gethostname(),
        platform=platform.platform(),
        git_revision=revision,
        git_dirty=(bool(status) if status is not None else None),
        versions=versions,
    )
    scheduler = {name: os.environ[name] for name in _SCHEDULER_ENV if name in os.environ}
    if scheduler:
        info["scheduler"] = scheduler
    if extra:
        info.update(extra)
    return info


# ------------------------------------------------------------------------------
# Rendering (text and JSON come from the same summary dict)
# ------------------------------------------------------------------------------
def _fmt_s(value):
    return "-" if value is None else f"{value:12.3f}"


def _fmt_bytes(value):
    if value is None:
        return "-"
    for unit, scale in (("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if abs(value) >= scale:
            return f"{value / scale:.2f} {unit}"
    return f"{value:.0f} B"


def _fmt_duration(seconds):
    """DAY:HR:MIN:SEC.ms, the historical display_timing_verbose format."""
    days = int(seconds // 86400)
    hours = int((seconds % 86400) // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    millis = int((seconds % 1) * 1000)
    return f"{days:02d}:{hours:02d}:{minutes:02d}:{secs:02d}.{millis:03d}"


def render_performance_report(summary: Mapping[str, Any]) -> str:
    width = 78
    meta = summary.get("metadata", {})
    ranks = summary.get("ranks", {})

    # The historical "[ICESEE] Performance Metrics" table, first and in its
    # original row order, so runs stay comparable with earlier output.
    header = f"[ICESEE] Performance Metrics ({ranks.get('world_size')} ranks)"
    lines = ["=" * width, f"{header:<58}(DAY:HR:MIN:SEC.ms)", "-" * width]
    no_events = set(summary.get("legacy_no_events", ()))
    for label, value in summary.get("legacy", {}).items():
        if label in no_events:
            shown = "no events"
        elif value is None:
            shown = "not measured"
        else:
            shown = _fmt_duration(value)
        lines.append(f"  {label:<44}{shown:>20}")
    lines.append("-" * width)

    lines.append("Run")
    for label, value in (
        ("Timestamp (UTC)", meta.get("timestamp_utc")),
        ("Hosts", ", ".join(ranks.get("hosts", [])[:4])
         + (" ..." if ranks.get("host_count", 0) > 4 else "")),
        ("Git revision", (meta.get("git_revision") or "-")[:12]
         + (" (dirty)" if meta.get("git_dirty") else "")),
        ("MPI ranks", ranks.get("world_size")),
    ):
        lines.append(f"  {label:<34}{value}")
    for key, value in summary.get("run", {}).items():
        lines.append(f"  {key:<34}{value}")
    versions = dict(meta.get("versions", {}))
    mpi_library = versions.pop("mpi_library", None)
    if versions:
        lines.append(f"  {'Versions':<34}" + ", ".join(f"{k} {v}" for k, v in versions.items()))
    if mpi_library:
        lines.append(f"  {'MPI library':<34}{mpi_library}")
    for key, value in meta.get("scheduler", {}).items():
        lines.append(f"  {key:<34}{value}")

    time_block = summary.get("time", {})
    elapsed = time_block.get("elapsed_s", {})
    lines.append("")
    lines.append(f"Timing [s]{'min':>36}{'mean':>12}{'max':>12}{'imbal.':>8}")
    lines.append(
        f"  {'Elapsed (per rank)':<34}{_fmt_s(elapsed.get('min'))}"
        f"{_fmt_s(elapsed.get('mean'))}{_fmt_s(elapsed.get('max'))}"
    )
    def _phase_line(label, entry):
        if entry.get("operations") == 0:
            return f"  {label:<34}{'no events':>36}"
        imbalance = entry.get("imbalance")
        line = (
            f"  {label:<34}{_fmt_s(entry['min'])}{_fmt_s(entry['mean'])}{_fmt_s(entry['max'])}"
            f"{'' if imbalance is None else f'{imbalance:8.2f}'}"
        )
        if "mean_s_per_operation" in entry:
            line += f"   ({entry['operations']} ops, {entry['mean_s_per_operation']:.4f} s/op)"
        elif "operations" in entry:
            line += f"   ({entry['operations']} ops)"
        return line

    phase_entries = time_block.get("phases", {})
    for name, entry in phase_entries.items():
        if entry.get("subset_of") in phase_entries:
            continue  # rendered under its parent
        lines.append(_phase_line(name, entry))
        for child, child_entry in phase_entries.items():
            if child_entry.get("subset_of") == name:
                lines.append(_phase_line("  of which " + child.split("/", 1)[1], child_entry))
    lines.append("  " + "-" * (width - 2))
    lines.append(f"  {'Wall-clock time (max over ranks)':<34}{_fmt_s(elapsed.get('max'))}")
    lines.append(f"  {'Computational time (sum over ranks)':<34}{_fmt_s(elapsed.get('sum'))}")
    lines.append("  Phase rows are accumulated time per phase; phases may overlap and do")
    lines.append("  not add up to the wall-clock time.")

    memory = summary.get("memory")
    if memory:
        rss = memory["peak_rss_bytes"]
        lines.append("")
        lines.append(f"Memory{'min':>40}{'mean':>12}{'max':>12}")
        lines.append(
            f"  {'Peak RSS per rank [GB]':<34}{rss['min'] / 1e9:12.3f}"
            f"{rss['mean'] / 1e9:12.3f}{rss['max'] / 1e9:12.3f}"
        )
        lines.append(f"  {'Sum of per-rank peak RSS [GB]':<34}{rss['sum'] / 1e9:12.3f}")

    io = summary.get("io", {})
    if io:
        lines.append("")
        lines.append("I/O (bytes summed over ranks; rate = bytes / slowest rank's time)")
        for name, entry in io.items():
            for direction, bytes_key, time_key, rate_key in (
                ("read", "bytes_read", "read_time_s", "read_bytes_per_s"),
                ("write", "bytes_written", "write_time_s", "write_bytes_per_s"),
            ):
                if bytes_key not in entry:
                    continue
                rate = entry.get(rate_key)
                lines.append(
                    f"  {name + ' ' + direction:<34}{_fmt_bytes(entry[bytes_key]):>12}"
                    f"{entry[time_key]['max']:12.3f} s"
                    f"{'' if rate is None else '   ' + _fmt_bytes(rate) + '/s'}"
                )

    for section, values in summary.get("metrics", {}).items():
        lines.append("")
        lines.append(section)
        for key, value in values.items():
            lines.append(f"  {key:<34}{value}")

    lines.append("=" * width)
    return "\n".join(lines)


def write_performance_json(summary: Mapping[str, Any], path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as handle:
        json.dump(summary, handle, indent=2, default=str)
    os.replace(tmp, path)
    return path


# ------------------------------------------------------------------------------
# Driver entry point (collective over ``comm``)
# ------------------------------------------------------------------------------
def emit_performance_report(
    comm,
    *,
    elapsed_s: float,
    phases: Mapping[str, float],
    counts: Optional[Mapping[str, int]] = None,
    output_dir=None,
    emit: Optional[Callable[[str], None]] = None,
) -> Optional[dict]:
    """Collective over ``comm``: every rank must call it once, with its own
    local timers. Rank 0 prints the summary and writes
    ``<output_dir>/performance.json``; returns the summary on rank 0 and None
    elsewhere. Adds no barrier beyond the single gather."""
    # Phases recorded inside the pipeline (record_phase) complete the
    # driver's own timers and are listed first, in the order they were
    # recorded; a phase the driver passes explicitly wins.
    counts = dict(counts or {})
    merged = {}
    for name, (seconds, operations) in _RECORDED_PHASES.items():
        if name not in phases:
            merged[name] = seconds
            counts.setdefault(name, operations)
    merged.update(phases)
    phases = merged
    record = {
        "elapsed_s": float(elapsed_s),
        "phases": {
            name: None if value is None else float(value) for name, value in phases.items()
        },
        "counts": counts,
        "peak_rss_bytes": peak_rss_bytes(),
        "io": _local_io_counters(),
        "host": socket.gethostname(),
    }
    records = comm.gather(record, root=0)
    if comm.Get_rank() != 0:
        return None

    summary = aggregate_rank_records(records, metadata=collect_reproducibility_metadata())
    text = render_performance_report(summary)
    if emit is None:
        print(text, flush=True)
    else:
        emit(text)
    if output_dir is not None:
        try:
            write_performance_json(summary, Path(output_dir) / "performance.json")
        except OSError as error:
            print(f"[ICESEE] could not write performance.json: {error}", flush=True)
    return summary
