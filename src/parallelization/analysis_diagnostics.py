"""Small per-field diagnostics at native analysis boundaries."""
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np

from .distributed_analysis import iter_local_ensemble_row_blocks

STAGES = ("forecast", "raw_analysis", "controlled_analysis", "balanced_analysis")


class NativeAnalysisDiagnostics:
    """Accumulate nodal moments from bounded rows, never a global ensemble.

    Only ensemble-slot zero contributes rows. Spatial ranks exchange small
    per-field summaries once at finish; ghost rows are absent from the layout.
    Optional truth is an owned-row snapshot at the same physical model time.
    It is used only for reporting and never passed to the analysis update.
    """

    def __init__(self, layout, number_of_members, topology, truth_owned=None, field_lower_bounds=None):
        self.layout, self.topology = layout, topology
        self.lower_bounds = dict(field_lower_bounds or {})
        if any(not np.isfinite(float(value)) for value in self.lower_bounds.values()):
            raise ValueError('diagnostic lower bounds must be finite')
        self.number_of_members = int(number_of_members)
        if self.number_of_members < 2:
            raise ValueError("analysis diagnostics require at least two members")
        if hasattr(layout, "iter_owned_global_intervals"):
            intervals = list(layout.iter_owned_global_intervals())
            self.fields = {name: rows for name, rows, _ in intervals}
            self.global_intervals = {name: list(global_rows) for name, _, global_rows in intervals}
        else:
            self.fields = {"state": slice(0, layout.owned_size)}
            self.global_intervals = {"state": [layout.owned_start, layout.owned_stop]}
        self.truth = None if truth_owned is None else np.asarray(truth_owned, dtype=float)
        if self.truth is not None and (self.truth.shape != (layout.owned_size,)
                                      or not np.isfinite(self.truth).all()):
            raise ValueError("diagnostic truth must be finite and match owned rows")
        self.errors = []
        self.hashers = {stage: {name: hashlib.sha256() for name in self.fields} for stage in STAGES}
        self.stats = {stage: {name: dict(rows=0, variance=0., error2=0.,
                          member_min=None, member_max=None, nonfinite=0,
                          owned_digest=self.hashers[stage][name].hexdigest(),
                          owned_global_interval=self.global_intervals[name])
                             for name in self.fields} for stage in STAGES}

    def observe(self, stage, rows, members):
        """Local-only callback; defer errors until every rank finishes its blocks."""
        if not self.topology.is_ensemble_root:
            return
        try:
            if stage not in self.stats or members.shape != (rows.stop-rows.start, self.number_of_members):
                raise ValueError("invalid diagnostic stage or member block shape")
            for name, owned in self.fields.items():
                begin, end = max(rows.start, owned.start), min(rows.stop, owned.stop)
                if end <= begin:
                    continue
                block = np.asarray(members[begin-rows.start:end-rows.start], dtype=float)
                item = self.stats[stage][name]
                self.hashers[stage][name].update(np.ascontiguousarray(block, dtype='<f8').tobytes())
                item["owned_digest"] = self.hashers[stage][name].hexdigest()
                valid = np.isfinite(block).all(axis=1)
                item["nonfinite"] += int((~valid).sum())
                block = block[valid]
                if name in self.lower_bounds:
                    below = block < float(self.lower_bounds[name])
                    item['below_lower_bound_entries'] = item.get('below_lower_bound_entries', 0)+int(below.sum())
                item["rows"] += int(len(block))
                if not len(block):
                    continue
                item["variance"] += float(np.var(block, axis=1, ddof=1).sum())
                if self.truth is not None:
                    item["error2"] += float(np.sum((block.mean(axis=1)-self.truth[begin:end][valid])**2))
                low, high = float(block.min()), float(block.max())
                item["member_min"] = low if item["member_min"] is None else min(item["member_min"], low)
                item["member_max"] = high if item["member_max"] is None else max(item["member_max"], high)
        except Exception as exc:
            self.errors.append(f"{type(exc).__name__}: {exc}")

    def observe_snapshot(self, stage, snapshot, *, row_chunk_size):
        for rows, members in iter_local_ensemble_row_blocks(
                snapshot.members, self.number_of_members, self.topology.ensemble_comm,
                row_chunk_size=row_chunk_size):
            self.observe(stage, rows, members)

    def finish(self):
        report, error = None, None
        if self.topology.is_ensemble_root:
            pieces = self.topology.spatial_comm.allgather(
                (self.stats, self.truth is not None, self.errors))
            errors = [message for _, _, messages in pieces for message in messages]
            if errors:
                error = "; ".join(errors)
            elif len({present for _, present, _ in pieces}) != 1:
                error = "diagnostic truth availability differs across spatial ranks"
            else:
                report = {"weighting": "uniform_owned_rows_not_area_weighted",
                          "number_of_members": self.number_of_members,
                          "truth_available": pieces[0][1], "stages": {}}
                for stage in STAGES:
                    report["stages"][stage] = {}
                    for name in self.fields:
                        parts = [stats[stage][name] for stats, _, _ in pieces]
                        count = sum(part["rows"] for part in parts)
                        variance = sum(part["variance"] for part in parts)
                        error2 = sum(part["error2"] for part in parts)
                        lows = [part["member_min"] for part in parts if part["member_min"] is not None]
                        highs = [part["member_max"] for part in parts if part["member_max"] is not None]
                        report["stages"][stage][name] = dict(
                            rows=count, rms_member_sd=float(np.sqrt(variance/count)) if count else None,
                            mean_rmse=float(np.sqrt(error2/count)) if count and pieces[0][1] else None,
                            member_min=min(lows) if lows else None, member_max=max(highs) if highs else None,
                            nonfinite_rows=sum(part["nonfinite"] for part in parts),
                            lower_bound=self.lower_bounds.get(name),
                            below_lower_bound_entries=(sum(part.get('below_lower_bound_entries', 0) for part in parts)
                                                       if name in self.lower_bounds else None),
                            owned_digests=[part["owned_digest"] for part in parts],
                            owned_global_intervals=[part["owned_global_interval"] for part in parts])
        report, error = self.topology.ensemble_comm.bcast((report, error), root=0)
        if error:
            raise RuntimeError(f"analysis diagnostics failed: {error}")
        return report


def prepare_analysis_diagnostics(layout, number_of_members, topology, *, truth_file=None,
                                 reference_rows=None, column=None, field_lower_bounds=None):
    """Read one reference snapshot on ensemble roots; broadcast setup errors."""
    truth, diagnostic, error = None, None, None
    try:
        if truth_file is not None and topology.is_ensemble_root:
            import h5py
            rows = np.asarray(reference_rows)
            if (rows.shape != (layout.owned_size,) or rows.dtype.kind not in "iu"
                    or len(np.unique(rows)) != len(rows)):
                raise ValueError("diagnostic reference rows must uniquely match owned state")
            with h5py.File(truth_file, "r") as handle:
                dataset = handle["true_state"]
                if column is None or not 0 <= column < dataset.shape[1]:
                    raise ValueError("diagnostic reference time column is unavailable")
                if rows.size and (rows.min() < 0 or rows.max() >= dataset.shape[0]):
                    raise ValueError("diagnostic reference rows leave the truth state")
                truth = np.empty(len(rows))
                order = np.argsort(rows)
                for begin in range(0, len(order), 4096):
                    indices = order[begin:begin+4096]
                    selected = rows[indices]
                    selection = (slice(int(selected[0]), int(selected[-1])+1)
                                 if np.all(np.diff(selected) == 1) else selected)
                    truth[indices] = dataset[selection, int(column)]
        diagnostic = NativeAnalysisDiagnostics(layout, number_of_members, topology, truth, field_lower_bounds)
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
    errors = [value for value in topology.world.allgather(error) if value]
    if errors:
        raise RuntimeError("analysis diagnostic setup failed: " + "; ".join(errors))
    return diagnostic


def publish_analysis_diagnostics(directory, report, *, topology, run_id, timestep, model_time):
    """Atomically replace a cycle's report; an interrupted cycle may be replayed."""
    error = None
    if topology.world_rank == 0:
        temporary = None
        try:
            directory = Path(directory)
            directory.mkdir(parents=True, exist_ok=True)
            target = directory / f"analysis_{int(timestep):08d}.json"
            if target.exists() and json.loads(target.read_text()).get("run_id") != run_id:
                raise ValueError("existing diagnostic report belongs to a different run")
            record = dict(report, run_id=run_id, timestep=int(timestep), model_time=float(model_time),
                          checkpoint_commit_required=True)
            encoded = json.dumps(record, allow_nan=False, indent=2) + "\n"
            fd, temporary = tempfile.mkstemp(prefix=".analysis-", dir=directory)
            with os.fdopen(fd, "w") as handle:
                handle.write(encoded)
            os.replace(temporary, target)
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
        finally:
            if temporary is not None and os.path.exists(temporary):
                os.unlink(temporary)
    error = topology.world.bcast(error, root=0)
    if error:
        raise RuntimeError("analysis diagnostic publication failed: " + error)
