"""Coordinated stops at committed cycle boundaries."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import signal
import tempfile
import time


def atomic_text(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix="." + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            stream.write(value)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def atomic_json(path, value):
    atomic_text(path, json.dumps(value, sort_keys=True, indent=2, allow_nan=False))


class CycleStopController:
    """Observe signals/deadlines only after a checkpoint has committed.

    The launcher must reserve enough wall time for the slowest complete cycle
    and checkpoint. A signal cannot interrupt or repair a hung model collective.
    """

    def __init__(self, options, comm, *, elapsed=0.0, clock=time.monotonic):
        self.comm, self.clock = comm, clock
        self.rank = int(comm.Get_rank())
        self.start = clock() - elapsed
        self.limit = options.get("run_walltime_seconds")
        self.cycles = options.get("run_segment_cycles")
        self.enabled = self.limit is not None or self.cycles is not None
        self.reason = None
        self.previous = {}
        if self.limit is not None:
            self.limit = float(self.limit)
            if not math.isfinite(self.limit) or self.limit <= 0:
                raise ValueError("run_walltime_seconds must be finite and positive")
        if self.cycles is not None:
            if isinstance(self.cycles, bool) or int(self.cycles) != self.cycles or self.cycles < 1:
                raise ValueError("run_segment_cycles must be a positive integer")
        if self.enabled and int(options.get("checkpoint_every", 1)) != 1:
            raise ValueError("graceful continuation requires checkpoint_every=1")

    def __enter__(self):
        if self.enabled:
            for sig in (signal.SIGTERM, signal.SIGUSR1):
                self.previous[sig] = signal.getsignal(sig)
                signal.signal(sig, self._signal)
        return self

    def _signal(self, signum, frame):
        self.reason = signal.Signals(signum).name

    def after_commit(self, cycles_completed):
        reason = self.reason
        if self.limit is not None and self.clock() - self.start >= self.limit:
            reason = reason or "walltime_budget"
        if self.cycles is not None and cycles_completed >= self.cycles:
            reason = reason or "cycle_budget"
        # A signal delivered to any task participates in the common decision.
        reasons = self.comm.allgather(reason)
        self.reason = next((r for r in reasons if r), None)
        return self.reason is not None

    def __exit__(self, *exc):
        for sig, previous in self.previous.items():
            signal.signal(sig, previous)
        self.previous.clear()
