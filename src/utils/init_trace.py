# ==============================================================================
# @des: Opt-in, rank-labelled initialization tracing and hang watchdog.
#
#   ICESEE_DEBUG_INIT=1      print one flushed line per ENTER/EXIT of the
#                            traced initialization/dispatch stages, labelled
#                            with world rank/size, host and any caller fields
#                            (group, model rank, communicator size, ...).
#   ICESEE_HANG_DUMP_S=<s>   every <s> seconds, dump every thread's Python
#                            stack to stderr (faulthandler), so a rank stuck
#                            in a collective shows exactly where. Defaults to
#                            900 s when ICESEE_DEBUG_INIT is set; 0 disables.
#
# Both are off by default and cost nothing when off.
# ==============================================================================

from __future__ import annotations

import faulthandler
import os
import socket
import sys
import time
from contextlib import contextmanager

_WATCHDOG_ARMED = False


def init_trace_enabled() -> bool:
    return os.environ.get("ICESEE_DEBUG_INIT", "").strip().lower() not in ("", "0", "false", "no")


def _world():
    try:
        from mpi4py import MPI

        if MPI.Is_initialized():
            return MPI.COMM_WORLD.Get_rank(), MPI.COMM_WORLD.Get_size()
    except Exception:
        pass
    return -1, -1


def _fmt_fields(fields) -> str:
    parts = []
    for key, value in fields.items():
        if hasattr(value, "Get_size") and hasattr(value, "Get_rank"):
            value = f"{value.Get_rank()}/{value.Get_size()}"
        parts.append(f"{key}={value}")
    return " ".join(parts)


def trace_init(message: str, **fields) -> None:
    """Print one rank-labelled, flushed line when ICESEE_DEBUG_INIT is set.

    A field whose value is an MPI communicator is shown as ``rank/size``.
    """

    if not init_trace_enabled():
        return
    rank, size = _world()
    extra = _fmt_fields(fields)
    print(
        f"[ICESEE-init w{rank}/{size} {socket.gethostname()} t={time.time():.1f}] "
        f"{message}{' ' + extra if extra else ''}",
        file=sys.stdout,
        flush=True,
    )


@contextmanager
def traced_init(stage: str, **fields):
    """ENTER/EXIT trace (with elapsed seconds) around one initialization stage."""

    if not init_trace_enabled():
        yield
        return
    trace_init(f"ENTER {stage}", **fields)
    start = time.perf_counter()
    yield
    trace_init(f"EXIT  {stage} ({time.perf_counter() - start:.2f}s)", **fields)


def arm_hang_watchdog() -> None:
    """Arm the periodic all-thread stack dump (idempotent; opt-in)."""

    global _WATCHDOG_ARMED
    if _WATCHDOG_ARMED:
        return
    raw = os.environ.get("ICESEE_HANG_DUMP_S")
    if raw is None:
        seconds = 900.0 if init_trace_enabled() else 0.0
    else:
        try:
            seconds = float(raw)
        except ValueError:
            seconds = 0.0
    if seconds <= 0:
        return
    faulthandler.dump_traceback_later(seconds, repeat=True, file=sys.stderr)
    _WATCHDOG_ARMED = True
    trace_init(f"hang watchdog armed: stack dump every {seconds:g}s")
