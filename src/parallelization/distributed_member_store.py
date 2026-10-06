# ==============================================================================
# @des: Inactive-member state backing for execution-mode-3 bounded-memory
# round scheduling (Option B: one active native member per ensemble slot,
# inactive round-assigned members retained as packed distributed-memory
# state rather than live native model objects).
# @date: 2026-09-26
# ==============================================================================
"""Inactive-member state storage, decoupled from the runtime layer that uses it.

The scientific/runtime layer (see ``distributed_streaming_runtime.py``) never
needs to know HOW an inactive member's packed owned-state array is kept
around between rounds -- only that it can ``put``, ``get``, and ``release``
it by member ID. This keeps "state ownership" (the runtime's job) separate
from "inactive-state backing" (this module's job), which is itself separate
from "durable checkpointing" (``distributed_checkpoint.py``, unrelated and
untouched by this module).

Only a memory-resident backend is implemented here (Option B, the first
authorized implementation). A future file/checkpoint-backed backend
(Option C) or an HPC-specific backend can implement the same
``InactiveMemberStore`` protocol without any runtime-layer change -- that is
the entire point of this abstraction existing as its own module.

Scale note (Ne up to ~1000, 37 GB/member): whole-array ``get``/``put`` is
sufficient for Option B (bounding NATIVE/PETSc overhead to one live member)
but does NOT by itself bound how much PACKED state must be addressable at
once -- a caller that calls ``get`` for every one of 1000 members before
using any of them would still need ~37 TB addressable. The optional
``get_rows``/``put_rows`` methods below exist so a future backend (and a
future streaming analysis primitive) can pull/push only a bounded
state-row slice for one member without ever materializing that member's
complete array -- the array-level ``get``/``put`` remain the only
REQUIRED methods; ``get_rows``/``put_rows`` are duck-typed optional
extensions (checked via ``getattr(store, "get_rows", None)``, exactly like
this codebase's existing optional-adapter-method convention), so a simple
backend (like ``MemoryInactiveMemberStore``) can still satisfy the full
protocol by implementing them trivially in terms of the whole array.
"""

from __future__ import annotations

import os
from typing import Any, Mapping, Protocol, runtime_checkable

import numpy as np


@runtime_checkable
class InactiveMemberStore(Protocol):
    """Small contract any inactive-member backing must satisfy."""

    def put(self, member_id: int, local_state: np.ndarray) -> None:
        """Store (replacing any previous value) this member's packed owned state."""

    def get(self, member_id: int) -> np.ndarray:
        """Return this member's packed owned state. Raises if never put()."""

    def release(self, member_id: int) -> None:
        """Discard this member's stored state (no-op if never put())."""

    def has(self, member_id: int) -> bool:
        """Whether this member currently has stored state."""

    def get_rows(self, member_id: int, row_slice: slice) -> np.ndarray:
        """Optional: return only ``row_slice`` of this member's packed
        state, without requiring the caller (or this store) to ever hold
        the complete array at once. A memory backend can only ever
        implement this in terms of the full array (no saving there); a
        chunked file/HDF5-hyperslab backend is where this actually bounds
        memory, by reading only the requested slice from disk."""

    def put_rows(self, member_id: int, row_slice: slice, values: np.ndarray) -> None:
        """Optional: write only ``row_slice`` of this member's packed
        state (see ``get_rows``)."""

    def get_ensemble_rows(
        self, member_ids: tuple[int, ...], row_slice: slice
    ) -> np.ndarray:
        """Optional bulk primitive: return ``row_slice`` for ALL of
        ``member_ids`` in one call, shape ``(rows, len(member_ids))``.

        Exists so a backend that can serve this natively (e.g. a chunked
        2D state-row x ensemble-member layout) does not pay Ne separate
        physical I/O operations for what is logically one request -- see
        this session's Gate 7 finding: a filesystem backend doing 1000
        independent small reads per row block scales terribly even though
        an in-memory backend doing the same loop is cheap. Callers should
        prefer this over looping ``get_rows`` when available; a backend
        that cannot serve it natively may simply loop internally (see
        ``MemoryInactiveMemberStore`` below) -- the DA layer never needs to
        know which case applies.
        """

    def put_ensemble_rows(
        self, member_ids: tuple[int, ...], row_slice: slice, block: np.ndarray
    ) -> None:
        """Optional bulk primitive: write an analyzed ``(rows,
        len(member_ids))`` block back across all of ``member_ids`` in one
        call (see ``get_ensemble_rows``)."""


class MemoryInactiveMemberStore:
    """Option B's backend: a plain dict of packed owned-state arrays.

    O(owned_size) bytes per stored member -- never a global-sized array,
    never a live native model object (no Firedrake/PETSc handles retained
    here at all). This is strictly cheaper than keeping a live native
    member resident, but does NOT eliminate the underlying O(rounds_per_slot
    * owned_size) storage requirement -- it only removes the (potentially
    much larger) native/PETSc overhead on top of the raw array. See the
    large-state benchmark's Option-B measurements for how much that
    actually saves in practice; this class intentionally does not
    editorialize about it.
    """

    def __init__(self) -> None:
        self._states: dict[int, np.ndarray] = {}

    def put(self, member_id: int, local_state: np.ndarray) -> None:
        self._states[int(member_id)] = np.ascontiguousarray(
            np.asarray(local_state, dtype=np.float64)
        )

    def get(self, member_id: int) -> np.ndarray:
        try:
            return self._states[int(member_id)]
        except KeyError as error:
            raise KeyError(
                f"no stored inactive state for member {member_id}"
            ) from error

    def release(self, member_id: int) -> None:
        self._states.pop(int(member_id), None)

    def has(self, member_id: int) -> bool:
        return int(member_id) in self._states

    def get_rows(self, member_id: int, row_slice: slice) -> np.ndarray:
        # Memory backend: no bounded-read benefit (the full array is
        # already resident), but implementing this trivially keeps the
        # protocol satisfiable so runtime code can call get_rows uniformly
        # regardless of backend.
        return self.get(member_id)[row_slice]

    def put_rows(self, member_id: int, row_slice: slice, values: np.ndarray) -> None:
        array = self._states[int(member_id)]
        array[row_slice] = np.asarray(values, dtype=np.float64)

    def get_ensemble_rows(
        self, member_ids: tuple[int, ...], row_slice: slice
    ) -> np.ndarray:
        # No native bulk benefit for a memory backend (every member's full
        # array is already resident) -- implemented as a loop purely to
        # satisfy the protocol uniformly for callers that prefer it.
        columns = [self.get_rows(mid, row_slice) for mid in member_ids]
        return np.stack(columns, axis=1)

    def put_ensemble_rows(
        self, member_ids: tuple[int, ...], row_slice: slice, block: np.ndarray
    ) -> None:
        block = np.asarray(block)
        for i, mid in enumerate(member_ids):
            self.put_rows(mid, row_slice, block[:, i])

    def total_bytes(self) -> int:
        """Diagnostic only: current total bytes held by this store."""

        return sum(array.nbytes for array in self._states.values())


class InstrumentedInactiveMemberStore:
    """Transparent instrumenting wrapper around any ``InactiveMemberStore``
    (Gate 11): counts calls and bytes moved per operation category,
    without touching any backend's own implementation. Works uniformly for
    ``MemoryInactiveMemberStore``, ``HDF5MemberMajorStore``, or any future
    backend -- wrap once, measure anything.

    Categories tracked separately (Gate 12's forecast-vs-analysis
    distinction): ``whole`` (get/put -- activation/deactivation traffic),
    ``row`` (get_rows/put_rows -- per-member row-block traffic, used only
    when a backend lacks the bulk primitive), ``bulk``
    (get_ensemble_rows/put_ensemble_rows -- analysis traffic, one call
    covering many members at once).
    """

    def __init__(self, inner: Any) -> None:
        self._inner = inner
        self.stats: dict[str, float] = {
            "whole_get_count": 0, "whole_put_count": 0,
            "whole_get_bytes": 0, "whole_put_bytes": 0,
            "row_get_count": 0, "row_put_count": 0,
            "row_get_bytes": 0, "row_put_bytes": 0,
            "bulk_get_count": 0, "bulk_put_count": 0,
            "bulk_get_bytes": 0, "bulk_put_bytes": 0,
            "read_time_s": 0.0, "write_time_s": 0.0,
        }
        from ICESEE.src.utils.performance import register_io_provider
        register_io_provider("member_store", self.io_counters)

    def io_counters(self) -> dict[str, float]:
        """Cumulative counters in the generic performance-report I/O keys."""
        stats = self.stats
        return {
            "bytes_read": stats["whole_get_bytes"] + stats["row_get_bytes"] + stats["bulk_get_bytes"],
            "bytes_written": stats["whole_put_bytes"] + stats["row_put_bytes"] + stats["bulk_put_bytes"],
            "reads": stats["whole_get_count"] + stats["row_get_count"] + stats["bulk_get_count"],
            "writes": stats["whole_put_count"] + stats["row_put_count"] + stats["bulk_put_count"],
            "read_time_s": stats["read_time_s"],
            "write_time_s": stats["write_time_s"],
        }

    def _time(self, key: str):
        import time
        return _Timer(self.stats, key)

    def put(self, member_id: int, local_state: np.ndarray) -> None:
        array = np.asarray(local_state)
        with self._time("write_time_s"):
            self._inner.put(member_id, array)
        self.stats["whole_put_count"] += 1
        self.stats["whole_put_bytes"] += array.nbytes

    def get(self, member_id: int) -> np.ndarray:
        with self._time("read_time_s"):
            result = self._inner.get(member_id)
        self.stats["whole_get_count"] += 1
        self.stats["whole_get_bytes"] += np.asarray(result).nbytes
        return result

    def release(self, member_id: int) -> None:
        self._inner.release(member_id)

    def has(self, member_id: int) -> bool:
        return self._inner.has(member_id)

    def get_rows(self, member_id: int, row_slice: slice) -> np.ndarray:
        with self._time("read_time_s"):
            result = self._inner.get_rows(member_id, row_slice)
        self.stats["row_get_count"] += 1
        self.stats["row_get_bytes"] += np.asarray(result).nbytes
        return result

    def put_rows(self, member_id: int, row_slice: slice, values: np.ndarray) -> None:
        array = np.asarray(values)
        with self._time("write_time_s"):
            self._inner.put_rows(member_id, row_slice, array)
        self.stats["row_put_count"] += 1
        self.stats["row_put_bytes"] += array.nbytes

    def get_ensemble_rows(self, member_ids: tuple[int, ...], row_slice: slice) -> np.ndarray:
        with self._time("read_time_s"):
            result = self._inner.get_ensemble_rows(member_ids, row_slice)
        self.stats["bulk_get_count"] += 1
        self.stats["bulk_get_bytes"] += np.asarray(result).nbytes
        return result

    def put_ensemble_rows(self, member_ids: tuple[int, ...], row_slice: slice, block: np.ndarray) -> None:
        array = np.asarray(block)
        with self._time("write_time_s"):
            self._inner.put_ensemble_rows(member_ids, row_slice, array)
        self.stats["bulk_put_count"] += 1
        self.stats["bulk_put_bytes"] += array.nbytes

    def __getattr__(self, name: str) -> Any:
        # Forward anything not explicitly instrumented above (e.g.
        # file_paths(), close(), total_bytes()) straight to the wrapped
        # backend, so this wrapper is a transparent drop-in.
        return getattr(self._inner, name)


class _Timer:
    """Tiny context manager: adds elapsed wall time to stats[key]."""

    def __init__(self, stats: dict[str, float], key: str) -> None:
        self._stats = stats
        self._key = key

    def __enter__(self):
        import time
        self._t0 = time.perf_counter()
        return self

    def __exit__(self, *exc):
        import time
        self._stats[self._key] += time.perf_counter() - self._t0
        return False


def build_inactive_member_store(
    icesee_kwargs: Mapping[str, Any],
    *,
    world_rank: int,
) -> Any:
    """Backend-selection factory -- the ONLY place in the generic Mode-3
    infrastructure that imports a concrete non-memory backend.

    Applications (Icepack, or any future model) never import
    ``HDF5MemberMajorStore``/``HDF5Chunked2DStore`` directly and never see
    which backend is active; they only see the ``InactiveMemberStore``
    protocol via ``StreamingNativeDistributedMemberPool``.

    Config keys (icesee_kwargs), following this codebase's existing
    generic-CLI-override convention rather than inventing new controls:

      ``member_store_backend``: ``"memory"`` (default) | ``"hdf5_member_major"``.
      ``member_store_root``: directory for a file-backed backend. If unset,
        defaults to ``<data_path>/_mode3_member_store/<run_id>`` -- a
        run-local, portable location consistent with this codebase's
        existing output conventions (never ``/tmp``, never a hardcoded
        HPC path). Unique per run via ``run_id``; unique per rank via the
        backend's own ``rank_{rank:06d}.h5`` naming (see
        ``distributed_member_store_hdf5.py``), so simultaneous jobs (using
        different ``run_id``s, as production runs already do for
        checkpoints) cannot collide.

    Default (``member_store_backend`` unset) is unchanged: a fresh
    ``MemoryInactiveMemberStore``, identical to every prior behavior.
    """

    backend = str(icesee_kwargs.get("member_store_backend", "memory")).strip().lower()
    if backend == "memory":
        store: Any = MemoryInactiveMemberStore()
    elif backend == "hdf5_member_major":
        from .distributed_member_store_hdf5 import HDF5MemberMajorStore

        run_id = str(icesee_kwargs.get("run_id", "mode3-run"))
        data_path = str(icesee_kwargs.get("data_path", "."))
        root = icesee_kwargs.get("member_store_root") or os.path.join(
            data_path, "_mode3_member_store", run_id
        )
        # owned_size is accepted for interface symmetry with a future
        # backend that needs it upfront (e.g. a pre-shaped 2D dataset) but
        # is NOT used by HDF5MemberMajorStore itself: each member's HDF5
        # dataset is created lazily, sized from whatever array its own
        # first put() call supplies.
        store = HDF5MemberMajorStore(str(root), rank=int(world_rank), owned_size=0)
    else:
        raise ValueError(
            f"unknown member_store_backend {backend!r}; expected 'memory' or "
            "'hdf5_member_major'"
        )
    # Instrumentation overhead is a handful of dict increments and one
    # perf_counter call per operation -- negligible, so this wrapper is
    # applied unconditionally (Gate 11) rather than gated behind yet
    # another opt-in flag.
    return InstrumentedInactiveMemberStore(store)
