"""Model-agnostic distributed field packing for execution mode 3.

Models retain ownership of their native vectors.  ICESEE sees only compact
owned arrays and stable global block metadata; it never needs to import the
model's mesh, finite-element, or linear-algebra package.
"""

from __future__ import annotations

from collections import OrderedDict
from typing import Any, Hashable, Protocol, Sequence, runtime_checkable

import numpy as np

from .distributed_adapter import (
    DistributedBlockStateLayout,
    DistributedStateBlockLayout,
)


@runtime_checkable
class DistributedFieldHandle(Protocol):
    """Small contract implemented by a model-native distributed field."""

    name: str
    global_size: int
    owned_start: int
    owned_stop: int

    def read_owned(self) -> np.ndarray:
        """Return this rank's owned values in stable block-local ordering."""

    def write_owned(self, values: np.ndarray) -> None:
        """Replace this rank's owned values without writing ghost entries."""

    def synchronize_ghosts(self) -> None:
        """Refresh model-native read-only halos after owned values change."""


def _validate_field(field: DistributedFieldHandle) -> None:
    name = str(field.name)
    size = int(field.global_size)
    start = int(field.owned_start)
    stop = int(field.owned_stop)
    if not name:
        raise ValueError("distributed field name must be nonempty")
    if size < 0 or not 0 <= start <= stop <= size:
        raise ValueError(f"distributed field {name!r} has invalid ownership")
    for callback in ("read_owned", "write_owned", "synchronize_ghosts"):
        if not callable(getattr(field, callback, None)):
            raise TypeError(
                f"distributed field {name!r} is missing callback {callback}"
            )


class DistributedFieldRegistry:
    """Ordered model fields represented as one compact local ICESEE state.

    Global rows remain variable-major and compatible with ICESEE's existing
    packed-state convention.  Local storage contains only each field's owned
    interval, so memory is proportional to local degrees of freedom.
    """

    def __init__(
        self,
        fields: Sequence[DistributedFieldHandle],
        *,
        layout_id: str,
    ) -> None:
        ordered: OrderedDict[str, DistributedFieldHandle] = OrderedDict()
        blocks: list[DistributedStateBlockLayout] = []
        global_offset = 0
        for field in fields:
            _validate_field(field)
            name = str(field.name)
            if name in ordered:
                raise ValueError(f"duplicate distributed field name: {name}")
            ordered[name] = field
            blocks.append(
                DistributedStateBlockLayout(
                    name=name,
                    global_offset=global_offset,
                    global_size=int(field.global_size),
                    owned_start=int(field.owned_start),
                    owned_stop=int(field.owned_stop),
                )
            )
            global_offset += int(field.global_size)
        if not ordered:
            raise ValueError("distributed field registry requires at least one field")
        self._fields = ordered
        self._layout = DistributedBlockStateLayout(
            tuple(blocks), layout_id=str(layout_id)
        )

    @property
    def layout(self) -> DistributedBlockStateLayout:
        return self._layout

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(self._fields)

    def field(self, name: str) -> DistributedFieldHandle:
        try:
            return self._fields[str(name)]
        except KeyError as error:
            raise KeyError(f"unknown distributed field: {name}") from error

    def pack_owned(self, *, dtype=None) -> np.ndarray:
        """Copy native owned fields into one compact local state vector."""

        pieces = []
        for name, field in self._fields.items():
            values = np.asarray(field.read_owned())
            expected = self.layout.block(name).owned_size
            if values.ndim != 1 or values.size != expected:
                raise ValueError(
                    f"distributed field {name!r} returned shape {values.shape}; "
                    f"expected ({expected},)"
                )
            pieces.append(np.asarray(values, dtype=dtype) if dtype else values)
        return np.ascontiguousarray(np.concatenate(pieces))

    def owned_global_rows(self) -> np.ndarray:
        """Return stable variable-major row numbers owned by this rank.

        The array is proportional to local ownership, never to the global
        model dimension.  It is primarily useful when partitioning a global
        identity observation operator among spatial ranks.
        """

        return self.layout.local_to_global_rows(
            np.arange(self.layout.owned_size, dtype=np.int64)
        )

    def observe_owned_rows(
        self,
        global_rows: np.ndarray,
        *,
        dtype=None,
    ) -> np.ndarray:
        """Read identity observations for globally numbered owned rows.

        Requested row ordering and duplicates are preserved.  Every row must
        belong to this registry's local ownership; callers must partition an
        observation stencil before invoking this method.  The lookup walks
        the small list of variable blocks and never constructs a global-sized
        index or packed model state.
        """

        raw_rows = np.asarray(global_rows)
        if raw_rows.ndim != 1:
            raise ValueError("global observation rows must be one-dimensional")
        if not np.issubdtype(raw_rows.dtype, np.integer):
            raise TypeError("global observation rows must use an integer dtype")
        rows = raw_rows.astype(np.int64, copy=False)
        if rows.size and (rows.min() < 0 or rows.max() >= self.layout.global_size):
            raise ValueError("global observation rows leave the distributed state")

        output_dtype = np.dtype(dtype) if dtype is not None else None
        result = None
        matched = np.zeros(rows.size, dtype=bool)
        for name, local, (global_start, global_stop) in (
            self.layout.iter_owned_global_intervals()
        ):
            selected = (rows >= global_start) & (rows < global_stop)
            if not np.any(selected):
                continue
            values = np.asarray(self._fields[name].read_owned())
            expected = local.stop - local.start
            if values.ndim != 1 or values.size != expected:
                raise ValueError(
                    f"distributed field {name!r} returned shape {values.shape}; "
                    f"expected ({expected},)"
                )
            if result is None:
                selected_dtype = output_dtype or values.dtype
                result = np.empty(rows.size, dtype=selected_dtype)
            result[selected] = values[rows[selected] - global_start]
            matched[selected] = True

        if not np.all(matched):
            missing = rows[~matched]
            preview = ", ".join(str(int(row)) for row in missing[:5])
            suffix = "..." if missing.size > 5 else ""
            raise ValueError(
                "observation rows are not owned by this spatial rank: "
                f"{preview}{suffix}"
            )
        if result is None:
            return np.empty(0, dtype=output_dtype or np.float64)
        return np.ascontiguousarray(result)

    def partition_owned_rows(
        self,
        global_rows: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Select observation rows owned by this spatial rank.

        Returns the positions in the caller's canonical observation array and
        the corresponding globally numbered state rows.  Ordering and
        duplicates are preserved.  This is deliberately separate from
        :meth:`observe_owned_rows`: routing an irregular observation stencil
        must not require packing a global state or allocating an array whose
        size is the global model dimension.
        """

        raw_rows = np.asarray(global_rows)
        if raw_rows.ndim != 1:
            raise ValueError("global observation rows must be one-dimensional")
        if not np.issubdtype(raw_rows.dtype, np.integer):
            raise TypeError("global observation rows must use an integer dtype")
        rows = raw_rows.astype(np.int64, copy=False)
        if rows.size and (rows.min() < 0 or rows.max() >= self.layout.global_size):
            raise ValueError("global observation rows leave the distributed state")

        owned = np.zeros(rows.size, dtype=bool)
        for _name, _local, (global_start, global_stop) in (
            self.layout.iter_owned_global_intervals()
        ):
            owned |= (rows >= global_start) & (rows < global_stop)
        positions = np.flatnonzero(owned).astype(np.int64, copy=False)
        return (
            np.ascontiguousarray(positions),
            np.ascontiguousarray(rows[positions]),
        )

    def unpack_owned(
        self,
        local_state: np.ndarray,
        *,
        synchronize: bool = True,
    ) -> None:
        """Write a compact analysis/forecast state into native owned fields."""

        state = np.asarray(local_state)
        if state.ndim != 1 or state.size != self.layout.owned_size:
            raise ValueError(
                f"local state must have shape ({self.layout.owned_size},)"
            )
        for name, field in self._fields.items():
            field.write_owned(np.ascontiguousarray(state[self.layout.local_slice(name)]))
        if synchronize:
            self.synchronize_ghosts()

    def synchronize_ghosts(self) -> None:
        """Refresh halos once per unique model-native backing object."""

        synchronized: set[Hashable] = set()
        for field in self._fields.values():
            key: Any = getattr(field, "synchronization_key", id(field))
            try:
                hash(key)
            except TypeError as error:
                raise TypeError(
                    "distributed field synchronization_key must be hashable"
                ) from error
            if key in synchronized:
                continue
            field.synchronize_ghosts()
            synchronized.add(key)
