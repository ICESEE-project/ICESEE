# Rank-local HDF5 inactive-member stores; serial files are private to each rank.
"""Rank-local HDF5 backends for inactive member state.

HDF5MemberMajorStore uses one dataset per member. HDF5Chunked2DStore uses one
chunked (owned_rows, members) dataset for bulk ensemble-row operations.
Each rank owns one serial HDF5 file; no parallel h5py driver is required.

Releasing entries does not shrink the HDF5 file. Long-lived stores may need
repacking or rotation. These temporary stores are not restart checkpoints."""

from __future__ import annotations

import os
from typing import Any

import h5py
import numpy as np


class HDF5MemberMajorStore:
    """Gate 3 prototype: one dataset per member, in one file per rank."""

    def __init__(self, directory: str, rank: int, owned_size: int, *, run_id: str = "member_major") -> None:
        os.makedirs(directory, exist_ok=True)
        self._path = os.path.join(directory, f"{run_id}_rank{int(rank):06d}.h5")
        self._owned_size = int(owned_size)
        self._file = h5py.File(self._path, "a")
        # Dataset-handle cache (Phase 4): file[name] re-resolves the HDF5
        # object each call even for an already-open dataset -- cheap, but
        # not free, and this store's access pattern re-reads/re-writes the
        # SAME small set of member datasets every timestep. Handles are
        # only ever used while self._file is open (cleared in release()
        # and never handed outside this class), so this never outlives
        # the file's own lifecycle.
        self._handles: dict[str, "h5py.Dataset"] = {}

    def flush(self) -> None:
        self._file.flush()

    def close(self) -> None:
        self._handles.clear()
        self._file.close()

    def __del__(self) -> None:
        try:
            self._handles.clear()
            self._file.close()
        except Exception:
            pass

    def _name(self, member_id: int) -> str:
        return f"member_{int(member_id)}"

    def _handle(self, name: str) -> "h5py.Dataset | None":
        handle = self._handles.get(name)
        if handle is not None:
            return handle
        if name not in self._file:
            return None
        handle = self._file[name]
        self._handles[name] = handle
        return handle

    def put(self, member_id: int, local_state: np.ndarray) -> None:
        name = self._name(member_id)
        array = np.ascontiguousarray(np.asarray(local_state, dtype=np.float64))
        handle = self._handle(name)
        if handle is not None:
            handle[:] = array
        else:
            self._handles[name] = self._file.create_dataset(name, data=array, dtype="f8")

    def get(self, member_id: int) -> np.ndarray:
        name = self._name(member_id)
        handle = self._handle(name)
        if handle is None:
            raise KeyError(f"no stored inactive state for member {member_id}")
        return np.asarray(handle[:])

    def release(self, member_id: int) -> None:
        name = self._name(member_id)
        if name in self._file:
            del self._file[name]
        self._handles.pop(name, None)

    def has(self, member_id: int) -> bool:
        return self._name(member_id) in self._file

    def get_rows(self, member_id: int, row_slice: slice) -> np.ndarray:
        name = self._name(member_id)
        handle = self._handle(name)
        if handle is None:
            raise KeyError(f"no stored inactive state for member {member_id}")
        return np.asarray(handle[row_slice])

    def put_rows(self, member_id: int, row_slice: slice, values: np.ndarray) -> None:
        name = self._name(member_id)
        handle = self._handle(name)
        if handle is None:
            raise KeyError(f"no stored inactive state for member {member_id}")
        handle[row_slice] = np.asarray(values, dtype=np.float64)

    def get_ensemble_rows(self, member_ids: tuple[int, ...], row_slice: slice) -> np.ndarray:
        # No native bulk primitive for member-major storage: one physical
        # HDF5 read per member (the exact "Ne-separate-read" cost Gate 7
        # asks to measure, not hide).
        columns = [self.get_rows(mid, row_slice) for mid in member_ids]
        return np.stack(columns, axis=1)

    def put_ensemble_rows(self, member_ids: tuple[int, ...], row_slice: slice, block: np.ndarray) -> None:
        block = np.asarray(block)
        for i, mid in enumerate(member_ids):
            self.put_rows(mid, row_slice, block[:, i])

    def file_paths(self) -> list[str]:
        return [self._path]


class HDF5Chunked2DStore:
    """Gate 4 prototype: one 2D (state-row x ensemble-member) chunked
    dataset per rank, in one file per rank."""

    def __init__(
        self,
        directory: str,
        rank: int,
        owned_size: int,
        number_of_members: int,
        *,
        state_chunk: int = 4096,
        member_chunk: int = 32,
        run_id: str = "chunked2d",
    ) -> None:
        os.makedirs(directory, exist_ok=True)
        self._path = os.path.join(directory, f"{run_id}_rank{int(rank):06d}.h5")
        self._owned_size = int(owned_size)
        self._number_of_members = int(number_of_members)
        self._file = h5py.File(self._path, "a")
        chunk_shape = (
            min(int(state_chunk), self._owned_size),
            min(int(member_chunk), self._number_of_members),
        )
        if "ensemble" not in self._file:
            self._file.create_dataset(
                "ensemble",
                shape=(self._owned_size, self._number_of_members),
                dtype="f8",
                chunks=chunk_shape,
                fillvalue=np.nan,
            )
        self._dataset = self._file["ensemble"]
        self._known_members: set[int] = set()

    def flush(self) -> None:
        self._file.flush()

    def close(self) -> None:
        self._file.close()

    def __del__(self) -> None:
        try:
            self._file.close()
        except Exception:
            pass

    def put(self, member_id: int, local_state: np.ndarray) -> None:
        self._dataset[:, int(member_id)] = np.asarray(local_state, dtype=np.float64)
        self._known_members.add(int(member_id))

    def get(self, member_id: int) -> np.ndarray:
        if int(member_id) not in self._known_members:
            raise KeyError(f"no stored inactive state for member {member_id}")
        return np.asarray(self._dataset[:, int(member_id)])

    def release(self, member_id: int) -> None:
        self._known_members.discard(int(member_id))

    def has(self, member_id: int) -> bool:
        return int(member_id) in self._known_members

    def get_rows(self, member_id: int, row_slice: slice) -> np.ndarray:
        if int(member_id) not in self._known_members:
            raise KeyError(f"no stored inactive state for member {member_id}")
        return np.asarray(self._dataset[row_slice, int(member_id)])

    def put_rows(self, member_id: int, row_slice: slice, values: np.ndarray) -> None:
        self._dataset[row_slice, int(member_id)] = np.asarray(values, dtype=np.float64)
        self._known_members.add(int(member_id))

    def get_ensemble_rows(self, member_ids: tuple[int, ...], row_slice: slice) -> np.ndarray:
        ids = list(member_ids)
        # Contiguous, sorted member IDs (the common case: all Ne members)
        # can be served as ONE plain hyperslab read; otherwise fall back to
        # h5py's fancy indexing along the member axis -- still ONE physical
        # HDF5 call, not a Python-level loop.
        if ids == list(range(ids[0], ids[0] + len(ids))):
            return np.asarray(self._dataset[row_slice, ids[0]:ids[0] + len(ids)])
        return np.asarray(self._dataset[row_slice, ids])

    def put_ensemble_rows(self, member_ids: tuple[int, ...], row_slice: slice, block: np.ndarray) -> None:
        ids = list(member_ids)
        block = np.asarray(block)
        if ids == list(range(ids[0], ids[0] + len(ids))):
            self._dataset[row_slice, ids[0]:ids[0] + len(ids)] = block
        else:
            self._dataset[row_slice, ids] = block
        self._known_members.update(ids)

    def file_paths(self) -> list[str]:
        return [self._path]
