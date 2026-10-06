# ==============================================================================
# @des: Tests for the production backend-selection factory
# (build_inactive_member_store) and the generic instrumenting wrapper
# (InstrumentedInactiveMemberStore) -- Gate 1/2/11/20 of the member-major
# production-wiring investigation. No MPI/Firedrake needed.
# ==============================================================================
from __future__ import annotations

import os

import numpy as np
import pytest

from ICESEE.src.parallelization.distributed_member_store import (
    InstrumentedInactiveMemberStore,
    MemoryInactiveMemberStore,
    build_inactive_member_store,
)
from ICESEE.src.parallelization.distributed_member_store_hdf5 import HDF5MemberMajorStore


def test_default_backend_is_memory():
    store = build_inactive_member_store({}, world_rank=0)
    # Unwrap the instrumentation wrapper to check the concrete backend.
    assert isinstance(store._inner, MemoryInactiveMemberStore)


def test_explicit_memory_backend():
    store = build_inactive_member_store({"member_store_backend": "memory"}, world_rank=0)
    assert isinstance(store._inner, MemoryInactiveMemberStore)


def test_unknown_backend_raises():
    with pytest.raises(ValueError, match="unknown member_store_backend"):
        build_inactive_member_store({"member_store_backend": "nonsense"}, world_rank=0)


def test_hdf5_backend_selection_and_path(tmp_path):
    icesee_kwargs = {
        "member_store_backend": "hdf5_member_major",
        "data_path": str(tmp_path),
        "run_id": "test-run-1",
    }
    store = build_inactive_member_store(icesee_kwargs, world_rank=3)
    assert isinstance(store._inner, HDF5MemberMajorStore)
    (path,) = store.file_paths()
    assert path == os.path.join(str(tmp_path), "_mode3_member_store", "test-run-1", "member_major_rank000003.h5")
    store.close()


def test_hdf5_backend_default_root_uses_data_path_and_run_id(tmp_path):
    icesee_kwargs = {
        "member_store_backend": "hdf5_member_major",
        "data_path": str(tmp_path),
        "run_id": "another-run",
    }
    store = build_inactive_member_store(icesee_kwargs, world_rank=0)
    (path,) = store.file_paths()
    assert "_mode3_member_store" in path
    assert "another-run" in path
    store.close()


def test_two_ranks_get_isolated_files(tmp_path):
    """Store path isolation (Gate 2/5): different world_rank -> different file."""
    icesee_kwargs = {
        "member_store_backend": "hdf5_member_major",
        "data_path": str(tmp_path),
        "run_id": "iso-test",
    }
    store0 = build_inactive_member_store(icesee_kwargs, world_rank=0)
    store1 = build_inactive_member_store(icesee_kwargs, world_rank=1)
    assert store0.file_paths() != store1.file_paths()
    store0.put(0, np.array([1.0, 2.0]))
    store1.put(0, np.array([9.0, 9.0]))
    # Each rank's own data must be unaffected by the other's writes.
    np.testing.assert_array_equal(store0.get(0), [1.0, 2.0])
    np.testing.assert_array_equal(store1.get(0), [9.0, 9.0])
    store0.close()
    store1.close()


def test_instrumented_wrapper_counts_whole_member_operations():
    inner = MemoryInactiveMemberStore()
    store = InstrumentedInactiveMemberStore(inner)
    store.put(0, np.array([1.0, 2.0, 3.0]))
    store.get(0)
    assert store.stats["whole_put_count"] == 1
    assert store.stats["whole_get_count"] == 1
    assert store.stats["whole_put_bytes"] == 24
    assert store.stats["whole_get_bytes"] == 24


def test_instrumented_wrapper_counts_row_and_bulk_operations():
    inner = MemoryInactiveMemberStore()
    store = InstrumentedInactiveMemberStore(inner)
    store.put(0, np.zeros(10))
    store.put(1, np.zeros(10))
    store.get_rows(0, slice(0, 5))
    store.put_rows(0, slice(0, 5), np.ones(5))
    store.get_ensemble_rows((0, 1), slice(0, 5))
    store.put_ensemble_rows((0, 1), slice(0, 5), np.ones((5, 2)))
    assert store.stats["row_get_count"] == 1
    assert store.stats["row_put_count"] == 1
    assert store.stats["bulk_get_count"] == 1
    assert store.stats["bulk_put_count"] == 1
    assert store.stats["bulk_get_bytes"] == 5 * 2 * 8
    assert store.stats["read_time_s"] >= 0.0
    assert store.stats["write_time_s"] >= 0.0


def test_instrumented_wrapper_forwards_unknown_attributes(tmp_path):
    inner = HDF5MemberMajorStore(str(tmp_path), rank=0, owned_size=0)
    store = InstrumentedInactiveMemberStore(inner)
    assert store.file_paths() == inner.file_paths()
    store.close()
