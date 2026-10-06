# ==============================================================================
# @des: Correctness tests for the two Option-C prototype file-backed
# InactiveMemberStore backends (Gates 3/4 of the store-backend
# investigation) against the reference MemoryInactiveMemberStore -- pure
# I/O correctness, no MPI/Firedrake needed.
# ==============================================================================
from __future__ import annotations

import numpy as np
import pytest

from ICESEE.src.parallelization.distributed_member_store import MemoryInactiveMemberStore
from ICESEE.src.parallelization.distributed_member_store_hdf5 import (
    HDF5Chunked2DStore,
    HDF5MemberMajorStore,
)


@pytest.fixture(params=["member_major", "chunked2d"])
def store(request, tmp_path):
    owned_size = 20
    number_of_members = 4
    if request.param == "member_major":
        s = HDF5MemberMajorStore(str(tmp_path), rank=0, owned_size=owned_size)
    else:
        s = HDF5Chunked2DStore(
            str(tmp_path), rank=0, owned_size=owned_size,
            number_of_members=number_of_members, state_chunk=8, member_chunk=2,
        )
    yield s
    s.close()


def test_put_get_round_trip(store):
    values = np.arange(20.0)
    store.put(3, values)
    assert store.has(3)
    np.testing.assert_array_equal(store.get(3), values)


def test_get_missing_member_raises(store):
    with pytest.raises(KeyError):
        store.get(99)


def test_get_rows_matches_full_array_slice(store):
    values = np.arange(20.0) * 2.0
    store.put(1, values)
    np.testing.assert_array_equal(store.get_rows(1, slice(4, 10)), values[4:10])


def test_put_rows_updates_only_the_slice(store):
    values = np.zeros(20)
    store.put(2, values)
    store.put_rows(2, slice(5, 8), np.array([9.0, 9.0, 9.0]))
    expected = np.zeros(20)
    expected[5:8] = 9.0
    np.testing.assert_array_equal(store.get(2), expected)


def test_release_removes_member(store):
    store.put(0, np.ones(20))
    assert store.has(0)
    store.release(0)
    assert not store.has(0)


def test_get_ensemble_rows_matches_reference_memory_store():
    """Both prototype backends' bulk primitive must agree exactly with a
    reference MemoryInactiveMemberStore's per-member loop, for both
    contiguous (0..N-1) and non-contiguous member ID sets."""
    owned_size = 24
    number_of_members = 6
    member_ids_contiguous = (0, 1, 2, 3, 4, 5)
    member_ids_sparse = (0, 2, 5)

    ref = MemoryInactiveMemberStore()
    data = {}
    for mid in range(number_of_members):
        arr = np.arange(owned_size, dtype=np.float64) + 100.0 * mid
        data[mid] = arr
        ref.put(mid, arr)

    for backend_name in ("member_major", "chunked2d"):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            if backend_name == "member_major":
                store = HDF5MemberMajorStore(tmp, rank=0, owned_size=owned_size)
            else:
                store = HDF5Chunked2DStore(
                    tmp, rank=0, owned_size=owned_size,
                    number_of_members=number_of_members, state_chunk=8, member_chunk=2,
                )
            for mid, arr in data.items():
                store.put(mid, arr)

            for member_ids in (member_ids_contiguous, member_ids_sparse):
                expected = ref.get_ensemble_rows(member_ids, slice(3, 10))
                actual = store.get_ensemble_rows(member_ids, slice(3, 10))
                np.testing.assert_array_equal(actual, expected, err_msg=backend_name)
            store.close()


def test_put_ensemble_rows_matches_reference_memory_store():
    owned_size = 16
    number_of_members = 4
    member_ids = (0, 1, 2, 3)
    block = np.arange(5 * 4, dtype=np.float64).reshape(5, 4)

    ref = MemoryInactiveMemberStore()
    for mid in range(number_of_members):
        ref.put(mid, np.zeros(owned_size))
    ref.put_ensemble_rows(member_ids, slice(2, 7), block)

    for backend_name in ("member_major", "chunked2d"):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            if backend_name == "member_major":
                store = HDF5MemberMajorStore(tmp, rank=0, owned_size=owned_size)
            else:
                store = HDF5Chunked2DStore(
                    tmp, rank=0, owned_size=owned_size,
                    number_of_members=number_of_members, state_chunk=8, member_chunk=2,
                )
            for mid in range(number_of_members):
                store.put(mid, np.zeros(owned_size))
            store.put_ensemble_rows(member_ids, slice(2, 7), block)
            for mid in range(number_of_members):
                np.testing.assert_array_equal(
                    store.get(mid), ref.get(mid), err_msg=f"{backend_name} member {mid}"
                )
            store.close()


def test_file_count_is_one_per_rank_regardless_of_member_count(tmp_path):
    owned_size = 10
    number_of_members = 50
    store_mm = HDF5MemberMajorStore(str(tmp_path / "mm"), rank=0, owned_size=owned_size)
    for mid in range(number_of_members):
        store_mm.put(mid, np.zeros(owned_size))
    assert len(store_mm.file_paths()) == 1
    store_mm.close()

    store_2d = HDF5Chunked2DStore(
        str(tmp_path / "2d"), rank=0, owned_size=owned_size,
        number_of_members=number_of_members, state_chunk=4, member_chunk=8,
    )
    for mid in range(number_of_members):
        store_2d.put(mid, np.zeros(owned_size))
    assert len(store_2d.file_paths()) == 1
    store_2d.close()
