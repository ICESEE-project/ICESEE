from __future__ import annotations

import numpy as np
import pytest

from src.parallelization.distributed_fields import DistributedFieldRegistry


class _Field:
    def __init__(self, name, values, *, global_size, owned_start, sync_key=None):
        self.name = name
        self.values = np.asarray(values, dtype=float)
        self.global_size = global_size
        self.owned_start = owned_start
        self.owned_stop = owned_start + self.values.size
        self.synchronization_key = id(self) if sync_key is None else sync_key
        self.sync_count = 0

    def read_owned(self):
        return self.values.copy()

    def write_owned(self, values):
        self.values[:] = values

    def synchronize_ghosts(self):
        self.sync_count += 1


def test_registry_packs_variable_major_owned_pieces_and_unpacks():
    h = _Field("h", [1, 2, 3], global_size=8, owned_start=2)
    u = _Field("u", [10, 11], global_size=5, owned_start=1)
    registry = DistributedFieldRegistry([h, u], layout_id="ice-native-v1")

    assert registry.names == ("h", "u")
    assert registry.layout.global_size == 13
    assert registry.layout.owned_size == 5
    assert registry.layout.block("h").global_owned_interval == (2, 5)
    assert registry.layout.block("u").global_owned_interval == (9, 11)
    np.testing.assert_array_equal(registry.pack_owned(), [1, 2, 3, 10, 11])

    registry.unpack_owned(np.asarray([4, 5, 6, 12, 13], dtype=float))
    np.testing.assert_array_equal(h.values, [4, 5, 6])
    np.testing.assert_array_equal(u.values, [12, 13])
    assert h.sync_count == 1
    assert u.sync_count == 1


def test_registry_synchronizes_shared_backing_once():
    u = _Field("u", [1, 2], global_size=4, owned_start=0, sync_key=77)
    v = _Field("v", [3, 4], global_size=4, owned_start=0, sync_key=77)
    registry = DistributedFieldRegistry([u, v], layout_id="velocity-v1")

    registry.unpack_owned(registry.pack_owned())
    assert u.sync_count + v.sync_count == 1


def test_registry_rejects_duplicate_names_and_bad_read_shape():
    field = _Field("h", [1, 2], global_size=2, owned_start=0)
    with pytest.raises(ValueError, match="duplicate"):
        DistributedFieldRegistry([field, field], layout_id="bad")

    registry = DistributedFieldRegistry([field], layout_id="bad-read")
    field.values = np.asarray([[1, 2]], dtype=float)
    with pytest.raises(ValueError, match="returned shape"):
        registry.pack_owned()


def test_registry_can_skip_ghost_synchronization():
    field = _Field("bed", [1], global_size=3, owned_start=1)
    registry = DistributedFieldRegistry([field], layout_id="bed-v1")
    registry.unpack_owned(np.asarray([9.0]), synchronize=False)
    assert field.values[0] == 9.0
    assert field.sync_count == 0


def test_registry_memory_scales_with_owned_not_global_entries():
    """A huge logical state must never cause a global NumPy allocation."""

    trillion = 10**12
    h = _Field("h", np.arange(7.0), global_size=trillion, owned_start=41)
    u = _Field("u", np.arange(5.0), global_size=trillion, owned_start=19)
    registry = DistributedFieldRegistry([h, u], layout_id="huge-logical-state-v1")

    packed = registry.pack_owned(dtype=np.float32)
    assert registry.layout.global_size == 2 * trillion
    assert registry.layout.owned_size == 12
    assert packed.shape == (12,)
    assert packed.nbytes == 12 * np.dtype(np.float32).itemsize


def test_registry_observes_global_rows_without_global_packing():
    h = _Field("h", [1, 2, 3], global_size=8, owned_start=2)
    u = _Field("u", [10, 11], global_size=5, owned_start=1)
    registry = DistributedFieldRegistry([h, u], layout_id="observed-v1")

    np.testing.assert_array_equal(
        registry.owned_global_rows(), [2, 3, 4, 9, 10]
    )
    np.testing.assert_array_equal(
        registry.observe_owned_rows(np.asarray([10, 2, 9, 2])),
        [11, 1, 10, 1],
    )

    positions, rows = registry.partition_owned_rows(
        np.asarray([0, 10, 2, 12, 9, 2], dtype=np.int64)
    )
    np.testing.assert_array_equal(positions, [1, 2, 4, 5])
    np.testing.assert_array_equal(rows, [10, 2, 9, 2])
    np.testing.assert_array_equal(
        registry.observe_owned_rows(rows), [11, 1, 10, 1]
    )


def test_registry_observation_rejects_remote_or_invalid_rows():
    field = _Field("h", [1, 2], global_size=8, owned_start=2)
    registry = DistributedFieldRegistry([field], layout_id="observed-errors-v1")

    with pytest.raises(ValueError, match="not owned"):
        registry.observe_owned_rows(np.asarray([1], dtype=np.int64))
    with pytest.raises(ValueError, match="leave the distributed state"):
        registry.observe_owned_rows(np.asarray([8], dtype=np.int64))
    with pytest.raises(TypeError, match="integer dtype"):
        registry.observe_owned_rows(np.asarray([2.0]))
    with pytest.raises(ValueError, match="one-dimensional"):
        registry.observe_owned_rows(np.asarray([[2]], dtype=np.int64))

    with pytest.raises(ValueError, match="leave the distributed state"):
        registry.partition_owned_rows(np.asarray([-1], dtype=np.int64))
    with pytest.raises(TypeError, match="integer dtype"):
        registry.partition_owned_rows(np.asarray([2.0]))
    with pytest.raises(ValueError, match="one-dimensional"):
        registry.partition_owned_rows(np.asarray([[2]], dtype=np.int64))


def test_registry_partition_owned_rows_never_scales_with_global_dimension():
    trillion = 10**12
    field = _Field("h", [1, 2, 3], global_size=trillion, owned_start=41)
    registry = DistributedFieldRegistry([field], layout_id="huge-observations-v1")

    positions, rows = registry.partition_owned_rows(
        np.asarray([0, 41, trillion - 1, 43], dtype=np.int64)
    )
    np.testing.assert_array_equal(positions, [1, 3])
    np.testing.assert_array_equal(rows, [41, 43])
    assert positions.nbytes + rows.nbytes == 4 * np.dtype(np.int64).itemsize
