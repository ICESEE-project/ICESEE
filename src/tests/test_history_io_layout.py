"""Check member/time storage access and atomic snapshot finalization."""
import h5py
import numpy as np
import pytest

from ICESEE.src.parallelization._parallel_i_o import parallel_write_full_ensemble_from_root
from ICESEE.src.utils.tools import build_vds


class RootComm:
    def Get_rank(self):
        return 0

    def bcast(self, value, root=0):
        return value

    def Barrier(self):
        pass


def test_history_chunks_isolate_each_member_and_time(tmp_path):
    initial = np.arange(21, dtype=np.float64).reshape(7, 3)
    advanced = initial + 0.125
    kwargs = dict(data_path=str(tmp_path), Nens=3, nt=2)
    for t, values in enumerate((initial, advanced)):
        parallel_write_full_ensemble_from_root(
            t, values.mean(axis=1), kwargs, values, RootComm()
        )
    with h5py.File(tmp_path / 'icesee_ensemble_data.h5', 'r') as f:
        assert f['ensemble'].chunks == (7, 1, 1)
        assert f['ensemble_mean'].chunks == (7, 1)
        for t, values in enumerate((initial, advanced)):
            np.testing.assert_array_equal(f['ensemble'][:, :, t], values)
            np.testing.assert_array_equal(f['ensemble_mean'][:, t], values.mean(axis=1))
            for member in range(3):
                np.testing.assert_array_equal(f['ensemble'][:, member, t], values[:, member])


def test_existing_contiguous_history_remains_compatible(tmp_path):
    values = np.arange(15, dtype=np.float64).reshape(5, 3)
    path = tmp_path / 'icesee_ensemble_data.h5'
    with h5py.File(path, 'w') as f:
        f.create_dataset('ensemble', (5, 3, 2), dtype='f8')
        f.create_dataset('ensemble_mean', (5, 2), dtype='f8')
    parallel_write_full_ensemble_from_root(
        1, values.mean(axis=1), dict(data_path=str(tmp_path), Nens=3, nt=1),
        values, RootComm(),
    )
    with h5py.File(path, 'r') as f:
        assert f['ensemble'].chunks is None
        np.testing.assert_array_equal(f['ensemble'][:, :, 1], values)


def _snapshots(tmp_path):
    arrays = []
    for t in range(3):
        arr = np.arange(35, dtype=np.float64).reshape(7, 5) + t / 8
        arr[2, 1] = np.nan
        with h5py.File(tmp_path / f'icesee_enkf_ens_{t:04d}.h5', 'w') as f:
            f.create_dataset('states', data=arr)
        arrays.append(arr)
    return arrays


def test_vds_finalization_reads_source_slabs_not_virtual_history(tmp_path, monkeypatch):
    arrays = _snapshots(tmp_path)
    original = h5py.Dataset.__getitem__
    selections = []

    def read(ds, selection):
        assert not ds.is_virtual, 'Finalization must not read through its VDS'
        if ds.name == '/states':
            assert len(selection) == 2
            assert selection[0].stop - selection[0].start <= 2
            selections.append(selection)
        return original(ds, selection)

    with monkeypatch.context() as patch:
        patch.setattr(h5py.Dataset, '__getitem__', read)
        out = build_vds(str(tmp_path), dset_name='states', row_chunk_size=2)
    assert len(selections) == 12
    with h5py.File(out, 'r') as f:
        assert f['ensemble'].is_virtual
        np.testing.assert_array_equal(f['states'][:], np.stack(arrays, axis=2))
        np.testing.assert_array_equal(
            f['ensemble_mean'][:], np.stack([np.nanmean(a, axis=1) for a in arrays], axis=1)
        )


def test_vds_failure_keeps_previous_output_and_sources(tmp_path):
    arrays = _snapshots(tmp_path)
    out = tmp_path / 'icesee_ensemble_data.h5'
    out.write_bytes(b'previous published output')
    bad = tmp_path / 'icesee_enkf_ens_0001.h5'
    with h5py.File(bad, 'w') as f:
        f.create_dataset('states', data=np.zeros((6, 5)))
    with pytest.raises(ValueError, match='Incompatible snapshot'):
        build_vds(str(tmp_path), dset_name='states', row_chunk_size=2)
    assert out.read_bytes() == b'previous published output'
    assert len(list(tmp_path.glob('icesee_enkf_ens_*.h5'))) == len(arrays)
    assert not list(tmp_path.glob('.*.tmp'))
