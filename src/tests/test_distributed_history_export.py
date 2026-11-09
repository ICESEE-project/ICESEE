"""Export real HDF5 checkpoint fixtures with reordered, spatially split fields."""
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import pytest

from ICESEE.src.parallelization.distributed_history_export import (
    export_ensemble_history, finalize_ensemble_history, prepare_reference_layout,
    save_reference_layout,
)


@pytest.fixture
def history(tmp_path):
    root = tmp_path / '_mode3_state_history'
    blocks = [dict(name='h', global_offset=0, global_size=4),
              dict(name='forcing', global_offset=4, global_size=2)]
    maps = {'h': np.array([2, 0, 3, 1]), 'forcing': np.array([1, 0])}
    reference = np.empty((6, 3, 4))
    for column in range(4):
        checkpoint = (root / 'initial_condition/checkpoint_00000000' if column == 0 else
                      root / 'steps' / f'checkpoint_{column - 1:08d}')
        checkpoint.mkdir(parents=True)
        native = np.arange(6)[:, None] * 100 + np.arange(3)[None, :] * 10 + column
        reference[:4, :, column] = native[maps['h']]
        reference[4:, :, column] = native[4 + maps['forcing']]
        shards = []
        for group, members in enumerate(([0, 2], [1])):
            for rank in range(2):
                owned = [dict(block, owned_start=rank * (block['global_size'] // 2),
                              owned_stop=(rank + 1) * (block['global_size'] // 2))
                         for block in blocks]
                name = f'shard_{group * 2 + rank:08d}.h5'
                with h5py.File(checkpoint / name, 'w') as handle:
                    for member in members:
                        for block in owned:
                            start = block['global_offset'] + block['owned_start']
                            stop = block['global_offset'] + block['owned_stop']
                            handle.create_dataset(f"member_{member}/{block['name']}",
                                                  data=native[start:stop, member].astype('f8'))
                shards.append(dict(file=name, member_ids=members, blocks=owned,
                                   layout_fingerprint=f'partition-{rank}'))
        manifest = dict(complete=True, run_id='run-1', layout_id='mesh-1',
                        global_size=6, number_of_members=3, dtype='<f8', blocks=blocks,
                        timestep=max(column - 1, 0), shards=shards,
                        metadata=dict(nt=3, next_step=column, completed_step=column - 1,
                                      model_time=column * 0.05, did_analysis=(column == 2)))
        (checkpoint / 'manifest.json').write_text(json.dumps(manifest))
    save_reference_layout(root, run_id='run-1', layout_id='mesh-1', row_maps=maps)
    return tmp_path, reference, maps


def test_materialized_history_is_reference_ordered_complete_and_standalone(history, monkeypatch):
    data, reference, _ = history
    original = h5py.Dataset.__getitem__
    reads = []

    def bounded_read(dataset, selection):
        if Path(dataset.file.filename).name.startswith('shard_'):
            assert isinstance(selection, np.ndarray)
            assert len(selection) <= 2
            reads.append(len(selection))
        return original(dataset, selection)

    with monkeypatch.context() as patch:
        patch.setattr(h5py.Dataset, '__getitem__', bounded_read)
        output = export_ensemble_history(data, row_chunk_size=2)
    assert reads
    shutil.rmtree(data / '_mode3_state_history')
    with h5py.File(output, 'r') as handle:
        assert handle.attrs['complete']
        assert handle.attrs['state_order'] == 'serial_reference'
        assert not handle['ensemble'].is_virtual
        assert handle['ensemble'].chunks == (2, 1, 1)
        np.testing.assert_array_equal(handle['ensemble'][:], reference)
        np.testing.assert_array_equal(handle['ensemble_mean'][:], reference.mean(axis=1))
        np.testing.assert_allclose(handle['t'][:], np.arange(4) * 0.05)
        np.testing.assert_array_equal(handle['checkpoint_steps'][:], [-1, 0, 1, 2])
        np.testing.assert_array_equal(handle['analysis_applied'][:], [False, False, True, False])


def test_pruned_snapshot_refuses_partial_history_and_preserves_previous_output(history):
    data, _, _ = history
    output = data / 'icesee_ensemble_data.h5'
    output.write_bytes(b'previous complete output')
    shutil.rmtree(data / '_mode3_state_history/steps/checkpoint_00000001')
    with pytest.raises(ValueError, match='missing history snapshot'):
        export_ensemble_history(data)
    assert output.read_bytes() == b'previous complete output'
    assert not list(data.glob('.icesee_ensemble_data.h5.*.tmp'))


def test_incomplete_run_is_not_published(history):
    data, _, _ = history
    shutil.rmtree(data / '_mode3_state_history/steps/checkpoint_00000002')
    with pytest.raises(ValueError, match='run is not complete'):
        export_ensemble_history(data)
    assert not (data / 'icesee_ensemble_data.h5').exists()


def test_corrupt_shard_failure_keeps_old_output_and_checkpoint_files(history):
    data, _, _ = history
    output = data / 'icesee_ensemble_data.h5'
    output.write_bytes(b'old output')
    source = data / '_mode3_state_history/steps/checkpoint_00000002/shard_00000000.h5'
    with h5py.File(source, 'a') as handle:
        del handle['member_0/h']
        handle['member_0/h'] = np.zeros(1)
    with pytest.raises(ValueError, match='shape/dtype'):
        export_ensemble_history(data, row_chunk_size=2)
    assert output.read_bytes() == b'old output'
    assert source.exists()
    assert not list(data.glob('.icesee_ensemble_data.h5.*.tmp'))


@pytest.mark.parametrize('problem', ['overlap', 'wrong_member', 'wrong_time', 'missing_shard'])
def test_invalid_snapshot_is_rejected(history, problem):
    data, _, _ = history
    path = data / '_mode3_state_history/steps/checkpoint_00000001/manifest.json'
    manifest = json.loads(path.read_text())
    if problem == 'overlap':
        manifest['shards'][1]['blocks'][0]['owned_start'] = 1
    elif problem == 'wrong_member':
        manifest['shards'][0]['member_ids'] = [8]
    elif problem == 'wrong_time':
        manifest['metadata']['next_step'] = 99
    else:
        manifest['shards'][0]['file'] = 'absent.h5'
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        export_ensemble_history(data)
    assert not (data / 'icesee_ensemble_data.h5').exists()


def test_reference_mapping_is_immutable_on_resume(history):
    data, _, maps = history
    root = data / '_mode3_state_history'
    save_reference_layout(root, run_id='run-1', layout_id='mesh-1', row_maps=maps)
    with pytest.raises(ValueError, match='numbering changed'):
        save_reference_layout(root, run_id='run-1', layout_id='mesh-1',
                              row_maps=dict(maps, h=np.arange(4)))
    with h5py.File(root / 'reference_layout.h5', 'r') as handle:
        np.testing.assert_array_equal(handle['fields/h'][:], maps['h'])


def test_export_rejects_invalid_reference_permutation(history):
    data, _, _ = history
    with h5py.File(data / '_mode3_state_history/reference_layout.h5', 'a') as handle:
        handle['fields/h'][:] = [0, 0, 1, 2]
    with pytest.raises(ValueError, match='invalid reference row permutation'):
        export_ensemble_history(data)


class Comm:
    def __init__(self, rank=0, broadcast=None, identities=None):
        self.rank, self.broadcast, self.identities = rank, broadcast, identities
        self.sent = None

    def Get_rank(self):
        return self.rank

    def bcast(self, value, root=0):
        self.sent = value
        return value if self.rank == 0 else self.broadcast

    def allgather(self, value):
        return self.identities or [value, value]


def test_export_error_is_broadcast_to_peers(tmp_path):
    root = Comm()
    with pytest.raises(RuntimeError, match='no committed step checkpoints'):
        finalize_ensemble_history(root, tmp_path)
    peer = Comm(rank=1, broadcast=root.sent)
    with pytest.raises(RuntimeError, match='no committed step checkpoints'):
        finalize_ensemble_history(peer, tmp_path)


def test_mapping_disagreement_fails_before_writing(tmp_path):
    with pytest.raises(ValueError, match='different reference/native'):
        prepare_reference_layout(Comm(identities=['a', 'b']), tmp_path,
                                 run_id='run', layout_id='layout', row_maps={'h': np.arange(3)})
    assert not (tmp_path / 'reference_layout.h5').exists()


def test_pig_notebook_loads_selected_full_history_without_configuration_import(history):
    import os
    from ICESEE.src.utils.tools import extract_datasets_from_h5

    data, reference, _ = history
    export_ensemble_history(data, row_chunk_size=2)
    with h5py.File(data / 'true-wrong-icepack.h5', 'w') as handle:
        handle['t'] = np.arange(4) * 0.05
        handle['obs_index'] = [2]
        handle['obs_max_time'] = [0.1]
        handle['run_mode'] = [0]
    notebook = (Path(__file__).resolve().parents[2] /
                'applications/icepack_model/examples/idealized_pig/read_results.ipynb')
    cells = json.loads(notebook.read_text())['cells']
    assert '_utility_imports' not in ''.join(cells[1]['source'])
    namespace = dict(os=os, h5py=h5py, np=np, results_dir=str(data),
                     extract_datasets_from_h5=extract_datasets_from_h5)
    try:
        exec(compile(''.join(cells[3]['source']), str(notebook), 'exec'), namespace)
        assert namespace['results_dir'] == str(data)
        assert isinstance(namespace['ensemble_vec_full'], h5py.Dataset)
        np.testing.assert_array_equal(namespace['ensemble_vec_full'][0, :, :], reference[0])
        np.testing.assert_array_equal(namespace['ensemble_vec_mean'][:], reference.mean(axis=1))
    finally:
        if '_ensemble_file' in namespace:
            namespace['_ensemble_file'].close()


def test_reference_mapping_refuses_fractional_indices(tmp_path):
    with pytest.raises(ValueError, match='must contain integers'):
        save_reference_layout(tmp_path, run_id='run', layout_id='mesh', row_maps={'h': [0.5, 1.5]})


def test_nonfinite_snapshot_is_not_published(history):
    data, _, _ = history
    source = data / '_mode3_state_history/steps/checkpoint_00000002/shard_00000000.h5'
    with h5py.File(source, 'a') as handle:
        handle['member_0/h'][0] = np.nan
    with pytest.raises(ValueError, match='nonfinite state'):
        export_ensemble_history(data)
    assert not (data / 'icesee_ensemble_data.h5').exists()
