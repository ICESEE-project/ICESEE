"""Materialize committed block checkpoints as a standalone plotting history.

Checkpoint state remains distributed during the run. Finalization reads one
reference-row chunk and member at a time, after every source file is closed.
The output contains no external links and never replaces restart checkpoints.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path
import tempfile

import h5py
import numpy as np


def _temporary(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f'.{path.name}.', suffix='.tmp', dir=path.parent)
    os.close(fd)
    return Path(name)


def _publish(temporary, target):
    with temporary.open('rb') as stream:
        os.fsync(stream.fileno())
    os.replace(temporary, target)
    fd = os.open(target.parent, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def save_reference_layout(history_root, *, run_id, layout_id, row_maps):
    """Persist reference-row -> native-field-row permutations for finalization.

    An existing mapping must agree, so resume cannot silently relabel retained
    snapshots. Field sizes and permutation coverage are checked again on export.
    """
    path = Path(history_root) / 'reference_layout.h5'
    maps = {}
    if not row_maps:
        raise ValueError('reference layout must describe at least one state field')
    for name, value in row_maps.items():
        mapping = np.asarray(value)
        if mapping.dtype.kind not in 'iu':
            raise ValueError('reference row permutations must contain integers')
        maps[str(name)] = mapping.astype(np.int64, copy=False)
    for name, mapping in maps.items():
        if mapping.ndim != 1 or not np.array_equal(np.sort(mapping), np.arange(mapping.size)):
            raise ValueError(f'{name}: reference layout must be a complete permutation')
    if path.exists():
        with h5py.File(path, 'r') as handle:
            if handle.attrs['run_id'] != run_id or handle.attrs['layout_id'] != layout_id:
                raise ValueError('reference layout belongs to another run/layout')
            if set(handle['fields']) != set(maps):
                raise ValueError('reference layout fields changed on resume')
            for name, mapping in maps.items():
                if not np.array_equal(handle[f'fields/{name}'][:], mapping):
                    raise ValueError('reference row numbering changed on resume')
        return path
    temporary = _temporary(path)
    try:
        with h5py.File(temporary, 'w') as handle:
            handle.attrs.update(run_id=run_id, layout_id=layout_id,
                                row_order='reference_to_native')
            for name, mapping in maps.items():
                handle.create_dataset(f'fields/{name}', data=mapping)
        _publish(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


def _inventory(history_root):
    final_paths = sorted((history_root / 'steps').glob('checkpoint_*/manifest.json'))
    complete = []
    for path in final_paths:
        manifest = json.loads(path.read_text())
        if manifest.get('complete'):
            complete.append((int(manifest['timestep']), path, manifest))
    if not complete:
        raise ValueError('no committed step checkpoints to export')
    _, _, final = max(complete, key=lambda item: item[0])
    nt = int(final['metadata']['nt'])
    if final['timestep'] != nt - 1 or final['metadata']['next_step'] != nt:
        raise ValueError('run is not complete; aggregate history was not published')
    paths = [history_root / 'initial_condition/checkpoint_00000000/manifest.json']
    paths += [history_root / 'steps' / f'checkpoint_{step:08d}/manifest.json'
              for step in range(nt)]
    snapshots = []
    signature = {key: final[key] for key in
                 ('run_id', 'layout_id', 'global_size', 'number_of_members', 'dtype', 'blocks')}
    blocks = signature['blocks']
    names, offset = set(), 0
    for block in blocks:
        name, size = block['name'], int(block['global_size'])
        if name in names or not name or '/' in name or size <= 0 or block['global_offset'] != offset:
            raise ValueError('checkpoint blocks must uniquely cover the packed state')
        names.add(name)
        offset += size
    if offset != int(signature['global_size']):
        raise ValueError('checkpoint block sizes disagree with global state size')
    nens = int(signature['number_of_members'])
    if nens < 1:
        raise ValueError('checkpoint ensemble must contain members')
    times = []
    for column, path in enumerate(paths):
        if not path.is_file():
            raise ValueError(f'missing history snapshot: {path.parent}; retain every timestep before running')
        manifest = json.loads(path.read_text())
        if not manifest.get('complete') or any(manifest.get(key) != value for key, value in signature.items()):
            raise ValueError(f'incomplete or incompatible checkpoint: {path.parent}')
        metadata = manifest['metadata']
        if any(metadata.get(key) != final['metadata'].get(key)
               for key in ('schedule_digest', 'error_mode')):
            raise ValueError(f'checkpoint scientific schedule differs: {path.parent}')
        if (int(metadata['nt']) != nt or int(metadata['next_step']) != column
                or int(metadata['completed_step']) != column - 1
                or (column and int(manifest['timestep']) != column - 1)):
            raise ValueError(f'checkpoint time indices disagree: {path.parent}')
        times.append(float(metadata['model_time']))
        coverage = {(member, block['name']): [] for member in range(nens) for block in blocks}
        sources = []
        for shard in manifest['shards']:
            source = (path.parent / shard['file']).resolve()
            if source.parent != path.parent.resolve() or not source.is_file():
                raise ValueError(f'missing or invalid shard path: {source}')
            sources.append((source, shard))
            for member in shard['member_ids']:
                for block in shard['blocks']:
                    key = (int(member), block['name'])
                    if key not in coverage:
                        raise ValueError('unexpected checkpoint member or field')
                    definition = next(item for item in blocks if item['name'] == block['name'])
                    if any(block[key] != definition[key] for key in ('global_offset', 'global_size')):
                        raise ValueError('shard block definition differs from manifest')
                    coverage[key].append((int(block['owned_start']), int(block['owned_stop'])))
        for (member, name), intervals in coverage.items():
            end = 0
            for start, stop in sorted(intervals):
                if start != end or stop < start:
                    raise ValueError(f'checkpoint ownership gap/overlap: member {member}, {name}')
                end = stop
            if end != next(item['global_size'] for item in blocks if item['name'] == name):
                raise ValueError(f'incomplete checkpoint ownership: member {member}, {name}')
        snapshots.append((manifest, sources))
    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError('checkpoint model times must be finite and increasing')
    return signature, snapshots, np.asarray(times)


def prepare_reference_layout(world, history_root, *, run_id, layout_id, row_maps):
    """Check group numbering agreement and persist the plotting permutation."""
    digest = hashlib.sha256()
    for name, mapping in sorted(row_maps.items()):
        digest.update(name.encode())
        digest.update(np.asarray(mapping, dtype='<i8').tobytes())
    identities = world.allgather(digest.hexdigest())
    if len(set(identities)) != 1:
        raise ValueError('model groups have different reference/native row mappings')
    error = None
    if int(world.Get_rank()) == 0:
        try:
            save_reference_layout(history_root, run_id=run_id, layout_id=layout_id,
                                  row_maps=row_maps)
        except Exception as exc:
            error = f'{type(exc).__name__}: {exc}'
    error = world.bcast(error, root=0)
    if error is not None:
        raise RuntimeError(f'cannot save reference row mapping: {error}')


def export_ensemble_history(data_path, *, row_chunk_size=8192):
    """Atomically create ``icesee_ensemble_data.h5`` from a complete native run.

    ``ensemble`` has (reference state, member, time) shape; ``ensemble_mean``
    has (reference state, time) shape. Column zero is the initial ensemble.
    Reference permutations are mandatory: native PETSc row numbers need not
    identify the same physical nodes as truth/observation row numbers.
    """
    data_path = Path(data_path)
    history_root = data_path / '_mode3_state_history'
    row_chunk_size = int(row_chunk_size)
    if row_chunk_size < 1:
        raise ValueError('row_chunk_size must be positive')
    signature, snapshots, times = _inventory(history_root)
    blocks = signature['blocks']
    maps = {}
    with h5py.File(history_root / 'reference_layout.h5', 'r') as handle:
        if any(handle.attrs[key] != signature[key] for key in ('run_id', 'layout_id')):
            raise ValueError('reference row mapping belongs to another run/layout')
        if set(handle['fields']) != {block['name'] for block in blocks}:
            raise ValueError('reference row mapping fields do not match checkpoints')
        for block in blocks:
            mapping = handle[f"fields/{block['name']}"][:]
            size = int(block['global_size'])
            if (mapping.dtype.kind not in 'iu' or mapping.shape != (size,)
                    or not np.array_equal(np.sort(mapping), np.arange(size))):
                raise ValueError('invalid reference row permutation')
            maps[block['name']] = mapping
    nd, nens = int(signature['global_size']), int(signature['number_of_members'])
    dtype = np.dtype(signature['dtype'])
    if dtype.kind != 'f':
        raise ValueError('plotting history requires floating point state')
    target = data_path / 'icesee_ensemble_data.h5'
    temporary = _temporary(target)
    try:
        with h5py.File(temporary, 'w') as output:
            output.attrs.update(format='icesee-ensemble-history-v1', execution_mode=3,
                                state_order='serial_reference', run_id=signature['run_id'],
                                layout_id=signature['layout_id'], complete=False,
                                history_stage='post_cycle_analysis_when_scheduled',
                                blocks=json.dumps(blocks))
            ensemble = output.create_dataset('ensemble', (nd, nens, len(times)), dtype=dtype,
                                             chunks=(min(row_chunk_size, nd), 1, 1))
            mean = output.create_dataset('ensemble_mean', (nd, len(times)), dtype=dtype,
                                         chunks=(min(row_chunk_size, nd), 1))
            output.create_dataset('t', data=times)
            output.create_dataset('checkpoint_steps', data=np.arange(-1, len(times) - 1))
            output.create_dataset('analysis_applied', data=np.array([
                bool(manifest['metadata'].get('did_analysis', False))
                for manifest, _ in snapshots], dtype=bool))
            for column, (_, sources) in enumerate(snapshots):
                with ExitStack() as stack:
                    handles = [(stack.enter_context(h5py.File(path, 'r')), shard) for path, shard in sources]
                    for block in blocks:
                        name, size, offset = block['name'], int(block['global_size']), int(block['global_offset'])
                        for start in range(0, size, row_chunk_size):
                            stop = min(start + row_chunk_size, size)
                            native_rows = maps[name][start:stop]
                            total = np.zeros(stop - start, dtype=np.float64)
                            for member in range(nens):
                                values = np.empty(stop - start, dtype=dtype)
                                for handle, shard in handles:
                                    if member not in shard['member_ids']:
                                        continue
                                    for owned in shard['blocks']:
                                        if owned['name'] != name:
                                            continue
                                        first, last = int(owned['owned_start']), int(owned['owned_stop'])
                                        positions = np.flatnonzero((native_rows >= first) & (native_rows < last))
                                        if not positions.size:
                                            continue
                                        local = native_rows[positions] - first
                                        order = np.argsort(local)
                                        dataset = handle[f'member_{member}/{name}']
                                        if dataset.shape != (last - first,) or dataset.dtype != dtype:
                                            raise ValueError('shard dataset shape/dtype disagrees with manifest')
                                        values[positions[order]] = dataset[local[order]]
                                if not np.isfinite(values).all():
                                    raise ValueError('checkpoint contains nonfinite state')
                                ensemble[offset + start:offset + stop, member, column] = values
                                total += values
                            mean[offset + start:offset + stop, column] = total / nens
                print(f'[ICESEE] ensemble history: {column + 1}/{len(times)} snapshots', flush=True)
            output.attrs['complete'] = True
            output.flush()
        _publish(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    return target


def finalize_ensemble_history(world, data_path, *, row_chunk_size=8192):
    """Export on world root and broadcast failure before any rank proceeds."""
    result, error = None, None
    if int(world.Get_rank()) == 0:
        try:
            result = str(export_ensemble_history(data_path, row_chunk_size=row_chunk_size))
        except Exception as exc:
            error = f'{type(exc).__name__}: {exc}'
    result, error = world.bcast((result, error), root=0)
    if error is not None:
        raise RuntimeError(f'cannot finalize native ensemble history: {error}')
    return Path(result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data_path', type=Path)
    parser.add_argument('--row-chunk-size', type=int, default=8192)
    args = parser.parse_args()
    print(export_ensemble_history(args.data_path, row_chunk_size=args.row_chunk_size))
