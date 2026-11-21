#!/usr/bin/env python3
"""Inspect a closed ICESEE sample history in bounded row/time slabs.

Checks array integrity and stored means, not physical skill or observation timing.
Use only after the simulation and all HDF5 writers have closed their files.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import h5py
import numpy as np


def inspect(path: Path, members: int, columns: int, rows: int | None,
            rtol: float, atol: float) -> dict:
    path = path.resolve(strict=True)
    max_error = 0.0
    minimum, maximum = float('inf'), float('-inf')
    virtual_sources = []
    with h5py.File(path, 'r') as file:
        ensemble, mean = file['ensemble'], file['ensemble_mean']
        if ensemble.ndim != 3:
            raise ValueError(f'Expected state/member/time axes; found {ensemble.shape}')
        n, ne, nt = ensemble.shape
        if ne != members or nt != columns or (rows is not None and n != rows):
            raise ValueError(f'Unexpected ensemble shape {ensemble.shape}; expected ({rows or "any"},{members},{columns})')
        if n < 1 or mean.shape != (n, nt):
            raise ValueError(f'Unexpected mean shape {mean.shape}')
        # Missing virtual backing files can silently yield fill values. Refuse
        # that condition before interpreting finite data as a complete history.
        for dataset in (ensemble, mean):
            if not dataset.is_virtual:
                continue
            for source in dataset.virtual_sources():
                name = source.file_name.decode() if isinstance(source.file_name, bytes) else source.file_name
                name = name.rstrip('\0')
                source_path = path if name == '.' else Path(name)
                if not source_path.is_absolute():
                    source_path = path.parent / source_path
                source_path = source_path.resolve(strict=True)
                dset = source.dset_name.decode() if isinstance(source.dset_name, bytes) else source.dset_name
                if source_path == path:
                    if dset not in file:
                        raise ValueError(f'Missing virtual source dataset: {dset}')
                else:
                    with h5py.File(source_path, 'r') as backing:
                        if dset not in backing:
                            raise ValueError(f'Missing virtual source dataset: {source_path.name}:{dset}')
                virtual_sources.append(source_path.name)
        for column in range(nt):
            for first in range(0, n, 256):
                block = ensemble[first:first+256, :, column]
                stored = mean[first:first+256, column]
                if not np.isfinite(block).all() or not np.isfinite(stored).all():
                    raise ValueError(f'Nonfinite history at column {column}, row block {first}')
                computed = np.mean(block, axis=1)
                np.testing.assert_allclose(stored, computed, rtol=rtol, atol=atol,
                                           err_msg=f'Mean mismatch at column {column}, row block {first}')
                max_error = max(max_error, float(np.max(np.abs(stored-computed))))
                minimum = min(minimum, float(np.min(block)))
                maximum = max(maximum, float(np.max(block)))
    return {'status': 'pass', 'file': path.name, 'ensemble_shape': [n, ne, nt],
            'finite': True, 'mean_consistent': True, 'mean_max_absolute_error': max_error,
            'packed_state_min': minimum, 'packed_state_max': maximum,
            'virtual_backing_files_checked': sorted(set(virtual_sources)),
            'scope': 'Closed-file numerical integrity only; physical units differ between packed blocks. Verify physical time, analyses, geometry and scientific skill separately.'}


def inspect_native(directory: Path, members: int, columns: int, rows: int | None) -> dict:
    """Check the whole-member file-backed native policy's state snapshots."""
    directory = directory.resolve(strict=True)
    n = rows
    minimum, maximum = float('inf'), float('-inf')
    for column in range(columns):
        path = directory / f'native_analysis_{column:04d}.h5'
        with h5py.File(path, 'r') as file:
            state = file['states']
            if n is None:
                n = state.shape[0] if state.ndim == 2 else 0
            if n < 1 or state.shape != (n, members):
                raise ValueError(f'Unexpected native state shape {state.shape} at column {column}')
            for first in range(0, n, 256):
                block = state[first:first+256, :]
                if not np.isfinite(block).all():
                    raise ValueError(f'Nonfinite native state at column {column}, row block {first}')
                minimum = min(minimum, float(np.min(block)))
                maximum = max(maximum, float(np.max(block)))
    return {'status': 'pass', 'layout': 'whole-member native_analysis_NNNN.h5',
            'ensemble_shape': [n, members, columns], 'finite': True,
            'mean_consistency': 'not applicable; this format does not store ensemble_mean',
            'packed_state_min': minimum, 'packed_state_max': maximum,
            'scope': 'Snapshot integrity only; verify native history manifests/shards, physical time, analyses and scientific skill separately.'}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('history', type=Path)
    parser.add_argument('--native-analysis-series', action='store_true',
                        help='history is a directory of whole-member native_analysis_NNNN.h5 snapshots')
    parser.add_argument('--members', type=int, required=True)
    parser.add_argument('--columns', type=int, required=True)
    parser.add_argument('--rows', type=int)
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--atol', type=float, default=1e-10)
    args = parser.parse_args()
    if args.members < 1 or args.columns < 1 or (args.rows is not None and args.rows < 1):
        parser.error('Expected dimensions must be positive')
    result = (inspect_native(args.history, args.members, args.columns, args.rows)
              if args.native_analysis_series else
              inspect(args.history, args.members, args.columns, args.rows, args.rtol, args.atol))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
