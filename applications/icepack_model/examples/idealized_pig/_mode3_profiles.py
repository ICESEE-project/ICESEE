"""Bounded flowline diagnostics for native Idealized PIG members."""
from pathlib import Path

import h5py
import numpy as np


def profile_columns(step, dt, save_steps):
    """Map completed model steps to requested physical-year columns."""
    if step == 0:
        return [0]
    return [i + 1 for i, year in enumerate(save_steps)
            if np.isclose(step * dt, year, rtol=0.0, atol=1e-9)]


def sample_fields(fields, points):
    """Evaluate physical points without extrapolated overlap values."""
    meshes = {}
    for field in fields:
        mesh = field.function_space().mesh()
        meshes[id(mesh)] = (mesh, mesh.tolerance)
    try:
        values = [field.at(points, dont_raise=True, tolerance=1e-10)
                  for field in fields]
    finally:
        for mesh, tolerance in meshes.values():
            mesh.tolerance = tolerance
    return values


def write_member_profile(member_id, h, s, bed, *, step, comm, icesee_kwargs):
    """Sample collectively; only the spatial root opens the serial HDF5 file.

    Profiles contain at most 1000 points, never a gathered model state.
    Missing scheduled columns are NaN. Existing columns survive restart;
    only initialization (step zero) creates a fresh diagnostic file.
    """
    save_steps = tuple(icesee_kwargs.get('save_steps') or ())
    if not save_steps:
        return
    columns = profile_columns(step, float(icesee_kwargs['dt']), save_steps)
    if not columns:
        return
    start = np.array([-1587750.0, -200000.0])
    end = np.array([-1610250.0, -500000.0])
    fractions = np.linspace(0.0, 1.0, 1000)
    points = start + fractions[:, None] * (end - start)
    # Firedrake's vector point evaluation is collective and returns None
    # for points outside the distributed mesh when dont_raise is enabled.
    values = sample_fields((h, s, bed), points)
    valid = np.array([all(v[i] is not None for v in values) for i in range(len(points))])
    points = points[valid]
    arrays = [np.asarray([v[i] for i in np.flatnonzero(valid)], dtype=float).reshape(-1)
              for v in values]
    error = None
    if comm.Get_rank() == 0:
        try:
            directory = Path(icesee_kwargs.get('data_path') or '_modelrun_datasets')
            directory.mkdir(parents=True, exist_ok=True)
            path = directory / f'hs_ensemble_profiles{int(member_id)}.h5'
            mode = 'w' if step == 0 else 'a'
            with h5py.File(path, mode) as handle:
                shape = (len(points), len(save_steps) + 1)
                if 'h_profiles' not in handle:
                    for name in ('h_profiles', 's_profiles'):
                        handle.create_dataset(name, shape, dtype='f8', fillvalue=np.nan)
                    handle.create_dataset('bed_values', data=arrays[2])
                    handle.create_dataset('valid_points', data=points)
                    handle.create_dataset('distances', data=np.linalg.norm(points - start, axis=1))
                    handle.create_dataset('profile_times', data=(0.0, *save_steps))
                    handle.create_dataset('written_steps', data=np.full(shape[1], -1, dtype='i8'))
                if (handle['h_profiles'].shape != shape
                    or not np.allclose(handle['valid_points'][:], points)
                    or not np.allclose(handle['profile_times'][:], (0.0, *save_steps))):
                    raise ValueError(f'Profile geometry or save schedule changed: {path}')
                for column in columns:
                    handle['h_profiles'][:, column] = arrays[0]
                    handle['s_profiles'][:, column] = arrays[1]
                    handle['written_steps'][column] = step
        except Exception as exc:
            error = f'{type(exc).__name__}: {exc}'
    error = comm.bcast(error, root=0)
    if error:
        raise RuntimeError(f'Mode-3 flowline output failed: {error}')
