"""Flowline I/O contracts without running the ice solver."""
import importlib.util
from pathlib import Path

import h5py
import numpy as np
import pytest

_path = Path(__file__).parents[2] / 'applications/icepack_model/examples/idealized_pig/_mode3_profiles.py'
_spec = importlib.util.spec_from_file_location('pig_mode3_profiles_test', _path)
profiles = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(profiles)


class Comm:
    def Get_rank(self):
        return 0

    def bcast(self, value, root=0):
        return value


class Field:
    def __init__(self, value):
        self.value = value
        self.calls = 0

    def function_space(self):
        from types import SimpleNamespace
        return SimpleNamespace(mesh=lambda: self)

    tolerance = 0.5

    def at(self, points, dont_raise=False, tolerance=None):
        assert dont_raise
        assert tolerance == 1e-10
        self.calls += 1
        return [self.value if i % 2 == 0 else None for i in range(len(points))]


def write(tmp_path, step, value=10, member_id=7, save_steps=(0.3, 0.5)):
    fields = [Field(value), Field(value + 1), Field(-100)]
    profiles.write_member_profile(member_id, *fields, step=step, comm=Comm(),
        icesee_kwargs={'dt': 0.1, 'save_steps': save_steps, 'data_path': str(tmp_path)})
    return fields


def test_initial_profile_and_scheduled_resume_preserve_columns(tmp_path):
    write(tmp_path, 0)
    write(tmp_path, 3, 20)
    write(tmp_path, 5, 30)
    path = tmp_path / 'hs_ensemble_profiles7.h5'
    with h5py.File(path) as f:
        assert f['h_profiles'].shape == (500, 3)
        np.testing.assert_allclose(f['h_profiles'][0], [10, 20, 30])
        np.testing.assert_array_equal(f['written_steps'][:], [0, 3, 5])
        np.testing.assert_allclose(f['profile_times'][:], [0, 0.3, 0.5])
    write(tmp_path, 3, 40)
    with h5py.File(path) as f:
        np.testing.assert_allclose(f['h_profiles'][0], [10, 40, 30])


def test_no_sampling_or_file_open_between_save_steps(tmp_path):
    assert all(f.calls == 0 for f in write(tmp_path, 1))
    assert not list(tmp_path.iterdir())


def test_resume_missing_file_does_not_invent_history(tmp_path):
    write(tmp_path, 5)
    with h5py.File(tmp_path / 'hs_ensemble_profiles7.h5') as f:
        assert np.isnan(f['h_profiles'][:, :2]).all()
        np.testing.assert_array_equal(f['written_steps'][:], [-1, -1, 5])


def test_members_have_separate_files(tmp_path):
    write(tmp_path, 0, member_id=2)
    write(tmp_path, 0, member_id=5)
    assert {p.name for p in tmp_path.iterdir()} == {'hs_ensemble_profiles2.h5', 'hs_ensemble_profiles5.h5'}


def test_resume_changed_schedule_is_rejected(tmp_path):
    write(tmp_path, 0)
    with pytest.raises(RuntimeError, match='schedule changed'):
        write(tmp_path, 3, save_steps=(0.3, 0.8))


def test_true_profile_file_includes_final_step_and_marks_future_years(tmp_path):
    """Execute the actual truth-generator body with a tiny deterministic solver."""
    import ast
    import os
    from types import SimpleNamespace
    source = _path.with_name('_icepack_enkf.py').read_text()
    function = next(node for node in ast.parse(source).body
                    if isinstance(node, ast.FunctionDef) and node.name == 'generate_true_state')
    class FakeField:
        def __init__(self, data):
            self.dat = SimpleNamespace(data_ro=np.asarray(data))
        def copy(self, deepcopy=False):
            return FakeField(self.dat.data_ro.copy())
    h, u, s, melt = FakeField([10., 11.]), FakeField([[1., 2.], [3., 4.]]), FakeField([20., 21.]), FakeField([0., 0.])
    mapping = {name: np.arange(i * 2, i * 2 + 2)
               for i, name in enumerate(('h', 'u', 'v', 's', 'basal_melt_field'))}
    env = {'np': np, 'os': os, 'h5py': h5py,
           'profile_columns': profiles.profile_columns,
           'icesee_get_index': lambda **kw: (None, mapping, None),
           'initial_flowline_profile': lambda kw: (np.array([[10., 11.]]), np.array([[20., 21.]]), np.zeros((2, 2)), [0., 1.], [-1., -2.]),
           'flowline_profile': lambda h, s, points: (h.dat.data_ro, s.dat.data_ro),
           'BasalMeltRate': lambda *args, **kw: (melt, 0.),
           'Icepack': lambda *args: (FakeField([12., 13.]), u, s, None, None),
           'EXPERIMENT_TRUE': 'true'}
    exec(compile(ast.Module(body=[function], type_ignores=[]), '<truth generator>', 'exec'), env)
    trajectory = np.zeros((10, 2))
    env['generate_true_state'](dt=0.1, nt=1, h0=h, u=u, s0=s, basal_melt_field=melt,
        statevec_true=trajectory, save_steps=[0.1, 34.], data_path=str(tmp_path), vec_inputs=list(mapping))
    with h5py.File(tmp_path / 'hs_profiles_true.h5') as f:
        np.testing.assert_allclose(f['h_profiles'][:, 1], [12., 13.])
        assert np.isnan(f['h_profiles'][:, 2]).all()
        np.testing.assert_array_equal(f['written_steps'][:], [0, 1, -1])
    np.testing.assert_allclose(trajectory[mapping['h'], 1], [12., 13.])
