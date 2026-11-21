"""Exercise the real basal callbacks with HDF5 and a stub external solver.

Only callback definitions are compiled, avoiding side-effecting application
configuration imports. No MATLAB/ISSM solve is performed by these unit tests.
"""
from pathlib import Path
import ast
import os
import types

import h5py
import numpy as np
import pytest


def _callbacks():
    source = (Path(__file__).resolve().parents[2] / "applications" /
              "issm_model/examples/basal_friction_variation/_issm_enkf.py")
    names = {"generate_true_state", "generate_nurged_state"}
    tree = ast.parse(source.read_text())
    selected = ast.Module(body=[node for node in tree.body
                                if isinstance(node, ast.FunctionDef) and node.name in names],
                          type_ignores=[])
    namespace = {"os": os, "np": np, "h5py": h5py}
    exec(compile(ast.fix_missing_locations(selected), str(source), "exec"), namespace)
    return namespace


@pytest.mark.parametrize("callback,state_key,file_stem", [
    ("generate_true_state", "statevec_true", "ensemble_true_state"),
    ("generate_nurged_state", "statevec_nurged", "ensemble_nurged_state"),
])
@pytest.mark.parametrize("disk_backed", [True, False])
@pytest.mark.parametrize("initial_state_only", [True, False])
def test_generation_restores_application_cwd_and_preserves_result(
        tmp_path, monkeypatch, callback, state_key, file_stem, disk_backed, initial_state_only):
    app = tmp_path / "application"
    work = tmp_path / "issm_example"
    data = app / "relative_data"
    data.mkdir(parents=True)
    work.mkdir()
    initial = np.array([1.0, 2.0])
    terminal = np.array([10.0, 20.0])
    with h5py.File(data / f"{file_stem}_0.h5", "w") as f:
        f["Thickness_0"] = initial[None, :]
        f["Thickness_1"] = terminal[None, :]
    with h5py.File(data / "mesh_idxy_0.h5", "w") as f:
        f["fric_x"] = [0.0, 1.0]
        f["fric_y"] = [0.0, 1.0]

    namespace = _callbacks()
    def external_solver(**kwargs):
        assert Path.cwd() == work
        expected_fname = ("initial_" if initial_state_only else "") + (
            "true_state.mat" if callback == "generate_true_state" else "nurged_state.mat")
        assert kwargs["fname"] == expected_fname
    namespace["ISSM_model"] = external_solver
    namespace["icesee_get_index"] = lambda **kw: (
        [], {"Thickness": np.arange(2)}, None)
    namespace["gs"] = types.SimpleNamespace(
        Gaussian=lambda **kw: object(),
        SRF=lambda model, seed: lambda xy: np.zeros(2),
    )
    monkeypatch.chdir(app)
    kwargs = dict(t=[0.0, 0.2], server=None, issm_examples_dir=str(work),
                  icesee_path=str(app), data_path="relative_data",
                  comm=types.SimpleNamespace(Get_rank=lambda: 0),
                  vec_inputs=["Thickness"], ens_id=0, nt=1, nd=2,
                  total_state_param_vars=1, range_friction=1.0,
                  sill_friction=1.0, nugget_friction=0.0,
                  initial_state_only=initial_state_only)
    with h5py.File(data / "target.h5", "w") as target:
        columns = 1 if initial_state_only else 2
        state = (target.create_dataset("state", shape=(2, columns))
                 if disk_backed else np.zeros((2, columns)))
        kwargs[state_key] = state
        result = namespace[callback](**kwargs)
        assert Path.cwd() == app
        assert Path("relative_data/target.h5").is_file()
        np.testing.assert_array_equal(state[:, 0], initial)
        if not initial_state_only:
            np.testing.assert_array_equal(state[:, 1], terminal)
        if disk_backed:
            assert result is None
        else:
            np.testing.assert_array_equal(result["Thickness"][:, 0], initial)
            if not initial_state_only:
                np.testing.assert_array_equal(result["Thickness"][:, 1], terminal)
