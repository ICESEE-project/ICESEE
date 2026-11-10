"""Retain the basal driver's requested inversion across an earlier deferral.

Execute its actual configuration update and the runner's actual activation
statements without importing the application or starting MATLAB.
"""
from pathlib import Path
import ast
import types
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]


def _driver_update(kwargs):
    p = ROOT / "applications/issm_model/examples/basal_friction_variation/run_da_issm.py"
    tree = ast.parse(p.read_text())
    updates = [node for node in ast.walk(tree) if isinstance(node, ast.Expr)
               and isinstance(node.value, ast.Call)
               and isinstance(node.value.func, ast.Attribute)
               and node.value.func.attr == "update"
               and node.value.args and isinstance(node.value.args[0], ast.Dict)
               and any(isinstance(key, ast.Constant) and key.value == "inversion_flag"
                       for key in node.value.args[0].keys)]
    assert len(updates) == 1
    namespace = dict(icesee_kwargs=kwargs, np=np,
                     socket=types.SimpleNamespace(gethostname=lambda: "testhost"),
                     icesee_cwd="application", issm_dir="installation", issm_examples_dir="work",
                     icesee_rank=0, icesee_size=1, ens_id=0)
    exec(compile(ast.fix_missing_locations(ast.Module(body=updates, type_ignores=[])),
                 str(p), "exec"), namespace)


def _activate_from_actual_runner(kwargs, time):
    p = ROOT / "src/run_model_da/icesee_da_partial_parallel.py"
    tree = ast.parse(p.read_text())
    names = {"inversion_enabled", "inversion_start_time", "inversion_flag"}
    statements = []
    for name in ["inversion_enabled", "inversion_start_time", "inversion_flag"]:
        matches = [node for node in ast.walk(tree) if isinstance(node, ast.Assign)
                   and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)
                   and node.targets[0].id == name]
        assert len(matches) == 1
        statements.extend(matches)
    namespace = dict(icesee_kwargs=kwargs, cycle_time=time)
    exec(compile(ast.fix_missing_locations(ast.Module(body=statements, type_ignores=[])),
                 str(p), "exec"), namespace)
    kwargs["inversion_flag"] = namespace["inversion_flag"]
    return namespace["inversion_flag"]


@pytest.mark.parametrize("requested,expected", [(True, [False, True, True]),
                                              (False, [False, False, False])])
def test_actual_driver_and_runner_preserve_delayed_inversion(requested, expected):
    kwargs = dict(Lx=640e3, Ly=80e3, nx=4, ny=2, timesteps_per_year=.2,
                  num_years=2., tinitial=0., steps=55, Nens=4,
                  inversion_flag=requested, inversion_start_time=2.)
    _driver_update(kwargs)
    actual = [_activate_from_actual_runner(kwargs, t) for t in [1., 2., 3.]]
    assert actual == expected
    assert kwargs.get("inversion_enabled", False) is requested
    assert kwargs["inversion_start_time"] == 2.
