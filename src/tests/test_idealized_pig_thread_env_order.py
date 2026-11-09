"""The Idealized PIG drivers pin one BLAS thread per rank before NumPy loads.

OpenBLAS reads OMP_NUM_THREADS once, when NumPy is imported. run_da_icepack.py
used to import NumPy first and set the variable afterwards, so under a bare
``mpirun -np 8 python run_da_icepack.py`` (no launcher exporting it) every
rank's NumPy BLAS still used all cores.
"""

import ast
from pathlib import Path

import pytest

PIG = Path(__file__).resolve().parents[2] / "applications/icepack_model/examples/idealized_pig"


def _first_line(tree, predicate):
    lines = [node.lineno for node in tree.body if predicate(node)]
    return min(lines) if lines else None


def _imports_numpy(node):
    if isinstance(node, ast.Import):
        return any(alias.name.split(".")[0] == "numpy" for alias in node.names)
    if isinstance(node, ast.ImportFrom):
        return (node.module or "").split(".")[0] == "numpy"
    return False


def _sets_omp_threads(node):
    if not isinstance(node, ast.Assign):
        return False
    target = node.targets[0]
    return (
        isinstance(target, ast.Subscript)
        and ast.unparse(target.value) == "os.environ"
        and ast.unparse(target.slice).strip("'\"") == "OMP_NUM_THREADS"
    )


@pytest.mark.parametrize(
    "script", ["run_da_icepack.py", "tools/profile_reference_trajectory.py"]
)
def test_omp_threads_set_before_numpy_import(script):
    tree = ast.parse((PIG / script).read_text())
    set_line = _first_line(tree, _sets_omp_threads)
    numpy_line = _first_line(tree, _imports_numpy)
    assert set_line is not None and numpy_line is not None
    assert set_line < numpy_line
