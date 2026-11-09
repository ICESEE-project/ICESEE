import ast
from pathlib import Path


def test_partial_runner_has_no_unused_dense_process_covariance():
    source = (Path(__file__).parents[1]/"run_model_da/icesee_da_partial_parallel.py").read_text()
    tree = ast.parse(source)
    allocations = []
    for node in ast.walk(tree):
        if not isinstance(node,ast.Call) or not isinstance(node.func,ast.Attribute):
            continue
        if node.func.attr not in {"zeros","eye","empty","ones"} or not node.args:
            continue
        size = node.args[0]
        if node.func.attr == "eye" and isinstance(size,ast.Name) and size.id == "nd":
            allocations.append(node.lineno)
        if isinstance(size,ast.Tuple) and len(size.elts)==2 and all(
            isinstance(x,ast.Name) and x.id=="nd" for x in size.elts):
            allocations.append(node.lineno)
    assert not allocations, f"unused quadratic state allocation at {allocations}"
