"""The results contract: an example's declaration of where its outputs are
and which result views read them.

Each example that supports post-processing ships a ``results.yaml`` next to
its ``read_results.ipynb``::

    schema: icesee-results/1
    model_name: lorenz            # the <model> in true-wrong-<model>.h5
    notebook: read_results.ipynb
    outputs:
      data_dir: _modelrun_datasets
      results_dir: _modelrun_datasets   # filter outputs when run_mode != 0
      ensemble: auto                    # or icesee_ensemble_data
    views:
      - lorenz.state_estimation
      - da.error_and_spread
    not_reproduced:                     # notebook cells deliberately left out
      - cell: 9
        reason: ...

An example may also declare itself unsupported (``status: unsupported`` with
a ``reason``) so a user interface can say why there are no results views
instead of guessing.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

SCHEMA = "icesee-results/1"
CONTRACT_FILE = "results.yaml"


class ContractError(Exception):
    pass


@dataclass(frozen=True)
class ResultsContract:
    model_name: str
    views: tuple = ()
    data_dir: str = "_modelrun_datasets"
    results_dir: str = "_modelrun_datasets"
    ensemble: str = "auto"
    notebook: str = "read_results.ipynb"
    status: str = "supported"
    reason: str = ""
    not_reproduced: tuple = field(default_factory=tuple)
    path: Path | None = None

    @property
    def supported(self) -> bool:
        return self.status == "supported" and bool(self.views)


def parse_contract(data: dict, *, path: Path | None = None) -> ResultsContract:
    """Validate a contract mapping; every problem is reported at once."""
    from .views import all_views

    errors = []
    if not isinstance(data, dict):
        raise ContractError(f"{path or CONTRACT_FILE}: not a mapping")
    if data.get("schema") != SCHEMA:
        errors.append(f"schema must be {SCHEMA!r}")
    model = data.get("model_name")
    if not isinstance(model, str) or not model.strip():
        errors.append("model_name is required")
    status = data.get("status", "supported")
    if status not in ("supported", "unsupported"):
        errors.append("status must be 'supported' or 'unsupported'")
    outputs = data.get("outputs") or {}
    if not isinstance(outputs, dict):
        errors.append("outputs must be a mapping")
        outputs = {}
    for key in ("data_dir", "results_dir"):
        value = outputs.get(key, "_modelrun_datasets")
        if not isinstance(value, str) or value.startswith("/") or ".." in Path(value).parts:
            errors.append(f"outputs.{key} must be a relative path inside the example")
    ensemble = outputs.get("ensemble", "auto")
    if ensemble not in ("auto", "icesee_ensemble_data"):
        errors.append("outputs.ensemble must be 'auto' or 'icesee_ensemble_data'")
    views = data.get("views") or []
    if not isinstance(views, list) or not all(isinstance(v, str) for v in views):
        errors.append("views must be a list of view ids")
        views = []
    known = all_views()
    unknown = [v for v in views if v not in known]
    if unknown:
        errors.append(f"unknown view id(s): {', '.join(unknown)}")
    if status == "supported" and not views:
        errors.append("a supported contract must list at least one view")
    if status == "unsupported" and not str(data.get("reason") or "").strip():
        errors.append("an unsupported contract must give a reason")
    not_reproduced = data.get("not_reproduced") or []
    if not isinstance(not_reproduced, list):
        errors.append("not_reproduced must be a list")
        not_reproduced = []
    if errors:
        raise ContractError(f"{path or CONTRACT_FILE}: " + "; ".join(errors))
    return ResultsContract(
        model_name=model.strip(), views=tuple(views),
        data_dir=outputs.get("data_dir", "_modelrun_datasets"),
        results_dir=outputs.get("results_dir", outputs.get("data_dir", "_modelrun_datasets")),
        ensemble=ensemble, notebook=str(data.get("notebook") or "read_results.ipynb"),
        status=status, reason=str(data.get("reason") or ""),
        not_reproduced=tuple(dict(item) for item in not_reproduced if isinstance(item, dict)),
        path=path,
    )


def load_contract(example_dir: Path | str) -> ResultsContract | None:
    """The contract in ``example_dir``, or ``None`` when the example has none."""
    import yaml

    path = Path(example_dir) / CONTRACT_FILE
    if not path.is_file():
        return None
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ContractError(f"{path}: not valid YAML ({exc})") from exc
    return parse_contract(data, path=path)
