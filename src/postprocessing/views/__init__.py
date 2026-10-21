"""Result views: the scientific figures and tables of the examples'
``read_results.ipynb`` notebooks, as functions.

A view is a :class:`ViewSpec` -- an id, a title, where it comes from in the
notebook, what it needs, its options -- plus a ``build`` function that turns
loaded outputs into a matplotlib figure or a table. Views never read files
themselves (that is :mod:`ICESEE.src.postprocessing.outputs`) and never save
anything (that is the caller's business), so the same function serves the
notebook, the launcher and the tests.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from ..outputs import DAOutputs, ResultsError


@dataclass(frozen=True)
class Option:
    """One user-adjustable view parameter. ``default`` may be a callable
    taking the :class:`ViewContext` (e.g. a grid point computed from the
    mesh dimensions, as the notebook does)."""

    name: str
    label: str
    kind: str                      # "int" | "float" | "choice" | "float_list"
    default: Any = None
    minimum: float | None = None
    maximum: Any = None            # number, or callable(ctx) -> number
    choices: tuple = ()
    help: str = ""


@dataclass
class ViewContext:
    outputs: DAOutputs
    params: dict | None = None
    data_dir: Path | None = None

    def section(self, name: str) -> dict:
        value = (self.params or {}).get(name)
        return value if isinstance(value, dict) else {}


@dataclass
class ViewResult:
    figure: Any = None             # matplotlib.figure.Figure
    table: dict | None = None      # {"columns": [...], "rows": [[...], ...], "caption": str}
    notes: list = field(default_factory=list)


@dataclass(frozen=True)
class ViewSpec:
    id: str
    title: str
    description: str
    source: str                    # where the notebook computes it
    build: Callable[[ViewContext, dict], ViewResult]
    kind: str = "figure"           # "figure" | "table"
    needs: tuple = ()              # "ensemble_members" | "dof_coordinates" | "params"
    options: tuple = ()
    rc: dict = field(default_factory=dict)
    supplementary: bool = False    # not part of the notebook; labelled as such

    def resolve_options(self, ctx: ViewContext, given: dict | None = None) -> dict:
        """Defaults filled in, user values coerced and range-checked."""
        given = dict(given or {})
        out = {}
        for opt in self.options:
            default = opt.default(ctx) if callable(opt.default) else opt.default
            value = given.pop(opt.name, None)
            if value is None or value == "":
                out[opt.name] = default
                continue
            out[opt.name] = _coerce(opt, value, ctx)
        if given:
            raise ResultsError("invalid", f"unknown option(s) for {self.id}: {', '.join(sorted(given))}")
        return out


def _coerce(opt: Option, value, ctx: ViewContext):
    try:
        if opt.kind == "int":
            value = int(value)
        elif opt.kind == "float":
            value = float(value)
        elif opt.kind == "float_list":
            if isinstance(value, str):
                value = [v for v in value.replace(",", " ").split() if v]
            value = [float(v) for v in value]
        elif opt.kind == "choice":
            if value not in opt.choices:
                raise ValueError(f"must be one of {', '.join(map(str, opt.choices))}")
    except (TypeError, ValueError) as exc:
        raise ResultsError("invalid", f"{opt.label}: {exc}") from exc
    maximum = opt.maximum(ctx) if callable(opt.maximum) else opt.maximum
    for item in (value if isinstance(value, list) else [value]):
        if opt.kind in ("int", "float", "float_list"):
            if opt.minimum is not None and item < opt.minimum:
                raise ResultsError("invalid", f"{opt.label} must be at least {opt.minimum}")
            if maximum is not None and item > maximum:
                raise ResultsError("invalid", f"{opt.label} must be at most {maximum}")
    return value


VIEWS: dict[str, ViewSpec] = {}
#: ids registered (or replaced) by a user's own view modules (load_user_views)
USER_VIEWS: set = set()
_LOADING_USER = [False]
USER_MODULE_PREFIX = "icesee_user_views"


def register(spec: ViewSpec) -> ViewSpec:
    VIEWS[spec.id] = spec
    if _LOADING_USER[0]:
        USER_VIEWS.add(spec.id)
    return spec


def get_view(view_id: str) -> ViewSpec:
    _load_builtin()
    try:
        return VIEWS[view_id]
    except KeyError:
        raise ResultsError("unsupported", f"no result view named {view_id!r}") from None


def all_views() -> dict[str, ViewSpec]:
    _load_builtin()
    return dict(VIEWS)


def _load_builtin() -> None:
    from . import common, flowline, icepack, issm, lorenz  # noqa: F401  (registration side effect)


def load_user_views(directory: Path | str) -> list[str]:
    """Import every ``*.py`` of ``directory`` (a user's copy of view modules)
    after the built-in views, so a view it registers under an existing id
    replaces the built-in one and a new id adds a view. Returns the ids the
    user's modules registered. The caller decides where this may run (the
    launcher runs it only inside a sandbox); this function only imports."""
    import importlib.util
    import sys

    _load_builtin()
    directory = Path(directory)
    _LOADING_USER[0] = True
    try:
        for path in sorted(directory.glob("*.py")):
            if path.name.startswith("_"):
                continue
            name = f"{USER_MODULE_PREFIX}.{path.stem}"
            spec = importlib.util.spec_from_file_location(name, path)
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
    finally:
        _LOADING_USER[0] = False
    return sorted(USER_VIEWS)
