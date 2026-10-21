"""Rendering result views to files, and the command-line entry point.

::

    python -m ICESEE.src.postprocessing.render \\
        --data-dir RUN/_modelrun_datasets --contract EXAMPLE/results.yaml \\
        --params RUN/params.yaml --view lorenz.state_estimation \\
        --options '{"h_indx": 242}' --png out.png --json out.json

The JSON written next to the image always describes the outcome -- the
view's title, resolved options, notes, an optional table, or the error
``kind`` and message -- so a caller never has to parse a traceback. Exit
status: 0 rendered, 2 outputs missing/incomplete, 3 outputs corrupt,
4 invalid request (unknown view, bad option), 5 the user's own view
modules (``--user-views``) failed, 1 anything else.
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from dataclasses import dataclass, field
from pathlib import Path

from .contract import ResultsContract, load_contract
from .outputs import DOF_COORDINATES, ResultsError, load_da_outputs, load_params
from .views import USER_VIEWS, ViewContext, ViewSpec, get_view, load_user_views

EXIT = {"missing": 2, "incomplete": 2, "corrupt": 3, "invalid": 4, "unsupported": 4, "user_code": 5}


@dataclass
class Rendered:
    view: str
    title: str
    source: str
    options: dict
    notes: list = field(default_factory=list)
    table: dict | None = None
    png: Path | None = None
    supplementary: bool = False
    customized: bool = False

    def to_json(self) -> dict:
        return {
            "ok": True, "view": self.view, "title": self.title, "source": self.source,
            "options": self.options, "notes": self.notes, "table": self.table,
            "png": str(self.png) if self.png else None, "supplementary": self.supplementary,
            "customized": self.customized,
        }


def open_context(data_dir: Path | str, contract: ResultsContract, *, filter_type: str | None = None,
                 params: dict | None = None, results_dir: Path | str | None = None) -> ViewContext:
    """Load a run's outputs for the views of ``contract``."""
    data_dir = Path(data_dir)
    if filter_type is None:
        filter_type = str(((params or {}).get("enkf-parameters") or {}).get("filter_type") or "EnKF")
    outputs = load_da_outputs(data_dir, model_name=contract.model_name, filter_type=filter_type,
                              results_dir=results_dir)
    if contract.ensemble == "icesee_ensemble_data" and outputs.run_mode != 0:
        from .outputs import ENSEMBLE, _read

        path = data_dir / ENSEMBLE
        outputs.ensemble_mean = _read(path, ("ensemble_mean",))["ensemble_mean"]
        outputs.ensemble_file, outputs.ensemble_key = path, "ensemble"
    return ViewContext(outputs=outputs, params=params, data_dir=data_dir)


def availability(ctx: ViewContext, spec: ViewSpec) -> tuple[bool, str]:
    """Whether ``spec`` can be drawn for this run, and if not, why."""
    for need in spec.needs:
        if need == "params" and not ctx.params:
            return False, "needs the run's params.yaml"
        if need == "dof_coordinates" and not (Path(ctx.data_dir or ctx.outputs.data_dir) / DOF_COORDINATES).is_file():
            return False, ("needs the mesh coordinates (dof_coordinates.h5), exported by the "
                           "post-processing step in the Icepack environment")
    if "ensemble_members" in spec.needs and ctx.outputs.ensemble_file is None:
        return True, "ensemble members were not saved; only the mean is drawn"
    return True, ""


def render_view(ctx: ViewContext, view_id: str, options: dict | None = None, *,
                png: Path | None = None, dpi: int = 110):
    """Build one view and (for figures) save it. Returns ``(Rendered, figure)``."""
    import matplotlib

    spec = get_view(view_id)
    ok, why = availability(ctx, spec)
    if not ok:
        raise ResultsError("missing", f"{spec.title}: {why}")
    resolved = spec.resolve_options(ctx, options)
    with matplotlib.rc_context(spec.rc or {}):
        result = spec.build(ctx, resolved)
        if result.figure is not None and png is not None:
            png = Path(png)
            png.parent.mkdir(parents=True, exist_ok=True)
            tmp = png.with_name(png.name + ".part")
            result.figure.savefig(tmp, dpi=dpi, bbox_inches="tight", format="png")
            tmp.replace(png)
    notes = list(result.notes)
    if why:
        notes.append(why)
    rendered = Rendered(view=spec.id, title=spec.title, source=spec.source,
                        options=_jsonable(resolved), notes=notes, table=result.table,
                        png=png if result.figure is not None else None,
                        supplementary=spec.supplementary, customized=spec.id in USER_VIEWS)
    return rendered, result.figure


def describe(ctx: ViewContext, contract: ResultsContract) -> dict:
    """The contract's views for this run: availability, why not, and each
    option with its resolved default -- what a user interface needs to
    offer the views without loading anything itself."""
    views = []
    for view_id in contract.views:
        spec = get_view(view_id)
        ok, why = availability(ctx, spec)
        options = []
        for opt in spec.options:
            try:
                default = opt.default(ctx) if callable(opt.default) else opt.default
                maximum = opt.maximum(ctx) if callable(opt.maximum) else opt.maximum
            except ResultsError as exc:
                ok, why, default, maximum = False, exc.message, None, None
            options.append({"name": opt.name, "label": opt.label, "kind": opt.kind,
                            "default": _jsonable(default), "minimum": opt.minimum,
                            "maximum": _jsonable(maximum), "choices": list(opt.choices), "help": opt.help})
        views.append({"id": spec.id, "title": spec.title, "description": spec.description,
                      "source": spec.source, "kind": spec.kind, "available": ok, "reason": why,
                      "supplementary": spec.supplementary, "customized": spec.id in USER_VIEWS,
                      "options": options})
    o = ctx.outputs
    return {
        "ok": True, "views": views, "warnings": list(o.warnings),
        "run": {"model": o.model_name, "filter": o.filter_type, "time_steps": o.n_times,
                "analysis_steps": int(o.ensemble_mean.shape[1]), "state_size": o.n_state,
                "ensemble_size": o.n_ens, "observation_times": int(o.obs_index.size),
                "t_start": float(o.t[0]) if o.t.size else None,
                "t_end": float(o.t[-1]) if o.t.size else None, "run_mode": o.run_mode},
        "not_reproduced": list(contract.not_reproduced), "notebook": contract.notebook,
    }


def figure_data(fig) -> list[dict]:
    """Every plotted series of a figure, for comparisons: per axes, its
    labels and each line's / collection's data and style."""
    import numpy as np

    out = []
    for i, ax in enumerate(fig.axes):
        entry = {"axes": i, "title": ax.get_title(), "xlabel": ax.get_xlabel(),
                 "ylabel": ax.get_ylabel(), "lines": [], "collections": []}
        for line in ax.get_lines():
            entry["lines"].append({
                "label": line.get_label(), "x": np.asarray(line.get_xdata(), dtype=float),
                "y": np.asarray(line.get_ydata(), dtype=float), "linestyle": line.get_linestyle(),
                "marker": line.get_marker(), "color": line.get_color(),
            })
        for coll in ax.collections:
            arr = coll.get_array()
            entry["collections"].append({
                "type": type(coll).__name__,
                "array": None if arr is None else np.asarray(arr, dtype=float),
                "clim": coll.get_clim() if arr is not None else None,
                "cmap": getattr(coll.get_cmap(), "name", None) if arr is not None else None,
            })
        out.append(entry)
    return out


def _jsonable(value):
    import numpy as np

    if isinstance(value, dict):
        return {k: _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m ICESEE.src.postprocessing.render")
    p.add_argument("--data-dir", help="the directory holding true-wrong-<model>.h5")
    p.add_argument("--run-dir", help="a run directory to search for the data directory instead")
    p.add_argument("--contract", required=True, help="results.yaml, or the example directory")
    p.add_argument("--params", help="the run's params.yaml")
    p.add_argument("--results-dir")
    p.add_argument("--view", help="the view to draw (omit with --describe)")
    p.add_argument("--describe", action="store_true", help="list the views available for this run")
    p.add_argument("--options", default="{}")
    p.add_argument("--png")
    p.add_argument("--json", required=True)
    p.add_argument("--dpi", type=int, default=110)
    p.add_argument("--user-views", help="a directory of the user's own view modules, loaded after the built-in ones")
    args = p.parse_args(argv)

    out_json = Path(args.json)
    out_json.parent.mkdir(parents=True, exist_ok=True)

    def fail(kind: str, message: str, code: int) -> int:
        out_json.write_text(json.dumps({"ok": False, "kind": kind, "message": message, "view": args.view}))
        return code

    try:
        import logging

        import matplotlib

        matplotlib.use("Agg")
        logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
        if args.user_views:
            _USER_FILES.extend(str(p) for p in Path(args.user_views).glob("*.py"))
            try:
                load_user_views(args.user_views)
            except Exception as exc:
                return fail("user_code", "your plotting code could not be loaded: " + _user_error(exc), 5)
        contract_path = Path(args.contract)
        contract = load_contract(contract_path if contract_path.is_dir() else contract_path.parent)
        if contract is None:
            return fail("unsupported", "this example declares no result views (no results.yaml)", 4)
        if not contract.supported:
            return fail("unsupported", contract.reason or "this example has no result views", 4)
        params = load_params(args.params) if args.params and Path(args.params).is_file() else None
        if not args.data_dir:
            if not args.run_dir:
                return fail("invalid", "give --data-dir or --run-dir", 4)
            from .outputs import TRUE_WRONG, find_data_dir

            found = find_data_dir(args.run_dir, model_name=contract.model_name, data_dir_name=contract.data_dir)
            if found is None:
                return fail("no_outputs", f"no {TRUE_WRONG.format(model=contract.model_name)} in the run", 2)
            args.data_dir = str(found)
        if args.describe:
            ctx = open_context(args.data_dir, contract, params=params, results_dir=args.results_dir)
            out_json.write_text(json.dumps({**describe(ctx, contract), "data_dir": str(args.data_dir)}))
            return 0
        if args.view not in contract.views:
            return fail("invalid", f"{args.view} is not a view of this example", 4)
        options = json.loads(args.options or "{}")
        ctx = open_context(args.data_dir, contract, params=params, results_dir=args.results_dir)
        rendered, _fig = render_view(ctx, args.view, options, png=Path(args.png) if args.png else None,
                                     dpi=args.dpi)
        out_json.write_text(json.dumps(rendered.to_json()))
        return 0
    except ResultsError as exc:
        return fail(exc.kind, exc.message, EXIT.get(exc.kind, 1))
    except json.JSONDecodeError as exc:
        return fail("invalid", f"options are not valid JSON ({exc})", 4)
    except Exception as exc:  # never a bare traceback to the caller
        traceback.print_exc(file=sys.stderr)
        if args.user_views and _user_frames(exc):
            return fail("user_code", "your plotting code failed: " + _user_error(exc), 5)
        return fail("error", f"{type(exc).__name__}: {exc}", 1)


def _user_frames(exc: BaseException) -> list:
    frames = traceback.extract_tb(exc.__traceback__)
    return [f for f in frames if f.filename in _USER_FILES]


_USER_FILES: list = []


def _user_error(exc: BaseException) -> str:
    """``file line N: Error: message`` pointing at the user's own code."""
    frames = _user_frames(exc)
    if isinstance(exc, SyntaxError) and exc.filename:
        return f"{Path(exc.filename).name} line {exc.lineno}: SyntaxError: {exc.msg}"
    where = f"{Path(frames[-1].filename).name} line {frames[-1].lineno}: " if frames else ""
    return f"{where}{type(exc).__name__}: {exc}"[:600]


if __name__ == "__main__":
    sys.exit(main())
