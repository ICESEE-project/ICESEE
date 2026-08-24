# ICESEE/config/_compat_metadata.py
from __future__ import annotations
from contextlib import nullcontext
from typing import Iterable, Optional

try:
    # stdlib (Py>=3.8)
    from importlib import resources as _ires
    from importlib import metadata as _im
except Exception:  # if ever needed on very old Pythons:
    import importlib_resources as _ires  # type: ignore
    import importlib_metadata as _im  # type: ignore

# -------- Versions / distributions --------
def version(dist_name: str) -> str:
    return _im.version(dist_name)

def distribution(dist_name: str):
    return _im.distribution(dist_name)

def distributions() -> Iterable:
    return _im.distributions()

# -------- Entry points (plugins) --------
def entry_points(group: Optional[str] = None, name: Optional[str] = None):
    eps = _im.entry_points()
    if hasattr(eps, "select"):  # Py>=3.10
        return eps.select(group=group, name=name) if (group or name) else eps
    if group:
        return [ep for ep in eps.get(group, []) if (name is None or ep.name == name)]
    return sum(eps.values(), [])

# -------- Package resources (data files) --------
def files(package):
    return _ires.files(package)

def open_text(package, resource):
    return _ires.open_text(package, resource)

def open_binary(package, resource):
    return _ires.open_binary(package, resource)

def path(package, resource):
    # yields a real filesystem path if needed (e.g., when bundled in wheels)
    target = files(package).joinpath(resource)
    return _ires.as_file(target) if hasattr(_ires, "as_file") else nullcontext(str(target))