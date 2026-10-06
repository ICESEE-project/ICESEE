# ==============================================================================
# @des: Architectural regression test (2026-09-28 reconciliation, item 2):
# the generic Mode-3 execution infrastructure -- distributed runtime,
# streaming runtime, member stores, scheduling/registry, projections --
# must never encode knowledge of any specific application's field names,
# state-vector bookkeeping flags, or physics. That knowledge belongs
# exclusively to the application's own adapter (e.g. applications/
# icepack_model/examples/idealized_pig/_icepack_native.py for Icepack).
#
# This is a plain source scan, not a behavioral test: it fails loudly (and
# names the offending file/line) the moment any generic runtime module
# starts referencing an Icepack-specific field name or bookkeeping flag,
# which is exactly the class of abstraction leak this reconciliation fixed
# in _icepack_native.py (a `joint_estimation`-keyed predicate that should
# have been driven by the already-generic `vec_inputs` instead).
# ==============================================================================
from __future__ import annotations

from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
_PARALLELIZATION_DIR = _REPO_ROOT / "ICESEE" / "src" / "parallelization"
if not _PARALLELIZATION_DIR.is_dir():
    _PARALLELIZATION_DIR = _REPO_ROOT / "src" / "parallelization"

# Every module that forms the generic Mode-3 execution path: topology,
# native/streaming runtime, member stores, field registry/adapter
# contracts, scheduling, checkpointing, and the pure-arithmetic
# projection formulas. Explicitly excludes nothing application-specific
# lives here by construction (see this file's own docstring).
_GENERIC_RUNTIME_FILES = sorted(
    p for p in _PARALLELIZATION_DIR.glob("distributed_*.py")
) + [_PARALLELIZATION_DIR / "mode3_projections.py"]

# Field names / bookkeeping flags that are legitimately known ONLY to a
# specific application's own adapter file -- never to generic execution
# infrastructure. This list is Icepack-focused (the application this
# reconciliation touched) but the test's shape generalizes to any
# application: if a future app's adapter needs a new field name, that
# name must never appear in these generic files either.
_FORBIDDEN_SUBSTRINGS = (
    "basal_melt_field",
    "joint_estimation",
    "wrong_basal_melt_field",
    "EXPERIMENT_TRUE",
    "EXPERIMENT_WRONG",
    "BasalMeltRate",
)


def test_generic_runtime_files_exist():
    assert len(_GENERIC_RUNTIME_FILES) >= 10, (
        f"expected to find the generic Mode-3 runtime modules under "
        f"{_PARALLELIZATION_DIR}, found only {_GENERIC_RUNTIME_FILES}"
    )


def test_generic_mode3_runtime_never_mentions_icepack_field_names():
    offenses = []
    for path in _GENERIC_RUNTIME_FILES:
        text = path.read_text()
        for lineno, line in enumerate(text.splitlines(), start=1):
            for token in _FORBIDDEN_SUBSTRINGS:
                if token in line:
                    offenses.append(f"{path.name}:{lineno}: contains {token!r}: {line.strip()}")
    assert not offenses, (
        "Generic Mode-3 runtime files must never reference Icepack-specific "
        "field names or bookkeeping flags -- that knowledge belongs only in "
        "the application adapter (e.g. _icepack_native.py). Offenses:\n"
        + "\n".join(offenses)
    )
