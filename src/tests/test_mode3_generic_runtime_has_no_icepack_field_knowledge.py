# Keep model-specific field names out of generic native runtime modules.
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

# Field names belong to application adapters, not generic execution modules.
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
