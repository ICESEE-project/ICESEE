# ==============================================================================
# @des: Regression tests for application-data isolation in the test suite.
#
#       config/_utility_imports.py removes and recreates data_path at import
#       time. Launched from an example directory without --data_path, it
#       therefore wipes that example's real _modelrun_datasets directory; this
#       is what emptied idealized_pig/_modelrun_datasets during test runs.
#       These tests run the real config import from an example-like directory
#       that already holds a dataset and check that an isolated --data_path
#       leaves it untouched, with a control showing the default path does not,
#       plus unit checks of the session guard in conftest.py.
# ==============================================================================
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

from ICESEE.src.tests.conftest import application_data_changes

_REPO_ROOT = Path(__file__).resolve().parents[2]
_LORENZ96_PARAMS = (
    _REPO_ROOT / "applications" / "lorenz_model" / "examples" / "lorenz96" / "params.yaml"
)


def _example_with_existing_dataset(tmp_path):
    example = tmp_path / "example"
    dataset = example / "_modelrun_datasets"
    dataset.mkdir(parents=True)
    shutil.copy(_LORENZ96_PARAMS, example / "params.yaml")
    (dataset / "existing_result.h5").write_bytes(b"precious")
    return example, dataset


def _import_config(cwd, *argv):
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{_REPO_ROOT.parent}:{_REPO_ROOT}"
    for name in ("OMPI_COMM_WORLD_RANK", "PMIX_RANK", "PMI_RANK", "MV2_COMM_WORLD_RANK"):
        env.pop(name, None)
    return subprocess.run(
        [sys.executable, "-c", "import ICESEE.config._utility_imports", *argv],
        cwd=str(cwd), env=env, capture_output=True, text=True, timeout=120,
    )


def test_isolated_data_path_preserves_existing_application_dataset(tmp_path):
    example, dataset = _example_with_existing_dataset(tmp_path)
    inode = os.stat(dataset).st_ino
    isolated = tmp_path / "isolated_data_path"

    result = _import_config(example, f"--data_path={isolated}")

    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    assert (dataset / "existing_result.h5").read_bytes() == b"precious"
    assert os.stat(dataset).st_ino == inode
    assert isolated.is_dir()


def test_default_data_path_is_cleaned_which_is_why_tests_must_isolate(tmp_path):
    # Control: production behavior cleans the configured data_path, so a test
    # that omits --data_path from an example directory deletes real results.
    example, dataset = _example_with_existing_dataset(tmp_path)

    result = _import_config(example)

    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    assert dataset.is_dir()
    assert not (dataset / "existing_result.h5").exists()


def test_session_guard_detects_deleted_recreated_and_emptied_directories():
    before = {
        "a/_modelrun_datasets": (1, ("x.h5",)),
        "b/_modelrun_datasets": (2, ()),
        "c/_modelrun_datasets": (3, ("y.h5",)),
        "d/_modelrun_datasets": (4, ("z.h5",)),
    }
    after = {
        "b/_modelrun_datasets": (20, ()),
        "c/_modelrun_datasets": (3, ()),
        "d/_modelrun_datasets": (4, ("z.h5", "new.h5")),
    }
    changes = application_data_changes(before, after)
    assert changes == [
        "a/_modelrun_datasets: deleted",
        "b/_modelrun_datasets: deleted and recreated",
        "c/_modelrun_datasets: lost ['y.h5']",
    ]
