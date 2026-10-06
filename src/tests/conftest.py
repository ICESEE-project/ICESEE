# ==============================================================================
# @des: Session-wide guard that the test suite never deletes or recreates an
#       application's working data directory.
#
#       config/_utility_imports.py cleans data_path at import time (rank 0
#       removes it, every rank recreates it). A test that imports an
#       application, or launches one of its scripts, from the example
#       directory without its own --data_path therefore wipes that example's
#       real _modelrun_datasets* directory. Tests must isolate data_path
#       (tmp_path / tempfile); this guard makes any regression fail the run.
#
#       The inode is compared as well as the listing, so an empty directory
#       that was removed and recreated is also detected. Subprocess-launched
#       workers are covered because only the filesystem is inspected.
# ==============================================================================
from __future__ import annotations

import os
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_EXAMPLES_GLOB = "applications/*/examples/*/_modelrun_datasets*"


def _snapshot():
    snapshot = {}
    for path in sorted(_REPO_ROOT.glob(_EXAMPLES_GLOB)):
        if not path.is_dir():
            continue
        try:
            snapshot[str(path.relative_to(_REPO_ROOT))] = (
                os.stat(path).st_ino,
                tuple(sorted(os.listdir(path))),
            )
        except OSError:
            continue
    return snapshot


def application_data_changes(before, after):
    """Return a description of every application data directory that was
    deleted, recreated, or lost entries between two snapshots."""
    changes = []
    for rel, (inode, entries) in before.items():
        if rel not in after:
            changes.append(f"{rel}: deleted")
            continue
        new_inode, new_entries = after[rel]
        if new_inode != inode:
            changes.append(f"{rel}: deleted and recreated")
        missing = sorted(set(entries) - set(new_entries))
        if missing:
            changes.append(f"{rel}: lost {missing[:5]}")
    return changes


def pytest_sessionstart(session):
    session.config._icesee_app_data_snapshot = _snapshot()


def pytest_sessionfinish(session, exitstatus):
    before = getattr(session.config, "_icesee_app_data_snapshot", {})
    changes = application_data_changes(before, _snapshot())
    if changes:
        reporter = session.config.pluginmanager.get_plugin("terminalreporter")
        message = (
            "Application data directories were modified by the test run "
            "(a test is not isolating data_path):\n  " + "\n  ".join(changes)
        )
        if reporter is not None:
            reporter.write_sep("=", "application data isolation failure", red=True)
            reporter.write_line(message)
        session.exitstatus = pytest.ExitCode.TESTS_FAILED
