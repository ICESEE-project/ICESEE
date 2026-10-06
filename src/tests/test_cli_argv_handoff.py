# ==============================================================================
# @des: ICESEE consumes its own command-line arguments and leaves only the
#       rest in sys.argv, so a library that reads the command line afterwards
#       (PETSc, when Firedrake is imported) neither receives ICESEE's options
#       nor warns about them as unused, while its own options still work.
# ==============================================================================
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_LORENZ96_PARAMS = (
    _REPO_ROOT / "applications" / "lorenz_model" / "examples" / "lorenz96" / "params.yaml"
)

# The loader parses sys.argv at import time; give it a configuration and an
# isolated data_path instead of pytest's own arguments.
_saved_argv = sys.argv
sys.argv = [sys.argv[0], "-F", str(_LORENZ96_PARAMS), "--data_path",
            tempfile.mkdtemp(prefix="icesee_cli_argv_test_")]
try:
    from ICESEE.config import _utility_imports as cli_module
finally:
    sys.argv = _saved_argv


def test_unclaimed_argv_keeps_only_non_icesee_tokens():
    # The tokens argparse leaves over (named options such as --Nens or
    # --verbose are consumed before this point).
    argv = ["--num_years", "2", "-snes_monitor", "--dt=0.05", "-log_view"]
    assert cli_module._unclaimed_argv(argv) == ["-snes_monitor", "-log_view"]


def test_unclaimed_argv_matches_the_override_parser_on_the_same_tokens():
    argv = ["--obs_start_time", "0.5", "--flag", "--other=1", "-pc_type", "lu"]
    overrides = cli_module._parse_generic_cli_overrides(argv)
    assert overrides == {"obs_start_time": "0.5", "flag": "true", "other": "1"}
    assert cli_module._unclaimed_argv(argv) == ["-pc_type", "lu"]


def test_petsc_sees_only_its_own_options_after_the_icesee_config(tmp_path):
    pytest.importorskip("petsc4py")
    example = tmp_path / "example"
    example.mkdir()
    shutil.copy(_LORENZ96_PARAMS, example / "params.yaml")
    script = (
        "import sys\n"
        "import ICESEE.config._utility_imports\n"
        "import petsc4py\n"
        "petsc4py.init(sys.argv)  # as Firedrake does on import\n"
        "from petsc4py import PETSc\n"
        "print('OPTIONS', sorted(PETSc.Options().getAll()))\n"
        "print('ARGV', sys.argv[1:])\n"
    )
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{_REPO_ROOT.parent}:{_REPO_ROOT}"
    # PETSc options that take a value go through PETSc's own PETSC_OPTIONS
    # (a bare value token on the ICESEE command line is read as ICESEE's
    # optional positional distribution mode).
    env["PETSC_OPTIONS"] = "-icesee_test_petsc_value 7"
    for name in ("OMPI_COMM_WORLD_RANK", "PMIX_RANK", "PMI_RANK", "MV2_COMM_WORLD_RANK"):
        env.pop(name, None)
    script_path = tmp_path / "probe.py"
    script_path.write_text(script)
    result = subprocess.run(
        [sys.executable, str(script_path), f"--data_path={tmp_path / 'data'}",
         "--Nens=4", "--num_years=2", "--obs_start_time=0.5", "-icesee_test_petsc_flag"],
        cwd=str(example), env=env, capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    output = result.stdout + result.stderr
    options_line = next(line for line in result.stdout.splitlines() if line.startswith("OPTIONS"))
    assert "icesee_test_petsc_flag" in options_line
    assert "icesee_test_petsc_value" in options_line
    for icesee_key in ("Nens", "num_years", "obs_start_time", "data_path"):
        assert icesee_key not in options_line, icesee_key
    assert "ARGV ['-icesee_test_petsc_flag']" in result.stdout
    # PETSc may still list options nothing consumed in this probe, but never
    # an ICESEE argument.
    unused = [line for line in output.splitlines() if line.startswith("Option left:")]
    for icesee_key in ("Nens", "num_years", "obs_start_time", "data_path"):
        assert not any(icesee_key in line for line in unused), unused
