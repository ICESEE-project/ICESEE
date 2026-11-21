"""Run a small configured Icepack application using the existing model hooks.

No solver/model implementation is copied. The established entrypoint keeps its
own working-directory resolution and application-scoped hook registration.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
import runpy
import sys


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True,
                        help="new dedicated absolute data directory")
    parser.add_argument("--mode", type=int, choices=(0, 1, 2, 3), default=2,
                        help="execution mode; native mode 3 uses one spatial rank per member")
    parser.add_argument("--check", action="store_true",
                        help="inspect merged YAML/entrypoint without model imports or output writes")
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    repository = here.parents[5]
    example = here.parents[1]
    config = here / "params.yaml"
    driver = example / "run_da_icepack.py"
    loader = repository / "config/config_loader.py"
    if not driver.is_file() or not loader.is_file():
        parser.error("keep this wrapper in the synthetic_ice_stream/tutorials/coarse_assimilation tree")
    if not args.output.is_absolute():
        parser.error("--output must be an absolute path to a new run directory")
    output = args.output.resolve()
    if output in {Path('/'), Path.home().resolve(), repository, example}:
        parser.error("output must be a dedicated run subdirectory")
    rank_hint = int(next((os.environ[key] for key in
        ("OMPI_COMM_WORLD_RANK", "PMIX_RANK", "PMI_RANK", "MV2_COMM_WORLD_RANK")
        if key in os.environ), "0"))
    # Only the root checks: a non-root can enter later after the application
    # has already created this directory. Existing-data refusal remains on
    # the root launch path and a root failure terminates the MPI command.
    if rank_hint == 0 and (output.exists() or output.is_symlink()):
        parser.error(f"output already exists; preserve it and choose another path: {output}")
    spec = importlib.util.spec_from_file_location("icesee_teaching_config_loader", loader)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    merged = module.load_yaml_to_dict(str(config))
    physical, modeling, ensemble = (merged[key] for key in
        ("physical-parameters", "modeling-parameters", "enkf-parameters"))
    if rank_hint == 0:
        print(json.dumps({
            "application": "synthetic_ice_stream/coarse_assimilation",
            "driver": str(driver), "configuration": str(config),
            "output": str(output), "execution_mode": args.mode,
            "mesh_elements": [physical["nx"], physical["ny"]],
            "degree": physical["degree"], "members": ensemble["Nens"],
            "forecast_steps": int(modeling["num_years"] * modeling["timesteps_per_year"]),
            "expected_serial_state_shape": [140, 4, 11],
            "expected_observation_state_columns": [5, 10],
            "requires_external_input_dataset": False,
            "check_only": args.check,
        }, indent=2), flush=True)
    if args.check:
        return
    for dependency in ("firedrake", "icepack"):
        if importlib.util.find_spec(dependency) is None:
            parser.error(f"{dependency} is unavailable; activate the approved Icepack/Firedrake environment")
    sys.path[:0] = [str(repository.parent), str(repository)]
    sys.argv = [str(driver), "-F", str(config),
                f"--execution_mode={args.mode}", f"--data_path={output}",
                "--ranks_per_model=1", "--model_nprocs=1"]
    runpy.run_path(str(driver), run_name="__main__")


if __name__ == "__main__":
    main()
