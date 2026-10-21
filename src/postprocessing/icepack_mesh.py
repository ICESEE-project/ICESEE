"""Export the Icepack mesh's degree-of-freedom coordinates for post-processing.

Icepack states are stored in the scalar function space's degree-of-freedom
ordering, which only Firedrake knows. Run this once, in the environment the
model ran in, after the run::

    python -m ICESEE.src.postprocessing.icepack_mesh params.yaml _modelrun_datasets

It rebuilds the model's mesh exactly as ``_icepack_model.py`` does
(``RectangleMesh(nx, ny, Lx, Ly, quadrilateral=True)``, ``CG(degree)``) and
writes ``dof_coordinates.h5`` (dataset ``coordinates``, shape ``(n, 2)``)
with ICESEE's own coordinate recipe (``icepack_utils/_coordinates.py``), so
maps can be drawn later without Firedrake.
"""
from __future__ import annotations

import sys
from pathlib import Path

from .outputs import DOF_COORDINATES


def export_dof_coordinates(params_file: Path | str, data_dir: Path | str) -> Path:
    import firedrake
    import h5py
    import yaml

    params = yaml.safe_load(Path(params_file).read_text(encoding="utf-8")) or {}
    phys = params.get("physical-parameters") or {}
    Lx, Ly = int(float(phys["Lx"])), int(float(phys["Ly"]))
    nx, ny = int(float(phys["nx"])), int(float(phys["ny"]))
    degree = int(float(phys.get("degree", 2)))
    from firedrake import COMM_SELF

    mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly, quadrilateral=True, comm=COMM_SELF)
    Q = firedrake.FunctionSpace(mesh, "CG", degree)
    space = firedrake.VectorFunctionSpace(mesh, Q.ufl_element())
    coords = firedrake.assemble(firedrake.interpolate(firedrake.SpatialCoordinate(mesh), space))
    out = Path(data_dir) / DOF_COORDINATES
    out.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out, "w") as f:
        ds = f.create_dataset("coordinates", data=coords.dat.data_ro.copy())
        ds.attrs.update({"Lx": Lx, "Ly": Ly, "nx": nx, "ny": ny, "degree": degree,
                         "mesh": "RectangleMesh(quadrilateral=True)", "space": f"CG{degree}"})
    return out


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: python -m ICESEE.src.postprocessing.icepack_mesh PARAMS_YAML DATA_DIR")
    print(export_dof_coordinates(sys.argv[1], sys.argv[2]))
