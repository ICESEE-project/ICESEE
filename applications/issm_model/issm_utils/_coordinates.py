"""Coordinate provider shared by ISSM application backends."""

import os

import h5py
import numpy as np


def get_issm_node_coordinates(model_kwargs):
    """Load physical ISSM nodal coordinates in state-vector ordering."""
    data_path = model_kwargs.get("data_path", "_modelrun_datasets")
    icesee_path = model_kwargs.get("icesee_path", os.getcwd())
    mesh_path = os.path.join(icesee_path, data_path, "mesh_idxy_0.h5")

    with h5py.File(mesh_path, "r") as mesh_file:
        x_coord = np.asarray(mesh_file["fric_x"][:], dtype=float).ravel()
        y_coord = np.asarray(mesh_file["fric_y"][:], dtype=float).ravel()

    if x_coord.size != y_coord.size:
        raise ValueError(
            f"ISSM coordinate sizes disagree in {mesh_path}: "
            f"x={x_coord.size}, y={y_coord.size}"
        )
    return np.column_stack((x_coord, y_coord))


def register_issm_coordinate_provider():
    """Register the shared provider under ICESEE's ISSM model key."""
    from ICESEE.src.utils.localization import register_coord_provider

    register_coord_provider("issm", get_issm_node_coordinates)

