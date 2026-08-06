"""Coordinate provider shared by Icepack application backends."""


def get_icepack_node_coordinates(model_kwargs):
    """Return physical coordinates in the scalar-space DOF ordering.

    Icepack state blocks in ICESEE use ``Q``'s nodal ordering.  The mesh can
    either be exposed explicitly by the application or recovered from ``Q``.
    Imports stay local so non-Icepack applications do not require Firedrake.
    """
    import firedrake

    Q = model_kwargs.get("Q")
    if Q is None:
        raise ValueError(
            "Icepack coordinate provider requires the scalar function space "
            "'Q' in model_kwargs"
        )
    mesh = model_kwargs.get("mesh")
    if mesh is None:
        mesh = Q.mesh()

    coordinate_space = firedrake.VectorFunctionSpace(mesh, Q.ufl_element())
    coordinate_field = firedrake.assemble(
        firedrake.interpolate(firedrake.SpatialCoordinate(mesh), coordinate_space)
    )
    return coordinate_field.dat.data_ro.copy()


def register_icepack_coordinate_provider():
    """Register the shared provider under ICESEE's Icepack model key."""
    from ICESEE.src.utils.localization import register_coord_provider

    register_coord_provider("icepack", get_icepack_node_coordinates)

