import firedrake

Lx, Ly = 50e2, 12e2
nx, ny = 24, 16
mesh = firedrake.RectangleMesh(nx, ny, Lx, Ly)

Q = firedrake.FunctionSpace(mesh, "CG", 2)
V = firedrake.VectorFunctionSpace(mesh, "CG", 2)

x, y = firedrake.SpatialCoordinate(mesh)

# the bedrock slopes down from 200m ABS at the inflow boundary to -400m at the terminus
b_in, b_out = 200, -400
b = firedrake.interpolate(b_in - (b_in - b_out) * x / Lx, Q)
