from modelfunc import argusToFiredrakeMesh
import firedrake
from modelfunc.myerror import myerror


def setupMesh(meshFile, degree=2, meshOversample=None, savegmsh=False,
              newMesh=None, comm=None):
    """
    Read argus mesh file and return mesh alongw ith function spaces

    Parameters
    ----------
    meshFile : str
        argus meshfile name
    degree : int, optional
        degree of function spaces, by default 2
    comm : mpi4py.MPI.Comm, optional
        Communicator to build the (discarded-unless-newMesh-is-None)
        Argus mesh on. See argusToFiredrakeMesh's own docstring for the
        Mode-3 hang this forwarding fixes -- previously silently dropped
        here, defaulting the call three levels down to COMM_WORLD
        regardless of the caller's actual intended scope.

    Returns
    -------
    mesh
        firedrake mesh
    Q, V
        firedrake scalar and vectory functions
    """
    # Input the mesh
    maxOversample = 4  # Arbitrary could be increased
    mesh1, opts = argusToFiredrakeMesh(meshFile, savegmsh=savegmsh, comm=comm)

    if meshOversample is not None and newMesh is None:
        numLevels = meshOversample - 1
        if numLevels < 0 or numLevels > (maxOversample-1):
            myerror(f'meshOverample={meshOversample} but  0 < '
                    'meshOversample < 4')
        mesh1 = firedrake.MeshHierarchy(mesh1, numLevels)[numLevels]
    # Create scalar and vector function spaces
    if newMesh is not None:
        mesh = newMesh
        #myerror("stop")
    else:
        mesh = mesh1
    # print(type(newMesh))
    # print(type(mesh1))
    Q = firedrake.FunctionSpace(mesh, family='CG', degree=degree)
    V = firedrake.VectorFunctionSpace(mesh, family='CG', degree=degree)
    return mesh, Q, V, opts
