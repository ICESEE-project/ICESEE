import modelfunc as m
import firedrake
import uuid
import os
import yaml


def argusToFiredrakeMesh(meshFile, savegmsh=False, comm=None):
    """Convert an argus mesh to firedrake
    Parameters
    ----------
    meshFile : string
        File name of an Argus mesh export
    savegmsh : bool
        Save gmsh file
    comm : mpi4py.MPI.Comm, optional
        Communicator to build the mesh on. Root-caused Mode-3 hang fix:
        firedrake.Mesh() defaults to COMM_WORLD when comm= is omitted
        (confirmed directly against the installed firedrake/mesh.py
        source, kwargs.get("comm", COMM_WORLD)). Every caller up the
        chain (initializeMesh) already receives a correctly-scoped
        communicator (topology.spatial_comm for Mode-3's per-ensemble-
        -group setup, a per-rank subcomm for the modes-0-2-style
        top-level call) -- it was simply never forwarded this deep, so
        this call silently became a real COMM_WORLD collective for every
        caller regardless of their own, different, intended scope.
        Confirmed via direct instrumentation (ICESEE_MODE3_TRACE=1) that
        this is the exact operation where world_size>1 Mode-3 runs hang:
        one rank/ensemble-group reaches this call while its COMM_WORLD
        "peer" is elsewhere in independent, unsynchronized work, so the
        implicit collective never completes.
    Return
    ----------
    mesh : firedrake mesh
    opts : dict
    """
    myMesh = m.argusMesh(meshFile)  # Input the mesh
    # Quick Check that file is a value argus mesh
    if '.exp' not in meshFile:
        m.myerror(f'Invalid mesh file [{meshFile}]: missing .exp')
    # create unique ide to avoid multiple jobs overwriting
    myId = f'.{uuid.uuid4().hex[:8]}.msh'
    gmshFile = meshFile.replace(".exp", myId)
    myMesh.toGmsh(gmshFile)  # Save in gmsh format
    _mesh_comm = comm if comm is not None else firedrake.COMM_WORLD
    mesh = firedrake.Mesh(gmshFile, comm=_mesh_comm)
    # delete gmsh file
    if not savegmsh:
        os.remove(gmshFile)
    else:
        os.rename(gmshFile, gmshFile.replace(myId, '.msh'))
    opts = readOpts(meshFile)
    # return mesh
    return mesh, opts


def readOpts(meshFile):
    """Read and opts file if there is one
    Parameters
    ----------
    meshFile : str
        opts file name
    Return
    ----------
    opts: dict
        opts data
    """
    optsFile = meshFile.replace('.exp', '.yaml')
    if not os.path.exists(optsFile):
        return {}
    with open(optsFile, 'r') as fp:
        try:
            opts = yaml.load(fp, Loader=yaml.FullLoader)
        except Exception:
            m.myerror(f'Error reading opts file: {optsFile}')
    if 'opts' in opts:
        opts = opts['opts']
    return opts
