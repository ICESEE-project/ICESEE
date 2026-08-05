import modelfunc as m
import firedrake
import uuid
import os
import yaml


def argusToFiredrakeMesh(meshFile, savegmsh=False):
    """Convert an argus mesh to firedrake
    Parameters
    ----------
    meshFile : string
        File name of an Argus mesh export
    savegmsh : bool
        Save gmsh file
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
    mesh = firedrake.Mesh(gmshFile)
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
