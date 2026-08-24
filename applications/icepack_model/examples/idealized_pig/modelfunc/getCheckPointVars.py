import firedrake
from modelfunc.myerror import myerror
import os
import warnings


def getCheckPointVars(checkFile, varNames, Q, t=None, mesh=None):
    """ Read a variable from a firedrake checkpoint file

    Parameters
    ----------
    checkFile : str
        checkfile name sans .h5
    varNames : str or list of str
        Names of variables to extract
    Q : firedrake function space
        firedrake function space can be a vector space, V, but not mixed
    Returns
    -------
    myVars: dict
        {'myVar':}
    """
    # Ensure a list since a single str is allowed
    if type(varNames) is not list:
        varNames = [varNames]
    if not os.path.exists(f'{checkFile}.h5'):
        myerror(f'getCheckPointVar: file {checkFile}.h5 does not exist')
    # open checkpoint
    try:
        # print('filter')
        myVars = getDumbCheckPointVars(checkFile, varNames, Q, t=t)
        # print('read dumb checkpointfile')
        return myVars
    except Exception:
        print('dumbcheckpoint failed, trying new checkpoint format')
    try:
        return getNewCheckPointVars(checkFile, varNames, Q, t=t, mesh=mesh)
    except Exception:
        myerror(f'Could not read new checkpoint file {checkFile}')
    return None


def getNewCheckPointVars(checkFile, varNames, Q, t=None, mesh=None):
    """ Read a variable from a firedrake checkpoint file

    Parameters
    ----------
    checkFile : str
        checkfile name sans .h5
    varNames : str or list of str
        Names of variables to extract
    Q : firedrake function space
        firedrake function space can be a vector space, V, but not mixed
    Returns
    -------
    myVars: dict
        {'myVar':}
    """
    myVars = {}
    if '.h5' not in checkFile:
        checkFile = f'{checkFile}.h5'
    # Open file and read variables
    with firedrake.CheckpointFile(checkFile, 'r') as chk:
        # Try mesh in file, otherwise use function space mesh
        if mesh is None:
            try:
                #mesh = chk.load_mesh()
                mesh = Q.mesh()
            except Exception:
                myerror(f"Could not load mesh from checkFile: {checkFile}")
        keywords = {}
        if t is not None:
            # note this only works for integer years
            #print('t', t)
            keywords['idx'] = int(t-1)
        #print(keywords)
        for varName in varNames:
            myVars[varName] = chk.load_function(mesh, varName, **keywords)
    return myVars


def getDumbCheckPointVars(checkFile, varNames, Q, t=None):
    """ Read a variable from a firedrake  dumb checkpoint file

    Parameters
    ----------
    checkFile : str
        checkfile name sans .h5
    varNames : str or list of str
        Names of variables to extract
    Q : firedrake function space
        firedrake function space can be a vector space, V, but not mixed
    Returns
    -------
    myVars: dict
        {'myVar':}
    """
    myVars = {}
    # Open file and read variables
    with warnings.catch_warnings(record=True):
        warnings.simplefilter('ignore', category=DeprecationWarning)
        with firedrake.DumbCheckpoint(checkFile, mode=firedrake.FILE_READ) as chk:
            if t is not None:
                # note this only works for integer years
                tt, ii = chk.get_timesteps()
                # print(tt[-1], ii[-1], len(tt))
                chk.set_timestep(t, idx=int(t-1))
            for varName in varNames:
                myVar = firedrake.Function(Q, name=varName)
                chk.load(myVar, name=varName)
                myVars[varName] = myVar
    return myVars