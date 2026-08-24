import rasterio
import icepack
import firedrake
import yaml
from modelfunc.myerror import myerror
from modelfunc import firedrakeSmooth, flotationMask, flotationHeight
from icepack.constants import ice_density as rhoI
import os

rhoW = rhoI * 1028./917.  # This ensures rhoW based on 1028


def getModelGeometry(geometryFile, Q, smooth=False, alpha=2e3, zFirn=0.,
                     rhoI=rhoI, rhoW=rhoW):
    """Load geometry data for model and create firedrake interpolators
    Parameters
    ----------
    geometryFile : str
        Path to a yaml file with bed, surface, thickness, and floatMask
    Q : firedrake function space
        function space
    smooth: bool, optional
        apply firedrakeSmooth to the result
    alpha : float, optional
        parameter that controls the amount of smoothing, which is approximately
        the smoothing lengthscale in m, by default 2e3
    zFirn : float, optional
        Correct elevation for firn thickness (m), by default 14 m
    rhoI : [type], optional
        [description], by default rhoI
    rhoW : [type], optional
        [description], by default rhoW
    Returns
    -------
    zb firedrake interp function
        bed elevation (m)
    s firedrake interp function
        surface elevation (m)
    h firedrake interp function
        ice thickness (m)
    floatMask firedrake interp function
        mask with 1 for floating 0 for grounded
    """
    # load geometry files
    try:
        with open(geometryFile) as fp:
            geom = yaml.load(fp, Loader=yaml.FullLoader)
    except Exception:
        myerror(f'Could not open geomtery file: {geometryFile}')
    # Load and convert to firedrake
    fd = {'bed': None, 'surface': None, 'thickness': None, 'floatMask': None}
    # Read and process data
    for myVar in geom:
        print(myVar, geom[myVar])
        fd[myVar] = getModelVarFromTiff(geom[myVar], Q)
        if smooth and alpha > 1 and myVar != 'floatMask':
            fd[myVar] = firedrakeSmooth(fd[myVar], alpha=alpha)
        if myVar == 'surface':
            fd[myVar] = icepack.interpolate(fd[myVar] - zFirn, Q)
    # If data are smoothed, regenerate a new mask from smoothed results.
    if smooth and alpha > 1:
        zF = flotationHeight(fd['bed'], Q, rhoI=rhoI, rhoW=rhoW)
        fd['floatMask'], g = flotationMask(fd['surface'], zF, Q, rhoI=rhoI,
                                           rhoW=rhoW)
    else:
        g = icepack.interpolate(fd['floatMask'] < 1, Q)
    # Don't allow negative values
    for myVar in ['surface', 'thickness']:
        fd[myVar] = icepack.interpolate(firedrake.max_value(10, fd[myVar]), Q)
    for myVar in geom:
        print(f'{myVar} min/max {fd[myVar].dat.data_ro.min():10.2f} '
              f'{fd[myVar].dat.data_ro.max():10.2f}')
    return fd['bed'], fd['surface'], fd['thickness'], fd['floatMask'], g


def getModelVarFromTiff(myTiff, Q):
    """Read a model variable from a tiff file using rasterio
    Parameters
    ----------
    myTiff : str
        tiff file with a scalar variable
    Q : firedrake function space
        function space
    Returns
    -------
    firedrake function
        Data from tiff
    """
    if not os.path.exists(myTiff):
        myerror(f'Geometry file {myTiff} does not exist')
    x = rasterio.open(myTiff)
    return icepack.interpolate(x, Q)
