import firedrake
from firedrake import inner
import icepack
import rasterio
import os
from modelfunc import myerror


def getModelVelocity(baseName, Q, V, minSigma=5, maxSigma=100):
    """Read in a tiff velocity data set and return
    firedrake interpolate functions.

    Parameters
    ----------
    baseName : str
        baseName should be of the form pattern.*.abc or pattern
        The wildcard (*) will be filled with the suffixes (vx, vy.)
        e.g.,pattern.vx.abc.tif, pattern.vy.abc.tif.
    Q : firedrake function space
        function space
    V : firedrake vector space
        vector space
    Returns
    -------
    uObs firedrake interp function on V
        velocity (m/yr)
    speed firedrake interp function on Q
        speed in (m)
    sigmaX firedrake interp function on Q
        vx error (m)
    sigmaY firedrake interp function on Q
        vy error (m)
    """
    # suffixes for products used
    suffixes = ['vx', 'vy', 'ex', 'ey']
    rasters = {}
    # prep baseName - baseName.*.xyz.tif or baseName.*
    if '*' not in baseName:
        baseName += '.*'
    if '.tif' not in baseName:
        baseName += '.tif'
    # read data
    for suffix in suffixes:
        myBand = baseName.replace('*', suffix)
        # print(myBand)
        if not os.path.exists(myBand):
            myerror(f'Velocity/error file - {myBand} - does not exist')
        rasters[suffix] = rasterio.open(myBand, 'r')
    # Firedrake interpolators
    uObs = icepack.interpolate((rasters['vx'], rasters['vy']), V)
    # force error to be at least 1 to avoid 0 or negatives.
    sigmaX = icepack.interpolate(rasters['ex'], Q)
    sigmaX = icepack.interpolate(firedrake.max_value(sigmaX, minSigma), Q)
    sigmaX = icepack.interpolate(firedrake.min_value(sigmaX, maxSigma), Q)
    sigmaY = icepack.interpolate(rasters['ey'], Q)
    sigmaY = icepack.interpolate(firedrake.max_value(sigmaY, minSigma), Q)
    sigmaY = icepack.interpolate(firedrake.min_value(sigmaY, maxSigma), Q)
    speed = icepack.interpolate(firedrake.sqrt(inner(uObs, uObs)), Q)
    # return results
    return uObs, speed, sigmaX, sigmaY
