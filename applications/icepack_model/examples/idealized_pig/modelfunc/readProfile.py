"""
Created on Fri Aug 14 08:44:51 2020
@author: ian
"""
import numpy as np
from modelfunc.myStrip import myStrip


def readProfile(profileFile, testVar):
    ''' Read a profile file with x y pairs
    # 2
    x0 y0
    x1 y1
    ...
    &
    With units of kilometers. Return results in meters along with distance
    along profile.

    Parameters
    ----------
    profileFile : str
        path of file with profile.
    testVar: firedrake function
        any function with meshdomain used to test if points in domain. Points
        not in domain are discarded.

    Returns
    -------
    xyProf : np array (npts, 2)
        xy array with points in m.
    d : np array (npts)
        Distance alogn profile in km.

    '''
    xyProf = myStrip(profileFile)  # read file
    # Compute distance along profile
    i = 0
    xyEdited = []
    xyLast = None
    for xy in xyProf:
        # use try in case profile has points outside of mesh
        xy[0] *= 1000.  # In meters
        xy[1] *= 1000.
        try:
            if i == 0:
                d = [0]
            else:
                _ = testVar((xy[0], xy[1]))  # Cause exception if outside mesh
                d.append(d[i-1] + np.sqrt((xy[0] - xyLast[0])**2 +
                                          (xy[1] - xyLast[1])**2))
            xyLast = xy
            xyEdited.append(xy)
            i += 1
        except Exception:
            continue
    d = np.array(d) * 0.001  # In km
    #
    xyProf = np.array(xyEdited)
    return xyProf, d
