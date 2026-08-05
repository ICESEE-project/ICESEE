import firedrake
from icepack.constants import ice_density as rhoI, water_density as rhoW

rhoW = rhoI * 1028./917.  # This ensures rhoW based on 1028


def flotationHeight(zb, Q, rhoI=rhoI, rhoW=rhoW):
    """Given bed elevation, determine height of flotation for function space Q.

    Parameters
    ----------
    zb  : firedrake interp function
        bed elevation (m)
    Q : firedrake function space
        function space
    rhoI : [type], optional
        [description], by default rhoI
    rhoW : [type], optional
        [description], by default rhoW
    Returns
    -------
    zF firedrake interp function
        Flotation height (m)
    """
    # computation for height above flotation
    zF = firedrake.interpolate(firedrake.max_value(-zb * (rhoW/rhoI-1), 0), Q)
    return zF
