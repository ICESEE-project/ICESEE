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
    #
    # firedrake.interpolate(expr, Q) (the free-function form) now returns a
    # lazy symbolic Interpolate object in the installed Firedrake/UFL
    # version, not an evaluated Function -- Function(Q).interpolate(expr)
    # is the eager-evaluation equivalent (same math, confirmed against
    # icepack.interpolate's own behavior, which already used this pattern
    # and was unaffected).
    zF = firedrake.Function(Q).interpolate(firedrake.max_value(-zb * (rhoW/rhoI-1), 0))
    return zF
