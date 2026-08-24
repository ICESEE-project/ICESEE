from icepack.constants import ice_density as rhoI, water_density as rhoW
from modelfunc.firedrakeSmooth import firedrakeSmooth
import icepack


def flotationMask(s, zF, Q, rhoI=rhoI, rhoW=rhoW, flotationThreshold=-1,
                  alpha=100):
    """
    Using flotation height, create masks for floating and grounded ice.
    Normally flotation would be when surface == flotationHeight, which implies
    flotationHeight = 0. Using this value, however, can lead to interior lakes
    with artificially high melt rates. Setting to small value (e.g., 1) stops
    melt until the bottom the ice is a bit above the bed (e.g., ~10m), avoids
    this problem while making little difference or no difference for steep
    grounding lines.

    Parameters
    ----------
    zF firedrake interp function
        Flotation height (m)
    Q : firedrake function space
        function space
    rhoI : [float], optional
        density of ice, by default rhoI
    rhoW : [float], optional
        density of seawater, by default rhoW
    flotationThresh, [float] optional
        threshold to determine floating applied to height above flotation
    alpha, [float] optional
        smoothing parameter to help avoid lakes
    Returns
    -------
    floating firedrake interp function
         ice shelf mask 1 floating, 0 grounded
    grounded firedrake interp function
        Grounded mask 1 grounded, 0 floating
    """
    # smooth to avoid isolated points dipping below flotation.
    zAbove = firedrakeSmooth(icepack.interpolate(s - zF, Q), alpha=100)
    floating = icepack.interpolate(zAbove < flotationThreshold, Q)
    grounded = icepack.interpolate(zAbove > flotationThreshold, Q)
    return floating, grounded
