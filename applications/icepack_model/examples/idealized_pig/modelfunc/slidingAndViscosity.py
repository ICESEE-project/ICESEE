import firedrake
import icepack
from icepack.constants import weertman_sliding_law as m


def viscosityTheta(u, h, grounded, floating, A, theta):
    """Combine existing A on grounded ice with theta^2 from floating ice.
    Parameters
    ----------
    u : firedrake function
        modelled velocity
    h : firedrake function
        model thickness
    grounded : firedrake function
        Mask with 1s for grounded 0 for floating.
    floating : firedrake function
        Mask with 1s for floating 0 for grounded.
    A : firedrake function
        A for grounded areas
    theta : firedrake function
        theta for A on floating ice as theta^2
    Returns
    -------
    viscosityFunction function
        Viscosity funciton for model
    """
    # A = grounded * A + floating * theta**2
    A0 = grounded * A + floating * firedrake.exp(theta)
    # A 0 = firedrake.max_value(firedrake.min_value(A0, 1.), 100.)
    return icepack.models.viscosity.viscosity_depth_averaged(velocity=u,
                                                             thickness=h,
                                                             fluidity=A0)


def weertmanFriction(velocity, grounded, beta, uThresh):
    """Friction coefficient for weertman sliding law
    Parameters
    ----------
    velocity : firedrake function
        modelled velocity
    grounded : firedrake function
        Mask with 1s for grounded 0 for floating.
    beta : firedrake function
        Beta for friction model with C=beta^2

    Returns
    -------
    frictionFunction function
        Bed friction model
    """
    C = grounded * beta**2
    return icepack.models.friction.bed_friction(velocity=velocity, friction=C)


# Note m is a global variable
def schoofFriction(velocity, grounded, beta, uThresh):
    """Friction coefficient for regularized coloumb (Schoof) sliding law
    Parameters
    ----------
    velocity : firedrake function
        modelled velocity
    grounded : firedrake function
        Mask with 1s for grounded 0 for floating.
    beta : firedrake function
        Beta for friction model with C=beta^2
    uThresh : float
        Velocity transition for sliding law
    Returns
    -------
    frictionFunction function
        Bed friction model
    """
    C = grounded * beta**2
    mExp = (1./m + 1.)
    U = firedrake.sqrt(firedrake.inner(velocity, velocity))
    return C * ((uThresh**mExp + U**mExp)**(m/(m+1.)) - uThresh)
