from firedrake import inner, grad, dx
import firedrake


def firedrakeSmooth(q0, alpha=2e3):
    """[summary]

    Parameters
    ----------
    q0 : firedrake function
        firedrake function to be smooth
    alpha : float, optional
        parameter that controls the amount of smoothing, which is
        approximately the smoothing lengthscale in m, by default 2e3
    Returns
    -------
    q firedrake interp function
        smoothed result
    """
    q = q0.copy(deepcopy=True)
    J = 0.5 * ((q - q0)**2 + alpha**2 * inner(grad(q), grad(q))) * dx
    F = firedrake.derivative(J, q)
    firedrake.solve(F == 0, q)
    return q
