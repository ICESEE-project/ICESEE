import icepack
import firedrake


def reduceNearGLBeta(s, sOrig, zF, grounded, Q, thresh, limit=False):
    """Compute beta reduction in area near grounding line where the height
    above floation is less than thresh.

    Parameters
    ----------
    s : firedrake function
        Current surface.
    sOrig : firedrake function
        Original surface.
    zF : firedrake function
        Flotation height.
    grounded : firedrake function
        grounde mask
    Q : firedrake function space
        scaler function space for model
    thresh : float
        Threshold to determine where to reduce beta (zAbove < thresh)
    """
    # compute original height above flotation for grounded region
    sAboveOrig = (sOrig - zF) * grounded
    # avoid negative/zero values
    sAboveOrig = firedrake.max_value(sAboveOrig, 0.0)
    # Current height above flotation with negative values zeroed out.
    sAbove = firedrake.max_value((s - zF) * grounded, 0)
    # mask so only areas less than thresh but grounded
    # This is the area below hT
    sMask = (sAbove < thresh) * grounded
    # print(f'{sAbove.dat.data_ro.min()}, {sAbove.dat.data_ro.max()} {thresh}')
    # Inverse mask (the area above hT)
    sMaskInv = sMask < 1
    # scale = fraction of original height above flotation
    # Use 5 to avoid potentially large ratio at small values.
    scaleBeta = sAbove / \
        firedrake.max_value(firedrake.min_value(thresh, sAboveOrig), 3)
    if limit:
        scaleBeta = firedrake.min_value(3, scaleBeta)
    scaleBeta = scaleBeta * sMask + sMaskInv
    # scaleBeta = icepack.interpolate(firedrake.min_value(scaleBeta,1.),Q)
    # sqrt so tau = scale * beta^2
    # scaleBeta = icepack.interpolate(firedrake.sqrt(scaleBeta) * grounded, Q)
    # Removed grounded above because grounded is always applied in friction
    scaleBeta = icepack.interpolate(firedrake.sqrt(scaleBeta), Q)
    # print(f'{scaleBeta.dat.data_ro.min()}, {scaleBeta.dat.data_ro.max()}')
    return scaleBeta
