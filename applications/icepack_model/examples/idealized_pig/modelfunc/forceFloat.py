import icepack
import firedrake
from modelfunc import flotationMask, flotationHeight
from icepack.constants import ice_density as rhoI


rhoW = rhoI * 1028./917.  # This ensures rhoW based on 1028


def forceFloat(floatingMask, bedElevation, sufaceElevation, Q, headRoom=50):
    """Adjust bed as needed to force flotation
    Parameters
    ----------
    floatingMask : firedrake function
        Mask where 1 indicates ice should be afloati
    bedElevation : firedrake function
        Original bed elevation
    sufaceElevation : firedrake function
        Surface elevation ice equivalent
    Q : firedrake function space
        Domain
    """
    nItMax = 20
    newBed = bedElevation.copy(deepcopy=True)
    dZ = 15
    i = 0
    print('start')
    while i < nItMax:
        zF = flotationHeight(newBed, Q, rhoI=rhoI, rhoW=rhoW)
        currentFlotation, _ = flotationMask(sufaceElevation, zF, Q, rhoI=rhoI,
                                            rhoW=rhoW)
        notFloating = firedrake.And(currentFlotation < firedrake.Constant(0.5),
                                    floatingMask > 0.5)
        area = firedrake.assemble(notFloating * firedrake.dx)
        newBed = icepack.interpolate(newBed - dZ * notFloating, Q)
        print(i, area/1e6)
        if area < 1:
            break
        i = i + 1
    newBed = icepack.interpolate(newBed - headRoom * floatingMask, Q)
    return newBed
