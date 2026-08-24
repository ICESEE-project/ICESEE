#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 09:10:16 2020

@author: ian
"""
import modelfunc as mf


def getInversionData(inversionCheckpointFile, Q, V, mesh=None):
    '''
    Read an inversion results file and return
    betaInv, thetaInv, Ainv, sInv, hInv, zbInv, floatingInv,
    groundedInv, uInv, uObsInv

    Parameters
    ----------
    inversionCheckpointFile : str
        Check point filename WITHOUT .h5 suffix
    Q : firedrake function space
        modal scalar function space.
    V : firedrake function space
        model vector function space.
    mesh: Mesh for new checkpoint function
    chk: Echo data to chk file

    Returns
    -------
    betaInv: firedrake
        Beta solution from inversion
    thetaInv: firedrake
        Theta solution from inversion
    Ainv: firedrake
        Initial A used in inversion
    sInv: firedrake
        Elevation used in inversion
    hInv: firedrake
        Thickness used in inversion
    zbInv: firedrake
        Bed elevation used in version
    floatingInv:
        floating masking used in inversion
    groundedInv: firedrake
        grounded mask used in inversion
    uInv: firedrake
        velocity solution from inversion
    uObsInv: firedrake
        observed velocity used to constrain inversion
    '''
    #
    inversionData = ['betaInv', 'thetaInv', 'AInv', 'sInv', 'hInv', 'zbInv',
                     'floatingInv', 'groundedInv']
    myVars = mf.getCheckPointVars(inversionCheckpointFile, inversionData, Q, mesh=mesh)
    myList = [myVars[myVar] for myVar in inversionData]
    vels = mf.getCheckPointVars(inversionCheckpointFile, ['uInv', 'uObsInv'],
                                V, mesh=mesh)
    myList += [vels['uInv'], vels['uObsInv']]
    #  
    # betaInv, thetaInv, Ainv,  sInv, hInv, zbInv, floatingInv,
    # groundedInv, uInv, uObsInv
    return myList[:]
