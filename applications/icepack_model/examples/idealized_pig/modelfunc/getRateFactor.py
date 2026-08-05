#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 10:23:27 2020

@author: ian
"""
import rasterio
import icepack


def getRateFactor(rateFile, Q):
    """Read rate factor as B and convert to A
    Parameters
    ----------
    rateFile : str
        File with rate factor data
    Q : firedrake function space
        function space
    Returns
    -------
    A : firedrake function
        A Glenns flow law parameter
    """
    Bras = rasterio.open(rateFile)
    B = icepack.interpolate(Bras, Q)
    convFactor = 1e-6 * (86400*365.25)**(-1./3.)
    A = icepack.interpolate((B * convFactor)**-3, Q)
    return A
