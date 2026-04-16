#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 16:22:57 2020

@author: ian
"""
import firedrake
import numpy as np


def velocityError(uO, uI, area, message=''):
    """
    Compute and print velocity error
    """
    deltaV = uO - uI
    vError = firedrake.inner(deltaV, deltaV)
    vErrorAvg = np.sqrt(firedrake.assemble(vError * firedrake.dx) / area)
    print(f'{message} v error {vErrorAvg:10.2f} (m/yr)')
    return vErrorAvg.item()
