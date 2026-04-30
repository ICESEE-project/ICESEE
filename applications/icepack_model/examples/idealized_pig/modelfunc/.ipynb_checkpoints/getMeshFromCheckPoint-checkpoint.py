#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 08:51:49 2023

@author: ian
"""
import firedrake
import warnings


def getMeshFromCheckPoint(checkFile):
    '''
    If a new checkpoint file, read the mesh from there.
    '''
    if '.h5' not in checkFile:
        checkFile = f'{checkFile}.h5'
    # print(checkFile)

    with firedrake.CheckpointFile(checkFile, 'r') as chk:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=DeprecationWarning)
            try:
                mesh = chk.load_mesh()
                return mesh
            except Exception:
                return None
