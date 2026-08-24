__all__ = ['argusMesh', 'argusToFiredrakeMesh', 'divMelt', 'firedrakeSmooth',
           'flotationHeight', 'flotationMask', 'forceFloat',
           'getCheckPointVars', 'getNewCheckPointVars',
           'getInversionData', 'getMeshFromCheckPoint',
           'getModelGeometry',
           'getRateFactor', 'getModelVarFromTiff', 'getModelVelocity',
           'gMesh', 'inputMeltParams', 'myerror', 'piecewiseWithDepth',
           'myStrip',
           'readModelParams', 'reduceNearGLBeta', 'readProfile',
           'schoofFriction', 'setupMesh',  'viscosityTheta', 'velocityError',
           'weertmanFriction']
from modelfunc.argusMesh import argusMesh
# from modelfunc.shpDomain import shpDomain
from modelfunc.argusToFiredrakeMesh import argusToFiredrakeMesh
from modelfunc.meltFunctions import divMelt
from modelfunc.firedrakeSmooth import firedrakeSmooth
from modelfunc.flotationHeight import flotationHeight
from modelfunc.flotationMask import flotationMask
from modelfunc.forceFloat import forceFloat
from modelfunc.getCheckPointVars import getCheckPointVars
from modelfunc.getCheckPointVars import getNewCheckPointVars
from modelfunc.getInversionData import getInversionData
from modelfunc.getMeshFromCheckPoint import getMeshFromCheckPoint
from modelfunc.getModelGeometry import getModelGeometry
from modelfunc.getModelGeometry import getModelVarFromTiff
from modelfunc.getModelVelocity import getModelVelocity
from modelfunc.getRateFactor import getRateFactor
from modelfunc.gMesh import gMesh
from modelfunc.meltFunctions import inputMeltParams
from modelfunc.meltFunctions import piecewiseWithDepth
from modelfunc.myerror import myerror
from modelfunc.readModelParams import readModelParams
from modelfunc.reduceNearGLBeta import reduceNearGLBeta
from modelfunc.readProfile import readProfile
from modelfunc.slidingAndViscosity import schoofFriction
from modelfunc.setupMesh import setupMesh
from modelfunc.myStrip import myStrip
from modelfunc.velocityError import velocityError
from modelfunc.slidingAndViscosity import viscosityTheta
from modelfunc.slidingAndViscosity import weertmanFriction
