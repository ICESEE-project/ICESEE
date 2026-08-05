#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 08:41:43 2020

@author: ian
"""

import shapefile
import numpy as np
import pygmsh
import matplotlib.pyplot as plt


class shpDomain:

    def __init__(self):
        self.geometry = pygmsh.built_in.Geometry()

    def addPoints(self, xyPoly):
        ''' add points to geometry and return list of points '''
        myPoints = []
        for xy in xyPoly:
            myPoints.append(self.geometry.add_point([xy[0], xy[1], 0],
                                                    lcar=10000.))
        return myPoints

    def checkPtExist(self, pt, priorPts):
        ''' Subsititute existing point if available '''
        for pPt in priorPts:
            d = (pt.x[0]-pPt.x[0])**2 + (pt.x[1]-pPt.x[1])**2
            if d < 0.001:
                return pPt
        return pt

    def addLines(self, myPoints, priorPoints):
        ''' Create a list of lines form points '''
        myLines = []
        for pt1, pt2 in zip(myPoints[0:-1], myPoints[1:]):
            pt1 = self.checkPtExist(pt1, priorPoints)
            pt2 = self.checkPtExist(pt2, priorPoints)
            myLines.append(self.geometry.add_line(pt1, pt2))
        return myLines

    def readDomainShapes(self, shapeFile, show=False, noPlane=False):
        ''' Read a domain poly line and add it to the node and element list '''
        # open shape file
        shape = shapefile.Reader(shapeFile)
        # Loop over each feature in shape file
        fig, ax = plt.subplots(1, 1)
        priorPts, allLines = [], []
        for feature in shape.shapeRecords():
            poly = feature.shape.__geo_interface__  # extract polyline feature
            xyPoly = np.array(poly['coordinates'])
            myPts = self.addPoints(xyPoly)  # Add points to geometry
            myLines = self.addLines(myPts, priorPts)  # Add lines to geometry
            allLines.extend(myLines)  # Build list of all lines
            self.geometry.add_physical(myLines)  # Add physical line
            priorPts.extend(myPts)  # Save points
            # Debug
            if show:
                x, y = zip(*list(xyPoly))
                plt.plot(x, y)
                for x1, y1, i in zip(x, y, range(len(x))):
                    plt.text(x1, y1, str(i))
        if show:
            plt.show()
        # Create lineLoop and plane
        lineLoop = self.geometry.add_line_loop(allLines)
        if not noPlane:
            planeSurface = self.geometry.add_plane_surface(lineLoop)
            self.geometry.add_physical(planeSurface)

    def writeGeometry(self, geoFile):
        ''' write to geojson file'''
        with open(geoFile, 'w') as fpGeo:
            fpGeo.write(self.geometry.get_code())
