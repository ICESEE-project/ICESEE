# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 14:33:23 2018

@author: ian
"""

import matplotlib.pyplot as plt
import numpy as np
from modelfunc.myerror import myerror


class argusMesh:
    ''' Class to read argusMeshes and then plot, or convert them to gmsh
    format. Currently it only works with with a single domain outline.  Will
    increment domain label by 100 to indicate a nodes flagged as shelf front.

    Need to define a set of labels in Argus for domains with multiple
    polygons. '''

    def __init__(self, expFile=None, meshType='Triangle'):
        self.nNodes, self.nElements = 0, 0
        self.nodes, self.elements = [[0, 0]], [[0, 0, 0]]  # 0 index not used
        self.nodeValues, self.elementValues = [[0]*4], [[0]*3]
        self.nNVal, self.nEVal = 0, 0
        self.nodeIndex, self.elIndex = [0], [0]
        self.nodeOnBoundary, self.elOnBoundary = [False], [False]
        self.shelfFrontNode = [False]
        self.segments = []
        self.nodeType = [0]
        if expFile is not None:
            self.readExp(expFile)

    def parseHeader(self, fpExp):
        ''' Parse header info from Arfua exp file. '''
        headerPieces = fpExp.readline().split()
        if len(headerPieces) != 4:
            myerror('header not found')
        # parse header
        nNewNodes = int(headerPieces[1])  # Number of nodes
        nNewElements = int(headerPieces[0])  # Number of elements
        self.nNVal = int(headerPieces[3])  # Number of cols of node data
        self.nEVal = int(headerPieces[2])  # Number of cols of el data
        self.nodeValues = [[0]*self.nNVal]  # Zero unused index 0 node vals
        self.elementValues = [[0]*self.nEVal]  # Zero unused index 0 el vals
        print(f'Header info {nNewNodes} {nNewElements} '
              f'{self.nNVal} {self.nEVal}')
        return nNewNodes, nNewElements

    def parseNodes(self, fpExp, nNewNodes):
        ''' Parse nodal values from argus exp file'''
        nodesRead = 0
        nodeIndex = len(self.nodeIndex)  # WIth zero filled, will start with 1
        print('parseNodes')
        for line in fpExp:
            if line.isspace():  # Skip blank lines
                continue
            pieces = line.split()
            if pieces[0] == 'N':
                # Append coordinates for a point
                self.nodes.append((float(pieces[2]), float(pieces[3])))
                # parse node data values
                self.nodeOnBoundary.append(bool(int(pieces[4])))
                self.shelfFrontNode.append(bool(int(pieces[6])))
                self.nodeValues.append([float(pieces[4+i])
                                        for i in range(self.nNVal)])
                self.nodeIndex.append(nodeIndex)  # Inct node index and append
                # I need to add a type field for multiple domain boundaries
                if len(pieces) <= 8:
                    self.nodeType.append(1)
                else:
                    self.nodeType.append(int(pieces[8]))
                nodeIndex += 1
                nodesRead += 1  # Increment node counter
                if nodesRead == nNewNodes:
                    self.nNodes += nodesRead
                    self.nodes = np.array(self.nodes)
                    self.nodeOnBoundary = np.array(self.nodeOnBoundary)
                    self.nodeType = np.array(self.nodeType)
                    # print(self.nodeType.shape, self.nodeOnBoundary.shape)
                    break

        self.uniqueNodeTypes = np.unique(self.nodeType)[1:]
        # print(self.uniqueNodeTypes)
        if nodesRead != self.nNodes:
            myerror(f'Expected {nNewNodes} nodes but only read {nodesRead}')
        self.nodeIndex = np.array(self.nodeIndex)
        return nodesRead

    def parseElements(self, fpExp, nNewElements):
        ''' Parse Elements from an Argus export file. '''
        # Loop through line
        elementsRead = 0
        elementIndex = len(self.elIndex)
        for line in fpExp:
            pieces = line.split()
            if pieces[0] == 'E':
                # Read indices with node increment if needed
                self.elements.append([int(pieces[i+2]) for i in range(3)])
                self.elementValues.append([int(pieces[i+5])
                                           for i in range(self.nEVal)])
                self.elIndex.append(elementIndex)
                self.elOnBoundary.append(bool(int(pieces[5])))
                elementsRead += 1
                elementIndex += 1
                if elementsRead == nNewElements:
                    self.nElements += elementsRead
                    self.elements = np.array(self.elements)
                    return elementsRead
        #
        if elementsRead != self.nEVal:
            myerror(
                f'Expected {nNewElements} elements but found {elementsRead}')

    def readExp(self, expFile):
        ''' Read Argus export (.exp) file. '''
        try:  # open file
            fpExp = open(expFile, 'r')
        except Exception:
            myerror(f'argusMesh.readExp: could not open meshfile {expFile} '
                    'for read')
        #
        nNewNodes, nNewElements = self.parseHeader(fpExp)
        nodesRead = self.parseNodes(fpExp, nNewNodes)
        elementsRead = self.parseElements(fpExp, nNewElements)
        print(nodesRead, self.nNodes, elementsRead, self.nElements)
        fpExp.close()
        return

    def findNearest(self, n0, bNodes, bNotUsed, nIndex, bIndex):
        ''' find nearest node to last node used. Used to find boundary nodes'''
        # Compute distances to n0 from unused nodes
        d = (bNodes[bNotUsed, 0] - n0[0])**2 + (bNodes[bNotUsed, 1] - n0[1])**2
        i = np.where(d == min(d))[0][0]  # Find closest node
        i0 = nIndex[bNotUsed][i]
        ii = bIndex[bNotUsed][i]  # map to boundary node index
        bNotUsed[ii] = False
        return self.nodeIndex[i0]

    def typeSegment(self, n1, n2):
        ''' Compute type for segment, use type, but increment by 100 for
        floating segments '''
        if self.shelfFrontNode[n1] or self.shelfFrontNode[n2]:
            return self.nodeType[n1] + 100  # Both flt type+100
        else:
            return self.nodeType[n1]  # Both grounded, use type

    def minDistance(self, nodes):
        '''
        For list of node indices, find pair with minimum distances
        '''
        minD = 1e9
        for ni in nodes:
            otherNodes = np.array([x for x in nodes if x != ni])
            d = np.array((self.nodes[otherNodes, 0] - self.nodes[ni, 0])**2 +
                         (self.nodes[otherNodes, 1] - self.nodes[ni, 1])**2)
            #
            if np.min(d) < minD:
                n1 = ni
                n2 = otherNodes[np.where(d == min(d))[0][0]]
                minD = np.min(d)
        return n1, n2

    def boundarySegments(self):
        self.segments = []
        self.timesUsed = np.zeros(self.nodes.shape[0])
        for segType in self.uniqueNodeTypes:
            self.boundarySegmentType(segType)
        #
        # Now close up contours where nodes have not been used in 2 segments
        while True:
            usedOnce = self.nodeIndex[self.timesUsed == 1]
            # Quit when run out of usedOnce
            if len(usedOnce) < 2:
                break
            n1, n2 = self.minDistance(usedOnce)
            self.timesUsed[n1] += 1
            self.timesUsed[n2] += 1
            segType = self.typeSegment(n1, n2)
            self.segments.append({'n1': n1, 'n2': n2, 'type': segType})
        #
        self.segments = sorted(self.segments, key=lambda d: d['type'])

    def boundarySegmentType(self, segType):
        ''' find line segments on boundaries. Need to modify to handle
        multiple  polygons'''
        # Boundary nodes
        segTypeNodes = np.logical_and(self.nodeType == segType,
                                      self.nodeOnBoundary)

        bNodes = np.copy(self.nodes[segTypeNodes])
        nIndex = np.array(self.nodeIndex)[segTypeNodes]  # Boundary ind
        bIndex = np.arange(0, bNodes.shape[0])  # Index into nIndex
        bNotUsed = np.ones(bNodes.shape[0], dtype=bool)  # Nodes not yet used
        bNotUsed[0] = False  # Start with first index
        n0 = nIndex[0]
        nLast = n0
        # Loop to find nearest node and form segments
        # (np.sum(bNotUsed))
        while np.sum(bNotUsed) > 0:
            n0 = self.findNearest(self.nodes[n0], bNodes, bNotUsed, nIndex,
                                  bIndex)
            segType = self.typeSegment(nLast, n0)
            self.segments.append({'n1': nLast, 'n2':  n0, 'type': segType})
            # Track number of times each node used in segment
            self.timesUsed[nLast] += 1
            self.timesUsed[n0] += 1
            nLast = n0

    def plotSegments(self, ax, segType=None):
        ''' Plot boundary segments '''
        self.computeIDColors()
        #
        self.boundarySegments()
        #
        for segment in self.segments:
            node1 = self.nodes[segment['n1']]
            node2 = self.nodes[segment['n2']]
            if segType == None or segment["type"] == segType:
                ax.plot([node1[0], node2[0]], [node1[1], node2[1]], '-',
                        color=self.cTable[f'{segment["type"]}'])

    def meshLim(self):
        ''' Return minx,maxx,miny,maxy of domain outline '''
        xmin, xmax = np.min(self.nodes[1:, 0]), np.max(self.nodes[1:, 0])
        ymin, ymax = np.min(self.nodes[1:, 1]), np.max(self.nodes[1:, 1])
        return xmin, xmax, ymin, ymax

    def computeIDColors(self):
        nUnique = len(self.uniqueNodeTypes)
        myColors = plt.cm.get_cmap('jet', nUnique*2)
        self.cTable = {'0': (0, 0, 0, 0)}
        for x, y in zip(self.uniqueNodeTypes, range(0, nUnique)):
            self.cTable[str(x)] = myColors(y)
            self.cTable[str(x+100)] = myColors(y + nUnique)
        #
        self.IDColors = np.array([self.cTable[f'{x}'] for x in self.nodeType])

    def plotMesh(self, nodes=True, elements=True, nodeOnBoundary=True,
                 shelfFrontNode=True, ax=None):
        ''' plot mesh  - very slow mostly for debugging'''
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        # Plot nodes

        if nodeOnBoundary:
            self.computeIDColors()
            for myID in self.uniqueNodeTypes:
                toPlot = np.logical_and(self.nodeOnBoundary,
                                        self.nodeType == myID)
                ax.scatter(self.nodes[toPlot, 0],
                           self.nodes[toPlot, 1], c=self.IDColors[toPlot],
                           marker='o', label=myID, s=5)
        if shelfFrontNode:
            ax.plot(self.nodes[self.shelfFrontNode, 0],
                    self.nodes[self.shelfFrontNode, 1], '.', color='m',
                    label='shelfFrontNode')
        if nodes:
            ax.plot(self.nodes[1:, 0], self.nodes[1:, 1], 'k.',
                    label='Nodes', markersize=2)
        label = 'ElOnBoundary'
        if elements:
            for myEl, onBoundary in zip(self.elements, self.elOnBoundary):
                if onBoundary:
                    xy = self.nodes[myEl, :]
                    ax.plot(xy[:, 0], xy[:, 1], 'g-', label=label)
                    label = None
        ax.legend()

    def gmshHeader(self, gmshFile):
        ''' Open gmsh file and write header'''
        try:
            fpGmsh = open(gmshFile, 'w')  # open file
        except Exception:
            myerror(f'error opening gmsh output file {gmshFile}')
        # Header info
        print('$MeshFormat', file=fpGmsh)
        print('2.2 0 8', file=fpGmsh)
        print('$EndMeshFormat', file=fpGmsh)
        return fpGmsh

    def writeGmshSegments(self, fpGmsh):
        ''' write segments as 1-d nodes'''
        # print(len(self.segments))
        nSegments = len(self.segments)
        for seg, i in zip(self.segments, range(1, nSegments+1)):
            print(f'{i:8} {1:3} {2:3} {seg["type"]:4} {1:6} '
                f'{seg["n1"]:6} {seg["n2"]:6}',
                file=fpGmsh)
            #print(f'{i:8} {1:3} {2:3} {0:3} {seg["type"]:4}  '
            #      f'{seg["n1"]:6} {seg["n2"]:6}',
                  
        return nSegments

    def writeGmshNodes(self, fpGmsh):
        ''' write node data '''
        print('$Nodes', file=fpGmsh)
        print(self.nNodes, file=fpGmsh)
        #
        for i in range(1, self.nNodes+1):
            print(f'{self.nodeIndex[i]:10} {self.nodes[i][0]:15.5f} '
                  f'{self.nodes[i][1]:15.5f} 0.0 ', file=fpGmsh)
        print('$EndNodes', file=fpGmsh)

    def writeGmshElements(self, fpGmsh):
        ''' write node data '''
        print('$Elements', file=fpGmsh)
        nSegs = len(self.segments)
        nElements = self.nElements + nSegs
        print(nElements, file=fpGmsh)
        nSegments = self.writeGmshSegments(fpGmsh)
        #print(nSegments, self.nElements, self.nElements+nSegments)
        #
        for i in range(1, self.nElements+1):
            # [rint(self.elements[i])
            print(f'{self.elIndex[i] + nSegments:8} {2:3} 2 0 99 '
                  f'{"".join([f"{x:10}" for x in self.elements[i]])}',
                  file=fpGmsh)
        print('$EndElements', file=fpGmsh)

    def toGmsh(self, gmshFile):
        ''' Write msh in gmsh format '''
        #
        self.boundarySegments()
        # init
        fpGmsh = self.gmshHeader(gmshFile)
        # Nodes
        self.writeGmshNodes(fpGmsh)
        #
        # Elements
        self.writeGmshElements(fpGmsh)
        #
        fpGmsh.close()
