# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 14:33:23 2018

@author: ian
"""
from modelfunc import myerror
import matplotlib.pyplot as plt
import numpy as np

elSize = {1: 2, 2: 3, 3: 4}


class gMesh:
    # nodes
    def __init__(self, expFile=None, meshType='Triangle'):
        self.nNodes, self.nElements = 0, 0
        self.nodeValues = []
        self.nodes, self.elements = [], []
        self.nNVal, self.nEVal = 0, 0
        self.nodeIndex, self.elIndex = [], []
        self.meshpy = None
        self.elementType = []
        self.tag, self.group = [], []
        self.meshtype = meshType
        if expFile is not None:
            self.readExp(expFile)

    def readGmshHeader(self, gmshFile):
        ''' read gmsh header, return pointer to start of nodes'''
        try:
            fpGmsh = open(gmshFile, 'r')
        except Exception:
            myerror(f'error opening gmsh input file {gmshFile}')
        for line in fpGmsh:
            if '$Nodes' in line:
                return fpGmsh
        myerror(f'Reached end of {gmshFile} with no $Nodes found')

    def readGmshNodes(self, fpGmsh):
        ''' read nodes from a gmsh file '''
        try:
            nNewNodes = int(fpGmsh.readline().strip())
        except Exception:
            myerror('readGmshNodes: cannot parse nLines')
        #
        nodeIndex, nodesRead = len(self.nodeIndex), 0
        for line in fpGmsh:
            nodeData = [eval(x) for x in line.split()]
            self.nodes.append(nodeData[1:3])
            self.nodeValues.append(nodeData[3])
            nodeIndex += 1  # Increment node index and append
            self.nodeIndex.append(nodeIndex)
            nodesRead += 1  # Increment node counter
            if nodesRead == nNewNodes:
                print(f'Read {nodesRead} nodes')
                self.nNodes += nodesRead
                break
        self.nodes = np.array(self.nodes)
        print(self.nodes.shape)
        if nodesRead != nNewNodes:
            myerror(f'readGmshNodes: expected {nNewNodes} nodes read '
                    f'{nodesRead}')
        if '$EndNodes' not in fpGmsh.readline():
            myerror(f'readGmshNodes: expected $EendNodes')

    def readGmshElements(self, fpGmsh):
        ''' read elements from a gmsh file '''
        try:
            print(f'{fpGmsh.readline()}')
            nNewElements = int(fpGmsh.readline().strip())
        except Exception:
            myerror('readGmshElements: cannot parse nLines')
        #
        elIndex, elRead = len(self.elIndex), 0
        for line in fpGmsh:
            elData = [eval(x) for x in line.split()]
            self.elements.append(elData[5:])
            self.elementType.append(elData[1])
            self.tag.append(elData[3])
            self.group.append(elData[4])
            elIndex += 1  # Increment el index and append
            self.elIndex.append(elIndex)
            elRead += 1  # Increment el counter
            if elRead == nNewElements:
                print(f'Read {elRead} elements')
                self.nElements += elRead
                break
        if elRead != nNewElements:
            myerror(f'readGmshElements: expected {nNewElements} nodes read '
                    f'{elRead}')
        if '$EndElements' not in fpGmsh.readline():
            myerror(f'readGmshNodes: expected $EndElements')

    def readGmsh(self, gmshFile):
        '''Read gmsh file'''
        fpGmsh = self.readGmshHeader(gmshFile)
        self.readGmshNodes(fpGmsh)
        self.readGmshElements(fpGmsh)

    def getNewTag(self):
        ''' Get a new incremental tag '''
        if len(self.tag) == 0:
            return 1  # First tag
        else:
            return self.tag[-1] + 1  # Increment last tag

    def gmshHeader(self, gmshFile):
        ''' Open gmsh file and write header'''
        # open file
        try:
            fpGmsh = open(gmshFile, 'w')
        except Exception:
            myerror(f'error opening gmsh output file {gmshFile}')
        # Header info
        print('$MeshFormat', file=fpGmsh)
        print('2.2 0 8', file=fpGmsh)
        print('$EndMeshFormat', file=fpGmsh)
        return fpGmsh

    def writeGmshNodes(self, fpGmsh):
        ''' write node data '''
        print('$Nodes', file=fpGmsh)
        print(self.nNodes, file=fpGmsh)
        #
        for i in range(self.nNodes):
            print(f'{self.nodeIndex[i]:10} {self.nodes[i][0]:15.5f} '
                  f'{self.nodes[i][1]:15.5f} 0.0 ', file=fpGmsh)
        print('$EndNodes', file=fpGmsh)

    def writeGmshElements(self, fpGmsh):
        ''' write node data '''
        print('$Elements', file=fpGmsh)
        print(self.nElements, file=fpGmsh)
        for i in range(self.nElements):
            # [rint(self.elements[i])
            if self.elementType[i] > 1:
                group = 1000
            else:
                group = self.elIndex[i]
            print(f'{self.elIndex[i]:8} {self.elementType[i]:3} 2 '
                  f'{self.tag[i]:4} {self.group[i]:6} '
                  f'{"".join([f"{x:10}" for x in self.elements[i]])}',
                  file=fpGmsh)
        print('$EndElements', file=fpGmsh)

    def toGmsh(self, gmshFile):
        ''' Write msh in gmsh format '''
        #
        # init
        fpGmsh = self.gmshHeader(gmshFile)
        # Nodes
        self.writeGmshNodes(fpGmsh)
        #
        # Elements
        self.writeGmshElements(fpGmsh)
        #
        fpGmsh.close()

    def plotMesh(self, nodes=True, elements=True):
        ''' plot mesh  - very slow mostly for debugging'''
        fig, ax = plt.subplots(1, 1)
        #
        i = 0
        ax.plot(self.nodes[:,0], self.nodes[:,1], 'k.')
        for myEl in self.elements:
            #print(self.nodes[myEl,:])
            kk=np.array(myEl)-1
            xy = self.nodes[kk, :]
            #print(xy.shape)
            ax.plot(xy[:,0], xy[:,1],'r-')
            i+=1
            if i == 1000 :
                 break
            #     print(myEl)
            #
            #     print(xy[:,0], xy[:,1])
        plt.show()