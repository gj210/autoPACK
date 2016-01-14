# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 16:34:28 2016

@author: ludo
"""
import numpy as np
from autopack.Ingredient import rotVectToVect,ApplyMatrix

class lipidsCG:
    def __init__(self):
        #bilayer is around 50A thick, coverage of beads around 12.5
        self.radii = [6.8,6.8]
        self.pcpalAxis=np.array([1,0,0])
        self.beadsoffset=[9,25]
        self.beads = np.array([[-25,0,0],[-9,0,0],[9,0,0],[25,0,0]])
        self.arrayOfBeads=[]
    def coverShape(self, vertices,normals):       
        self.arrayOfBeads=[]
        #foreach vertices place a outer beads align to normal and inner beads align to -normal
        for i in range(len(vertices)):
            v=vertices[i]
            vn=normals[i]
            rotMat = np.array( rotVectToVect(self.pcpalAxis, vn ), 'f')
            #apply rot mat to beadposition
            new_beads = ApplyMatrix(rotMat,self.beads).tolist()
            self.arrayOfBeads.append(new_beads)
            
            