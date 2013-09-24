# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 11:35:38 2013

@author: ludo
"""

#remove all colliding atoms/residues with given geometries.
#two options : kdtree/raycast
import sys
#sys.path.append("/Users/ludo/anaconda/envs/p2.6/lib/python2.6/site-packages/")
sys.path.append("/Users/ludo/anaconda/lib/python2.7/site-packages/")
import numpy
import scipy
import math
#import upy
#helper = upy.getHelperClass()()
def displaySphere(coords,rotate=False):
    helper = self.autoPACK.apgui.helper
    from DejaVu.Spheres import Spheres
    spheres = Spheres("lipids",centers=coords,visible=1,inheritMaterial = False,)
    self.GUI.VIEWER.AddObject(spheres)    
    if rotate : helper.rotateObj(spheres,[math.pi/2.0,0,0])
    return spheres

def displaySphereCenter(center):
    helper = self.autoPACK.apgui.helper
    from DejaVu.Spheres import Spheres
    spheres = Spheres("center",centers=center,radii=[10.,]*len(center),visible=1,inheritMaterial = False,)
    self.GUI.VIEWER.AddObject(spheres)    
    #helper.rotateObj(spheres,[math.pi/2.0,0,0])
    return spheres
 
   
def collide_panda(coordinates,mesh):
    #testcollisio between sphere and mesh rigid_body
    collide  = False
    return collide
    
def collide_rapid(coordinates,mesh):
    #testcollisio between sphere mesh and mesh body
    collide  = False
    return collide
    
def collide_raycast(coordinates,mesh):  
    #collide when distance < atom.radius, or inside
    collide  = False
    return collide
    
def collide_kdtree(coordinates,radii,mesh):
    #collide when distance < atom.radius
    collide  = False
    return collide

def collide_bhtree(residues,matrix,mesh):
    #object tree
    from bhtree import bhtreelib
    v=helper.ApplyMatrix(mesh.getVertices(),matrix.reshape(4,4))
    print len(v),v[0],len(coordinates) 
    tree = bhtreelib.BHtree(tuple(v), None, 10)
    #lipids positon
    c=[]
    at=[]
    res=[]
    for r in residues :
        for a in r :
            c.append(a.get_coord())
            at.append(a)
            res.append(r)
    closest = tree.closestPointsArray(numpy.array(c).tolist(), 30, 0)
    indices=numpy.nonzero(numpy.greater(closest,0))
    indices_close = numpy.take(closest,indices)    
#    collide  = False   
    return indices_close[0],at,res

def getTheMol(filename):
#    from MolKit.protein import Protein
    from MolKit.pdbParser import PdbParser    
    parser = PdbParser(filename=filename)
    mol = parser.parse()
#    mol.read( filename, PdbParser() )
    return mol    
    
def getStructure(name,filename):
    #faster in shell ?
    from Bio.PDB import PDBParser
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(name,filename)
    return structure
    
def getClosestKDtree(coords,pos):
    from scipy.spatial import cKDTree
    tree = cKDTree(coords)
    index = tree.query_ball_point(pos,77.0)
    return index
    
def getClosestBHtree(coords,mpos):
    from bhtree import bhtreelib
    tree=bhtreelib.BHtree( coords, None, 10)
    index=[]
    for i in range(len(mpos)):
        pos = mpos[i]
        result = numpy.zeros(len(coords) ).astype('i')
        nb = tree.closePoints((pos[0],pos[1],pos[2]), 77.0, result )
        index.append(result[:nb])
    #index = [result[:tree.closePoints((pos[0],pos[1],pos[2]), 77.0, result )] for pos in mpos]
    return index
    
def getResidueToCheck(mol,coords, pos):
    remove=[]
    k=0
    index = getClosestBHtree(coords,pos)
#    print (len(index)," atoms found close ")
    residues = [a.parent for a in mol.get_atoms()]
#    rtoremove=[]
#    for i in range(len(index)) : rtoremove.append(numpy.take(residues,index[i]))
    c=[]    
    close_residues=[list(set(numpy.take(residues,index[i]))) for i in range(len(index))]
    c.extend(close_residues)    
    return c,index,residues

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def keep(r,residues_list):
    for i in range(len(residues_list)):
       if (r in residues_list[i]) :
           return False
    return True
    
def removeAround(mol,remove):
    coords2=[]
    toremove = []
    for r in mol.get_residues():
        c=[a.get_coord() for a in r]
        if keep(r,remove) : 
            coords2.extend(c)
        else :
            toremove.extend(c)
    c2=unique_rows(coords2)
    #transform
    mat  = helper.eulerToMatrix([-math.pi/2.0,0,0])
    c2 = helper.ApplyMatrix(c2,mat)    
    displaySphere(c2)
    return c2,toremove


#mpos=[ [-13.18,49.488,-13.182],[-13.17383,49.54086,20.8162],[223.684066,8.164540,-210.3868255], [421.026733398,-43.89,-36.7927],]
#object/kdtree ?
#rotation ?
#s=parser.get_structure("lipid","/Users/ludo/DEV/plane_noH.pdb")
mol = getStructure("lipids","/Users/ludo/DEV/plane_noH.pdb")
#mol = self.Mols[0]
coords = [a.get_coord() for a in mol.get_atoms()]
#need to ApplyMatrix
helper = self.autoPACK.apgui.helper
mat  = helper.eulerToMatrix([-math.pi/2.0,0,0])
rcoords = helper.ApplyMatrix(coords,mat)
#pc = helper.PointCloudObject("cloud",vertices=coords)
h=self.autoPACK.apgui.histoVol.values()[0]
liste_molecules = h.organelles[0].surfaceRecipe.ingredients
#DejaVu
meshes = [m.mesh.children[0] for m in liste_molecules]
transformations = [m.mesh.instanceMatricesFortran for m in liste_molecules]
trans=[transformations[1][0].reshape(4,4),transformations[0][2].reshape(4,4),transformations[0][1].reshape(4,4),transformations[0][0].reshape(4,4)]
meshs=[meshes[1],meshes[0],meshes[0],meshes[0]]
mpos = [tr[3,:3].tolist() for tr in trans]
remove,index,allresidues = getResidueToCheck(mol,rcoords.tolist(),mpos)
c2,toremove=removeAround(mol,remove)
r=[]
rove=[]
for i in range(len(h.organelles[0].molecules)): 
    indice_to_remove,at,res = collide_bhtree(remove[i],trans[i],meshs[i])
    #exclueat = numpy.array(at)[indice_to_remove]
    rove.append(list(set(res)))
    r.append([indice_to_remove,at])
c2,toremove=removeAround(mol,rove)