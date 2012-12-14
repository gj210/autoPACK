# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 23:39:56 2012

@author: -
"""
from math import sqrt
import sys
sys.path.append(".")
sys.path.append("/Users/ludo/DEV/MGLTOOLS/mgl32/MGLToolsPckgs/")
from bhtree import bhtreelib
import AutoFillClean    
from AutoFillClean.pandautil import PandaUtil
pud = PandaUtil()
from AutoFillClean.HistoVol import Grid  
from AutoFillClean.ray import vlen, vdiff, vcross

def tetrahedron(radius):
    """
    Create the mesh data of a tetrahedron of a given radius
    
    @type  radius: float
    @param radius: radius of the embeding sphere
    
    @rtype:   array
    @return:  vertex,face, face normal of the tetrahedron        
    """
    
    faceIndices = ( (0,2,3), (3,2,1), #bot
                    (3, 1, 4), 
                    (3, 4, 0), 
                    (2, 4, 1), 
                    (2, 0, 4), 
                    )
    a = 1 / sqrt(3)
    faceNormals = ( ( a,  a, -a),
                    (-a,  a, -a), 
                    ( a, -a, -a),
                    (-a, -a, -a), )        
    diameter = radius * 2
    width  = diameter
    height = diameter
    depth  = diameter     
    boundingRadius = radius
    surfaceArea = (sqrt(3) / 2) * (radius ** 2)
    volume = (4/3) * (radius ** 3)
    _corners = (
        (-radius, 0.0, 0.0), (radius, 0.0, 0.0),
        (0.0, radius, 0.0), (0.0, -radius, 0.0),
        (0.0, 0.0, radius), )#(0.0, 0.0, radius))
    return _corners,faceIndices,faceNormals
spacing = 5.0      
grid = Grid()
boundingBox = grid.boundingBox = [[-10.0,-10.0,-10.],[10.0,10.0,10.]]
grid.gridSpacing = spacing# = self.smallestProteinSize*1.1547  # 2/sqrt(3)????
#t=time()
#helper.progressBar(label="BuildGRid")
grid.gridVolume,grid.nbGridPoints = grid.computeGridNumberOfPoint(boundingBox,spacing)
grid.create3DPointLookup()
nbPoints = grid.gridVolume
grid.gridPtId = [0]*nbPoints
xl,yl,zl = boundingBox[0]
xr,yr,zr = boundingBox[1]
# distToClosestSurf is set to self.diag initially
grid.diag = diag = vlen( vdiff((xr,yr,zr), (xl,yl,zl) ) )
#grid.distToClosestSurf = [diag]*nbPoints        
#distances = grid.distToClosestSurf
#idarray = grid.gridPtId
#diag = grid.diag
grdPos = grid.masterGridPositions
insidePoints = []
#surfacePoints = vertices
NPT=len(grdPos)
rads = [spacing,]*NPT
iPtList=[]
#helper.progressBar(label="BuildWorldAndNode")
#t=time()
#def addSphere(r,pos,i):
#    node = pud.addSingleSphereRB(r,name=str(i))
#    node.setPos(pos[0],pos[1],pos[2])
#    return node
#nodes =[addSphere(rads[i],grdPos[i],i) for i in range(NPT)  ]
#node = pud.addMultiSphereRB(rads,grdPos)
#print ("time sphere ",time()-t)
#t=time()
#add the mesh
vertices, faces, fn = tetrahedron(10.0)
meshnode = pud.addMeshRB(vertices, faces)
#print ("time mesh",time()-t)
#computeCollisionTest
#t=time()     
srfPts=vertices
from panda3d.core import Mat4,Vec3,Point3
bht =  bhtreelib.BHtree(tuple(vertices), None, 10)
returnNullIfFail = 0
closest = bht.closestPointsArray(tuple(grdPos), diag, returnNullIfFail)
hit=[]
for i in range(NPT):
    start = Point3(grdPos[i][0],grdPos[i][1],grdPos[i][2])
    end = Point3(srfPts[closest[i]][0]*diag,srfPts[closest[i]][1]*diag,srfPts[closest[i]][2]*diag)
    res = pud.world.rayTestAll(start, end)
    if res.hasHits() :
        print res.getNumHits()
        hit.append(res)
#meshcontacts = pud.world.contactTest(meshnode.node())
#N=meshcontacts.getNumContacts()
#for ct in meshcontacts.getContacts():
#    m=ct.getManifoldPoint ()
#    d = m.getDistance ()
#    print ct.getNode0().getName(), ct.getNode1().getName(),d
#    i = eval(ct.getNode0().getName())
#    if i not in iPtList :
#        insidePoints.append(grdPos[i])
#        iPtList.append(i)
#
##meshcontact = pud.world.contactTestPair(meshnode.node(), node.node())
#r=[ pud.world.contactTestPair(meshnode.node(), n.node()) for n in nodes]  
#print len(insidePoints),NPT
##print ("N",N,NPT)
#print ("time contact",time()-t)
#return insidePoints, surfacePoints
