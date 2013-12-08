# -*- coding: utf-8 -*-
"""
###############################################################################
#
# autoPACK Authors: Graham T. Johnson, Mostafa Al-Alusi, Ludovic Autin, Michel Sanner
#   Based on COFFEE Script developed by Graham Johnson between 2005 and 2010 
#   with assistance from Mostafa Al-Alusi in 2009 and periodic input 
#   from Arthur Olson's Molecular Graphics Lab
#
# HistoVol.py Authors: Graham Johnson & Michel Sanner with editing/enhancement from Ludovic Autin
#
# Translation to Python initiated March 1, 2010 by Michel Sanner with Graham Johnson
#
# Class restructuring and organization: Michel Sanner
#
# Copyright: Graham Johnson ©2010
#
# This file "HistoVol.py" is part of autoPACK, cellPACK.
#    
#    autoPACK is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    autoPACK is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with autoPACK (See "CopyingGNUGPL" in the installation.
#    If not, see <http://www.gnu.org/licenses/>.
#
#
###############################################################################
@author: Ludovic Autin, Graham Johnson,  & Michel Sanner
"""
import numpy
from autopack.ldSequence import cHaltonSequence3
from scipy import spatial
from math import ceil
from math import floor

class Grid:
    """
    The Grid class
    ==========================
    This class handle the use of grid to control the packing. The grid keep information of
    3d positions, distances, freePoints and inside/surface points from organelles.
    NOTE : thi class could be completly replace if openvdb is wrapped to python.
    """
    def __init__(self,boundingBox=([0,0,0], [.1,.1,.1]),space=1,setup=True):
        #a grid is attached to an environement
        self.boundingBox=boundingBox
        # this list provides the id of the component this grid points belongs
        # to. The id is an integer where 0 is the Histological Volume, and +i is
        # the surface of compartment i and -i is the interior of compartment i
        # in the list self. compartments
        self.gridPtId = []
        # will be a list of indices into 3D of histovol
        # of points that have not yet been used by the fill algorithm
        # entries are removed from this list as grid points are used up
        # during hte fill. This list is used to pick points randomly during
        # the fill
        self.freePoints = []
        self.nbFreePoints = 0
        # this list evolves in parallel with self.freePoints and provides
        # the distance to the closest surface (either an already placed
        # object (or an compartment surface NOT IMPLEMENTED)
        self.distToClosestSurf = []
        self.distToClosestSurf_store = []
        
        self.diag=None
        self.gridSpacing = space*1.1547
        self.nbGridPoints = None
        self.nbSurfacePoints = 0
        self.gridVolume = 0 # will be the toatl number of grid points
        # list of (x,y,z) for each grid point (x index moving fastest)
        self.masterGridPositions = []
        self._x = None
        self._y = None
        self._z = None
        
        #this are specific for each compartment
        self.aInteriorGrids = []
        self.aSurfaceGrids = []
        #bhtree
        self.surfPtsBht=None
        self.ijkPtIndice = []
        self.filename=None          #used for storing before fill so no need rebuild
        self.result_filename=None   #used after fill to store result
        self.tree = None
        self.encapsulatingGrid = 1
        if setup :
            self.setup(boundingBox,space)
        #use np.roll to have periodic condition
        #what about collision ?
            
    def setup(self,boundingBox,space):
        self.gridSpacing = space*1.1547
        self.boundingBox = boundingBox
        #self.gridVolume,self.nbGridPoints=self.computeGridNumberOfPoint(boundingBox,self.gridSpacing)
        self.create3DPointLookup()
#        self.create3DPointLookupCover()
        self.getDiagonal()
        self.nbSurfacePoints = 0
        print ("$$$$$$", self.gridVolume,self.gridSpacing)
        self.gridPtId = numpy.zeros(self.gridVolume,'i')#[0]*nbPoints
        #self.distToClosestSurf = [self.diag]*self.gridVolume#surface point too?
        self.distToClosestSurf = numpy.ones(self.gridVolume)*self.diag#(self.distToClosestSurf)
        self.freePoints = list(range(self.gridVolume))
        self.nbFreePoints =len(self.freePoints)
        print ("$$$$$$$$  1 ",boundingBox,space,self.gridSpacing,len(self.gridPtId))
#        self.create3DPointLookup()#whatdoI do it twice ?
        print("$$$$$$$$  gridVolume = nbPoints = ", self.gridVolume, 
                          " grid.nbGridPoints = ", self.nbGridPoints,
                          "gridPId = ",len(self.gridPtId),
                        "self.nbFreePoints =",self.nbFreePoints)

    def reset(self,):
        #reset the  distToClosestSurf and the freePoints
        #boundingBox shoud be the same otherwise why keeping the grid
        self.distToClosestSurf = numpy.array([self.diag]*len(self.distToClosestSurf))#surface point too?
        self.freePoints = list(range(len(self.freePoints)))
        self.nbFreePoints =len(self.freePoints)
        
    def removeFreePoint(self,pti):
        tmp = self.freePoints[self.nbFreePoints] #last one
        self.freePoints[self.nbFreePoints] = pti
        self.freePoints[pti] = tmp
        self.nbFreePoints -= 1        

# Very dangerous to manipulate the grids... lets solve this problem much earlier in the setup with the new PseudoCode
#    def updateDistances(self, histoVol ,insidePoints, freePoints,
#                        nbFreePoints ):
#        verbose = histoVol.verbose
#        nbPts = len(insidePoints)
#        for pt in insidePoints:  #Reversing is not necessary if you use the correct Swapping GJ Aug 17,2012
#            try :
#                # New system replaced by Graham on Aug 18, 2012
#                nbFreePoints -= 1  
#                vKill = freePoints[pt]
#                vLastFree = freePoints[nbFreePoints]
#                freePoints[vKill] = vLastFree
#                freePoints[vLastFree] = vKill
#            except :
#                pass 
#            
#        return nbFreePoints,freePoints
#
#    def removeFreePointdeque(self,pti):
#        self.freePoints.remove(pti)

    def getDiagonal(self, boundingBox=None):
        if boundingBox is None :
            boundingBox = self.boundingBox 
        self.diag = numpy.sqrt((numpy.array(boundingBox[0])-numpy.array(boundingBox[1]))**2).sum()
        return self.diag
        
    def create3DPointLookup_loop(self, boundingBox=None):
        """
        Fill the orthogonal bounding box described by two global corners
        with an array of points spaces pGridSpacing apart.:
        """
        if boundingBox is None :
            boundingBox= self.boundingBox
        xl,yl,zl = boundingBox[0]
        xr,yr,zr = boundingBox[1]

        nx,ny,nz = self.nbGridPoints
        pointArrayRaw = numpy.zeros( (nx*ny*nz, 3), 'f')
        #self.ijkPtIndice = numpy.zeros( (nx*ny*nz, 3), 'i')#this is unused 
        #try :
        self.ijkPtIndice = numpy.ndindex(nx,ny,nz)
        space = self.gridSpacing
        # Vector for lower left broken into real of only the z coord.
        i = 0
        for zi in xrange(nz):
            for yi in xrange(ny):
                for xi in xrange(nx):
                    pointArrayRaw[i] = (xl+xi*space, yl+yi*space, zl+zi*space)
                    #self.ijkPtIndice[i] = (xi,yi,zi)
                    #print ("add i",i,xi,yi,zi,nx,ny,nz)
                    i+=1
        self.masterGridPositions = pointArrayRaw

    def create3DPointLookup(self, boundingBox=None):
        """
        Fill the orthogonal bounding box described by two global corners
        with an array of points spaces pGridSpacing apart. Optimized version using
        numpy broadcasting
        """
        if boundingBox is None :
            boundingBox= self.boundingBox
        space = self.gridSpacing
        # we want the diagonal of the voxel, not the diagonal of the plane, so the second 1.1547 is was incorrect
        environmentBoxEqualFillBox = True
        #np.linspace(2.0, 3.0, num=5)
        if environmentBoxEqualFillBox: #environment.environmentBoxEqualFillBox:
            self._x = x = numpy.arange(boundingBox[0][0], boundingBox[1][0], space)#*1.1547) gridspacing is already multiplied by 1.1547
            self._y = y = numpy.arange(boundingBox[0][1], boundingBox[1][1], space)#*1.1547)
            self._z = z = numpy.arange(boundingBox[0][2], boundingBox[1][2], space)#*1.1547)
        else:
            self._x = x = numpy.arange(boundingBox[0][0], boundingBox[1][0] + space, space)#*1.1547) gridspacing is already multiplied by 1.1547
            self._y = y = numpy.arange(boundingBox[0][1], boundingBox[1][1] + space, space)#*1.1547)
            self._z = z = numpy.arange(boundingBox[0][2], boundingBox[1][2] + space, space)#*1.1547)
#        self._x = x = numpy.ogrid(boundingBox[0][0]: boundingBox[1][0]: space*1.1547j)
#        self._y = y = numpy.ogrid(boundingBox[0][1]: boundingBox[1][1]: space*1.1547j)
#        self._z = z = numpy.ogrid(boundingBox[0][2]: boundingBox[1][2]: space*1.1547j)
        #nx,ny,nz = self.nbGridPoints
        nx = len(x) # sizes must be +1 or the right, top, and back edges don't get any points using this numpy.arange method
        ny = len(y)
        nz = len(z)
        # Dec 5 2013, we need to confirm that the getPointsInBox function is also using +1, or potential neighbors will be missed
        # This used to be fine, but it may have changed?
        
        self.nbGridPoints = [nx,ny,nz]
        self.gridVolume = nx*ny*nz
        self.ijkPtIndice = numpy.ndindex(nx,ny,nz)
        #this is 60% faster than the for loop
#        self.masterGridPositions = numpy.array(list(numpy.broadcast(*numpy.ix_(x, y, z))))
#        self.masterGridPositions = numpy.vstack(numpy.meshgrid(x,y,z)).reshape(3,-1).T
        self.masterGridPositions = numpy.vstack(numpy.meshgrid(x,y,z,copy=False)).reshape(3,-1).T
        #this ay be faster but dont kno the implication
#        np.vstack((ndmesh(x_p,y_p,z_p,copy=False))).reshape(3,-1).T
        #from http://stackoverflow.com/questions/18253210/creating-a-numpy-array-of-3d-coordinates-from-three-1d-arrays

    def create3DPointLookupCover(self, boundingBox=None):
        """
        Fill the orthogonal bounding box described by two global corners
        with an array of points spaces pGridSpacing apart. Optimized version using
        numpy broadcasting
        """
        if boundingBox is None :
            boundingBox= self.boundingBox
        space = self.gridSpacing
        S =  numpy.array(boundingBox[1])-numpy.array(boundingBox[0])
        NX,NY,NZ = numpy.around(S/self.gridSpacing)
        if NX == 0 : NX=1
        if NY == 0 : NY=1
        if NZ == 0 : NZ=1
        print NX,NY,NZ
        # we want the diagonal of the voxel, not the diagonal of the plane, so the second 1.1547 is was incorrect
        environmentBoxEqualFillBox = True
        #np.linspace(2.0, 3.0, num=5)
        if environmentBoxEqualFillBox: #environment.environmentBoxEqualFillBox:
            self._x = x = numpy.linspace(boundingBox[0][0], boundingBox[1][0], int(NX))#*1.1547) gridspacing is already multiplied by 1.1547
            self._y = y = numpy.linspace(boundingBox[0][1], boundingBox[1][1], int(NY))#*1.1547)
            self._z = z = numpy.linspace(boundingBox[0][2], boundingBox[1][2], int(NZ))#*1.1547)
        else:
            self._x = x = numpy.arange(boundingBox[0][0], boundingBox[1][0] + space, space)#*1.1547) gridspacing is already multiplied by 1.1547
            self._y = y = numpy.arange(boundingBox[0][1], boundingBox[1][1] + space, space)#*1.1547)
            self._z = z = numpy.arange(boundingBox[0][2], boundingBox[1][2] + space, space)#*1.1547)
        xyz = numpy.meshgrid(x,y,z,copy=False)
        nx = len(x) # sizes must be +1 or the right, top, and back edges don't get any points using this numpy.arange method
        ny = len(y)
        nz = len(z)
        self.gridSpacing = x[1]-x[0]
        self.nbGridPoints = [nx,ny,nz]
        self.gridVolume = nx*ny*nz
        self.ijkPtIndice = numpy.ndindex(nx,ny,nz)
        self.masterGridPositions = numpy.vstack(xyz).reshape(3,-1).T

    def getPointCompartmentId(self,point,ray=1):
        #check if point inside on of the compartments
        #surface point ?
        ncomp = len(self.histoVol.compartments)
        if ncomp:
            comp = ncomp
            for i in range(ncomp):                
                inside = checkPointInside_rapid( self.histoVol.compartments[i],
                                            point,self.histoVol.grid.diag,ray=ray)
                if inside :
                    return -(i+1)
                #comp=comp-1
            #the point is not inside , is it on the surface ? ie distance to surface < X?
            for i in range(ncomp):                
                distance,nb = self.histoVol.compartments[i].OGsrfPtsBht.query(point)
                if distance < 10.0 :
                    return i+1
        return 0

    def getClosestGridPoint(self,pt3d):
        if self.tree is None :
            self.tree = spatial.cKDTree(self.masterGridPositions, leafsize=10)
        distance,nb = self.tree.query(pt3d)#len of ingr posed so far
        return distance,nb

    def getPointFrom3D(self, pt3d):
        """
        get point number from 3d coordinates
        """
        x, y, z = pt3d  # Continuous 3D point to be discretized
        spacing1 = 1./self.gridSpacing  # Grid spacing = diagonal of the voxel determined by smalled packing radius
        NX, NY, NZ = self.nbGridPoints  # vector = [length, height, depth] of grid, units = gridPoints
        OX, OY, OZ = self.boundingBox[0] # origin of fill grid
        # Algebra gives nearest gridPoint ID to pt3D
        i = min( NX-1, max( 0, round((x-OX)*spacing1)))
        j = min( NY-1, max( 0, round((y-OY)*spacing1)))
        k = min( NZ-1, max( 0, round((z-OZ)*spacing1)))
        return int(k*NX*NY + j*NX + i)

    def getIJK(self,ptInd):
        """
        get i,j,k (3d) indices from u (1d)
        """
        #ptInd = k*(sizex)*(sizey)+j*(sizex)+i;#want i,j,k
        return self.ijkPtIndice[ptInd]


    def checkPointInside(self,pt3d,dist=None,jitter=[1,1,1]):
        """
        Check if the given 3d points is inside the grid
        """        
        O = numpy.array(self.boundingBox[0])
        E = numpy.array(self.boundingBox[1])
        P = numpy.array(pt3d)
        test1 = P < O 
        test2 =  P > E
        if True in test1 or True in test2:
            #outside
            print "outside"
            return False
        else :
            if dist is not None:
#                print "ok distance ",dist,P
                #distance to closest wall
                d1 = (P - O)*jitter
                s1=min(x for x in d1[d1 != 0] if x != 0)
#                print d1,s1,s1<=dist
                #s1 = numpy.sum(d1*d1)
                d2 = (E - P)*jitter
                s2=min(x for x in d2[d2 != 0] if x != 0)
#                print d2,s2,s2<=dist
                #s2 = numpy.sum(d2*d2)
                if s1 <= dist or s2 <=dist:
#                    print ("too close")
                    return False 
            return True
        
    def getCenter(self):
        """
        Get the center of the grid
        """        
        center=[0.,0.,0.]
        for i in range(3):
            center[i]=(self.boundingBox[0][i]+self.boundingBox[1][i])/2.
        return center
        
    def getRadius(self):
        """
        Get the radius the grid
        """        
        d = numpy.array(self.boundingBox[0]) - numpy.array(self.boundingBox[1])
        s = numpy.sum(d*d)
        return math.sqrt(s)

    def getPointsInSphere(self, bb, pt, radius,addSP=True,info=False):
        #diag = numpy.sqrt((numpy.array(bb[0])-numpy.array(bb[1]))**2).sum()
        diag = abs(bb[0][0]-bb[1][0])/2.0
        if self.tree is None :
            self.tree = spatial.cKDTree(self.masterGridPositions, leafsize=10)
        return self.tree.query_ball_point(pt,diag)#len of ingr posed so far
  
    def getPointsInCubeFillBB(self, bb, pt, radius,addSP=True,info=False):
        """
        Return all grid points indices inside the given bouding box.
        NOTE : need to fix with grid build with numpy arrange
        """        
        #return self.getPointsInSphere(bb, pt, radius,addSP=addSP,info=info)
        spacing1 = 1./(self.gridSpacing)
        NX, NY, NZ = self.nbGridPoints
        OX, OY, OZ = self.boundingBox[0] # origin of fill grid-> bottom lef corner not origin
        ox, oy, oz = bb[0]
        ex, ey, ez = bb[1]

        i0 = int(max(0, floor((ox-OX)*spacing1)))
        i1 = int(min(NX, int((ex-OX)*spacing1))+1)

        j0 = int(max(0, floor((oy-OY)*spacing1)))
        j1 = int(min(NY, int((ey-OY)*spacing1))+1)

        k0 = int(max(0, floor((oz-OZ)*spacing1)))
        k1 = int(min(NZ, int((ez-OZ)*spacing1))+1)

        i0 = int(min( NX-1, max( 0, round((ox-OX)*spacing1))))
        j0 = int(min( NY-1, max( 0, round((oy-OY)*spacing1))))
        k0 = int(min( NZ-1, max( 0, round((oz-OZ)*spacing1))))
        i1 = int(min( NX, max( 0, round((ex-OX)*spacing1))))
        j1 = int(min( NY, max( 0, round((ey-OY)*spacing1))))
        k1 = int(min( NZ, max( 0, round((ez-OZ)*spacing1))))

        if NZ == 1 :
            k0=0
            k1=1
        elif NY == 1:
            j0=0
            j1=1
        elif NX == 1:
            i0=0
            i1=1

        ptIndices=[]
#        print "Ob",self.boundingBox[0] 
#        print "bb",bb[0]
#        print "ijk",i0,j0,k0,i,j,k
#        print "coord",self._x[i0],self._y[j0],self._z[k0],self._x[i],self._y[j],self._z[k]
#        print i0,i1,j0,j1,k0,k1
#        x=self._x[i0:i1]
#        y=self._y[j0:j1]
#        z=self._z[k0:k1]
#        positions = numpy.vstack(numpy.meshgrid(x,y,z,copy=False)).reshape(3,-1).T
#        distances=spatial.distance.cdist(self.masterGridPositions,positions)
#        ptIndices = numpy.nonzero(distances == 0)[0]
        pid=numpy.mgrid[i0:i1,j0:j1,k0:k1]
        ijk=numpy.vstack(pid).reshape(3,-1).T
        #in case 2D, meaning one of the dimension is 1
        if NZ == 1 :
            ptIndices = [p[2]+p[1]+NX*p[0] for p in ijk]
        elif NY == 1:
            ptIndices = [p[2]+p[1]+NX*p[0] for p in ijk]
        elif NX == 1:
            ptIndices = [p[2]+NY*p[1]+p[0] for p in ijk]
        else :
            ptIndices = [p[2]+NY*p[1]+NX*NY*p[0] for p in ijk]
#        print "coordi",ptIndices[0],self.masterGridPositions[ptIndices[0]]
       # add surface points
        if addSP and self.nbSurfacePoints != 0:
            result = numpy.zeros( (self.nbSurfacePoints,), 'i')
            nb = self.surfPtsBht.closePoints(tuple(pt), radius, result )
#            nb = self.surfPtsBht.query(tuple(pt),k=self.nbSurfacePoints)
            dimx, dimy, dimz = self.nbGridPoints
            ptIndices.extend(list(map(lambda x, length=self.gridVolume:x+length,
                             result[:nb])) )
        return ptIndices

      
    def getPointsInCube(self, bb, pt, radius,addSP=True,info=False):
        """
        Return all grid points indices inside the given bouding box.
        NOTE : need to fix with grid build with numpy arrange
        """        
        #return self.getPointsInCubeFillBB(bb, pt, radius,addSP=addSP,info=info)
        #return self.getPointsInSphere(bb, pt, radius,addSP=addSP,info=info)
        spacing1 = 1./self.gridSpacing
        
        NX, NY, NZ = self.nbGridPoints
        OX, OY, OZ = self.boundingBox[0] # origin of fill grid-> bottom lef corner not origin
        ox, oy, oz = bb[0]
        ex, ey, ez = bb[1]
#        print("getPointsInCube bb[0] = ",bb[0])
#        print("getPointsInCube bb[1] = ",bb[1])
        #        i0 = max(0, int((ox-OX)*spacing1)+1)
        i0 = int(max(0, floor((ox-OX)*spacing1)))
        i1 = int(min(NX, int((ex-OX)*spacing1)+1))
        #        j0 = max(0, int((oy-OY)*spacing1)+1)
        j0 = int(max(0, floor((oy-OY)*spacing1)))
        j1 = int(min(NY, int((ey-OY)*spacing1)+1))
        #        k0 = max(0, int((oz-OZ)*spacing1)+1)
        k0 = int(max(0, floor((oz-OZ)*spacing1)))
        k1 = int(min(NZ, int((ez-OZ)*spacing1)+1))
#        print("oz-OZ = ", oz-OZ)
#        print("((oz-OZ)*spacing1) = ", ((oz-OZ)*spacing1))
#        print("int((oz-OZ)*spacing1)) = ", int((oz-OZ)*spacing1))
#        print("int((oz-OZ)*spacing1)+1 = ", int((oz-OZ)*spacing1)+1)
#        print("floor(oz-OZ)*spacing1) = ", floor((oz-OZ)*spacing1))
#        print("i0= ", i0, ", i1= ", i1,", j0= ", j0,", j1= ",j1,", k0= ", k0, ", k1= ", k1)

        zPlaneLength = NX*NY

        ptIndices = []
        for z in range(int(k0),int(k1)):
            offz = z*zPlaneLength
            for y in range(int(j0),int(j1)):
                off = y*NX + offz
                #ptIndices.extend(numpy.arange(i0,i1)+off)
                for x in range(int(i0),int(i1)):
                    ptIndices.append( x + off)
#                    print("position of point ",x+off," = ", self.masterGridPositions[x+off])

        # add surface points
        if addSP and self.nbSurfacePoints != 0:
            result = numpy.zeros( (self.nbSurfacePoints,), 'i')
            nb = self.surfPtsBht.closePoints(tuple(pt), radius, result )
#            nb = self.surfPtsBht.query(tuple(pt),k=self.nbSurfacePoints)
            dimx, dimy, dimz = self.nbGridPoints
            ptIndices.extend(list(map(lambda x, length=self.gridVolume:x+length,
                             result[:nb])) )
        return ptIndices

    def computeGridNumberOfPoint(self,boundingBox,space):
        """
        Return the grid size : total number of point and number of point per axes
        """        
        xl,yl,zl = boundingBox[0]
        xr,yr,zr = boundingBox[1]
        encapsulatingGrid = self.encapsulatingGrid  
        #Graham Added on Oct17 to allow for truly 2D grid for test fills... may break everything!        
        
        nx = int(ceil((xr-xl)/space))+encapsulatingGrid
        ny = int(ceil((yr-yl)/space))+encapsulatingGrid
        nz = int(ceil((zr-zl)/space))+encapsulatingGrid
        return nx*ny*nz,(nx,ny,nz)

    def set_surfPtsBht(self,verts):
        from bhtree import bhtreelib
        self.surfPtsBht = None
        if verts :
            self.surfPtsBht = bhtreelib.BHtree( verts, None, 10)

    def set_surfPtscht(self,verts):
        from scipy import spatial
        self.surfPtsBht = None
        if verts :
            self.surfPtsBht = spatial.cKDTree( verts, leafsize=10)

    def computeExteriorVolume(self,compartments=None,space=None,fbox_bb=None):
        # compute exterior volume, totaVolume without compartments volume 
        unitVol = self.gridSpacing**3
        totalVolume = self.gridVolume*unitVol
        if fbox_bb is not None :
            V,nbG = self.computeGridNumberOfPoint(fbox_bb,space)
            totalVolume = V*unitVol
        if compartments is not None :
            for o in compartments:
                #totalVolume -= o.surfaceVolume
                totalVolume -= o.interiorVolume
        return totalVolume

    def computeVolume(self,space=None,fbox_bb=None):
        # compute exterior volume, totaVolume without compartments volume 
        unitVol = self.gridSpacing**3
        totalVolume = self.gridVolume*unitVol
        if fbox_bb is not None :
            V,nbG = self.computeGridNumberOfPoint(fbox_bb,space)
            totalVolume = V*unitVol
        return totalVolume
        
#==============================================================================
# TO DO File IO 
#==============================================================================
    def save(self):
        pass
    def restore(self):
        pass

#dont forget to use spatial.distance.cdist
class HaltonGrid(Grid):
    def __init__(self,boundingBox=([0,0,0], [.1,.1,.1]),space=1):
        Grid.__init__(self, boundingBox=boundingBox,space=space)
        self.haltonseq = cHaltonSequence3()
        self.tree = None
        
    def getPointFrom3D(self, pt3d):
        """
        get point number from 3d coordinates
        """
        #actually look at closest point using ckdtree
        #nb = self.tree.query_ball_point(point,cutoff)
        distance,nb = self.tree.query(pt3d,1)#len of ingr posed so far
        return nb

    def getScale(self, boundingBox=None):
        if boundingBox is None :
            boundingBox = self.boundingBox
        xl,yl,zl = boundingBox[0]#0
        xr,yr,zr = boundingBox[1]#1
        encapsulatingGrid = self.encapsulatingGrid  #Graham Added on Oct17 to allow for truly 2D grid for test fills... may break everything!
        txyz = numpy.array(boundingBox[0])
        scalexyz = numpy.array(boundingBox[1]) - txyz        
        return scalexyz,txyz
        
    def create3DPointLookup(self, boundingBox=None):
        #we overwrite here by using the halton sequence dimension 5 to get
        #the coordinate
        self.haltonseq.reset()
        pointArrayRaw = numpy.zeros( (nx*ny*nz, 3), 'f')
        self.ijkPtIndice = numpy.zeros( (nx*ny*nz, 3), 'i')
        nx,ny,nz = self.nbGridPoints
        scalexyz,txyz = self.getScale()
        i = 0
        for zi in range(nz):
            for yi in range(ny):
                for xi in range(nx):
                    self.haltonseq.inc()
                    pointArrayRaw[i] = numpy.array([self.haltonseq.mX, 
                                                    self.haltonseq.mY, 
                                                    self.haltonseq.mZ])
                    self.ijkPtIndice[i] = (xi,yi,zi)
                    i+=1
        #scale the value from the halton(0...1) to grid boundin box
        self.masterGridPositions = pointArrayRaw*scalexyz+txyz        
        self.tree = spatial.cKDTree(self.masterGridPositions, leafsize=10)
        
