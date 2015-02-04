# -*- coding: utf-8 -*-
"""
# autoPACK Authors: Graham T. Johnson, Mostafa Al-Alusi, Ludovic Autin, Michel Sanner
#   Based on COFFEE Script developed by Graham Johnson between 2005 and 2010 
#   with assistance from Mostafa Al-Alusi in 2009 and periodic input 
#   from Arthur Olson's Molecular Graphics Lab
#
# Compartment.py Authors: Graham Johnson & Michel Sanner with editing/enhancement from Ludovic Autin
#
# Translation to Python initiated March 1, 2010 by Michel Sanner with Graham Johnson
#
# Class restructuring and organization: Michel Sanner
#
# Copyright: Graham Johnson ©2010
#
# This file "Compartment.py" is part of autoPACK, cellPACK.
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
###############################################################################
@author: Graham Johnson, Ludovic Autin, & Michel Sanner

# Hybrid version merged from Graham's Sept 6, 2011 and Ludo's April 2012
#version on May 16, 2012, remerged on July 5, 2012 with thesis versions
"""

# Hybrid version merged from Graham's Sept 2011 and Ludo's April 2012 version on May 16, 2012
# Updated with Sept 16, 2011 thesis versions on July 5, 2012

# TODO: Describe Organelle class here at high level

# TODO: Graham and Ludovic implemented a 2D density function to obtain target numbers for
#   filling surfaces.  This should be formalized and named something other than molarity
#   or molarity should be converted to a 2D value behind the scenes.
# IDEA: We should offer the user an option to override molarity with a specific
#   number, e.g., "I want to place 3 1xyz.pdb files in compartment A" rather than
#   forcing them to calculate- "I need to place 0.00071M of 1xyz.pdb to get 3 of them
#   in an compartment A of volume=V."

## IDEAS

## randomly select recipe and then randomly select free point in set of free
## points corresponding to this recipe would allow giving surface more
## chances to get filled

## NOTE changing smallest molecule radius changes grid spacing and invalidates
##      arrays saved to file
import sys
import numpy.oldnumeric as N
import numpy
import pickle
import weakref
import pdb
from time import time,sleep
from .ray import vlen, vdiff, vcross
import math
from math import sqrt,ceil
import os
try :
    import urllib.request as urllib# , urllib.parse, urllib.error
except :
    import urllib

#from .Ingredient import Ingredient
from .Recipe import Recipe
#print "import AutoFill"
import autopack
from autopack import checkURL
from RAPID import RAPIDlib

AFDIR = autopack.__path__[0]

try :
    helper = autopack.helper
except :
    helper = None
print ("helper is "+str(helper))

try :
    import panda3d
    print ("Should have Panda3D now because panda3d = ", panda3d)
    
    from panda3d.core import Mat3,Mat4,Vec3,Point3
    from panda3d.core import TransformState
    from panda3d.core import BitMask32
    from panda3d.bullet import BulletSphereShape,BulletBoxShape,BulletCylinderShape
    #        from panda3d.bullet import BulletUpAxis
    from panda3d.bullet import BulletRigidBodyNode
    from panda3d.ode import OdeBody, OdeMass
    from panda3d.ode import OdeSphereGeom
    from panda3d.core import NodePath
    print ("Got Panda3D Except")
except :
    panda3d = None
    print ("Failed to get Panda, because panda3d = ", panda3d)

#from autopack import intersect_RayTriangle as iRT
#from autopack.Environment import Grid   
if sys.version > "3.0.0":
    xrange = range
    
class CompartmentList:
    """
    The CompartmentList class
    ==========================
    Handle a list of compartments.
    """
    
    def __init__(self):
        # list of compartments inside this compartment
        self.compartments = []

        # point to parent compartment or Environment
        self.parent = None

    def addCompartment(self, compartment):
        """add a new compartment to the list"""
        assert compartment.parent == None
        assert isinstance(compartment, Compartment)
        self.compartments.append(compartment)
        compartment.parent = self

class  Compartment(CompartmentList):
    """
    The Compartment class
    ==========================
    This class represents a sub volume delimited by a polyhedral
    surface. Compartment can be nested
    """

    def __init__(self, name, vertices, faces, vnormals, **kw):
        CompartmentList.__init__(self)
        #print ("compartment init",name,kw)
        self.name = name
        self.center= None
        self.vertices = vertices
        self.faces = faces
        self.vnormals = vnormals
        self.fnormals = None
        self.area = 0.0
        if "fnormals" in kw :
            self.fnormals = kw["fnormals"]
        self.mesh = None
        self.rbnode = None
        self.gname= ""
        self.ghost =False
        self.bb = None
        self.diag = 9999.9
        if "ghost" in kw :
            self.ghost = kw["ghost"]
        self.ref_obj = None
        if "ref_obj" in kw :
            self.ref_obj = kw["ref_obj"]
        if vertices == None :
            if "filename" in kw :
                self.faces,self.vertices,self.vnormals = self.getMesh(filename=kw["filename"])
                print (self.vertices[0])
                self.filename=kw["filename"]
                self.ref_obj = name
                #print ("mesh",self.name,self.filename)
        self.encapsulatingRadius = 9999.9  
        if self.vertices is not None and len(self.vertices):
            #can be dae/fbx file, object name that have to be in the scene or dejaVu indexedpolygon file
            self.bb = self.getBoundingBox()
            v=numpy.array(self.vertices,'f')
            l=numpy.sqrt((v*v).sum(axis=1))
            self.encapsulatingRadius = max(l)   

        self.checkinside = True
        self.representation = None
        self.representation_file = None
        if "object_name" in kw :
            if kw["object_name"] is not None:
                print ("rep",kw["object_name"],kw["object_filename"])
                self.representation = kw["object_name"]
                self.representation_file =  kw["object_filename"]
                self.getMesh(filename=self.representation_file,rep=self.representation)
        self.innerRecipe = None
        self.surfaceRecipe = None
        self.surfaceVolume = 0.0
        self.interiorVolume = 0.0
        self.rapid_model = None
        # list of grid point indices inside organelle
        self.insidePoints = None
        # list of grid point indices on compartment surface
        self.surfacePoints = None
        self.surfacePointsNormals = {} # will be point index:normal
        
        self.number = None # will be set to an integer when this compartment
                           # is added to a Environment. Positivefor surface pts
                           # negative for interior points
#        self.parent = None
        self.molecules = [] 
        # list of ( (x,y,z), rotation, ingredient) triplet generated by fill 
        self.overwriteSurfacePts = True 
        # do we discretize surface point per edges
        self.highresVertices = None 
        # if a highres vertices is provided this give the surface point, 
        #not the one provides
        self.isBox = False #the compartment shape is a box, no need 
        #to compute inside points.
        if "isBox" in kw :
            self.isBox = kw["isBox"]
        self.isOrthogonalBoudingBox = None #the compartment shape is a box, no need
        #to compute inside points.
        if "isOrthogonalBoudingBox" in kw :
            self.isOrthogonalBoudingBox = kw["isOrthogonalBoudingBox"]

    def reset(self):
        """reset the inner compartment data, surface and inner points"""
        # list of grid point indices inside compartment
        self.insidePoints = None
        # list of grid point indices on compartment surface
        self.surfacePoints = None
        self.surfacePointsNormals = {} # will be point index:normal
        #self.molecules = [] 
        

    def getDejaVuMesh(self, filename, geomname):
        """
        Create a DejaVu polygon mesh object from a filename 
    
        @type  filename: string
        @param filename: the name of the input file
        @type  geomname: string
        @param geomname: the name of the ouput geometry
    
        @rtype:   DejaVu.IndexedPolygons
        @return:  the created dejavu mesh  
        """
        #filename or URL
        from DejaVu.IndexedPolygons import IndexedPolygonsFromFile
        #seems not ok...when they came from c4d ... some transformation are not occuring.
        print ("dejavu mesh", filename)        
        geom = IndexedPolygonsFromFile(filename, 'mesh_%s'%self.name)
#        if helper is not None:
#            helper.rotateObj(geom,[0.0,-math.pi/2.0,0.0])
        return geom

    def rapid_model(self):
        rapid_model = RAPIDlib.RAPID_model()
        rapid_model.addTriangles(numpy.array(self.vertices,'f'), numpy.array(self.faces,'i'))
        return rapid_model       
 
    def create_rapid_model(self):
        self.rapid_model = RAPIDlib.RAPID_model()
        #need triangle and vertices
        #faces,vertices,vnormals = helper.DecomposeMesh(self.mesh,
        #                   edit=False,copy=False,tri=True,transform=True) 
        self.rapid_model.addTriangles(numpy.array(self.vertices,'f'), numpy.array(self.faces,'i'))
                   
    def get_rapid_model(self):
        if (self.rapid_model is None ):
            self.create_rapid_model()
        return self.rapid_model 

    def addMeshRB(self,):
        from panda3d.bullet import BulletRigidBodyNode
        inodenp = self.parent.worldNP.attachNewNode(BulletRigidBodyNode(self.name))
        inodenp.node().setMass(1.0)
        helper = autopack.helper
        from panda3d.core import GeomVertexFormat,GeomVertexWriter,GeomVertexData,Geom,GeomTriangles
        from panda3d.core import GeomVertexReader
        from panda3d.bullet import BulletTriangleMesh,BulletTriangleMeshShape,BulletConvexHullShape
        #step 1) create GeomVertexData and add vertex information
        format=GeomVertexFormat.getV3()
        vdata=GeomVertexData("vertices", format, Geom.UHStatic)        
        vertexWriter=GeomVertexWriter(vdata, "vertex")
        [vertexWriter.addData3f(v[0],v[1],v[2]) for v in self.vertices]

        #step 2) make primitives and assign vertices to them
        tris=GeomTriangles(Geom.UHStatic)
        [self.setGeomFaces(tris,face) for face in self.faces]

        #step 3) make a Geom object to hold the primitives
        geom=Geom(vdata)
        geom.addPrimitive(tris)
        #step 4) create the bullet mesh and node
#        if ingr.convex_hull:
#            shape = BulletConvexHullShape()
#            shape.add_geom(geom)
#        else :
        mesh = BulletTriangleMesh()
        mesh.addGeom(geom)    
        shape = BulletTriangleMeshShape(mesh, dynamic=False)#BulletConvexHullShape            
        print ("shape ok",shape)
        #inodenp = self.worldNP.attachNewNode(BulletRigidBodyNode(ingr.name))
        #inodenp.node().setMass(1.0)
        inodenp.node().addShape(shape)#,TransformState.makePos(Point3(0, 0, 0)))#, pMat)#TransformState.makePos(Point3(jtrans[0],jtrans[1],jtrans[2])))#rotation ?
        return inodenp

    def create_rbnode(self, ):
        # Sphere
        if panda3d is None :
            return None
        trans=[0,0,0]
        rotMat = mat = numpy.identity(4)
        mat = mat.transpose().reshape((16,))
        mat3x3 =Mat3(mat[0],mat[1],mat[2],
                     mat[4],mat[5],mat[6],
                     mat[8],mat[9],mat[10])
        pmat   = Mat4(mat[0],mat[1],mat[2],mat[3],
                                           mat[4],mat[5],mat[6],mat[7],
                                           mat[8],mat[9],mat[10],mat[11],
                                           trans[0],trans[1],trans[2],mat[15],)
        pMat = TransformState.makeMat(pmat)
        if self.parent.panda_solver == "ode":
            pMat = mat3x3
        shape = None
        inodenp = self.addMeshRB(pMat,trans,rotMat)
        if self.panda_solver == "bullet":
            inodenp.setCollideMask(BitMask32.allOn())
            inodenp.node().setAngularDamping(1.0)
            inodenp.node().setLinearDamping(1.0)
            inodenp.setMat(pmat)
            self.parent.world.attachRigidBody(inodenp.node())
            inodenp = inodenp.node()
        elif self.panda_solver == "ode":
            inodenp.setCollideBits(BitMask32(0x00000002))
            inodenp.setCategoryBits(BitMask32(0x00000001))
            #boxGeom.setBody(boxBody)
        self.parent.rb_panda.append(inodenp)
        #self.moveRBnode(inodenp.node(), trans, rotMat)
        return inodenp
        
    def get_rb_model(self):
        if self.rbnode is None :
            self.rbnode=self.create_rbnode()
        return self.rbnode
       
    def getMesh(self,filename=None,rep=None,**kw):
        """
        Retrieve the compartment 3d representaton from the given filename
        
        @type  filename: string
        @param filename: the name of the input file
        @type  rep: string
        @param rep: the name of the input file for the representation
        """
        geom = None
        gname = self.name
        helper = autopack.helper
        if helper is not None :
            parent=helper.getObject("autopackHider")
            if parent is None:
                parent = helper.newEmpty("autopackHider")
        if rep is not None :
            gname =rep 
            if helper is not None :
                parent=helper.getObject('O%s'%self.name)
        print ("compartments ",self.name,filename,gname,rep)
        #identify extension
        name = filename.split("/")[-1]
        fileName, fileExtension = os.path.splitext(name)
        if fileExtension is '' :  
            tmpFileName1 =autopack.retrieveFile(filename+".indpolface",cache="geometries")
            tmpFileName2 =autopack.retrieveFile(filename+".indpolvert",cache="geometries")
            filename = os.path.splitext(tmpFileName1)[0]
        else :
            filename =autopack.retrieveFile(filename,cache="geometries")  
        if filename is None :
            print ("problem with getMesh of compartm",self.name,fileName, fileExtension,rep)
            return            
        if not os.path.isfile(filename) and fileExtension != '' :
            print ("problem with "+filename)
            return
        print ("filename compartements is",filename,type(filename))     
        fileName, fileExtension = os.path.splitext(filename)
        print('found fileName '+fileName+' fileExtension '+fileExtension)
        if fileExtension.lower() == ".fbx":
            print ("read withHelper",filename)
            #use the host helper if any to read
            if helper is not None:#neeed the helper
                helper.read(filename)
                geom = helper.getObject(gname)
                #reparent to the fill parent
                if helper.host == "3dsmax" or helper.host.find("blender") != -1:
                    helper.resetTransformation(geom)#remove rotation and scale from importing
                if helper.host != "c4d" and rep == None  and helper.host != "softimage":
                    #need to rotate the transform that carry the shape
                    helper.rotateObj(geom,[0.,-math.pi/2.0,0.0])
                if helper.host =="softimage"  :
                    helper.rotateObj(geom,[0.0,-math.pi/2.0,0.0],primitive=True)#need to rotate the primitive                    
#                p=helper.getObject("autopackHider")
#                if p is None:
#                    p = helper.newEmpty("autopackHider")
#                    helper.toggleDisplay(p,False)
                helper.reParent(geom,parent)
#                return geom
#            return None
        elif fileExtension == ".dae":
            #use the host helper if any to read
            if helper is None :
                from upy.dejavuTk.dejavuHelper import dejavuHelper
                #need to get the mesh directly. Only possible if dae or dejavu format
                #get the dejavu heper but without the View, and in nogui mode
                h =  dejavuHelper(vi="nogui")
                dgeoms = h.read(filename)
                v,vn,f = dgeoms.values()[0]["mesh"]
                vn = self.getVertexNormals(v,f)                
                return f,v,vn
            else : #if helper is not None:#neeed the helper
                if helper.host == "dejavu" and helper.nogui:
                    dgeoms = helper.read(filename)
                    v,vn,f = dgeoms.values()[0]["mesh"]
#                    print v[0],vn[0]
#                    vn = self.getVertexNormals(numpy.array(v),f[:])     
                    self.mesh = helper.createsNmesh(gname,v,None,f)[0]
#                    print v[0],vn[0]
                    return f,v,vn
                helper.read(filename)
                geom = helper.getObject(gname)
                print ("should have read...",gname,geom,parent)
#                helper.update()
                if helper.host == "3dsmax" or helper.host.find("blender") != -1:
                    helper.resetTransformation(geom)#remove rotation and scale from importing
#                if helper.host != "c4d"  and helper.host != "dejavu"  and helper.host != "softimage":
                    #need to rotate the transform that carry the shape
#                    helper.rotateObj(geom,[0.0,-math.pi/2.0,0.0])#wayfront as well euler angle
                if helper.host =="softimage"  :
                    helper.rotateObj(geom,[0.0,-math.pi/2.0,0.0],primitive=True)#need to rotate the primitive                    
#                if helper.host =="c4d"  :
                    #remove the normaltag if exist
#                    helper.removeNormalTag(geom)
#                p=helper.getObject("AutoFillHider")
#                if p is None:
#                    p = helper.newEmpty("autopackHider")
#                    helper.toggleDisplay(p,False)
                print ("reparent ",geom,parent)
                helper.reParent(geom,parent)            
#                return geom
#            return None
        elif fileExtension is '' :
            geom =  self.getDejaVuMesh(filename, gname)
        else :#speficif host file
            if helper is not None:#neeed the helper
                helper.read(filename)
                geom = helper.getObject(gname)
#                p=helper.getObject("autopackHider")
#                if p is None:
#                    p = helper.newEmpty("autopackHider")
#                    helper.toggleDisplay(p,False)
                helper.reParent(geom,parent) 
        if rep is None:
            if helper.host.find("blender") == -1 :
                helper.toggleDisplay(parent,False) 
            elif helper.host == "dejavu":
                helper.toggleDisplay(geom,False)                 
        else :
            return None 
        self.gname = gname                    
        if geom is not None and fileExtension != '' and not self.ghost:
            faces,vertices,vnormals = helper.DecomposeMesh(geom,
                                edit=False,copy=False,tri=True,transform=True)
            return faces,vertices,vnormals
        else :
            if self.ghost:
                return [],[],[]
            faces = geom.getFaces()
            vertices = geom.getVertices()
            vnormals = geom.getVNormals()       
            return faces,vertices,vnormals
        
    def setMesh(self, filename=None,vertices=None, 
                faces=None, vnormals=None, **kw )   :
        """
        Set the 3d mesh from the given filename or the given mesh data (v,f,n)
        
        @type  filename: string
        @param filename: the name of the input file
        @type  vertices: array
        @param vertices: mesh vertices or None
        @type  faces: array
        @param faces: mesh faces or None
        @type  vnormals: array
        @param vnormals: mesh vnormals or None
        """                    
        if vertices is None and filename is not None:
            self.faces,self.vertices,self.vnormals = self.getMesh(filename)
        else :
            self.vertices = vertices
            self.faces = faces
            self.vnormals = vnormals
        if "fnormals" in kw :
            self.fnormals = kw["fnormals"]
        self.mesh = None
        self.ref_obj = filename 
        self.bb = self.getBoundingBox()
        v=numpy.array(self.vertices,'f')
        l=numpy.sqrt((v*v).sum(axis=1))
        self.encapsulatingRadius = max(l)   

    def saveGridToFile(self, f):
        """Save insidePoints and surfacePoints to file"""
        #print 'surface', len(self.surfacePoints)
        pickle.dump(self.insidePoints, f)
        #print 'interior', len(self.surfacePoints)
        pickle.dump(self.surfacePoints, f)
        pickle.dump(self.surfacePointsNormals, f)
        pickle.dump(self.surfacePointsCoords, f)
        
    def readGridFromFile(self, f):
        """read insidePoints and surfacePoints from file"""
        self.insidePoints = insidePoints = pickle.load(f)
        self.surfacePoints = surfacePoints = pickle.load(f)
        self.surfacePointsNormals = surfacePointsNormals = pickle.load(f)
        self.surfacePointsCoords = surfacePointsCoords = pickle.load(f)
        return surfacePoints, insidePoints, surfacePointsNormals,surfacePointsCoords
        #surfacePointscoords = pickle.load(f)
        #return surfacePoints, insidePoints, surfacePointsNormals, \
        #       surfacePointscoords
        
    def setNumber(self, num):
        """set compartment uniq id"""
        self.number = num

    def setInnerRecipe(self, recipe):
        """set the inner recipe that define the ingredient to pack inside"""
        assert self.number is not None
        assert isinstance(recipe, Recipe)
        self.innerRecipe = recipe
        self.innerRecipe.number= self.number
        recipe.compartment = self#weakref.ref(self)
        for ingr in recipe.ingredients:
            ingr.compNum = -self.number
            if hasattr(ingr,"compMask"):
                if not ingr.compMask :
                    ingr.compMask=[ingr.compNum]
                    
    def setSurfaceRecipe(self, recipe):
        """set the inner recipe that define the ingredient to pack at the surface"""
        assert self.number is not None
        assert isinstance(recipe, Recipe)
        self.surfaceRecipe = recipe
        self.surfaceRecipe.number= self.number
        recipe.compartment = self#weakref.ref(self)
        for ingr in recipe.ingredients:
            ingr.compNum = self.number
            

    def getCenter(self):
        """get the center of the mesh (vertices barycenter)"""
        if self.center is None :
            coords = numpy.array(self.vertices)#self.allAtoms.coords
            center = sum(coords)/(len(coords)*1.0)
            center = list(center)
            for i in range(3):
                center[i] = round(center[i], 4)
    #        print "center =", center
            self.center = center
    
    def getRadius(self):
        """get the radius as the distance between vertices center and bottom left bouding box"""
        import math
        d = self.center - self.bb[0]
        s = numpy.sum(d*d)
        return math.sqrt(s)
        
    def getBoundingBox(self):
        """get the bounding box"""
        from autopack.Environment import vlen,vdiff
        mini = numpy.min(self.vertices, 0)
        maxi = numpy.max(self.vertices, 0)
        xl,yl,zl = mini
        xr,yr,zr = maxi
        self.diag = vlen( vdiff((xr,yr,zr), (xl,yl,zl) ) )    
        return (mini, maxi)

    def getSizeXYZ(self):
        """get the size per axe"""
        sizexyz=[0,0,0]
        for i in range(3):
            sizexyz[i] = self.bb[1][i] - self.bb[0][i]
        return sizexyz
        
#    def checkPointInsideBBold(self,pt3d,dist=None):
#        """check if the given 3d coordinate is inside the compartment bouding box"""        
#        O = numpy.array(self.bb[0])
#        E = numpy.array(self.bb[1])
#        P = numpy.array(pt3d)
#        test1 = P < O 
#        test2 =  P > E
#        if True in test1 or True in test2:
#            #outside
#            return False
#        else :
#            if dist is not None:
#                d1 = P - O
#                s1 = numpy.sum(d1*d1)
#                d2 = E - P
#                s2 = numpy.sum(d2*d2)
#                if s1 <= dist or s2 <=dist:
#                    return False 
#            return True
 
    def checkPointInsideBB(self,pt3d,dist=None):
        """check if the given 3d coordinate is inside the compartment bouding box"""
        O = numpy.array(self.bb[0])
        E = numpy.array(self.bb[1])
        P = numpy.array(pt3d)
        
        #a point is inside is  < min and > maxi etc.. 
        test1 = P < O 
        test2 =  P > E
        if True in test1 or True in test2:
            #outside
            return False
        else :
            if dist is not None:
                d1 = P - O
                s1 = numpy.sum(d1*d1)
                d2 = E - P
                s2 = numpy.sum(d2*d2)
                if s1 <= dist or s2 <=dist:
                    return False 
            return True
 
    def inBox(self, box):
        """
        check if bounding box of this compartment fits inside the give box
        returns true or false and the extended bounding box if this compartment
        did not fit
        """
        if self.ghost:
            return False, None
        bb = self.bb
        xm,ym,zm = box[0]
        xM,yM,zM = box[1]
        # padding 50 shows problem
        padding = 0.
        
        newBB = [box[0][:], box[1][:]]
        fits = True

        if xm > bb[0][0] - padding:
            newBB[0][0] = bb[0][0] - padding
            fits = False

        if ym > bb[0][1] - padding:
            newBB[0][1] = bb[0][1] - padding
            fits = False

        if zm > bb[0][2] - padding:
            newBB[0][2] = bb[0][2] - padding
            fits = False

        if xM < bb[1][0] + padding:
            newBB[1][0] = bb[1][0] + padding
            fits = False

        if yM < bb[1][1] + padding:
            newBB[1][1] = bb[1][1] + padding
            fits = False

        if zM < bb[1][2] + padding:
            newBB[1][2] = bb[1][2] + padding
            fits = False

        return fits, newBB

    def inGrid(self, point, fillBB):
        """
        check if bounding box of this compartment fits inside the give box
        returns true or false and the extended bounding box if this compartment
        did not fit
        """
        mini, maxi = fillBB
        mx, my, mz = mini
        Mx, My, Mz = maxi
        x,y,z = point
        if (x>=mx and x<=Mx and y>=my and y<=My and z>=mz and z<=Mz):
            return True
        else :
            return False

    def getMinMaxProteinSize(self):
        """retrieve minimum and maximum ingredient size for inner and surface recipe ingredients"""
        #for compartment in self.compartments:
        #    mini, maxi = compartment.getSmallestProteinSize(size)
        mini1 = mini2 = 9999999.
        maxi1 = maxi2 = 0.
        if self.surfaceRecipe:
            mini1, maxi1 = self.surfaceRecipe.getMinMaxProteinSize()
        if self.innerRecipe:
            mini2, maxi2 = self.innerRecipe.getMinMaxProteinSize()
        return min(mini1, mini2), max(maxi1, maxi2)

    def getVertexNormals(self, vertices, faces):
        vnormals = vertices[:]
        face = [[0,0,0],[0,0,0],[0,0,0]]
        v = [[0,0,0],[0,0,0],[0,0,0]]        
        for f in faces:
            for i in range(3) :
                face [i] = vertices[f[i]]
            for i in range(3) :
                v[0][i] = face[1][i]-face[0][i]
                v[1][i] = face[2][i]-face[0][i]                
            normal = vcross(v[0],v[1])
            n = vlen(normal)
            if n == 0. :
                n1=1.
            else :
                n1 = 1./n
            for i in range(3) :
                vnormals[f[i]] = [normal[0]*n1, normal[1]*n1, normal[2]*n1]
        return vnormals #areas added by Graham

    def getFaceNormals(self, vertices, faces,fillBB=None):
        """compute the face normal of the compartment mesh"""
        normals = []
        areas = [] #added by Graham
        face = [[0,0,0],[0,0,0],[0,0,0]]
        v = [[0,0,0],[0,0,0],[0,0,0]]        
        for f in faces:
            for i in range(3) :
                face [i] = vertices[f[i]]
            for i in range(3) :
                v[0][i] = face[1][i]-face[0][i]
                v[1][i] = face[2][i]-face[0][i]                
            normal = vcross(v[0],v[1])
            n = vlen(normal)
            if n == 0. :
                n1=1.
            else :
                n1 = 1./n
            normals.append( (normal[0]*n1, normal[1]*n1, normal[2]*n1) )
            if fillBB is not None:
                if self.inGrid(vertices[f[0]],fillBB) and self.inGrid(vertices[f[0]],fillBB) and self.inGrid(vertices[f[0]],fillBB):
                    areas.append(0.5*vlen(normal)) #added by Graham
        self.area = sum(areas)
        return normals, areas #areas added by Graham

    def getInterpolatedNormal(self, pt, tri):
        """compute an interpolated normal for te given triangle at the given point"""
        v1,v2,v3 = self.faces[tri]
        verts = self.vertices
        d1 = vlen(vdiff(pt, verts[v1]))
        d2 = vlen(vdiff(pt, verts[v2]))
        d3 = vlen(vdiff(pt, verts[v3]))
        sumlen1 = d1+d2+d3
        w1 = sumlen1/d1
        w2 = sumlen1/d2
        w3 = sumlen1/d3
        n1 = self.vnormals[v1]
        n2 = self.vnormals[v2]
        n3 = self.vnormals[v3]
        norm = ( (n1[0]*w1 + n2[0]*w2 +n3[0]*w3),
                 (n1[1]*w1 + n2[1]*w2 +n3[1]*w3),
                 (n1[2]*w1 + n2[2]*w2 +n3[2]*w3) )
        l1 = 1./vlen(norm)
        return ( norm[0]*l1, norm[1]*l1, norm[2]*l1 )

    def createSurfacePoints(self, maxl=20):
        """
        create points inside edges and faces with max distance between then maxl
        creates self.surfacePoints and self.surfacePointsNormals
        """
        vertices = self.vertices
        faces = self.faces
        vnormals = self.vnormals
        
        points = list(vertices)[:]
        normals = list(vnormals)[:]

        # create points in edges
        edges = {}
        for fn, tri in enumerate(faces):
            s1,s2 = tri[0],tri[1]
            if (s2, s1) in edges:
                edges[(s2,s1)].append(fn)
            else:
                edges[(s1,s2)] = [fn]

            s1,s2 = tri[1],tri[2]
            if (s2, s1) in edges:
                edges[(s2,s1)].append(fn)
            else:
                edges[(s1,s2)] = [fn]

            s1,s2 = tri[2],tri[0]
            if (s2, s1) in edges:
                edges[(s2,s1)].append(fn)
            else:
                edges[(s1,s2)] = [fn]

        lengths = list(map(len, list(edges.values())))
        assert max(lengths)==2
        assert min(lengths)==2

        for edge, faceInd in list(edges.items()):
            s1, s2 = edge
            p1 = vertices[s1]
            p2 = vertices[s2]
            v1 = vdiff(p2, p1) # p1->p2
            l1 = vlen(v1)
            if l1 <= maxl: continue

            # compute number of points
            nbp1 = int(l1 / maxl)
            if nbp1<1: continue

            # compute interval size to spread the points
            dl1 = l1/(nbp1+1)

            # compute interval vector
            dx1 = dl1*v1[0]/l1
            dy1 = dl1*v1[1]/l1
            dz1 = dl1*v1[2]/l1
            x, y, z = p1
            nx1, ny1, nz1 = vnormals[s1]
            nx2, ny2, nz2 = vnormals[s2]
            edgeNorm = ( (nx1+nx2)*.5,  (ny1+ny2)*.5,  (nz1+nz2)*.5 )
            for i in range(1, nbp1+1):
                points.append( (x + i*dx1, y + i*dy1, z + i*dz1) )
                normals.append( edgeNorm )

        for fn,t in enumerate(faces):
            #if t[0]==16 and t[1]==6 and t[2]==11:
            #    pdb.set_trace()
            pa = vertices[t[0]]
            pb = vertices[t[1]]
            pc = vertices[t[2]]

            va = vdiff(pb, pa) # p1->p2
            la = vlen(va)
            if la <= maxl: continue

            vb = vdiff(pc, pb) # p2->p3
            lb = vlen(vb)
            if lb <= maxl: continue

            vc = vdiff(pa, pc) # p3->p1
            lc = vlen(vc)
            if lc <= maxl: continue

            #if fn==0:
            #    pdb.set_trace()
            # pick shortest edge to be second vector
            if la<=lb and la<=lc:
                p1 = pc
                p2 = pa
                p3 = pb
                v1 = vc
                l1 = lc
                v2 = va
                l2 = la
                v3 = vb

            if lb<=la and lb<=lc:
                p1 = pa
                p2 = pb
                p3 = pc
                v1 = va
                l1 = la
                v2 = vb
                l2 = lb
                v3 = vc

            if lc<=lb and lc<=la:
                p1 = pb
                p2 = pc
                p3 = pa
                v1 = vb
                l1 = lb
                v2 = vc
                l2 = lc
                v3 = va

            lengthRatio = l2/l1

            nbp1 = int(l1 / maxl)
            if nbp1<1: continue

            dl1 = l1/(nbp1+1)
            #if dl1<15:
            #    pdb.set_trace()
            #print l1, nbp1, dl1, lengthRatio
            dx1 = dl1*v1[0]/l1
            dy1 = dl1*v1[1]/l1
            dz1 = dl1*v1[2]/l1
            x,y,z = p1
            fn = vcross(v1, (-v3[0], -v3[1], -v3[2]) )
            fnl = 1.0/vlen(fn)
            faceNorm = ( (fn[0]*fnl, fn[1]*fnl, fn[2]*fnl) )

            for i in range(1, nbp1+1):
                l2c = (i*dl1)*lengthRatio
                nbp2 = int(l2c/maxl)
#                percentage = (i*dl1)/l1
                #nbp2 = int(l2*lengthRatio*percentage/maxl)
                if nbp2<1: continue
                #dl2 = l2*percentage/(nbp2+1)
                dl2 = l2c/(nbp2+1)
                #print '   ',i, percentage, dl1, l2c, dl2, nbp2, l2
                #if dl2<15:
                #    pdb.set_trace()

                dx2 = dl2*v2[0]/l2
                dy2 = dl2*v2[1]/l2
                dz2 = dl2*v2[2]/l2
                for j in range(1, nbp2+1):
                    points.append( (
                        x + i*dx1 + j*dx2, y + i*dy1 + j*dy2, z + i*dz1 + j*dz2) )
                    normals.append(faceNorm)

        self.ogsurfacePoints = points
        self.ogsurfacePointsNormals = normals

    def one_rapid_ray(self,pt1,pt2,diag):
        #return number of triangle /triangle contact
        helper = autopack.helper
        rm = self.get_rapid_model()
        v1=numpy.array(pt1)
        direction = helper.unit_vector(pt2-pt1)*diag
        v2=v1+direction
        if sum(v1) == 0.0 :
            v3 = v2+numpy.array([0.,1.,0.])
        else :
            v3=v2+helper.unit_vector(numpy.cross(v1,v2))
        f=[0,1,2]
        ray_model = RAPIDlib.RAPID_model()
        ray_model.addTriangles(numpy.array([v1,v2,v3],'f'), numpy.array([f,],'i'))
        RAPIDlib.RAPID_Collide_scaled(numpy.identity(3), numpy.array([0.0,0.0,0.0],'f'), 
                                      1.0,rm, numpy.identity(3), 
                                    numpy.array([0.0,0.0,0.0],'f'), 1.0,ray_model,
                                      RAPIDlib.cvar.RAPID_ALL_CONTACTS);
        #could display it ?
        #print numpy.array([v1,v2,v3],'f')
        return RAPIDlib.cvar.RAPID_num_contacts
        

    def checkPointInside_rapid(self,point,diag,ray=1):
        #we want to be sure to cover the organelle
        if diag < self.diag :
            diag = self.diag
        inside = False
        v1=numpy.array(point)
        self.getCenter()
        count1 = self.one_rapid_ray(v1,v1+numpy.array([0.,0.0,1.1]),diag)
        r= ((count1 % 2) == 1) #inside ?
        #we need 2 out of 3 ?
        if ray == 3 :    
            count2 = self.one_rapid_ray(v1, numpy.array(self.center) ,diag )
            if (count2 % 2) == 1 and r :
                return True
            count3 = self.one_rapid_ray(v1, v1+numpy.array([0.0,1.1,0.]),diag )
            if (count3 % 2) == 1 and r :
                return True
            return False
#            c = count1+count2+count3
#            if c == 3 : #1 1 1 
#                r = True
#            elif c == 4 : # 1 2 1
#                r = True
#            elif c == 5 : #2 2 1 or 3 1 1 or 3 2 0
#                r = False
#            elif c == 6 : # 2 2 2 or 3 3 0 or ...
#                r = False
#            if r < 5 :
#                r = True
#            else :
#                r = False
#            if r :
#               if (count2 % 2) == 1 and (count3 % 2) == 1 :
#                   r=True
#               else : 
#                   r=False
        if r : # odd inside
            inside = True
        return inside
        
    def checkPointInside(self,point,diag,ray=1):
        inside = False
        insideBB  = self.checkPointInsideBB(point)#cutoff?
        r=False
        if insideBB:
            helper = autopack.helper
            if helper.host == "dejavu":
                return insideBB
            geom =   helper.getObject(self.gname)      
            if geom is None :
                self.gname = '%s_Mesh'%self.name            
                geom = helper.getObject(self.gname)  
            center = helper.getTranslation( geom )
            intersect, count = helper.raycast(geom, point, center, diag, count = True )
            r= ((count % 2) == 1)
            if ray == 3 :    
                intersect2, count2 = helper.raycast(geom, point, point+[0.,0.,1.1], diag, count = True )
                intersect3, count3 = helper.raycast(geom, point, point+[0.,1.1,0.], diag, count = True )
                if r :
                   if (count2 % 2) == 1 and (count3 % 2) == 1 :
                       r=True
                   else : 
                       r=False
        if r : # odd inside
            inside = True
        return inside

    def BuildGrid(self, env):
        if self.isOrthogonalBoudingBox==1:
            self.prepare_buildgrid_box(env)
        if env.innerGridMethod =="sdf" and self.isOrthogonalBoudingBox!=1: # A fillSelection can now be a mesh too... it can use either of these methods
            a, b = self.BuildGrid_utsdf(env) # to make the outer most selection from the master and then the compartment
        elif env.innerGridMethod == "bhtree" and self.isOrthogonalBoudingBox!=1:  # surfaces and interiors will be subtracted from it as normal!
            a, b = self.BuildGrid_bhtree(env)
        elif env.innerGridMethod == "jordan" and self.isOrthogonalBoudingBox!=1:  # surfaces and interiors will be subtracted from it as normal!
            a, b = self.BuildGrid_jordan(env)
        elif env.innerGridMethod == "jordan3" and self.isOrthogonalBoudingBox!=1:  # surfaces and interiors will be subtracted from it as normal!
            a, b = self.BuildGrid_jordan(env,ray=3)
        elif env.innerGridMethod == "pyray" and self.isOrthogonalBoudingBox!=1:  # surfaces and interiors will be subtracted from it as normal!
            a, b = self.BuildGrid_pyray(env) 
        elif env.innerGridMethod == 'kevin' and self.isOrthogonalBoudingBox != 1:
            a, b = self.BuildGrid_kevin(env)   
        return a,b 

    def prepare_buildgrid_box(self,env):
        a = env.grid.getPointsInCube(self.bb, None, None) #This is the highspeed shortcut for inside points! and no surface! that gets used if the fillSelection is an orthogonal box and there are no other compartments.
        b = []
        env.grid.gridPtId[a] = -self.number
        self.surfacePointsCoords = None
        bb0x, bb0y, bb0z = self.bb[0]
        bb1x, bb1y, bb1z = self.bb[1]
        AreaXplane = (bb1y-bb0y)*(bb1z-bb0z)
        AreaYplane = (bb1x-bb0x)*(bb1z-bb0z)
        AreaZplane = (bb1y-bb0y)*(bb1x-bb0x)
        vSurfaceArea = abs(AreaXplane)*2+abs(AreaYplane)*2+abs(AreaZplane)*2
        if autopack.verbose :
            print("vSurfaceArea = ", vSurfaceArea)
        self.insidePoints = a
        self.surfacePoints = []
        self.surfacePointsCoords = []
        self.surfacePointsNormals = []
        if autopack.verbose :
            print(' %d inside pts, %d tot grid pts, %d master grid'%( len(a),len(a), len(self.grid.masterGridPositions)))
        self.computeVolumeAndSetNbMol(env, b, a,areas=vSurfaceArea)               
        #print("I've built a grid in the compartment test with no surface", a)
        if autopack.verbose :
            print("The size of the grid I build = ", len(a))
        return a,b,vSurfaceArea

    def BuildGrid_box(self, env,vSurfaceArea):
        nbGridPoints = len(env.grid.masterGridPositions)
        insidePoints = env.grid.getPointsInCube(self.bb, None, 
                                            None,addSP=False)     
        for p in insidePoints : 
            env.grid.gridPtId[p]=-self.number
        if autopack.verbose :
            print('is BOX Total time', time()-t0)
        surfPtsBB,surfPtsBBNorms = self.getSurfaceBB(self.ogsurfacePoints,env)
        srfPts = surfPtsBB
        surfacePoints,surfacePointsNormals  = self.extendGridArrays(nbGridPoints,
                                    srfPts, surfPtsBBNorms, histoVol)
        self.insidePoints = insidePoints
        self.surfacePoints = surfacePoints
        self.surfacePointsCoords = surfPtsBB
        self.surfacePointsNormals = surfacePointsNormals
        if autopack.verbose :
            print('%s surface pts, %d inside pts, %d tot grid pts, %d master grid'%(
            len(self.surfacePoints), len(self.insidePoints),
            nbGridPoints, len(env.grid.masterGridPositions)))
    
        self.computeVolumeAndSetNbMol(env, self.surfacePoints,
                                      self.insidePoints,areas=vSurfaceArea)    
      
    def BuildGrid_kevin(self, env, superFine = False, visualize = False):
        """
        Takes a polyhedron, and builds a grid. In this grid:
            - Projects the polyhedron to the grid.
            - Determines which points are inside/outside the polyhedron
            - Determines point's distance to the polyhedron.
        Usage:
        Select desired object, run this program with desired grid spacing resolution. By default,
        radius is 2 * gridSpacing. superFine provides the option doing a super leakproof test when
        determining which points are inside or outside. This usually not necessary, because the
        built-in algorithm has no known leakage cases. It is simply there as a safeguard.
        """

        def f_ray_intersect_polyhedron(pRayStartPos, pRayEndPos, faces,vertices, pTruncateToSegment):
            """This function returns TRUE if a ray intersects a triangle.
            It also calculates and returns the UV coordinates of said colision as part of the intersection test,
            Makes sure that we are working with arrays

            * TAKES IN pRayStartPos as LEFT-HANDED COORDINATE (Z,Y,X)
            * TAKES IN pRayEndPos AS LEFT-HANDED COORDINATE (Z,Y,X)
            * faces is a list of vertex indices that defines the polyhedron
            * vertices is the list of global coordinates that faces refers to
            * pTruncateToSegment decides whether or not the segment terminates at the end position, or keeps on
              going forever
            """
            pRayStartPos = numpy.array(pRayStartPos)
            pRayEndPos = numpy.array(pRayEndPos)

            vLineSlope = pRayEndPos - pRayStartPos #This line segment defines an infinite line to test for intersection
            vPolyhedronPos = numpy.array((0,0,0))
            vTriPoints = vertices
            vTriPolys = faces
            vBackface = None
            vBackfaceFinal = None
            maxDistance = 99999
            vHitCount = 0
            
            # Globalize the polyhedron
            # Present in original coffee script, but unnecessary because DecomposeMesh gives global coords
            #for i in range(len(vTriPoints)):
            #   vTriPoints[i] = vTriPoints[i] + vPolyhedronPos

            vEpsilon = 0.00001
            vBreakj = False
            vCollidePos = None
            
            # Walk through each polygon in a polyhedron
            for testingFace in vTriPolys:
                #  Loop through all the polygons in an input polyhedron
                #vQuadrangle = 1  #  Default says polygon is a quadrangle.
                #vLoopLimit = 2  #  Default k will loop through polygon assuming its a quad.
                #if (vTriPolys[j+3] == vTriPolys[j+2])  #  Test to see if quad is actually just a triangle.
                #   {
                vQuadrangle = 0  #  Current polygon is not a quad, it's a triangle.
                vLoopLimit = 1  #  Set k loop to only cycle one time.
                
                for k in range(vLoopLimit):
                    vTriPt0 = numpy.array(vTriPoints[testingFace[0]])  # Always get the first point of a quad/tri
                    vTriPt1 = numpy.array(vTriPoints[testingFace[1]])  # Get point 1 for a tri and a quad's first pass, but skip for a quad's second pass
                    vTriPt2 = numpy.array(vTriPoints[testingFace[2]])  # Get point 2 for a tri and a quad's first pass, but get point 3 only for a quad on its second pass.
            
                    vE1 = vTriPt1 - vTriPt0  #  Get the first edge as a vector.
                    vE2 = vTriPt2 - vTriPt0  #  Get the second edge.
                    h = vcross(vLineSlope, vE2)
                    
                    a = f_dot_product(vE1,h)  #  Get the projection of h onto vE1.
                    if a > -vEpsilon and a < vEpsilon:
                        continue# If the ray is parallel to the plane then it does not intersect it, i.e, a = 0 +/- given rounding slope.
                        #  If the polygon is a quadrangle, test the other triangle that comprises it.

                    F = 1.0/a       
                    s = pRayStartPos - vTriPt0  #  Get the vector from the origin of the triangle to the ray's origin.
                    u = F * f_dot_product(s,h)
                    if u < 0.0 or u > 1.0:
                        continue
                        #/* Break if its outside of the triangle, but try the other triangle if in a quad.
                        #U is described as u = : start of vE1 = 0.0,  to the end of vE1 = 1.0 as a percentage.  
                        #If the value of the U coordinate is outside the range of values inside the triangle,
                        #then the ray has intersected the plane outside the triangle.*/

                    q = vcross(s, vE1)
                    v = F * f_dot_product(vLineSlope,q)
                    if v <0.0 or u+v > 1.0:
                        continue  
                        #/*  Break if outside of the triangles v range.
                        #If the value of the V coordinate is outside the range of values inside the triangle,
                        #then the ray has intersected the plane outside the triangle.
                        #U + V cannot exceed 1.0 or the point is not in the triangle.           
                        #If you imagine the triangle as half a square this makes sense.  U=1 V=1 would be  in the 
                        #lower left hand corner which would be in the second triangle making up the square.*/

                    vCollidePos = vTriPt0 + u*vE1 + v*vE2  #  This is the global collision position.
                    assert len(vCollidePos) == 3

                    #  The ray is hitting a triangle, now test to see if its a triangle hit by the ray.
                    vBackface = False
                    if f_dot_product(vLineSlope, vCollidePos - pRayStartPos) > 0:  #  This truncates our infinite line to a ray pointing from start THROUGH end positions.
                        vHitCount += 1
                        #print(testingFace)
                        if pTruncateToSegment and vlen(vLineSlope) < vlen(vCollidePos - pRayStartPos):
                            print('broken')
                            break # This truncates our ray to a line segment from start to end positions.

                        d = squaredTwoPointDistance(pRayStartPos,vCollidePos)
                        if d >= maxDistance:
                            continue
                        if a < vEpsilon:  #  Test to see if the triangle hit is a backface.
                            #set master grid to organelle->getname inside
                            vBackface = True
                            #  This stuff is specific to our Point inside goals.
                            vBreakj = True  #  To see if a point is inside, I can stop at the first backface hit.
                            #break
                        else:
                            vBreakj = True
                        vBackfaceFinal = vBackface
                        maxDistance = d
            #vBackfaceFinal = vBackfaces[distances.index(min(distances))]
            return vHitCount,vBackfaceFinal
        
        def makeGrid(lowerLeft,upperRight,spacing, inclusive = True):
            """
            Given lower left back corner, an upper right front corner, and a grid spacing, construct a list of
            coordinates that representsly an evenly spaced grid spanning the entire volume. Returned coordinates
            preserve the right/lefthandedness of the input coordinates, making this function work equally well
            in both 3D coordinate systems.
            """
            pointCoords = []
            xSize,ySize,zSize = upperRight[0] - lowerLeft[0], upperRight[1] - lowerLeft[1], upperRight[2] - lowerLeft[2]
            if inclusive:
                nx,ny,nz = xSize / spacing + 1, ySize / spacing + 1, zSize / spacing + 1
            else:
                nx,ny,nz = xSize / spacing, ySize / spacing, zSize / spacing
            nx,ny,nz = round(nx),round(ny),round(nz)
            x = numpy.linspace(lowerLeft[0],upperRight[0],nx)
            y = numpy.linspace(lowerLeft[1],upperRight[1],ny)
            z = numpy.linspace(lowerLeft[2],upperRight[2],nz)
            for a in x:
                for b in y:
                    for c in z:
                        pointCoords.append((a,b,c))
            # Order of above had to change to correctly produce a viable BBOX
            return pointCoords, (nx,ny,nz)

        def makeMarchingCube(gridSpacing,r):
            """
            Create a numpy array that represents the precomputed distances to each point
            for the cube of points surrounding our center point.
            """
            def _pythagorean(*edgeLengths):
                from math import sqrt
                total = 0
                for length in edgeLengths:
                    total += length * length
                distance = sqrt(total)
                return distance
            from math import ceil
            pointsForRadius = ceil(r/gridSpacing) # Number of grid points required to represent our radius, rounded up
            pointsInEdge = 2 * pointsForRadius + 1 # Number of points in one edge of our cube

            center = pointsForRadius # The index if our center point
            cube = numpy.zeros(shape=(pointsInEdge,pointsInEdge,pointsInEdge))
            distX = numpy.zeros(shape=(pointsInEdge,pointsInEdge,pointsInEdge))
            distY = numpy.zeros(shape=(pointsInEdge,pointsInEdge,pointsInEdge))
            distZ = numpy.zeros(shape=(pointsInEdge,pointsInEdge,pointsInEdge))
            for a in range(pointsInEdge):
                lenX = a - center
                for b in range(pointsInEdge):
                    lenY = b - center
                    for c in range(pointsInEdge):
                        lenZ = c - center
                        cube[a][b][c] = _pythagorean(lenX,lenY,lenZ) * gridSpacing
                        distX[a][b][c] = lenX
                        distY[a][b][c] = lenY
                        distZ[a][b][c] = lenZ
            return cube, distX, distY, distZ
        
        def getPointFrom3D(pt3d, gridSpacing, gridPtsPerEdge, bb):
            """
            Starts from orthogonal grid. Has a lower left, upper right coordinate from bounding box.
            Given lower left and known spacing of grid, we can provide any vector in 3D space and it will
            algebraically return index of closest grid point.

            Given a RIGHT-HANDED fine coordinate, it will return the index of the closest point in the LEFT-HANDED coarse grid.
            """
            x, y, z = pt3d  # Continuous 3D point to be discretized
            spacing1 = 1./gridSpacing  # Grid spacing = diagonal of the voxel determined by smalled packing radius
            NZ, NY, NX = gridPtsPerEdge  # vector = [length, height, depth] of grid, units = gridPoints
            
            # Added to fix some weird edge behavior
            #NZ, NY, NX = round(NZ),round(NY),round(NX)
            #print(NX,NY,NZ)
            OZ, OY, OX = bb[0] # origin of Pack grid. OX and OZ were switched to compensate for coordinate systems
            #print(OX,OY,OZ)
            
            # Algebra gives nearest gridPoint ID to pt3D
            i = min( NX-1, max( 0, round((x-OX)*spacing1)))
            j = min( NY-1, max( 0, round((y-OY)*spacing1)))
            k = min( NZ-1, max( 0, round((z-OZ)*spacing1)))
            result = k*NX*NY + j*NX + i
            return int(result)

        def squaredTwoPointDistance(point1,point2):
            """Computes and returns the distance between two points in 3D space"""
            dists = [numpy.float64(point1[i] - point2[i]) for i in range(len(point1))]
            distsSquare = [x**2 for x in dists]
            dist = sum(distsSquare)
            return dist

        def findPointsCenter(*args):
            # Average down the column, such that we're averaging across all measurements in one dimension
            center = numpy.mean(args[0], axis = 0)
            return center

        def vcross(v1,v2):
            """
            Returns the cross-product of two 3D vectors.
            """
            x1,y1,z1 = v1
            x2,y2,z2 = v2
            return (y1*z2-y2*z1, z1*x2-z2*x1, x1*y2-x2*y1)

        def vlen(v1):
            """
            Returns the length of a 3D vector.
            """
            x1,y1,z1 = v1
            return sqrt(x1*x1 + y1*y1 + z1*z1)

        def f_dot_product(vector1,vector2):
            """
            Return the dot product of two 3D vectors.
            """
            dottedVectors = [vector1[i] * vector2[i] for i in range(len(vector1))]
            return sum(dottedVectors)

        def visualizeAtomArray(name = 'Test', vertexCoordinates = [(0,0,0)], sphereRadius = 5):
            """
            Given coordinates, visualize them in Cinema 4D using an Atom Array Object.
            """
            atom = c4d.BaseObject(c4d.Oatomarray)
            atom.SetName(name)
            atom[1000] = 3 #radius cylinder
            atom[1001] = sphereRadius #radius sphere
            atom[1002] = 5 #subdivision
            helper.addObjectToScene(helper.getCurrentScene(),atom,parent=None)
            
            parent = atom
            coords = vertexCoordinates
            nface = 0
            obj= c4d.PolygonObject(len(coords), nface)
            obj.SetName(name + '_child')
            cd4vertices = map(helper.FromVec,coords)
            map(obj.SetPoint,range(len(coords)),cd4vertices)
            helper.addObjectToScene(helper.getCurrentScene(),obj,parent=parent)

        class gridPoint:
            def __init__(self,i,globalC,isPolyhedron):
                self.index = int(i)
                self.isOutside = None
                self.minDistance = 99999 # Only store a number here if within certain distance from polyhedron
                self.representsPolyhedron = isPolyhedron
                self.closeFaces = []
                self.closestFaceIndex = 0
                self.testedEndpoint = None
                self.allDistances = [] # Stores a tuple list of distances to all points. (point,distance) = (5,2.5)
                self.globalCoord = numpy.array(globalC) # Stores the global coordinate associated with this point

        if self.ghost : return
        from time import time
        # Copied from jordan
        t0 = t1 = time()

        # Get the object, and decompose it (standalone version)
        # object_target = helper.getCurrentSelection()[0]
        # faces,verticesLH,vnormals,faceNormals = helper.DecomposeMesh(object_target,edit=False,copy=False,tri=True,transform=True,fn=True)
        # This has been replaced with the below, upon integration with Compartment.py
        verticesLH = self.vertices
        faces = self.faces
        normalList2, areas = self.getFaceNormals(verticesLH, faces,fillBB=env.fillBB)
        vSurfaceArea = sum(areas)

        distances = env.grid.distToClosestSurf
        idarray = env.grid.gridPtId
        diag = env.grid.diag

        t1 = time()

        # Number is the compartment ID of our current polyhedron. This will be important later on in the idarray magic.
        number = self.number

        # Inferred from reading Compartment.py
        # This is incorrect. THis seems to return the diagnonal distance, not the straight
        # gridSpacing = env.grid.gridSpacing

        startTime = time()
        from math import ceil
        # Setup a bounding box for our polyhedron. Deprecated upon integration into Compartment.py
        # b = setupBBOX(object_target)
        # corners = b.getBoundingBox()
        corners = self.bb

        # Create grid points to fill in our bbox with grid points at regular distances
        pointsRH = env.grid.masterGridPositions #This returns right-handed coords. We need left-handed coords, so we'll flip them.
        points = [x[::-1] for x in pointsRH]
        gridPtsPerEdgeRH = env.grid.nbGridPoints # Also right-handed. So we'll need to flip these as well.
        gridPtsPerEdge = gridPtsPerEdgeRH[::-1]
        gridSpacing = max(points[1] - points[0])
        radius = gridSpacing * 1

        # Create a gridPoint object for every single point we have in our grid. This pre-allocates all memory, so that future
        # computations are faster.
        gridPoints = []
        i = 0
        for point in points:
            gridPoints.append(gridPoint(i,point,isPolyhedron = False))
            i += 1
        # assert len(gridPoints) == len(points)

        # Make a precomputed distance cube, where 
        distanceCube,distX,distY,distZ = makeMarchingCube(gridSpacing,radius)
        # Flatten and combine these arrays. This makes it easier to iterate over.
        distanceCubeF,distXF,distYF,distZF = distanceCube.flatten(),distX.flatten(),distY.flatten(),distZ.flatten()
        zippedNumbers = zip(distanceCubeF,distXF,distYF,distZF)

        NZ,NY,NX = gridPtsPerEdge
        OZ, OY, OX = self.bb[0][::-1] # It looks like the bb is in right-handed coords, so we'll have to reverse it
        spacing1 = 1./gridSpacing # Inverse of the spacing. We compute this here, so we don't have to recompute it repeatedly
        allCoordinates = [] # Tracker for all the fine coordiantes that we have interpolated for the faces of the polyhedron
        # pointsToTestInsideOutside = set()
        # Walk through the faces, projecting each to the grid and marking immediate neighbors so we can test said 
        # neighbors for inside/outside later.
        for face in faces:
            # Dist is a small function that calculates the distance between two points.
            # dist = lambda x,y: vlen([y[i] - x[i] for i in range(len(x))])
            # Get the vertex coordinates and conver to numpy arrays...
            triCoords = [numpy.array(verticesLH[i]) for i in face]
            thisFaceFineCoords = list(triCoords)
            allCoordinates.extend(triCoords)
            # We will use these u/v vectors to interpolate points that reside on the face
            pos = triCoords[0]
            u = triCoords[1] - pos
            v = triCoords[2] - pos
            # Smetimes the hypotenuse isn't fully represented. To remedy this, we will use an additional w vector
            # to interpolate points on the hypotenuse
            w = triCoords[2] - triCoords[1]

            # If either u or v is greater than the grid spacing, then we need to subdivide it
            # We will use ceil: if we have a u of length 16, and grid spacing of 5, then we want
            # a u at 0, 5, 10, 15 which is [0, 1, 2, 3] * gridSpacing. We need ceil to make this happen.
            
            # This is a much more efficient solution than the previous method of artificially increasing
            # the densith of the polygon mesh by 4, then mapping each of those small points to the grid.
            # However, it is also leaks more, for (as of yet) unknown reasons.
            
            # It was discovered that by using the default gridspacing, some faces will produce leakage. One solution is to
            # temporarily use a somewhat denser gridspacing to interpolate, and then project back onto our original spacing.
            # We'll decrease the gridspacing by 25% (so that it's 75% of the original). This seems to fix the issue without
            # any appreciable increase in computation time.
            gridSpacingTempFine = gridSpacing / 3
            # Determine the number of grid spacing-sized points we can fit on each vector.
            # Minimum is one because range(1) gives us [0]
            uSubunits, vSubunits, wSubunits = 1, 1, 1
            if vlen(u) > gridSpacingTempFine:
                uSubunits = ceil(vlen(u)/gridSpacingTempFine) + 1
            if vlen(v) > gridSpacingTempFine:
                vSubunits = ceil(vlen(v)/gridSpacingTempFine) + 1
            if vlen(w) > gridSpacingTempFine:
                wSubunits = ceil(vlen(w)/gridSpacingTempFine) + 1
            # Because we have observed leakage, maybe we want to try trying a denser interpolation, using numpy's linspace?
            for uSub in range(uSubunits):
                percentU = uSub * gridSpacingTempFine / vlen(u)
                percentU = min(percentU, 1.0) # Make sure that we have not stepped outside of our original u vector
                # h represents the height of the hypotenuse at this u. Naturally, we cannot go past the hypotenuse, so this will be
                # our upper bound.
                h = percentU * u + (1 - percentU) * v
                for vSub in range(vSubunits):
                    percentV = vSub * gridSpacingTempFine / vlen(v)
                    percentV = min(percentV, 1.0) # Make sure that we have not stepped oustide of our original v vector.
                    interpolatedPoint = percentU * u + percentV * v
                    # The original if: statement asks if the distance from the origin to the interpolated point is less than
                    # the distance from the origin to the hypotenuse point, as such:
                    # if vlen(interpolatedPoint) < vlen(h):
                    # Wouldn't it be a better idea to measure distance to the u position instead? This is implemented below.
                    if (vlen(interpolatedPoint - percentU * u) < vlen(h - percentU * u)):
                        allCoordinates.append(interpolatedPoint + pos)
                        thisFaceFineCoords.append(interpolatedPoint + pos)
                    else:
                        break
            # Because the above only interpolates the face, and may not completely capture the hypotenuse, let's separately
            # interpolate points on the hypotenuse (the w vector). This way we can be 100% sure that our final solution is 
            # completely leak-proof
            for wSub in range(wSubunits):
                # Apply the same proceudre we did above for u/v, just for w (for hypotenuse interpolation)
                percentW = wSub * gridSpacingTempFine / vlen(w)
                percentW = min(percentW, 1.0)
                interpolatedPoint = percentW * w
                allCoordinates.append(interpolatedPoint + triCoords[1])
                thisFaceFineCoords.append(interpolatedPoint + triCoords[1])
            projectedIndices = set()
            # Once we have interpolated the face, let's project each interpolated point to the grid.
            for coord in thisFaceFineCoords:
                projectedPointIndex = getPointFrom3D(coord[::-1],gridSpacing,gridPtsPerEdge,self.bb)
                projectedIndices.add(projectedPointIndex)

            # Walk through each grid point that our face spans, gather its closest neighbors, and annotate the neighbors with
            # minimum distance and closest faces, and flag them for testing inside/outside later.
            for P in list(projectedIndices):
                # Get the point object corresponding to the index, and set its polyhedron attribute to true
                g = gridPoints[P]
                g.representsPolyhedron = True
                # Get the coordinates of the point, and convert them to grid units
                zTemp,yTemp,xTemp = g.globalCoord
                i,j,k = round((xTemp-OX)*spacing1), round((yTemp-OY)*spacing1), round((zTemp-OZ)*spacing1)
                # Let's step through our distance cube, and assign faces/closest distances to each
                for d,x,y,z in zippedNumbers:
                    # Get the grid indices for the point we're considering, and skip it if we're stepping oustide the boundaries
                    newI, newJ, newK = i + x, j + y, k + z
                    if newI < 0 or newI > (NX-1) or newJ < 0 or newJ > (NY-1) or newK < 0 or newK > (NZ-1):
                        continue
                    # Get the point index that this coordinate corresponds to.
                    desiredPointIndex = int(round(newK*NX*NY + newJ*NX + newI))
                    # if desiredPointIndex == 54199:
                    #     print('HIT with face ' + str(face) + ' and polygon point ' + str(P))
                    desiredPoint = gridPoints[desiredPointIndex]
                    if desiredPoint.representsPolyhedron == True:
                        continue
                    # Add the current face to the its list of closest faces
                    if face not in desiredPoint.closeFaces:
                        desiredPoint.closeFaces.append(face)
                    # Add the distance to the point's list of distances, and overwrite minimum distance if appropriate
                    desiredPoint.allDistances.append((v,d))
                    if d < desiredPoint.minDistance:
                        desiredPoint.minDistance = d
                    # Later down the road, we want to test as few points as possible for inside/outside. Therefore,
                    # we will only test points that are 
                    # if abs(x) <= 1 and abs(y) <= 1 and abs(z) <= 1:
                    #     pointsToTestInsideOutside.add(desiredPointIndex)
        # visualizeAtomArray('fineCoords',allCoordinates,5)
        # projectedIndices = [x.index for x in gridPoints if x.representsPolyhedron == True]
        # projectedIndices = set()
        # for coord in allCoordinates:
        #     projectedPointIndex = getPointFrom3D(coord[::-1],gridSpacing,gridPtsPerEdge,b)
        #     projectedIndices.add(projectedPointIndex)
        # visualizeAtomArray('projectedCoords',[points[i] for i in projectedIndices],5)
        
        # ogsurfacePoints represents the off-grid surface points that have been interpolated to fit our grid.
        # However, allCoordinates is in left-handed coordinates. We'll have to flip it since it seems that
        # the self.ogsurfacePoints works in right-handed coordinates.
        if self.isBox :
            self.overwriteSurfacePts = True
        if self.overwriteSurfacePts:
            self.ogsurfacePoints = self.vertices[:]
            self.ogsurfacePointsNormals = self.vnormals[:]
        else :
            # self.createSurfacePoints(maxl=env.grid.gridSpacing)
            self.ogsurfacePoints = [x[::-1] for x in allCoordinates]
        print('******************************************************************')
        print(allCoordinates)
        return None

        srfPts = self.ogsurfacePoints
        timeFinishProjection = time()
        print('Projecting polyhedron to grid took ' + str(timeFinishProjection - startTime) + ' seconds.')
        
        # Let's start filling in inside outside. Here's the general algorithm:
        # Walk through all the points in our grid. Once we encounter a point that has closest faces, 
        # then we know we need to test it for inside/outside. Once we test that for inside/outside, we
        # fill in all previous points with that same inside outisde property. To account for the possible
        # situation that there is a surface that is only partially bound by the bbox, then we need to
        # reset the insideOutsideTracker every time we have a change in more than 1 of the 3 coordinates
        # because that indicates we're starting a new row/column of points.

        # This tracks our inside/outside.
        isOutsideTracker = None
        # This tracks the points that we've iterated over which we do not know if inside/outside. Remember
        # to reset every time we find an inside/outside.
        emptyPointIndicies = []
        mismatchCounter = 0
        for g in gridPoints:        
            # Check if we've started a new line. If so, then we reset everything.
            # This test should precede all other test, because we don't want old knowldge
            # to carry over to the new line, since we don't know if the polygon is only partially encapsulated by the bounding box.
            if g.index > 0: # We can't check the first element, so we can skip it. 
                coordDiff = g.globalCoord - gridPoints[g.index - 1].globalCoord
                coordDiffNonzero = [x != 0 for x in coordDiff]
                if sum(coordDiffNonzero) > 1:
                    # assert len(emptyPointIndicies) == 0 # When starting a new line, we shouldn't have any unknowns from the previous line
                    isOutsideTracker = None
                    emptyPointIndicies = []

            # There's no point testing inside/outside for points that are on the surface. So if
            # we get to one, we just skip over it.
            if g.representsPolyhedron == True:
                g.isOutside = None
                continue

            if len(g.closeFaces) == 0:
                # If it's not close to any faces, and we don't know if this row is inside/outside, then
                # we have to wait till later to figure it out
                if isOutsideTracker == None:
                    emptyPointIndicies.append(g.index)
                # However, if we do know , we can just use the previous one to fill
                else:
                    g.isOutside = isOutsideTracker        
            # If there are close faces attached to it, then we need to test it for inside/outside.
            else:
                # Find centroid of all the vertices of all the close faces. This will be our endpoint
                # when casting a ray for collision testing. 
                uniquePoints = []
                # This takes just the first face and projects to the center of it.
                # [uniquePoints.append(x) for x in g.closeFaces[0] if x not in uniquePoints]
                [uniquePoints.append(x) for x in g.closeFaces[g.closestFaceIndex] if x not in uniquePoints]
                uniquePointsCoords = verticesLH[uniquePoints]
                endPoint = findPointsCenter(uniquePointsCoords)
                g.testedEndpoint = endPoint

                # Draw a ray to that point, and see if we hit a backface or not
                numHits, thisBackFace = f_ray_intersect_polyhedron(g.globalCoord,endPoint,g.closeFaces,verticesLH,False)
                
                # We can check the face as well if we want to be super precise. If they dont' agree, we then check against the entire polyhedron.
                if superFine == True:
                    if len(g.closeFaces) > 1:
                        uniquePoints2 = []
                        [uniquePoints2.append(x) for x in g.closeFaces[1] if x not in uniquePoints2]
                        uniquePointsCoords2 = verticesLH[uniquePoints2]
                        endPoint2 = findPointsCenter(uniquePointsCoords2)
                        numHits2, thisBackFace2 = f_ray_intersect_polyhedron(g.globalCoord,endPoint2,g.closeFaces,verticesLH,False)
                    if len(g.closeFaces) == 1 or thisBackFace != thisBackFace2:
                        mismatchCounter += 1
                        numHits, thisBackFace = f_ray_intersect_polyhedron(g.globalCoord,numpy.array([0.0,0.0,0.0]),faces,verticesLH,False)
                
                # Fill in inside outside attribute for this point, as pRayStartPos, pRayEndPos, faces, vertices, pTruncateToSegmentll as for any points before it
                g.isOutside = not thisBackFace
                isOutsideTracker = not thisBackFace
                for i in emptyPointIndicies:
                    gridPoints[i].isOutside = isOutsideTracker
                # Because we have filled in all the unknowns, we can reset that counter.
                emptyPointIndicies = []

        # Final pass through for quality/sanity checks.
        for g in gridPoints:
            if g.representsPolyhedron == True:
                assert g.isOutside == None
                continue
            else:
                if g.isOutside == None:
                    g.isOutside = True
        print('Flood filling grid inside/outside took ' + str(time() - timeFinishProjection) + ' seconds.')
        insidePoints = [g.index for g in gridPoints if g.isOutside == False]
        # Carries out the idarray magic.
        for pt in insidePoints:
            idarray.itemset(pt, -number)
        outsidePoints = [g.index for g in gridPoints if g.isOutside == True]
        surfacePoints = [g.index for g in gridPoints if g.representsPolyhedron == True]

        if visualize == True:
            visualizeAtomArray('insidePoints',[points[i] for i in insidePoints],0.5)
            visualizeAtomArray('outsidePoints',[points[i] for i in outsidePoints],0.5)
            visualizeAtomArray('surfacePoints',[points[i] for i in surfacePoints],0.5)
        if superFine == True:
            print('Superfine was ' + str(superFine) + ' and there were ' + str(mismatchCounter) + ' mismatches.')
        print('Grid construction took ' + str(time() - startTime) + ' seconds for ' + str(len(faces)) + ' faces and ' + str(len(gridPoints)) + ' points.')
        
        nbGridPoints = len(env.grid.masterGridPositions)
        surfPtsBB,surfPtsBBNorms = self.getSurfaceBB(srfPts,env)
        srfPts = surfPtsBB

        ex = True
        surfacePoints,surfacePointsNormals = self.extendGridArrays(nbGridPoints,
                                            srfPts,surfPtsBBNorms,
                                            env,extended=ex)

        self.insidePoints = insidePoints
        self.surfacePoints = surfacePoints
        self.surfacePointsCoords = surfPtsBB
        self.surfacePointsNormals = surfacePointsNormals
        self.computeVolumeAndSetNbMol(env, self.surfacePoints,self.insidePoints,areas=vSurfaceArea)

        with open('C:/Users/Kevin/Mesoscope/Self Learning/testing.txt', 'w+') as output:
            # pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
            for attribute in dir(self):
                attrValue = str(getattr(self, attribute))
                output.write('\n' + attribute + attrValue + '\n')
        output.close()

        return self.insidePoints, self.surfacePoints

        # print('My output:')
        # print(surfacePoints)
        # print(len(surfacePoints))
        # print('Jordan Output')
        # test = self.BuildGrid_jordan(env)
        # print(test[1])
        # print(len(test[1]))
        # print("ARE THEY THE SAME")
        # print(set(surfacePoints) == test[1])
        # print("THIS IS THE INTERSECTION")
        # print(set(surfacePoints).intersection(test[1]))
        # print(len(set(surfacePoints).intersection(test[1])))
        # print("END TEST")
        # return None
        # return self.insidePoints, self.surfacePoints


    #Jordan Curve Theorem
    def BuildGrid_jordan(self, env,ray=1):
        """Build the compartment grid ie surface and inside point using jordan theorem and host raycast"""
        # create surface points 
        if self.ghost : return
        t0 = t1 = time()        
        if self.isBox :
            self.overwriteSurfacePts = True
        if self.overwriteSurfacePts:
            self.ogsurfacePoints = self.vertices[:]
            self.ogsurfacePointsNormals = self.vnormals[:]
        else :
            self.createSurfacePoints(maxl=env.grid.gridSpacing)
        
        # Graham Sum the SurfaceArea for each polyhedron
        vertices = self.vertices  #NEED to make these limited to selection box, not whole compartment
        faces = self.faces #         Should be able to use self.ogsurfacePoints and collect faces too from above
        normalList2,areas = self.getFaceNormals(vertices, faces,fillBB=env.fillBB)
        vSurfaceArea = sum(areas)
        
        # build a BHTree for the vertices
        if self.isBox :
            self.BuildGrid_box(env,vSurfaceArea)
            return self.insidePoints, self.surfacePoints
        
        if autopack.verbose :
            print('time to create surface points', 
                  time()-t1, len(self.ogsurfacePoints))

        distances = env.grid.distToClosestSurf
        idarray = env.grid.gridPtId
        diag = env.grid.diag
        if autopack.verbose :
            print ("distance ",len(distances),
                   "idarray ",len(idarray))
        t1 = time()

        #build BHTree for off grid surface points
        #or scipy ckdtree ?
#        from bhtree import bhtreelib
        from scipy import spatial
        srfPts = self.ogsurfacePoints
#        self.OGsrfPtsBht = bht =  bhtreelib.BHtree(tuple(srfPts), None, 10)
        self.OGsrfPtsBht = bht = spatial.cKDTree(tuple(srfPts),leafsize=10) 
        #res = numpy.zeros(len(srfPts),'f')
        #dist2 = numpy.zeros(len(srfPts),'f')

        number = self.number
        #ogNormals = self.ogsurfacePointsNormals
        insidePoints = []

        # find closest off grid surface point for each grid point 
        #FIXME sould be diag of compartment BB inside fillBB
        grdPos = env.grid.masterGridPositions
#        returnNullIfFail = 0
        if autopack.verbose :
            print ("compartment build grid jordan",
                   diag," nb points in grid ",len(grdPos))#[],None
#        closest = bht.closestPointsArray(tuple(grdPos), diag, 
#                                         returnNullIfFail)
        closest = bht.query(tuple(grdPos))#return both indices and distances

        helper = autopack.helper
#        geom = helper.getObject(self.gname)      
#        if geom is None :
#            self.gname = '%s_Mesh'%self.name            
#            geom = helper.getObject(self.gname)
#        if geom is not None :
#            center = helper.getTranslation( geom )
#        else :
#            self.getCenter()        
#            center = self.center        

        self.closestId = closest[1]
        new_distances =  closest[0]
        mask = distances[:len(grdPos)] > new_distances
        nindices  = numpy.nonzero(mask)
        distances[nindices] = new_distances[nindices]
        #now check if point inside 
        
        helper.resetProgressBar()
        for ptInd in xrange(len(grdPos)):#len(grdPos)):
            coord = [grdPos.item((ptInd,0)),grdPos.item((ptInd,1)),grdPos.item((ptInd,2))]
            insideBB  = self.checkPointInsideBB(coord,dist=new_distances.item(ptInd))
            r=False
            if insideBB:
                r=self.checkPointInside_rapid(coord,diag,ray=ray)
            if r : # odd inside
                insidePoints.append(ptInd)
                idarray.itemset(ptInd,-number)
                #idarray[ptInd] = -number
            p=(ptInd/float(len(grdPos)))*100.0 
            if (ptInd % 100 ) == 0 :
                helper.progressBar(progress=int(p),label=str(ptInd)+"/"+str(len(grdPos))+" inside "+str(r))
                if autopack.verbose:
                    print (str(ptInd)+"/"+str(len(grdPos))+" inside "+str(r))
        if autopack.verbose:
            print('time to update distance field and idarray', time()-t1)
        
        t1 = time()
        nbGridPoints = len(env.grid.masterGridPositions)
        
        surfPtsBB,surfPtsBBNorms = self.getSurfaceBB(srfPts,env)
        srfPts = surfPtsBB
        if autopack.verbose:
            print ("compare length id distances",nbGridPoints,len(idarray),len(distances),
                   (nbGridPoints == len(idarray)),
                        (nbGridPoints == len(distances)))
        ex = True#True if nbGridPoints == len(idarray) else False
        surfacePoints,surfacePointsNormals  = self.extendGridArrays(nbGridPoints,
                                            srfPts,surfPtsBBNorms,
                                            env,extended=ex)


        insidePoints = insidePoints
        if autopack.verbose:
            print('time to extend arrays', time()-t1)
            print('Total time', time()-t0)

        self.insidePoints = insidePoints
        self.surfacePoints = surfacePoints
        self.surfacePointsCoords = surfPtsBB
        self.surfacePointsNormals = surfacePointsNormals
        if autopack.verbose:
            print('%s surface pts, %d inside pts, %d tot grid pts, %d master grid'%(
                len(self.surfacePoints), len(self.insidePoints),
                nbGridPoints, len(env.grid.masterGridPositions)))

        self.computeVolumeAndSetNbMol(env, self.surfacePoints,
                                      self.insidePoints,areas=vSurfaceArea)
        return self.insidePoints, self.surfacePoints


    def BuildGrid_pyray(self, histoVol,ray = 1):
        """Build the compartment grid ie surface and inside point using bhtree"""
        # create surface points 
        if self.ghost : return
        t0 = t1 = time()        
        if self.isBox :
            self.overwriteSurfacePts = True
        if self.overwriteSurfacePts:
            self.ogsurfacePoints = self.vertices[:]
            self.ogsurfacePointsNormals = self.vnormals[:]
        else :
            self.createSurfacePoints(maxl=histoVol.grid.gridSpacing)
        
        # Graham Sum the SurfaceArea for each polyhedron
        vertices = self.vertices  #NEED to make these limited to selection box, not whole compartment
        faces = self.faces #         Should be able to use self.ogsurfacePoints and collect faces too from above
        normalList2,areas = self.getFaceNormals(vertices, faces,fillBB=histoVol.fillBB)
        vSurfaceArea = sum(areas)
        #for gnum in range(len(normalList2)):
        #    vSurfaceArea = vSurfaceArea + areas[gnum]
        
        # build a BHTree for the vertices
        if self.isBox :
            nbGridPoints = len(histoVol.grid.masterGridPositions)
            insidePoints = histoVol.grid.getPointsInCube(self.bb, None, 
                                                None,addSP=False)     
            for p in insidePoints : histoVol.grid.gridPtId[p]=-self.number
            print('is BOX Total time', time()-t0)
            surfPtsBB,surfPtsBBNorms = self.getSurfaceBB(self.ogsurfacePoints,histoVol)
            srfPts = surfPtsBB
            surfacePoints,surfacePointsNormals  = self.extendGridArrays(nbGridPoints,srfPts,
                        surfPtsBBNorms,histoVol)
            self.insidePoints = insidePoints
            self.surfacePoints = surfacePoints
            self.surfacePointsCoords = surfPtsBB
            self.surfacePointsNormals = surfacePointsNormals
            print('%s surface pts, %d inside pts, %d tot grid pts, %d master grid'%(
                len(self.surfacePoints), len(self.insidePoints),
                nbGridPoints, len(histoVol.grid.masterGridPositions)))
        
            self.computeVolumeAndSetNbMol(histoVol, self.surfacePoints,
                                          self.insidePoints,areas=vSurfaceArea)    
            print('time to create surface points', time()-t1, len(self.ogsurfacePoints))
            return self.insidePoints, self.surfacePoints
        
        print('time to create surface points', time()-t1, len(self.ogsurfacePoints))

        distances = histoVol.grid.distToClosestSurf
        idarray = histoVol.grid.gridPtId
        diag = histoVol.grid.diag

        t1 = time()

        #build BHTree for off grid surface points
        from bhtree import bhtreelib
        from scipy import spatial
        srfPts = self.ogsurfacePoints
        #??why the bhtree behave like this
        self.OGsrfPtsBht = bht =  bhtreelib.BHtree(tuple(srfPts), None, 10)
        ctree = spatial.cKDTree(srfPts, leafsize=10)
        res = numpy.zeros(len(srfPts),'f')
        dist2 = numpy.zeros(len(srfPts),'f')
#        print "nspt",len(srfPts)
#        print srfPts
        number = self.number
        ogNormals = self.ogsurfacePointsNormals
        insidePoints = []

        # find closest off grid surface point for each grid point 
        #FIXME sould be diag of compartment BB inside fillBB
        grdPos = histoVol.grid.masterGridPositions
        returnNullIfFail = 0
        print ("compartment build grid ",grdPos,"XX",diag,"XX",len(grdPos))#[],None
#        closest = bht.closestPointsArray(tuple(grdPos), diag, returnNullIfFail)
        new_distance,nb = ctree.query(grdPos,1)
        wh = numpy.greater(distances,new_distance)

        self.closestId = nb
        numpy.copyto(distances,new_distance,where=wh)
        #how to get closest triangle ? may help
        print ("distance have been updated")
        helper = autopack.helper
        geom = self.mesh
        center = helper.getCenter(srfPts)
        print ("mesh center",center)
        for ptInd in range(len(grdPos)):
            inside = False
            insideBB  = self.checkPointInsideBB(grdPos[ptInd],dist=distances[ptInd])
            #print ("insideBB",ptInd,insideBB)
            r=False
            if insideBB:
                #should use an optional direction for the ray, which will help for unclosed surface....
                intersect, count = helper.raycast(geom, grdPos[ptInd], center, diag, count = True ,
                                                  vertices=self.vertices, faces=self.faces)
                #intersect, count = helper.raycast(geom, grdPos[ptInd], grdPos[ptInd]+[0.,1.,0.], diag, count = True )
                r= ((count % 2) == 1)
                if ray == 3 :
                    intersect2, count2 = helper.raycast(geom, grdPos[ptInd], grdPos[ptInd]+[0.,0.,1.1], diag, count = True,
                                                        vertices=self.vertices, faces=self.faces)
                    center = helper.rotatePoint(helper.ToVec(center),[0.,0.,0.],[1.0,0.0,0.0,math.radians(33.0)])
                    intersect3, count3 = helper.raycast(geom, grdPos[ptInd], grdPos[ptInd]+[0.,1.1,0.], diag, count = True,
                                                        vertices=self.vertices, faces=self.faces)
                    #intersect3, count3 = helper.raycast(geom, grdPos[ptInd], center, diag, count = True )#grdPos[ptInd]+[0.,1.1,0.]
                    if r :
                       if (count2 % 2) == 1 and (count3 % 2) == 1 :
                           r=True
                       else : 
                           r=False
            if r : # odd inside
                inside = True
                if inside :
                    insidePoints.append(ptInd)
                    idarray[ptInd] = -number
            p=(ptInd/float(len(grdPos)))*100.0
            helper.progressBar(progress=int(p),label=str(ptInd)+"/"+str(len(grdPos))+" inside "+str(inside))
        print('time to update distance field and idarray', time()-t1)

        t1 = time()
        nbGridPoints = len(histoVol.grid.masterGridPositions)
        
        surfPtsBB,surfPtsBBNorms = self.getSurfaceBB(srfPts,histoVol)
        srfPts = surfPtsBB
        print ("compare length id distances",nbGridPoints,(nbGridPoints == len(idarray)),(nbGridPoints == len(distances)))
        ex = True#True if nbGridPoints == len(idarray) else False
        #back to list type
        #histoVol.grid.distToClosestSurf = histoVol.grid.distToClosestSurf.tolist()
        surfacePoints,surfacePointsNormals  = self.extendGridArrays(nbGridPoints,
                                            srfPts,surfPtsBBNorms,histoVol,extended=ex)

        insidePoints = insidePoints
        print('time to extend arrays', time()-t1)

        print('Total time', time()-t0)

        self.insidePoints = insidePoints
        self.surfacePoints = surfacePoints
        self.surfacePointsCoords = surfPtsBB
        self.surfacePointsNormals = surfacePointsNormals
        print('%s surface pts, %d inside pts, %d tot grid pts, %d master grid'%(
            len(self.surfacePoints), len(self.insidePoints),
            nbGridPoints, len(histoVol.grid.masterGridPositions)))

        self.computeVolumeAndSetNbMol(histoVol, self.surfacePoints,
                                      self.insidePoints,areas=vSurfaceArea)
        bhtreelib.freeBHtree(bht)
        return self.insidePoints, self.surfacePoints
        
    def BuildGrid_bhtree(self, histoVol):
        """Build the compartment grid ie surface and inside point using bhtree"""
        # create surface points 
        if self.ghost : return
        t0 = t1 = time()        
        if self.isBox :
            self.overwriteSurfacePts = True
        if self.overwriteSurfacePts:
            self.ogsurfacePoints = self.vertices[:]
            self.ogsurfacePointsNormals = self.vnormals[:]
        else :
            self.createSurfacePoints(maxl=histoVol.grid.gridSpacing)
        
        # Graham Sum the SurfaceArea for each polyhedron
        vertices = self.vertices  #NEED to make these limited to selection box, not whole compartment
        faces = self.faces #         Should be able to use self.ogsurfacePoints and collect faces too from above
        normalList2,areas = self.getFaceNormals(vertices, faces,fillBB=histoVol.fillBB)
        vSurfaceArea = sum(areas)
        #for gnum in range(len(normalList2)):
        #    vSurfaceArea = vSurfaceArea + areas[gnum]
        
        # build a BHTree for the vertices
        if self.isBox :
            nbGridPoints = len(histoVol.grid.masterGridPositions)
            insidePoints = histoVol.grid.getPointsInCube(self.bb, None, 
                                                None,addSP=False)     
            for p in insidePoints : histoVol.grid.gridPtId[p]=-self.number
            print('is BOX Total time', time()-t0)
            surfPtsBB,surfPtsBBNorms = self.getSurfaceBB(self.ogsurfacePoints,histoVol)
            srfPts = surfPtsBB
            surfacePoints,surfacePointsNormals  = self.extendGridArrays(nbGridPoints,srfPts,
                        surfPtsBBNorms,histoVol)
            self.insidePoints = insidePoints
            self.surfacePoints = surfacePoints
            self.surfacePointsCoords = surfPtsBB
            self.surfacePointsNormals = surfacePointsNormals
            print('%s surface pts, %d inside pts, %d tot grid pts, %d master grid'%(
                len(self.surfacePoints), len(self.insidePoints),
                nbGridPoints, len(histoVol.grid.masterGridPositions)))
        
            self.computeVolumeAndSetNbMol(histoVol, self.surfacePoints,
                                          self.insidePoints,areas=vSurfaceArea)    
            print('time to create surface points', time()-t1, len(self.ogsurfacePoints))
            return self.insidePoints, self.surfacePoints
        
        print('time to create surface points', time()-t1, len(self.ogsurfacePoints))

        distances = histoVol.grid.distToClosestSurf
        idarray = histoVol.grid.gridPtId
        diag = histoVol.grid.diag

        t1 = time()

        #build BHTree for off grid surface points
        from bhtree import bhtreelib
        srfPts = self.ogsurfacePoints
        #??why the bhtree behave like this
        self.OGsrfPtsBht = bht =  bhtreelib.BHtree(tuple(srfPts), None, 10)
        res = numpy.zeros(len(srfPts),'f')
        dist2 = numpy.zeros(len(srfPts),'f')
#        print "nspt",len(srfPts)
#        print srfPts
        number = self.number
        ogNormals = self.ogsurfacePointsNormals
        insidePoints = []

        # find closest off grid surface point for each grid point 
        #FIXME sould be diag of compartment BB inside fillBB
        grdPos = histoVol.grid.masterGridPositions
        returnNullIfFail = 0
        print ("compartment build grid ",grdPos,"XX",diag,"XX",len(grdPos))#[],None
        closest = bht.closestPointsArray(tuple(grdPos), diag, returnNullIfFail)
        
#        print len(closest),diag,closest[0]
#        print closest
#        bhtreelib.freeBHtree(bht)
        self.closestId = closest
#        print len(self.closestId)
#would it be faster using c4d vector ? or hots system?
##        import c4d
#        c4d.StatusSetBar(0)
        for ptInd in range(len(grdPos)):#len(grdPos)):
            # find closest OGsurfacepoint
            gx, gy, gz = grdPos[ptInd]
            sptInd = closest[ptInd]
            if closest[ptInd]==-1:
                print("ouhoua, closest OGsurfacePoint = -1")
                #pdb.set_trace()#???  Can be used to debug with http://docs.python.org/library/pdb.html
            if sptInd < len(srfPts):
                sx, sy, sz = srfPts[sptInd]
                d = sqrt( (gx-sx)*(gx-sx) + (gy-sy)*(gy-sy) +
                          (gz-sz)*(gz-sz))
            else :
                try :
                    n=bht.closePointsDist2(tuple(grdPos[ptInd]),diag,res,dist2)
                    d = min(dist2[0:n])
                    sptInd = res[tuple(dist2).index(d)]
                except :
                    #this is quite long
                    delta = numpy.array(srfPts)-numpy.array(grdPos[ptInd])
                    delta *= delta
                    distA = numpy.sqrt( delta.sum(1) )
                    d = min(distA)
                    sptInd = list(distA).index(d)
                sx, sy, sz = srfPts[sptInd]
#            target = afvi.vi.getObject("test")
#            if target is None :
#                target = afvi.vi.Sphere("test",radius=10.0,color=[0,0,1])[0]
#            afvi.vi.setTranslation(target,pos=srfPts[sptInd])
#            target2 = afvi.vi.getObject("test2")
#            if target2 is None :
#                target2 = afvi.vi.Sphere("test2",radius=10.0,color=[0,0,1])[0]
#            afvi.vi.changeObjColorMat(target2,[0,0,1])
#            afvi.vi.setTranslation(target2,pos=grdPos[ptInd])
#            afvi.vi.update()
            # update distance field
            # we should not reompute this ...
            #if ptInd < len(distances)-1:  # Oct. 20, 2012 Graham turned this if off because this dist override is necessary in
            if distances[ptInd]>d: distances[ptInd] = d  # case a diffent surface ends up being closer in the linear walk through the grid

            # check if ptInd in inside
            nx, ny, nz = numpy.array(ogNormals[sptInd])
#            target3 = afvi.vi.getObject("test3")
#            if target3 is None :
#                target3 = afvi.vi.Sphere("test3",radius=10.0,color=[0,0,1])[0]
#            afvi.vi.changeObjColorMat(target3,[0,0,1])
#            afvi.vi.setTranslation(target3,pos=numpy.array(srfPts[sptInd])+numpy.array(ogNormals[sptInd])*10.)
#            afvi.vi.update()
            
            # check on what side of the surface point the grid point is
            vx,vy,vz = (gx-sx, gy-sy, gz-sz)
            dot = vx*nx + vy*ny + vz*nz
            if dot <= 0: # inside
                #and the point is actually inside the mesh bounding box
                inside = True
                if self.checkinside :
                    inside  = self.checkPointInsideBB(grdPos[ptInd],dist=d)
                #this is not working for a plane, or any unclosed compartment...
                if inside :
                    if ptInd < len(idarray)-1:   #Oct 20, 2012 Graham asks: why do we do this if test? not in old code
                        idarray[ptInd] = -number
                    insidePoints.append(ptInd)
#                if target2 is not None :
#                    afvi.vi.changeObjColorMat(target2,[1,0,0])
            #sleep(0.01)
#            c4d.StatusSetBar(int((ptInd/len(grdPos)*100)))
        print('time to update distance field and idarray', time()-t1)
        
        t1 = time()
        nbGridPoints = len(histoVol.grid.masterGridPositions)
        
        surfPtsBB,surfPtsBBNorms = self.getSurfaceBB(srfPts,histoVol)
        srfPts = surfPtsBB
        print ("compare length id distances",nbGridPoints,(nbGridPoints == len(idarray)),(nbGridPoints == len(distances)))
        ex = True#True if nbGridPoints == len(idarray) else False
        surfacePoints,surfacePointsNormals  = self.extendGridArrays(nbGridPoints,
                                            srfPts,surfPtsBBNorms,histoVol,extended=ex)

        insidePoints = insidePoints
        print('time to extend arrays', time()-t1)

        print('Total time', time()-t0)

        self.insidePoints = insidePoints
        self.surfacePoints = surfacePoints
        self.surfacePointsCoords = surfPtsBB
        self.surfacePointsNormals = surfacePointsNormals
        print('%s surface pts, %d inside pts, %d tot grid pts, %d master grid'%(
            len(self.surfacePoints), len(self.insidePoints),
            nbGridPoints, len(histoVol.grid.masterGridPositions)))

        self.computeVolumeAndSetNbMol(histoVol, self.surfacePoints,
                                      self.insidePoints,areas=vSurfaceArea)
#        bhtreelib.freeBHtree(bht)
        return self.insidePoints, self.surfacePoints

    def extendGridArrays(self,nbGridPoints,srfPts,surfPtsBBNorms,histoVol,extended=True):
        """Extend the environment grd using the compartment point"""
        if extended  :      
            length = len(srfPts)
            pointArrayRaw = numpy.zeros( (nbGridPoints + length, 3), 'f')
            pointArrayRaw[:nbGridPoints] = histoVol.grid.masterGridPositions
            pointArrayRaw[nbGridPoints:] = srfPts
            self.surfacePointscoords = srfPts
            histoVol.grid.nbSurfacePoints += length
            histoVol.grid.masterGridPositions = pointArrayRaw
            if type(histoVol.grid.distToClosestSurf) == numpy.ndarray :
                #histoVol.grid.distToClosestSurf = numpy.append(histoVol.grid.distToClosestSurf,numpy.array([histoVol.grid.diag,]*length ))
                distCS = numpy.ones(length)*histoVol.grid.diag                
                histoVol.grid.distToClosestSurf = numpy.hstack((histoVol.grid.distToClosestSurf,distCS))
            else :
                histoVol.grid.distToClosestSurf.extend( [histoVol.grid.diag,]*length )
            ptId = numpy.ones(length,'i')*self.number#surface point
            histoVol.grid.gridPtId=numpy.hstack((histoVol.grid.gridPtId,ptId))
            #histoVol.grid.gridPtId=numpy.append(numpy.array(histoVol.grid.gridPtId), [self.number]*length ,axis=0)#surface point ID
            histoVol.grid.freePoints = numpy.arange(nbGridPoints+length)
            surfacePoints = list(range(nbGridPoints, nbGridPoints+length))
            #histoVol.grid.freePoints.extend(surfacePoints)
    
            surfacePointsNormals = {}
            for i, n in enumerate(surfPtsBBNorms):
                surfacePointsNormals[nbGridPoints + i] = n
        else :
            length = len(srfPts)
            pointArrayRaw = histoVol.grid.masterGridPositions
            self.surfacePointscoords = srfPts
            histoVol.grid.nbSurfacePoints += length
            surfacePoints = list(range(nbGridPoints-length, nbGridPoints))
            surfacePointsNormals = {}
            for i, n in enumerate(surfPtsBBNorms):
                surfacePointsNormals[nbGridPoints-length + i] = n            
        return surfacePoints,surfacePointsNormals
    
    def getSurfaceBB(self,srfPts,histoVol):
        """get the bounding box from the environment grid that encapsulated the mesh"""
        surfPtsBB = []
        surfPtsBBNorms  = []
        mini, maxi = histoVol.fillBB
        mx, my, mz = mini
        Mx, My, Mz = maxi
        ogNorms = self.ogsurfacePointsNormals
        for i,p in enumerate(srfPts):
            x,y,z = p
            if (x>=mx and x<=Mx and y>=my and y<=My and z>=mz and z<=Mz):
                surfPtsBB.append(p)
                surfPtsBBNorms.append(ogNorms[i])        
        if self.highresVertices is not None :
            srfPts = self.highresVertices   
            surfPtsBB = []     
            for i,p in enumerate(srfPts):
                x,y,z = p
                if (x>=mx and x<=Mx and y>=my and y<=My and z>=mz and z<=Mz):
                    surfPtsBB.append(p)
        print('surf points going from to', len(srfPts), len(surfPtsBB))
        srfPts = surfPtsBB
        return surfPtsBB,surfPtsBBNorms

#    def BuildGrid_break(self, histoVol):
#        # create surface points
#        t0 = t1 = time()
#        self.createSurfacePoints(maxl=histoVol.grid.gridSpacing)
#        print("Creating surface points and preparing to sum the surface area, TODO- limit to selection box")
#        # Graham Sum the SurfaceArea for each polyhedron
#        vertices = self.vertices  #NEED to make these limited to selection box, not whole compartment
#        faces = self.faces #         Should be able to use self.ogsurfacePoints and collect faces too from above
#        normalList2,areas = self.getFaceNormals(vertices, faces,fillBB=histoVol.fillBB)
#        vSurfaceArea = sum(areas)
#        #for gnum in range(len(normalList2)):
#        #    vSurfaceArea = vSurfaceArea + areas[gnum]
#        print('Graham says Surface Area is %d' %(vSurfaceArea))
#        print('Graham says the last triangle area is is %d' %(areas[-1]))
#        #print '%d surface points %.2f unitVol'%(len(surfacePoints), unitVol)
#        
#        # build a BHTree for the vertices
#        
#        print('time to create surface points', time()-t1, len(self.ogsurfacePoints))
#        
#        distances = histoVol.grid.distToClosestSurf
#        idarray = histoVol.grid.gridPtId
#        diag = histoVol.grid.diag
#
#        t1 = time()
#
#        #build BHTree for off grid surface points
#        from bhtree import bhtreelib
#        srfPts = self.ogsurfacePoints
#        bht = self.OGsrfPtsBht = bhtreelib.BHtree( srfPts, None, 10)
#
#        number = self.number
#        ogNormals = self.ogsurfacePointsNormals
#        insidePoints = []
#        #self.vnormals = ogNormals = normalList2
#        # find closest off grid surface point for each grid point 
#        #FIXME sould be diag of compartment BB inside fillBB
#        grdPos = histoVol.grid.masterGridPositions
#        returnNullIfFail = 0
#        closest = bht.closestPointsArray(grdPos, diag, returnNullIfFail)
##        def distanceLoop(ptInd,distances,grdPos,closest,srfPts,ogNormals,idarray,insidePoints,number):
##            # find closest OGsurfacepoint
##            gx, gy, gz = grdPos[ptInd]
##            sptInd = closest[ptInd]
##            if closest[ptInd]==-1:
##                pdb.set_trace()
##            sx, sy, sz = srfPts[sptInd]
##
##            # update distance field
##            d = sqrt( (gx-sx)*(gx-sx) + (gy-sy)*(gy-sy) +
##                      (gz-sz)*(gz-sz))
##            if distances[ptInd]>d: distances[ptInd] = d
##
##            # check if ptInd in inside
##            nx, ny, nz = ogNormals[sptInd]
##            # check on what side of the surface point the grid point is
##            vx,vy,vz = (gx-sx, gy-sy, gz-sz)
##            dot = vx*nx + vy*ny + vz*nz
##            if dot < 0: # inside
##                idarray[ptInd] = -number
##                insidePoints.append(ptInd)
#            
#        
#        #[distanceLoop(x,distances,grdPos,closest,srfPts,ogNormals,idarray,insidePoints,number) for x in xrange(len(grdPos))]
#        for ptInd in range(len(grdPos)):
#
#            # find closest OGsurfacepoint
#            gx, gy, gz = grdPos[ptInd]
#            sptInd = closest[ptInd]
#            if closest[ptInd]==-1:
#                pdb.set_trace()
#            sx, sy, sz = srfPts[sptInd]
#
#            # update distance field
#            #measure distance between the grid Point and the surface point
#            d = sqrt( (gx-sx)*(gx-sx) + (gy-sy)*(gy-sy) +
#                      (gz-sz)*(gz-sz))
#            if distances[ptInd] > d : 
#                distances[ptInd] = d
#
#            # check if ptInd in inside, and look at the normal at this points
#            nx, ny, nz = ogNormals[sptInd]
#            # check on what side of the surface point the grid point is
#            vx,vy,vz = (gx-sx, gy-sy, gz-sz)
#            dot = vx*nx + vy*ny + vz*nz
#            if dot < 0: # inside
#                idarray[ptInd] = -number
#                insidePoints.append(ptInd)
#
#        print('time to update distance field and idarray', time()-t1)
#
#        t1 = time()
#        nbGridPoints = len(histoVol.grid.masterGridPositions)
#
#        surfPtsBB = []
#        surfPtsBBNorms  = []
#        mini, maxi = histoVol.fillBB
#        mx, my, mz = mini
#        Mx, My, Mz = maxi
#        ogNorms = self.ogsurfacePointsNormals
#        for i,p in enumerate(srfPts):
#            x,y,z = p
#            if (x>=mx and x<=Mx and y>=my and y<=My and z>=mz and z<=Mz):
#                surfPtsBB.append(p)
#                surfPtsBBNorms.append(ogNorms[i])
#
#        print('surf points going from to', len(srfPts), len(surfPtsBB))
#        srfPts = surfPtsBB
#        length = len(srfPts)
#
#        pointArrayRaw = numpy.zeros( (nbGridPoints + length, 3), 'f')
#        pointArrayRaw[:nbGridPoints] = histoVol.grid.masterGridPositions
#        pointArrayRaw[nbGridPoints:] = srfPts
#        self.surfacePointsCoords = srfPts #surfacePointscoords ?
#        histoVol.grid.nbSurfacePoints += length
#        histoVol.grid.masterGridPositions = pointArrayRaw
#        histoVol.grid.distToClosestSurf.extend( [histoVol.grid.diag]*length )
#        
#        histoVol.grid.gridPtId.extend( [number]*length )
#        surfacePoints = list(range(nbGridPoints, nbGridPoints+length))
#        histoVol.grid.freePoints.extend(surfacePoints)
#
#        surfacePointsNormals = {}
#        for i, n in enumerate(surfPtsBBNorms):
#            surfacePointsNormals[nbGridPoints + i] = n
#
#        insidePoints = insidePoints
#        print('time to extend arrays', time()-t1)
#
#        print('Total time', time()-t0)
#
#        self.insidePoints = insidePoints
#        self.surfacePoints = surfacePoints
#        self.surfacePointsCoords = surfPtsBB
#        self.surfacePointsNormals = surfacePointsNormals
#        print('%s surface pts, %d inside pts, %d tot grid pts, %d master grid'%(
#            len(self.surfacePoints), len(self.insidePoints),
#            nbGridPoints, len(histoVol.grid.masterGridPositions)))
#
#        self.computeVolumeAndSetNbMol(histoVol, self.surfacePoints,
#                                      self.insidePoints,areas=vSurfaceArea)
#        return self.insidePoints, self.surfacePoints

    
    def BuildGridEnviroOnly(self, histoVol, location=None):
        """Build the compartment grid ie surface and inside only environment"""
        # create surface points
        t0 = t1 = time()
        self.createSurfacePoints(maxl=histoVol.grid.gridSpacing)
        
        # Graham Sum the SurfaceArea for each polyhedron
        vertices = self.vertices  #NEED to make these limited to selection box, not whole compartment
        faces = self.faces #         Should be able to use self.ogsurfacePoints and collect faces too from above
        normalList2,areas = self.getFaceNormals(vertices, faces,fillBB=histoVol.fillBB)
        vSurfaceArea = sum(areas)
        #for gnum in range(len(normalList2)):
        #    vSurfaceArea = vSurfaceArea + areas[gnum]
        #print 'Graham says Surface Area is %d' %(vSurfaceArea)
        #print 'Graham says the last triangle area is is %d' %(areas[-1])
        #print '%d surface points %.2f unitVol'%(len(surfacePoints), unitVol)
        
        # build a BHTree for the vertices
        
        print('time to create surface points', time()-t1, len(self.ogsurfacePoints))

        
        distances = histoVol.grid.distToClosestSurf
        idarray = histoVol.grid.gridPtId
#        diag = histoVol.grid.diag

        t1 = time()

        #build BHTree for off grid surface points
#        from bhtree import bhtreelib
        srfPts = self.ogsurfacePoints
#        bht = self.OGsrfPtsBht = bhtreelib.BHtree( srfPts, None, 10)

        number = self.number
        ogNormals = self.ogsurfacePointsNormals
        insidePoints = []

        # find closest off grid surface point for each grid point 
        #FIXME sould be diag of compartment BB inside fillBB
        grdPos = histoVol.grid.masterGridPositions
#        returnNullIfFail = 0
        closest = []#bht.closestPointsArray(grdPos, diag, returnNullIfFail)
        def distanceLoop(ptInd,distances,grdPos,closest,srfPts,ogNormals,idarray,insidePoints,number):
            # find closest OGsurfacepoint
            gx, gy, gz = grdPos[ptInd]
            sptInd = closest[ptInd]
            if closest[ptInd]==-1:
                pdb.set_trace()
            sx, sy, sz = srfPts[sptInd]

            # update distance field
            d = sqrt( (gx-sx)*(gx-sx) + (gy-sy)*(gy-sy) +
                      (gz-sz)*(gz-sz))
            if distances[ptInd]>d: distances[ptInd] = d

            # check if ptInd in inside
            nx, ny, nz = ogNormals[sptInd]
            # check on what side of the surface point the grid point is
            vx,vy,vz = (gx-sx, gy-sy, gz-sz)
            dot = vx*nx + vy*ny + vz*nz
            if dot < 0: # inside
                idarray[ptInd] = -number
                insidePoints.append(ptInd)
            
        if location is None :
            re=[distanceLoop(x, distances, grdPos, closest, srfPts,
                              ogNormals, idarray, insidePoints,
                               number) for x in range(len(grdPos))]
        else :
            insidePoints = list(range(len(grdPos)))
            for ptInd in range(len(grdPos)):
                distances[ptInd] = 99999.
                idarray[ptInd] = location

#        for ptInd in xrange(len(grdPos)):

            # find closest OGsurfacepoint
#            gx, gy, gz = grdPos[ptInd]
#            sptInd = closest[ptInd]
#            if closest[ptInd]==-1:
#                pdb.set_trace()
#            sx, sy, sz = srfPts[sptInd]

            # update distance field
#            d = sqrt( (gx-sx)*(gx-sx) + (gy-sy)*(gy-sy) +
#                      (gz-sz)*(gz-sz))
#            if distances[ptInd]>d: distances[ptInd] = d

            # check if ptInd in inside
#            nx, ny, nz = ogNormals[sptInd]
            # check on what side of the surface point the grid point is
#            vx,vy,vz = (gx-sx, gy-sy, gz-sz)
#            dot = vx*nx + vy*ny + vz*nz
#            if dot < 0: # inside
#                idarray[ptInd] = -number
#                insidePoints.append(ptInd)

        print('time to update distance field and idarray', time()-t1)

        t1 = time()
        nbGridPoints = len(histoVol.grid.masterGridPositions)

        surfPtsBB = []
        surfPtsBBNorms  = []
        mini, maxi = histoVol.fillBB
        mx, my, mz = mini
        Mx, My, Mz = maxi
        ogNorms = self.ogsurfacePointsNormals
        for i,p in enumerate(srfPts):
            x,y,z = p
            if (x>=mx and x<=Mx and y>=my and y<=My and z>=mz and z<=Mz):
                surfPtsBB.append(p)
                surfPtsBBNorms.append(ogNorms[i])

        print('surf points going from to', len(srfPts), len(surfPtsBB))
        srfPts = surfPtsBB
        length = len(srfPts)

        pointArrayRaw = numpy.zeros( (nbGridPoints + length, 3), 'f')
        pointArrayRaw[:nbGridPoints] = histoVol.grid.masterGridPositions
        pointArrayRaw[nbGridPoints:] = srfPts
        self.surfacePointsCoords = srfPts #surfacePointscoords ?
        histoVol.grid.nbSurfacePoints += length
        histoVol.grid.masterGridPositions = pointArrayRaw
        histoVol.grid.distToClosestSurf.extend( [histoVol.grid.diag]*length )
        
        histoVol.grid.gridPtId.extend( [number]*length )
        surfacePoints = list(range(nbGridPoints, nbGridPoints+length))
        histoVol.grid.freePoints.extend(surfacePoints)

        surfacePointsNormals = {}
        for i, n in enumerate(surfPtsBBNorms):
            surfacePointsNormals[nbGridPoints + i] = n

        insidePoints = insidePoints
        print('time to extend arrays', time()-t1)

        print('Total time', time()-t0)

        self.insidePoints = insidePoints
        self.surfacePoints = surfacePoints
        self.surfacePointsCoords = surfPtsBB
        self.surfacePointsNormals = surfacePointsNormals
        print('%s surface pts, %d inside pts, %d tot grid pts, %d master grid'%(
            len(self.surfacePoints), len(self.insidePoints),
            nbGridPoints, len(histoVol.grid.masterGridPositions)))

        self.computeVolumeAndSetNbMol(histoVol, self.surfacePoints,
                                      self.insidePoints,areas=vSurfaceArea)
        return self.insidePoints, self.surfacePoints

    def BuildGrid_OrthogonalBox(self, histoVol):
        """Build the compartment grid ie surface and inside point using a box"""
        t0 = time()
        self.ogsurfacePoints = None
        self.ogsurfacePointsNormals = None
#        vertices = None
#        faces = None
        normalList2,areas = None # in future, get area of box from corner points
        vSurfaceArea = None#sum(areas)
#        bbox = self.bb
#        xmin = bbox[0][0]; ymin = bbox[0][1]; zmin = bbox[0][2]
#        xmax = bbox[1][0]; ymax = bbox[1][1]; zmax = bbox[1][2]
#        sizex = self.getSizeXYZ()
#        gboundingBox = histoVol.grid.boundingBox
#        gspacing = histoVol.grid.gridSpacing
#        
#        from UTpackages.UTsdf import utsdf
#        #can be 16,32,64,128,256,512,1024
#        #        if spacing not in [16,32,64,128,256,512,1024]:
#        #            spacing = self.find_nearest(numpy.array([16,32,64,128,256,512,1024]),spacing)
#        # compute SDF
#        dim=16
#        dim1=dim+1
#        print ("ok2 dim ",dim)
#        size = dim1*dim1*dim1
#        from UTpackages.UTsdf import utsdf
#        verts = N.array(self.vertices,dtype='f')
#        
#        tris = N.array(self.faces,dtype="int")
#        utsdf.setParameters(dim,0,1,[0,0,0,0,0,0])#size, bool isNormalFlip, bool insideZero,bufferArr
#        surfacePoints = srfPts = self.vertices
#        print ("ok grid points")
#        datap = utsdf.computeSDF(N.ascontiguousarray(verts, dtype=N.float32),N.ascontiguousarray(tris, dtype=N.int32))
#        print ("ok computeSDF")
#        data = utsdf.createNumArr(datap,size)
#        volarr = data[:]
#        volarr.shape = (dim1, dim1, dim1)
#        volarr = numpy.ascontiguousarray(numpy.transpose(volarr), 'f')
#        
#        # get grid points distances to compartment surface
#        from Volume.Operators.trilinterp import trilinterp
#        invstep =(1./(sizex[0]/dim), 1./(sizex[1]/dim), 1./(sizex[2]/dim))
#        origin = self.bb[0]
#        distFromSurf = trilinterp(histoVol.grid.masterGridPositions,
#                                  volarr, invstep, origin)
#        
#        
#        # save SDF
#        #        self.sdfData = volarr
#        #        self.sdfOrigin = origin
#        #        self.sdfGridSpacing = (gSizeX, gSizeY, gSizeZ)
#        #        self.sdfDims = (dimx, dimy, dimz)
#        
#        ## update histoVol.distToClosestSurf
#        distance = histoVol.grid.distToClosestSurf
#        for i,d in enumerate(distFromSurf):
#            if distance[i] > d:
#                distance[i] = d
        
        # loop over fill box grid points and build the idarray
        # identify inside and surface points and update the distance field
        number = self.number
        insidePoints = []
        surfacePoints = []
        allNormals = {}
        idarray = histoVol.grid.gridPtId
        #surfaceCutOff = histoVol.gridSpacing*.5
        #print 'BBBBBBBBBBBBBB', surfaceCutOff, min(distFromSurf), max(distFromSurf)
        #print 'We should get', len(filter(lambda x:fabs(x)<surfaceCutOff, distance))
    
        #import pdb
        #pdb.set_trace()
        indice = numpy.nonzero(numpy.less(distance,0.0))
        pointinside = numpy.take(histoVol.grid.masterGridPositions,indice,0)[0]
        #        print (len(indice[0]),indice,len(pointinside))
        if len(indice) == 1  and len(indice[0]) != 1:
            indice = indice[0]
        if len(pointinside) == 1  and len(pointinside[0]) != 1:
            pointinside = pointinside[0]
        histoVol.grid.gridPtId[indice] = -self.number
        print ("sdf pointID N ",self.number,len(histoVol.grid.gridPtId[indice]),histoVol.grid.gridPtId[indice])
        t1 = time()
        nbGridPoints = len(histoVol.grid.masterGridPositions)

        surfPtsBB,surfPtsBBNorms = self.getSurfaceBB(srfPts,histoVol)
        srfPts = surfPtsBB
        surfacePoints,surfacePointsNormals  = self.extendGridArrays(nbGridPoints,
                                                                    srfPts,surfPtsBBNorms,histoVol)
        
        print (len(histoVol.grid.gridPtId[indice]),histoVol.grid.gridPtId[indice])
        insidePoints = pointinside
        print('time to extend arrays', time()-t1)
        
        print('Total time', time()-t0)
        
        self.insidePoints = insidePoints
        self.surfacePoints = surfacePoints
        self.surfacePointsCoords = surfPtsBB
        self.surfacePointsNormals = surfacePointsNormals
        print('%s surface pts, %d inside pts, %d tot grid pts, %d master grid'%(
                                                                                len(self.surfacePoints), len(self.insidePoints),
                                                                                nbGridPoints, len(histoVol.grid.masterGridPositions)))
        self.computeVolumeAndSetNbMol(histoVol, self.surfacePoints,
                                      self.insidePoints,areas=vSurfaceArea)    
        return insidePoints, surfacePoints



    def BuildGrid_utsdf(self, histoVol):        
        """
        Build the compartment grid ie surface and inside point using signed distance fields
        from the UT package
        """
        t0 = time()
        self.ogsurfacePoints = self.vertices[:]
        self.ogsurfacePointsNormals = self.vnormals[:]
        vertices = self.vertices
        faces = self.faces
        normalList2,areas = self.getFaceNormals(vertices, faces,fillBB=histoVol.fillBB)
        vSurfaceArea = sum(areas)
#        labels = numpy.ones(len(faces), 'i')
        
        # FIXME .. dimensions on SDF should addapt to compartment size
        bbox = self.bb
        xmin = bbox[0][0]; ymin = bbox[0][1]; zmin = bbox[0][2]
        xmax = bbox[1][0]; ymax = bbox[1][1]; zmax = bbox[1][2]
        sizex = self.getSizeXYZ()
        gboundingBox = histoVol.grid.boundingBox
        gspacing = histoVol.grid.gridSpacing
        
        from UTpackages.UTsdf import utsdf
        #can be 16,32,64,128,256,512,1024
#        if spacing not in [16,32,64,128,256,512,1024]:
#            spacing = self.find_nearest(numpy.array([16,32,64,128,256,512,1024]),spacing)
        # compute SDF
        dim=16
        dim1=dim+1        
        print ("ok2 dim ",dim)
        size = dim1*dim1*dim1
        from UTpackages.UTsdf import utsdf        
        verts = N.array(self.vertices,dtype='f')
        
        tris = N.array(self.faces,dtype="int")
        utsdf.setParameters(dim,0,1,[0,0,0,0,0,0])#size, bool isNormalFlip, bool insideZero,bufferArr
        surfacePoints = srfPts = self.vertices
        print ("ok grid points")
        datap = utsdf.computeSDF(N.ascontiguousarray(verts, dtype=N.float32),N.ascontiguousarray(tris, dtype=N.int32))
        print ("ok computeSDF")
        data = utsdf.createNumArr(datap,size)
        volarr = data[:]
        volarr.shape = (dim1, dim1, dim1)
        volarr = numpy.ascontiguousarray(numpy.transpose(volarr), 'f')
        
        # get grid points distances to compartment surface
        from Volume.Operators.trilinterp import trilinterp
        invstep =(1./(sizex[0]/dim), 1./(sizex[1]/dim), 1./(sizex[2]/dim))
        origin = self.bb[0]
        distFromSurf = trilinterp(histoVol.grid.masterGridPositions,
                                  volarr, invstep, origin)

        
        # save SDF
#        self.sdfData = volarr
#        self.sdfOrigin = origin
#        self.sdfGridSpacing = (gSizeX, gSizeY, gSizeZ)
#        self.sdfDims = (dimx, dimy, dimz)

        ## update histoVol.distToClosestSurf
        distance = histoVol.grid.distToClosestSurf
        for i,d in enumerate(distFromSurf):
            if distance[i] > d:
                distance[i] = d
        
        # loop over fill box grid points and build the idarray
        # identify inside and surface points and update the distance field
        number = self.number
        insidePoints = []
        surfacePoints = []
        allNormals = {}
        idarray = histoVol.grid.gridPtId
        #surfaceCutOff = histoVol.gridSpacing*.5
        #print 'BBBBBBBBBBBBBB', surfaceCutOff, min(distFromSurf), max(distFromSurf)
        #print 'We should get', len(filter(lambda x:fabs(x)<surfaceCutOff, distance))

        #import pdb
        #pdb.set_trace()
        indice = numpy.nonzero(numpy.less(distance,0.0))
        pointinside = numpy.take(histoVol.grid.masterGridPositions,indice,0)[0]
#        print (len(indice[0]),indice,len(pointinside))
        if len(indice) == 1  and len(indice[0]) != 1:
            indice = indice[0]
        if len(pointinside) == 1  and len(pointinside[0]) != 1:
            pointinside = pointinside[0]
        histoVol.grid.gridPtId[indice] = -self.number
        print ("sdf pointID N ",self.number,len(histoVol.grid.gridPtId[indice]),histoVol.grid.gridPtId[indice])
        t1 = time()
        nbGridPoints = len(histoVol.grid.masterGridPositions)
        
        surfPtsBB,surfPtsBBNorms = self.getSurfaceBB(srfPts,histoVol)
        srfPts = surfPtsBB
        surfacePoints,surfacePointsNormals  = self.extendGridArrays(nbGridPoints,
                                            srfPts,surfPtsBBNorms,histoVol)

        print (len(histoVol.grid.gridPtId[indice]),histoVol.grid.gridPtId[indice])                
        insidePoints = pointinside
        print('time to extend arrays', time()-t1)

        print('Total time', time()-t0)

        self.insidePoints = insidePoints
        self.surfacePoints = surfacePoints
        self.surfacePointsCoords = surfPtsBB
        self.surfacePointsNormals = surfacePointsNormals
        print('%s surface pts, %d inside pts, %d tot grid pts, %d master grid'%(
            len(self.surfacePoints), len(self.insidePoints),
            nbGridPoints, len(histoVol.grid.masterGridPositions)))
        self.computeVolumeAndSetNbMol(histoVol, self.surfacePoints,
                                      self.insidePoints,areas=vSurfaceArea)    
        return insidePoints, surfacePoints
        
    def get_bbox(self, vert_list, BB_SCALE = 0.0):
        """get bounding box for the given list of vertices"""
        from multisdf import multisdf
        multisdf.cvar.BB_SCALE = BB_SCALE
        HUGE = 999999

        bbox = []
        x_min = HUGE; x_max = -HUGE
        y_min = HUGE; y_max = -HUGE
        z_min = HUGE; z_max = -HUGE
        for i in range(len(vert_list)):
            p = vert_list[i]
            #check x-span
            if p[0] < x_min: 
                x_min = p[0]
            if p[0] > x_max: 
                x_max = p[0]
            #check y-span
            if p[1] < y_min: 
                y_min = p[1]
            if p[1] > y_max: 
                y_max = p[1]
            #check z-span
            if p[2] < z_min: 
                z_min = p[2]
            if p[2] > z_max: 
                z_max = p[2]

        bbox.append(x_min - BB_SCALE*(x_max-x_min))
        bbox.append(y_min - BB_SCALE*(y_max-y_min))
        bbox.append(z_min - BB_SCALE*(z_max-z_min))

        bbox.append(x_max + BB_SCALE*(x_max-x_min))
        bbox.append(y_max + BB_SCALE*(y_max-y_min))
        bbox.append(z_max + BB_SCALE*(z_max-z_min))
        return bbox

        
    def BuildGrid_multisdf(self, histoVol):
        """Build the compartment grid ie surface and inside point using multisdf"""
        vertices = self.vertices
        faces = self.faces
        labels = numpy.ones(len(faces), 'i')
        
        # FIXME .. dimensions on SDF should addapt to compartment size
        bbox = self.get_bbox(vertices)
        xmin = bbox[0]; ymin = bbox[1]; zmin = bbox[2]
        xmax = bbox[3]; ymax = bbox[4]; zmax = bbox[5]

        # compute SDF
        from multisdf import multisdf
        gridSpacing = 30.
        dimx = int( (xmax-xmin)/gridSpacing ) + 1
        dimy = int( (ymax-ymin)/gridSpacing ) + 1
        dimz = int( (zmax-zmin)/gridSpacing ) + 1

        gSizeX = (xmax-xmin)/(dimx-1) 
        gSizeY = (ymax-ymin)/(dimy-1)
        gSizeZ = (zmax-zmin)/(dimz-1)

        print('SDF grid size', dimx,dimy,dimz,gSizeX, gSizeY, gSizeZ)

        mind = -1000.
        maxd = 1000.
        datap = multisdf.computeSDF(vertices, faces, labels, dimx, dimy, dimz,
                                    maxd, mind)
        grid_size  = dimx*dimy*dimz
        volarr = multisdf.createNumArr(datap, grid_size)
        volarr.shape = (dimz, dimy, dimx)
        volarr = numpy.ascontiguousarray(numpy.transpose(volarr), 'f')

        # get grid points distances to compartment surface
        from Volume.Operators.trilinterp import trilinterp
        invstep =(1./gridSpacing, 1./gridSpacing, 1./gridSpacing)
        origin = (xmin, ymin, zmin)
        distFromSurf = trilinterp(histoVol.masterGridPositions,
                                  volarr, invstep, origin)

        # save SDF
        self.sdfData = volarr
        self.sdfOrigin = origin
        self.sdfGridSpacing = (gSizeX, gSizeY, gSizeZ)
        self.sdfDims = (dimx, dimy, dimz)

        ## update histoVol.distToClosestSurf
        distance = histoVol.distToClosestSurf
        for i,d in enumerate(distFromSurf):
            if distance[i] > d:
                distance[i] = d
        
        # loop over fill box grid points and build the idarray
        # identify inside and surface points and update the distance field
        number = self.number
        insidePoints = []
        surfacePoints = []
        allNormals = {}
        idarray = histoVol.gridPtId
        #surfaceCutOff = histoVol.gridSpacing*.5
        #print 'BBBBBBBBBBBBBB', surfaceCutOff, min(distFromSurf), max(distFromSurf)
        #print 'We should get', len(filter(lambda x:fabs(x)<surfaceCutOff, distance))

        #import pdb
        #pdb.set_trace()
        
        for i, d in enumerate(distance):

            # identify surface and interior points
            # there is a problem with SDF putting large negative values
            # for inside points. For now we pick all negative != mind as
            # surface points
            if d>0:
                continue
            elif d<mind:
                surfacePoints.append(i)
                idarray[i] = number
                allNormals[i] = (1,0,0)
            else:
                insidePoints.append(i)
                idarray[i] = -number

        self.computeVolumeAndSetNbMol(histoVol, surfacePoints, insidePoints)

        self.insidePoints = insidePoints
        self.surfacePoints = surfacePoints
        self.surfacePointsNormals = allNormals
        print('AAAAAAAAAAAA', len(surfacePoints))
        
        return insidePoints, surfacePoints


    def getSurfacePoint( self, p1, p2, w1, w2):
        """compute point between p1 and p2 with weight w1 and w2"""
        x1,y1,z1 = p1
        x2,y2,z2 = p2
#        totalWeight = w1+w2
        ratio = w1/(w1+w2)
        vec = (x2-x1, y2-y1, z2-z1)
        return x1 + ratio*vec[0], y1 + ratio*vec[1], z1 + ratio*vec[2]

    def estimateVolume(self,unitVol=None,hBB=None):
        """ 
        get the volume of the compartment
        v: A pointer to the array of vertices
        // i: A pointer to the array of indices
        // n: The number of indices (multiple of 3)
        // This function uses Gauss's Theorem to calculate the volume of a body
        // enclosed by a set of triangles. The triangle set must form a closed
        // surface in R3 and be outward facing. Outward facing triangles index
        // their vertices in a counterclockwise order where the x-axis points
        // left, the y-axis point up and the z-axis points toward you (rhs).
        // from http://www.gamedev.net/page/resources/_/technical/game-programming/area-and-volume-calculations-r2247
        """
        if self.interiorVolume is None or self.interiorVolume == 0.0:
            v = self.vertices
            i = self.faces
            n = len(self.faces)
            volume = 0.0
            for j in range(n):#(j = 0; j < n; j+=3)
                v1 = v[i[j][0]]
                v2 = v[i[j][1]]
                v3 = v[i[j][2]]    
                volume += ((v2[1]-v1[1])*(v3[2]-v1[2])-(v2[2]-v1[2])*(v3[1]-v1[1]))*(v1[0]+v2[0]+v3[0])    
            self.interiorVolume =  volume / 6.0
        if self.surfaceVolume is None or self.surfaceVolume == 0.0:
            if hBB is not None :
                normalList2,areas = self.getFaceNormals(self.vertices,self.faces,fillBB=hBB)
                self.surfaceVolume =  sum(areas)            
            elif unitVol != None :
                self.surfaceVolume = len(self.vertices) * unitVol

    def computeVolumeAndSetNbMol(self, histoVol, surfacePoints, insidePoints,
                                 areas=None):
        """
        Compute volume of surface and interior
        set 'nbMol' in each ingredient of both recipes
        """        
        unitVol = histoVol.grid.gridSpacing**3
        if surfacePoints !=None:
            print('%d surface points %.2f unitVol'%(len(surfacePoints), unitVol))
        #FIXME .. should be surface per surface point instead of unitVol
            self.surfaceVolume = len(surfacePoints)*unitVol 
        area = False
        if areas is not None :
            self.surfaceVolume = areas
            area = True
        self.interiorVolume = len(insidePoints)*unitVol
        if self.surfaceVolume != None:
            print('%d surface volume %.2f interior volume'%(self.surfaceVolume, self.interiorVolume))
        print('%.2f interior volume'%(self.interiorVolume))


        # compute number of molecules and save in recipes
        rs = self.surfaceRecipe
        if rs:
            volume = self.surfaceVolume
            rs.setCount(volume,area=area)

        ri = self.innerRecipe
        if ri:
            volume = self.interiorVolume
            a = ri.setCount(volume)
            print ("number of molecules for Special Cube = ", a, ", because interiorVolume = ", volume)

    def getFacesNfromV(self,vindice,ext=0):
        """
        Retrieve the face normal from the indice of a vertice
        """                
        f=[]        
        for i,af in enumerate(self.faces) :
            if vindice in af :
                if ext :
                    for vi in af :
                        if vi != vindice :
                            ff=self.getFacesNfromV(vi)
                            f.extend(ff)
                else :
                    f.append(self.fnormals[i])
        return f

    def getVNfromF(self,i):
        """
        Retrieve the vertice normal from the indice of a vertice
        """                
        self.normals
        fi=[]        
        for k,af in enumerate(self.faces) :
            if i in af :
                for j in af :
                    if j not in fi :
                        fi.append(j)
        n=[]
        for ind in fi:
            n.append(self.normals[ind])
        return n

    def create3DPointLookup(self, nbGridPoints,gridSpacing,dim,boundingBox=None):
        """ 
        Fill the orthogonal bounding box described by two global corners
         with an array of points spaces pGridSpacing apart. Duplicate from grid class
        """
        if boundingBox is None :
            boundingBox= self.bb
        xl,yl,zl = boundingBox[0]
        xr,yr,zr = boundingBox[1]

        nx,ny,nz = nbGridPoints
        pointArrayRaw = numpy.zeros( (nx*ny*nz, 3), 'f')
        ijkPtIndice = numpy.zeros( (nx*ny*nz, 3), 'i')
        space = gridSpacing
        size = self.getSizeXYZ()
        # Vector for lower left broken into real of only the z coord.
        i = 0
        for zi in range(nz):
            for yi in range(ny):
                for xi in range(nx):
                    pointArrayRaw[i] = (xl+xi*(size[0]/dim), yl+yi*(size[1]/dim), zl+zi*(size[2]/dim))
                    ijkPtIndice[i] = (xi,yi,zi)
                    i+=1
        return ijkPtIndice,pointArrayRaw

    def find_nearest(self,array,value):
        """find nearest point indice of value in array using numpy"""
        idx=(numpy.abs(array-value)).argmin()
        return array[idx]


    def getSurfaceInnerPoints_sdf(self,boundingBox,spacing,display = True,useFix=False):
        """
        Only compute the inner point. No grid.
        This is independant from the packing. Help build ingredient sphere tree and representation 
        """                
        print ("beforea import" )        
        print ("ok1" )
        from UTpackages.UTsdf import utsdf
        #can be 16,32,64,128,256,512,1024
        if spacing not in [16,32,64,128,256,512,1024]:
            spacing = self.find_nearest(numpy.array([16,32,64,128,256,512,1024]),spacing)
        dim = spacing
        dim1=dim+1
        
        print ("ok2 ",dim1,dim)
        size = dim1*dim1*dim1
         #can be 16,32,64,128,256,512,1024
        verts = N.array(self.vertices,dtype='f')
        tris = N.array(self.faces,dtype="int")
        utsdf.setParameters(int(dim),0,1,[0,0,0,0,0,0])#size, bool isNormalFlip, bool insideZero,bufferArr
        print ("ok3")

        #spacing = length / 64
        sizes=self.getSizeXYZ()
        L = max(sizes)        
        spacing = L/dim# = self.smallestProteinSize*1.1547  # 2/sqrt(3)????
#        helper.progressBar(label="BuildGRid")        
#        grid.gridVolume,grid.nbGridPoints = grid.computeGridNumberOfPoint(boundingBox,spacing)
        xl,yl,zl = boundingBox[0]
        xr,yr,zr = boundingBox[1]
        
        print (len(verts),len(tris),type(verts),type(tris),verts[0],tris[0])
        surfacePoints = srfPts = self.vertices
        print ("ok grid points")
#        verts = N.ascontiguousarray(verts, dtype='f')
#        tris = N.ascontiguousarray(tris, dtype=N.int32)
#        verts = N.ascontiguousarray(verts, dtype='f')
#        tris = N.ascontiguousarray(tris, dtype=N.int32)

        datap = utsdf.computeSDF(verts,tris)
        #datap = utsdf.computeSDF(verts,tris)
        print ("ok computeSDF")
        data = utsdf.createNumArr(datap,size)

        nbGridPoints = [dim1,dim1,dim1]
        ijkPtIndice,pointArrayRaw = self.create3DPointLookup(nbGridPoints,spacing,dim)
        print ("ok grid",len(data),size)
        nbPoints = len(pointArrayRaw)
        print ("n pts",nbPoints)
        gridPtId = [0]*nbPoints
        grdPos = pointArrayRaw        
        indice = numpy.nonzero(numpy.less(data,0.0))
        pointinside = numpy.take(grdPos,indice,0)
        #need to update the surface. need to create a aligned grid
        return pointinside[0],self.vertices

    def getSurfaceInnerPoints_jordan(self,boundingBox,spacing,display = True,useFix=False,ray=1):
        """
        Only compute the inner point. No grid.
        This is independant from the packing. Help build ingredient sphere tree and representation 
        """                
        from autopack.Environment import Grid        
        grid = Grid()
        grid.boundingBox = boundingBox
        grid.gridSpacing = spacing# = self.smallestProteinSize*1.1547  # 2/sqrt(3)????
        helper.progressBar(label="BuildGRid")        
        grid.gridVolume,grid.nbGridPoints = grid.computeGridNumberOfPoint(boundingBox,spacing)
        grid.create3DPointLookup()
        nbPoints = grid.gridVolume
        grid.gridPtId = [0]*nbPoints
        xl,yl,zl = boundingBox[0]
        xr,yr,zr = boundingBox[1]
        # distToClosestSurf is set to self.diag initially
        grid.diag = diag = vlen( vdiff((xr,yr,zr), (xl,yl,zl) ) )
        grid.distToClosestSurf = [diag]*nbPoints        
        distances = grid.distToClosestSurf
        idarray = grid.gridPtId
        diag = grid.diag
        
        from bhtree import bhtreelib
        self.ogsurfacePoints = self.vertices[:]
        self.ogsurfacePointsNormals = self.vnormals[:]#helper.FixNormals(self.vertices,self.faces,self.vnormals,fn=self.fnormals)
        mat = helper.getTransformation(self.ref_obj)
        surfacePoints = srfPts = self.ogsurfacePoints
        self.OGsrfPtsBht = bht =  bhtreelib.BHtree(tuple(srfPts), None, 10)

        res = numpy.zeros(len(srfPts),'f')
        dist2 = numpy.zeros(len(srfPts),'f')
        number = self.number
        insidePoints = []
        grdPos = grid.masterGridPositions
        returnNullIfFail = 0
        t1=time()
        center = helper.getTranslation( self.ref_obj )
        helper.resetProgressBar()
        if display :
            cyl2 = helper.oneCylinder("ray",[0.,0.,0.],[1.0,1.0,1.0],radius=20.0)
            if ray == 3 :
                cyl1 = helper.oneCylinder("ray2",[0.,0.,0.],[1.0,1.0,1.0],radius=20.0)
                cyl3 = helper.oneCylinder("ray3",[0.,0.,0.],[1.0,1.0,1.0],radius=20.0) 
                helper.changeObjColorMat(cyl1,(1.,1.,1.))
                helper.changeObjColorMat(cyl3,(1.,1.,1.))
        for ptInd in xrange(len(grdPos)):#len(grdPos)):
            inside = False
            t2=time()
            gx, gy, gz = grdPos[ptInd]
            if display :
                helper.updateOneCylinder("ray",grdPos[ptInd],grdPos[ptInd]+(numpy.array([1.1,0.,0.])*spacing*10.0),radius=10.0)
                helper.changeObjColorMat(cyl2,(1.,1.,1.))
                helper.update()
            intersect, count = helper.raycast(self.ref_obj, grdPos[ptInd], grdPos[ptInd]+[1.1,0.,0.], diag, count = True )
            r= ((count % 2) == 1)
            if ray == 3 :
                if display :
                    helper.updateOneCylinder("ray2",grdPos[ptInd],grdPos[ptInd]+(numpy.array([0.,0.0,1.1])*spacing*10.0),radius=10.0)
                    helper.changeObjColorMat(cyl1,(1.,1.,1.))
                    helper.update()
                intersect2, count2 = helper.raycast(self.ref_obj, grdPos[ptInd], grdPos[ptInd]+[0.,0.,1.1], diag, count = True )
                center = helper.rotatePoint(helper.ToVec(center),[0.,0.,0.],[1.0,0.0,0.0,math.radians(33.0)])
                if display :
                    helper.updateOneCylinder("ray3",grdPos[ptInd],grdPos[ptInd]+(numpy.array([0.,1.1,0.])*spacing*10.0),radius=10.0)
                    helper.changeObjColorMat(cyl3,(1.,1.,1.))
                    helper.update()
                intersect3, count3 = helper.raycast(self.ref_obj, grdPos[ptInd], grdPos[ptInd]+[0.,1.1,0.], diag, count = True )
                if r :
                   if (count2 % 2) == 1 and (count3 % 2) == 1 :
                       r=True
                   else : 
                       r=False
            if r : # odd inside
                inside = True
                if inside :
                    if display :
                        helper.changeObjColorMat(cyl2,(1.,0.,0.))
                        if ray == 3 :
                            helper.changeObjColorMat(cyl1,(1.,0.,0.))
                            helper.changeObjColorMat(cyl3,(1.,0.,0.))    
                    #idarray[ptInd] = -number
                    insidePoints.append(grdPos[ptInd]) 
            p=(ptInd/float(len(grdPos)))*100.0
            helper.progressBar(progress=int(p),label=str(ptInd)+"/"+str(len(grdPos))+" inside "+str(inside))
        print('total time', time()-t1)
        return insidePoints, surfacePoints

         
    def getSurfaceInnerPoints_sdf_interpolate(self,boundingBox,spacing,display = True,useFix=False):
        """
        Only compute the inner point. No grid.
        This is independant from the packing. Help build ingredient sphere tree and representation 
        """                
        from autopack.Environment import Grid        
        grid = Grid()
        grid.boundingBox = boundingBox
        grid.gridSpacing = spacing# = self.smallestProteinSize*1.1547  # 2/sqrt(3)????
        helper.progressBar(label="BuildGRid")        
        grid.gridVolume,grid.nbGridPoints = grid.computeGridNumberOfPoint(boundingBox,spacing)
        grid.create3DPointLookup()
        nbPoints = grid.gridVolume
        grid.gridPtId = [0]*nbPoints
        xl,yl,zl = boundingBox[0]
        xr,yr,zr = boundingBox[1]
        # distToClosestSurf is set to self.diag initially
        grid.diag = diag = vlen( vdiff((xr,yr,zr), (xl,yl,zl) ) )
        grid.distToClosestSurf = [diag]*nbPoints        
        distances = grid.distToClosestSurf
        idarray = grid.gridPtId
        diag = grid.diag
        dim=16
        dim1=dim+1        
        print ("ok2")
        size = dim1*dim1*dim1
        from UTpackages.UTsdf import utsdf        
        verts = N.array(self.vertices,dtype='f')
        tris = N.array(self.faces,dtype="int")
        utsdf.setParameters(dim,0,1,[0,0,0,0,0,0])#size, bool isNormalFlip, bool insideZero,bufferArr
        surfacePoints = srfPts = self.vertices
        print ("ok grid points")
        #datap = utsdf.computeSDF(N.ascontiguousarray(verts, dtype=N.float32),N.ascontiguousarray(tris, dtype=N.int32))
        datap = utsdf.computeSDF(verts,tris) #noncontiguous?
        print ("ok computeSDF ",len(verts),len(tris))
        data = utsdf.createNumArr(datap,size)
        volarr = data[:]
        volarr.shape = (dim1, dim1, dim1)
        volarr = numpy.ascontiguousarray(numpy.transpose(volarr), 'f')
        
        # get grid points distances to compartment surface
        from Volume.Operators.trilinterp import trilinterp
        sizex = self.getSizeXYZ()
        invstep =(1./(sizex[0]/dim), 1./(sizex[0]/dim), 1./(sizex[0]/dim))
        origin = self.bb[0]
        distFromSurf = trilinterp(grid.masterGridPositions,
                                  volarr, invstep, origin)
        ## update histoVol.distToClosestSurf
        distance = grid.distToClosestSurf
        for i,d in enumerate(distFromSurf):
            if distance[i] > d:
                distance[i] = d
        
        number = self.number
        insidePoints = []
        surfacePoints = []
        allNormals = {}
#        idarray = histoVol.gridPtId
        indice = numpy.nonzero(numpy.less(distance,0.0))
        pointinside = numpy.take(grid.masterGridPositions,indice,0)
        #need to update the surface. need to create a aligned grid
        return pointinside[0],self.vertices

             
         
    def getSurfaceInnerPoints(self,boundingBox,spacing,display = True,useFix=False):
        """
        Only compute the inner point. No grid.
        This is independant from the packing. Help build ingredient sphere tree and representation 
        """                
        from autopack.Environment import Grid        
        grid = Grid()
        grid.boundingBox = boundingBox
        grid.gridSpacing = spacing# = self.smallestProteinSize*1.1547  # 2/sqrt(3)????
        helper.progressBar(label="BuildGRid")        
        grid.gridVolume,grid.nbGridPoints = grid.computeGridNumberOfPoint(boundingBox,spacing)
        grid.create3DPointLookup()
        nbPoints = grid.gridVolume
        grid.gridPtId = [0]*nbPoints
        xl,yl,zl = boundingBox[0]
        xr,yr,zr = boundingBox[1]
        # distToClosestSurf is set to self.diag initially
        grid.diag = diag = vlen( vdiff((xr,yr,zr), (xl,yl,zl) ) )
        grid.distToClosestSurf = [diag]*nbPoints        
        distances = grid.distToClosestSurf
        idarray = grid.gridPtId
        diag = grid.diag
        
        from bhtree import bhtreelib
        self.ogsurfacePoints = self.vertices[:]
        self.ogsurfacePointsNormals = self.vnormals[:]#helper.FixNormals(self.vertices,self.faces,self.vnormals,fn=self.fnormals)
        mat = helper.getTransformation(self.ref_obj)
            #c4dmat = poly.GetMg()
            #mat,imat = self.c4dMat2numpy(c4dmat)
        self.normals= helper.FixNormals(self.vertices,self.faces,self.vnormals,fn=self.fnormals)
        self.ogsurfacePointsNormals = helper.ApplyMatrix(numpy.array(self.normals),helper.ToMat(mat))
#        faces = self.faces[:]
#        self.createSurfacePoints(maxl=grid.gridSpacing/2.0)
        surfacePoints = srfPts = self.ogsurfacePoints
        print (len(self.ogsurfacePointsNormals),self.ogsurfacePointsNormals)
        self.OGsrfPtsBht = bht =  bhtreelib.BHtree(tuple(srfPts), None, 10)

        res = numpy.zeros(len(srfPts),'f')
        dist2 = numpy.zeros(len(srfPts),'f')

        number = self.number
        ogNormals = numpy.array(self.ogsurfacePointsNormals)
        insidePoints = []

        # find closest off grid surface point for each grid point 
        #FIXME sould be diag of compartment BB inside fillBB
        grdPos = grid.masterGridPositions
        returnNullIfFail = 0
        closest = bht.closestPointsArray(tuple(grdPos), diag, returnNullIfFail)#diag is  cutoff ? meanin max distance ?
        
        self.closestId = closest
        t1=time()
        helper.resetProgressBar()
#        helper.progressBar(label="checking point %d" % point)
#       what abou intractive display ?
        if display :
            sph = helper.Sphere("gPts",res=10,radius=20.0)[0]
            sph2 = helper.Sphere("sPts",res=10,radius=20.0)[0]
            cylN = helper.oneCylinder("normal",[0.,0.,0.],[1.0,1.0,1.0],radius=20.0)
            cyl2 = helper.oneCylinder("V",[0.,0.,0.],[1.0,1.0,1.0],radius=20.0)
            helper.changeObjColorMat(sph2,(0.,0.,1.))
        for ptInd in xrange(len(grdPos)):#len(grdPos)):
            # find closest OGsurfacepoint
            if display :
                helper.changeObjColorMat(sph,(1.,1.,1.))
                helper.changeObjColorMat(cylN,(1.,0.,0.))
            inside = False
            t2=time()
            gx, gy, gz = grdPos[ptInd]
            sptInd = closest[ptInd]#this is a vertices 
            if display :
                helper.setTranslation(sph,grdPos[ptInd])
                helper.setTranslation(sph2,srfPts[sptInd])
#            helper.update()
            if closest[ptInd]==-1:
                print("ouhoua, closest OGsurfacePoint = -1")
                pdb.set_trace()
            if sptInd < len(srfPts):
                sx, sy, sz = srfPts[sptInd]
                d = sqrt( (gx-sx)*(gx-sx) + (gy-sy)*(gy-sy) +
                          (gz-sz)*(gz-sz))
            else :
                t0=times()
#                try :
                n=bht.closePointsDist2(tuple(grdPos[ptInd]),diag,res,dist2)#wthis is not working  
                d = min(dist2[0:n])
                sptInd = res[tuple(dist2).index(d)]
#                except :
                    #this is quite long
                    #what about C4d/host
#                delta = ((numpy.array(srfPts)-numpy.array(grdPos[ptInd]))**2).sum(axis=1)  # compute distances
#                ndx = delta.argsort() # indirect sort 
#                d = delta[ndx[0]]
#                sptInd = ndx[0]
#                    delta = numpy.array(srfPts)-numpy.array(grdPos[ptInd])
#                    delta *= delta
#                    distA = numpy.sqrt( delta.sum(1) )
#                    d = min(distA)
#                print('distance time', time()-t0)
                sptInd = list(distA).index(d)
                sx, sy, sz = srfPts[sptInd]
            if distances[ptInd]>d: distances[ptInd] = d

            if self.fnormals is not None and useFix:
                #too slow
                facesN = self.getFacesNfromV(sptInd,ext=1)
                #now lets get all fnormals and averge them
                n = nx, ny, nz = numpy.average( numpy.array(facesN),0)
#            print (faces)
        
            # check if ptInd in inside
            else :
                n = nx, ny, nz = numpy.array(ogNormals[sptInd])
#            vRayCollidePos = iRT.f_ray_intersect_polyhedron(numpy.array(grdPos[ptInd]), numpy.array(srfPts[sptInd]), self.ref_obj, 0,point = ptInd);
#            if (vRayCollidePos %  2):
#                print ("inside")
#                inside = True
#                idarray[ptInd] = -number
#                insidePoints.append(grdPos[ptInd]) 
#            vnpos = numpy.array(npost[sptInd])
            facesN = self.getVNfromF(sptInd)
            d1 = helper.measure_distance(numpy.array(grdPos[ptInd]),numpy.array(srfPts[sptInd])+(n*0.00001))
            d2 = helper.measure_distance(numpy.array(grdPos[ptInd]),numpy.array(srfPts[sptInd]))
            print ("gridpont distance from surf normal %0.10f from surf  %0.10f closer to snormal %s" % (d1,d2,str(d1 < d2)))
#             check on what side of the surface point the grid point is
            vptos = numpy.array(srfPts[sptInd]) - numpy.array(grdPos[ptInd])
            if display : 
#                helper.updateOneCylinder("normal",[0.,0.,0.],(n*spacing),radius=1.0)#srfPts[sptInd],numpy.array(srfPts[sptInd])+(n*spacing*10.0),radius=10.0)
#                helper.updateOneCylinder("V",[0.,0,0.],vptos,radius=1.0)#srfPts[sptInd],numpy.array(srfPts[sptInd])+(v*spacing*10.0),radius=10.0)
                helper.updateOneCylinder("normal",srfPts[sptInd],numpy.array(srfPts[sptInd])+(n*spacing*10.0),radius=10.0)
                helper.updateOneCylinder("V",srfPts[sptInd],numpy.array(srfPts[sptInd])+(vptos*spacing*10.0),radius=10.0)
                helper.update()
            dots=[]
            vptos=helper.normalize(vptos)
            for fn in facesN :
                dot = numpy.dot(vptos,fn)
                dots.append(dot)
                if display : 
                    helper.updateOneCylinder("normal",srfPts[sptInd],numpy.array(srfPts[sptInd])+(fn*spacing*10.0),radius=10.0)
                    helper.update()
            gr=numpy.greater(dots, 0.0)
#            print dots
#            print gr
            include = True
            if True in gr and False in gr:
                include = False 
            dot = numpy.dot(vptos,n)#project vptos on n -1 0 1 
            vx,vy,vz = (gx-sx, gy-sy, gz-sz)
            dot2 = vx*nx + vy*ny + vz*nz
            a = helper.angle_between_vectors(vptos,n)
            print (dot,dot2,math.degrees(a),include)
#            if math.degrees(a) > 250. :#and math.degrees(a) <= 271 :
#                print (dot,dot2,a,math.degrees(a))
#            if a > (math.pi/2.)+0.1 and a < (math.pi+(math.pi/2.)):#<= gave he outside points 
            if dot > 0 and a < math.pi/2.0 and include :#and d1 > d2 :#and dot < (-1.*10E-5): # inside ?
                print("INSIDE",dot,dot2,a,math.degrees(a))  
                
#                print (grdPos[ptInd],srfPts[sptInd])
#                print ("point ",v," normal ",n)
#                print ("inside",dot,dot2,a,math.degrees(a),helper.vector_norm(v),helper.vector_norm(n))
                #and the point is actually inside the mesh bounding box
                inside = True
#                if self.checkinside :
#                    inside  = self.checkPointInsideBB(grdPos[ptInd])
                #this is not working for a plane, or any unclosed compartment...
                if inside :
                    idarray[ptInd] = -number
                    insidePoints.append(grdPos[ptInd]) 
                if display : 
                    helper.changeObjColorMat(sph,(1.,0.,0.))
                    helper.update()                    
                    res = helper.drawQuestion(title="Inside?",question="%0.2f %0.2f %0.2f %0.2f %s" % (d1,d2,a,math.degrees(a),str(inside)))
                    if not res :
                        return insidePoints, surfacePoints
#                sleep(5.0)
                
            p=(ptInd/float(len(grdPos)))*100.0
            helper.progressBar(progress=int(p),label=str(ptInd)+"/"+str(len(grdPos))+" inside "+str(inside))

#            print('inside time', time()-t2)
        print('total time', time()-t1)
        return insidePoints, surfacePoints

    def getSurfaceInnerPointsPandaRay(self,boundingBox,spacing,display = True,useFix=False):
        """
        Only compute the inner point. No grid.
        This is independant from the packing. Help build ingredient sphere tree and representation 
        """                
        #should use the ray and see if it gave better reslt
        from autopack.pandautil import PandaUtil
        pud = PandaUtil()
        from autopack.Environment import Grid        
        grid = Grid()
        grid.boundingBox = boundingBox
        grid.gridSpacing = spacing# = self.smallestProteinSize*1.1547  # 2/sqrt(3)????
        t=time()
        helper.progressBar(label="BuildGRid")
        grid.gridVolume,grid.nbGridPoints = grid.computeGridNumberOfPoint(boundingBox,spacing)
        grid.create3DPointLookup()
        nbPoints = grid.gridVolume
        grid.gridPtId = [0]*nbPoints
        xl,yl,zl = boundingBox[0]
        xr,yr,zr = boundingBox[1]
        # distToClosestSurf is set to self.diag initially
#        grid.diag = diag = vlen( vdiff((xr,yr,zr), (xl,yl,zl) ) )
#        grid.distToClosestSurf = [diag]*nbPoints        
#        distances = grid.distToClosestSurf
#        idarray = grid.gridPtId
#        diag = grid.diag
        grdPos = grid.masterGridPositions
        insidePoints = []
        surfacePoints = self.vertices
        NPT=len(grdPos)
        meshnode = pud.addMeshRB(self.vertices, self.faces)
        #then sed ray from pointgrid to closest surface oint and see if collide ?
        #distance ? dot ? angle 
        grid.diag = diag = vlen( vdiff((xr,yr,zr), (xl,yl,zl) ) )
        grid.distToClosestSurf = [diag]*nbPoints        
        distances = grid.distToClosestSurf
        idarray = grid.gridPtId
        diag = grid.diag
        
        from bhtree import bhtreelib
        self.ogsurfacePoints = self.vertices[:]
        self.ogsurfacePointsNormals = helper.FixNormals(self.vertices,self.faces,self.vnormals,fn=self.fnormals)
#        faces = self.faces[:]
#        self.createSurfacePoints(maxl=grid.gridSpacing)
        surfacePoints = srfPts = self.ogsurfacePoints
        print (len(self.ogsurfacePointsNormals),self.ogsurfacePointsNormals)
        self.OGsrfPtsBht = bht =  bhtreelib.BHtree(tuple(srfPts), None, 10)

        res = numpy.zeros(len(srfPts),'f')
        dist2 = numpy.zeros(len(srfPts),'f')

        number = self.number
        ogNormals = numpy.array(self.ogsurfacePointsNormals)
        insidePoints = []

        # find closest off grid surface point for each grid point 
        #FIXME sould be diag of compartment BB inside fillBB
        grdPos = grid.masterGridPositions
        returnNullIfFail = 0
        closest = bht.closestPointsArray(tuple(grdPos), diag, returnNullIfFail)#diag is  cutoff ? meanin max distance ?
        
        self.closestId = closest
        t1=time()
        helper.resetProgressBar()
#        helper.progressBar(label="checking point %d" % point)
#       what abou intractive display ?
        if display :
            sph = helper.Sphere("gPts",res=10,radius=20.0)[0]
            sph2 = helper.Sphere("sPts",res=10,radius=20.0)[0]
            sph3 = helper.Sphere("hitPos",res=10,radius=20.0)[0]
            cylN = helper.oneCylinder("normal",[0.,0.,0.],[1.0,1.0,1.0],radius=20.0)
            cyl2 = helper.oneCylinder("V",[0.,0.,0.],[1.0,1.0,1.0],radius=20.0)
            helper.changeObjColorMat(sph2,(0.,0.,1.))
        for ptInd in xrange(len(grdPos)):#len(grdPos)):
            inside = False
            sptInd=closest[ptInd]
            v =- numpy.array(grdPos[ptInd]) + numpy.array(srfPts[closest[ptInd]])
            an = nx, ny, nz = numpy.array(ogNormals[sptInd])
#            start = Point3(grdPos[i][0],grdPos[i][1],grdPos[i][2])
            if display :
                helper.setTranslation(sph,grdPos[ptInd])
                helper.setTranslation(sph2,srfPts[closest[ptInd]])
                helper.update()  
#            end = Point3(srfPts[closest[i]][0]*diag,srfPts[closest[i]][1]*diag,srfPts[closest[i]][2]*diag)
            #raycats and see what it it on the mesh
            #or result = world.sweepTestClosest(shape, tsFrom, tsTo, penetration)
            res = pud.rayCast(grdPos[ptInd],(numpy.array(grdPos[ptInd])+v)*99999,closest=True)#world.rayTestAll(start, end)
            #can we get the number of hit?
            if res.hasHit():
                h = res
#                hit=res.getHits() 
##                for h in hit :
#                if len(hit):
#                h = hit[0]
                n = numpy.array(h.getHitNormal())
                a = helper.angle_between_vectors(v,n)
                dot = numpy.dot(v,n)
                dot2 = numpy.dot(an,v)
                a2 = helper.angle_between_vectors(-v,an)
                print ("hit with ",a,math.degrees(a),a2,math.degrees(a2),dot,dot2)
                if display : 
                    helper.setTranslation(sph3,numpy.array(h.getHitPos()))
                    helper.updateOneCylinder("normal",srfPts[sptInd],numpy.array(srfPts[sptInd])+(n*spacing*10.0),radius=10.0)
                    helper.updateOneCylinder("V",grdPos[ptInd],numpy.array(grdPos[ptInd])+(v),radius=10.0)
                    helper.update()
#                    if dot < 0 :#and dot < (-1.*10E-5): # inside ?
                if dot < 0.0 and dot2 < 0.0 :#a2 < (math.pi/2.)+0.1 and a > (math.pi/2.):# and a < (math.pi/2.) :#and a > (math.pi+(math.pi/2.)):
                    print("INSIDE",dot,a,math.degrees(a))
                    inside = True
                    if inside :
                        idarray[ptInd] = -number
                        insidePoints.append(grdPos[ptInd]) 
                    if display : 
                        helper.changeObjColorMat(sph,(1.,0.,0.))
                        helper.update()                    
                        res = helper.drawQuestion(title="Inside?",question="%0.2f %0.2f %0.2f %0.2f %s" % (dot,dot2,a,math.degrees(a),str(inside)))
                        if not res :
                            return insidePoints, surfacePoints
            p=(ptInd/float(len(grdPos)))*100.0
            helper.progressBar(progress=int(p),label=str(ptInd)+"/"+str(len(grdPos))+" inside "+str(inside))
                   
        return insidePoints, surfacePoints
        
    def getSurfaceInnerPointsPanda(self,boundingBox,spacing,display = True,useFix=False):
        """
        Only compute the inner point. No grid.
        This is independant from the packing. Help build ingredient sphere tree and representation 
        """                
        #work for small object
        from autopack.pandautil import PandaUtil
        pud = PandaUtil()
        from autopack.Environment import Grid        
        grid = Grid()
        grid.boundingBox = boundingBox
        grid.gridSpacing = spacing# = self.smallestProteinSize*1.1547  # 2/sqrt(3)????
        t=time()
        helper.progressBar(label="BuildGRid")
        grid.gridVolume,grid.nbGridPoints = grid.computeGridNumberOfPoint(boundingBox,spacing)
        grid.create3DPointLookup()
        nbPoints = grid.gridVolume
        grid.gridPtId = [0]*nbPoints
        xl,yl,zl = boundingBox[0]
        xr,yr,zr = boundingBox[1]
        # distToClosestSurf is set to self.diag initially
#        grid.diag = diag = vlen( vdiff((xr,yr,zr), (xl,yl,zl) ) )
#        grid.distToClosestSurf = [diag]*nbPoints        
#        distances = grid.distToClosestSurf
#        idarray = grid.gridPtId
#        diag = grid.diag
        grdPos = grid.masterGridPositions
        insidePoints = []
        surfacePoints = self.vertices
        NPT=len(grdPos)
        rads = [spacing,]*NPT
        helper.progressBar(label="BuildWorldAndNode")
        t=time()
        def addSphere(r,pos,i):
            node = pud.addSingleSphereRB(r,name=str(i))
            node.setPos(pos[0],pos[1],pos[2])
            helper.progressBar(progress=int((i/float(NPT))*100.0),label=str(i)+"/"+str(NPT))           
            return node
        nodes =[addSphere(rads[i],grdPos[i],i) for i in range(NPT)  ]
#        node = pud.addMultiSphereRB(rads,grdPos)
        helper.progressBar(label="OK SPHERE %0.2f" % (time()-t))# ("time sphere ",time()-t)
        t=time()
        #add the mesh
        meshnode = pud.addMeshRB(self.vertices, self.faces)
        helper.progressBar(label="OK MESH %0.2f" % (time()-t))#
        #computeCollisionTest
        t=time()   
        iPtList=[]
        meshcontacts = pud.world.contactTest(meshnode.node())
        N=meshcontacts.getNumContacts()
        for ct in meshcontacts.getContacts():
            m=ct.getManifoldPoint ()
            d = m.getDistance ()
            print (ct.getNode0().getName(), ct.getNode1().getName(),d)
            i = eval(ct.getNode0().getName())
            if i not in iPtList :
                insidePoints.append(grdPos[i])
                iPtList.append(i)
        print ("N",len(insidePoints),NPT)
        print ("time contact",time()-t)
        return insidePoints, surfacePoints

    def printFillInfo(self):
        """ print som info about the compartment and its recipe"""
        print('compartment %d'%self.number)
        r = self.surfaceRecipe
        if r is not None:
            print('    surface recipe:')
            r.printFillInfo('        ')

        r = self.innerRecipe
        if r is not None:
            print('    interior recipe:')
            r.printFillInfo('        ')
        
