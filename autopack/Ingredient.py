# -*- coding: utf-8 -*-
"""
############################################################################
#
# autoPACK Authors: Graham T. Johnson, Mostafa Al-Alusi, Ludovic Autin, 
#   and Michel Sanner
#   Based on COFFEE Script developed by Graham Johnson 
#    between 2005 and 2010 
#   with assistance from Mostafa Al-Alusi in 2009 and periodic input 
#   from Arthur Olson's Molecular Graphics Lab
#
# Ingredient.py Authors: Graham Johnson & Michel Sanner with 
#  editing/enhancement from Ludovic Autin
#
# Translation to Python initiated March 1, 2010 by Michel Sanner 
#  with Graham Johnson
#
# Class restructuring and organization: Michel Sanner
#
# Copyright: Graham Johnson ©2010
#
# This file "Ingredient.py" is part of autoPACK, cellPACK.
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
############################################################################
@author: Graham Johnson, Ludovic Autin, & Michel Sanner


# Hybrid version merged from Graham's Sept 2011 and Ludo's April 2012 
# version on May 16, 2012
# Updated with Correct Sept 25, 2011 thesis version on July 5, 2012

# TODO: Describe Ingredient class here at high level
"""
try :
    from scipy import spatial
except :
    autopack = None
import numpy
from numpy import matrix
#, weakref
from math import sqrt, pi,sin, cos, asin
from bhtree import bhtreelib
from random import uniform, gauss,random
from time import time,sleep
import math
from .ray import vlen, vdiff, vcross
from RAPID import RAPIDlib
from autopack.transformation import euler_from_matrix,matrixToEuler,euler_matrix,superimposition_matrix,rotation_matrix,affine_matrix_from_points
from autopack.transformation import quaternion_from_matrix
#RAPID require a uniq mesh. not an empty or an instance
#need to combine the vertices and the build the rapid model

#combining panda place with rapid place is not working properly


import sys

try :
    import urllib.request as urllib# , urllib.parse, urllib.error
except :
    import urllib

import autopack
from autopack import checkURL

AFDIR = autopack.__path__[0]#working dir ?
verbose = autopack.verbose
try :
    helper = autopack.helper
except :
    helper = None
print ("helper in Ingredient is "+str(helper))
helper = autopack.helper
reporthook = None
if helper is not None:        
    reporthook=helper.reporthook

#import numpy.oldnumeric as N
degtorad = pi/180.
KWDS = {   
                        "molarity":{"type":"float","name":"molarity","default":0,"value":0,"min":0,"max":500,"description":"molarity"}, 
                        "nbMol":{"type":"int","name":"nbMol","default":0,"value":0,"min":0,"max":50000,"description":"nbMol"},
                        "overwrite_nbMol_value":{"type":"int","name":"overwrite_nbMol_value","default":0,"value":0,"min":0,"max":50000,"description":"nbMol"},                        
                        "encapsulatingRadius":{"type":"float","name":"encapsulatingRadius","default":5,"value":5,"min":0,"max":500,"description":"encapsulatingRadius"}, 
                        "radii":{"type":"float"}, 
                        "positions":{"type":"vector"}, 
                        "positions2":{"type":"vector"},
                        "sphereFile":{"type":"string"}, 
                        "packingPriority":{"type":"float"}, 
                        "name":{"type":"string"}, 
                        "pdb":{"type":"string"}, 
                        "color":{"type":"vector"},
                        "principalVector":{"type":"vector"},
                        "meshFile":{"type":"string"},
                        "meshName":{"type":"string"},
                        "use_mesh_rb":{"name":"use_mesh_rb","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"use mesh for collision"},                             
                        "coordsystem":{"name":"coordsystem","type":"string","value":"left","default":"left","description":"coordinate system of the files"},
#                        "meshObject":{"type":"string"},
                        "Type":{"type":"string"},
                        "jitterMax":{"name":"jitterMax","value":[1.,1.,1.],"default":[1.,1.,1.],"min":0,"max":1,"type":"vector","description":"jitterMax"},
                        "nbJitter":{"name":"nbJitter","value":5,"default":5,"type":"int","min":0,"max":50,"description":"nbJitter"},
                        "perturbAxisAmplitude":{"name":"perturbAxisAmplitude","value":0.1,"default":0.1,"min":0,"max":1,"type":"float","description":"perturbAxisAmplitude"},
                        "useRotAxis":{"name":"useRotAxis","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"useRotAxis"},                             
                        "rotAxis":{"name":"rotAxis","value":[0.,0.,0.],"default":[0.,0.,0.],"min":0,"max":1,"type":"vector","description":"rotAxis"},
                        "rotRange":{"name":"rotRange","value":6.2831,"default":6.2831,"min":0,"max":12,"type":"float","description":"rotRange"},
                        
                        "useOrientBias":{"name":"useOrientBias","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"useOrientBias"},

                        "orientBiasRotRangeMin":{"name":"orientBiasRotRange","value":-pi,"default":-pi,"min":-pi,"max":pi,"type":"float","description":"orientBiasRotRangeMin"},
                        "orientBiasRotRangeMax":{"name":"orientBiasRotRange","value":pi,"default":pi,"min":-pi,"max":pi,"type":"float","description":"orientBiasRotRangeMax"},
                        
                        "rejectionThreshold":{"name":"rejectionThreshold","value":30,"default":30,"type":"float","min":0,"max":10000,"description":"rejectionThreshold"},
                        

                        "principalVector":{"name":"principalVector","value":[0.,0.,0.],"default":[0.,0.,0.],"min":-1,"max":1,"type":"vector","description":"principalVector"},
                        "cutoff_boundary":{"name":"cutoff_boundary","value":1.0,"default":1.0,"min":0.,"max":50.,"type":"float","description":"cutoff_boundary"},
                        "cutoff_surface":{"name":"cutoff_surface","value":5.0,"default":5.0,"min":0.,"max":50.,"type":"float","description":"cutoff_surface"},
                        "placeType":{"name":"placeType","value":"jitter","values":autopack.LISTPLACEMETHOD,"min":0.,"max":0.,
                                        "default":"jitter","type":"liste","description":"placeType"},
                        "packingMode":{"name":"packingMode","type":"string"},                     
                        "useLength":{"name":"useLength","type":"float"},
                        "length":{"name":"length","type":"float"},
                        "uLength":{"name":"uLength","type":"float"},
                        "closed":{"name":"closed","type":"bool"},
                        "biased":{"name":"biased","type":"float"},
                        "marge":{"name":"marge","type":"float"},
                        "constraintMarge":{"name":"constraintMarge","type":"bool"},
                        "orientation":{"name":"orientation","type":"vector"},
                        "partners_name":{"name":"partners_name","type":"liste_string"},
                        "excluded_partners_name":{"name":"excluded_partners_name","type":"liste_string", "value":"[]"},
                        "partners_position":{"name":"partners_position","type":"liste_float", "value":"[]"},  
                        "proba_not_binding":{"name":"proba_not_binding","type":"float", "value":"0.5"},                      
                        "walkingMode":{"name":"walkingMode","type":"string"},
                        "gradient":{"name":"gradient","value":"","values":[],"min":0.,"max":0.,
                                        "default":"jitter","type":"liste","description":"gradient name to use if histo.use_gradient"},
                        "isAttractor":{"name":"isAttractor","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"isAttractor"},
                        "weight":{"name":"weight","value":0.2,"default":0.2,"min":0.,"max":50.,"type":"float","description":"weight"},
                        "proba_binding":{"name":"proba_binding","value":0.5,"default":0.5,"min":0.,"max":1.0,"type":"float","description":"proba_binding"},
                        "proba_not_binding":{"name":"proba_not_binding","value":0.5,"default":0.5,"min":0.0,"max":1.,"type":"float","description":"proba_not_binding"},
                        "compMask":{"name":"compMask","value":"0","values":"0","min":0.,"max":0.,"default":'0',"type":"string","description":"allowed compartments"},
                        "use_rbsphere":{"name":"use_rbsphere","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"use sphere instead of cylinder for collision"},                             
                        "properties":{"name":"properties","value":{},"default":{},"min":0.,"max":1.0,"type":"dic","description":"properties"},        
                    }
#should use transform.py instead
def getNormedVectorOnes(a):
 n= a/numpy.linalg.norm(a)
 return numpy.round(n)

def getNormedVectorU(a):
 return a/numpy.linalg.norm(a)
    
def getNormedVector(a,b):
 return (b-a)/numpy.linalg.norm(b-a)

def getDihedral(a,b,c,d):
 v1 = getNormedVector(a, b)
 v2 = getNormedVector(b, c)
 v3 = getNormedVector(c, d)
 v1v2 = numpy.cross(v1,v2)
 v2v3 = numpy.cross(v2,v3)
 return angle_between_vectors(v1v2,v2v3)

def angle_between_vectors(v0, v1, directed=True, axis=0):
    """Return angle between vectors.
    If directed is False, the input vectors are interpreted as undirected axes,
    i.e. the maximum angle is pi/2.
    >>> a = angle_between_vectors([1, -2, 3], [-1, 2, -3])
    >>> numpy.allclose(a, math.pi)
    True
    >>> a = angle_between_vectors([1, -2, 3], [-1, 2, -3], directed=False)
    >>> numpy.allclose(a, 0)
    True
    >>> v0 = [[2, 0, 0, 2], [0, 2, 0, 2], [0, 0, 2, 2]]
    >>> v1 = [[3], [0], [0]]
    >>> a = angle_between_vectors(v0, v1)
    >>> numpy.allclose(a, [0, 1.5708, 1.5708, 0.95532])
    True
    >>> v0 = [[2, 0, 0], [2, 0, 0], [0, 2, 0], [2, 0, 0]]
    >>> v1 = [[0, 3, 0], [0, 0, 3], [0, 0, 3], [3, 3, 3]]
    >>> a = angle_between_vectors(v0, v1, axis=1)
    >>> numpy.allclose(a, [1.5708, 1.5708, 1.5708, 0.95532])
    True
    """
    v0 = numpy.array(v0, dtype=numpy.float64, copy=False)
    v1 = numpy.array(v1, dtype=numpy.float64, copy=False)
    dot = numpy.sum(v0 * v1, axis=axis)
    dot /= vector_norm(v0, axis=axis) * vector_norm(v1, axis=axis)
    return numpy.arccos(dot if directed else numpy.fabs(dot))
    
def vector_norm(data, axis=None, out=None):
    """Return length, i.e. Euclidean norm, of ndarray along axis.
    >>> v = numpy.random.random(3)
    >>> n = vector_norm(v)
    >>> numpy.allclose(n, numpy.linalg.norm(v))
    True
    >>> v = numpy.random.rand(6, 5, 3)
    >>> n = vector_norm(v, axis=-1)
    >>> numpy.allclose(n, numpy.sqrt(numpy.sum(v*v, axis=2)))
    True
    >>> n = vector_norm(v, axis=1)
    >>> numpy.allclose(n, numpy.sqrt(numpy.sum(v*v, axis=1)))
    True
    >>> v = numpy.random.rand(5, 4, 3)
    >>> n = numpy.empty((5, 3))
    >>> vector_norm(v, axis=1, out=n)
    >>> numpy.allclose(n, numpy.sqrt(numpy.sum(v*v, axis=1)))
    True
    >>> vector_norm([])
    0.0
    >>> vector_norm([1])
    1.0
    """
    data = numpy.array(data, dtype=numpy.float64, copy=True)
    if out is None:
        if data.ndim == 1:
            return math.sqrt(numpy.dot(data, data))
        data *= data
        out = numpy.atleast_1d(numpy.sum(data, axis=axis))
        numpy.sqrt(out, out)
        return out
    else:
        data *= data
        numpy.sum(data, axis=axis, out=out)
        numpy.sqrt(out, out)

def SphereHalton(n, p2,marge=math.pi):
    p=0.0
    t=0.0
    st=0.0
    phi=0.0
    phirad=0.0
    ip=0.0
    k=0
    kk=0
#    pos=0
    result=[]
    a=0
    for k in range(n):  
        #for (k=0, pos=0 ; k<n ; k++)
        t = 0;   
        p=0.5
        kk=k
#        for (p=0.5, kk=k ; kk ; p*=0.5, kk>>=1)
        while (kk):            
            if (kk & 1): #// kk mod 2 == 1
                t += p
            kk>>=1
            p*=0.5
        t = 2.0 * t - 1.0; #// map from [0,1] to [-1,1]
        st = math.sqrt(1.0-t*t);
        phi = 0;
        ip = 1.0/p2; #// inverse of p2
#        for (p=ip, kk=k ; kk ; p*=ip, kk/=p2) #// kk = (int)(kk/p2)
        p=ip
        kk=k        
        while(kk):
            a=kk % p2
            if (a):
                phi += a * p;
            kk/=p2
            p*=ip
        phirad = phi * 4.0 * marge; #// map from [0,0.5] to [0, 2 pi)
        px = st * math.cos(phirad);
        py = st * math.sin(phirad);
        pz = t;
        result.append([px,py,pz])
    return result

#from mglutil.math.rotax import rotax, rotVectToVect
def rotax( a, b, tau, transpose=1 ):
    """
    Build 4x4 matrix of clockwise rotation about axis a-->b
    by angle tau (radians).
    a and b are sequences of 3 floats each
    Result is a homogenous 4x4 transformation matrix.
    NOTE: This has been changed by Brian, 8/30/01: rotax now returns
    the rotation matrix, _not_ the transpose. This is to get
    consistency across rotax, mat_to_quat and the classes in
    transformation.py
    when transpose is 1 (default) a C-style rotation matrix is returned
    i.e. to be used is the following way Mx (opposite of OpenGL style which
    is using the FORTRAN style)
    """

    assert len(a) == 3
    assert len(b) == 3
    if tau <= -2*pi or tau >= 2*pi:
        tau = tau%(2*pi)

    ct = cos(tau)
    ct1 = 1.0 - ct
    st = sin(tau)

    # Compute unit vector v in the direction of a-->b. If a-->b has length
    # zero, assume v = (1,1,1)/sqrt(3).

    v = [b[0]-a[0], b[1]-a[1], b[2]-a[2]]
    s = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]
    if s > 0.0:
        s = sqrt(s)
        v = [v[0]/s, v[1]/s, v[2]/s]
    else:
        val = sqrt(1.0/3.0)
        v = (val, val, val)

    rot = numpy.zeros( (4,4), 'f' )
    # Compute 3x3 rotation matrix

    v2 = [v[0]*v[0], v[1]*v[1], v[2]*v[2]]
    v3 = [(1.0-v2[0])*ct, (1.0-v2[1])*ct, (1.0-v2[2])*ct]
    rot[0][0]=v2[0]+v3[0]
    rot[1][1]=v2[1]+v3[1]
    rot[2][2]=v2[2]+v3[2]
    rot[3][3] = 1.0;

    v2 = [v[0]*st, v[1]*st, v[2]*st]
    rot[1][0]=v[0]*v[1] * ct1-v2[2]
    rot[2][1]=v[1]*v[2] * ct1-v2[0]
    rot[0][2]=v[2]*v[0] * ct1-v2[1]
    rot[0][1]=v[0]*v[1] * ct1+v2[2]
    rot[1][2]=v[1]*v[2] * ct1+v2[0]
    rot[2][0]=v[2]*v[0] * ct1+v2[1]

    # add translation
    for i in (0,1,2):
        rot[3][i] = a[i]
    for j in (0,1,2):
        rot[3][i] = rot[3][i]-rot[j][i]*a[j]
    rot[i][3]=0.0

    if transpose:
        return rot
    else:
        return numpy.transpose(rot)




def rotVectToVect(vect1, vect2, i=None):
    """returns a 4x4 transformation that will align vect1 with vect2
vect1 and vect2 can be any vector (non-normalized)
"""
    v1x, v1y, v1z = vect1
    v2x, v2y, v2z = vect2
    
    # normalize input vectors
    norm = 1.0/sqrt(v1x*v1x + v1y*v1y + v1z*v1z )
    v1x *= norm
    v1y *= norm
    v1z *= norm    
    norm = 1.0/sqrt(v2x*v2x + v2y*v2y + v2z*v2z )
    v2x *= norm
    v2y *= norm
    v2z *= norm
    
    # compute cross product and rotation axis
    cx = v1y*v2z - v1z*v2y
    cy = v1z*v2x - v1x*v2z
    cz = v1x*v2y - v1y*v2x

    # normalize
    nc = sqrt(cx*cx + cy*cy + cz*cz)
    if nc==0.0:
        return [ [1., 0., 0., 0.],
                 [0., 1., 0., 0.],
                 [0., 0., 1., 0.],
                 [0., 0., 0., 1.] ]

    cx /= nc
    cy /= nc
    cz /= nc
    
    # compute angle of rotation
    if nc<0.0:
        if i is not None:
            print ('truncating nc on step:', i, nc)
        nc=0.0
    elif nc>1.0:
        if i is not None:
            print ('truncating nc on step:', i, nc)
        nc=1.0
        
    alpha = asin(nc)
    if (v1x*v2x + v1y*v2y + v1z*v2z) < 0.0:
        alpha = pi - alpha

    # rotate about nc by alpha
    # Compute 3x3 rotation matrix

    ct = cos(alpha)
    ct1 = 1.0 - ct
    st = sin(alpha)
    
    rot = [ [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.] ]


    rv2x, rv2y, rv2z = cx*cx, cy*cy, cz*cz
    rv3x, rv3y, rv3z = (1.0-rv2x)*ct, (1.0-rv2y)*ct, (1.0-rv2z)*ct
    rot[0][0] = rv2x + rv3x
    rot[1][1] = rv2y + rv3y
    rot[2][2] = rv2z + rv3z
    rot[3][3] = 1.0;

    rv4x, rv4y, rv4z = cx*st, cy*st, cz*st
    rot[0][1] = cx * cy * ct1 - rv4z
    rot[1][2] = cy * cz * ct1 - rv4x
    rot[2][0] = cz * cx * ct1 - rv4y
    rot[1][0] = cx * cy * ct1 + rv4z
    rot[2][1] = cy * cz * ct1 + rv4x
    rot[0][2] = cz * cx * ct1 + rv4y

    return rot

def ApplyMatrix(coords,mat):
    """
    Apply the 4x4 transformation matrix to the given list of 3d points.

    @type  coords: array
    @param coords: the list of point to transform.
    @type  mat: 4x4array
    @param mat: the matrix to apply to the 3d points

    @rtype:   array
    @return:  the transformed list of 3d points
    """

    #4x4matrix"
    mat = numpy.array(mat)
    coords = numpy.array(coords)
    one = numpy.ones( (coords.shape[0], 1), coords.dtype.char )
    c = numpy.concatenate( (coords, one), 1 )
    return numpy.dot(c, numpy.transpose(mat))[:, :3]

def bullet_checkCollision_mp(world, node1, node2):
#    world = 
#    node1 = histoVol.callFunction(self.histoVol.addRB,(self, jtrans, rotMatj,),{"rtype":self.Type},)
#    node2 = histoVol.callFunction(self.histoVol.addRB,(self, jtrans, rotMatj,),{"rtype":self.Type},)
    return world.contactTestPair(node1, node2).getNumContacts() > 0

#(v1,f1,numpy.array(rotMatj[:3,:3],'f'),numpy.array(jtrans,'f'),v2,f2,liste_nodes[n][2],liste_nodes[n][1],liste_nodes[n][3].name) 
def rapid_checkCollision_rmp(liste_input):
    from RAPID import RAPIDlib
    node1 = RAPIDlib.RAPID_model()
    node1.addTriangles(numpy.array(liste_input[0][0],'f'), numpy.array(liste_input[0][1],'i'))
    node2={}
    for inp in liste_input :
        if inp[-1] not in node2:
            node2[inp[-1]] = RAPIDlib.RAPID_model()
            node2[inp[-1]].addTriangles(numpy.array(inp[4],'f'), numpy.array(inp[5],'i'))
        RAPIDlib.RAPID_Collide_scaled(inp[2],inp[3], 1.0,node1, inp[6], inp[7], 1.0,
                                  node2[inp[-1]],RAPIDlib.cvar.RAPID_FIRST_CONTACT);                  
        if RAPIDlib.cvar.RAPID_num_contacts != 0 :
            return True
    return False
    
def rapid_checkCollision_mp(v1,f1,rot1,trans1,v2,f2,rot2,trans2):
    from RAPID import RAPIDlib
    node1 = RAPIDlib.RAPID_model()
    node1.addTriangles(numpy.array(v1,'f'), numpy.array(f1,'i'))
    node2 = RAPIDlib.RAPID_model()
    node2.addTriangles(numpy.array(v2,'f'), numpy.array(f2,'i'))
    RAPIDlib.RAPID_Collide_scaled(rot1,trans1, 1.0,node1, rot2, trans2, 1.0,
                                  node2,RAPIDlib.cvar.RAPID_FIRST_CONTACT);
#    print ("Num box tests: %d" % RAPIDlib.cvar.RAPID_num_box_tests)
#    print ("Num contact pairs: %d" % RAPIDlib.cvar.RAPID_num_contacts)
    #data3 = RAPIDlib.RAPID_Get_All_Pairs()
    return RAPIDlib.cvar.RAPID_num_contacts != 0
    
class Partner:
    def __init__(self,ingr,weight=0.0,properties=None):
        if type(ingr) is str :
            self.name = ingr
        else :
            self.name = ingr.name
        self.ingr = ingr
        self.weight = weight
        self.properties = {}
        self.distExpression = None
        if properties is not None:
            self.properties = properties
            
    def addProperties(self,name, value):
        self.properties[name] = value

    def getProperties(self,name):
        if name in self.properties:
            #if name == "pt1":
            #    return [0,0,0]
            #if name == "pt2":
            #    return [0,0,0]
            return self.properties[name]
        else :
            return None        

    def distanceFunction(self,d,expression=None,function=None):
        #default functino that can be overwrite or 
        #can provide an experssion which 1/d or 1/d^2 or d^2etc.w*expression
        #can provide directly a function that take as
        # arguments the w and the distance
        if expression is not None:
            val = self.weight * expression(d)
        elif function is not None :
            val = function(self.weight,d)
        else :
            val= self.weight*1.0/d
        return val

class IngredientInstanceDrop:
    def __init__(self, ptId, position, rotation, ingredient,rb=None):
        self.ptId = ptId
        self.position = position
        self.rotation = rotation
        self.ingredient = ingredient
        self.rigid_body = rb
        self.name = ingredient.name+str(ptId)
        x,y,z=position
        rad = ingredient.encapsulatingRadius
        self.bb=( [x-rad, y-rad, z-rad], [x+rad, y+rad, z+rad] )
        #maybe get bb from mesh if any ?
        if self.ingredient.mesh is not None :
            self.bb = autopack.helper.getBoundingBox(self.ingredient.mesh)
            for i in range(3):
                self.bb[0][i]=self.bb[0][i]+self.position[i]
                self.bb[1][i]=self.bb[1][i]+self.position[i]
                
class Agent:
    def __init__(self, name, concentration, packingMode='close', 
                 placeType="jitter", **kw):
        self.name=name
        self.concentration = concentration
        self.partners={}
        self.excluded_partners={}
        #the partner position is the local position 
        self.partners_position=[]
        if "partners_position" in kw :
            self.partners_position=kw["partners_position"]         
        self.partners_name=[]
        if "partners_name" in kw :
            self.partners_name=kw["partners_name"]
            if not self.partners_position:
                for i in self.partners_name:
                   self.partners_position.append([numpy.identity(4)])
        self.excluded_partners_name=[]
        if "excluded_partners_name" in kw :
            self.excluded_partners_name=kw["excluded_partners_name"]  
        assert packingMode in ['random', 'close', 'closePartner',
                               'randomPartner', 'gradient','hexatile','squaretile',
                               'triangletile']
        self.packingMode = packingMode

        #assert placeType in ['jitter', 'spring','rigid-body']
        self.placeType = placeType        
        self.mesh_3d = None
        self.isAttractor = False
        if "isAttractor" in kw:
            self.isAttractor=kw["isAttractor"]
        self.weight=0.2               #use for affinity ie partner.weight
        if "weight" in kw:
            self.weight=kw["weight"]
        self.proba_not_binding = 0.5  #chance to actually not bind
        if "proba_not_binding" in kw:
            self.proba_not_binding = kw["proba_not_binding"]
        self.proba_binding = 0.5      
        if "proba_binding" in kw:
            self.proba_binding = kw["proba_binding"]        
        self.force_random = False     #avoid any binding
        if "force_random" in kw:
            self.force_random = kw["force_random"]
        self.distFunction = None
        if "distFunction" in kw:
            self.distFunction = kw["distFunction"]
        self.distExpression = None
        if "distExpression" in kw:
            self.distExpression=kw["distExpression"]
        self.overwrite_distFunc = False     #overWrite
        if "overwrite_distFunc" in kw:
            self.overwrite_distFunc = kw["overwrite_distFunc"]        
        #chance to actually bind to any partner
        self.gradient=""
        if "gradient" in kw :
           self.gradient=kw["gradient"] 
        self.cb = None
        self.radii = None
        self.recipe = None #weak ref to recipe
        self.tilling = None
        
    def getProbaBinding(self,val =None):
        #get a value between 0.0 and 1.0and return the weight and success ?
        if val is None :
            val = random()
        if self.cb is not None :
            return self.cb(val)
        if val <= self.weight :
            return True,val
        else :
            return False,val

    def getPartnerweight(self,name):
        print("Deprecated use self.weight")
        partner = self.getPartner(name)
        w = partner.getProperties("weight")
        if w is not None :
            return w
            
    def getPartnersName(self):
        return list(self.partners.keys())
        
    def getPartner(self,name):
        if name in self.partners:
            return self.partners[name]
        else :
            return None
        
    def addPartner(self, ingr, weight=0.0, properties=None):
        self.partners[ingr.name] = Partner(ingr, weight=weight, 
                                properties=properties)
        return self.partners[ingr.name]
        
    def getExcludedPartnersName(self):
        return list(self.excluded_partners.keys())

    def getExcludedPartner(self, name):
        if name in self.excluded_partners:
            return self.excluded_partners[name]
        else :
            return None

    def addExcludedPartner(self, name, properties=None):
        self.excluded_partners[name] = Partner(name, properties=properties)
         

    def sortPartner(self, listeP=None):
        if listeP is None :
            listeP=[]
            for i,ingr in list(self.partners.keys()):
                listeP.append([i,ingr])
        #extract ing name unic
        listeIngrInstance={}
        for i,ingr in listeP:
            if ingr.name not in listeIngrInstance:
                listeIngrInstance[ingr.name]=[ingr.weight,[]]
            listeIngrInstance[ingr.name][1].append(i)
        #sort according ingrediant binding weight (proba to bind)
        sortedListe = sorted(list(listeIngrInstance.items()), 
                                 key=lambda elem: elem[1][0])   
        #sortedListe is [ingr,(weight,(instances indices))]
        # sort by weight/min->max
        #wIngrList = []
        #for i,ingr in listeP:
            #need to sort by ingr.weight
        #    wIngrList.append([i,ingr,ingr.weight])
        #sortedListe = sorted(wIngrList, key=lambda elem: elem[2])   # sort by weight/min->max
#        print sortedListe
        return sortedListe

    def weightListByDistance(self, listePartner):
        probaArray=[]
        w=0.
        for i,part,dist in listePartner:
            if self.overwrite_distFunc:
                wd = part.weight
            else:
                wd = part.distanceFunction(dist, 
                                           expression=part.distExpression)
#            print "calc ",dist, wd
            probaArray.append(wd)
            w = w+wd
        #probaArray.append(self.proba_not_binding)
        #w=w+self.proba_not_binding
        return probaArray,w

    def getProbaArray(self, weightD, total):
        probaArray=[]
        final=0.
        for w in weightD:
            p = w/total
#            print "norma ",w,total,p
            final = final + p
            probaArray.append(final)
        probaArray[-1]=1.0
        return probaArray

    def getSubWeighted(self,weights):
        """
        From http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/
        This method is about twice as fast as the binary-search technique, 
        although it has the same complexity overall. Building the temporary 
        list of totals turns out to be a major part of the functions runtime.
        This approach has another interesting property. If we manage to sort 
        the weights in descending order before passing them to 
        weighted_choice_sub, it will run even faster since the random 
        call returns a uniformly distributed value and larger chunks of 
        the total weight will be skipped in the beginning.
        """
        rnd = random() * sum(weights)
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                return i,rnd
                        
    def pickPartner(self, mingrs, listePartner, currentPos=[0,0,0]):
        #listePartner is (i,partner,d)
        #wieght using the distance function
#        print "len",len(listePartner)
        weightD,total = self.weightListByDistance(listePartner)
#        print "w", weightD,total
        i,b = self.getSubWeighted(weightD)
#        probaArray = self.getProbaArray(weightD,total)
##        print "p",probaArray
#        probaArray=numpy.array(probaArray)
#        #where is random in probaArray->index->ingr
#        b = random() 
#        test = b < probaArray
#        i = test.tolist().index(True)
##        print "proba",i,test,(len(probaArray)-1)    
#        if i == (len(probaArray)-1) :
#            #no binding due to proba not binding....
#            print ("no binding due to proba")
#            return None,b
        ing_indice = listePartner[i][0]#i,part,dist
        ing = mingrs[2][ing_indice]#[2]
        print ("binding to "+ing.name)
        targetPoint= mingrs[0][ing_indice]#[0]     
        if self.compNum > 0 :
#            organelle = self.histoVol.compartments[abs(self.compNum)-1]
#            dist,ind = organelle.OGsrfPtsBht.query(targetPoint)
            targetPoint=self.histoVol.grid.getClosestFreeGridPoint(targetPoint,
                                        compId=self.compNum,ball=(ing.encapsulatingRadius+self.encapsulatingRadius)*2.0,
                                        distance=self.encapsulatingRadius*1.5)
            print ("target point free tree is ",targetPoint,self.encapsulatingRadius,ing.encapsulatingRadius)
        else :                
            #get closestFreePoint using freePoint and masterGridPosition
            #if self.placeType == "rigid-body" or self.placeType == "jitter":
                #the new point is actually tPt -normalise(tPt-current)*radius
            print("tP",ing_indice,ing.name,targetPoint,ing.radii[0][0])
            print("cP",currentPos)
            #what I need it the closest free point from the target ingredient
            v=numpy.array(targetPoint) - numpy.array(currentPos)
            s = numpy.sum(v*v)
            factor = ((v/math.sqrt(s)) * (ing.encapsulatingRadius+self.encapsulatingRadius))#encapsulating radus ?
            targetPoint =  numpy.array(targetPoint) - factor    
        print("tPa",targetPoint)
        return targetPoint,b
        
    def pickPartner_old(self, mingrs,listePartner, currentPos=[0,0,0]):
        #pick the highest weighted partner
        #pick one instance of this ingrediant 
        #(distance or random or density nb of instance) 
        #roll a dice to decide if bind or not
        #problem where put the cb function
        #listePartner is (i,partner,d)
        sorted_listePartner=self.sortPartner(listePartner)
        #sortedListe is [ingrP,(weight,(instances indices))]
        binding = False
        targetPoint = None
#        pickedIngr = None     
#        found = False
#        safetycutoff=10 #10 roll dice, save most prob
        #do we take the one with highest weight of binding, 
        #or do we choose with a dice
#        p = random()
#        for bindingIngr in sorted_listePartner:
#            if p < part[1][0] :
#                break
        bindingIngr = sorted_listePartner[0]
        #pick the instance random, or distance
        i = self.pickPartnerInstance(bindingIngr, mingrs, 
                                     currentPos=currentPos)
        #roll a dice to see if we bind
        b=random()
        if b < self.proba_binding:
            binding = True
        if binding:
            ing = mingrs[i][2]
            targetPoint= mingrs[i][0]            
            if self.placeType == "rigid-body":
                #the new point is actually 
                #tPt -normalise(tPt-current)*radius
                v=numpy.array(targetPoint) - numpy.array(currentPos)
                s = numpy.sum(v*v)
                factor = ((v/s)* ing.radii[0][0])
                targetPoint = numpy.array(targetPoint) - factor
        return targetPoint,b

    def pickPartnerInstance(self, bindingIngr, mingrs, currentPos=None):
        #bindingIngr is ingr,(weight,(instances indices))
#        print "bindingIngr ",bindingIngr,bindingIngr[1]
        if currentPos is None : #random mode
            picked_I = random()*len(bindingIngr[1][1])
            i = bindingIngr[1][1][picked_I]
        else : #pick closest one
            mind=99999999.9
            i=0
            for ind in bindingIngr[1][1]:
                v=numpy.array(mingrs[ind][0]) - numpy.array(currentPos)
                d = numpy.sum(v*v)
                if d < mind:
                    mind = d
                    i = ind
        return i
                

#the ingrediant should derive from a class of Agent
class Ingredient(Agent):
    """
    Base class for Ingredients that can be added to a Recipe.
    Ingredients provide:
        - a molarity used to compute how many to place
        - a generic density value
        - a unit associated with the density value
        - a jitter amplitude vector specifying by how much the jittering
        algorithm can move fro the grid position.
        - a number of jitter attempts
        - an optional color used to draw the ingredient default (white)
        - an optional name
        - an optional pdb ID
        - an optional packing priority. If omited the priority will be based
        on the radius with larger radii first
        ham here: (-)packingPriority object will pack from high to low one at a time
        (+)packingPriority will be weighted by assigned priority value
        (0)packignPriority will be weighted by complexity and appended to what is left
        of the (+) values
        - an optional princial vector used to align the ingredient
        - recipe will be a weakref to the Recipe this Ingredient belongs to
        - compNum is th compartment number (0 for cyto plasm, positive for compartment
        surface and negative compartment interior
        - Attributes used by the filling algorithm:
        - nbMol counts the number of palced ingredients during a fill
        - counter is the target numbr of ingredients to place
        - completion is the ratio of placxed/target
        - rejectionCounter is used to eliminate ingredients after too many failed
        attempts

    """
    def __init__(self, molarity=0.0, radii=None, positions=None, positions2=None,
                 sphereFile=None, packingPriority=0, name=None, pdb='????', 
                 color=None, nbJitter=5, jitterMax=(1,1,1),
                 perturbAxisAmplitude = 0.1, principalVector=(1,0,0),
                 meshFile=None, meshName=None,packingMode='random',placeType="jitter",
                 meshObject=None,nbMol=0,Type="MultiSphere",**kw):
        Agent.__init__(self, name, molarity, packingMode=packingMode, 
                       placeType=placeType, **kw)
        self.molarity = molarity
        self.packingPriority = packingPriority
        print (packingPriority,self.packingPriority)
        if name == None:
            name = "%f"% molarity
        print ("CREATE INGREDIENT",str(name),("rejectionThreshold" in kw))
        self.name = str(name)
        self.o_name = str(name)
        self.Type = Type
        self.pdb = pdb        #pmv ?
        self.transform_sources = None
        self.source=None
        #should deal with source of the object
        if "source" in kw :
            sources = kw["source"].keys()
            self.source= kw["source"]
            if "pdb" in sources :
                self.pdb =  kw["source"]["pdb"]
                self.transform_sources= kw["source"]["transform"]
        else :
           self.source={"pdb":self.pdb,"transform":{"center":True}} 
           self.transform_sources= {"transform":{"center":True}}
        self.color = color    # color used for sphere display
        if self.color == "None":
            self.color = None
        self.modelType='Spheres'
        self.rRot=[]
        self.tTrans=[]
        self.htrans=[]
        self.moving = None
        self.moving_geom = None
        self.rb_nodes = []#store rbnode. no more than X ?
        self.bullet_nodes =[None,None]#try only store 2, and move them when needd
        self.limit_nb_nodes = 50
        self.vi = autopack.helper
        self.minRadius = 0
        self.encapsulatingRadius = 0
        self.maxLevel = 1
        self.is_previous = False
        self.vertices=[]
        self.faces=[]
        #self._place = self.place
        children = []
        self.sphereFile = None
        if sphereFile is not None and str(sphereFile) != "None":
            sphereFileo = autopack.retrieveFile(sphereFile,cache="collisionTrees") 
            fileName, fileExtension = os.path.splitext(sphereFile)
            print ("sphereTree",sphereFileo)
            if sphereFileo is not None :
                self.sphereFile=sphereFile
                sphereFile=sphereFileo
                if fileExtension == ".mstr" : #BD_BOX format
                    data = numpy.loadtxt(sphereFileo,converters = {0: lambda s: 0})
                    positions = data[:,1:4]
                    radii = data[:,4]
                    self.minRadius = min(radii)
                    #np.apply_along_axis(np.linalg.norm, 1, c)
                    self.encapsulatingRadius = max(numpy.sqrt(numpy.einsum('ij,ij->i',positions,positions)))#shoud be max distance
                    positions=[positions]
                    radii=[radii]
                else :
                    rm, rM, positions, radii, children = self.getSpheres(sphereFileo)
                    if not len(radii):
                        self.minRadius = 1.0
                        self.encapsulatingRadius = 1.0
                    else :
                        # minRadius is used to compute grid spacing. It represents the
                        # smallest radius around the anchor point(i.e. 
                        # the point where the
                        # ingredient is dropped that needs to be free
                        self.minRadius = rm
                        # encapsulatingRadius is the radius of the sphere 
                        # centered at 0,0,0
                        # and encapsulate the ingredient
                        self.encapsulatingRadius = rM
                
        if positions is None or positions[0] is None or positions[0][0] is None:#[0][0]
            positions = [[[0,0,0]]]
            if radii is not None :    
                self.minRadius = [radii[0]]
                self.encapsulatingRadius = max(radii[0])
        else :
            if radii is not None :
                self.minRadius = min(radii[0])
                self.encapsulatingRadius = max(radii[0])
#        print "sphereFile",sphereFile
#        print "positions",positions,len(positions)
#        print "rad",radii,len(radii)
        if radii is not None and positions is not None:
            for r,c in zip(radii, positions):
                assert len(r)==len(c)
        
        if radii is not None :
            self.maxLevel = len(radii)-1
        if radii is None :
            radii = [[0]]
        self.radii = radii
        self.positions = positions
        self.positions2 = positions2
        self.children = children
        self.rbnode =  {} #keep the rbnode if any
        self.collisionLevel = self.maxLevel 
        # first level used for collision detection
        self.jitterMax = jitterMax 
        # (1,1,1) means 1/2 grid spacing in all directions
        self.nbJitter = nbJitter 
        # number of jitter attemps for translation

        self.perturbAxisAmplitude = perturbAxisAmplitude
        
        self.principalVector = principalVector

        self.recipe = None # will be set when added to a recipe
        self.compNum = None 
        self.compId_accepted=[]#if this list is defined, point picked outise the list are rejected
        #should be self.compNum per default
        # will be set when recipe is added to HistoVol 
        #added to a compartment
        self.overwrite_nbMol = False
        self.overwrite_nbMol_value = nbMol
        self.nbMol =  nbMol
        self.vol_nbmol=0
        # used by fill() to count placed molecules,overwrite if !=0
#        if nbMol != 0:
#            self.overwrite_nbMol = True
#            self.overwrite_nbMol_value = nMol
#            self.nbMol = nMol
        self.counter = 0      # target number of molecules for a fill
        self.completion = 0.0 # ratio of counter/nbMol
        self.rejectionCounter = 0
        self.verts=None
        self.rad=None
        self.rapid_model = None
        if self.encapsulatingRadius <= 0.0 or self.encapsulatingRadius < max(self.radii[0]):
            self.encapsulatingRadius = max(self.radii[0])#
        #TODO : geometry : 3d object or procedural from PDB
        #TODO : usekeyword resolution->options dictionary of res :
        #TODO : {"simple":{"cms":{"parameters":{"gridres":12}},
        #TODO :            "obj":{"parameters":{"name":"","filename":""}}
        #TODO :            }
        #TODO : "med":{"method":"cms","parameters":{"gridres":30}}
        #TODO : "high":{"method":"msms","parameters":{"gridres":30}}
        #TODO : etc...
        self.coordsystem="left"
        if "coordsystem" in kw:
            self.coordsystem = kw["coordsystem"]
        self.rejectionThreshold = 30
        if "rejectionThreshold" in kw:
            print ("rejectionThreshold",kw["rejectionThreshold"])
            self.rejectionThreshold = kw["rejectionThreshold"]

        #get the collision mesh
        self.meshFile = None
        self.meshName = meshName
        self.mesh = None
        self.meshObject= None
        if meshFile is not None:
            print ("OK, meshFile is not none, it is = ",meshFile,self.name,self.coordsystem)
            gname = self.name
            if self.meshName is not None :
                gname = self.meshName
            self.mesh = self.getMesh(meshFile, gname)#self.name)
            print ("OK got",self.mesh)
            if self.mesh is None :
                #display a message ?
                print ("no geometrie for ingredient " + self.name)
            #should we reparent it ?
            self.meshFile = meshFile
        elif meshObject is not None:
            self.mesh = meshObject
        if "encapsulatingRadius" in kw:
            #we force the encapsulatingRadius
            if autopack.helper.host != "3dsmax":
                self.encapsulatingRadius = kw["encapsulatingRadius"]
        if self.mesh is not None :
           self.getEncapsulatingRadius()
        #need to build the basic shape if one provided
        self.use_mesh_rb = False
        self.current_resolution="Low"#should come from data
        self.available_resolution=["Low","Med","High"]#0,1,2
        self.resolution_dictionary = {"Low":"","Med":"","High":""}
        if "resolution_dictionary" in kw :
            if kw["resolution_dictionary"] is not None:
                self.resolution_dictionary = kw["resolution_dictionary"]

        #how to get the geom of different res?
        self.representation = None
        self.representation_file = None

        self.useRotAxis = False    
        if "useRotAxis" in kw:
            self.useRotAxis = kw["useRotAxis"]
        self.rotAxis = None
        if "rotAxis" in kw:
            self.rotAxis = kw["rotAxis"]
            #this could define the biased
        self.rotRange = 6.2831
        if "rotRange" in kw:
            self.rotRange = kw["rotRange"]
                
                
        self.useOrientBias = False
        if "useOrientBias" in kw:
            self.useOrientBias = kw["useOrientBias"]
                
        self.orientBiasRotRangeMin = -pi
        if "orientBiasRotRangeMin" in kw:
            self.orientBiasRotRangeMin = kw["orientBiasRotRangeMin"]
                
        self.orientBiasRotRangeMax = -pi
        if "orientBiasRotRangeMax" in kw:
            self.orientBiasRotRangeMax = kw["orientBiasRotRangeMax"]

        #cutoff are used for picking point far from surface and boundary
        self.cutoff_boundary = None#self.encapsulatingRadius
        self.cutoff_surface = float(self.encapsulatingRadius)
        if "cutoff_boundary" in kw:
            self.cutoff_boundary = kw["cutoff_boundary"]
        if "cutoff_surface" in kw:
            if kw["cutoff_surface"] != 0.0 :
                self.cutoff_surface = float(kw["cutoff_surface"])
        self.properties ={}#four tout
        if "properties" in kw:
            self.properties = kw["properties"]
        
        self.compareCompartment = False
        self.compareCompartmentTolerance = 0
        self.compareCompartmentThreshold = 0.0
        
        self.updateOwnFreePts = False #work for rer python not ??
        self.haveBeenRejected = False
        
        self.distances_temp=[]
        self.centT = None #transformed position

        self.results =[]
#        if self.mesh is not None :
#            self.getData()


        #add tiling property ? as any ingredient coud tile as hexagon. It is just the packing type
        self.KWDS = { 
                        "overwrite_nbMol_value":{"type":"int","name":"overwrite_nbMol_value","default":0,"value":0,"min":0,"max":50000,"description":"nbMol"},
                        "molarity":{"type":"float","name":"molarity","default":0,"value":0,"min":0,"max":500,"description":"molarity"}, 
                        "nbMol":{"type":"int","name":"nbMol","default":0,"value":0,"min":0,"max":50000,"description":"nbMol"},
                        "encapsulatingRadius":{"type":"float","name":"encapsulatingRadius","default":5,"value":5,"min":0,"max":500,"description":"encapsulatingRadius"}, 
                        "radii":{"type":"float"}, 
                        "positions":{}, "positions2":{},
                        "sphereFile":{"type":"string"}, 
                        "packingPriority":{"type":"float"}, 
                        "name":{"type":"string"}, 
                        "pdb":{"type":"string"}, 
                        "source":{},
                        "color":{"type":"vector"},"principalVector":{"type":"vector"},
                        "meshFile":{"type":"string"}, 
                        "meshName":{"type":"string"}, 
                        "coordsystem":{"name":"coordsystem","type":"string","value":"left","default":"left","description":"coordinate system of the files"},
#                        "meshObject":{"type":"string"},
                        "principalVector":{"name":"principalVector","value":[0.,0.,0.],"default":[0.,0.,0.],"min":-1,"max":1,"type":"vector","description":"principalVector"},
                        "Type":{"type":"string"},
                        "jitterMax":{"name":"jitterMax","value":[1.,1.,1.],"default":[1.,1.,1.],"min":0,"max":1,"type":"vector","description":"jitterMax"},
                        "nbJitter":{"name":"nbJitter","value":5,"default":5,"type":"int","min":0,"max":50,"description":"nbJitter"},
                        "perturbAxisAmplitude":{"name":"perturbAxisAmplitude","value":0.1,"default":0.1,"min":0,"max":1,"type":"float","description":"perturbAxisAmplitude"},
                        "useRotAxis":{"name":"useRotAxis","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"useRotAxis"},                             
                        "rotAxis":{"name":"rotAxis","value":[0.,0.,0.],"default":[0.,0.,0.],"min":0,"max":1,"type":"vector","description":"rotAxis"},
                        "rotRange":{"name":"rotRange","value":6.2831,"default":6.2831,"min":0,"max":12,"type":"float","description":"rotRange"},
                        
                        
                        "useOrientBias":{"name":"useOrientBias","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"useOrientBias"},
                        
                        "orientBiasRotRangeMin":{"name":"orientBiasRotRange","value":-pi,"default":-pi,"min":-pi,"max":pi,"type":"float","description":"orientBiasRotRangeMin"},
                        "orientBiasRotRangeMax":{"name":"orientBiasRotRange","value":pi,"default":pi,"min":-pi,"max":pi,"type":"float","description":"orientBiasRotRangeMax"},

                        
                        
                        "cutoff_boundary":{"name":"cutoff_boundary","value":1.0,"default":1.0,"min":0.,"max":50.,"type":"float","description":"cutoff_boundary"},
                        "cutoff_surface":{"name":"cutoff_surface","value":5.0,"default":5.0,"min":0.,"max":50.,"type":"float","description":"cutoff_surface"},
                        "placeType":{"name":"placeType","value":"jitter","values":autopack.LISTPLACEMETHOD,"min":0.,"max":0.,
                                        "default":"jitter","type":"liste","description":"placeType"},
                        "use_mesh_rb":{"name":"use_mesh_rb","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"use mesh for collision"},                             
                        "rejectionThreshold":{"name":"rejectionThreshold","value":30,"default":30,"type":"float","min":0,"max":10000,"description":"rejectionThreshold"},
                        "partners_name":{"name":"partners_name","type":"liste_string", "value":"[]"},
                        "excluded_partners_name":{"name":"excluded_partners_name","type":"liste_string", "value":"[]"},
                        "partners_position":{"name":"partners_position","type":"liste_float", "value":"[]"},
                        "packingMode":{"name":"packingMode","value":"random","values":['random', 'close', 'closePartner',
                               'randomPartner', 'gradient','hexatile','squaretile','triangletile'],"min":0.,"max":0.,"default":'random',"type":"liste","description":"packingMode"},
                        "gradient":{"name":"gradient","value":"","values":[],"min":0.,"max":0.,
                                        "default":"jitter","type":"liste","description":"gradient name to use if histo.use_gradient"},
                        "proba_not_binding":{"name":"proba_not_binding","type":"float", "value":"0.5"},
                        "isAttractor":{"name":"isAttractor","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"isAttractor"},
                        "weight":{"name":"weight","value":0.2,"default":0.2,"min":0.,"max":50.,"type":"float","description":"weight"},
                        "proba_binding":{"name":"proba_binding","value":0.5,"default":0.5,"min":0.,"max":1.0,"type":"float","description":"proba_binding"},
                        "properties":{"name":"properties","value":{},"default":{},"min":0.,"max":1.0,"type":"dic","description":"properties"},
                        
                        }
        self.OPTIONS = {
                        "molarity":{"type":"float","name":"molarity","default":0,"value":0,"min":0,"max":500,"description":"molarity"}, 
                        "nbMol":{"type":"int","name":"nbMol","default":0,"value":0,"min":0,"max":50000,"description":"nbMol"},
                        "overwrite_nbMol_value":{"type":"int","name":"overwrite_nbMol_value","default":0,"value":0,"min":0,"max":50000,"description":"overwrite_nbMol_value"},
                        "radii":{}, 
                        "encapsulatingRadius":{"type":"float","name":"encapsulatingRadius","default":5,"value":5,"min":0,"max":500,"description":"encapsulatingRadius"}, 
                        "positions":{}, 
                        "positions2":{},
                        "sphereFile":{}, 
                        "packingPriority":{}, 
                        "name":{}, 
                        "pdb":{}, 
                        "source":{},
                        "color":{},
                        "principalVector":{"name":"principalVector","value":[0.,0.,0.],"default":[0.,0.,0.],"min":-1,"max":1,"type":"vector","description":"principalVector"},
                        "meshFile":{}, 
                        "meshObject":{},
                        "coordsystem":{"name":"coordsystem","type":"string","value":"right","default":"right","description":"coordinate system of the files"},
                        "use_mesh_rb":{"name":"use_mesh_rb","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"use mesh for collision"},                             
                        "rejectionThreshold":{"name":"rejectionThreshold","value":30,"default":30,"type":"float","min":0,"max":10000,"description":"rejectionThreshold"},
                        "jitterMax":{"name":"jitterMax","value":[1.,1.,1.],"default":[1.,1.,1.],"min":0,"max":1,"type":"vector","description":"jitterMax"},
                        "nbJitter":{"name":"nbJitter","value":5,"default":5,"type":"int","min":0,"max":50,"description":"nbJitter"},
                        "perturbAxisAmplitude":{"name":"perturbAxisAmplitude","value":0.1,"default":0.1,"min":0,"max":1,"type":"float","description":"perturbAxisAmplitude"},
#                         "principalVector":{"name":"principalVector","value":9999999,"default":99999999,"type":"vector_norm","description":"principalVector"},
                        "useRotAxis":{"name":"useRotAxis","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"useRotAxis"},                             
                        "rotAxis":{"name":"rotAxis","value":[0.,0.,0.],"default":[0.,0.,0.],"min":0,"max":1,"type":"vector","description":"rotAxis"},
                        "rotRange":{"name":"rotRange","value":6.2831,"default":6.2831,"min":0,"max":12,"type":"float","description":"rotRange"},
                        
                        
                        "useOrientBias":{"name":"useOrientBias","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"useOrientBias"},
                        
                        "orientBiasRotRangeMin":{"name":"orientBiasRotRange","value":-pi,"default":-pi,"min":-pi,"max":pi,"type":"float","description":"orientBiasRotRangeMin"},
                        "orientBiasRotRangeMax":{"name":"orientBiasRotRange","value":pi,"default":pi,"min":-pi,"max":pi,"type":"float","description":"orientBiasRotRangeMax"},

                        
                        "packingMode":{"name":"packingMode","value":"random","values":['random', 'close', 'closePartner',
                               'randomPartner', 'gradient','hexatile','squaretile','triangletile'],"min":0.,"max":0.,"default":'random',"type":"liste","description":"packingMode"},
                        "placeType":{"name":"placeType","value":"jitter","values":autopack.LISTPLACEMETHOD,"min":0.,"max":0.,
                                        "default":"jitter","type":"liste","description":"placeType"},
                        "gradient":{"name":"gradient","value":"","values":[],"min":0.,"max":0.,
                                        "default":"jitter","type":"liste","description":"gradient name to use if histo.use_gradient"},
                        "isAttractor":{"name":"isAttractor","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"isAttractor"},
                        "weight":{"name":"weight","value":0.2,"default":0.2,"min":0.,"max":50.,"type":"float","description":"weight"},
                        "proba_binding":{"name":"proba_binding","value":0.5,"default":0.5,"min":0.,"max":1.0,"type":"float","description":"proba_binding"},
                        "proba_not_binding":{"name":"proba_not_binding","value":0.5,"default":0.5,"min":0.0,"max":1.,"type":"float","description":"proba_not_binding"},
                        "cutoff_boundary":{"name":"cutoff_boundary","value":1.0,"default":1.0,"min":0.,"max":50.,"type":"float","description":"cutoff_boundary"},
                        "cutoff_surface":{"name":"cutoff_surface","value":5.0,"default":5.0,"min":0.,"max":50.,"type":"float","description":"cutoff_surface"},
                        "compareCompartment":{"name":"compareCompartment","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"compareCompartment"},
                        "compareCompartmentTolerance":{"name":"compareCompartmentTolerance","value":0.0,"default":0.0,"min":0.,"max":1.0,"type":"float","description":"compareCompartmentTolerance"},
                        "compareCompartmentThreshold":{"name":"compareCompartmentThreshold","value":0.0,"default":0.0,"min":0.,"max":1.0,"type":"float","description":"compareCompartmentThreshold"},
                        "partners_name":{"name":"partners_name","type":"liste_string", "value":"[]"},
                        "excluded_partners_name":{"name":"excluded_partners_name","type":"liste_string", "value":"[]"},
                        "partners_position":{"name":"partners_position","type":"liste_float", "value":"[]"},
                        "properties":{"name":"properties","value":{},"default":{},"min":0.,"max":1.0,"type":"dic","description":"properties"},
                        }

    def setTilling(self,comp):
        if self.packingMode == 'hexatile' :
            from autopack.hexagonTile import tileHexaIngredient
            #self.histoVol attached to compratmentsmes
            self.tilling = tileHexaIngredient(self,comp,
                                         self.encapsulatingRadius,
                                         init_seed=self.histoVol.seed_used)
        elif self.packingMode == 'squaretile' :
            from autopack.hexagonTile import tileSquareIngredient
            #self.histoVol attached to compratmentsmes
            self.tilling = tileSquareIngredient(self,comp,
                                         self.encapsulatingRadius,
                                         init_seed=self.histoVol.seed_used)
        elif self.packingMode == 'triangletile' :
            from autopack.hexagonTile import tileTriangleIngredient
            #self.histoVol attached to compratmentsmes
            self.tilling = tileTriangleIngredient(self,comp,
                                         self.encapsulatingRadius,
                                         init_seed=self.histoVol.seed_used)


    def DecomposeMesh(self,poly,edit=True,copy=False,tri=True,transform=True) :         
        helper = autopack.helper
#        if hasattr(m,"getFaces"):#DejaVu object
#            faces = m.getFaces()
#            vertices = m.getVertices()
#            vnormals = m.getVNormals()       
#        else :
        m=None
        if helper.host == "dejavu" :
            m = helper.getMesh(poly)
            tr=False
        else :
            m = helper.getMesh(helper.getName(poly))
            tr=True
        print ("Decompose Mesh ingredient",helper.getName(poly),m)
        #what about empty, hyerarchi, should merged all the data?
        faces,vertices,vnormals = helper.DecomposeMesh(m,
                       edit=edit,copy=copy,tri=tri,transform=tr) 
        return faces,vertices,vnormals

    def getSpheres(self,sphereFile):
        """
        get spherical approximation of shape
        """
        # file format is space separated
        # float:Rmin float:Rmax
        # int:number of levels
        # int: number of spheres in first level
        # x y z r i j k ...# first sphere in first level and 0-based indices
                           # of spheres in next level covererd by this sphere
        # ...
        # int: number of spheres in second level
        f = open(sphereFile)
        datao = f.readlines()
        f.close()
        
        # strip comments
        data = [x for x in datao if x[0]!='#' and len(x)>1]
    
        rmin, rmax = list(map(float, data[0].split()))
        nblevels = int(data[1])
        radii = []
        centers = []
        children = []
        line = 2
        for level in range(nblevels):
            rl = []
            cl = []
            ch = []
            nbs = int(data[line])
            line += 1
            for n in range(nbs):
                w = data[line].split()
                x,y,z,r = list(map(float, w[:4]))
                if level<nblevels-1: # get sub spheres indices
                    ch.append( list(map(int, w[4:])) )
                cl.append( (x,y,z) )
                rl.append( r )
                line += 1
            centers.append(cl)
            radii.append(rl)
            children.append(ch)
    
        # we ignore the hierarchy for now
        return rmin, rmax, centers, radii, children 
    

    def rejectOnce(self,rbnode,moving,afvi):
        if rbnode : 
            self.histoVol.callFunction(self.histoVol.delRB,(rbnode,))
        if afvi is not None and moving is not None :
            afvi.vi.deleteObject(moving)
        self.haveBeenRejected = True
        self.rejectionCounter += 1
        if self.rejectionCounter >= self.rejectionThreshold: #Graham set this to 6000 for figure 13b (Results Fig 3 Test1) otehrwise it fails to fill small guys
            #if verbose :
            print('PREMATURE ENDING of ingredient', self.name)
            self.completion = 1.0


    def addRBsegment(self,pt1,pt2):
        #ovewrite by grow ingredient
        pass
                        
    def Set(self,**kw):
        self.nbMol = 0   
        if "nbMol" in kw :
            nbMol = kw["nbMol"]
#            if nbMol != 0:
#                self.overwrite_nbMol = True
#                self.overwrite_nbMol_value = nbMol
#                self.nbMol = nbMol
#            else :
#                self.overwrite_nbMol =False
            self.overwrite_nbMol_value = nbMol
            #self.nbMol = nbMol
        if "molarity" in kw :
            self.molarity  = kw["molarity"]  
        if "priority" in kw :
            self.packingPriority = kw["priority"]
        if "packingMode" in kw :
            self.packingMode = kw["packingMode"]
        if "compMask" in kw :
            if type(kw["compMask"]) is str :
                self.compMask = eval(kw["compMask"])    
            else :
                self.compMask = kw["compMask"]

    def getEncapsulatingRadius(self,mesh=None):
        if self.vertices is None or not len(self.vertices) :
            if self.mesh :
                helper = autopack.helper
                if helper.host == "3dsmax" : 
                    return
                if mesh is None :
                    mesh = self.mesh
                print ("getEncapsulatingRadius ",self.mesh,mesh )
                self.faces,self.vertices,vnormals = self.DecomposeMesh(mesh,
                                   edit=True,copy=False,tri=True) 
        #print ("create the triangle",len(faces))
        #encapsulating radius ?
        v=numpy.array(self.vertices,'f')
        l=numpy.sqrt((v*v).sum(axis=1))
        r = float(max(l))+15.0
        print "self.encapsulatingRadius ",self.encapsulatingRadius,r
        self.encapsulatingRadius = r
#        if r != self.encapsulatingRadius:
#            self.encapsulatingRadius = r

    def getData(self):
        if self.vertices is None or not len(self.vertices) :
            if self.mesh :
                helper = autopack.helper
                if helper.host == "3dsmax" : 
                    return
                print ("getData ",self.mesh )
                self.faces,self.vertices,vnormals = self.DecomposeMesh(self.mesh,
                                   edit=True,copy=False,tri=True) 
        

    def rapid_model(self):
        rapid_model = RAPIDlib.RAPID_model()
        self.getData()
        if len(self.vertices) :
            rapid_model.addTriangles(numpy.array(self.vertices,'f'), numpy.array(self.faces,'i'))
        return rapid_model       
            
    def create_rapid_model(self):
        self.rapid_model = RAPIDlib.RAPID_model()
        #need triangle and vertices
        self.getData()
        if len(self.vertices) :
            self.rapid_model.addTriangles(numpy.array(self.vertices,'f'), numpy.array(self.faces,'i'))
                   
    def get_rapid_model(self):
        if (self.rapid_model is None ):
            #print ("get rapid model, create it")
            self.create_rapid_model()
            #print ("OK")
        return self.rapid_model 

    def get_rb_model(self,alt=False):
        ret = 0        
        if alt :
            ret=1
        if (self.bullet_nodes[ret] is None ):
            self.bullet_nodes[ret] = self.histoVol.addRB(self,[0.,0.,0.], numpy.identity(4),rtype=self.Type)
        return self.bullet_nodes[ret]

    def getMesh(self, filename, geomname):
        """
        Create a mesh representation from a filename for the ingredient
    
        @type  filename: string
        @param filename: the name of the input file
        @type  geomname: string
        @param geomname: the name of the ouput geometry
    
        @rtype:   DejaVu.IndexedPolygons/HostObjec
        @return:  the created mesh  
        """
        #depending the extension of the filename, can be eitherdejaVu file, fbx or wavefront
        #no extension is DejaVu
        helper = autopack.helper
        reporthook = None
        if helper is not None:        
            reporthook=helper.reporthook
#        print('TODO: getMesh need safety check for no internet connection')
#        print ("helper in Ingredient is "+str(helper))
        #should wetry to see if it already exist inthescene 
        if helper is not None and not helper.nogui:
            o = helper.getObject(geomname)
            print ("retrieve ",geomname,o)
            if o is not None :
                return o
        #identify extension
        name =   filename.split("/")[-1]
        fileName, fileExtension = os.path.splitext(name)
        print ("retrieve ",filename,fileExtension)
        if fileExtension is '' :  
            tmpFileName1 =autopack.retrieveFile(filename+".indpolface",cache="geometries")
            tmpFileName2 =autopack.retrieveFile(filename+".indpolvert",cache="geometries")
            filename = os.path.splitext(tmpFileName1)[0]
        else :
            filename =autopack.retrieveFile(filename,cache="geometries") 
        if filename is None :
            print ("problem")
            return None
        if not os.path.isfile(filename) and fileExtension != '' :
            print ("problem with "+filename,fileExtension)
            return None
        fileName, fileExtension = os.path.splitext(filename)
        print('found fileName '+fileName+' fileExtension '+fileExtension)
        if fileExtension.lower() == ".fbx" :
#            print ("read fbx withHelper",filename,helper,autopack.helper)
            #use the host helper if any to read
            if helper is not None:#neeed the helper
#                print "read "+filename
                helper.read(filename)
#                print "try to get the object "+geomname
                geom = helper.getObject(geomname)
                print ("geom ",geom,geomname,helper.getName(geom))
                #reparent to the fill parent
                if helper.host == "3dsmax" or helper.host.find("blender") != -1 :
                    helper.resetTransformation(geom)#remove rotation and scale from importing
                    #helper.rotateObj(geom,[0.0,0.0,-math.pi/2.0])
                    #m = geom.GetNodeTM()
                    #m.PreRotateY(-math.pi/2.0)
                    #geom.SetNodeTM(m)
                if helper.host != "c4d" and self.coordsystem == "left" and helper.host != "softimage":
                    #need to rotate the transform that carry the shape
                    helper.rotateObj(geom,[0.,-math.pi/2.0,0.0])
                if helper.host =="softimage" and self.coordsystem == "left" :
                    helper.rotateObj(geom,[0.0,-math.pi/2.0,0.0],primitive=True)#need to rotate the primitive                    
                if helper.host == "c4d" and self.coordsystem == "right":
                    helper.resetTransformation(geom)
                    helper.rotateObj(geom,[0.0,math.pi/2.0,math.pi/2.0],primitive=True)
#                    oldv = self.principalVector[:]
#                    self.principalVector = [oldv[2],oldv[1],oldv[0]]
                p=helper.getObject("autopackHider")
                if p is None:
                    p = helper.newEmpty("autopackHider")
                    if helper.host.find("blender") == -1 :
                        helper.toggleDisplay(p,False)
                helper.reParent(geom,p)
                return geom
            return None
        elif fileExtension == ".dae":
            print ("read dae withHelper",filename,helper,autopack.helper)
            #use the host helper if any to read
            if helper is None :
                from upy.dejavuTk.dejavuHelper import dejavuHelper
                #need to get the mesh directly. Only possible if dae or dejavu format
                #get the dejavu heper but without the View, and in nogui mode
                h =  dejavuHelper(vi="nogui")
                dgeoms = h.read(filename)
                #v,vn,f = dgeoms.values()[0]["mesh"]
                self.vertices,self.vnormals,self.faces = helper.combineDaeMeshData(dgeoms.values())
                self.vnormals=helper.normal_array(self.vertices,numpy.array(self.faces))
                geom = h.createsNmesh(geomname,self.vertices,None,self.faces)[0]
                return geom
            else : #if helper is not None:#neeed the helper
                if helper.host == "dejavu" and helper.nogui:
                    dgeoms = helper.read(filename)
                    v,vn,f = dgeoms.values()[0]["mesh"]
#                    vn = self.getVertexNormals(v,f)     
                    self.vertices,self.vnormals,self.faces = helper.combineDaeMeshData(dgeoms.values())
#                    print (self.name,len(self.vertices))
                    self.vnormals=helper.normal_array(self.vertices,numpy.array(self.faces))
                    geom = helper.createsNmesh(geomname,self.vertices,self.vnormals,self.faces)[0]
                    return geom
                else :
                    if helper.host != "dejavu":
                        coll=True
                        try :
                            import collada
                        except :
                            print ("no collada")
                            coll = False
                        if coll :
                            from upy.dejavuTk.dejavuHelper import dejavuHelper
                            #need to get the mesh directly. Only possible if dae or dejavu format
                            #get the dejavu heper but without the View, and in nogui mode
                            h =  dejavuHelper(vi="nogui")
                            dgeoms = h.read(filename)
                            #should combine both
                            self.vertices,vnormals,self.faces = h.combineDaeMeshData(dgeoms.values())#dgeoms.values()[0]["mesh"] 
                            self.vnormals=helper.normal_array(self.vertices,numpy.array(self.faces))
                helper.read(filename)
#                helper.update()
                geom = helper.getObject(geomname)
                print ("should have read...",geomname,geom)
                #rotate ?
                if helper.host == "3dsmax" :#or helper.host.find("blender") != -1:
                    helper.resetTransformation(geom)#remove rotation and scale from importing??maybe not?
                if helper.host.find("blender") != -1 :
                    helper.resetTransformation(geom)
#                    if self.coordsystem == "left" :
#                        mA = helper.rotation_matrix(-math.pi/2.0,[1.0,0.0,0.0])
#                        mB = helper.rotation_matrix(math.pi/2.0,[0.0,0.0,1.0])
#                        m=matrix(mA)*matrix(mB)
#                        helper.setObjectMatrix(geom,matrice=m)
#                if helper.host != "c4d"  and helper.host != "dejavu" and self.coordsystem == "left" and helper.host != "softimage" and helper.host.find("blender") == -1:
                    #what about softimage
                    #need to rotate the transform that carry the shape, maya ? or not ?
#                    helper.rotateObj(geom,[0.0,-math.pi/2.0,0.0])#wayfront as well euler angle
                    #swicth the axe?
#                    oldv = self.principalVector[:]
#                    self.principalVector = [oldv[2],oldv[1],oldv[0]]
                if helper.host =="softimage" and self.coordsystem == "left" :
                    helper.rotateObj(geom,[0.0,-math.pi/2.0,0.0],primitive=True)#need to rotate the primitive
                if helper.host == "c4d" and self.coordsystem == "right":
                    helper.resetTransformation(geom)
                    helper.rotateObj(geom,[0.0,math.pi/2.0,math.pi/2.0],primitive=True)
                p=helper.getObject("autopackHider")
                if p is None:
                    p = helper.newEmpty("autopackHider")
                    if helper.host.find("blender") == -1 :
                        helper.toggleDisplay(p,False)
                helper.reParent(geom,p)            
                return geom
            return None
        elif fileExtension is '' :
            return self.getDejaVuMesh(filename, geomname)
        else :#host specific file
            if helper is not None:#neeed the helper
                helper.read(filename)# doesnt get the regular file ? conver state to object
                geom = helper.getObject(geomname)
                print ("should have read...",geomname,geom)
                p=helper.getObject("autopackHider")
                if p is None:
                    p = helper.newEmpty("autopackHider")
                    if helper.host.find("blender") == -1 :helper.toggleDisplay(p,False)
                helper.reParent(geom,p)            
                return geom
            return None
            
        
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
#        print ("dejavu mesh", filename)        
        geom = IndexedPolygonsFromFile(filename, 'mesh_%s'%self.pdb)
#        if helper is not None:
#            if helper.host != "maya" :
#                helper.rotateObj(geom,[0.0,-math.pi/2.0,0.0])
        return geom


    def jitterPosition(self, position, spacing, normal = None):
        """
        position are the 3d coordiantes of the grid point
        spacing is the grid spacing
        this will jitter gauss(0., 0.3) * Ingredient.jitterMax
        """
        if self.compNum > 0:
            vx, vy, vz = v1 = self.principalVector
            #surfacePointsNormals problem here
            v2 = normal
            try :
                rotMat = numpy.array( rotVectToVect(v1, v2 ), 'f')
            except :
                print('PROBLEM ', self.name)
                rotMat = numpy.identity(4)
            
        jx, jy, jz = self.jitterMax
        dx = jx*spacing*uniform(-1.0, 1.0)  #This needs to use the same rejection if outside of the sphere that the uniform cartesian jitters have.  Shoiuld use oneJitter instead?
        dy = jy*spacing*uniform(-1.0, 1.0)
        dz = jz*spacing*uniform(-1.0, 1.0)
#        d2 = dx*dx + dy*dy + dz*dz
#        if d2 < jitter2:
        if self.compNum > 0: # jitter less among normal
            dx, dy, dz, dum = numpy.dot(rotMat, (dx,dy,dz,0))
        position[0] += dx
        position[1] += dy
        position[2] += dz
        return position


    def getMaxJitter(self, spacing):
        return max(self.jitterMax)*spacing

##     def checkCollisions(self, pt, radius, pointsInCube, gridPointsCoords,
##                         distance, verbose):
##         insidePoints = {}
##         newDistPoints = {}
##         x, y, z = pt
##         for pt in pointsInCube:
##             x1,y1,z1 = gridPointsCoords[pt]
##             dist = sqrt((x1-x)*(x1-x) + (y1-y)*(y1-y) + (z1-z)*(z1-z))
##             d = dist-radius
##             if dist < radius:  # point is inside dropped sphere
##                 if distance[pt]<-0.0001:
##                     return None, None, pt
##                 else:
##                     if insidePoints.has_key(pt):
##                         if d < insidePoints[pt]:
##                             insidePoints[pt] = d
##                     else:
##                         insidePoints[pt] = d

##             else: # update distance is smaller that existing distance
##                 if d < distance[pt]:
##                     if newDistPoints.has_key(pt):
##                         if d < newDistPoints[pt]:
##                             newDistPoints[pt] = d
##                     else:
##                         newDistPoints[pt] = d

##                     if verbose > 5:
##                         print 'point ',pt, 'going from ', \
##                               distance[pt],'to', dist

##         return insidePoints, newDistPoints, None

    def swap(self,d, n):
        d.rotate(-n)
        d.popleft()
        d.rotate(n)

    def deleteblist(self,d, n):
        del d[n]

    def getDistancesCube(self,jtrans, rotMat,gridPointsCoords, distance, grid):
        radii = self.radii        
        insidePoints = {}
        newDistPoints = {}
        cent1T = self.transformPoints(jtrans, rotMat, self.positions[0])[0]#bb1
        cent2T = self.transformPoints(jtrans, rotMat, self.positions2[0])[0]#bb2
        center = self.transformPoints(jtrans, rotMat, [self.center,])[0]
#        cylNum = 0
#        for radc, p1, p2 in zip(radii, cent1T, cent2T):
        x1, y1, z1 = cent1T
        x2, y2, z2 = cent2T
        vx, vy, vz =  (x2-x1, y2-y1, z2-z1)
        lengthsq = vx*vx + vy*vy + vz*vz
        l = sqrt( lengthsq )
        cx, cy, cz = posc = center#x1+vx*.5, y1+vy*.5, z1+vz*.5
        radt = l/2. + self.encapsulatingRadius
        
        bb = [cent2T,cent1T]#self.correctBB(p1,p2,radc)
        x,y,z = posc
        bb = ( [x-radt, y-radt, z-radt], [x+radt, y+radt, z+radt] )
#        print ("pointsInCube",bb,posc,radt)        
        pointsInGridCube = grid.getPointsInCube(bb, posc, radt)
        
        # check for collisions with cylinder    
        pd = numpy.take(gridPointsCoords,pointsInGridCube,0)-center
        
        delta = pd.copy()        
        delta *= delta
        distA = numpy.sqrt( delta.sum(1) )
        
        m = numpy.matrix(numpy.array(rotMat).reshape(4,4))#
        mat = m.I
        #need to apply inverse mat to pd
        rpd = ApplyMatrix(pd,mat)
        #need to check if these point are inside the cube using the dimension of the cube
        #numpy.fabs        
        res = numpy.less_equal(numpy.fabs(rpd),numpy.array(radii[0])/2.)
        if len(res) :
            c=numpy.average(res,1)#.astype(int)
            d=numpy.equal(c,1.)
            ptinsideCube = numpy.nonzero(d)[0]
        else :
            ptinsideCube = []
        for pti in range(len(pointsInGridCube)):#ptinsideCube:#inside point but have been already computed during the check collision...?
            pt = pointsInGridCube[pti]
            if pt in insidePoints: continue
            dist = distA[pti]
            d = dist-self.encapsulatingRadius #should be distance to the cube, but will use approximation
            if pti in ptinsideCube:# dist < radt:  # point is inside dropped sphere
                if pt in insidePoints:
                    if d < insidePoints[pt]:
                        insidePoints[pt] = d
                else:
                    insidePoints[pt] = d
            elif d < distance[pt]: # point in region of influence
                if pt in newDistPoints:
                    if d < newDistPoints[pt]:
                        newDistPoints[pt] = d
                else:
                    newDistPoints[pt] = d
        return insidePoints,newDistPoints

    def getDistancesSphere(self,jtrans, rotMatj,gridPointsCoords, distance,level,dpad):
        self.centT = centT = self.transformPoints(jtrans, rotMatj, self.positions[level])
        centT = self.centT#self.transformPoints(jtrans, rotMatj, self.positions[-1])
        insidePoints = {}
        newDistPoints = {}
        for radc, posc in zip(self.radii[-1], centT):
            rad = radc + dpad
            x,y,z = posc
            #this have already be done in the checkCollision why doing it again
            bb = ( [x-rad, y-rad, z-rad], [x+rad, y+rad, z+rad] )
            pointsInCube = self.histoVol.callFunction(self.histoVol.grid.getPointsInCube,
                                                 (bb, posc, rad))
#
            delta = numpy.take(gridPointsCoords,pointsInCube,0)-posc
            delta *= delta
            distA = numpy.sqrt( delta.sum(1) )
            ptsInSphere = numpy.nonzero(numpy.less_equal(distA, rad))[0]
#            ptsInSphere =self.histoVol.grid.getPointsInSphere(self, bb, pt, radius,addSP=True,info=False)
            for pti in ptsInSphere:
                pt = pointsInCube[pti]
                if pt in insidePoints: continue
                dist = distA[pti]
                d = dist-radc
                if dist < radc:  # point is inside dropped sphere
                    if pt in insidePoints:
                        if d < insidePoints[pt]:
                            insidePoints[pt] = d
                    else:
                        insidePoints[pt] = d
                elif d < distance[pt]: # point in region of influence
                    if pt in newDistPoints:
                        if d < newDistPoints[pt]:
                            newDistPoints[pt] = d
                    else:
                        newDistPoints[pt] = d
        return insidePoints,newDistPoints

    def getDistancesCylinders(self,jtrans, rotMatj,gridPointsCoords, distance,dpad):
        insidePoints = {}
        newDistPoints = {}
        cent1T = self.transformPoints(jtrans, rotMatj, self.positions[-1])
        cent2T = self.transformPoints(jtrans, rotMatj, self.positions2[-1])

        for radc, p1, p2 in zip(self.radii[-1], cent1T, cent2T):

            x1, y1, z1 = p1
            x2, y2, z2 = p2
            vx, vy, vz = vect = (x2-x1, y2-y1, z2-z1)
            lengthsq = vx*vx + vy*vy + vz*vz
            l = sqrt( lengthsq )
            cx, cy, cz = posc = x1+vx*.5, y1+vy*.5, z1+vz*.5
            radt = l + radc + dpad
            bb = ( [cx-radt, cy-radt, cz-radt], [cx+radt, cy+radt, cz+radt] )
            pointsInCube = self.histoVol.callFunction(self.histoVol.grid.getPointsInCube,
                                                 (bb, posc, radt))

            pd = numpy.take(gridPointsCoords,pointsInCube,0) - p1
            dotp = numpy.dot(pd, vect)
            rad2 = radc*radc
            d2toP1 = numpy.sum(pd*pd, 1)
            dsq = d2toP1 - dotp*dotp/lengthsq

            pd2 = numpy.take(gridPointsCoords,pointsInCube,0) - p2
            d2toP2 = numpy.sum(pd2*pd2, 1)

            for pti, pt in enumerate(pointsInCube):
                if pt in insidePoints: continue

                if dotp[pti] < 0.0: # outside 1st cap
                    d = sqrt(d2toP1[pti])
                    if d < distance[pt]: # point in region of influence
                        if pt in newDistPoints:
                            if d < newDistPoints[pt]:
                                newDistPoints[pt] = d
                        else:
                            newDistPoints[pt] = d
                elif dotp[pti] > lengthsq:
                    d = sqrt(d2toP2[pti])
                    if d < distance[pt]: # point in region of influence
                        if pt in newDistPoints:
                            if d < newDistPoints[pt]:
                                newDistPoints[pt] = d
                        else:
                            newDistPoints[pt] = d
                else:
                    d = sqrt(dsq[pti])-radc
                    if d < 0.:  # point is inside dropped sphere
                        if pt in insidePoints:
                            if d < insidePoints[pt]:
                                insidePoints[pt] = d
                        else:
                            insidePoints[pt] = d
        return insidePoints,newDistPoints

    def getDistances(self,jtrans, rotMatj,gridPointsCoords, distance,dpad):
        insidePoints = {}
        newDistPoints = {}
        if self.modelType=='Spheres':
            insidePoints,newDistPoints=self.getDistancesSphere(jtrans, rotMatj,
                           gridPointsCoords, distance, self.collisionLevel,dpad)
        elif self.modelType=='Cylinders':
            insidePoints,newDistPoints=self.getDistancesCylinders(jtrans, rotMatj,
                                      gridPointsCoords, distance,dpad)
        elif self.modelType=='Cube':
            insidePoints,newDistPoints=self.getDistancesCube(jtrans, rotMatj,
                             gridPointsCoords, distance, self.histoVol.grid)
        return insidePoints,newDistPoints   
                
    def updateDistances(self, histoVol,insidePoints, newDistPoints, freePoints,
                        nbFreePoints, distance, masterGridPositions, verbose):
#        print("*************updating Distances")
        verbose = histoVol.verbose
        t1 = time()
        distChanges = {}
        self.nbPts = len(insidePoints)
#        print("nbPts = len(insidePoints) = ", self.nbPts)
#        print("nbFreePoints = ", nbFreePoints)
#        print("lenFreePoints = ", len(freePoints))
#        fptCount = 0
#        for val1 in freePoints:
#            print("FreePoint[",fptCount,"] = ", freePoints[fptCount]," = val1 = ", val1)
#            fptCount += 1
        for pt,dist in list(insidePoints.items()):  #Reversing is not necessary if you use the correct Swapping GJ Aug 17,2012
            #        for pt,dist in reversed(list(insidePoints.items()) ):  # GJ notes (August 17, 2012): Critical to reverse 
            # or points that need to get masked towards the end get used incorrectly as valid points during forward swap
            # Reversing the walk through the masked points cures this! 
            # swap reverse point at ptIndr with last free one
            # pt is the grid point indice not the freePt indice
#            try :
#                fi = freePoints.index(pt)#too slow
#            except :
#                pass
            try :
                # New system replaced by Graham on Aug 18, 2012
                nbFreePoints -= 1  
                vKill = freePoints[pt]
                vLastFree = freePoints[nbFreePoints]
                freePoints[vKill] = vLastFree
                freePoints[vLastFree] = vKill
                # End New replaced by Graham on Aug 18, 2012
            # Start OLD system replaced by Graham on Aug 18, 2012. This has subtle problems of improper swapping
            #                tmp = freePoints[nbFreePoints-1] #last one
            #                #freePoints.remove(pt)
            #                freePoints[nbFreePoints-1] = pt
            #                tmpDebug2 = freePoints[freePoints[pt]]
            ##                freePoints[nbFreePoints-1] = freePoints[pt]
            #                freePoints[freePoints[pt]] = tmp
            #                nbFreePoints -= 1
            # End OLD system replaced by Graham on Aug 18, 2012. This has subtle problems of improper swapping
            except :
#                print (pt, "not in freeePoints********************************")
                pass 
    
        #Turn on these printlines if there is a problem with incorrect points showing in display points        
#            print("*************pt = masterGridPointValue = ", pt)
#            print("nbFreePointAfter = ", nbFreePoints)    
#            print("vKill = ", vKill)
#            print("vLastFree = ", vLastFree)
#            print("freePoints[vKill] = ", freePoints[vKill])
#            print("freePoints[vLastFree] = ", freePoints[vLastFree])
#            print("pt = masterGridPointValue = ", pt)
#            print("freePoints[nbFreePoints-1] = ", freePoints[nbFreePoints])
#            print("freePoints[pt] = ", freePoints[pt])
                
            
            distChanges[pt] = (masterGridPositions[pt],
                               distance[pt], dist)
            distance[pt] = dist
#            print("distance[pt] = ", distance[pt])
#            print("distChanges[pt] = ", distChanges[pt])

                #self.updateDistanceForOther(pt=pt)
        if verbose >4:       
            print("update freepoints loop",time()-t1)
        if verbose >5:
            print(nbFreePoints)
        t2 = time()
        for pt,dist in list(newDistPoints.items()):
            if pt not in insidePoints:
                distChanges[pt] = (masterGridPositions[pt],
                                   distance[pt], dist)
                distance[pt] = dist                    
            #self.updateDistanceForOther(pt=pt,dist=dist)
        #hack c4d particle
        if verbose==0.4: 
            print("update distance loop",time()-t2)
        #bitmaps.ShowBitmap(bmp)
        #bmp.Save(name,c4d.FILTER_TIF)
        
            
        return nbFreePoints
        
#   DEPRECRATED FUNCTION
#    def updateDistanceForOther(self,pt=None,dist=None):
#        #return
#        r = self.recipe()#stored_recipe it is  weakref
#        for ingr in r.ingredients:
#            if hasattr(ingr,"allIngrPts") and ingr.updateOwnFreePts:
#                if dist is not None and dist  < ingr.cut :
#                    try :
#                        n = histoVol.lmethod.swap_value(ingr.allIngrPts,pt)
#                        ingr.NI = n
#                        #ingr.allIngrPts.remove(pt)
#                        #del ingr.allIngrPts[pt]
#                    except :
#                        pass
#                else :
#                    try :
#                        n = histoVol.lmethod.swap_value(ingr.allIngrPts,pt)
#                        ingr.NI = n
#                        #ingr.allIngrPts.remove(pt)
#                        #del ingr.allIngrPts[pt]
#                    except :
#                        pass
                

    def perturbAxis(self, amplitude):
        # modify axis using gaussian distribution but clamp
        # at amplitutde
        x,y,z = self.principalVector
        stddev = amplitude *.5
        dx = gauss(0., stddev)
        if dx>amplitude: dx = amplitude
        elif dx<-amplitude: dx = -amplitude
        dy = gauss(0., stddev)
        if dy>amplitude: dy = amplitude
        elif dy<-amplitude: dy = -amplitude
        dz = gauss(0., stddev)
        if dz>amplitude: dz = amplitude
        elif dz<-amplitude: dz = -amplitude
        #if self.name=='2bg9 ION CHANNEL/RECEPTOR':
        #    print 'FFFFFFFFFFFFF AXIS', x+dx,y+dy,z+dz
        return (x+dx,y+dy,z+dz)


    def transformPoints(self, trans, rot, points):
#        if helper is not None :
#            rot[:3,:3] = trans
#            return helper.ApplyMatrix(points,rot)
        tx,ty,tz = trans
        pos = []
        for xs,ys,zs in points:
            x = rot[0][0]*xs + rot[0][1]*ys + rot[0][2]*zs + tx
            y = rot[1][0]*xs + rot[1][1]*ys + rot[1][2]*zs + ty
            z = rot[2][0]*xs + rot[2][1]*ys + rot[2][2]*zs + tz
            pos.append( [x,y,z] )
#        print ("transform to",points,rot,trans, pos)
#        print (autopack.helper.ApplyMatrix(points,rot)+numpy.array(trans))
        return numpy.array(pos)

    def alignRotation(self,jtrans):
        # for surface points we compute the rotation which
        # aligns the principalVector with the surface normal
        if self.compNum == 0 :
            compartment = self.histoVol
        else :
            compartment = self.histoVol.compartments[abs(self.compNum)-1]
        vx, vy, vz = v1 = self.principalVector
        #surfacePointsNormals problem here
        gradient_center = self.histoVol.gradients[self.gradient].direction
        v2 = numpy.array(gradient_center) - numpy.array(jtrans)
        try :
            rotMat = numpy.array( rotVectToVect(v1, v2 ), 'f')
        except :
            print('PROBLEM ', self.name)
            rotMat = numpy.identity(4)
        return rotMat

    def getAxisRotation(self, rot):
        """
        combines a rotation about axis to incoming rot.
        rot aligns the principalVector with the surface normal
        rot aligns the principalVector with the biased diretion
        """
        if self.perturbAxisAmplitude!=0.0:
            axis = self.perturbAxis(self.perturbAxisAmplitude)
        else:
            axis = self.principalVector
        tau = uniform(-pi, pi)
        rrot = rotax( (0,0,0), axis, tau, transpose=1 )
        rot = numpy.dot(rot, rrot)
        return rot

    def getBiasedRotation(self, rot,weight=None):
        """
        combines a rotation about axis to incoming rot
        """
        if self.perturbAxisAmplitude!=0.0:
            axis = self.perturbAxis(self.perturbAxisAmplitude)
        else:
            axis = self.rotAxis
        #-30,+30 ?
        if weight is not None :
           tau = uniform(-pi*weight, pi*weight)#(-pi, pi)
        else : 
            tau = gauss(self.orientBiasRotRangeMin, self.orientBiasRotRangeMax)#(-pi, pi)
        rrot = rotax( (0,0,0), self.rotAxis, tau, transpose=1 )
        rot = numpy.dot(rot, rrot)
        return rot

    def correctBB(self,p1,p2,radc):
        #unprecised
        x1,y1,z1=p1
        x2,y2,z2=p2
#        bb = ( [x1-radc, y1-radc, z1-radc], [x2+radc, y2+radc, z2+radc] )
        mini=[]
        maxi=[]
        for i in range(3):
            mini.append(min(p1[i],p2[i])-radc)
            maxi.append(max(p1[i],p2[i])+radc)
        return numpy.array([numpy.array(mini).flatten(),numpy.array(maxi).flatten()])
        #precised:
#        kx=sqrt(((A.Y-B.Y)^2+(A.Z-B.Z)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))
#        ky=sqrt(((A.X-B.X)^2+(A.Z-B.Z)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))
#        kz=sqrt(((A.X-B.X)^2+(A.Y-B.Y)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))

    def checkDistSurface(self,point,cutoff):
        if not hasattr(self,"histoVol") :
            return False
        if self.compNum == 0 :
            compartment = self.histoVol
        else :
            compartment = self.histoVol.compartments[abs(self.compNum)-1]
        compNum = self.compNum
#        print "compNum ",compNum
        if compNum < 0 :
            sfpts = compartment.surfacePointsCoords
            delta = numpy.array(sfpts)-numpy.array(point)
            delta *= delta
            distA = numpy.sqrt( delta.sum(1) )
#            print len(distA)
            test = distA < cutoff
            if True in test:
                return True
        elif compNum == 0 :
            for o in self.histoVol.compartments:
                sfpts = o.surfacePointsCoords
                delta = numpy.array(sfpts)-numpy.array(point)
                delta *= delta
                distA = numpy.sqrt( delta.sum(1) )
#                print len(distA)
                test = distA < cutoff
                if True in test:
                    return True
        return False
        
    def getListCompFromMask(self,cId,ptsInSphere):
        #cID ie [-2,-1,-2,0...], ptsinsph = [519,300,etc]
        current = self.compNum
        if current < 0 : #inside
            mask = ["self"] #authorize in and surf
            ins=[i for i,x in enumerate(cId) if x == current]
            #surf=[i for i,x in enumerate(cId) if x == -current]
            liste = ins#+surf
        if current > 0 :#surface
            mask = ["self","neg"] #authorize in and surf and extra but not ther compartment
            ins=[i for i,x in enumerate(cId) if x == current]
            surf=[i for i,x in enumerate(cId) if x == -current]
            extra=[i for i,x in enumerate(cId) if x < 0]
            liste = ins+surf+extra         
        elif current == 0 :#extracellular
            mask = ["self"]
            liste=[i for i,x in enumerate(cId) if x == current]
        return liste

    def isInGoodComp(self,pId,nbs=None):
        #cID ie [-2,-1,-2,0...], ptsinsph = [519,300,etc]
        current = self.compNum
        cId = self.histoVol.grid.gridPtId[pId]
        if current <= 0 : #inside
            if current != cId :
                return False
            return True
        if current > 0 :#surface
            if current != cId and -current != cId :  
                return False
            return True
        return False

    def compareCompartmentPrimitive(self,level,jtrans, rotMatj, 
                                    gridPointsCoords, distance):
        collisionComp = False                                    
        if self.modelType=='Spheres':
            collisionComp = self.histoVol.callFunction(self.checkSphCompart,(
                self.positions[level], self.radii[level], jtrans, rotMatj,
                level, gridPointsCoords, distance, self.histoVol))
        elif self.modelType=='Cylinders':
            collisionComp = self.histoVol.callFunction(self.checkCylCompart,(
                self.positions[level], self.positions2[level],
                self.radii[level], jtrans, rotMatj, gridPointsCoords,
                distance, histoVol))
        elif self.modelType=='Cube':
            collisionComp = self.histoVol.callFunction(self.checkCubeCompart,(
                self.positions[0], self.positions2[0], self.radii,
                jtrans, rotMatj, gridPointsCoords,
                distance, histoVol))
        return collisionComp
        
    def checkCompartmentAlternative(self,ptsId,histoVol,nbs=None):
        compIds = numpy.take(histoVol.grid.gridPtId,ptsId,0)  
#        print "compId in listPtId",compIds
        if self.compNum <= 0 :
            wrongPt = [ cid for cid in compIds if cid != self.compNum ]
            if len(wrongPt):
#                print wrongPt
                return True
        return False
        
    def checkCompartment(self,ptsInSphere,nbs=None):
        trigger = False
#        print ("checkCompartment using",len(ptsInSphere))
#        print (ptsInSphere)
        if self.compareCompartment:
            cId = numpy.take(self.histoVol.grid.gridPtId,ptsInSphere,0)#shoud be the same ?
            if nbs != None:
                #print ("cId ",cId,ptsInSphere)
                if self.compNum <= 0 and nbs != 0 :
                    return trigger,True                               
            L = self.getListCompFromMask(cId,ptsInSphere)
            
            #print ("liste",L)
            if len(cId) <= 1 :
                return trigger,True
            p = float(len(L))/float(len(cId))#ratio accepted compId / totalCompId-> want 1.0
            if p < self.compareCompartmentTolerance:
                #print ("the ratio is ",p, " threshold is ",self.compareCompartmentThreshold," and tolerance is ",self.compareCompartmentTolerance)
                trigger = True
                return trigger,True
            #threshold
            if self.compareCompartmentThreshold != 0.0 and \
                p < self.compareCompartmentThreshold:
                    return trigger,True
                        #reject the ingr
        return trigger,False

    def checkCylCollisions(self, centers1, centers2, radii, jtrans, rotMat,
                           gridPointsCoords, distance, histoVol):
        """
        Check cylinders for collision
        """
#        print "#######################"
#        print jtrans
#        print rotMat
        cent1T = self.transformPoints(jtrans, rotMat, centers1)
        cent2T = self.transformPoints(jtrans, rotMat, centers2)

        cylNum = 0
        for radc, p1, p2 in zip(radii, cent1T, cent2T):
            if histoVol.runTimeDisplay > 1:
                name = "cyl"
                cyl = self.vi.getObject("cyl")
                if cyl is None:
                    cyl=self.vi.oneCylinder(name,p1,p2,
                                            color=(1.,1.,1.),
                                            radius=radc)
#                    self.vi.updateTubeMesh(cyl,cradius=radc)
                else :
                    self.vi.updateOneCylinder(cyl,p1,p2,radius=radc)
                self.vi.changeObjColorMat(cyl,(1.,1.,1.))
                name = "sph1"
                sph1 = self.vi.getObject("sph1")
                if sph1 is None:
                    sph1=self.vi.Sphere(name,radius=radc*2.)[0]
                self.vi.setTranslation(sph1,p1)
                name = "sph2"
                sph2 = self.vi.getObject("sph2")
                if sph2 is None:
                    sph2=self.vi.Sphere(name,radius=radc*2.)[0]
                self.vi.setTranslation(sph2,p2)
 
                self.vi.update()
            x1, y1, z1 = p1
            x2, y2, z2 = p2
            vx, vy, vz = vect = (x2-x1, y2-y1, z2-z1)
            lengthsq = vx*vx + vy*vy + vz*vz
            l = sqrt( lengthsq )
            cx, cy, cz = posc = x1+vx*.5, y1+vy*.5, z1+vz*.5
            radt = l + radc
            
            bb = self.correctBB(p1,p2,radc)
#            bb = self.correctBB(posc,posc,radt)
            if histoVol.runTimeDisplay > 1:
                box = self.vi.getObject("collBox")
                if box is None:
                    box = self.vi.Box('collBox', cornerPoints=bb,visible=1)
                else :
#                    self.vi.toggleDisplay(box,True)
                    self.vi.updateBox(box,cornerPoints=bb)
                    self.vi.update()
#                 sleep(1.0)
            pointsInCube = histoVol.grid.getPointsInCube(bb, posc, radt,info=True)
            
            # check for collisions with cylinder            
            pd = numpy.take(gridPointsCoords,pointsInCube,0)-p1
            dotp = numpy.dot(pd, vect)
            rad2 = radc*radc
            dsq = numpy.sum(pd*pd, 1) - dotp*dotp/lengthsq

            ptsWithinCaps = numpy.nonzero( numpy.logical_and(
               numpy.greater_equal(dotp, 0.), numpy.less_equal(dotp, lengthsq)))
#            if not len(ptsWithinCaps[0]):
#                print "no point inside the geom?"
#                return False
            if self.compareCompartment:
                ptsInSphereId = numpy.take(pointsInCube,ptsWithinCaps[0],0)
                compIdsSphere = numpy.take(histoVol.grid.gridPtId,ptsInSphereId,0)  
#                print "compId",compIdsSphere
                if self.compNum <= 0 :
                    wrongPt = [ cid for cid in compIdsSphere if cid != self.compNum ]
                    if len(wrongPt):
#                        print wrongPt
                        return True                
#            trigger, res = self.checkCompartment(numpy.take(pointsInCube,ptsWithinCaps[0],0),nbs=nbs)
#            print ("checkCompartment result",trigger, res)
#            if res :
                #reject
#                return True
            
            for pti in ptsWithinCaps[0]:
                pt = pointsInCube[pti]
                dist = dsq[pti]
                if dist > rad2: continue # outside radius
                elif distance[pt]<-0.0001 :#or trigger: # pt is inside cylinder
                    #changeObjColorMat
                    if histoVol.runTimeDisplay > 1:
                        self.vi.changeObjColorMat(cyl,(1.,0.,0.))
                        self.vi.update()
#                        sleep(1.0)
                    #reject
                    return True
            cylNum += 1
        return False

    def checkCylCompart(self, centers1, centers2, radii, jtrans, rotMat,
                           gridPointsCoords, distance, histoVol):
        """
        Check cylinders for collision
        """
#        print "#######################"
#        print jtrans
#        print rotMat
        cent1T = self.transformPoints(jtrans, rotMat, centers1)
        cent2T = self.transformPoints(jtrans, rotMat, centers2)

        cylNum = 0
        for radc, p1, p2 in zip(radii, cent1T, cent2T):
            x1, y1, z1 = p1
            x2, y2, z2 = p2
            vx, vy, vz = vect = (x2-x1, y2-y1, z2-z1)
            lengthsq = vx*vx + vy*vy + vz*vz
            l = sqrt( lengthsq )
            cx, cy, cz = posc = x1+vx*.5, y1+vy*.5, z1+vz*.5
            radt = l + radc
            
            bb = self.correctBB(p1,p2,radc)
            pointsInCube = histoVol.grid.getPointsInCube(bb, posc, radt,info=True)
            
            # check for collisions with cylinder            
            pd = numpy.take(gridPointsCoords,pointsInCube,0)-p1
            dotp = numpy.dot(pd, vect)
#            rad2 = radc*radc
#            dsq = numpy.sum(pd*pd, 1) - dotp*dotp/lengthsq
            ptsWithinCaps = numpy.nonzero( numpy.logical_and(
               numpy.greater_equal(dotp, 0.), numpy.less_equal(dotp, lengthsq)))

            ptsInSphereId = numpy.take(pointsInCube,ptsWithinCaps[0],0)
            compIdsSphere = numpy.take(histoVol.grid.gridPtId,ptsInSphereId,0)  
            if self.compNum <= 0 :
                    wrongPt = [ cid for cid in compIdsSphere if cid != self.compNum ]
                    if len(wrongPt):
#                        print wrongPt
                        return True
            cylNum += 1
        return False
        
    def checkSphCollisions(self, centers, radii, jtrans, rotMat, level,
                        gridPointsCoords, distance, histoVol):
        """
        Check spheres for collision
        """
        #wouldnt be faster to do sphere-sphere distance test ? than points/points from the grid
        self.centT = centT = self.transformPoints(jtrans, rotMat, centers)#this should be jtrans
#        print "sphCollision",centT,radii
        sphNum = 0
        self.distances_temp=[]
        if self.compareCompartment:
            listeCpmNum=[]
        for radc, posc in zip(radii, centT):
            r=[]
            x,y,z = posc
            bb = ( [x-radc, y-radc, z-radc], [x+radc, y+radc, z+radc] )
#            if histoVol.runTimeDisplay:
#                box = self.vi.getObject("collBox")
#                if box is None:
#                    box = self.vi.Box('collBox', cornerPoints=bb,visible=1)
#                else :
##                    self.vi.toggleDisplay(box,True)
#                    self.vi.updateBox(box,cornerPoints=bb)
#                    self.vi.update()
            pointsInCube = histoVol.grid.getPointsInCube(bb, posc, radc,info=True)#indices
            r.append(pointsInCube)
            
#            print("boundingBox forPointsInCube = ", bb)
#            print "boudnig",bb,len(pointsInCube)
            # check for collisions
            delta = numpy.take(gridPointsCoords,pointsInCube,0)-posc
            delta *= delta
            distA = numpy.sqrt( delta.sum(1) )
            ptsInSphere = numpy.nonzero(numpy.less_equal(distA, radc))[0]
            r.append(ptsInSphere)
            r.append(distA)
            self.distances_temp.append(r)

            ptsInSphereId = numpy.take(pointsInCube,ptsInSphere,0)
            if self.compareCompartment:
                compIdsSphere = numpy.take(histoVol.grid.gridPtId,ptsInSphereId,0)  
#                print "compId in sphere",compIdsSphere
                wrongPt = [ cid for cid in compIdsSphere if cid == 99999 ]                
                if len(wrongPt):
                    return True
                if self.compNum <= 0 :
                    wrongPt = [ cid for cid in compIdsSphere if cid != self.compNum ]
                    if len(wrongPt):
#                        print wrongPt
                        return True
#            trigger, res = self.checkCompartment(ptsInSphereId,nbs=nbs)            
##            print ("checkCompartment result trigger and res",trigger, res)
#            if res :
#                return True
            for pti in ptsInSphere:
                pt = pointsInCube[pti]
                dist = distA[pti]
                d = dist-radc
                #print dist,d,radc,distance[pt]
#                if dist < radc:  # point is inside dropped sphere
#                if self.compareCompartment:
#                    if not self.isInGoodComp(pt):
#                        return True
                if distance[pt]<-0.0001:# or trigger:#trigger mean different compId
#                    print 'Col level:%d  d:%.1f  distance:%.1f'%(level, d, distance[pt]),
                    #return True
                    if level < self.maxLevel:
                        nxtLevelSpheres = self.positions[level+1]
                        nxtLevelRadii = self.radii[level+1]
                        # get sphere that are children of this one
                        ccenters = []
                        cradii = []
                        for sphInd in self.children[level][sphNum]:
                            ccenters.append( nxtLevelSpheres[sphInd] )
                            cradii.append( nxtLevelRadii[sphInd] )
                        collision = self.checkSphCollisions(
                            ccenters, cradii, jtrans, rotMat, level+1,
                            gridPointsCoords, distance, histoVol)
                        if not collision and level>0:
                            #import pdb
                            #pdb.set_trace()
                            #print("returning notCollision and level>0 as Collision")
                            return collision
                        #print("returning regular collision")
                        return collision
                    else:
                        #print("in Collision, but returning True")
                        return True
                    # FIXME DEBUG INFO
                    if d+distance[pt] < histoVol.maxColl:
                        histoVol.maxColl = d+distance[pt]
                            #print("in collision histovol.maxColl if")
                    return True
                        #print("End of collision for pt = ", pt)
            sphNum += 1
                #print("collision returning False")
        
        return False

    def checkSphCompart(self, centers, radii, jtrans, rotMat, level,
                        gridPointsCoords, distance, histoVol):
        """
        Check spheres for collision
        TODO improve the testwhen grid stepSize is larger that size of the ingredient 
        """
        print ("OK sphere compartment checking",self.compNum)
        centT = self.transformPoints(jtrans, rotMat, centers)#this should be jtrans
#        print "sphCollision",centT,radii
        sphNum = 0
#        self.distances_temp=[]
#        if self.compareCompartment:
#            listeCpmNum=[]
        for radc, posc in zip(radii, centT):
#            r=[]
            x,y,z = posc
            bb = ( [x-radc, y-radc, z-radc], [x+radc, y+radc, z+radc] )
            pointsInCube = histoVol.grid.getPointsInCube(bb, posc, radc,info=True)#indices
#            r.append(pointsInCube)
            
            delta = numpy.take(gridPointsCoords,pointsInCube,0)-posc
            delta *= delta
            distA = numpy.sqrt( delta.sum(1) )
            ptsInSphere = numpy.nonzero(numpy.less_equal(distA, radc))[0]
            ptsInSphereId = numpy.take(pointsInCube,ptsInSphere,0)
            compIdsSphere = numpy.take(histoVol.grid.gridPtId,ptsInSphereId,0)  
            print (len(compIdsSphere),compIdsSphere)
            if self.compNum <= 0 :
                wrongPt = [ cid for cid in compIdsSphere if cid != self.compNum ]
                if len(wrongPt):
                    print ("OK false compartment",len(wrongPt))
                    return True
        return False

    def checkCubeCollisions(self, centers1, centers2, radii, jtrans, rotMat,
                           gridPointsCoords, distance, histoVol):
        """
        Check cube for collision
        centers1 and centers2 should be the cornerPoints ?
        can also use the center plus size (radii), or the position/position2
        """
        cent1T = self.transformPoints(jtrans, rotMat, centers1)[0]#bb1
        cent2T = self.transformPoints(jtrans, rotMat, centers2)[0]#bb2
        center = self.transformPoints(jtrans, rotMat, [self.center,])[0]
        
        cylNum = 0
#        for radc, p1, p2 in zip(radii, cent1T, cent2T):
        x1, y1, z1 = cent1T
        x2, y2, z2 = cent2T
        vx, vy, vz =  (x2-x1, y2-y1, z2-z1)
        lengthsq = vx*vx + vy*vy + vz*vz
        l = sqrt( lengthsq )
        cx, cy, cz = posc = center#x1+vx*.5, y1+vy*.5, z1+vz*.5
        radt = l/2. + self.encapsulatingRadius
        x,y,z = posc
        bb = ( [x-radt, y-radt, z-radt], [x+radt, y+radt, z+radt] )
        
#        bb = [cent2T,cent1T]#self.correctBB(p1,p2,radc)
#            bb = self.correctBB(posc,posc,radt)
        if histoVol.runTimeDisplay :#> 1:
            print ("collBox",bb)
            box = self.vi.getObject("collBox")
            if box is None:
                box = self.vi.Box('collBox', 
                cornerPoints=bb,
#                center=center, 
#                size = [radt,radt,radt],
                visible=1)# cornerPoints=bb,visible=1)
            else :
#                    self.vi.toggleDisplay(box,True)
                self.vi.updateBox(box,
                cornerPoints=bb,
#                center=center, 
#                size = [radt,radt,radt],
                )#cornerPoints=bb)
            self.vi.update()
#                 sleep(1.0)
#        print ("pointsInCube",bb,posc,radt)        
        pointsInCube = histoVol.grid.getPointsInCube(bb, posc, radt)
        
        # check for collisions with cylinder    
#        print   ("pointsInCube",pointsInCube)      
        pd = numpy.take(gridPointsCoords,pointsInCube,0)-center
#        print ("Cube",rotMat)
        m = numpy.matrix(numpy.array(rotMat).reshape(4,4))#
        mat = m.I
        #need to apply inverse mat to pd
        rpd = ApplyMatrix(pd,mat)
#        print (pd)
#        print (rpd)
        #need to check if these point are inside the cube using the dimension of the cube
#        res = numpy.less(rpd,radii[0])#size ?
        res = numpy.less_equal(numpy.fabs(rpd),numpy.array(radii[0])/2.)
#        print (res)
        c=numpy.average(res,1)#.astype(int)
#        print ("average",c)
        d=numpy.equal(c,1.)
#        print (d)
        ptinside = numpy.nonzero(d)[0]
#        print (ptinside)
#        if not len(ptinside):
#            print "no point inside the geom?"
#            return False
        if self.compareCompartment:
            ptinsideId = numpy.take(pointsInCube,ptinside,0)
            compIdsSphere = numpy.take(histoVol.grid.gridPtId,ptinsideId,0)  
#            print "compId",compIdsSphere
            if self.compNum <= 0 :
                wrongPt = [ cid for cid in compIdsSphere if cid != self.compNum ]
                if len(wrongPt):
#                    print wrongPt
                    return True                
#            
#        trigger, res = self.checkCompartment(numpy.take(pointsInCube,ptinside,0))
#        print ("checkCompartment result",trigger, res)
#        if res :
#            return True
#            for pt in pointsInCube:
        for pti in ptinside:
            pt = pointsInCube[pti]
#            print pt,distance[pt]
#                dist = dsq[pti]
#                if dist > rad2: continue # outside radius
            if distance[pt]<-0.0001 :#or trigger : # pt is inside cylinder
                #changeObjColorMat
#                    if histoVol.runTimeDisplay > 1:
#                        self.vi.changeObjColorMat(cyl,(1.,0.,0.))
#                        self.vi.update()
#                        sleep(1.0)
                return True

#            cylNum += 1
        return False

    def checkCubeCompart(self, centers1, centers2, radii, jtrans, rotMat,
                           gridPointsCoords, distance, histoVol):
        """
        Check cube for collision
        centers1 and centers2 should be the cornerPoints ?
        can also use the center plus size (radii), or the position/position2
        """
        cent1T = self.transformPoints(jtrans, rotMat, centers1)[0]#bb1
        cent2T = self.transformPoints(jtrans, rotMat, centers2)[0]#bb2
        center = self.transformPoints(jtrans, rotMat, [self.center,])[0]
        
        cylNum = 0
#        for radc, p1, p2 in zip(radii, cent1T, cent2T):
        x1, y1, z1 = cent1T
        x2, y2, z2 = cent2T
        vx, vy, vz = vect = (x2-x1, y2-y1, z2-z1)
        lengthsq = vx*vx + vy*vy + vz*vz
        l = sqrt( lengthsq )
        cx, cy, cz = posc = center#x1+vx*.5, y1+vy*.5, z1+vz*.5
        radt = l/2. + self.encapsulatingRadius
        x,y,z = posc
        bb = ( [x-radt, y-radt, z-radt], [x+radt, y+radt, z+radt] )
        
        pointsInCube = histoVol.grid.getPointsInCube(bb, posc, radt)
        
        pd = numpy.take(gridPointsCoords,pointsInCube,0)-center
        m = numpy.matrix(numpy.array(rotMat).reshape(4,4))#
        mat = m.I
        rpd = ApplyMatrix(pd,mat)
        res = numpy.less_equal(numpy.fabs(rpd),numpy.array(radii[0])/2.)
        c=numpy.average(res,1)#.astype(int)
        d=numpy.equal(c,1.)
        ptinside = numpy.nonzero(d)[0]
        ptinsideId = numpy.take(pointsInCube,ptinside,0)
        compIdsSphere = numpy.take(histoVol.grid.gridPtId,ptinsideId,0)  
#        print "compId",compIdsSphere
        if self.compNum <= 0 :
            wrongPt = [ cid for cid in compIdsSphere if cid != self.compNum ]
            if len(wrongPt):
#                print wrongPt
                return True                
        return False

    def checkPointComp(self,point):
        #if grid too sparse this will not work.
        #ptID = self.histoVol.grid.getPointFrom3D(point)
        cID = self.histoVol.getPointCompartmentId(point)
        #dist,ptID = self.histoVol.grid.getClosestGridPoint(point)
        #cID = self.histoVol.grid.gridPtId[ptID]              
        if self.compNum == 0 :
            organelle = self.histoVol
        else :
            organelle = self.histoVol.compartments[abs(self.compNum)-1]

        if self.compNum > 0 : #surface ingredient
            #r=compartment.checkPointInside_rapid(point,self.histoVol.grid.diag,ray=3)
            if self.Type == "Grow":
                #need a list of accepted compNum
                check = False
                if len(self.compMask):
                    if cID not in self.compMask:
                        check =  False
                    else :
                        check =  True
                else :
                    check =  True
#                if cID > 0 : #surface point look at surface cutoff 
#                    if dist < self.cutoff_surface :
#                        check = False
#                    else :
#                        check = True #grid probably too sparse, need to check where we are 
                return check
#        for i,o in self.histoVol.compartments:
#        if self.compNum != cID:
#            return False
#        else :
#            return True
        if self.compNum < 0 :#     
            inside = organelle.checkPointInside_rapid(point,self.histoVol.grid.diag,ray=3)
            if inside :#and cID < 0:
                return True
            else :
                return False
#            if inside and self.compNum >=0 :
#                return False
#            if not inside and self.compNum < 0 :
#                return False                
        if self.compNum == 0 : #shouldnt be in any compartments
            for o in self.histoVol.compartments:
                inside = o.checkPointInside_rapid(point,self.histoVol.grid.diag,ray=3)
                if inside :
                    return False
        if self.compNum != cID:
            return False
        else :
            return True


    def checkPointSurface(self,point,cutoff):
        if not hasattr(self,"histoVol") :
            return False
        if self.compNum == 0 :
            organelle = self.histoVol
        else :
            organelle = self.histoVol.compartments[abs(self.compNum)-1]
        compNum = self.compNum
        for o in self.histoVol.compartments:
            if self.compNum > 0 and o.name == organelle.name:
                continue
            #add some jitter to the cutoff ?
            #closest = o.OGsrfPtsBht.closestPointsArray(tuple(numpy.array([point,])), cutoff, 0)#default cutoff is 0.0
#            res = o.OGsrfPtsBht.closestPointsArrayDist2(tuple(numpy.array([point,])),self.histoVol.grid.diag*2.0)
            res = o.OGsrfPtsBht.query(tuple(numpy.array([point,])))
            if len(res) == 2 :
                d = res[0][0]
                #pt=res[1][0]
                if autopack.verbose :
                    print ("distance is ",d,cutoff,res,o.name,organelle.name)#d can be wrond for some reason,
                #d = autopack.helper.measure_distance(point,o.vertices[pt])
                if d < cutoff :
                    return True
                if compNum < 0 and o.name == organelle.name :
                    inside = o.checkPointInside(numpy.array(point),self.histoVol.grid.diag)
                    if autopack.verbose :
                        print("inside ? ",inside) 
                    if not inside : 
                        return True
                    
    #                        #print closest        
#                if closest[0] != -1 :
#                    return True
#                else :
#                    return False 
#                sfpts = o.surfacePointsCoords
#                delta = numpy.array(sfpts)-numpy.array(point)
#                delta *= delta
#                distA = numpy.sqrt( delta.sum(1) )
##                print len(distA)
#                test = distA < cutoff
#                if True in test:
#                    return True

        return False

    def testPoint(self,newPt):
        inComp = True
        closeS = False
        inside = self.histoVol.grid.checkPointInside(newPt,dist=self.cutoff_boundary,jitter=getNormedVectorOnes(self.jitterMax))
        #print ("testPoint",newPt,inside,self.jitterMax)    
        if inside :
            inComp = self.checkPointComp(newPt)
            if inComp :
                #check how far from surface ?
                closeS = self.checkPointSurface(newPt,cutoff=self.cutoff_surface)
        if autopack.verbose :
            print ("test",self.name,newPt,not inside, closeS,not inComp,
                   (not inside or closeS or not inComp))
        return not inside or closeS or not inComp
        
    def oneJitter(self,spacing,trans,rotMat):
#        spacing = histoVol.smallestProteinSize
        jx, jy, jz = self.jitterMax
        jitter = self.getMaxJitter(spacing)
        jitter2 = jitter * jitter
        compNum = self.compNum
        tx, ty, tz = trans
        verbose=False
        if jitter2 > 0.0:
            found = False
            while not found:
#                    dx = jx*jitter*gauss(0., 0.3)
#                    dy = jy*jitter*gauss(0., 0.3)
#                    dz = jz*jitter*gauss(0., 0.3)
                dx = jx*jitter*uniform(-1.0, 1.0) # These should be from -1 to 1, not from -0.5 to 0.5
                dy = jy*jitter*uniform(-1.0, 1.0)
                dz = jz*jitter*uniform(-1.0, 1.0)
                d2 = dx*dx + dy*dy + dz*dz
                if d2 < jitter2:
                    if compNum > 0: # jitter less among normal
                        #if self.name=='2uuh C4 SYNTHASE':
                        #    import pdb
                        #    pdb.set_trace()
                        dx, dy, dz, dum = numpy.dot(rotMat, (dx,dy,dz,0))
                    jtrans = (tx+dx, ty+dy, tz+dz)
                    found = True
                else:
                    if verbose :
                        print('JITTER REJECTED', d2, jitter2)
        else:
            jtrans = trans
            dx = dy = dz = 0.0
            # randomize rotation about axis
        if compNum>0:
            rotMatj = self.getAxisRotation(rotMat)
        else:
            if self.useRotAxis :# is not None :
                if sum(self.rotAxis) == 0.0 :
                    rotMatj=numpy.identity(4)
                else :
                    rotMatj=self.vi.rotation_matrix(random()*self.rotRange,self.rotAxis)
            else :
                rotMatj = rotMat.copy()
        return jtrans,rotMatj
        
    def getInsidePoints(self,grid,gridPointsCoords,dpad,distance,
                       centT=None,jtrans=None, rotMatj=None):
        insidePoints={}
        newDistPoints={}
        if self.modelType=='Spheres':
            for radc, posc in zip(self.radii[-1], centT):
             
                rad = radc + dpad
                x,y,z = posc
                
                bb = ( [x-rad, y-rad, z-rad], [x+rad, y+rad, z+rad] )
                #print ("pointsInCube",bb, posc, rad,radc,dpad)
                pointsInCube = grid.getPointsInCube(bb, posc, rad)

                delta = numpy.take(gridPointsCoords,pointsInCube,0)-posc
                delta *= delta
                distA = numpy.sqrt( delta.sum(1) )
                ptsInSphere = numpy.nonzero(numpy.less_equal(distA, rad))[0]

                for pti in ptsInSphere:
                    pt = pointsInCube[pti]
                    if pt in insidePoints: continue
                    dist = distA[pti]
                    d = dist-radc
                    if dist < radc:  # point is inside dropped sphere
                        if pt in insidePoints:
                            if d < insidePoints[pt]:
                                insidePoints[pt] = d
                        else:
                            insidePoints[pt] = d
                    elif d < distance[pt]: # point in region of influence
                        if pt in newDistPoints:
                            if d < newDistPoints[pt]:
                                newDistPoints[pt] = d
                        else:
                            newDistPoints[pt] = d
        elif self.modelType=='Cylinders':
#            print ("cent transformed \n",rotMatj,self.positions[-1],self.positions2[-1])            
            cent1T = self.transformPoints(jtrans, rotMatj, self.positions[-1])
            cent2T = self.transformPoints(jtrans, rotMatj, self.positions2[-1])
#            print ("cent transformed \n",cent1T,cent2T,rotMatj,self.positions[-1],self.positions2[-1])
            for radc, p1, p2 in zip(self.radii[-1], cent1T, cent2T):
                x1, y1, z1 = p1
                x2, y2, z2 = p2
                vx, vy, vz = vect = (x2-x1, y2-y1, z2-z1)
                lengthsq = vx*vx + vy*vy + vz*vz
                l = sqrt( lengthsq )
                cx, cy, cz = posc = x1+vx*.5, y1+vy*.5, z1+vz*.5
                radt = l + radc + dpad
                #bb = ( [cx-radt, cy-radt, cz-radt], [cx+radt, cy+radt, cz+radt] )
#                print (p1,p2,posc,radc)
                bb = self.correctBB(posc,posc,radt)#p1,p2,radc
#                print (bb,rotMatj)
                pointsInCube = grid.getPointsInCube(bb, posc, radt)
#                print ("correct BB ",bb," posc ",posc, " radt ",radt," poInCube",max(pointsInCube))
                if hasattr(self,"histoVol") and self.histoVol.runTimeDisplay > 1:
                    box = self.vi.getObject("insidePtBox")
                    if box is None:
                        box = self.vi.Box('insidePtBox', cornerPoints=bb,visible=1)
                    else :
#                        self.vi.toggleDisplay(box,False)
                        self.vi.updateBox(box,cornerPoints=bb)
                        self.vi.update()
                    sleep(1.)
                pd = numpy.take(gridPointsCoords,pointsInCube,0) - p1
                dotp = numpy.dot(pd, vect)
                rad2 = radc*radc
                d2toP1 = numpy.sum(pd*pd, 1)
                dsq = d2toP1 - dotp*dotp/lengthsq

                pd2 = numpy.take(gridPointsCoords,pointsInCube,0) - p2
                d2toP2 = numpy.sum(pd2*pd2, 1)

                for pti, pt in enumerate(pointsInCube):
                    if pt in insidePoints: continue

                    if dotp[pti] < 0.0: # outside 1st cap
                        d = sqrt(d2toP1[pti])
                        if d < distance[pt]: # point in region of influence
                            if pt in newDistPoints:
                                if d < newDistPoints[pt]:
                                    newDistPoints[pt] = d
                            else:
                                newDistPoints[pt] = d
                    elif dotp[pti] > lengthsq:
                        d = sqrt(d2toP2[pti])
                        if d < distance[pt]: # point in region of influence
                            if pt in newDistPoints:
                                if d < newDistPoints[pt]:
                                    newDistPoints[pt] = d
                            else:
                                newDistPoints[pt] = d
                    else:
                        d = sqrt(dsq[pti])-radc
                        if d < 0.:  # point is inside dropped sphere
                            if pt in insidePoints:
                                if d < insidePoints[pt]:
                                    insidePoints[pt] = d
                            else:
                                insidePoints[pt] = d
#                print ("ok",len(pointsInCube))
        elif self.modelType=='Cube': 
            insidePoints,newDistPoints = self.getDistancesCube(jtrans, rotMatj,gridPointsCoords, distance, grid)
        return insidePoints,newDistPoints

    def getIngredientsInBox(self,histoVol,jtrans,rotMat,compartment,afvi):
        if histoVol.windowsSize_overwrite :
            rad = histoVol.windowsSize
        else :
#            rad = self.minRadius*2.0# + histoVol.largestProteinSize + \
                #histoVol.smallestProteinSize + histoVol.windowsSize
            rad = self.minRadius + histoVol.largestProteinSize + \
                histoVol.smallestProteinSize + histoVol.windowsSize
        x,y,z = jtrans
        bb = ( [x-rad, y-rad, z-rad], [x+rad, y+rad, z+rad] )
        if self.modelType == "Cylinders":
            cent1T = self.transformPoints(jtrans, rotMat, self.positions[self.maxLevel])
            cent2T = self.transformPoints(jtrans, rotMat, self.positions2[self.maxLevel])
            bbs=[]
            for radc, p1, p2 in zip(self.radii, cent1T, cent2T):            
                bb = self.correctBB(p1,p2,radc)
                bbs.append(bb)
            #get min and max from all bbs
            maxBB = [0,0,0]
            minBB = [9999,9999,9999]
            for bb in bbs:
                for i in range(3):
                    if bb[0][i] < minBB[i]:
                        minBB[i] =bb[0][i]
                    if bb[1][i] > maxBB[i]:
                        maxBB[i] = bb[1][i]
                    if bb[1][i] < minBB[i]:
                        minBB[i] = bb[1][i]
                    if bb[0][i] > maxBB[i]:
                        maxBB[i] = bb[0][i]
            bb = [minBB,maxBB]
        if histoVol.runTimeDisplay > 1:
            box = self.vi.getObject("partBox")
            if box is None:
                box = self.vi.Box('partBox', cornerPoints=bb,visible=1)
            else :
                self.vi.toggleDisplay(box,True)
                self.vi.updateBox(box,cornerPoints=bb)
                self.vi.update()
#            sleep(1.0)
        pointsInCube = histoVol.grid.getPointsInCube(bb, jtrans, rad)
        #should we got all ingre from all recipes?
        #can use the kdtree for it...
        #maybe just add the surface if its not already the surface
        mingrs = [m for m in organelle.molecules if m[3] in pointsInCube]
        return mingrs

    def getIngredientsInTree(self,close_indice):
        if len(self.histoVol.rIngr):    
            ingrs= [self.histoVol.rIngr[i] for i in close_indice["indices"]]
            return [numpy.asarray(self.histoVol.rTrans)[close_indice["indices"]],
                    numpy.asarray(self.histoVol.rRot)[close_indice["indices"]],
                    ingrs,
                    close_indice["distances"]]
        else :
            return []
            
    def getListePartners(self,histoVol,jtrans,rotMat,organelle,afvi,
                         close_indice=None):
        nb_ingredients = 0        
        if close_indice is None :
            mingrs = self.getIngredientsInBox(histoVol,jtrans,rotMat,organelle,afvi)
            nb_ingredients = len(mingrs)
            mingrs = zip(*mingrs)
        else :
            mingrs = self.getIngredientsInTree(close_indice)
            nb_ingredients = len(close_indice)
        listePartner = []
        weight=0.
        if not len(mingrs) or not len(mingrs[2]) :
            print ("no close ingredient found")
            return [],[]       
        else :
            print ("nb close ingredient",self.name,len(mingrs), len(mingrs[2]),
                   len(close_indice),nb_ingredients)   
        listePartner =[]
        for i in range(len(mingrs[2])):
            ing = mingrs[2][i]
            t = mingrs[0][i]
#            print ("test "+ing.name,ing.o_name,ing.isAttractor,self.partners_name)
            if self.packingMode=="closePartner":
                if ing.o_name in self.partners_name :
#                    print ("is a partner of"+self.name)
                    listePartner.append([i,self.partners[ing.name],mingrs[3][i]])
#                                         autopack.helper.measure_distance(jtrans,mingrs[0][i])])
            if ing.isAttractor :#and self.compNum <= 0: #always attract! or rol a dice ?sself.excluded_partners.has_key(name)               
                if ing.name not in self.partners_name and self.name not in ing.excluded_partners_name \
                and ing.name not in self.excluded_partners_name :
                    print ("shoul attract "+self.name)
                    part = self.getPartner(ing.name)
                    if part is None :
                        part = self.addPartner(ing,weight=ing.weight)
                    if ing.distExpression is not None:
                        part.distExpression = ing.distExpression
                    #print "new Partner", part,part.name,part.weight
                    d=afvi.vi.measure_distance(jtrans,t)
                    listePartner.append([i,part,d])
        if not listePartner:
            print ("no partner found in close ingredient",self.packingMode)
            return [],[] 
        else : 
            print (len(listePartner)," partner found in close ingredient")
            return mingrs,listePartner

    def getTransform(self):      
        tTrans = self.vi.ToVec(self.vi.getTranslation(self.moving))
        rRot = self.vi.getMatRotation(self.moving)
        self.htrans.append(tTrans)
        avg = numpy.average(numpy.array(self.htrans))
        d=self.vi.measure_distance(tTrans,avg)
        #print "during",d,tTrans
        if d < 5.0:
#            print("during",d,tTrans)#,rRot
            return True
        else :
            return False

    def checkDistance(self,liste_nodes,point,cutoff):
        for node in liste_nodes:
            rTrans,rRot=self.histoVol.getRotTransRB(node)
            d=self.vi.measure_distance(rTrans,point)
            print ("checkDistance",d,d<cutoff)        

    def get_rapid_nodes(self,close_indice,curentpt,removelast=False,prevpoint=None):
        if self.compNum == 0 :
            organelle = self.histoVol
        else :
            organelle = self.histoVol.compartments[abs(self.compNum)-1]
        nodes = []
#        ingrCounter={}
#        a=numpy.asarray(self.histoVol.rTrans)[close_indice["indices"]]
#        b=numpy.array([curentpt,])
#        distances=spatial.distance.cdist(a,b)
        distances=close_indice["distances"]#spatial.distance.cdist(a,b)#close_indice["distance"]
        #print ("retrieve ",len(close_indice["indices"]),close_indice["indices"],distances)
        for nid,n in enumerate(close_indice["indices"]):
            if n == -1 :
                continue
            if n == len(close_indice["indices"]):
                continue
            if n >= len(self.histoVol.rIngr):
                continue
            if distances[nid] == 0.0 : continue
            ingr= self.histoVol.rIngr[n]
            jtrans=self.histoVol.rTrans[n]
            rotMat=self.histoVol.rRot[n]
            #print (self.name+" is close to "+ingr.name,jtrans,curentpt)
            if prevpoint != None :
                #print distances[nid],
                #if prevpoint == jtrans : continue
                d=self.vi.measure_distance(numpy.array(jtrans),numpy.array(prevpoint))
                if d==0.0:#distances[nid] == 0 : #same point
                    #print ("continue d=0",numpy.array(jtrans),numpy.array(prevpoint),d)
                    continue
            if self.Type == "Grow":
                #shouldnt we use sphere instead
                if self.name == ingr.name :
                    #dont want last n-2  point?
                    c = len(self.histoVol.rIngr)
#                    print jtrans,curentpt
#                    print ("whats ",n,nid,c,(n==c) or n==(c-1) or  (n==c-2),
#                                (nid==c) or nid==(c-1) or  (nid==c-2))
#                    raw_input()
                    if (n==c) or n==(c-1) :#or  (n==(c-2)):
                        continue
            if ingr.name in self.partners and self.Type == "Grow":
                #for now just do nothing
                c = len(self.histoVol.rIngr)
#                    print jtrans,curentpt
#                print ("whats ",n,nid,c,(n==c) or n==(c-1) or  (n==c-2),
#                            (nid==c) or nid==(c-1) or  (nid==c-2))
#                    raw_input()
                if (n==c) or n==(c-1) or (n==c-2):
                    continue                              
            if self.name in ingr.partners and ingr.Type=="Grow":
                c = len(self.histoVol.rIngr)
                if (n==c) or n==(c-1) or (n==c-2):
                    continue                
            #else :
            #   print (self.name+" is close to "+ingr.name,jtrans,curentpt)            
            if distances[nid] > (ingr.encapsulatingRadius+self.encapsulatingRadius)*self.histoVol.scaleER:
                #print (distances[nid][0],ingr.encapsulatingRadius+self.encapsulatingRadius)                
                continue
            node = ingr.get_rapid_model()
            #distance ? should be < ingrencapsRadius+self.encradius
            nodes.append([node,numpy.array(jtrans),numpy.array(rotMat[:3,:3],'f'),ingr])
        #append organelle rb nodes
        node = None    
        for o in self.histoVol.compartments:
            node = None
            if self.Type != "Grow":
                if self.compNum > 0 and o.name == organelle.name:
                    continue
            node = o.get_rapid_model() 
            if node is not None : 
                nodes.append([node,numpy.zeros((3), 'f'),numpy.identity(3, 'f'),o])
#        print len(nodes),nodes
        self.histoVol.nodes = nodes
        return nodes
        
    def get_rbNodes(self,close_indice,currentpt,removelast=False,prevpoint=None,
                     getInfo=False):
        #move around the rbnode and return it
        #self.histoVol.loopThroughIngr( self.histoVol.reset_rbnode )   
        if self.compNum == 0 :
            organelle = self.histoVol
        else :
            organelle = self.histoVol.compartments[abs(self.compNum)-1]
        nodes = []
        ingrCounter={}
#        a=numpy.asarray(self.histoVol.rTrans)[close_indice["indices"]]
#        b=numpy.array([currentpt,])
        distances=close_indice["distances"]#spatial.distance.cdist(a,b)#close_indice["distance"]
        #print ("retrieve ",len(close_indice["indices"]),close_indice["indices"],distances)
        for nid,n in enumerate(close_indice["indices"]):
            if n == -1 :
                continue
#            if n == len(close_indice["indices"]):
#                continue
            if n >= len(self.histoVol.rIngr):
                continue
            ingr= self.histoVol.rIngr[n]
            if len(distances):
                if distances[nid] == 0.0 : continue
                if distances[nid] > (ingr.encapsulatingRadius+self.encapsulatingRadius)*self.histoVol.scaleER:
                    continue
#            print ("distance",nid,n,distances[nid],(ingr.encapsulatingRadius+self.encapsulatingRadius)*self.histoVol.scaleER)
#            print ("get_rbNodes",nid,n,len(self.histoVol.rTrans)-1,distances[nid][0]) 
#            print self.name+" is close to "+ingr.name
            jtrans=self.histoVol.rTrans[n]
            rotMat=self.histoVol.rRot[n]
            if prevpoint != None :
                #if prevpoint == jtrans : continue
                d=self.vi.measure_distance(numpy.array(jtrans),numpy.array(prevpoint))
                if d == 0 : #same point
                    continue
            if self.Type == "Grow":
                if self.name == ingr.name :
                    c = len(self.histoVol.rIngr)
                    if (n==c) or n==(c-1) :#or  (n==(c-2)):
                        continue
            if ingr.name in self.partners and self.Type == "Grow":
                c = len(self.histoVol.rIngr)
                if (n==c) or n==(c-1) :#or (n==c-2):
                    continue   
            if self.name in ingr.partners and ingr.Type=="Grow":
                c = len(self.histoVol.rIngr)
                if (n==c) or n==(c-1) :#or (n==c-2):
                    continue                
#            if self.packingMode == 'hexatile' :
#                #no self collition for testing
#                if self.name == ingr.name :
#                    continue
            rbnode = ingr.get_rb_model(alt=(ingr.name==self.name))
            if getInfo :
                nodes.append([rbnode, jtrans, rotMat,ingr])
            else :
                nodes.append(rbnode)
#            print "get",ingr.name,self.name,rbnode,distances[nid],(ingr.encapsulatingRadius+self.encapsulatingRadius)
        #append organelle rb nodes
        for o in self.histoVol.compartments:
            if self.compNum > 0 and o.name == organelle.name:
                #this i notworking for growing ingredient like hair.
                #should had after second segments
                if self.Type != "Grow":
                    continue
                else :
                    #whats the current length 
                    if len(self.results) <= 1 :
                        continue
            orbnode = o.get_rb_model()        
            if orbnode is not None : 
                #test distance to surface ? 
                res = o.OGsrfPtsBht.query(tuple(numpy.array([currentpt,])))
                if len(res) == 2 :
                    d = res[0][0]
                    if d < self.encapsulatingRadius :
                        if not getInfo :
                            nodes.append(orbnode)
                        else :
                            nodes.append([orbnode, [0,0,0], numpy.identity(4),o])
#        if self.compNum < 0 or self.compNum == 0 :     
#            for o in self.histoVol.compartments:
#                if o.rbnode is not None : 
#                    if not getInfo :
#                        nodes.append(o.rbnode)
#        print ("GetNode ",len(nodes),nodes)
        self.histoVol.nodes = nodes
        return nodes

    def getClosePairIngredient(self,point,histoVol,cutoff=10.0): 
        R={"indices":[],"distances":[]}
        radius = [ingr.encapsulatingRadius for ingr in self.histoVol.rIngr]
        radius.append(self.encapsulatingRadius)
        pos = self.histoVol.rTrans[:]#).tolist()
        pos.append([point[0],point[1],point[2]])
        ind = len(pos)-1
        bht=bhtreelib.BHtree( pos, radius, 10)
        # find all pairs for which the distance is less than 1.1
        # times the sum of the radii
        pairs = bht.closePointsPairsInTree(1.) 
        for p in pairs :
            if p[0] == ind :
                R["indices"].append(p[1])
            elif p[1] == ind :
                R["indices"].append(p[0])
        bhtreelib.freeBHtree(bht)
        print ("getClosePairIngredient ",R)
        print ("all pairs ",pairs)
        print ("query was ind ",ind)
        return R
       
    def getClosestIngredient(self,point,histoVol,cutoff=10.0):
        #may have to rebuild the whale tree every time we add a point
        #grab the current result
        #set the bhtree
        #get closest ClosePoints() 
#        raw_input()
#        return self.getClosePairIngredient(point,histoVol,cutoff=cutoff) 
        R={"indices":[],"distances":[]}
        result = numpy.zeros(histoVol.totalNbIngr ).astype('i')        
        nb=0
        if histoVol.close_ingr_bhtree is not None :
            if histoVol.treemode == "bhtree":# "cKDTree"
                nb = histoVol.close_ingr_bhtree.closePoints((point[0],point[1],point[2]), cutoff, result )
                R["indices"] = result[:nb]
                return R
            else :
                #request kdtree
                nb=[]
                if len(histoVol.rTrans) >= 1 :
#                    nb = histoVol.close_ingr_bhtree.query_ball_point(point,cutoff)
#                else :#use the general query, how many we want
                    distance,nb = histoVol.close_ingr_bhtree.query(point,len(histoVol.rTrans),distance_upper_bound=cutoff)#len of ingr posed so far
                    if len(histoVol.rTrans) == 1 :
                        distance=[distance]
                        nb=[nb]
                    R["indices"] = nb
                    R["distances"] = distance#sorted by distance short -> long
                return R
        else :
            return R
#        closest = histoVol.close_ingr_bhtree.closestPointsArray(tuple(numpy.array([point,])), cutoff, 0)#returnNullIfFail
#        print ("getClosestIngredient",closest,cutoff )       
#        return closest

    def update_data_tree(self,jtrans,rotMatj,ptInd=0,pt1=None,pt2=None):
        #self.histoVol.static.append(rbnode)
        #self.histoVol.moving = None
        self.histoVol.nb_ingredient+=1
        self.histoVol.rTrans.append(numpy.array(jtrans).flatten())
        self.histoVol.rRot.append(numpy.array(rotMatj))#rotMatj 
        self.histoVol.rIngr.append(self)
        if pt1 is not None :
            self.histoVol.result.append([ [numpy.array(pt1).flatten(),
                                numpy.array(pt2).flatten()], rotMatj, 
                            self, ptInd ])
        else :
            self.histoVol.result.append([ jtrans, rotMatj, 
                            self, ptInd ])
        if self.histoVol.treemode == "bhtree":# "cKDTree"
            if len(self.histoVol.rTrans) >= 1 : bhtreelib.freeBHtree(self.histoVol.close_ingr_bhtree)
            self.histoVol.close_ingr_bhtree=bhtreelib.BHtree( self.histoVol.rTrans, None, 10)
        else :
            if len(self.histoVol.rTrans) >= 1 :
                self.histoVol.close_ingr_bhtree= spatial.cKDTree(self.histoVol.rTrans, leafsize=10)

    def reject(self,):
        # got rejected
        self.haveBeenRejected = True
        self.rejectionCounter += 1
        if verbose :
                print('Failed ingr:%s rejections:%d'%(
                self.name, self.rejectionCounter))
        if self.rejectionCounter >= self.rejectionThreshold: #Graham set this to 6000 for figure 13b (Results Fig 3 Test1) otehrwise it fails to fill small guys
                #if verbose :
                print('PREMATURE ENDING of ingredient', self.name)
                self.completion = 1.0
        
    def place(self,histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,usePP,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None):
        success = False
        #print self.placeType
        self.vi = autopack.helper
#        if histoVol.afviewer != None: 
#            self.vi = histoVol.afviewer.vi
        self.histoVol=histoVol        
        #grow doesnt use panda.......but could use all the geom produce by the grow as rb
        if self.placeType == "jitter" or self.Type == "Grow" or self.Type == "Actine":
            success, nbFreePoints = self.jitter_place(histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None,usePP=usePP)
        elif self.placeType == "spring" or self.placeType == "rigid-body":
            success, nbFreePoints = self.rigid_place(histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None,usePP=usePP)
        elif self.placeType == "pandaDev":
            success, nbFreePoints = self.pandaBullet_place_dev(histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None,usePP=usePP)
        elif self.placeType == "pandaBullet":
            success, nbFreePoints = self.pandaBullet_place(histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None,usePP=usePP)
        elif self.placeType == "pandaBulletRelax" or self.placeType == "pandaBulletSpring":
            success, nbFreePoints = self.pandaBullet_relax(histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None)
        elif self.placeType == "RAPID":
            success, nbFreePoints = self.rapid_place(histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None,usePP=usePP)
        #freePoints and distance are changed .. return them and see if they aredifferent
#        print ("multiprocessor",usePP)
#        if usePP :
#            return success, nbFreePoints, freePoints ,distance, histoVol.molecules         
#        else :
        return success, nbFreePoints
 
    def place_mp(self,histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None):
        if self.compNum == 0 :
            organelle = histoVol
        else :
            organelle = histoVol.compartments[abs(self.compNum)-1]
        success = False
        #print self.placeType
        self.vi = autopack.helper
#        if histoVol.afviewer != None: 
#            self.vi = histoVol.afviewer.vi
        self.histoVol=histoVol        

        if self.placeType == "jitter" or self.Type == "Grow" or self.Type == "Actine":
            success, insidePts,newDistancePts = self.jitter_place(histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None,usePP=True)
        #at this points, molecules and some other variable have been changed
        return success, self, insidePts, newDistancePts , histoVol ,organelle.molecules        
        
    def rigid_place(self, histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None):
        """
        drop the ingredient on grid point ptInd
        """
        #print "rigid",self.placeType
        self.vi = histoVol.afviewer.vi
        afvi = histoVol.afviewer
        windowsSize = histoVol.windowsSize
        simulationTimes = histoVol.simulationTimes
        runTimeDisplay = histoVol.runTimeDisplay
        springOptions = histoVol.springOptions
        self.histoVol = histoVol
        rejectionCount = 0
        spacing = histoVol.smallestProteinSize
        jx, jy, jz = self.jitterMax
        jitter = self.getMaxJitter(spacing)
        jitter2 = jitter * jitter

        if self.compNum == 0 :
            compartment = self.histoVol
        else :
            compartment = self.histoVol.compartments[abs(self.compNum)-1]
            #this is hisotVol for cytoplasme
        compNum = self.compNum
        radius = self.minRadius

        gridPointsCoords = histoVol.grid.masterGridPositions

        # compute rotation matrix rotMat
        if compNum>0:
            # for surface points we compute the rotation which
            # aligns the principalVector with the surface normal
            vx, vy, vz = v1 = self.principalVector
            v2 = compartment.surfacePointsNormals[ptInd]
            try :
                rotMat = numpy.array( rotVectToVect(v1, v2 ), 'f')
            except :
                print('PROBLEM ', self.name)
                rotMat = numpy.identity(4)
        else:
            if self.useRotAxis :
                if sum(self.rotAxis) == 0.0 :
                    rotMat=numpy.identity(4)
                else :
                    rotMat=afvi.vi.rotation_matrix(random()*self.rotRange,self.rotAxis)
            else :
                rotMat=histoVol.randomRot.get()

        if verbose :
            pass#print('%s nbs:%2d'%(self.pdb, len(self.positions[0])), end=' ')

        # jitter position loop
        jitterList = []
        collD1 = []
        collD2 = []

        trans = gridPointsCoords[ptInd] # drop point
        gridDropPoint = trans
        jtrans,rotMatj = self.oneJitter(spacing,trans,rotMat)
        
        ok = False
        #here should go the simulation
        #1- we build the ingrediant if not already and place the ingrediant at jtrans, rotMatj
        moving=None
        static=[]
        target=None
        targetPoint = jtrans
#        import c4d
        #c4d.documents.RunAnimation(c4d.documents.GetActiveDocument(), True)

        if self.mesh:
            if hasattr(self,"mesh_3d"):
                #create an instance of mesh3d and place it
                name = self.name + str(ptInd)
                if self.mesh_3d is None :
                    self.moving= moving = afvi.vi.Sphere(name,radius=self.radii[0][0],
                                                    parent=afvi.movingMesh)[0]
                    afvi.vi.setTranslation(moving,pos=jtrans)
                else :
                    #why the GetDown?
                    self.moving= moving = afvi.vi.newInstance(name,self.mesh_3d,#.GetDown()
                                                matrice=rotMatj,
                                                location=jtrans, parent = afvi.movingMesh)
        #2- get the neighboring object from ptInd
        mingrs,listePartner=self.getListePartners(histoVol,jtrans,rotMat,compartment,afvi)
        for i,elem in enumerate(mingrs):
            ing = elem[2]
            t = elem[0]
            r = elem[1]
            ind = elem[3]
            #print "neighbour",ing.name
            if hasattr(ing,"mesh_3d"):
                #create an instance of mesh3d and place it
                name = ing.name + str(ind)
                if ing.mesh_3d is None :
                    ipoly = afvi.vi.Sphere(name,radius=self.radii[0][0],parent=afvi.staticMesh)[0]
                    afvi.vi.setTranslation(ipoly,pos=t)
                else :
                    ipoly =afvi.vi.newInstance(name,ing.mesh_3d,matrice=r,#.GetDown()
                           location=t, parent = afvi.staticMesh)
                static.append(ipoly)
            elif isinstance(ing,GrowIngrediant):
                name = ing.name + str(ind)
                ipoly =afvi.vi.newInstance(name,afvi.orgaToMasterGeom[ing],
                                           parent = afvi.staticMesh)
                static.append(ipoly)
            
        if listePartner : #self.packingMode=="closePartner":
            if verbose:
                print("len listePartner = ", len(listePartner))
            if not self.force_random:
                targetPoint,weight = self.pickPartner(mingrs,listePartner,currentPos=jtrans)
                if targetPoint is None :
                    targetPoint = jtrans
            else :
                targetPoint = jtrans
#        print "targetPt",len(targetPoint),targetPoint   
        #setup the target position
        if self.placeType == "spring":
            afvi.vi.setRigidBody(afvi.movingMesh,**histoVol.dynamicOptions["spring"])
            #target can be partner position?
            target = afvi.vi.getObject("target"+name)
            if target is None :
                target = afvi.vi.Sphere("target"+name,radius=5.0)[0]
            afvi.vi.setTranslation(target,pos=targetPoint)
            afvi.vi.addObjectToScene(None,target)
            #3- we setup the spring (using the sphere position empty)
            spring = afvi.vi.getObject("afspring")
            if spring is None :
                spring = afvi.vi.createSpring("afspring",targetA=moving,tragetB=target,**springOptions)
            else :
                afvi.vi.updateSpring(spring,targetA=moving,tragetB=target,**springOptions)
        else :
            #before assigning should get outside thge object 
            afvi.vi.setRigidBody(afvi.movingMesh,**histoVol.dynamicOptions["moving"])
            afvi.vi.setTranslation(self.moving,pos=targetPoint)
        afvi.vi.setRigidBody(afvi.staticMesh,**histoVol.dynamicOptions["static"])
        #4- we run the simulation
        #c4d.documents.RunAnimation(c4d.documents.GetActiveDocument(), False,True)
        #if runTimeDisplay :
        afvi.vi.update()
#        rTrans = afvi.vi.ToVec(afvi.vi.getTranslation(moving))
#        rRot = afvi.vi.getMatRotation(moving)

        #print afvi.vi.ToVec(moving.GetAllPoints()[0])
        #afvi.vi.animationStart(duration = simulationTimes)
        #afvi.vi.update()
        afvi.vi.frameAdvanced(duration = simulationTimes,display = runTimeDisplay)#,
                              #cb=self.getTransfo)
#        moving=afvi.vi.makeEditable(moving,copy=False)                            
        #5- we get the resuling transofrmation matrix and decompose ->rTrans rRot
        #if runTimeDisplay :
        afvi.vi.update()
        rTrans = afvi.vi.ToVec(afvi.vi.getTranslation(moving))
        rRot = afvi.vi.getMatRotation(moving)
#        M=moving.GetMg()
        #print afvi.vi.ToVec(moving.GetAllPoints()[0])

#        print("OK AFTER",rTrans)#,rRot
#        print("save",self.tTrans)#,self.rRot
        #6- clean and delete everything except the spring
        afvi.vi.deleteObject(moving)        
        afvi.vi.deleteObject(target)
        for o in static:
            afvi.vi.deleteObject(o)
        ok = True
        jtrans = rTrans[:]
        rotMatj = rRot[:]
        jitterPos = 1
        if ok :
            ## get inside points and update distance
            ##
            # use best sperical approcimation
#            print(">>?",self.name,jtrans)
            centT = self.transformPoints(jtrans, rotMatj, self.positions[-1])

            insidePoints = {}
            newDistPoints = {}
            insidePoints,newDistPoints = self.getInsidePoints(histoVol.grid,
                                gridPointsCoords,dpad,distance,centT=centT,
                                jtrans=jtrans, 
                                rotMatj=rotMatj)

            
            # update free points
            nbFreePoints = self.updateDistances(histoVol,insidePoints, newDistPoints, 
                        freePoints, nbFreePoints, distance, 
                        histoVol.grid.masterGridPositions, verbose)

            # save dropped ingredient
            
            compartment.molecules.append([ jtrans, rotMatj, self, ptInd ])
            histoVol.order[ptInd]=histoVol.lastrank
            histoVol.lastrank+=1
            histoVol.nb_ingredient+=1
#            raw_input()
            histoVol.rTrans.append(jtrans)
            histoVol.result.append([ jtrans, rotMatj, self, ptInd ])
            histoVol.rRot.append(rotMatj)
            histoVol.rIngr.append(self)            
#            histoVol.close_ingr_bhtree.InsertRBHPoint((jtrans[0],jtrans[1],jtrans[2]),radius,None,histoVol.nb_ingredient)
#            histoVol.close_ingr_bhtree.MoveRBHPoint(histoVol.nb_ingredient,jtrans,0)
            
            self.rRot.append(rotMatj)
            self.tTrans.append(jtrans)
            # add one to molecule counter for this ingredient
            self.counter += 1
            self.completion = float(self.counter)/self.nbMol

            if jitterPos>0:
                histoVol.successfullJitter.append(
                    (self, jitterList, collD1, collD2) )  
            if verbose :
                print('Success nbfp:%d %d/%d dpad %.2f'%(
                nbFreePoints, self.counter, self.nbMol, dpad))
            if self.name=='in  inside':
                histoVol.jitterVectors.append( (trans, jtrans) )

            success = True
            self.rejectionCounter = 0
            
        else: # got rejected
            success = False
            histoVol.failedJitter.append(
                (self, jitterList, collD1, collD2) )

            distance[ptInd] = max(0, distance[ptInd]*0.9)
            self.rejectionCounter += 1
            if verbose :
                print('Failed ingr:%s rejections:%d'%(
                self.name, self.rejectionCounter))
            if self.rejectionCounter >= 30:
                if verbose :
                    print('PREMATURE ENDING of ingredient', self.name)
                self.completion = 1.0
        return success, nbFreePoints
    

    def jitter_place(self, histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None,drop=True,usePP=False):
        """
        drop the ingredient on grid point ptInd
        """
#        print("JitterPlace1****************************************************************")
        verbose = histoVol.verbose
        afvi = histoVol.afviewer
        rejectionCount = 0
        spacing = histoVol.smallestProteinSize
        jx, jy, jz = self.jitterMax
        jitter = histoVol.callFunction(self.getMaxJitter,(spacing,))
        jitter2 = jitter * jitter

        if self.compNum == 0 :
            compartment = histoVol
        else :
            compartment = histoVol.compartments[abs(self.compNum)-1]
        compNum = self.compNum
        radius = self.minRadius
        runTimeDisplay = histoVol.runTimeDisplay
        
        gridPointsCoords = histoVol.masterGridPositions
        
        #test force dpad to 1
 #       dpad = 1  Dangerous- this is specifically what is breaking all of the other Sphere/Cylinder test files
        # compute rotation matrix rotMat
        if compNum>0 :
            # for surface points we compute the rotation which
            # aligns the principalVector with the surface normal
            vx, vy, vz = v1 = self.principalVector
                
#            try :
            v2 = compartment.surfacePointsNormals[ptInd]
            rotMat = numpy.array( rotVectToVect(v1, v2 ), 'f')
#            except :
#                print('############ PROBLEM ', self.name, ptInd,len(compartment.surfacePointsNormals))
#                rotMat = numpy.identity(4)
        else:
            #this is where we could apply biased rotatio ie gradient/attractor
            if self.useRotAxis :
#                angle = random()*6.2831#math.radians(random()*360.)#random()*pi*2.
#                print "angle",angle,math.degrees(angle)
#                direction = self.rotAxis
                if sum(self.rotAxis) == 0.0 :
                    rotMat=numpy.identity(4)
                elif self.useOrientBias and self.packingMode =="gradient":#you need a gradient here
                    rotMat=self.alignRotation(gridPointsCoords[ptInd] )
                else :
                    rotMat=afvi.vi.rotation_matrix(random()*self.rotRange,self.rotAxis)
            # for other points we get a random rotation
            else :
                rotMat=histoVol.randomRot.get()

#        if verbose :
#            pass#print('%s nbs:%2d'%(self.pdb, len(self.positions[0])), end=' ')

        # jitter position loop
        jitterList = []
        collD1 = []
        collD2 = []

        trans = gridPointsCoords[ptInd] # drop point, surface points.
        targetPoint = trans
        moving = None
        if runTimeDisplay and self.mesh:
            if hasattr(self,"mesh_3d"):
                #create an instance of mesh3d and place it
                name = self.name + str(ptInd)
                moving = afvi.vi.getObject(name)
                if moving is None :
                    if self.mesh_3d is None :
                        moving = afvi.vi.Sphere(name,radius=self.radii[0][0],
                                                        parent=afvi.staticMesh)[0]
                        afvi.vi.setTranslation(moving,pos=targetPoint)
                    else :
                        moving=  afvi.vi.newInstance(name,self.mesh_3d,#.GetDown(),
                                                    matrice=rotMat,
                                                    location=targetPoint, parent = afvi.staticMesh)
#                else :   #Graham turned off this unnecessary update
                    #afvi.vi.setTranslation(moving,pos=targetPoint)#rot?
#                    mat = rotMat.copy() #Graham turned off this unnecessary update
#                    mat[:3, 3] = targetPoint #Graham turned off this unnecessary update
#                    afvi.vi.setObjectMatrix(moving,mat,transpose=True) #Graham turned off this unnecessary update
#                afvi.vi.update()  #Graham turned off this unnecessary update
#            print ('Must check for collision here before accepting final position')
        #do we get the list of neighbours first > and give a different trans...closer to the partner
        #we should look up for an available ptID around the picked partner if any
        #getListPartner
        if histoVol.ingrLookForNeighbours:
            mingrs,listePartner=self.getListePartners(histoVol,trans,rotMat,compartment,afvi)
            #if liste:pickPartner
            if listePartner : #self.packingMode=="closePartner":
#                print "ok partner",len(listePartner)
                if not self.force_random:
                    targetPoint,weight = self.pickPartner(mingrs,listePartner,currentPos=trans)
                    if targetPoint is None :
                        targetPoint = trans
                    else :#maybe get the ptid that can have it
                        #find a newpoint here?
                        x,y,z = targetPoint
                        rad = self.radii[0][0]*2.#why by 2.0
                        bb = ( [x-rad, y-rad, z-rad], [x+rad, y+rad, z+rad] )
                        pointsInCube = histoVol.grid.getPointsInCube(bb, targetPoint,rad )
                        #is one of this point can receive the current ingredient
                        cut  = rad-jitter
                        for pt in pointsInCube:
                            d = distance[pt]
                            if d>=cut:
                                #lets just take the first one
                                targetPoint = gridPointsCoords[pt]
                                break
                else :
                    targetPoint = trans
            #if partner:pickNewPoit like in fill3
        tx, ty, tz = jtrans = targetPoint
        gridDropPoint = targetPoint        
        #we may increase the jitter, or pick from xyz->Id free for its radius
        
        #jitter loop
        t1 = time()
        for jitterPos in range(self.nbJitter):   #This expensive Gauusian rejection system should not be the default should it?
            # jitter points location
            if jitter2 > 0.0:
                found = False
                while not found:
#                    dx = jx*jitter*gauss(0., 0.3)
#                    dy = jy*jitter*gauss(0., 0.3)
#                    dz = jz*jitter*gauss(0., 0.3)
                    dx = jx*jitter*uniform(-1.0, 1.0) # These should be from -1 to 1, not from -0.5 to 0.5
                    dy = jy*jitter*uniform(-1.0, 1.0)
                    dz = jz*jitter*uniform(-1.0, 1.0)
                    d2 = dx*dx + dy*dy + dz*dz
                    if d2 < jitter2:
                        if compNum > 0: # jitter less among normal
                            #if self.name=='2uuh C4 SYNTHASE':
                            #    import pdb
                            #    pdb.set_trace()
                            dx, dy, dz, dum = numpy.dot(rotMat, (dx,dy,dz,0))
                        jtrans = (tx+dx, ty+dy, tz+dz)
                        found = True
                    else:
                        if verbose :
                            print('JITTER REJECTED', d2, jitter2)
#                    if runTimeDisplay and moving is not None :  #Graham turned off this unnecessary update
#                        afvi.vi.setTranslation(moving,pos=jtrans)   #Graham turned off this unnecessary update
#                        afvi.vi.update()  #Graham turned off this unnecessary update
#                        print "ok moving"
            else:
                jtrans = targetPoint
                dx = dy = dz = 0.0
           
            histoVol.totnbJitter += 1
            histoVol.jitterLength += dx*dx + dy*dy + dz*dz  #Why is this expensive line needed?
            jitterList.append( (dx,dy,dz) )
            
            #print 'j%d %.2f,%.2f,%.2f,'%(jitterPos,tx, ty, tz),
#            if verbose :
#                print('j%d'%jitterPos)# end=' '
            
            # loop over all spheres representing ingredient
            modSphNum = 1
            if sphGeom is not None:
                modCent = []
                modRad = []

            ## check for collisions 
            ## 
            level = self.collisionLevel

            # randomize rotation about axis
            if compNum>0:
                rotMatj = self.getAxisRotation(rotMat)
            else:
                if self.useRotAxis :
                    if sum(self.rotAxis) == 0.0 :
                        rotMatj=numpy.identity(4)  #Graham Oct 16,2012 Turned on always rotate below as default.  If you want no rotation
                                                   #set useRotAxis = 1 and set rotAxis = 0, 0, 0 for that ingredient
                    elif self.useOrientBias and self.packingMode =="gradient":
#                        rotMatj = self.getAxisRotation(rotMat)
                        rotMatj = self.getBiasedRotation(rotMat,weight=None)
#                            weight = 1.0 - self.histoVol.gradients[self.gradient].weight[ptInd])
                    else :
                        #should we align to this rotAxis ?
                        rotMatj=afvi.vi.rotation_matrix(random()*self.rotRange,self.rotAxis)
                else :
                    rotMatj=histoVol.randomRot.get()  #Graham turned this back on to replace rotMat.copy() so ing rotate each time
#                    rotMatj = rotMat.copy()
            if runTimeDisplay and moving is not None :
#                print "ok rot copy"
                mat = rotMatj.copy()
                mat[:3, 3] = jtrans
                afvi.vi.setObjectMatrix(moving,mat,transpose=True)
#                afvi.vi.setTranslation(moving,pos=jtrans)
                afvi.vi.update()
            collision = True
            #perodicity check
#            periodic_pos = self.histoVol.grid.getPositionPeridocity(jtrans,
#                    getNormedVectorOnes(self.jitterMax),self.encapsulatingRadius)                
##            histoVol.callFunction(histoVol.moveRBnode,(rbnode, jtrans, rotMatj,))
#            perdiodic_collision = False
            r=[False]
#            if periodic_pos is not None and self.packingMode !="gradient" :
##                print ("OK Periodicity ",len(periodic_pos),periodic_pos)
#                for p in periodic_pos :
#                    perdiodic_collision = self.collision_jitter(p, rotMatj,
#                                    level, gridPointsCoords, distance, histoVol)
#                    r.extend([perdiodic_collision])
#                    if  True in r : break
#                    if runTimeDisplay and moving is not None :
#                        mat = rotMatj.copy()
#                        mat[:3, 3] = jtrans
#                        afvi.vi.setObjectMatrix(moving,mat,transpose=True)
#                        afvi.vi.update()
            test=False#self.testPoint(jtrans)
            if not test and not ( True in r) :                    
                collision = self.collision_jitter(jtrans, rotMatj,
                                    level, gridPointsCoords, distance, histoVol)
            if not collision:
                break # break out of jitter pos loop
        if verbose: 
            print("jitter loop ",time()-t1)
        if not collision and not  ( True in r) :
            ## get inside points and update distance
            ##
            # use best sperical approcimation
            centT = self.centT#self.transformPoints(jtrans, rotMatj, self.positions[-1])

            insidePoints = {}
            newDistPoints = {}
            t3=time()
            #should be replace by self.getPointInside
            if self.modelType=='Spheres':
#                for pointsInCube,ptsInSphere in self.distances_temp:
                for radc, posc in zip(self.radii[-1], centT):
#                for i,radc in enumerate(self.radii[-1]):
#                    pointsInCube,ptsInSphere,distA = self.distances_temp[i]
                    rad = radc + dpad
#                    print ("sphere ingr",rad,dpad,radc,posc)
                    x,y,z = posc
                    #this have already be done in the checkCollision why doing it again
                    bb = ( [x-rad, y-rad, z-rad], [x+rad, y+rad, z+rad] )
                    pointsInCube = histoVol.callFunction(histoVol.grid.getPointsInCube,
                                                         (bb, posc, rad))
              #      print ("sphere ingr",len(pointsInCube))
                    delta = numpy.take(gridPointsCoords,pointsInCube,0)-posc
                    delta *= delta
                    distA = numpy.sqrt( delta.sum(1) )
                    ptsInSphere = numpy.nonzero(numpy.less_equal(distA, rad))[0]

                    for pti in ptsInSphere:
                        pt = pointsInCube[pti]
                        if pt in insidePoints: continue
                        dist = distA[pti]
                        d = dist-radc
                        if dist < radc:  # point is inside dropped sphere
                            if pt in insidePoints:
                                if d < insidePoints[pt]:
                                    insidePoints[pt] = d
                            else:
                                insidePoints[pt] = d
                        elif d < distance[pt]: # point in region of influence
                            if pt in newDistPoints:
                                if d < newDistPoints[pt]:
                                    newDistPoints[pt] = d
                            else:
                                newDistPoints[pt] = d

            elif self.modelType=='Cylinders':
                cent1T = self.transformPoints(jtrans, rotMatj, self.positions[-1])
                cent2T = self.transformPoints(jtrans, rotMatj, self.positions2[-1])

                for radc, p1, p2 in zip(self.radii[-1], cent1T, cent2T):

                    x1, y1, z1 = p1
                    x2, y2, z2 = p2
                    vx, vy, vz = vect = (x2-x1, y2-y1, z2-z1)
                    lengthsq = vx*vx + vy*vy + vz*vz
                    l = sqrt( lengthsq )
                    cx, cy, cz = posc = x1+vx*.5, y1+vy*.5, z1+vz*.5
                    radt = l + radc + dpad
                    bb = ( [cx-radt, cy-radt, cz-radt], [cx+radt, cy+radt, cz+radt] )
                    pointsInCube = histoVol.callFunction(histoVol.grid.getPointsInCube,
                                                         (bb, posc, radt))

                    pd = numpy.take(gridPointsCoords,pointsInCube,0) - p1
                    dotp = numpy.dot(pd, vect)
                    rad2 = radc*radc
                    d2toP1 = numpy.sum(pd*pd, 1)
                    dsq = d2toP1 - dotp*dotp/lengthsq

                    pd2 = numpy.take(gridPointsCoords,pointsInCube,0) - p2
                    d2toP2 = numpy.sum(pd2*pd2, 1)

                    for pti, pt in enumerate(pointsInCube):
                        if pt in insidePoints: continue

                        if dotp[pti] < 0.0: # outside 1st cap
                            d = sqrt(d2toP1[pti])
                            if d < distance[pt]: # point in region of influence
                                if pt in newDistPoints:
                                    if d < newDistPoints[pt]:
                                        newDistPoints[pt] = d
                                else:
                                    newDistPoints[pt] = d
                        elif dotp[pti] > lengthsq:
                            d = sqrt(d2toP2[pti])
                            if d < distance[pt]: # point in region of influence
                                if pt in newDistPoints:
                                    if d < newDistPoints[pt]:
                                        newDistPoints[pt] = d
                                else:
                                    newDistPoints[pt] = d
                        else:
                            d = sqrt(dsq[pti])-radc
                            if d < 0.:  # point is inside dropped sphere
                                if pt in insidePoints:
                                    if d < insidePoints[pt]:
                                        insidePoints[pt] = d
                                else:
                                    insidePoints[pt] = d
            elif self.modelType=='Cube':
                insidePoints,newDistPoints=self.getDistancesCube(jtrans, rotMatj,gridPointsCoords, distance, histoVol)
            
            # save dropped ingredient
            if verbose:
                print("compute distance loop ",time()-t3)

            if drop:
                #print "ok drop",compartment.name,self.name
                compartment.molecules.append([ jtrans, rotMatj, self, ptInd ])
                histoVol.order[ptInd]=histoVol.lastrank
                histoVol.lastrank+=1
                
#                histoVol.close_ingr_bhtree.MoveRBHPoint(histoVol.nb_ingredient,jtrans,0)
                histoVol.nb_ingredient+=1
#                histoVol.close_ingr_bhtree.InsertRBHPoint((jtrans[0],jtrans[1],jtrans[2]),radius,None,histoVol.nb_ingredient)
                if False :#usePP
                    self.counter += 1
                    self.completion = float(self.counter)/float(self.nbMol)
        
                    if jitterPos>0:
                        histoVol.successfullJitter.append(
                            (self, jitterList, collD1, collD2) )
                       
                    #if verbose :
                    print('Success nbfp:%d %d/%d dpad %.2f'%(
                        nbFreePoints, self.counter, self.nbMol, dpad))
                    if self.name=='in  inside':
                        histoVol.jitterVectors.append( (trans, jtrans) )
        
                    success = True
                    self.rejectionCounter = 0
                    return True, insidePoints, newDistPoints
            # update free points
            if verbose:
                print ("updating distances and verbose = ", verbose)
            timeUpDistLoopStart=time()
            nbFreePoints = histoVol.callFunction(self.updateDistances,(histoVol,insidePoints,
                                                newDistPoints,freePoints,nbFreePoints, distance,
                                                histoVol.masterGridPositions, verbose))
            histoVol.timeUpDistLoopTotal+=time()-timeUpDistLoopStart
            

            if sphGeom is not None:
                for po1, ra1 in zip(modCent, modRad):
                    sphCenters.append(po1)
                    sphRadii.append(ra1)
                    sphColors.append(self.color)

            if labDistGeom is not None:
                verts = []
                labels = []
                colors = []
                #for po1, d1,d2 in distChanges.values():
                fpts = freePoints
                for i in range(nbFreePoints):
                    pt = fpts[i]
                    verts.append(histoVol.masterGridPositions[pt])
                    labels.append( "%.2f"%distance[pt])                    
#                for pt in freePoints[:nbFreePoints]:
#                    verts.append(histoVol.masterGridPositions[pt])
#                    labels.append( "%.2f"%distance[pt])
                labDistGeom.Set(vertices=verts, labels=labels)
                #materials=colors, inheritMaterial=0)

            # add one to molecule counter for this ingredient
            self.counter += 1
            self.completion = float(self.counter)/float(self.nbMol)

            if jitterPos>0:
                histoVol.successfullJitter.append(
                    (self, jitterList, collD1, collD2) )
               
            if verbose :
                print('Success nbfp:%d %d/%d dpad %.2f'%(
                nbFreePoints, self.counter, self.nbMol, dpad))
            if self.name=='in  inside':
                histoVol.jitterVectors.append( (trans, jtrans) )

            success = True
            self.rejectionCounter = 0
            
        else: # got rejected
            if runTimeDisplay and moving is not None :
                afvi.vi.deleteObject(moving)
            success = False
            self.haveBeenRejected = True
            histoVol.failedJitter.append(
                (self, jitterList, collD1, collD2) )

            distance[ptInd] = max(0, distance[ptInd]*0.9)# ???
            self.rejectionCounter += 1
            if verbose :
                print('Failed ingr:%s rejections:%d'%(
                self.name, self.rejectionCounter))
            if self.rejectionCounter >= self.rejectionThreshold: #Graham set this to 6000 for figure 13b (Results Fig 3 Test1) otehrwise it fails to fill small guys
                #if verbose :
                print('PREMATURE ENDING of ingredient', self.name)
                self.completion = 1.0
#            if usePP :
#                return False,None,None
        if sphGeom is not None:
            sphGeom.Set(vertices=sphCenters, radii=sphRadii,
                        materials=sphColors)
            sphGeom.viewer.OneRedraw()
            sphGeom.viewer.update()

        if drop :
            return success, nbFreePoints
        else :
            return success, nbFreePoints,jtrans, rotMatj

    def lookForNeighbours(self,trans,rotMat,organelle,afvi,distance,closest_indice=None):
        mingrs,listePartner=self.getListePartners(self.histoVol,trans,rotMat,
                                organelle,afvi,close_indice=closest_indice)
        targetPoint = trans
        found = False
        if listePartner : #self.packingMode=="closePartner":
            print ("partner found")
            if not self.force_random:
                targetPoint,weight = self.pickPartner(mingrs,listePartner,currentPos=trans)
                if targetPoint is None :
                    targetPoint = trans
                else :#maybe get the ptid that can have it
                    found = True
                    if self.compNum > 0 :
                        d,i=organelle.OGsrfPtsBht.query(targetPoint)
                        vx, vy, vz = v1 = self.principalVector
                        #surfacePointsNormals problem here
                        v2 = organelle.ogsurfacePointsNormals[i]
                        try :
                            rotMat = numpy.array( rotVectToVect(v1, v2 ), 'f')
                        except :
                            print('PROBLEM ', self.name)
                            rotMat = numpy.identity(4)
                    #find a newpoint here?
                    return targetPoint,rotMat,found
                    x,y,z = targetPoint
                    rad = self.radii[0][0]*2.
                    bb = ( [x-rad, y-rad, z-rad], [x+rad, y+rad, z+rad] )
                    pointsInCube = self.histoVol.grid.getPointsInCube(bb, targetPoint,rad )
                    #is one of this point can receive the current ingredient
                    cut  = rad-self.histoVol.smallestProteinSize
                    for pt in pointsInCube:
                        #get the closest point to object?
                        d = distance[pt]
                        if d>=cut and self.compNum == self.histoVol.grid.gridPtId[pt]:
                            #lets just take the first one
                            targetPoint = self.histoVol.masterGridPositions[pt]
                            break
            else :
                targetPoint = trans 
        else :
            print ("no partner found")
        return targetPoint,rotMat,found

    def pandaBullet_collision(self,pos,rot,rbnode,getnodes=False):
        r=[False]
        liste_nodes=[]
        if len(self.histoVol.rTrans) == 0 : r=[False]
        else :
            closesbody_indice = self.getClosestIngredient(pos,self.histoVol,cutoff=self.histoVol.largestProteinSize+self.encapsulatingRadius*2.0)#vself.radii[0][0]*2.0
            if len(closesbody_indice["indices"]) == 0: r =[False]         #closesbody_indice[0] == -1            
            else : 
                print ("get RB ",len(closesbody_indice["indices"]))
                if rbnode is None :
                    rbnode = self.get_rb_model()
                    self.histoVol.moveRBnode(rbnode,pos, rot )
                    print ("get RB for ",self.name,pos,rot)
                liste_nodes = self.get_rbNodes(closesbody_indice,pos,getInfo=True)
                print ("test collision against  ",len(liste_nodes))
                for node in liste_nodes:
                    print ("collision test with ",node)
                    self.histoVol.moveRBnode(node[0], node[1], node[2])  #Pb here ? 
                    col = (self.histoVol.world.contactTestPair(rbnode, node[0]).getNumContacts() > 0 )
                    r=[col]
                    if col :
                        break
        if getnodes:
            return  True in r,liste_nodes
        else :
            return True in r

    def pandaBullet_place(self, histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None,drop=True,usePP=False):
        """
        drop the ingredient on grid point ptInd
        """
        histoVol.setupPanda()
        afvi = histoVol.afviewer
        rejectionCount = 0
        spacing = histoVol.grid.gridSpacing#/1.1547#  histoVol.smallestProteinSize
        jx, jy, jz = self.jitterMax
        jitter = spacing#histoVol.callFunction(self.getMaxJitter,(spacing,))
        jitter2 = jitter * jitter
        
        if self.compNum == 0 :
            organelle = histoVol
        else :
            organelle = histoVol.compartments[abs(self.compNum)-1]
        compartment = organelle
        compNum = self.compNum
        radius = self.minRadius
        runTimeDisplay = histoVol.runTimeDisplay
        
        gridPointsCoords = histoVol.masterGridPositions
        periodic_pos = None
        # compute rotation matrix rotMat
        if compNum>0 :
            # for surface points we compute the rotation which
            # aligns the principalVector with the surface normal
            # no noise ?
            vx, vy, vz = v1 = self.principalVector
            #surfacePointsNormals problem here
            v2 = organelle.surfacePointsNormals[ptInd]
            try :
                rotMat = numpy.array( rotVectToVect(v1, v2 ), 'f')
            except :
                print('PROBLEM ', self.name)
                rotMat = numpy.identity(4)
        else:
#            self.useOrientBias = True
            if self.useRotAxis :
#                angle = random()*6.2831#math.radians(random()*360.)#random()*pi*2.
#                print "angle",angle,math.degrees(angle)
#                direction = self.rotAxis
                if sum(self.rotAxis) == 0.0 :
                    rotMat=numpy.identity(4)
                elif self.useOrientBias and self.packingMode =="gradient":#you need a gradient here
                    rotMat=self.alignRotation(gridPointsCoords[ptInd] )
#                    print("rotAxis = ", self.rotAxis)
#                    print("rotRange = ", self.rotRange)
#                    print("rotMat aligned = ", rotMatAligned)
#                    print("**************************************************Rotate GRADIENT ON **************************")
#                    rotMatRand = autopack.helper.rotation_matrix(random()*self.rotRange,self.rotAxis)
#                    print("rotMat random = ", rotMatRand)
#                    
#                    print("rotMat biased = ", rotMat)
                else :
                    rotMat=autopack.helper.rotation_matrix(random()*self.rotRange,self.rotAxis)
            # for other points we get a random rotation
            else :
                rotMat=histoVol.randomRot.get()

#        if verbose :
#            pass#print('%s nbs:%2d'%(self.pdb, len(self.positions[0])), end=' ')

        # jitter position loop
        jitterList = []
        collD1 = []
        collD2 = []

        trans = gridPointsCoords[ptInd] # drop point, surface points.
        targetPoint = trans
        moving = None
        if runTimeDisplay and self.mesh:
            if hasattr(self,"mesh_3d"):
                #create an instance of mesh3d and place it
                name = self.name + str(ptInd)
                moving = afvi.vi.getObject(name)
                if moving is None :
                    if self.mesh_3d is None :
                        moving = afvi.vi.Sphere(name,radius=self.radii[0][0],
                                                        parent=afvi.staticMesh)[0]
                        afvi.vi.setTranslation(moving,pos=targetPoint)
                    else :
                        moving=  afvi.vi.newInstance(name,self.mesh_3d,#.GetDown(),
                                                    matrice=rotMat,
                                                    location=targetPoint, parent = afvi.staticMesh)
                else :   
                    #afvi.vi.setTranslation(moving,pos=targetPoint)#rot?
                    mat = rotMat.copy()
                    mat[:3, 3] = targetPoint
                    afvi.vi.setObjectMatrix(moving,mat,transpose=True)
#                afvi.vi.update()  #Graham Killed this unneccesarry redraw
#            print ('Must check for collision here before accepting final position')
        #do we get the list of neighbours first > and give a different trans...closer to the partner
        #we should look up for an available ptID around the picked partner if any
        #getListCloseIngredient
        #should se a distance_of_influence ? or self.histoVol.largestProteinSize+self.encapsulatingRadius*2.0
        #or the grid diagonal
        #we need to change here in case tilling, the pos,rot ade deduced fromte tilling.
        if self.packingMode[-4:] == 'tile' :
            if self.tilling == None :
                self.setTilling(compartment)
            if self.counter!=0 : 
                #pick the next Hexa pos/rot.
                t, r = self.tilling.getNextHexaPosRot()
                if len(t) :
                    trans = t
                    rotMat = r
                    targetPoint = trans
                    if runTimeDisplay and self.mesh:
                        mat = rotMat.copy()
                        mat[:3, 3] = targetPoint
                        afvi.vi.setObjectMatrix(moving,mat,transpose=True)
                        afvi.vi.update()                                                   
                else :
                    self.reject()                
                    return False, nbFreePoints#,targetPoint, rotMat  
            else :
                self.tilling.init_seed(histoVol.seed_used)
        elif histoVol.ingrLookForNeighbours and self.packingMode == "closePartner":
            bind = True
            print ("look for ingredient",trans)
            #roll a dice about proba_not_binding
            if self.proba_not_binding != 0 :#between 0 and 1
                b=random()
                if b <= self.proba_not_binding :
                    bind = False
            if bind :
                closesbody_indice = self.getClosestIngredient(trans,self.histoVol,
                    cutoff=self.histoVol.grid.diag)#vself.radii[0][0]*2.0
                #return R[indice] and distance R["distances"] 
                targetPoint,rotMat,found = self.lookForNeighbours(trans,rotMat,organelle,afvi,distance,
                                                     closest_indice=closesbody_indice)
                if not found and self.counter!=0:
                    self.reject()                
                    return False, nbFreePoints#,targetPoint, rotMat       
    
                #if partner:pickNewPoit like in fill3
                if runTimeDisplay and self.mesh:
                    mat = rotMat.copy()
                    mat[:3, 3] = targetPoint
                    afvi.vi.setObjectMatrix(moving,mat,transpose=True)
                    afvi.vi.update()                                                   
        tx, ty, tz = jtrans = targetPoint
        gridDropPoint = targetPoint        
        #we may increase the jitter, or pick from xyz->Id free for its radius
        #create the rb only once and not at ever jitter
        #rbnode = histoVol.callFunction(self.histoVol.addRB,(self, jtrans, rotMat,),{"rtype":self.Type},)
        #jitter loop
        t1 = time()
        collision2 = False
        for jitterPos in range(self.nbJitter):  #  This expensive Gauusian rejection system should not be the default should it?
            collision2 = False
            # jitter points location
#            print ("jitter ",jitterPos,self.nbJitter)
            if jitter2 > 0.0:
                found = False
                while not found:
#                    dx = jx*jitter*gauss(0., 0.3)
#                    dy = jy*jitter*gauss(0., 0.3)
#                    dz = jz*jitter*gauss(0., 0.3)
                    dx = jx*jitter*uniform(-1.0, 1.0)
                    dy = jy*jitter*uniform(-1.0, 1.0)
                    dz = jz*jitter*uniform(-1.0, 1.0)
                    d2 = dx*dx + dy*dy + dz*dz
                    if d2 < jitter2:
                        if compNum > 0: # jitter less among normal
                            #if self.name=='2uuh C4 SYNTHASE':
                            #    import pdb
                            #    pdb.set_trace()
                            dx, dy, dz, dum = numpy.dot(rotMat, (dx,dy,dz,0))
                        jtrans = (tx+dx, ty+dy, tz+dz)
                        found = True
#                    else:  #Graham Killed thse unnecessary updates
#                        if verbose :  #Graham Killed thse unnecessary updates
#                            print('JITTER REJECTED', d2, jitter2)  #Graham Killed thse unnecessary updates
#                    if runTimeDisplay and moving is not None :  #Graham Killed thse unnecessary updates
#                        afvi.vi.setTranslation(moving,pos=jtrans)   
#                        afvi.vi.update()  
#                        print "ok moving"
            else:
                jtrans = targetPoint
                dx = dy = dz = 0.0
           
            histoVol.totnbJitter += 1
            histoVol.jitterLength += dx*dx + dy*dy + dz*dz
            jitterList.append( (dx,dy,dz) )
            
            #print 'j%d %.2f,%.2f,%.2f,'%(jitterPos,tx, ty, tz),
#            if verbose :
#                print('j%d'%jitterPos)# end=' '
            
            # loop over all spheres representing ingredient
            modSphNum = 1
            if sphGeom is not None:
                modCent = []
                modRad = []

            ## check for collisions 
            ## 
            level = self.collisionLevel

            # randomize rotation about axis
            if compNum>0:
                rotMatj = self.getAxisRotation(rotMat)
            else:
                if self.useRotAxis :
                    if sum(self.rotAxis) == 0.0 :
                        rotMatj=numpy.identity(4)
                    elif self.useOrientBias and self.packingMode =="gradient":
#                        rotMatj = self.getAxisRotation(rotMat)
                        rotMatj = self.getBiasedRotation(rotMat,weight = None)
#                            weight = 1.0 - self.histoVol.gradients[self.gradient].weight[ptInd])
                    else :
                        rotMatj=autopack.helper.rotation_matrix(random()*self.rotRange,self.rotAxis)
                else :
                    rotMatj=histoVol.randomRot.get()  #Graham turned this back on to replace rotMat.copy() so ing rotate each time
#                    rotMatj = rotMat.copy()
            if self.packingMode[-4:] == 'tile' :
                if self.counter==0 : 
                    euler = euler_from_matrix(rotMat)
                    print (euler)
#                    rotMat = euler_matrix(euler[0],euler[1],0).transpose()
                jtrans = targetPoint
                rotMatj =  rotMat[:]#self.tilling.getNextHexaPosRot()
            if runTimeDisplay and moving is not None :
#                print "ok rot copy"
                mat = rotMatj.copy()
                mat[:3, 3] = jtrans
                afvi.vi.setObjectMatrix(moving,mat,transpose=True)
#                afvi.vi.setTranslation(moving,pos=jtrans)
                afvi.vi.update()
#            closeS = self.checkPointSurface(jtrans,cutoff=float(self.cutoff_surface))
#            if closeS :
#                print ("ok reject once")
#                self.rejectOnce(None,moving,afvi)
#                continue
            r=[False]
            rbnode = self.get_rb_model()
#            print ("periodicity ?",getNormedVectorU(self.jitterMax),self.encapsulatingRadius)
            #cutoff ?
            periodic_pos = self.histoVol.grid.getPositionPeridocity(jtrans,
                    getNormedVectorOnes(self.jitterMax),self.encapsulatingRadius)                
            histoVol.callFunction(histoVol.moveRBnode,(rbnode, jtrans, rotMatj,))
            perdiodic_collision = False
            if periodic_pos is not None and self.packingMode !="gradient" :
                print ("OK Periodicity ",len(periodic_pos),periodic_pos)
                for p in periodic_pos :
                    histoVol.callFunction(histoVol.moveRBnode,(rbnode, p, rotMatj,))
                    perdiodic_collision = self.pandaBullet_collision(p,rotMatj,rbnode)                
                    r.extend([perdiodic_collision])
                    if  True in r : break
                    histoVol.callFunction(histoVol.moveRBnode,(rbnode, jtrans, rotMatj,))              
                    rbnode2 = self.get_rb_model(alt=True)
                    self.histoVol.moveRBnode(rbnode2, p, rotMatj)  #Pb here ? 
                    col = (self.histoVol.world.contactTestPair(rbnode, rbnode2).getNumContacts() > 0 )
#                    perdiodic_collision.extend([col])
                    r.extend([col])# = True in perdiodic_collision
                    if runTimeDisplay and moving is not None :
                        mat = rotMatj.copy()
                        mat[:3, 3] = jtrans
                        afvi.vi.setObjectMatrix(moving,mat,transpose=True)
                        afvi.vi.update()
                    if True in r : break
#            rbnode = histoVol.callFunction(self.histoVol.addRB,(self, jtrans, rotMatj,),{"rtype":self.Type},)
#            rbnode = histoVol.callFunction(self.histoVol.addRB,(self, jtrans, rotMatj,),{"rtype":self.Type},            
            t=time()   
            test=self.testPoint(jtrans)
#            print ("testPoints",r,test,jtrans)
            #       checkif rb collide 
#            r=[ (self.histoVol.world.contactTestPair(rbnode, n).getNumContacts() > 0 ) for n in self.histoVol.static]  
            if not test and not ( True in r):
                #check for close surface ?
                #closeS = self.checkPointSurface(newPt,cutoff=self.cutoff_surface)
#                raw_input()
#                if self.histoVol.close_ingr_bhtree.FreePts.NumPts == 0 : r=[False]
                if len(self.histoVol.rTrans) == 0 : r=[False]
                else :
#                    print("getClosestIngredient",jtrans)
                    closesbody_indice = self.getClosestIngredient(jtrans,self.histoVol,cutoff=self.histoVol.largestProteinSize+self.encapsulatingRadius*2.0)#vself.radii[0][0]*2.0
#                    print ("len(closesbody_indice) ",len(closesbody_indice["indices"]),str(self.histoVol.largestProteinSize+self.encapsulatingRadius) )                    
                    if len(closesbody_indice["indices"]) == 0: r =[False]         #closesbody_indice[0] == -1            
                    else : 
                        liste_nodes = self.get_rbNodes(closesbody_indice,jtrans,getInfo=True)
#                        print ("len(liste_nodes) ",len(liste_nodes) )
                        if usePP :
                            #use self.grab_cb and self.pp_server
                            ## Divide the task or just submit job
                            n=0
                            self.histoVol.grab_cb.reset()#can't pickle bullet world 
                            print ("usePP ",autopack.ncpus)
                            for i in range(len(liste_nodes)/autopack.ncpus):
                                for c in range(autopack.ncpus):
                                    print ("submit job",i,c,n)
                                    self.histoVol.pp_server.submit(checkCollision_mp,
                                                (self.histoVol.world,rbnode, liste_nodes[n]),
                                            callback=self.histoVol.grab_cb.grab)
#                                    self.histoVol.pp_server.submit(self.histoVol.world.contactTestPair, 
#                                                          (rbnode, liste_nodes[n]), 
#                                            callback=self.histoVol.grab_cb.grab)
                                    n+=1
                                self.histoVol.pp_server.wait()
                                r.extend(self.histoVol.grab_cb.collision[:])
                                if True in r :
                                    break
                        else :
                            #why will it be not woking with organelle ? 
                            #tranformation prolem ?
                            for node in liste_nodes:
                                self.histoVol.moveRBnode(node[0], node[1], node[2])  #Pb here ? 
#                                print (node[0],node[1])
                                col = (self.histoVol.world.contactTestPair(rbnode, node[0]).getNumContacts() > 0 )
#                                print ("collision ?",self.name,node[3].name,col,autopack.helper.measure_distance(jtrans,node[1]))
                                r=[col]
                                if col :
                                    break
                        #r=[ (self.histoVol.world.contactTestPair(rbnode, node).getNumContacts() > 0 ) for node in liste_nodes]# in closesbody_indice if n != len(self.histoVol.rTrans)-1]  #except last one  that should be last drop fragment                        
#                        r=[ (self.histoVol.world.contactTestPair(rbnode, self.histoVol.static[n]).getNumContacts() > 0 ) for n in closesbody_indice if n < len(self.histoVol.static)]  #except last one  that should be last drop fragment
#            else :
#                result2 = self.histoVol.world.contactTest(rbnode)
#                r = [( result2.getNumContacts() > 0),]    
#            print ("contact All ",(True in r), time()-t, len(r),self.encapsulatingRadius*2.,len(self.histoVol.rTrans))                 
#            t=time()     
            collision2=( True in r)
            collisionComp = False
#            print("collide??",collision2,r,test)
#            print ("contactTestPair",collision2,time()-t)
#            print ("contact Pair ",collision, r,self.histoVol.static) #gave nothing ???
            #need to check compartment too
            if not collision2 and not test :# and not collision2:
                print ("no collision")
                if self.compareCompartment:
                    collisionComp = self.compareCompartmentPrimitive(level,
                     jtrans, rotMatj,gridPointsCoords, distance)
                if not collisionComp :
                    #self.update_data_tree(jtrans,rotMatj,ptInd=ptInd)?
                    self.histoVol.static.append(rbnode)
                    self.histoVol.moving = None
                    self.histoVol.rTrans.append(jtrans)
                    self.histoVol.rRot.append(rotMatj)
                    self.histoVol.rIngr.append(self)
                    self.histoVol.result.append([ jtrans, rotMatj, self, ptInd ])
    #                histoVol.close_ingr_bhtree.MoveRBHPoint(histoVol.nb_ingredient,(jtrans[0],jtrans[1],jtrans[2]),1)
#                    self.histoVol.nb_ingredient+=1
                    if periodic_pos is not None and self.packingMode !="gradient" :
                        for p in periodic_pos :
                            self.histoVol.rTrans.append(p)
                            self.histoVol.rRot.append(rotMatj)
                            self.histoVol.rIngr.append(self)
                            self.histoVol.result.append([ p, rotMatj, self, ptInd ])
                            #self.histoVol.nb_ingredient+=1                        
                    if self.histoVol.treemode == "bhtree":# "cKDTree"
                        if len(self.histoVol.rTrans) >= 1 : bhtreelib.freeBHtree(self.histoVol.close_ingr_bhtree)
                        if len(self.histoVol.rTrans) : self.histoVol.close_ingr_bhtree=bhtreelib.BHtree( self.histoVol.rTrans, None, 10)
                    else :
                        #rebuild kdtree
                        if len(self.histoVol.rTrans) >= 1 : self.histoVol.close_ingr_bhtree= spatial.cKDTree(self.histoVol.rTrans, leafsize=10)
                    #self.histoVol.callFunction(self.histoVol.delRB,(rbnode,))                    
                    break # break out of jitter pos loop
                else :
#                    print ("collision")
                    collision2 = collisionComp
#            else :
#                histoVol.callFunction(self.histoVol.delRB,(rbnode,))
                #remove the node
#        if verbose: 
#        print("jitter loop ",time()-t1)
#            self.histoVol.callFunction(self.histoVol.delRB,(rbnode,)) 
        if not collision2 and not test:# and not collision2:
#            print("jtrans for NotCollision= ", jtrans)
            drop = True    
            ## get inside points and update distance
            ##
            # use best sperical approcimation   

            insidePoints = {}
            newDistPoints = {}
            t3=time()
            insidePoints,newDistPoints=self.getDistances(jtrans, rotMatj,
                                             gridPointsCoords, distance,dpad)
            # save dropped ingredient
            if verbose:
                print("compute distance loop ",time()-t3)
            nbFreePoints = histoVol.callFunction(self.updateDistances,(histoVol,insidePoints,
                                                            newDistPoints,
                                                freePoints,nbFreePoints, distance,
                                                histoVol.masterGridPositions, verbose))
            if periodic_pos is not None and self.packingMode !="gradient" :
                for p in periodic_pos :
                    insidePoints,newDistPoints=self.getDistances(p, rotMatj,
                                                     gridPointsCoords, distance,dpad)
                    nbFreePoints = histoVol.callFunction(self.updateDistances,(histoVol,insidePoints,
                                                                    newDistPoints,
                                                        freePoints,nbFreePoints, distance,
                                                        histoVol.masterGridPositions, verbose))
                
            if drop:
                #print "ok drop",organelle.name,self.name
                organelle.molecules.append([ jtrans, rotMatj, self, ptInd ])
                histoVol.order[ptInd]=histoVol.lastrank
                histoVol.lastrank+=1
                if self.packingMode[-4:] == 'tile' :
                    nexthexa = self.tilling.dropTile(self.tilling.idc,
                                            self.tilling.edge_id,jtrans,rotMatj)
                    print ("drop next hexa",nexthexa.name,self.tilling.idc,self.tilling.edge_id)
                    #('drop next hexa', 'hexa_0_', 0, '')
#                if periodic_pos is not None and self.packingMode !="gradient" :
#                    for p in periodic_pos :
#                        organelle.molecules.append([ p, rotMatj, self, -1 ])

            # add one to molecule counter for this ingredient
            self.counter += 1
            self.completion = float(self.counter)/float(self.nbMol)
#            if periodic_pos is not None :
#                for p in periodic_pos :
#                    self.counter += 1
#                    self.completion = float(self.counter)/float(self.nbMol)

            if jitterPos>0:
                histoVol.successfullJitter.append(
                    (self, jitterList, collD1, collD2) )
               
            if verbose :
                print('Success nbfp:%d %d/%d dpad %.2f'%(
                nbFreePoints, self.counter, self.nbMol, dpad))
            if self.name=='in  inside':#<<??<<
                histoVol.jitterVectors.append( (trans, jtrans) )

            success = True
            self.rejectionCounter = 0
            
        else: # got rejected
#            histoVol.callFunction(histoVol.delRB,(rbnode,))
            if runTimeDisplay and moving is not None :
                afvi.vi.deleteObject(moving)
            success = False
            self.haveBeenRejected = True
            if self.packingMode[-4:] == 'tile' :
                if self.tilling.start.nvisit[self.tilling.edge_id] >= 2 :
                    self.tilling.start.free_pos[self.tilling.edge_id]=0

#            histoVol.failedJitter.append(
#                (self, jitterList, collD1, collD2) )#<??

#            distance[ptInd] = max(0, distance[ptInd]*0.9)# ???
            distance[ptInd] = distance[ptInd]-1.0#?
            self.rejectionCounter += 1
            if verbose :
                print('Failed ingr:%s rejections:%d'%(
                self.name, self.rejectionCounter))
            if self.rejectionCounter >= self.rejectionThreshold: #Graham set this to 6000 for figure 13b (Results Fig 3 Test1) otehrwise it fails to fill small guys
                #if verbose :
                print('PREMATURE ENDING of ingredient', self.name)
                self.completion = 1.0
        if drop :
            return success, nbFreePoints
        else :
            return success, nbFreePoints,jtrans, rotMatj
            
    def pandaBullet_place_dev(self, histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None,drop=True):
        """
        drop the ingredient on grid point ptInd
        """
        histoVol.setupPanda()#do I need this everytime?
        histoVol.setupOctree()#do I need this everytime?
        
        if histoVol.panda_solver == "bullet":
            collideFunc = self.histoVol.world.contactTestPair
        elif histoVol.panda_solver == "ode":
            from panda3d.ode import OdeUtil
            collideFunc = OdeUtil.collide
        afvi = histoVol.afviewer
        rejectionCount = 0
        spacing = histoVol.smallestProteinSize
        jx, jy, jz = self.jitterMax
        jitter = histoVol.callFunction(self.getMaxJitter,(spacing,))
        jitter2 = jitter * jitter
        
        if self.compNum == 0 :
            compartment = histoVol
        else :
            compartment = histoVol.compartments[abs(self.compNum)-1]
        compNum = self.compNum
        radius = self.minRadius
        runTimeDisplay = histoVol.runTimeDisplay
        
        gridPointsCoords = histoVol.masterGridPositions

        # compute rotation matrix rotMat
        if compNum>0 :
            # for surface points we compute the rotation which
            # aligns the principalVector with the surface normal
            vx, vy, vz = v1 = self.principalVector
            v2 = compartment.surfacePointsNormals[ptInd]#1000 and it is a dictionary ?
            try :
                rotMat = numpy.array( rotVectToVect(v1, v2 ), 'f')
            except :
                print('PROBLEM ', self.name)
                rotMat = numpy.identity(4)
        else:
            if self.useRotAxis :
#                angle = random()*6.2831#math.radians(random()*360.)#random()*pi*2.
#                print "angle",angle,math.degrees(angle)
#                direction = self.rotAxis
                if sum(self.rotAxis) == 0.0 :
                    rotMat=numpy.identity(4)
                else :
                    rotMat=afvi.vi.rotation_matrix(random()*self.rotRange,self.rotAxis)
            # for other points we get a random rotation
            else :
                rotMat=histoVol.randomRot.get()

#        if verbose :
#            pass#print('%s nbs:%2d'%(self.pdb, len(self.positions[0])), end=' ')

        # jitter position loop
        jitterList = []
        collD1 = []
        collD2 = []

        trans = gridPointsCoords[ptInd] # drop point, surface points.
        targetPoint = trans
        moving = None
        if runTimeDisplay and self.mesh:
            if hasattr(self,"mesh_3d"):
                #create an instance of mesh3d and place it
                name = self.name + str(ptInd)
                moving = afvi.vi.getObject(name)
                if moving is None :
                    if self.mesh_3d is None :
                        moving = afvi.vi.Sphere(name,radius=self.radii[0][0],
                                                        parent=afvi.staticMesh)[0]
                        afvi.vi.setTranslation(moving,pos=targetPoint)
                    else :
                        moving=  afvi.vi.newInstance(name,self.mesh_3d,#.GetDown(),
                                                    matrice=rotMat,
                                                    location=targetPoint, parent = afvi.staticMesh)
                else :   
                    #afvi.vi.setTranslation(moving,pos=targetPoint)#rot?
                    mat = rotMat.copy()
                    mat[:3, 3] = targetPoint
                    afvi.vi.setObjectMatrix(moving,mat,transpose=True)
#                afvi.vi.update()  #Graham Killed this unneccesarry redraw
#            print ('Must check for collision here before accepting final position')
        #do we get the list of neighbours first > and give a different trans...closer to the partner
        #we should look up for an available ptID around the picked partner if any
        #getListPartner
        if histoVol.ingrLookForNeighbours:
            mingrs,listePartner=self.getListePartners(histoVol,trans,rotMat,compartment,afvi)
            #if liste:pickPartner
            if listePartner : #self.packingMode=="closePartner":
#                print "ok partner",len(listePartner)
                if not self.force_random:
                    targetPoint,weight = self.pickPartner(mingrs,listePartner,currentPos=trans)
                    if targetPoint is None :
                        targetPoint = trans
                    else :#maybe get the ptid that can have it
                        #find a newpoint here?
                        x,y,z = targetPoint
                        rad = self.radii[0][0]*2.
                        bb = ( [x-rad, y-rad, z-rad], [x+rad, y+rad, z+rad] )
                        pointsInCube = histoVol.grid.getPointsInCube(bb, targetPoint,rad )
                        #is one of this point can receive the current ingredient
                        cut  = rad-jitter
                        for pt in pointsInCube:
                            d = distance[pt]
                            if d>=cut:
                                #lets just take the first one
                                targetPoint = gridPointsCoords[pt]
                                break
                else :
                    targetPoint = trans
            #if partner:pickNewPoit like in fill3
        tx, ty, tz = jtrans = targetPoint
        gridDropPoint = targetPoint        
        #we may increase the jitter, or pick from xyz->Id free for its radius
        #create the rb only once and not at ever jitter
        rbnode = histoVol.callFunction(self.histoVol.addRB,(self, jtrans, rotMat,),{"rtype":self.Type},)
        #should I create the rbnode on the fly ?
        #ningr_rb = histoVol.callFunction(self.getNeighboursInBox,(histoVol,trans,rotMat,organelle,afvi),{"rb":True})
#        x,y,z=jtrans
#        rad = self.encapsulatingRadius
#        bb=( [x-rad, y-rad, z-rad], [x+rad, y+rad, z+rad] )
        dropedObject = IngredientInstanceDrop(ptInd, jtrans, 
                                                      rotMat, self,rb=rbnode)
        nodes = histoVol.octree.findContainingNodes(dropedObject, histoVol.octree.root)
        ningr_rb = nodes[0].objects
  
        #jitter loop
        t1 = time()
        for jitterPos in range(self.nbJitter):  #  This expensive Gauusian rejection system should not be the default should it?
            # jitter points location
            if jitter2 > 0.0:
                found = False
                while not found:
#                    dx = jx*jitter*gauss(0., 0.3)
#                    dy = jy*jitter*gauss(0., 0.3)
#                    dz = jz*jitter*gauss(0., 0.3)
                    dx = jx*jitter*uniform(-1.0, 1.0)
                    dy = jy*jitter*uniform(-1.0, 1.0)
                    dz = jz*jitter*uniform(-1.0, 1.0)
                    d2 = dx*dx + dy*dy + dz*dz
                    if d2 < jitter2:
                        if compNum > 0: # jitter less among normal
                            #if self.name=='2uuh C4 SYNTHASE':
                            #    import pdb
                            #    pdb.set_trace()
                            dx, dy, dz, dum = numpy.dot(rotMat, (dx,dy,dz,0))
                        jtrans = (tx+dx, ty+dy, tz+dz)
                        found = True
#                    else:  #Graham Killed thse unnecessary updates
#                        if verbose :  #Graham Killed thse unnecessary updates
#                            print('JITTER REJECTED', d2, jitter2)  #Graham Killed thse unnecessary updates
#                    if runTimeDisplay and moving is not None :  #Graham Killed thse unnecessary updates
#                        afvi.vi.setTranslation(moving,pos=jtrans)   
#                        afvi.vi.update()  
#                        print "ok moving"
            else:
                jtrans = targetPoint
                dx = dy = dz = 0.0
           
            histoVol.totnbJitter += 1
            histoVol.jitterLength += dx*dx + dy*dy + dz*dz
            jitterList.append( (dx,dy,dz) )
            
            #print 'j%d %.2f,%.2f,%.2f,'%(jitterPos,tx, ty, tz),
#            if verbose :
#                print('j%d'%jitterPos)# end=' '
            
            # loop over all spheres representing ingredient
            modSphNum = 1
            if sphGeom is not None:
                modCent = []
                modRad = []

            ## check for collisions 
            ## 
            level = self.collisionLevel

            # randomize rotation about axis
            if compNum>0:
                rotMatj = self.getAxisRotation(rotMat)
            else:
                if self.useRotAxis :
                    if sum(self.rotAxis) == 0.0 :
                        rotMatj=numpy.identity(4)
                    else :
                        rotMatj=afvi.vi.rotation_matrix(random()*self.rotRange,self.rotAxis)
                else :
                    rotMatj=histoVol.randomRot.get()  #Graham turned this back on to replace rotMat.copy() so ing rotate each time
#                    rotMatj = rotMat.copy()
            if runTimeDisplay and moving is not None :
#                print "ok rot copy"
                mat = rotMatj.copy()
                mat[:3, 3] = jtrans
                afvi.vi.setObjectMatrix(moving,mat,transpose=True)
#                afvi.vi.setTranslation(moving,pos=jtrans)
                afvi.vi.update()
                
#            rbnode = histoVol.callFunction(self.histoVol.addRB,(self, jtrans, rotMatj,),{"rtype":self.Type},)
            histoVol.callFunction(histoVol.moveRBnode,(rbnode, jtrans, rotMatj,))
            t=time()
            #       checkif rb collide 
#            result2 = self.histoVol.world.contactTest(rbnode)
#            collision = ( result2.getNumContacts() > 0)    
#            print ("contact All ",collision, time()-t, result2.getNumContacts())                 
#            t=time()     
#            ningr_rb = self.getNeighboursInBox(histoVol,trans,rotMat,organelle,afvi,rb=True)
            r=[False]
#            result = self.histoVol.world.contactTest(rbnode).getNumContacts() > 0
#            print ("contactTest find ",result)
#            if not result :
            t1=time()
            #ningr_rb = histoVol.octree.findPosition(histoVol.octree.root, jtrans)
            #ningr_rb = histoVol.octree.findPosition(histoVol.octree.root, jtrans)
            if ningr_rb is not None and len(ningr_rb):
                #ode is just contact
#                print ("get ",len(ningr_rb))
                r=[ (collideFunc(rbnode, n.rigid_body).getNumContacts() > 0 ) for n in ningr_rb]
            #r=[ (self.histoVol.world.contactTestPair(rbnode, n).getNumContacts() > 0 ) for n in self.histoVol.static]  
#                print ("contactTestPair find ",len(r),( True in r))
            #print(("time contactPair", time()-t1))
            collision2=( True in r)
            collisionComp = False

            #print ("contactTestPair",collision2,time()-t,len(r))
            #print ("contact Pair ", len(r),len(ningr_rb)) #gave nothing ???
            #need to check compartment too
                #Graham here:  If this is less expensive (compareCompartment less exp than mesh collision r=) we should do it first. Feb 28, 2013
            if not collision2 :# and not collision2:
                if self.compareCompartment:
                    if self.modelType=='Spheres':
        #                print("running jitter number ", histoVol.totnbJitter, " on Spheres for pt = ", ptInd)
        #                print("jtrans = ", jtrans)
                        collisionComp = histoVol.callFunction(self.checkSphCompart,(
                            self.positions[level], self.radii[level], jtrans, rotMatj,
                            level, gridPointsCoords, distance, histoVol))
        #                print("jitter collision = ", collision, " for pt = ", ptInd, " with jtrans = ", jtrans)
                    elif self.modelType=='Cylinders':
                        collisionComp = histoVol.callFunction(self.checkCylCompart,(
                            self.positions[level], self.positions2[level],
                            self.radii[level], jtrans, rotMatj, gridPointsCoords,
                            distance, histoVol))
                    elif self.modelType=='Cube':
                        collisionComp = histoVol.callFunction(self.checkCubeCompart,(
                            self.positions[0], self.positions2[0], self.radii,
                            jtrans, rotMatj, gridPointsCoords,
                            distance, histoVol))
                if not collisionComp :
                    #self.rbnode[ptInd] = rbnode
                    self.histoVol.static.append(rbnode)
                    self.histoVol.moving = None
#                    if len(self.histoVol.rTrans) > 1 : bhtreelib.freeBHtree(self.histoVol.close_ingr_bhtree)
#                    self.histoVol.close_ingr_bhtree=bhtreelib.BHtree( self.histoVol.rTrans, None, 10)
                    if self.histoVol.treemode == "bhtree":# "cKDTree"
                        if len(self.histoVol.rTrans) > 1 : bhtreelib.freeBHtree(self.histoVol.close_ingr_bhtree)
                        self.histoVol.close_ingr_bhtree=bhtreelib.BHtree( self.histoVol.rTrans, None, 10)
                    else :
                        #rebuild kdtree
                        if len(self.histoVol.rTrans) > 1 :self.histoVol.close_ingr_bhtree= spatial.cKDTree(self.histoVol.rTrans, leafsize=10)

                      #add to the octree
                    break # break out of jitter pos loop
                else :
                    collision2 = collisionComp
#            else :
#                histoVol.callFunction(self.histoVol.delRB,(rbnode,))
                #remove the node
#        if verbose: 
#        print("jitter loop ",time()-t1)
        if not collision2 :# and not collision2:
#            print("jtrans for NotCollision= ", jtrans)
            
            ## get inside points and update distance
            ##
            # use best sperical approcimation   

            insidePoints = {}
            newDistPoints = {}
            t3=time()
            #should be replace by self.getPointInside
            if self.modelType=='Spheres':
                self.centT = centT = self.transformPoints(jtrans, rotMatj, self.positions[level])
                centT = self.centT#self.transformPoints(jtrans, rotMatj, self.positions[-1])
#                for pointsInCube,ptsInSphere in self.distances_temp:
                for radc, posc in zip(self.radii[-1], centT):
#                for i,radc in enumerate(self.radii[-1]):
#                    pointsInCube,ptsInSphere,distA = self.distances_temp[i]
#                 
                    rad = radc + dpad
                    x,y,z = posc
#                    #this have already be done in the checkCollision why doing it again
                    bb = ( [x-rad, y-rad, z-rad], [x+rad, y+rad, z+rad] )
                    pointsInCube = histoVol.callFunction(histoVol.grid.getPointsInCube,
                                                         (bb, posc, rad))
#
                    delta = numpy.take(gridPointsCoords,pointsInCube,0)-posc
                    delta *= delta
                    distA = numpy.sqrt( delta.sum(1) )
                    ptsInSphere = numpy.nonzero(numpy.less_equal(distA, rad))[0]

                    for pti in ptsInSphere:
                        pt = pointsInCube[pti]
                        if pt in insidePoints: continue
                        dist = distA[pti]
                        d = dist-radc
                        if dist < radc:  # point is inside dropped sphere
                            if pt in insidePoints:
                                if d < insidePoints[pt]:
                                    insidePoints[pt] = d
                            else:
                                insidePoints[pt] = d
                        elif d < distance[pt]: # point in region of influence
                            if pt in newDistPoints:
                                if d < newDistPoints[pt]:
                                    newDistPoints[pt] = d
                            else:
                                newDistPoints[pt] = d

            elif self.modelType=='Cylinders':
                cent1T = self.transformPoints(jtrans, rotMatj, self.positions[-1])
                cent2T = self.transformPoints(jtrans, rotMatj, self.positions2[-1])

                for radc, p1, p2 in zip(self.radii[-1], cent1T, cent2T):

                    x1, y1, z1 = p1
                    x2, y2, z2 = p2
                    vx, vy, vz = vect = (x2-x1, y2-y1, z2-z1)
                    lengthsq = vx*vx + vy*vy + vz*vz
                    l = sqrt( lengthsq )
                    cx, cy, cz = posc = x1+vx*.5, y1+vy*.5, z1+vz*.5
                    radt = l + radc + dpad
                    bb = ( [cx-radt, cy-radt, cz-radt], [cx+radt, cy+radt, cz+radt] )
                    pointsInCube = histoVol.callFunction(histoVol.grid.getPointsInCube,
                                                         (bb, posc, radt))

                    pd = numpy.take(gridPointsCoords,pointsInCube,0) - p1
                    dotp = numpy.dot(pd, vect)
                    rad2 = radc*radc
                    d2toP1 = numpy.sum(pd*pd, 1)
                    dsq = d2toP1 - dotp*dotp/lengthsq

                    pd2 = numpy.take(gridPointsCoords,pointsInCube,0) - p2
                    d2toP2 = numpy.sum(pd2*pd2, 1)

                    for pti, pt in enumerate(pointsInCube):
                        if pt in insidePoints: continue

                        if dotp[pti] < 0.0: # outside 1st cap
                            d = sqrt(d2toP1[pti])
                            if d < distance[pt]: # point in region of influence
                                if pt in newDistPoints:
                                    if d < newDistPoints[pt]:
                                        newDistPoints[pt] = d
                                else:
                                    newDistPoints[pt] = d
                        elif dotp[pti] > lengthsq:
                            d = sqrt(d2toP2[pti])
                            if d < distance[pt]: # point in region of influence
                                if pt in newDistPoints:
                                    if d < newDistPoints[pt]:
                                        newDistPoints[pt] = d
                                else:
                                    newDistPoints[pt] = d
                        else:
                            d = sqrt(dsq[pti])-radc
                            if d < 0.:  # point is inside dropped sphere
                                if pt in insidePoints:
                                    if d < insidePoints[pt]:
                                        insidePoints[pt] = d
                                else:
                                    insidePoints[pt] = d
            elif self.modelType=='Cube':
                insidePoints,newDistPoints=self.getDistancesCube(jtrans, rotMatj,gridPointsCoords, distance, histoVol)
            
#            print("Graham.Mira........................insidePoints=", insidePoints,"\n while ptInd = ", ptInd)
            insidePoints[ptInd]=-0.1
#            print("Graham.Mira........AFTER...........insidePoints=", insidePoints,"\n while ptInd = ", ptInd)

            # save dropped ingredient
            if verbose:
                print("compute distance loop ",time()-t3)
            if drop:
                dropedObject = IngredientInstanceDrop(ptInd, jtrans, 
                                                      rotMatj, self,rb=rbnode)
                #r= self.encapsulatingRadius
                r = self.minRadius + histoVol.largestProteinSize + \
                        histoVol.smallestProteinSize + histoVol.windowsSize
#                histoVol.octree.insertNode(histoVol.octree.root, r, 
#                                           histoVol.octree.root, dropedObject) 
                histoVol.octree.insertNode(dropedObject,[histoVol.octree.root]) 
                #print "ok drop",organelle.name,self.name
                organelle.molecules.append([ jtrans, rotMatj, self, ptInd ])
                histoVol.order[ptInd]=histoVol.lastrank
                histoVol.lastrank+=1
#                histoVol.close_ingr_bhtree.MoveRBHPoint(histoVol.nb_ingredient,jtrans,0)
                histoVol.nb_ingredient+=1
#                histoVol.close_ingr_bhtree.InsertRBHPoint((jtrans[0],jtrans[1],jtrans[2]),radius,None,histoVol.nb_ingredient)
                
            # update free points
            nbFreePoints = histoVol.callFunction(self.updateDistances,(histoVol,insidePoints,
                                                            newDistPoints,
                                                freePoints,nbFreePoints, distance,
                                                histoVol.masterGridPositions, verbose))
            
#            distChanges = {}
#            for pt,dist in insidePoints.items():
#                # swap point at ptIndr with last free one
#                try:
#                    ind = freePoints.index(pt)
#                    tmp = freePoints[nbFreePoints] #last one
#                    freePoints[nbFreePoints] = pt
#                    freePoints[ind] = tmp
#                    nbFreePoints -= 1
#                except ValueError: # pt not in list of free points
#                    pass
#                distChanges[pt] = (histoVol.masterGridPositions[pt],
#                                   distance[pt], dist)
#                distance[pt] = dist
#            print "update freepoints loop ",time()-t4
#            t5=time()
#            # update distances
#            for pt,dist in newDistPoints.items():
#                if not insidePoints.has_key(pt):
#                    distChanges[pt] = (histoVol.masterGridPositions[pt],
#                                       distance[pt], dist)
#                    distance[pt] = dist
#            print "update distances loop ",time()-t5
            
            if sphGeom is not None:
                for po1, ra1 in zip(modCent, modRad):
                    sphCenters.append(po1)
                    sphRadii.append(ra1)
                    sphColors.append(self.color)

            if labDistGeom is not None:
                verts = []
                labels = []
                colors = []
                #for po1, d1,d2 in distChanges.values():
                fpts = freePoints
                for i in range(nbFreePoints):
                    pt = fpts[i]
                    verts.append(histoVol.masterGridPositions[pt])
                    labels.append( "%.2f"%distance[pt])                    
#                for pt in freePoints[:nbFreePoints]:
#                    verts.append(histoVol.masterGridPositions[pt])
#                    labels.append( "%.2f"%distance[pt])
                labDistGeom.Set(vertices=verts, labels=labels)
                #materials=colors, inheritMaterial=0)

            # add one to molecule counter for this ingredient
            self.counter += 1
            self.completion = float(self.counter)/float(self.nbMol)

            if jitterPos>0:
                histoVol.successfullJitter.append(
                    (self, jitterList, collD1, collD2) )
               
            if verbose :
                print('Success nbfp:%d %d/%d dpad %.2f'%(
                nbFreePoints, self.counter, self.nbMol, dpad))
            if self.name=='in  inside':
                histoVol.jitterVectors.append( (trans, jtrans) )

            success = True
            self.rejectionCounter = 0
#            histoVol.callFunction(histoVol.delRB,(rbnode,))
        else: # got rejected
            #self.rbnode = None
            histoVol.callFunction(histoVol.delRB,(rbnode,))
            if runTimeDisplay and moving is not None :
                afvi.vi.deleteObject(moving)
            success = False
            self.haveBeenRejected = True
            histoVol.failedJitter.append(
                (self, jitterList, collD1, collD2) )

            distance[ptInd] = max(0, distance[ptInd]*0.9)# ???
            self.rejectionCounter += 1
            if verbose :
                print('Failed ingr:%s rejections:%d'%(
                self.name, self.rejectionCounter))
            if self.rejectionCounter >= self.rejectionThreshold: #Graham set this to 6000 for figure 13b (Results Fig 3 Test1) otehrwise it fails to fill small guys
                #if verbose :
                print('PREMATURE ENDING of ingredient', self.name)
                self.completion = 1.0
#        [histoVol.delRB(n) for n in ningr_rb]
        if sphGeom is not None:
            sphGeom.Set(vertices=sphCenters, radii=sphRadii,
                        materials=sphColors)
            sphGeom.viewer.OneRedraw()
            sphGeom.viewer.update()

        if drop :
            return success, nbFreePoints
        else :
            return success, nbFreePoints,jtrans, rotMatj

    def rapid_place(self, histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None,drop=True,usePP=False):
        """
        drop the ingredient on grid point ptInd
        """
        #histoVol.setupPanda()
        afvi = histoVol.afviewer
        rejectionCount = 0
        spacing = histoVol.grid.gridSpacing#histoVol.smallestProteinSize#or the gridSpace?
        jx, jy, jz = self.jitterMax
        jitter = spacing#histoVol.callFunction(self.getMaxJitter,(spacing,))
        jitter2 = jitter * jitter
        
        if self.compNum == 0 :
            organelle = histoVol
        else :
            organelle = histoVol.compartments[abs(self.compNum)-1]
        compNum = self.compNum
        radius = self.minRadius
        runTimeDisplay = histoVol.runTimeDisplay
        
        gridPointsCoords = histoVol.masterGridPositions

        # compute rotation matrix rotMat
        if compNum>0 :
            # for surface points we compute the rotation which
            # aligns the principalVector with the surface normal
            vx, vy, vz = v1 = self.principalVector
            #surfacePointsNormals problem here
            v2 = organelle.surfacePointsNormals[ptInd]
            try :
                rotMat = numpy.array( rotVectToVect(v1, v2 ), 'f')
            except :
                print('PROBLEM ', self.name)
                rotMat = numpy.identity(4)
        else:
            if self.useRotAxis :
#                angle = random()*6.2831#math.radians(random()*360.)#random()*pi*2.
#                print "angle",angle,math.degrees(angle)
#                direction = self.rotAxis
                if sum(self.rotAxis) == 0.0 :
                    rotMat=numpy.identity(4)
                elif self.useOrientBias and self.packingMode =="gradient":#you need a gradient here
                    rotMat=self.alignRotation(gridPointsCoords[ptInd] )
                else :
                    rotMat=afvi.vi.rotation_matrix(random()*self.rotRange,self.rotAxis)
            # for other points we get a random rotation
            else :
                rotMat=histoVol.randomRot.get()

#        if verbose :
#            pass#print('%s nbs:%2d'%(self.pdb, len(self.positions[0])), end=' ')

        # jitter position loop
        jitterList = []
        collD1 = []
        collD2 = []

        trans = gridPointsCoords[ptInd] # drop point, surface points.
        targetPoint = trans
        moving = None
        if runTimeDisplay and self.mesh:
            if hasattr(self,"mesh_3d"):
                #create an instance of mesh3d and place it
                name = self.name + str(ptInd)
                moving = afvi.vi.getObject(name)
                if moving is None :
                    if self.mesh_3d is None :
                        moving = afvi.vi.Sphere(name,radius=self.radii[0][0],
                                                        parent=afvi.staticMesh)[0]
                        afvi.vi.setTranslation(moving,pos=targetPoint)
                    else :
                        moving=  afvi.vi.newInstance(name,self.mesh_3d,#.GetDown(),
                                                    matrice=rotMat,
                                                    location=targetPoint, parent = afvi.staticMesh)
                else :   
                    #afvi.vi.setTranslation(moving,pos=targetPoint)#rot?
                    mat = rotMat.copy()
                    mat[:3, 3] = targetPoint
                    afvi.vi.setObjectMatrix(moving,mat,transpose=True)
#                afvi.vi.update()  #Graham Killed this unneccesarry redraw
#            print ('Must check for collision here before accepting final position')
        #do we get the list of neighbours first > and give a different trans...closer to the partner
        #we should look up for an available ptID around the picked partner if any
        #getListPartner
        if histoVol.ingrLookForNeighbours:
            closesbody_indice = self.getClosestIngredient(trans,self.histoVol,
                cutoff=self.histoVol.grid.diag)#vself.radii[0][0]*2.0
            targetPoint,rotMat = self.lookForNeighbours(trans,rotMat,organelle,afvi,distance,
                                                 closest_indice=closesbody_indice)
        tx, ty, tz = jtrans = targetPoint
        gridDropPoint = targetPoint        
        #we may increase the jitter, or pick from xyz->Id free for its radius
        #create the rb only once and not at ever jitter
        #rbnode = histoVol.callFunction(self.histoVol.addRB,(self, jtrans, rotMat,),{"rtype":self.Type},)
        #jitter loop
        t1 = time()
        collision2 = False
        for jitterPos in range(self.nbJitter):  #  This expensive Gauusian rejection system should not be the default should it?
            collision2 = False
            # jitter points location
            if jitter2 > 0.0:
                found = False
                while not found:
#                    dx = jx*jitter*gauss(0., 0.3)
#                    dy = jy*jitter*gauss(0., 0.3)
#                    dz = jz*jitter*gauss(0., 0.3)
                    dx = jx*jitter*uniform(-1.0, 1.0)
                    dy = jy*jitter*uniform(-1.0, 1.0)
                    dz = jz*jitter*uniform(-1.0, 1.0)
                    d2 = dx*dx + dy*dy + dz*dz
                    if d2 < jitter2:
                        if compNum > 0: # jitter less among normal
                            #if self.name=='2uuh C4 SYNTHASE':
                            #    import pdb
                            #    pdb.set_trace()
                            dx, dy, dz, dum = numpy.dot(rotMat, (dx,dy,dz,0))
                        jtrans = (tx+dx, ty+dy, tz+dz)
                        found = True
#                    else:  #Graham Killed thse unnecessary updates
#                        if verbose :  #Graham Killed thse unnecessary updates
#                            print('JITTER REJECTED', d2, jitter2)  #Graham Killed thse unnecessary updates
#                    if runTimeDisplay and moving is not None :  #Graham Killed thse unnecessary updates
#                        afvi.vi.setTranslation(moving,pos=jtrans)   
#                        afvi.vi.update()  
#                        print "ok moving"
            else:
                jtrans = targetPoint
                dx = dy = dz = 0.0
           
            histoVol.totnbJitter += 1
            histoVol.jitterLength += dx*dx + dy*dy + dz*dz
            jitterList.append( (dx,dy,dz) )
            
            #print 'j%d %.2f,%.2f,%.2f,'%(jitterPos,tx, ty, tz),
#            if verbose :
#                print('j%d'%jitterPos)# end=' '
            
            # loop over all spheres representing ingredient
            modSphNum = 1
            if sphGeom is not None:
                modCent = []
                modRad = []

            ## check for collisions 
            ## 
            level = self.collisionLevel

            # randomize rotation about axis
            if compNum>0:
                rotMatj = self.getAxisRotation(rotMat)
            else:
                if self.useRotAxis :
                    if sum(self.rotAxis) == 0.0 :
                        rotMatj=numpy.identity(4)
                    elif self.useOrientBias and self.packingMode =="gradient":
#                        rotMatj = self.getAxisRotation(rotMat)
                        rotMatj = self.getBiasedRotation(rotMat,weight = None)
#                            weight = 1.0 - self.histoVol.gradients[self.gradient].weight[ptInd])
                    else :
                        rotMatj=afvi.vi.rotation_matrix(random()*self.rotRange,self.rotAxis)
                else :
                    rotMatj=histoVol.randomRot.get()  #Graham turned this back on to replace rotMat.copy() so ing rotate each time
#                    rotMatj = rotMat.copy()
            if runTimeDisplay and moving is not None :
#                print "ok rot copy"
                mat = rotMatj.copy()
                mat[:3, 3] = jtrans
                afvi.vi.setObjectMatrix(moving,mat,transpose=True)
#                afvi.vi.setTranslation(moving,pos=jtrans)
                afvi.vi.update()
            
#            closeS = self.checkPointSurface(jtrans,cutoff=float(self.cutoff_surface))
#            if closeS :
#                self.rejectOnce(None,moving,afvi)
#                continue
            r=False#self.testPoint(jtrans)
            if not r :
                collision2,liste_nodes = self.collision_rapid(jtrans,rotMatj,usePP=usePP)
                collisionComp = False
            if autopack.verbose:
                print ("collision",collision2,r)
#            print ("contact Pair ",collision, r,self.histoVol.static) #gave nothing ???
            #need to check compartment too
            if not collision2 and not r:# and not collision2:
                if self.compareCompartment:
                    if self.modelType=='Spheres':
        #                print("running jitter number ", histoVol.totnbJitter, " on Spheres for pt = ", ptInd)
        #                print("jtrans = ", jtrans)
                        collisionComp = histoVol.callFunction(self.checkSphCompart,(
                            self.positions[level], self.radii[level], jtrans, rotMatj,
                            level, gridPointsCoords, distance, histoVol))
        #                print("jitter collision = ", collision, " for pt = ", ptInd, " with jtrans = ", jtrans)
                    elif self.modelType=='Cylinders':
                        collisionComp = histoVol.callFunction(self.checkCylCompart,(
                            self.positions[level], self.positions2[level],
                            self.radii[level], jtrans, rotMatj, gridPointsCoords,
                            distance, histoVol))
                    elif self.modelType=='Cube':
                        collisionComp = histoVol.callFunction(self.checkCubeCompart,(
                            self.positions[0], self.positions2[0], self.radii,
                            jtrans, rotMatj, gridPointsCoords,
                            distance, histoVol))
                if not collisionComp :
                    #update_data_tree
                    self.update_data_tree(jtrans,rotMatj,ptInd=ptInd)
#                    self.histoVol.static.append(rbnode)
#                    self.histoVol.moving = None
#                    self.histoVol.rTrans.append(jtrans)
#                    self.histoVol.rRot.append(rotMatj)
#                    self.histoVol.rIngr.append(self)
#                    self.histoVol.result.append([ jtrans, rotMatj, self, ptInd ])
#    #                histoVol.close_ingr_bhtree.MoveRBHPoint(histoVol.nb_ingredient,(jtrans[0],jtrans[1],jtrans[2]),1)
#                    self.histoVol.nb_ingredient+=1
##                    self.histoVol.close_ingr_bhtree.InsertRBHPoint((jtrans[0],jtrans[1],jtrans[2]),radius,None,histoVol.nb_ingredient)
#                    #update bhree
##                    if len(self.histoVol.rTrans) > 1 : bhtreelib.freeBHtree(self.histoVol.close_ingr_bhtree)
##                    if len(self.histoVol.rTrans) : self.histoVol.close_ingr_bhtree=bhtreelib.BHtree( self.histoVol.rTrans, None, 10)
#                    if self.histoVol.treemode == "bhtree":# "cKDTree"
#                        if len(self.histoVol.rTrans) > 1 : bhtreelib.freeBHtree(self.histoVol.close_ingr_bhtree)
#                        if len(self.histoVol.rTrans) : self.histoVol.close_ingr_bhtree=bhtreelib.BHtree( self.histoVol.rTrans, None, 10)
#                    else :
#                        #rebuild kdtree
#                        if len(self.histoVol.rTrans) > 1 :self.histoVol.close_ingr_bhtree= spatial.cKDTree(self.histoVol.rTrans, leafsize=10)
#                    #self.histoVol.callFunction(self.histoVol.delRB,(rbnode,))                    
                    break # break out of jitter pos loop
                else :
                    collision2 = collisionComp
#            else :
#                histoVol.callFunction(self.histoVol.delRB,(rbnode,))
                #remove the node
#        if verbose: 
#        print("jitter loop ",time()-t1)
#            self.histoVol.callFunction(self.histoVol.delRB,(rbnode,)) 
        if not collision2 and not r:# and not collision2:
            insidePoints = {}
            newDistPoints = {}
            t3=time()
            #should be replace by self.getPointInside
            if self.modelType=='Spheres':
                self.centT = centT = self.transformPoints(jtrans, rotMatj, self.positions[level])
                centT = self.centT#self.transformPoints(jtrans, rotMatj, self.positions[-1])
#                for pointsInCube,ptsInSphere in self.distances_temp:
                for radc, posc in zip(self.radii[-1], centT):
#                for i,radc in enumerate(self.radii[-1]):
#                    pointsInCube,ptsInSphere,distA = self.distances_temp[i]
#                 
                    rad = radc + dpad
                    x,y,z = posc
#                    #this have already be done in the checkCollision why doing it again
                    bb = ( [x-rad, y-rad, z-rad], [x+rad, y+rad, z+rad] )
                    pointsInCube = histoVol.callFunction(histoVol.grid.getPointsInCube,
                                                         (bb, posc, rad))
#
                    delta = numpy.take(gridPointsCoords,pointsInCube,0)-posc
                    delta *= delta
                    distA = numpy.sqrt( delta.sum(1) )
                    ptsInSphere = numpy.nonzero(numpy.less_equal(distA, rad))[0]

                    for pti in ptsInSphere:
                        pt = pointsInCube[pti]
                        if pt in insidePoints: continue
                        dist = distA[pti]
                        d = dist-radc
                        if dist < radc:  # point is inside dropped sphere
                            if pt in insidePoints:
                                if d < insidePoints[pt]:
                                    insidePoints[pt] = d
                            else:
                                insidePoints[pt] = d
                        elif d < distance[pt]: # point in region of influence
                            if pt in newDistPoints:
                                if d < newDistPoints[pt]:
                                    newDistPoints[pt] = d
                            else:
                                newDistPoints[pt] = d

            elif self.modelType=='Cylinders':
                cent1T = self.transformPoints(jtrans, rotMatj, self.positions[-1])
                cent2T = self.transformPoints(jtrans, rotMatj, self.positions2[-1])

                for radc, p1, p2 in zip(self.radii[-1], cent1T, cent2T):

                    x1, y1, z1 = p1
                    x2, y2, z2 = p2
                    vx, vy, vz = vect = (x2-x1, y2-y1, z2-z1)
                    lengthsq = vx*vx + vy*vy + vz*vz
                    l = sqrt( lengthsq )
                    cx, cy, cz = posc = x1+vx*.5, y1+vy*.5, z1+vz*.5
                    radt = l + radc + dpad
                    bb = ( [cx-radt, cy-radt, cz-radt], [cx+radt, cy+radt, cz+radt] )
                    pointsInCube = histoVol.callFunction(histoVol.grid.getPointsInCube,
                                                         (bb, posc, radt))

                    pd = numpy.take(gridPointsCoords,pointsInCube,0) - p1
                    dotp = numpy.dot(pd, vect)
                    rad2 = radc*radc
                    d2toP1 = numpy.sum(pd*pd, 1)
                    dsq = d2toP1 - dotp*dotp/lengthsq

                    pd2 = numpy.take(gridPointsCoords,pointsInCube,0) - p2
                    d2toP2 = numpy.sum(pd2*pd2, 1)

                    for pti, pt in enumerate(pointsInCube):
                        if pt in insidePoints: continue

                        if dotp[pti] < 0.0: # outside 1st cap
                            d = sqrt(d2toP1[pti])
                            if d < distance[pt]: # point in region of influence
                                if pt in newDistPoints:
                                    if d < newDistPoints[pt]:
                                        newDistPoints[pt] = d
                                else:
                                    newDistPoints[pt] = d
                        elif dotp[pti] > lengthsq:
                            d = sqrt(d2toP2[pti])
                            if d < distance[pt]: # point in region of influence
                                if pt in newDistPoints:
                                    if d < newDistPoints[pt]:
                                        newDistPoints[pt] = d
                                else:
                                    newDistPoints[pt] = d
                        else:
                            d = sqrt(dsq[pti])-radc
                            if d < 0.:  # point is inside dropped sphere
                                if pt in insidePoints:
                                    if d < insidePoints[pt]:
                                        insidePoints[pt] = d
                                else:
                                    insidePoints[pt] = d
            elif self.modelType=='Cube':
                insidePoints,newDistPoints=self.getDistancesCube(jtrans, rotMatj,gridPointsCoords, distance, histoVol)
                
            # save dropped ingredient
            if verbose:
                print("compute distance loop ",time()-t3)
            if drop:
                #print "ok drop",organelle.name,self.name
                organelle.molecules.append([ jtrans, rotMatj, self, ptInd ])
                histoVol.order[ptInd]=histoVol.lastrank
                histoVol.lastrank+=1
#                histoVol.rTrans.append(jtrans)
##                histoVol.close_ingr_bhtree.MoveRBHPoint(histoVol.nb_ingredient,(jtrans[0],jtrans[1],jtrans[2]),1)
#                histoVol.nb_ingredient+=1
##                histoVol.close_ingr_bhtree.InsertRBHPoint((jtrans[0],jtrans[1],jtrans[2]),radius,None,histoVol.nb_ingredient)
#                print ("RBHTree", histoVol.close_ingr_bhtree.FreePts.NumPts) 
            # update free points
            nbFreePoints = histoVol.callFunction(self.updateDistances,(histoVol,insidePoints,
                                                            newDistPoints,
                                                freePoints,nbFreePoints, distance,
                                                histoVol.masterGridPositions, verbose))
            
#            distChanges = {}
#            for pt,dist in insidePoints.items():
#                # swap point at ptIndr with last free one
#                try:
#                    ind = freePoints.index(pt)
#                    tmp = freePoints[nbFreePoints] #last one
#                    freePoints[nbFreePoints] = pt
#                    freePoints[ind] = tmp
#                    nbFreePoints -= 1
#                except ValueError: # pt not in list of free points
#                    pass
#                distChanges[pt] = (histoVol.masterGridPositions[pt],
#                                   distance[pt], dist)
#                distance[pt] = dist
#            print "update freepoints loop ",time()-t4
#            t5=time()
#            # update distances
#            for pt,dist in newDistPoints.items():
#                if not insidePoints.has_key(pt):
#                    distChanges[pt] = (histoVol.masterGridPositions[pt],
#                                       distance[pt], dist)
#                    distance[pt] = dist
#            print "update distances loop ",time()-t5
            
            if sphGeom is not None:
                for po1, ra1 in zip(modCent, modRad):
                    sphCenters.append(po1)
                    sphRadii.append(ra1)
                    sphColors.append(self.color)

            if labDistGeom is not None:
                verts = []
                labels = []
                colors = []
                #for po1, d1,d2 in distChanges.values():
                fpts = freePoints
                for i in range(nbFreePoints):
                    pt = fpts[i]
                    verts.append(histoVol.masterGridPositions[pt])
                    labels.append( "%.2f"%distance[pt])                    
#                for pt in freePoints[:nbFreePoints]:
#                    verts.append(histoVol.masterGridPositions[pt])
#                    labels.append( "%.2f"%distance[pt])
                labDistGeom.Set(vertices=verts, labels=labels)
                #materials=colors, inheritMaterial=0)

            # add one to molecule counter for this ingredient
            self.counter += 1
            self.completion = float(self.counter)/float(self.nbMol)

            if jitterPos>0:
                histoVol.successfullJitter.append(
                    (self, jitterList, collD1, collD2) )
               
            if verbose :
                print('Success nbfp:%d %d/%d dpad %.2f'%(
                nbFreePoints, self.counter, self.nbMol, dpad))
            if self.name=='in  inside':
                histoVol.jitterVectors.append( (trans, jtrans) )

            success = True
            self.rejectionCounter = 0
            
        else: # got rejected
            #histoVol.callFunction(histoVol.delRB,(rbnode,))
            if runTimeDisplay and moving is not None :
                afvi.vi.deleteObject(moving)
            success = False
            self.haveBeenRejected = True
            histoVol.failedJitter.append(
                (self, jitterList, collD1, collD2) )

            if (self.encapsulatingRadius - self.histoVol.smallestProteinSize) == 0:
                distance[ptInd] = distance[ptInd] - 10.0
#            distance[ptInd] = distance[ptInd]-5.0# ???
            else :
                distance[ptInd] = max(0, distance[ptInd]*0.9)# ???
                
            print ("distance reduce ",distance[ptInd] )
            self.rejectionCounter += 1
            if verbose :
                print('Failed ingr:%s rejections:%d'%(
                self.name, self.rejectionCounter))
            if self.rejectionCounter >= self.rejectionThreshold: #Graham set this to 6000 for figure 13b (Results Fig 3 Test1) otehrwise it fails to fill small guys
                #if verbose :
                print('PREMATURE ENDING of ingredient', self.name)
                self.completion = 1.0

        if sphGeom is not None:
            sphGeom.Set(vertices=sphCenters, radii=sphRadii,
                        materials=sphColors)
            sphGeom.viewer.OneRedraw()
            sphGeom.viewer.update()

        if drop :
            return success, nbFreePoints
        else :
            return success, nbFreePoints,jtrans, rotMatj

    def collision_jitter(self,jtrans, rotMatj,
                level, gridPointsCoords, distance, histoVol):
        collision=False
        if self.modelType=='Spheres':
#                print("running jitter number ", histoVol.totnbJitter, " on Spheres for pt = ", ptInd)
#                print("jtrans = ", jtrans)
            collision = histoVol.callFunction(self.checkSphCollisions,(
                self.positions[level], self.radii[level], jtrans, rotMatj,
                level, gridPointsCoords, distance, histoVol))
#                print("jitter collision = ", collision, " for pt = ", ptInd, " with jtrans = ", jtrans)
        elif self.modelType=='Cylinders':
            collision = histoVol.callFunction(self.checkCylCollisions,(
                self.positions[level], self.positions2[level],
                self.radii[level], jtrans, rotMatj, gridPointsCoords,
                 distance, histoVol))
        elif self.modelType=='Cube':
            collision = histoVol.callFunction(self.checkCubeCollisions,(
                self.positions[0], self.positions2[0], self.radii,
                jtrans, rotMatj, gridPointsCoords,
                distance, histoVol))  
        return collision

    def collision_rapid(self,jtrans,rotMatj,cutoff=None,usePP=False,point=None,
                        prevpoint=None,liste_nodes=None):
        r=[False]
        #liste_nodes=[]
        if cutoff is None :
            cutoff = self.histoVol.largestProteinSize+self.encapsulatingRadius*2.0
        rbnode = self.get_rapid_model()  
        #liste_nodes=None          
        if len(self.histoVol.rTrans) == 0 : 
            print ("no history rTrans")
            r=[False]
        else :
            if liste_nodes is None :
                if point is not None :
                    closesbody_indice = self.getClosestIngredient(point,self.histoVol,cutoff=cutoff)
                else  :               
                    closesbody_indice = self.getClosestIngredient(jtrans,self.histoVol,cutoff=cutoff)
                if len(closesbody_indice) == 0:
                    #print ("no closesbody_indice")
                    r =[False]         #closesbody_indice[0] == -1            
                else : 
                    #print ("found closesbody_indice",len(closesbody_indice))
                    liste_nodes = self.get_rapid_nodes(closesbody_indice,jtrans,#
                                                       prevpoint=prevpoint)
            #print ("test against",len(liste_nodes),"nodes")
            r=[]
            if usePP :
                n=0
                self.histoVol.grab_cb.reset()
                inputp={}
                for c in range(autopack.ncpus): inputp[c]=[]
                while n < len(liste_nodes):
                    for c in range(autopack.ncpus):
                        if n == len(liste_nodes) : break
                        v1=self.vertices
                        f1=self.faces
                        v2=liste_nodes[n][3].vertices
                        f2=liste_nodes[n][3].faces  
                        inp=(v1,f1,numpy.array(rotMatj[:3,:3],'f'),numpy.array(jtrans,'f'),
                         v2,f2,liste_nodes[n][2],liste_nodes[n][1],liste_nodes[n][3].name)                           
                        inputp[c].append(inp)
                        n+=1
                jobs=[]                                                                                
                for c in range(autopack.ncpus):
                    if not len(inputp[c]) : continue                                
                    j=self.histoVol.pp_server.submit(rapid_checkCollision_rmp,(inputp[c],),
                                callback=self.histoVol.grab_cb.grab,
                                modules=("numpy",))
                    jobs.append(j)
                self.histoVol.pp_server.wait()
                r.extend(self.histoVol.grab_cb.collision[:])                             
            else :    
#                    print jtrans, rotMatj                                        
                for node,trans,rot,ingr in liste_nodes:
#                        print node,trans,rot,ingr
                    RAPIDlib.RAPID_Collide_scaled(numpy.array(rotMatj[:3,:3],'f'), 
                                                  numpy.array(jtrans,'f'), 1.0,
                                        rbnode, rot, trans, 1.0,node,
                                        RAPIDlib.cvar.RAPID_FIRST_CONTACT);
                    collision2 = RAPIDlib.cvar.RAPID_num_contacts!=0
                    r.append(collision2)
                    if collision2 :                        
                        break
        return True in r,liste_nodes

    def pandaBullet_relax(self, histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None,drop=True):
        """
        drop the ingredient on grid point ptInd
        """
        histoVol.setupPanda()
        self.vi = histoVol.afviewer.vi
        afvi = histoVol.afviewer
        windowsSize = histoVol.windowsSize
        simulationTimes = histoVol.simulationTimes
        runTimeDisplay = histoVol.runTimeDisplay
        springOptions = histoVol.springOptions
        self.histoVol = histoVol
        rejectionCount = 0
        spacing = histoVol.smallestProteinSize
        jx, jy, jz = self.jitterMax
        jitter = self.getMaxJitter(spacing)
        jitter2 = jitter * jitter

        if self.compNum == 0 :
            compartment = self.histoVol
        else :
            compartment = self.histoVol.compartments[abs(self.compNum)-1]
            #this is hisotVol for cytoplasme
        compNum = self.compNum
        radius = self.minRadius

        gridPointsCoords = histoVol.grid.masterGridPositions

        # compute rotation matrix rotMat
        if compNum>0:
            # for surface points we compute the rotation which
            # aligns the principalVector with the surface normal
            vx, vy, vz = v1 = self.principalVector
            v2 = compartment.surfacePointsNormals[ptInd]
            try :
                rotMat = numpy.array( rotVectToVect(v1, v2 ), 'f')
            except :
                print('PROBLEM ', self.name)
                rotMat = numpy.identity(4)
        else:
            if self.useRotAxis :
                if sum(self.rotAxis) == 0.0 :
                    rotMat=numpy.identity(4)
                else :
                    rotMat=afvi.vi.rotation_matrix(random()*self.rotRange,self.rotAxis)
            else :
                rotMat=histoVol.randomRot.get()

        if verbose :
            pass#print('%s nbs:%2d'%(self.pdb, len(self.positions[0])), end=' ')

        # jitter position loop
        jitterList = []
        collD1 = []
        collD2 = []

        trans = gridPointsCoords[ptInd] # drop point
#        print ("ptID ",ptInd," coord ",trans)
        gridDropPoint = trans
        jtrans,rotMatj = self.oneJitter(spacing,trans,rotMat)
#        print ("jtrans ",jtrans)
        ok = False
        #here should go the simulation
        #1- we build the ingrediant if not already and place the ingrediant at jtrans, rotMatj
        moving=None
        static=[]
        target=None
        targetPoint = jtrans
#        import c4d
        #c4d.documents.RunAnimation(c4d.documents.GetActiveDocument(), True)

        if runTimeDisplay and self.mesh:
            if hasattr(self,"mesh_3d"):
                #create an instance of mesh3d and place it
                name = self.name + str(ptInd)
                if self.mesh_3d is None :
                    self.moving_geom= afvi.vi.Sphere(name,radius=self.radii[0][0],
                                                    parent=afvi.movingMesh)[0]
                    afvi.vi.setTranslation(self.moving_geom,pos=jtrans)
                else :
                    self.moving_geom= afvi.vi.newInstance(name,self.mesh_3d,
                                                matrice=rotMatj,
                                                location=jtrans, parent = afvi.movingMesh)
        #2- get the neighboring object from ptInd
        if histoVol.ingrLookForNeighbours:
            mingrs,listePartner=self.getListePartners(histoVol,jtrans,rotMat,compartment,afvi)
            for i,elem in enumerate(mingrs):
                ing = elem[2]
                t = elem[0]
                r = elem[1]
                ind = elem[3]
                #print "neighbour",ing.name
                if hasattr(ing,"mesh_3d"):
                    #create an instance of mesh3d and place it
                    name = ing.name + str(ind)
                    if ing.mesh_3d is None :
                        ipoly = afvi.vi.Sphere(name,radius=self.radii[0][0],parent=afvi.staticMesh)[0]
                        afvi.vi.setTranslation(ipoly,pos=t)
                    else :
                        ipoly =afvi.vi.newInstance(name,ing.mesh_3d,matrice=r,
                               location=t, parent = afvi.staticMesh)
#                    static.append(ipoly)
                elif isinstance(ing,GrowIngrediant):
                    name = ing.name + str(ind)
                    ipoly =afvi.vi.newInstance(name,afvi.orgaToMasterGeom[ing],
                                               parent = afvi.staticMesh)
#                    static.append(ipoly)
            
            if listePartner : #self.packingMode=="closePartner":
                if verbose:
                    print("len listePartner = ", len(listePartner))
                if not self.force_random:
                    targetPoint,weight = self.pickPartner(mingrs,listePartner,currentPos=jtrans)
                    if targetPoint is None :
                        targetPoint = jtrans
                else :
                    targetPoint = jtrans
#        print "targetPt",len(targetPoint),targetPoint   
#       should be panda util
#        add the rigid body
        self.histoVol.moving = rbnode = self.histoVol.callFunction(self.histoVol.addRB,(self, jtrans, rotMat,),{"rtype":self.Type},)
        self.histoVol.callFunction(self.histoVol.moveRBnode,(rbnode, jtrans, rotMatj,))
        #run he simulation for simulationTimes
#        afvi.vi.frameAdvanced(duration = simulationTimes,display = runTimeDisplay)#,
        histoVol.callFunction(self.histoVol.runBullet,(self,simulationTimes, runTimeDisplay,))
                              #cb=self.getTransfo)
        rTrans,rRot = self.histoVol.getRotTransRB(rbnode)
        #5- we get the resuling transofrmation matrix and decompose ->rTrans rRot
        #use
        #r=[ (self.histoVol.world.contactTestPair(rbnode, n).getNumContacts() > 0 ) for n in self.histoVol.static]
        self.histoVol.static.append(rbnode)
#        rbnode.setMass(0.0)#make it non dynmaics
#        rTrans = afvi.vi.ToVec(afvi.vi.getTranslation(moving))#getPos
#        rRot = afvi.vi.getMatRotation(moving)                   #getRot or getMat?
        #continue, we assume the object is droped test for contact ? and continueuntil there is  some ?
        ok = True
        jtrans = rTrans[:]
        rotMatj = rRot[:]
        jitterPos = 1
        if ok :
            level = 0
            ## get inside points and update distance
            ##
            # use best sperical approcimation
#            print(">>?",self.name,jtrans)
            insidePoints = {}
            newDistPoints = {}
            t3=time()
            #should be replace by self.getPointInside
            if self.modelType=='Spheres':
                self.centT = centT = self.transformPoints(jtrans, rotMatj, self.positions[level])
                centT = self.centT#self.transformPoints(jtrans, rotMatj, self.positions[-1])
#                for pointsInCube,ptsInSphere in self.distances_temp:
                for radc, posc in zip(self.radii[-1], centT):
#                for i,radc in enumerate(self.radii[-1]):
#                    pointsInCube,ptsInSphere,distA = self.distances_temp[i]
#                 
                    rad = radc + dpad
                    x,y,z = posc
#                    #this have already be done in the checkCollision why doing it again
                    bb = ( [x-rad, y-rad, z-rad], [x+rad, y+rad, z+rad] )
                    pointsInCube = histoVol.callFunction(histoVol.grid.getPointsInCube,
                                                         (bb, posc, rad))
#
                    delta = numpy.take(gridPointsCoords,pointsInCube,0)-posc
                    delta *= delta
                    distA = numpy.sqrt( delta.sum(1) )
                    ptsInSphere = numpy.nonzero(numpy.less_equal(distA, rad))[0]

                    for pti in ptsInSphere:
                        pt = pointsInCube[pti]
                        if pt in insidePoints: continue
                        dist = distA[pti]
                        d = dist-radc
                        if dist < radc:  # point is inside dropped sphere
                            if pt in insidePoints:
                                if d < insidePoints[pt]:
                                    insidePoints[pt] = d
                            else:
                                insidePoints[pt] = d
                        elif d < distance[pt]: # point in region of influence
                            if pt in newDistPoints:
                                if d < newDistPoints[pt]:
                                    newDistPoints[pt] = d
                            else:
                                newDistPoints[pt] = d

            elif self.modelType=='Cylinders':
                cent1T = self.transformPoints(jtrans, rotMatj, self.positions[-1])
                cent2T = self.transformPoints(jtrans, rotMatj, self.positions2[-1])

                for radc, p1, p2 in zip(self.radii[-1], cent1T, cent2T):

                    x1, y1, z1 = p1
                    x2, y2, z2 = p2
                    vx, vy, vz = vect = (x2-x1, y2-y1, z2-z1)
                    lengthsq = vx*vx + vy*vy + vz*vz
                    l = sqrt( lengthsq )
                    cx, cy, cz = posc = x1+vx*.5, y1+vy*.5, z1+vz*.5
                    radt = l + radc + dpad
                    bb = ( [cx-radt, cy-radt, cz-radt], [cx+radt, cy+radt, cz+radt] )
                    pointsInCube = histoVol.callFunction(histoVol.grid.getPointsInCube,
                                                         (bb, posc, radt))

                    pd = numpy.take(gridPointsCoords,pointsInCube,0) - p1
                    dotp = numpy.dot(pd, vect)
                    rad2 = radc*radc
                    d2toP1 = numpy.sum(pd*pd, 1)
                    dsq = d2toP1 - dotp*dotp/lengthsq

                    pd2 = numpy.take(gridPointsCoords,pointsInCube,0) - p2
                    d2toP2 = numpy.sum(pd2*pd2, 1)

                    for pti, pt in enumerate(pointsInCube):
                        if pt in insidePoints: continue

                        if dotp[pti] < 0.0: # outside 1st cap
                            d = sqrt(d2toP1[pti])
                            if d < distance[pt]: # point in region of influence
                                if pt in newDistPoints:
                                    if d < newDistPoints[pt]:
                                        newDistPoints[pt] = d
                                else:
                                    newDistPoints[pt] = d
                        elif dotp[pti] > lengthsq:
                            d = sqrt(d2toP2[pti])
                            if d < distance[pt]: # point in region of influence
                                if pt in newDistPoints:
                                    if d < newDistPoints[pt]:
                                        newDistPoints[pt] = d
                                else:
                                    newDistPoints[pt] = d
                        else:
                            d = sqrt(dsq[pti])-radc
                            if d < 0.:  # point is inside dropped sphere
                                if pt in insidePoints:
                                    if d < insidePoints[pt]:
                                        insidePoints[pt] = d
                                else:
                                    insidePoints[pt] = d
            elif self.modelType=='Cube':
                insidePoints,newDistPoints=self.getDistancesCube(jtrans, rotMatj,gridPointsCoords, distance, histoVol)
            # update free points
            organelle.molecules.append([ jtrans, rotMatj, self, ptInd ])
            histoVol.order[ptInd]=histoVol.lastrank
            histoVol.lastrank+=1    
#            histoVol.close_ingr_bhtree.MoveRBHPoint(histoVol.nb_ingredient,jtrans,0)
            histoVol.nb_ingredient+=1            
#            histoVol.close_ingr_bhtree.InsertRBHPoint((jtrans[0],jtrans[1],jtrans[2]),radius,None,histoVol.nb_ingredient)
            nbFreePoints = self.updateDistances(histoVol,insidePoints, newDistPoints, 
                        freePoints, nbFreePoints, distance, 
                        histoVol.grid.masterGridPositions, verbose)

            # save dropped ingredient
            
        
            self.rRot.append(rotMatj)
            self.tTrans.append(jtrans)
            # add one to molecule counter for this ingredient
            self.counter += 1
            self.completion = float(self.counter)/self.nbMol

            if jitterPos>0:
                histoVol.successfullJitter.append(
                    (self, jitterList, collD1, collD2) )  
            if verbose :
                print('Success nbfp:%d %d/%d dpad %.2f'%(
                nbFreePoints, self.counter, self.nbMol, dpad))
            if self.name=='in  inside':
                histoVol.jitterVectors.append( (trans, jtrans) )

            success = True
            self.rejectionCounter = 0
            
        else: # got rejected
            success = False
            histoVol.failedJitter.append(
                (self, jitterList, collD1, collD2) )

            distance[ptInd] = max(0, distance[ptInd]*0.9)
            self.rejectionCounter += 1
            if verbose :
                print('Failed ingr:%s rejections:%d'%(
                self.name, self.rejectionCounter))
            if self.rejectionCounter >= 30:
                if verbose :
                    print('PREMATURE ENDING of ingredient', self.name)
                self.completion = 1.0
        if drop :
            return success, nbFreePoints
        else :
            return success, nbFreePoints,jtrans, rotMatj
                
class SingleSphereIngr(Ingredient):
    """
    This Ingredient is represented by a single sphere
    and either a single radius, or a list of radii and offset vectors
    for each sphere representing the ingredient
    """
    def __init__(self, molarity=0.0, radius=None, position=None,sphereFile=None,
                 packingPriority=0, name=None, pdb=None,
                 color=None, nbJitter=5, jitterMax=(1,1,1),
                 perturbAxisAmplitude = 0.1,
                 principalVector=(1,0,0), meshFile=None, packingMode='random',
                 placeType="jitter",Type="SingleSphere",
                 meshObject=None,nbMol=0,**kw):

        Ingredient.__init__(self, molarity=molarity, radii=[[radius],], positions=[[position]], #positions2=None,
                 sphereFile=sphereFile, packingPriority=packingPriority, name=name, pdb=pdb, 
                 color=color, nbJitter=nbJitter, jitterMax=jitterMax,
                 perturbAxisAmplitude = perturbAxisAmplitude, principalVector=principalVector,
                 meshFile=meshFile, packingMode=packingMode,placeType=placeType,
                 meshObject=meshObject,nbMol=nbMol,Type=Type,**kw)

        if name == None:
            name = "%5.2f_%f"% (radius,molarity)
        self.name = name
        self.singleSphere = True
        self.minRadius = self.radii[0][0]
        self.encapsulatingRadius = radius
        #make a sphere ?->rapid ?
        if self.mesh is None and autopack.helper is not None  :
            if not autopack.helper.nogui :
#            if not autopack.helper.nogui :
                #build a cylinder and make it length uLength, radius radii[0]
                #this mesh is used bu RAPID for collision
                p=autopack.helper.getObject("autopackHider")
                if p is None:
                    p = autopack.helper.newEmpty("autopackHider")
                    if autopack.helper.host.find("blender") == -1 :
                        autopack.helper.toggleDisplay(p,False)
                self.mesh = autopack.helper.Sphere(self.name+"_basic",
                                radius=self.radii[0][0],color=self.color,
                                parent=p,res=24)[0]
            else :
                self.mesh = autopack.helper.unitSphere(self.name+"_basic",5,
                                radius=self.radii[0][0])[0]
                self.getData()
        #should do that for all ingredient type
        if self.representation is None and not hasattr(self.mesh,"getFaces"):#this is not working with dejavu
            #and should go in the graphics.
            if not autopack.helper.nogui :
                self.representation = autopack.helper.Sphere(self.name+"_rep",
                                radius=self.radii[0][0],color=self.color,
                                parent=self.mesh,res=24)[0]
            else :
                self.representation = autopack.helper.Icosahedron(self.name+"_rep",
                                radius=self.radii[0][0])[0]
            
        
        
        
class SingleCubeIngr(Ingredient):
    """
    This Ingredient is represented by a single cube
    """
    def __init__(self, molarity=0.0, radii=None, 
                 positions=[[[0,0,0],[0,0,0],[0,0,0],]], 
                 positions2=[[[0,0,0]]],
                 sphereFile=None,
                 packingPriority=0, name=None, pdb=None,
                 color=None, nbJitter=5, jitterMax=(1,1,1),
                 perturbAxisAmplitude = 0.1,
                 principalVector=(1,0,0), meshFile=None, packingMode='random',
                 placeType="jitter",Type="SingleCube",
                 meshObject=None,nbMol=0,**kw):

        Ingredient.__init__(self, molarity=molarity, radii=radii, positions=positions, #positions2=None,
                 sphereFile=sphereFile, packingPriority=packingPriority, name=name, pdb=pdb, 
                 color=color, nbJitter=nbJitter, jitterMax=jitterMax,
                 perturbAxisAmplitude = perturbAxisAmplitude, principalVector=principalVector,
                 meshFile=meshFile, packingMode=packingMode,placeType=placeType,
                 meshObject=meshObject,nbMol=nbMol,Type=Type,**kw)

        if name == None:
            name = "%5.2f_%f"% (radii[0][0],molarity)
        self.name = name
        self.singleSphere = False
        self.modelType = "Cube"
        self.collisionLevel = 0
        self.encapsulatingRadius = max(radii[0]) #should the sphere that encapsulated the cube
        self.center = [0.,0.,0.]
        radii = numpy.array(radii)
        self.minRadius =  min(radii[0]/2.0) #should have three radii sizex,sizey,sizez 
        self.encapsulatingRadius = self.maxRadius = math.sqrt(max(radii[0]/2.0)*max(radii[0]/2.0)+min(radii[0]/2.0)*min(radii[0]/2.0))
        self.bb=[-radii[0]/2.,radii[0]/2.]
        self.positions = [[-radii[0]/2.]]
        self.positions2 = [[radii[0]/2.]]
#        if positions2 is not None and positions is not None:
#            self.bb=[positions[0],positions2[0]]
#            x1, y1, z1 = self.bb[0]
#            x2, y2, z2 = self.bb[0]
#            vx, vy, vz = vect = (x2-x1, y2-y1, z2-z1)
#            self.center = x1+vx*.5, y1+vy*.5, z1+vz*.5
#            self.positions = positions
#            self.positions2 = positions2
        #position and position2 are the cornerPoints of the cube
        self.radii = radii


        
class MultiSphereIngr(Ingredient):
    """
    This Ingredient is represented by a collection of spheres specified by radii
    and positions. The principal Vector will be used to align the ingredient
    """
    def __init__(self, molarity=0.0, radii=None, positions=None, sphereFile=None,
                 packingPriority=0, name=None,
                 pdb=None, color=None, nbJitter=5, jitterMax=(1,1,1),
                 perturbAxisAmplitude = 0.1,Type="MultiSphere",
                 principalVector=(1,0,0), meshFile=None, packingMode='random',
                 placeType="jitter",
                 meshObject=None,nbMol=0,**kw):

        Ingredient.__init__(self, molarity=molarity, radii=radii, positions=positions, #positions2=None,
                 sphereFile=sphereFile, packingPriority=packingPriority, name=name, pdb=pdb, 
                 color=color, nbJitter=nbJitter, jitterMax=jitterMax,
                 perturbAxisAmplitude = perturbAxisAmplitude, principalVector=principalVector,
                 meshFile=meshFile, packingMode=packingMode,placeType=placeType,
                 meshObject=meshObject,nbMol=nbMol,Type=Type,**kw)

        if name == None:
            name = "%s_%f"% (str(radii),molarity)
        self.name = name
        self.singleSphere = False
#        print ("done  MultiSphereIngr")
#        raw_input()


class MultiCylindersIngr(Ingredient):
    """
    This Ingredient is represented by a collection of cylinder specified by
    radii, positions and positions2.
    The principal Vector will be used to align the ingredient
    """
    def __init__(self, molarity=0.0, radii=None, positions=None, positions2=None,
                 sphereFile=None, packingPriority=0, name=None,
                 pdb=None, color=None, nbJitter=5, jitterMax=(1,1,1),
                 perturbAxisAmplitude = 0.1,
                 principalVector=(1,0,0), meshFile=None, packingMode='random',
                 placeType="jitter",Type="MultiCylinder",
                 meshObject=None,nbMol=0,**kw):

        Ingredient.__init__(self, molarity=molarity, radii=radii, positions=positions, positions2=positions2,
                 sphereFile=sphereFile, packingPriority=packingPriority, name=name, pdb=pdb, 
                 color=color, nbJitter=nbJitter, jitterMax=jitterMax,
                 perturbAxisAmplitude = perturbAxisAmplitude, principalVector=principalVector,
                 meshFile=meshFile, packingMode=packingMode,placeType=placeType,
                 meshObject=meshObject,nbMol=nbMol,Type=Type,**kw)

        if name == None:
            name = "%s_%f"% (str(radii),molarity)
        self.name = name
        self.singleSphere = False
        self.modelType = "Cylinders"
        self.collisionLevel = 0
        self.minRadius = radii[0][0]
#        self.encapsulatingRadius = radii[0][0]  #Graham worry: 9/8/11 This is incorrect... shoudl be max(radii[0]) or radii[0][1] 
#        self.encapsulatingRadius = radii[0][0]#nope should be  half length ?
        self.length= 1.0
        
#        print('encapsulating Radius is probably wrong- fix the array')
        self.useLength = False
        if "useLength" in kw:
            self.useLength = kw["useLength"]
        if positions2 is not None and positions is not None:
            #shoulde the overall length of the object from bottom to top
            bb = self.getBigBB()
            d = numpy.array(bb[1])-numpy.array(bb[0])
            s = numpy.sum(d*d)
            self.length  = math.sqrt(s ) #diagonal
#        if self.mesh is None and autopack.helper is not None :
#            #build a cylinder and make it length uLength, radius radii[0]
#            self.mesh = autopack.helper.Cylinder(self.name+"_basic",radius=self.radii[0][0],
#                                       length=self.uLength,parent="autopackHider")[0]
        if self.mesh is None and autopack.helper is not None  :
            p=None
            if not autopack.helper.nogui :
                #build a cylinder and make it length uLength, radius radii[0]
                #this mesh is used bu RAPID for collision
                p=autopack.helper.getObject("autopackHider")
                if p is None:
                    p = autopack.helper.newEmpty("autopackHider")
                    if autopack.helper.host.find("blender") == -1 :
                        autopack.helper.toggleDisplay(p,False)
#                self.mesh = autopack.helper.Cylinder(self.name+"_basic",
#                                radius=self.radii[0][0]*1.24, length=self.uLength,
#                                res= 5, parent="autopackHider",axis="+X")[0]
            length = 1
            if positions2 is not None and positions is not None:
                d = numpy.array(positions2[0][0])-numpy.array(positions[0][0])
                s = numpy.sum(d*d)
                length  = math.sqrt(s ) #diagonal
            self.mesh = autopack.helper.Cylinder(self.name+"_basic",
                                radius=self.radii[0][0]*1.24, length=length,
                                res= 5, parent="autopackHider",axis=self.principalVector)[0]                
#            self.mesh = autopack.helper.oneCylinder(self.name+"_basic",
#                                self.positions[0][0],self.positions2[0][0],
#                                radius=self.radii[0][0]*1.24,
#                                parent = p,color=self.color)
#            self.getData()

        self.KWDS["useLength"]={}

    def getBigBB(self):
        #one level for cylinder
        bbs=[]
        for radc, p1, p2 in zip(self.radii[0], self.positions[0], self.positions2[0]):            
            bb = self.correctBB(p1,p2,radc)
            bbs.append(bb)
        #get min and max from all bbs
        maxBB = [0,0,0]
        minBB = [9999,9999,9999]
        for bb in bbs:
            for i in range(3):
                if bb[0][i] < minBB[i]:
                    minBB[i] =bb[0][i]
                if bb[1][i] > maxBB[i]:
                    maxBB[i] = bb[1][i]
                if bb[1][i] < minBB[i]:
                    minBB[i] = bb[1][i]
                if bb[0][i] > maxBB[i]:
                    maxBB[i] = bb[0][i]
        bb = [minBB,maxBB]
        return bb        

#do we need this class?
class GrowIngrediant(MultiCylindersIngr):
    def __init__(self, molarity=0.0, radii=None, positions=None, positions2=None,
                 sphereFile=None, packingPriority=0, name=None,
                 pdb=None, color=None, nbJitter=5, jitterMax=(1,1,1),
                 perturbAxisAmplitude = 0.1,length = 10.,closed = False,
                 modelType="Cylinders",biased=1.0,
                 principalVector=(1,0,0), meshFile=None, packingMode='random',
                 placeType="jitter",marge=20.0,meshObject=None,orientation = (1,0,0),
                 nbMol=0,Type="Grow",walkingMode="sphere",**kw):

        MultiCylindersIngr.__init__(self, molarity=molarity, radii=radii, positions=positions, positions2=positions2,
                 sphereFile=sphereFile, packingPriority=packingPriority, name=name, pdb=pdb, 
                 color=color, nbJitter=nbJitter, jitterMax=jitterMax,
                 perturbAxisAmplitude = perturbAxisAmplitude, principalVector=principalVector,
                 meshFile=meshFile, packingMode=packingMode,placeType=placeType,
                 meshObject=meshObject,nbMol=nbMol,Type=Type,**kw)
        if name == None:
            name = "%s_%f"% (str(radii),molarity)
        self.name = name
        self.singleSphere = False
        self.modelType = modelType
        self.collisionLevel = 0
        self.minRadius = radii[0][0]
#        self.encapsulatingRadius = radii[0][0] #Graham worry: 9/8/11 This is incorrect... shoudl be max(radii[0]) or radii[0][1] 
        self.encapsulatingRadius = radii[0][0]
#        print('encapsulating Radius is probably wrong- fix the array')
        self.marge = marge
        self.length = length
        self.closed = closed
        self.nbCurve=0
        self.results=[]
        self.start_positions=[]
        self.listePtLinear=[]
        self.listePtCurve=[] #snake
        self.Ptis=[]    #snakePts
        self.currentLength=0.#snakelength
        self.direction=None#direction of growing
        #can be place either using grid point/jittering or dynamics
#        self.uLength = 0. #length of the cylinder or averall radius for sphere, this is the lenght of one unit
        self.uLength = 10
        if "uLength" in kw :
            self.uLength = kw["uLength"]
        else :
            v,u = self.vi.measure_distance(self.positions,self.positions2,vec=True)
            self.uLength = abs(u)
        self.encapsulatingRadius = self.uLength/2.0
        self.unitNumberF = 0 #number of unit pose so far forward
        self.unitNumberR = 0 #number of unit pose so far reverse
        self.orientation = orientation
        self.seedOnPlus = True   #The filament should continue to grow on its (+) end
        self.seedOnMinus = False #The filamen should continue to grow on its (-) end.
#        if self.compNum > 0 :
#            self.seedOnMinus = False
        self.vector=[0.,0.,0.]
        self.biased = biased
        self.absoluteLengthMax =99999999.9999#(default value is infinite or some safety number like 1billion)
        self.probableLengthEquation = None 
        #(this can be a number or an equation, e.g., every actin should grow to 
        #10 units long, or this actin fiber seed instance should grow to (random*10)^2
        #actually its a function of self.uLength like self.uLength * 10. or *(random*10)^2
        self.ingGrowthMatrix = numpy.identity(4)
        #(ultimately, we'll build a database for these (David Goodsell has a few examples), 
        #but users should be able to put in their own, so for a test for now, lets say we'll 
        #grow one unit for a singleSphereIng r=60 along its X as [55,0,0;1,0,0;0,1,0;0,0,1]
        self.ingGrowthWobbleFormula = None
        #(this could be a rotation matrix to make a helix, or just a formula, 
        #like the precession algorithm Michel and I already put in 
        #for surface ingredients.
        self.constraintMarge = False
        self.cutoff_boundary = 1.0
        self.cutoff_surface = 5.0
        self.useHalton = True        
        self.use_rbsphere = False #use sphere instead of cylinder or collision with bullet
        if "useHalton" in kw :
            self.useHalton = kw["useHalton"]
        if "constraintMarge" in kw :
            self.constraintMarge = kw["constraintMarge"]
        if "cutoff_boundary" in kw :
            self.cutoff_boundary = kw["cutoff_boundary"]
        if "cutoff_surface" in kw :
            self.cutoff_surface = kw["cutoff_surface"]
        if "use_rbsphere" in kw :
            self.use_rbsphere = kw["use_rbsphere"]
        if "encapsulatingRadius" in kw :
            self.use_rbsphere = kw["encapsulatingRadius"]
        #mesh object representing one uLength? or a certain length
        self.unitParent = None
        self.unitParentLength = 0.
        self.walkingMode = walkingMode #["sphere","lattice"]
        self.walkingType = "stepbystep" #or atonce
        self.compMask = []
        #default should be compId
        if "compMask" in kw :
            if type(kw["compMask"]) is str :
                self.compMask = eval(kw["compMask"])    
            else :
                self.compMask = kw["compMask"]
        #create a simple geom if none pass?        
        #self.compMask=[]
        if self.mesh is None and autopack.helper is not None  :
            p=None
            if not autopack.helper.nogui :
                #build a cylinder and make it length uLength, radius radii[0]
                #this mesh is used bu RAPID for collision
                p=autopack.helper.getObject("autopackHider")
                if p is None:
                    p = autopack.helper.newEmpty("autopackHider")
                    if autopack.helper.host.find("blender") == -1 :
                        autopack.helper.toggleDisplay(p,False)
#                self.mesh = autopack.helper.Cylinder(self.name+"_basic",
#                                radius=self.radii[0][0]*1.24, length=self.uLength,
#                                res= 5, parent="autopackHider",axis="+X")[0]
#            else :
            #is axis working ?
            self.mesh = autopack.helper.Cylinder(self.name+"_basic",
                                radius=self.radii[0][0]*1.24, length=self.uLength,
                                res= 32, parent="autopackHider",axis=self.orientation)[0]                
            if autopack.helper.nogui : self.getData()
        self.sphere_points_nb = 50000    
        self.sphere_points = numpy.array(SphereHalton(self.sphere_points_nb,5))
        self.sphere_points_mask = numpy.ones(self.sphere_points_nb,'i')
        self.sphere_points_masked = None

        #need to define the binder/modifier. This is different from partner
        #Every nth place alternative repesentation
        #name proba is this for ingredient in general ?
        self.alternates_names = []
        if "alternates_names" in kw :
            self.alternates_names = kw["alternates_names"]
        self.alternates_proba = []
        if "alternates_proba" in kw :
            self.alternates_proba = kw["alternates_proba"]
        self.alternates_weight = []
        if "alternates_weight" in kw :
            self.alternates_weight = kw["alternates_weight"]
        self.prev_alt=None
        self.prev_was_alternate=False
        self.prev_alt_pt=None        
        self.mini_interval = 2
        self.alternate_interval = 0
        #keep record of point Id that are bound to alternate and change the 
        #representation according.
        self.safetycutoff = 10
        
        self.KWDS["length"]={}
        self.KWDS["closed"]={}
        self.KWDS["uLength"]={}        
        self.KWDS["biased"]={}
        self.KWDS["marge"]={}
        self.KWDS["orientation"]={}
        self.KWDS["walkingMode"]={}
        self.KWDS["constraintMarge"]={}
        self.KWDS["useHalton"]={}
        self.KWDS["compMask"]={}
        self.KWDS["use_rbsphere"]={}
        
        self.OPTIONS["length"]={"name":"length","value":self.length,"default":self.length,"type":"float","min":0,"max":10000,"description":"snake total length"}
        self.OPTIONS["uLength"]={"name":"uLength","value":self.uLength,"default":self.uLength,"type":"float","min":0,"max":10000,"description":"snake unit length"}
        self.OPTIONS["closed"]={"name":"closed","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"closed snake"}                          
        self.OPTIONS["biased"]={"name":"biased","value":0.0,"default":0.0,"type":"float","min":0,"max":10,"description":"snake biased"}
        self.OPTIONS["marge"]={"name":"marge","value":10.0,"default":10.0,"type":"float","min":0,"max":10000,"description":"snake angle marge"}
        self.OPTIONS["constraintMarge"]={"name":"constraintMarge","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"lock the marge"}                          
        self.OPTIONS["orientation"]={"name":"orientation","value":[0.,1.,0.],"default":[0.,1.,0.],"min":0,"max":1,"type":"vector","description":"snake orientation"}
        self.OPTIONS["walkingMode"]={"name":"walkingMode","value":"random","values":['sphere', 'lattice'],"min":0.,"max":0.,"default":'sphere',"type":"liste","description":"snake mode"}
        self.OPTIONS["useHalton"]={"name":"useHalton","value":True,"default":True,"type":"bool","min":0.,"max":0.,"description":"use spherica halton distribution"}                          
        self.OPTIONS["compMask"]={"name":"compMask","value":"0","values":"0","min":0.,"max":0.,"default":'0',"type":"string","description":"allowed compartments"}
        self.OPTIONS["use_rbsphere"]={"name":"use_rbsphere","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"use sphere instead of cylinder wih bullet"}                          
        
    def resetSphereDistribution(self):
        #given a radius, create the sphere distribution
        self.sphere_points = SphereHalton(self.sphere_points_nb,5)
        self.sphere_points_mask = numpy.ones(self.sphere_points_nb,'i')
        
    def getNextPoint(self):
        #pick a random point from the sphere point distribution
        pointsmask=numpy.nonzero(self.sphere_points_mask)[0]
        if not len(pointsmask):
            print ("no sphere point available from mask")
            return None
        ptIndr = int(uniform(0.0,1.0)*len(pointsmask))
        sp_pt_indice = pointsmask[ptIndr]
        np = numpy.array(self.sphere_points[sp_pt_indice])*numpy.array(self.jitterMax) 
        return numpy.array(self.vi.unit_vector(np)) *self.uLength   #biased by jitterMax ?                      

    def mask_sphere_points_boundary(self,pt,boundingBox = None):
        if boundingBox is None :
            boundingBox = self.histoVol.fillBB
        pts =  (numpy.array(self.sphere_points)*self.uLength)+pt
        points_mask = numpy.nonzero(self.sphere_points_mask)[0]
        if len(points_mask):
            mask = [not self.testPoint(pt) for pt in pts[points_mask]]
            if len(mask) :
                self.sphere_points_mask[points_mask]=numpy.logical_and(mask,
                                self.sphere_points_mask[points_mask])

    def mask_sphere_points_ingredients(self,pt,listeclosest):
        listeclosest=[elem for elem in listeclosest if not isinstance(elem[3],autopack.Compartment.Compartment)]
        points_mask = numpy.nonzero(self.sphere_points_mask)[0]
        if len(listeclosest) and len(points_mask):
            points = numpy.array(listeclosest)[:,1]
            ingrs = numpy.array(listeclosest)[:,3]
            radius = [float(ingr.encapsulatingRadius) for ingr in ingrs]
            #this distance is between unit vector and 3d points...
            #translate and scale the spheres points
            sp=numpy.array(self.sphere_points, dtype=numpy.float64, copy=False)
            dp = sp[points_mask]*self.uLength+pt
            pts = numpy.array(points.tolist(), dtype=numpy.float64, copy=False)
            #distance between sphere point and ingredient positions
            distances=spatial.distance.cdist(pts,dp)
            #empty mask for the point
            mask = numpy.nonzero(numpy.ones(len(dp)))[0]
            #mas cumulative ?for 
            for i in range(len(distances)):
                 #if distance is >= to ingredient encapsulatingRadius we keep the point
                 m=numpy.greater_equal(distances[i],radius[i])   
                 mask = numpy.logical_and(mask,m)
            #ponts to keep
            self.sphere_points_mask[points_mask]=numpy.logical_and(mask,
                                self.sphere_points_mask[points_mask])
            #indices = numpy.nonzero(mask)[0]#indice f
            #i = self.sphere_points_mask[indices]
            #self.sphere_points_mask = i  
                                
    def mask_sphere_points_dihedral(self,v1,v2,marge_out,marge_diedral):
        points_mask = numpy.nonzero(self.sphere_points_mask)[0]
        a=angle_between_vectors(self.vi.unit_vector(v2),self.sphere_points[points_mask], axis=1)
        b=angle_between_vectors(self.vi.unit_vector(v1),self.sphere_points[points_mask], axis=1)
        if type(marge_out) is float :
            marge_out = [0.0,marge_out]
        if type(marge_diedral) is float :
            marge_diedral = [0.0,marge_diedral]
        mask1 = numpy.logical_and(a<math.radians(marge_out[1]), a > math.radians(marge_out[0]))
        mask2 = numpy.logical_and(b<math.radians(marge_diedral[1]), b > math.radians(marge_diedral[0]))
        mask3 = numpy.logical_and(mask1,mask2)
        self.sphere_points_mask[points_mask]=numpy.logical_and(mask3,
                                self.sphere_points_mask[points_mask])
                                
    def mask_sphere_points_angle(self,v,marge_in):
        #mask first with angle
        #adjust the points to current transfomation? or normalize current vector ?
        points_mask = numpy.nonzero(self.sphere_points_mask)[0]
        a=angle_between_vectors(self.vi.unit_vector(v),self.sphere_points[points_mask], axis=1)
        if type(marge_in) is float :
            mask = numpy.less (a,math.radians(marge_in))
        else :
            mask = numpy.logical_and(a<math.radians(marge_in[1]), a > math.radians(marge_in[0]))                   
#        print len(mask),marge_in,v
#        self.sphere_points_mask = mask#numpy.nonzero(mask)[0]#points to keep
        self.sphere_points_mask[points_mask]=numpy.logical_and(mask,
                                self.sphere_points_mask[points_mask])

    def mask_sphere_points(self,v,pt,marge,listeclosest,cutoff,pv=None,marge_diedral=None):
        #self.sphere_points_mask=numpy.ones(10000,'i') 
        print ("mask_sphere_points ",marge_diedral,marge,len(listeclosest))
        if marge_diedral is not None :
            self.mask_sphere_points_dihedral(pv,v,marge,marge_diedral)        
        else :            
            self.mask_sphere_points_angle(v,marge)
        #storethe mask point
        sphere_points_mask_copy = numpy.copy(self.sphere_points_mask)
        self.mask_sphere_points_ingredients(pt,listeclosest)
        if not len( numpy.nonzero(self.sphere_points_mask)[0] ):
            self.sphere_points_mask=numpy.copy(sphere_points_mask_copy)
        #self.mask_sphere_points_boundary(pt)              
                
    def reset(self):
        self.nbCurve=0
        self.results=[]
        self.listePtLinear=[]
        self.listePtCurve=[] #snake
        self.Ptis=[]    #snakePts
        self.currentLength=0.#snakelength        
        #update he cylinder ?
        
    def getNextPtIndCyl(self,jtrans, rotMatj,freePoints,histoVol):
#        print jtrans, rotMatj
        cent2T = self.transformPoints(jtrans, rotMatj, self.positions[-1])
        jx, jy, jz = self.jitterMax
        jitter = self.getMaxJitter(histoVol.smallestProteinSize)
        jitter2 = jitter * jitter
        if len(cent2T) == 1 :
            cent2T=cent2T[0]
        tx, ty, tz = cent2T
        dx = jx*jitter*gauss(0., 0.3)  #This is an incorrect jitter use the uniform random with sphereical rejection
        dy = jy*jitter*gauss(0., 0.3)
        dz = jz*jitter*gauss(0., 0.3)
#        print "d",dx,dy,dz
        nexPt = (tx+dx, ty+dy, tz+dz)
        #where is this point in the grid
        #ptInd = histoVol.grid.getPointFrom3D(nexPt)
        t,r = self.oneJitter(histoVol.smallestProteinSize,cent2T,rotMatj)
        dist,ptInd = histoVol.grid.getClosestGridPoint(t)
        dv = numpy.array(nexPt) - numpy.array(cent2T)
        d = numpy.sum(dv*dv)
        return ptInd,dv,sqrt(d)

    def getJtransRot_r(self,pt1,pt2,length = None):
        if length is None :
            length = self.uLength
        #print "input is ",pt1,pt2,self.orientation
        v = numpy.array(pt2) - numpy.array(pt1)
        pmx = rotVectToVect(numpy.array(self.orientation) * length,v, i=None)
        return numpy.array(pmx),numpy.array(pt1)+numpy.array(v)/2.#.transpose()jtrans      

    def getJtransRot(self,pt1,pt2):
#        print "input is ",pt1,pt2
#        v = numpy.array(pt1) - numpy.array(pt2)
#        pmx = rotVectToVect(v,numpy.array(self.orientation) * self.uLength, i=None)
#        return  numpy.array(pmx),numpy.array(pt1)+numpy.array(v)/2.#.transpose()jtrans      
        v,d = self.vi.measure_distance(pt1,pt2,vec=True)
        length, mat = autopack.helper.getTubePropertiesMatrix(pt1,pt2)
        return  numpy.array(mat),numpy.array(pt1)+numpy.array(v)/2.#.transpose()jtrans    
        
        
        #Start jtrans section that is new since Sept 8, 2011 version
        n = numpy.array(pt1) - numpy.array(pt2)
        #normalize the vector n
        nn=self.vi.unit_vector(n) #why normalize ?
#        print ("getJtranseRot for")
#        print (nn)
#        print (numpy.array(self.orientation))
        
        #get axis of rotation between the plane normal and Z
        v1 = nn
        v2 = numpy.array([0.,0.,1.0])#self.orientation) #[0.,0.,1.0]
        cr = numpy.cross(v1,v2)
        axis = self.vi.unit_vector(cr)
        
        #get the angle between the plane normal and Z
        angle = self.vi.angle_between_vectors(v2,v1)
        #get the rotation matrix between plane normal and Z
        print (axis,angle)

        mx = self.vi.rotation_matrix(-angle, axis)#.transpose()-angle ?
        jtrans = n
        #End jtrans section that is new since Sept 8, 2011 version       
        matrix = mx.transpose()# self.vi.ToMat(mx).transpose()#Why ToMat here ?
        rotMatj = matrix.reshape((4,4))
        return rotMatj.transpose(),numpy.array(pt1)+numpy.array(v)/2.#.transpose()jtrans

    def walkLatticeSurface(self,pts,distance,histoVol, size,mask,marge=999.0,
                    checkcollision=True,saw=True):

        o = self.histoVol.compartments[abs(self.compNum)-1]
        sfpts = o.surfacePointsCoords
        
        found = False
        attempted = 0
        safetycutoff = 200
        if self.runTimeDisplay:
            name = "walking"+self.name
            sp = self.vi.getObject(name)
            if sp is None :
                sp=self.vi.Sphere(name,radius=10.0)[0]

            self.vi.update()        
        while not found :
            #print attempted
            if attempted > safetycutoff:
                return None,False
            newPtId = int(random()*len(sfpts))
            v = sfpts[newPtId]#histoVol.grid.masterGridPositions[newPtId]
#            print newPtId,v,len(pointsInCube)
            if self.runTimeDisplay :
                self.vi.setTranslation(sp,self.vi.FromVec(v))
                self.vi.update()            
            if saw : #check if already taken, but didnt prevent cross
                if pointsInCube[newPtId] in self.Ptis:
                    attempted +=1
                    continue
            angle=self.vi.angle_between_vectors(numpy.array(posc),numpy.array(v))
            v,d = self.vi.measure_distance(numpy.array(posc),numpy.array(v),vec=True)
#            print angle,abs(math.degrees(angle)),marge,d
            if abs(math.degrees(angle)) <= marge :
                closeS = self.checkPointSurface(v,cutoff=self.cutoff_surface)
                inComp = self.checkPointComp(v)
                if closeS or not inComp :#or d > self.uLength:
                    #print "closeS or not good comp or too long"
                    attempted +=1
                    continue
                if checkcollision:
                    #print "checkColl"
                    m = numpy.identity(4)
                    collision = self.checkSphCollisions([v,],[float(self.uLength)*1.,], 
                                            [0.,0.,0.], m, 0,
                                            histoVol.grid.masterGridPositions,
                                            distance, 
                                            histoVol)
                    if not collision :
                        found = True
                        self.Ptis.append(pointsInCube[newPtId])
                        return v,True
                    else : #increment the range
                        #print "collision"
                        if not self.constraintMarge :
                            if marge >=180 :
                                return None,False
                            marge+=1
                        attempted +=1
                        continue
                found = True
                self.Ptis.append(pointsInCube[newPtId])
                return v,True
            else :
                attempted+=1
                continue
        

    def walkLattice(self,pts,distance,histoVol, size,marge=999.0,
                    checkcollision=True,saw=True):
        #take the next random point in the windows size +1 +2
        #extended = histoVol.getPointsInCube()
        cx,cy,cz = posc = pts#histoVol.grid.masterGridPositions[ptId]
        step = histoVol.grid.gridSpacing*size
        bb = ( [cx-step, cy-step, cz-step], [cx+step, cy+step, cz+step] )
#        print pts,step,bb
        if self.runTimeDisplay > 1:
            box = self.vi.getObject("collBox")
            if box is None:
                box = self.vi.Box('collBox', cornerPoints=bb,visible=1)
            else :
#                    self.vi.toggleDisplay(box,True)
                self.vi.updateBox(box,cornerPoints=bb)
                self.vi.update()
        pointsInCube = histoVol.grid.getPointsInCube(bb, posc, step,addSP=False)
        pointsInCubeCoords=numpy.take(histoVol.grid.masterGridPositions,pointsInCube,0)
#        print pointsInCube
        #take a random point from it OR use gradient info OR constrain by angle
        found = False
        attempted = 0
        safetycutoff = 200
        if self.runTimeDisplay:
            name = "walking"+self.name
            sp = self.vi.getObject(name)
            if sp is None :
                sp=self.vi.Sphere(name,radius=10.0)[0]
            namep = "latticePoints"
            pts= self.vi.getObject(namep)
            if pts is None :
                pts = self.vi.Points(namep)
            pts.Set(vertices = pointsInCubeCoords)
            #sp.SetAbsPos(self.vi.FromVec(startingPoint))
            self.vi.update()        
        while not found :
            #print attempted
            if attempted > safetycutoff:
                return None,False
            newPtId = int(random()*len(pointsInCube))
            v = pointsInCubeCoords[newPtId]#histoVol.grid.masterGridPositions[newPtId]
#            print newPtId,v,len(pointsInCube)
            if self.runTimeDisplay :
                self.vi.setTranslation(sp,self.vi.FromVec(v))
                self.vi.update()            
            if saw : #check if already taken, but didnt prevent cross
                if pointsInCube[newPtId] in self.Ptis:
                    attempted +=1
                    continue
            angle=self.vi.angle_between_vectors(numpy.array(posc),numpy.array(v))
            v,d = self.vi.measure_distance(numpy.array(posc),numpy.array(v),vec=True)
#            print angle,abs(math.degrees(angle)),marge,d
            if abs(math.degrees(angle)) <= marge :
                closeS = self.checkPointSurface(v,cutoff=self.cutoff_surface)
                inComp = self.checkPointComp(v)
                if closeS or not inComp :#or d > self.uLength:
                    #print "closeS or not good comp or too long"
                    attempted +=1
                    continue
                if checkcollision:
                    #print "checkColl"
                    m = numpy.identity(4)
                    collision = self.checkSphCollisions([v,],[float(self.uLength)*1.,], 
                                            [0.,0.,0.], m, 0,
                                            histoVol.grid.masterGridPositions,
                                            distance, 
                                            histoVol)
                    if not collision :
                        found = True
                        self.Ptis.append(pointsInCube[newPtId])
                        return v,True
                    else : #increment the range
                        #print "collision"
                        if not self.constraintMarge :
                            if marge >=180 :
                                return None,False
                            marge+=1
                        attempted +=1
                        continue
                found = True
                self.Ptis.append(pointsInCube[newPtId])
                return v,True
            else :
                attempted+=1
                continue
        
    def walkSphere(self,pt1,pt2,distance,histoVol,marge = 90.0,checkcollision=True):
        """ use a random point on a sphere of radius uLength, and useCylinder collision on the grid """ 
        v,d = self.vi.measure_distance(pt1,pt2,vec=True)
        found = False
        attempted = 0
        pt = [0.,0.,0.]
        angle=0.
        safetycutoff = 10000
        if self.constraintMarge:
            safetycutoff = 200
        if self.runTimeDisplay:
            name = "walking"+self.name
            sp = self.vi.getObject(name)
            if sp is None :
                sp=self.vi.Sphere(name,radius=2.0)[0]
            #sp.SetAbsPos(self.vi.FromVec(startingPoint))
            self.vi.update()
        while not found :
            #main loop thattryto found the next point (similar to jitter)
            if attempted >= safetycutoff:
                return None,False#numpy.array(pt2).flatten()+numpy.array(pt),False
            print ("length is ",self.uLength)
            pt = self.vi.randpoint_onsphere(self.uLength)#*numpy.array(self.jitterMax)
            #the new position is the previous point (pt2) plus the random point
            newPt = numpy.array(pt2).flatten()+numpy.array(pt)
            if self.runTimeDisplay  >= 2:
                self.vi.setTranslation(sp,newPt)
                self.vi.update()
            #compute the angle between the previous direction (pt1->pt2) and the new random one (pt)
            angle=self.vi.angle_between_vectors(numpy.array(v),numpy.array(pt))
#            print angle,math.degrees(angle),marge
#            print self.histoVol.grid.boundingBox #left-bottom to right-top
            #first test angle less than the constraint angle
            if abs(math.degrees(angle)) <= marge :
#                print "ok"
                #check if in bounding box
                inside = histoVol.grid.checkPointInside(newPt,dist=self.cutoff_boundary,jitter=self.jitterMax)
                closeS = self.checkPointSurface(newPt,cutoff=self.cutoff_surface)
                inComp = self.checkPointComp(newPt)
#                print "inside,closeS ",inside,closeS
                if not inside or closeS or not inComp:
#                    print "oustide"
                    if not self.constraintMarge :
                        if marge >=175 :
                            return None,False
                        marge+=1
                    else :
                        attempted +=1
                    continue
                #optionally check for collision
                if checkcollision:
#                    print self.modelType
                    if self.modelType=='Cylinders':
                        #outise is consider as collision...?
#                        rotMatj,jtrans=self.getJtransRot(numpy.array(pt2).flatten(),newPt)
                        m=[[ 1.,  0.,  0.,  0.],
                           [ 0.,  1.,  0.,  0.],
                           [ 0.,  0.,  1.,  0.],
                           [ 0.,  0.,  0.,  1.]]
##                        print rotMatj,jtrans
#                        print "before collide"
#                        collision = self.checkSphCollisions([newPt,],[float(self.uLength)*1.,], 
#                                            [0.,0.,0.], m, 0,
#                                            histoVol.grid.masterGridPositions,
#                                            distance, 
#                                            histoVol)
                        #use panda ?    
                        collision = self.checkCylCollisions(
                            [numpy.array(pt2).flatten()], [newPt],
                            self.radii[-1], [0.,0.,0.], m, 
                            histoVol.grid.masterGridPositions,
                            distance, 
                            histoVol)
#                        print "collision",collision
                        if not collision :
                            found = True
                            return numpy.array(pt2).flatten()+numpy.array(pt),True
                        else : #increment the range
                            if not self.constraintMarge :
                                if marge >=180 :
                                    return None,False
                                marge+=1
                            else :
                                attempted +=1
                            continue
                else :
                    found = True
                    return numpy.array(pt2).flatten()+numpy.array(pt),True
#                attempted += 1
            else :
                attempted += 1
                continue
            attempted += 1
        return numpy.array(pt2).flatten()+numpy.array(pt),True

    def getInterpolatedSphere(self,pt1,pt2):
        v,d = self.vi.measure_distance(pt1,pt2,vec=True)
#        d=self.uLength
        nbSphere = int (d / self.minRadius)
        sps = numpy.arange(0,d,self.minRadius*2)
        r=[] 
        p=[]
        pt1= numpy.array(pt1)
        pt2= numpy.array(pt2)
        vn= numpy.array(v)/numpy.linalg.norm(numpy.array(v))#normalized
#        print pt1,pt2,d,v,vn,sps,self.minRadius
        p.append(pt1)
        r.append(self.minRadius)
        for i,sp in enumerate(sps[1:]) :
            r.append(self.minRadius)
            p.append(pt1+(vn*sp))
        p.append(pt2)
        r.append(self.minRadius)
#        print ("get ",len(p),p)
#        print ("get ",len(r),r)
        return [r,p]  
        
    def addRBsegment(self,pt1,pt2,nodeid=""):
        #build a rigid body of multisphere along pt1topt2
        r,p = self.getInterpolatedSphere(pt1,pt2)
        print ("pos len",len( p)," ",len(r))
        inodenp=self.histoVol.multiSphereRB(self.name+nodeid,p,r)
        print ("node build ",inodenp)
        inodenp.setCollideMask(self.histoVol.BitMask32.allOn())
        inodenp.node().setAngularDamping(1.0)
        inodenp.node().setLinearDamping(1.0)
        print ("attach node to world")
        #inodenp.setMat(pmat)
        self.histoVol.world.attachRigidBody(inodenp.node())
        print ("node attached to world")
        inodenp = inodenp.node()
        print ("add msphere ",inodenp)
        self.histoVol.rb_panda.append(inodenp)
        return inodenp

    def walkSpherePanda(self,pt1,pt2,distance,histoVol,marge = 90.0,
                        checkcollision=True,usePP=False):
        """ use a random point on a sphere of radius uLength, and useCylinder collision on the grid """ 
#        print ("walkSpherePanda ",pt1,pt2  )      
        v,d = self.vi.measure_distance(pt1,pt2,vec=True)
        found = False
        attempted = 0
        pt = [0.,0.,0.]
        angle=0.
        safetycutoff  = self.rejectionThreshold
        mask=None
        if self.constraintMarge:
            safetycutoff  = self.rejectionThreshold
        if self.runTimeDisplay:
            name = "walking"+self.name
            sp = self.vi.getObject(name)
            if sp is None :
                sp=self.vi.Sphere(name,radius=2.0)[0]
            #sp.SetAbsPos(self.vi.FromVec(startingPoint))
            self.vi.update()
        liste_nodes=[]
        cutoff=self.histoVol.largestProteinSize+self.uLength
        closesbody_indice = self.getClosestIngredient(pt2,self.histoVol,cutoff=cutoff)
        liste_nodes = self.get_rbNodes(closesbody_indice,pt2,prevpoint=pt1,getInfo=True)  
        alternate,ia = self.pick_alternate()
        print ("pick alternate",alternate,ia,self.prev_alt)
        if self.prev_was_alternate :
            alternate = None
        p_alternate = None#self.partners[alternate]#self.histoVol.getIngrFromNameInRecipe(alternate,self.recipe )
        #if p_alternate.getProperties("bend"):
        nextPoint= None#p_alternate.getProperties("pt2")
        marge_in= None#p_alternate.getProperties("marge_in")            
        dihedral= None#p_alternate.getProperties("diehdral")#for next point
        length= None#p_alternate.getProperties("length")#for next point
        #prepare halton if needed
        newPt = None
        newPts=[]
        if self.prev_alt is not None :#and self.prev_alt_pt is not None:
            p_alternate = self.partners[self.prev_alt]
            dihedral= p_alternate.getProperties("diehdral")
            nextPoint= p_alternate.getProperties("pt2")
            marge_in= p_alternate.getProperties("marge_out")  #marge out ? 
            if dihedral is not None :
                self.mask_sphere_points(v,pt1,marge_in,liste_nodes,0,
                                        pv=self.prev_vec,marge_diedral=dihedral)
            alternate = None
        t1 = time()
        if self.prev_alt is not None :
            test = True
        elif alternate and not self.prev_was_alternate :
            #next point shouldnt be an alternate
            p_alternate = self.partners[alternate]#self.histoVol.getIngrFromNameInRecipe(alternate,self.recipe )
            #if p_alternate.getProperties("bend"):
            nextPoint=p_alternate.getProperties("pt2")
            marge_in=p_alternate.getProperties("marge_in")            
            dihedral=p_alternate.getProperties("diehdral")#for next point
            length=p_alternate.getProperties("length")#for next point
            if marge_in is not None :
                self.mask_sphere_points(v,pt2,marge_in,liste_nodes,0,pv=None,marge_diedral=None)
            else :
                self.mask_sphere_points(v,pt2,marge+2.0,liste_nodes,0) 
            self.prev_was_alternate = True
        elif self.useHalton: 
            self.mask_sphere_points(v,pt2,marge+2.0,liste_nodes,0)
        if histoVol.runTimeDisplay:
            points_mask=numpy.nonzero(self.sphere_points_mask)[0]
            verts=self.sphere_points[points_mask]*self.uLength+pt2            
            name = "Hcloud"+self.name
            pc = self.vi.getObject(name)
            if pc is None :
                pc=self.vi.PointCloudObject(name,vertices=verts)[0]
            else :
                self.vi.updateMesh(pc,vertices=verts)
            self.vi.update()
        
        while not found :            
            print ("attempt ", attempted, marge)
            #main loop thattryto found the next point (similar to jitter)
            if attempted >= safetycutoff:
                print ("break too much attempt ", attempted, safetycutoff)
                return None,False#numpy.array(pt2).flatten()+numpy.array(pt),False
            #pt = numpy.array(self.vi.randpoint_onsphere(self.uLength,biased=(uniform(0.,1.0)*marge)/360.0))*numpy.array([1,1,0])
            test=True
            newPt = None
            if self.prev_alt is not None :
                if dihedral is not None :
                    newPt = self.pickAlternateHalton(pt1,pt2,None)
                elif nextPoint is not None :                                              
                    newPt= self.prev_alt_pt
                if newPt is None :
                    print ("no sphere points available with prev_alt",dihedral,nextPoint)
                    self.prev_alt = None
                    return None,False
                    attempted+=1#?
                    continue
            elif alternate :
                if marge_in is not None :
                    newPt = self.pickAlternateHalton(pt1,pt2,length)
                    if newPt is None :
                        print ("no sphere points available with marge_in")
                        return None,False
                    jtrans,rotMatj = self.get_alternate_position(alternate,ia,v,pt2,newPt)
                elif nextPoint is not None :
                    newPts,jtrans,rotMatj = self.place_alternate(alternate,ia,v,pt1,pt2)    
                    #found = self.check_alternate_collision(pt2,newPts,jtrans,rotMatj)
                    newPt = newPts[0]
                    #attempted should be to safetycutoff
                    attempted=safetycutoff
                    if newPt is None :
                        print ("no  points available place_alternate")
                        return None,False
                else : #no constraint we just place alternate relatively to the given halton new points
                    newPt = self.pickAlternateHalton(pt1,pt2,length)
                    if newPt is None :
                        print ("no sphere points available with marge_in")
                        return None,False
                    jtrans,rotMatj = self.get_alternate_position(alternate,ia,v,pt2,newPt)                        
            elif self.useHalton:
                newPt = self.pickHalton(pt1,pt2)
                if newPt is None :
                    print ("no sphere points available with marge_in")
                    return None,False
            else :
                newPt = self.pickRandomSphere(pt1,pt2,marge,v)                 
            if histoVol.runTimeDisplay :
                self.vi.setTranslation(sp,newPt)
                self.vi.update()
            print ("picked point",newPt)
            if newPt is None :
                print ("no  points available")
                return None,False
            r=[False]
            test =  self.testPoint(newPt)
            print ("testPoint",test,self.constraintMarge,marge,attempted,self.rejectionThreshold)
            if test :
                if not self.constraintMarge :
                    if marge >=175 :
                        attempted +=1
                        continue
                        #print ("no second point not constraintMarge 1 ", marge)
                        #self.prev_alt = None
                        #return None,False
                    if attempted % (self.rejectionThreshold/3) == 0 and not alternate:
                        marge+=1
                        attempted = 0
                        #need to recompute the mask
                        if not alternate and self.useHalton and self.prev_alt is None:
                            self.sphere_points_mask=numpy.ones(self.sphere_points_nb,'i') 
                            self.mask_sphere_points(v,pt2,marge+2.0,liste_nodes,0)
                            if histoVol.runTimeDisplay:
                                points_mask=numpy.nonzero(self.sphere_points_mask)[0]
                                verts=self.sphere_points[points_mask]*self.uLength+pt2            
                                name = "Hcloud"+self.name
                                pc = self.vi.getObject(name)
                                if pc is None :
                                    pc=self.vi.PointCloudObject(name,vertices=verts)[0]
                                else :
                                    self.vi.updateMesh(pc,vertices=verts)
                                self.vi.update()
                attempted +=1
                print ("rejected boundary")
                continue
            if checkcollision:
                collision=False
                cutoff=histoVol.largestProteinSize+self.uLength
                if not alternate :
                    prev=None
                    if len(histoVol.rTrans) > 2:
                        prev = histoVol.rTrans[-1]
                
                    a=numpy.array(newPt)-numpy.array(pt2).flatten()
                    b=numpy.array(pt2).flatten()+a
                    #this s where we use panda
                    rotMatj=[[ 1.,  0.,  0.,  0.],
                       [ 0.,  1.,  0.,  0.],
                       [ 0.,  0.,  1.,  0.],
                       [ 0.,  0.,  0.,  1.]]
                    jtrans = [0.,0.,0.]
                    #move it or generate it inplace
#                    oldpos1=self.positions
#                    oldpos2=self.positions2
#                    self.positions=[[numpy.array(pt2).flatten()],]
#                    self.positions2=[[newPt],]
#                    if self.use_rbsphere :
#                        print ("new RB ")
#                        rbnode = self.addRBsegment(numpy.array(pt2).flatten(),newPt)
#                    else :
#                        rbnode = histoVol.callFunction(histoVol.addRB,(self, numpy.array(jtrans), numpy.array(rotMatj),),{"rtype":self.Type},)#cylinder
#                    #histoVol.callFunction(histoVol.moveRBnode,(rbnode, jtrans, rotMatj,))
                    #if inside organelle check for collision with it ?
#                    self.positions=oldpos1
#                    self.positions2=oldpos2
                    rotMatj,jtrans=self.getJtransRot(numpy.array(pt2).flatten(),newPt)
                    rbnode = self.get_rb_model()
                    self.histoVol.callFunction(self.histoVol.moveRBnode,(rbnode, jtrans, rotMatj,))
                    if len(self.histoVol.rTrans) == 0 : r=[False]
                    else :
                        prev=None
                        if len(self.histoVol.rTrans) > 2 :
                            prev=self.histoVol.rTrans[-1]
                        #print("getClosestIngredient",b/2.)
                        closesbody_indice = self.getClosestIngredient(newPt,histoVol,cutoff=cutoff)#vself.radii[0][0]*2.0
                        if len(closesbody_indice) == 0: 
                            print ("No CloseBody")                            
                            r =[False]         #closesbody_indice[0] == -1            
                        else : 
    #                            print("get_rbNodes",closesbody_indice)
                            print ("collision get RB ",len(closesbody_indice))
                            liste_nodes = self.get_rbNodes(closesbody_indice,
                                            jtrans,prevpoint=prev,getInfo=True)
                            if usePP :
                                #use self.grab_cb and self.pp_server
                                ## Divide the task or just submit job
                                n=0
                                self.histoVol.grab_cb.reset()
                                for i in range(len(liste_nodes)/autopack.ncpus):
                                    for c in range(autopack.ncpus):
                                        self.histoVol.pp_server.submit(self.histoVol.world.contactTestPair, 
                                                              (rbnode, liste_nodes[n][0]), 
                                                callback=self.histoVol.grab_cb.grab)
                                        n+=1
                                    self.histoVol.pp_server.wait()
                                    r.extend(self.histoVol.grab_cb.collision[:])
                                    if True in r :
                                        break
                            else :
                                for node in liste_nodes:
                                    print ("collision test with ",node)
                                    self.histoVol.moveRBnode(node[0], node[1], node[2])
                                    col = (self.histoVol.world.contactTestPair(rbnode, node[0]).getNumContacts() > 0 )
                                    print ("collision? ",col)
                                    r=[col]
                                    if col :
                                        break
                    collision=( True in r)
                    if not collision :
                        self.alternate_interval+=1
                        if self.alternate_interval >= self.mini_interval:
                            self.prev_was_alternate = False
                        self.prev_alt = None
                        self.prev_vec = None
                        self.update_data_tree(numpy.array(pt2).flatten(),rotMatj,pt1=pt2,pt2=newPt)#jtrans
#                        histoVol.callFunction(histoVol.delRB,(rbnode,))
                        return newPt,True

                else :
                    print ("alternate collision")
                    rotMatj1,jtrans1=self.getJtransRot_r(numpy.array(pt2).flatten(),newPt) 
                    #collision,liste_nodes = self.collision_rapid(jtrans1,rotMatj1,cutoff=cutoff,usePP=usePP,point=newPt)
                    #the collision shouldnt look for previous cylinder   
                    collision_alternate,liste_nodes=self.partners[alternate].ingr.pandaBullet_collision(jtrans,rotMatj,None,getnodes=True)                      
                    collision = collision_alternate#(collision or collision_alternate)
#                        print "collision",collision,collision_alternate,len(liste_nodes)
                    if not collision :
                        #what about point in curve and result
                        #self.update_data_tree(jtrans1,rotMatj1,pt1=pt2,pt2=newPt)
                        #self.update_data_tree(jtrans1,rotMatj1,pt1=newPt,pt2=newPts[1])
                        self.partners[alternate].ingr.update_data_tree(jtrans,rotMatj)
                        self.compartment.molecules.append([ jtrans, rotMatj.transpose(), self.partners[alternate].ingr, 0 ])     
                        newv,d1 = self.vi.measure_distance(pt2,newPt,vec=True)
                        #v,d2 = self.vi.measure_distance(newPt,newPts[1],vec=True)
                        #self.currentLength += d1
                        if dihedral is not None :
                            self.prev_alt=alternate
                        self.prev_vec = v
                        if nextPoint is not None and dihedral is None:                                              
                            self.prev_alt_pt = newPts[1]
                        #treat special case of starting other ingr
                        start_ingr_name = self.partners[alternate].getProperties("st_ingr")
                        if start_ingr_name is not None :
                            #add new starting positions
                            start_ingr = self.histoVol.getIngrFromName(start_ingr_name)
                            matrice = numpy.array(rotMatj)#.transpose()
                            matrice[3,:3] = jtrans
                            newPts = self.get_alternate_starting_point(numpy.array(pt2).flatten(),newPt,alternate)
                            start_ingr.start_positions.append([newPts[0],newPts[1]])
                            start_ingr.nbMol+=1
                            #add a mol
                        #we need to store 
                        self.alternate_interval = 0    
                        return newPt,True
                    else : 
                        self.prev_alt_pt = None                    #print (" collide ?",collision)
                if collision : #increment the range
                    if not self.constraintMarge :
                        if marge >=180 : #pi
                            attempted +=1
                            continue
                            #print ("no second point not constraintMarge 2 ", marge)
                            #return None,False
#                            print ("upate marge because collision ", marge)
                        if attempted %  (self.rejectionThreshold/3) == 0 and not alternate:
                            marge+=1
                            attempted = 0 
                            #need to recompute the mask
                            if not alternate and self.useHalton and self.prev_alt is None:
                                self.sphere_points_mask=numpy.ones(self.sphere_points_nb,'i') 
                                self.mask_sphere_points(v,pt2,marge+2.0,liste_nodes,0)
                                if self.runTimeDisplay:
                                    points_mask=numpy.nonzero(self.sphere_points_mask)[0]
                                    v=self.sphere_points[points_mask]*self.uLength+pt2            
                                    name = "Hcloud"+self.name
                                    sp = self.vi.getObject(name)
                                    if sp is None :
                                        pc=self.vi.PointCloudObject("bbpoint",vertices=v)[0]
                                    else :
                                        self.vi.updateMesh(pc,vertices=v)
                                    self.vi.update()
                        else :
                            attempted +=1
                    else :
                        attempted +=1
                    print ("rejected collision")
                    continue
            else :
                found = True
#                histoVol.callFunction(histoVol.delRB,(rbnode,))
                return numpy.array(pt2).flatten()+numpy.array(pt),True
            print ("end loop add attempt ",attempted)
            attempted += 1
#        histoVol.callFunction(histoVol.delRB,(rbnode,))
        return numpy.array(pt2).flatten()+numpy.array(pt),True        
        
    def walkSpherePandaOLD(self,pt1,pt2,distance,histoVol,marge = 90.0,
                        checkcollision=True,usePP=False):
        """ use a random point on a sphere of radius uLength, and useCylinder collision on the grid """ 
#        print ("walkSpherePanda ",pt1,pt2  )      
        v,d = self.vi.measure_distance(pt1,pt2,vec=True)
        found = False
        attempted = 0
        pt = [0.,0.,0.]
        angle=0.
        safetycutoff = 1000
        mask=None
        if self.constraintMarge:
            safetycutoff = 200
        if self.runTimeDisplay:
            name = "walking"+self.name
            sp = self.vi.getObject(name)
            if sp is None :
                sp=self.vi.Sphere(name,radius=2.0)[0]
            #sp.SetAbsPos(self.vi.FromVec(startingPoint))
            self.vi.update()
        liste_nodes=[]
        while not found :            
            print ("attempt ", attempted, marge)
            #main loop thattryto found the next point (similar to jitter)
            if attempted >= safetycutoff:
                print ("break too much attempt ", attempted, safetycutoff)
                return None,False#numpy.array(pt2).flatten()+numpy.array(pt),False
            #pt = numpy.array(self.vi.randpoint_onsphere(self.uLength,biased=(uniform(0.,1.0)*marge)/360.0))*numpy.array([1,1,0])
            test=False
            if self.useHalton:
                self.mask_sphere_points(v,pt2,marge+2.0,liste_nodes,0)
                p = self.getNextPoint()
                if p is None :
                    print ("no sphere points available")
                    return None,False#numpy.array(pt2).flatten()+numpy.array(pt),False                    
                p=numpy.array(p)*self.uLength
                pt = numpy.array(p)*numpy.array(self.jitterMax)#?
                newPt = numpy.array(pt2).flatten()+numpy.array(pt)
                if self.runTimeDisplay  >= 2:
                    self.vi.setTranslation(sp,newPt)
                    self.vi.update()
                test = True
            else :
                p = self.vi.advance_randpoint_onsphere(self.uLength,
                                                   marge=math.radians(marge),
                                                    vector=v)
                pt = numpy.array(p)*numpy.array(self.jitterMax)#?
                #the new position is the previous point (pt2) plus the random point
                newPt = numpy.array(pt2).flatten()+numpy.array(pt)
                if self.runTimeDisplay  >= 2:
                    self.vi.setTranslation(sp,newPt)
                    self.vi.update()
                #compute the angle between the previous direction (pt1->pt2) and the new random one (pt)
                angle=self.vi.angle_between_vectors(numpy.array(v),numpy.array(pt))
                test= abs(math.degrees(angle)) <= marge+2.0
            if test :
                r=[False]
#                print "ok"
                #check if in bounding box
                inComp = True
                closeS = False
                inside = histoVol.grid.checkPointInside(newPt,dist=self.cutoff_boundary,jitter=self.jitterMax)
                if inside :
                    inComp = self.checkPointComp(newPt)
                    if inComp :
                        #check how far from surface ?
                        closeS = self.checkPointSurface(newPt,cutoff=self.cutoff_surface)
                if not inside or closeS or not inComp:
                    print ("inside,closeS ",not inside,closeS,not inComp,newPt,marge)
                    if not self.constraintMarge :
                        if marge >=175 :
                            print ("no second point not constraintMarge 1 ", marge)
                            return None,False
                        #print ("increase marge because inside,closeS ",inside,closeS,inComp,newPt,marge)
                        marge+=1
                    else :
                        print ("inside,closeS ",inside,closeS,inComp,newPt,marge)
                        attempted +=1
                    continue
                #optionally check for collision
                if checkcollision:
                    a=numpy.array(newPt)-numpy.array(pt2).flatten()
                    b=numpy.array(pt2).flatten()+a
                    #this s where we use panda
                    rotMatj=[[ 1.,  0.,  0.,  0.],
                       [ 0.,  1.,  0.,  0.],
                       [ 0.,  0.,  1.,  0.],
                       [ 0.,  0.,  0.,  1.]]
                    jtrans = [0.,0.,0.]
                    #move it or generate it inplace
                    oldpos1=self.positions
                    oldpos2=self.positions2
                    self.positions=[[numpy.array(pt2).flatten()],]
                    self.positions2=[[newPt],]
                    if self.use_rbsphere :
                        print ("new RB ")
                        rbnode = self.addRBsegment(numpy.array(pt2).flatten(),newPt)
                    else :
                        rbnode = histoVol.callFunction(histoVol.addRB,(self, numpy.array(jtrans), numpy.array(rotMatj),),{"rtype":self.Type},)#cylinder
                    #histoVol.callFunction(histoVol.moveRBnode,(rbnode, jtrans, rotMatj,))
                    #if inside organelle check for collision with it ?
                    self.positions=oldpos1
                    self.positions2=oldpos2
                    #check collision using bullet
                    #get closest object?
                    cutoff=histoVol.largestProteinSize+self.uLength
                    if len(self.histoVol.rTrans) == 0 : r=[False]
                    else :
                        prev=None
                        if len(self.histoVol.rTrans) > 2 :
                            prev=self.histoVol.rTrans[-1]
                        #print("getClosestIngredient",b/2.)
                        closesbody_indice = self.getClosestIngredient(newPt,histoVol,cutoff=cutoff)#vself.radii[0][0]*2.0
                        if len(closesbody_indice) == 0: r =[False]         #closesbody_indice[0] == -1            
                        else : 
#                            print("get_rbNodes",closesbody_indice)
                            #print ("collision RB ",len(closesbody_indice))
                            liste_nodes = self.get_rbNodes(closesbody_indice,
                                            jtrans,prevpoint=prev,getInfo=True)
                            if usePP :
                                #use self.grab_cb and self.pp_server
                                ## Divide the task or just submit job
                                n=0
                                self.histoVol.grab_cb.reset()
                                for i in range(len(liste_nodes)/autopack.ncpus):
                                    for c in range(autopack.ncpus):
                                        self.histoVol.pp_server.submit(self.histoVol.world.contactTestPair, 
                                                              (rbnode, liste_nodes[n][0]), 
                                                callback=self.histoVol.grab_cb.grab)
                                        n+=1
                                    self.histoVol.pp_server.wait()
                                    r.extend(self.histoVol.grab_cb.collision[:])
                                    if True in r :
                                        break
                            else :
                                for node in liste_nodes:
                                    col = (self.histoVol.world.contactTestPair(rbnode, node[0]).getNumContacts() > 0 )
                                    r=[col]
                                    if col :
                                        break
                                #r=[ (self.histoVol.world.contactTestPair(rbnode, node).getNumContacts() > 0 ) for node in liste_nodes]# in closesbody_indice if n != len(self.histoVol.rTrans)-1]  #except last one  that should be last drop fragment                    
                    #closesbody_indice = self.getClosestIngredient(newPt,histoVol,cutoff=self.uLength)
                    #r=[ (histoVol.world.contactTestPair(rbnode, self.histoVol.static[n]).getNumContacts() > 0 ) for n in closesbody_indice]  #except last one  that should be last drop fragment
                    #result2 = histoVol.world.contactTest(rbnode)
                    #r = [( result2.getNumContacts() > 0),]    
#                    print ("contact All ",(True in r))                 
                    collision=( True in r)
                    #print (" collide ?",collision)
                    if not collision :
                        #print angle,math.degrees(angle),marge
                        histoVol.static.append(rbnode)
                        histoVol.moving = None
                        found = True
#                        histoVol.close_ingr_bhtree.MoveRBHPoint(histoVol.nb_ingredient,jtrans,0)
                        histoVol.nb_ingredient+=1
#                        r,j=self.getJtransRot(numpy.array(pt2).flatten(),newPt)
                        histoVol.rTrans.append(numpy.array(pt2).flatten())
#                        m=autopack.helper.getTubePropertiesMatrix(numpy.array(pt2).flatten(),newPt)[1]
                        histoVol.rRot.append(numpy.array(rotMatj))#rotMatj r 
                        histoVol.rIngr.append(self)
                        histoVol.result.append([ [numpy.array(pt2).flatten(),newPt], rotMatj, self, 0 ])
                        histoVol.callFunction(histoVol.delRB,(rbnode,))
                        #histoVol.close_ingr_bhtree.InsertRBHPoint((jtrans[0],jtrans[1],jtrans[2]),radius,None,histoVol.nb_ingredient)
#                        print ("update bhtree")
                        if histoVol.treemode == "bhtree":# "cKDTree"
                            if len(histoVol.rTrans) > 1 : bhtreelib.freeBHtree(histoVol.close_ingr_bhtree)
                            histoVol.close_ingr_bhtree=bhtreelib.BHtree( histoVol.rTrans, None, 10)
                        else :
                            #rebuild kdtree
                            if len(self.histoVol.rTrans) > 1 :histoVol.close_ingr_bhtree= spatial.cKDTree(histoVol.rTrans, leafsize=10)
#                        print ("bhtree updated")
                        return numpy.array(pt2).flatten()+numpy.array(pt),True
                    else : #increment the range
                        if not self.constraintMarge :
                            if marge >=180 : #pi
                                #print ("no second point not constraintMarge 2 ", marge)
                                return None,False
                            #print ("upate marge because collision ", marge)
                            marge+=1
                        else :
                            #print ("collision")
                            attempted +=1
                        continue
                else :
                    found = True
                    print ("found !")
                    histoVol.callFunction(histoVol.delRB,(rbnode,))
                    return numpy.array(pt2).flatten()+numpy.array(pt),True
#                attempted += 1
            else :
                print ("not in the marge ",abs(math.degrees(angle)),marge)
                attempted += 1
                continue
#            print ("end loop add attempt ",attempted)
            attempted += 1
        histoVol.callFunction(histoVol.delRB,(rbnode,))
        return numpy.array(pt2).flatten()+numpy.array(pt),True
        

    def walkSphereRAPIDold(self,pt1,pt2,distance,histoVol,marge = 90.0,checkcollision=True,usePP=False):
        """ use a random point on a sphere of radius uLength, and useCylinder collision on the grid """ 
#        print ("walkSphereRapid ",pt1,pt2  )      
        v,d = self.vi.measure_distance(pt1,pt2,vec=True)
        found = False
        attempted = 0
        pt = [0.,0.,0.]
        angle=0.
        safetycutoff = 50
        if self.constraintMarge:
            safetycutoff = 50
        if self.runTimeDisplay:
            name = "walking"+self.name
            sp = self.vi.getObject(name)
            if sp is None :
                sp=self.vi.Sphere(name,radius=2.0)[0]
            #sp.SetAbsPos(self.vi.FromVec(startingPoint))
            self.vi.update()
        #do we use the cylindr or the alternate / partner
        liste_nodes=[]     
        cutoff=self.histoVol.largestProteinSize+self.uLength
        closesbody_indice = self.getClosestIngredient(pt2,self.histoVol,cutoff=cutoff)
        liste_nodes = self.get_rapid_nodes(closesbody_indice,pt2,prevpoint=pt1)  
        #mask using boundary and ingredient 
        #self.mask_sphere_points_boundary(pt2)
        #self.mask_sphere_points_ingredients(pt2,liste_nodes)
        #mask_sphere_points_start = self.sphere_points_mask[:]
        while not found :
            self.sphere_points_mask = numpy.ones(10000,'i') #mask_sphere_points_start[:]
            dihedral= None
            nextPoint = None
            #liste_nodes=[]
            print ("attempt ", attempted, marge)
            #main loop thattryto found the next point (similar to jitter)
            if attempted >= safetycutoff:   
                print ("break too much attempt ", attempted, safetycutoff)
                return None,False#numpy.array(pt2).flatten()+numpy.array(pt),False
            #pt = numpy.array(self.vi.randpoint_onsphere(self.uLength,biased=(uniform(0.,1.0)*marge)/360.0))*numpy.array([1,1,0])
            test=False
            newPt = None
            alternate,ia = self.pick_alternate()
            print ("pick alternate",alternate,ia,self.prev_alt)
            #thats the ame of the ingedient used as alternate
            if self.prev_alt is not None :#and self.prev_alt_pt is not None:
#                print ("prev_alt",self.prev_alt,self.prev_alt_pt)
                newPt = self.prev_alt_pt
                test=True     
                if self.prev_alt_pt is not None :
                    self.prev_alt_pt=None 
                alternate = None
            t1 = time()
            if newPt is not None :
                test = True
#            elif autopack.helper.measure_distance(pt1,pt2) == 0.0:
#                return None,False
            elif alternate :
#                print ("try to place alernate",alternate,ia) 
                p_alternate = self.partners[alternate]#self.histoVol.getIngrFromNameInRecipe(alternate,self.recipe )
                #if p_alternate.getProperties("bend"):
                nextPoint=p_alternate.getProperties("pt2")
                marge_in=p_alternate.getProperties("marge_in")            
                dihedral=p_alternate.getProperties("diehdral")#for next point
                length=p_alternate.getProperties("length")#for next point
                if marge_in is not None :
                    self.mask_sphere_points(v,pt2,marge_in,liste_nodes,0,pv=None,marge_diedral=None)
                    p = self.getNextPoint()
                    if p is None :
                        print ("no sphere points available with marge_in")
                        #try again ?
                        attempted +=1
                        continue
                        #return None,False
                    if length is not None:
                        p=(p/self.uLength)*length
                    newPt = numpy.array(pt2).flatten()+numpy.array(p)
                    jtrans,rotMatj = self.get_alternate_position(alternate,ia,v,pt2,newPt)
                else :
                    newPts,jtrans,rotMatj = self.place_alternate(alternate,ia,v,pt1,pt2)    
                    #found = self.check_alternate_collision(pt2,newPts,jtrans,rotMatj)
                    newPt = newPts[0]
                test=True
            elif self.useHalton: 
#                print ("halton point")
                self.mask_sphere_points(v,pt2,marge+2.0,liste_nodes,0)
                p = self.getNextPoint()
                if p is None :
                    print ("no sphere points available with halton",pt2,v)
                    marge+=1
                    attempted +=1
                    continue
                    #return None,False#numpy.array(pt2).flatten()+numpy.array(pt),False                    
                p=numpy.array(p)#*self.uLength
                pt = numpy.array(p)#*numpy.array(self.jitterMax)#?
                newPt = numpy.array(pt2).flatten()+numpy.array(pt)
#                print ("pick halton",pt2,p,newPt)
                test = True
            else :
                p = self.vi.advance_randpoint_onsphere(self.uLength,
                                                       marge=math.radians(marge),
                                                       vector=v)
                pt = numpy.array(p)*numpy.array(self.jitterMax)
                #the new position is the previous point (pt2) plus the random point
                newPt = numpy.array(pt2).flatten()+numpy.array(pt)
                #compute the angle between the previous direction (pt1->pt2) and the new random one (pt)
                angle=self.vi.angle_between_vectors(numpy.array(v),numpy.array(pt))
                test= abs(math.degrees(angle)) <= marge+2.0
            if self.runTimeDisplay  >= 2:
                self.vi.setTranslation(sp,newPt)
                self.vi.update()
            print ("time to pick point ",time()-t1)
            if test :
                r=[False]
                test =  self.testPoint(newPt)
#                if alternate or self.prev_alt is not None :
#                    print ("do nothing")
#                    test2 = self.testPoint(newPts[1])
                    #collide 
                    #print ("problem2")
#                    if test :#r test2:
#                        print ("problem")
#                        return None,False
#                else :
                if test :
    #                    print "inside,closeS ",not inside,closeS,not inComp,newPt,marge
                        if not self.constraintMarge :
                            if marge >=175 :
                                print ("no second point not constraintMarge 1 ", marge)
                                return None,False
    #                        print ("increase marge because inside,closeS ",inside,closeS,inComp,newPt,marge)
                            marge+=1
                        else :
#                             print "newAttempt"
#                            print ("inside,closeS ",inside,closeS,inComp,newPt,marge)
                            attempted +=1
                        print ("rejected boundary")
                        continue
                #optionally check for collision
                if checkcollision:  
#                    print ("check collision",alternate)                     
                    collision=False
                    cutoff=histoVol.largestProteinSize+self.uLength
                    if not alternate :
                        prev=None
                        if len(histoVol.rTrans) > 2:
                            prev = histoVol.rTrans[-1]
                        a=numpy.array(newPt)-numpy.array(pt2).flatten()
                        b=numpy.array(pt2).flatten()+a
                        #this s where we use panda
                        rotMatj=[[ 1.,  0.,  0.,  0.],
                           [ 0.,  1.,  0.,  0.],
                           [ 0.,  0.,  1.,  0.],
                           [ 0.,  0.,  0.,  1.]]
                        jtrans = [0.,0.,0.]
    #                    rbnode = self.get_rapid_model()
                        rotMatj,jtrans=self.getJtransRot(numpy.array(pt2).flatten(),newPt)                    
                        collision,liste_nodes = self.collision_rapid(jtrans,rotMatj,
                                        cutoff=cutoff,usePP=usePP,point=newPt,
                                        prevpoint=prev)
#                        print "collision",collision,cutoff
#                        raw_input()
                        if not collision :
                            self.prev_alt = None
                            self.update_data_tree(jtrans,rotMatj,pt1=pt2,pt2=newPt)
                            return newPt,True
                    else :
                        rotMatj1,jtrans1=self.getJtransRot_r(numpy.array(pt2).flatten(),newPt) 
                        collision,liste_nodes = self.collision_rapid(jtrans1,rotMatj1,cutoff=cutoff,usePP=usePP,point=newPt)
                        #the collision shouldnt look for previous cylinder                        
                        collision_alternate,liste_nodes = self.partners[alternate].ingr.collision_rapid(jtrans,rotMatj,usePP=usePP,liste_nodes=liste_nodes)
                        collision = collision_alternate#(collision or collision_alternate)
#                        print "collision",collision,collision_alternate,len(liste_nodes)
                        if not collision :
                            #what about point in curve and result
                            #self.update_data_tree(jtrans1,rotMatj1,pt1=pt2,pt2=newPt)
                            #self.update_data_tree(jtrans1,rotMatj1,pt1=newPt,pt2=newPts[1])
                            self.partners[alternate].ingr.update_data_tree(jtrans,rotMatj)
                            self.compartment.molecules.append([ jtrans, rotMatj.transpose(), self.partners[alternate].ingr, 0 ])     
                            newv,d1 = self.vi.measure_distance(pt2,newPt,vec=True)
                            #v,d2 = self.vi.measure_distance(newPt,newPts[1],vec=True)
                            #self.currentLength += d1
                            self.prev_alt=alternate
                            if dihedral is not None :
                                self.mask_sphere_points(newv,pt2,marge_in,liste_nodes,0,pv=v,marge_diedral=dihedral)
                                p = self.getNextPoint()
                                self.prev_alt_pt = numpy.array(newPt).flatten()+numpy.array(pt)
                            elif nextPoint is not None :                                              
                                self.prev_alt_pt = newPts[1]
                            #treat special case of starting other ingr
                            start_ingr_name = self.partners[alternate].getProperties("st_ingr")
                            if start_ingr_name is not None :
                                #add new starting positions
                                start_ingr = self.histoVol.getIngrFromName(start_ingr_name)
                                matrice = numpy.array(rotMatj)#.transpose()
                                matrice[3,:3] = jtrans
                                newPts = self.get_alternate_starting_point(matrice,alternate)
                                start_ingr.start_positions.append([newPts[0],newPts[1]])
                            return newPt,True
                        else : 
                            self.prev_alt_pt = None
                    if collision : #increment the range
                        if not self.constraintMarge :
                            if marge >=180 : #pi
                                print ("no second point not constraintMarge 2 ", marge)
                                return None,False
#                            print ("upate marge because collision ", marge)
                            marge+=1
                        else :
#                            print ("collision")
                            attempted +=1
                        print ("rejected collision")
                        continue
                else :
                    found = True
#                    print ("found !")
                    return numpy.array(pt2).flatten()+numpy.array(pt),True
#                attempted += 1
            else :
                print ("not in the marge ",abs(math.degrees(angle)),marge)
                attempted += 1
                continue
            print ("end loop add attempt ",attempted)
            attempted += 1
        return numpy.array(pt2).flatten()+numpy.array(pt),True

    def pickHalton(self,pt1,pt2):
        p = self.getNextPoint()
        if p is None :
            return None                
        p=numpy.array(p)#*self.uLength
        pt = numpy.array(p)#*numpy.array(self.jitterMax)#?
        return numpy.array(pt2).flatten()+numpy.array(pt)

    def pickRandomSphere(self,pt1,pt2,marge,v):
        p = self.vi.advance_randpoint_onsphere(self.uLength,
                                               marge=math.radians(marge),
                                               vector=v)
        pt = numpy.array(p)*numpy.array(self.jitterMax)
        #the new position is the previous point (pt2) plus the random point
        newPt = numpy.array(pt2).flatten()+numpy.array(pt)
        #compute the angle between the previous direction (pt1->pt2) and the new random one (pt)
#        angle=self.vi.angle_between_vectors(numpy.array(v),numpy.array(pt))
#        test= abs(math.degrees(angle)) <= marge+2.0
        return newPt
        
    def pickAlternateHalton(self,pt1,pt2,length):
        p = self.getNextPoint()
        if p is None :
            return None
            #return None,False
        if length is not None:
            p=(p/self.uLength)*length
        newPt = numpy.array(pt2).flatten()+numpy.array(p)
        return newPt
       
    def walkSphereRAPID(self,pt1,pt2,distance,histoVol,marge = 90.0,checkcollision=True,usePP=False):
        """ use a random point on a sphere of radius uLength, and useCylinder collision on the grid """ 
#        print ("walkSphereRapid ",pt1,pt2  )      
        v,d = self.vi.measure_distance(pt1,pt2,vec=True)
        found = False
        attempted = 0
        pt = [0.,0.,0.]
        angle=0.
        safetycutoff = self.rejectionThreshold#angle  / 360
        mask = None
        sp=None
        pc=None
        if self.constraintMarge:
            safetycutoff = self.rejectionThreshold
        if histoVol.runTimeDisplay:
            name = "walking"+self.name
            sp = self.vi.getObject(name)
            if sp is None :
                sp=self.vi.Sphere(name,radius=2.0)[0]
            #sp.SetAbsPos(self.vi.FromVec(startingPoint))
            self.vi.update()
        #do we use the cylinder or the alternate / partner
        liste_nodes=[]     
        cutoff=self.histoVol.largestProteinSize+self.uLength
        closesbody_indice = self.getClosestIngredient(pt2,self.histoVol,cutoff=cutoff)
        liste_nodes = self.get_rapid_nodes(closesbody_indice,pt2,prevpoint=pt1)  
        #mask using boundary and ingredient 
        #self.mask_sphere_points_boundary(pt2)
        #self.mask_sphere_points_ingredients(pt2,liste_nodes)
        #mask_sphere_points_start = self.sphere_points_mask[:]
        alternate,ia = self.pick_alternate()
        print ("pick alternate",alternate,ia,self.prev_alt)
        if self.prev_was_alternate :
            alternate = None
        p_alternate = None#self.partners[alternate]#self.histoVol.getIngrFromNameInRecipe(alternate,self.recipe )
        #if p_alternate.getProperties("bend"):
        nextPoint= None#p_alternate.getProperties("pt2")
        marge_in= None#p_alternate.getProperties("marge_in")            
        dihedral= None#p_alternate.getProperties("diehdral")#for next point
        length= None#p_alternate.getProperties("length")#for next point
        #prepare halton if needed
        newPt = None
        newPts=[]
        #thats the name of the ingedient used as alternate
        if self.prev_alt is not None :#and self.prev_alt_pt is not None:
            p_alternate = self.partners[self.prev_alt]
            dihedral= p_alternate.getProperties("diehdral")
            nextPoint= p_alternate.getProperties("pt2")
            marge_in= p_alternate.getProperties("marge_out")  #marge out ? 
            if dihedral is not None :
                self.mask_sphere_points(v,pt1,marge_in,liste_nodes,0,
                                        pv=self.prev_vec,marge_diedral=dihedral)
            alternate = None
        t1 = time()
        if self.prev_alt is not None :
            test = True
        elif alternate and not self.prev_was_alternate :
            #next point shouldnt be an alternate
            p_alternate = self.partners[alternate]#self.histoVol.getIngrFromNameInRecipe(alternate,self.recipe )
            #if p_alternate.getProperties("bend"):
            nextPoint=p_alternate.getProperties("pt2")
            marge_in=p_alternate.getProperties("marge_in")            
            dihedral=p_alternate.getProperties("diehdral")#for next point
            length=p_alternate.getProperties("length")#for next point
            if marge_in is not None :
                self.mask_sphere_points(v,pt2,marge_in,liste_nodes,0,pv=None,marge_diedral=None)
            else :
                self.mask_sphere_points(v,pt2,marge+2.0,liste_nodes,0) 
            self.prev_was_alternate = True
        elif self.useHalton: 
            self.mask_sphere_points(v,pt2,marge+2.0,liste_nodes,0)
        if histoVol.runTimeDisplay:
            points_mask=numpy.nonzero(self.sphere_points_mask)[0]
            verts=self.sphere_points[points_mask]*self.uLength+pt2            
            name = "Hcloud"+self.name
            pc = self.vi.getObject(name)
            if pc is None :
                pc=self.vi.PointCloudObject(name,vertices=verts)[0]
            else :
                self.vi.updateMesh(pc,vertices=verts)
            self.vi.update()

        while not found :
            #try to drop at newPoint, 
#            self.sphere_points_mask = numpy.ones(10000,'i') #mask_sphere_points_start[:]
            #dihedral= None
            #nextPoint = None
            #liste_nodes=[]
            print ("attempt ", attempted, marge)
            #main loop thattryto found the next point (similar to jitter)
            if attempted > safetycutoff:   
                print ("break too much attempt ", attempted, safetycutoff)
                return None,False#numpy.array(pt2).flatten()+numpy.array(pt),False
            #pt = numpy.array(self.vi.randpoint_onsphere(self.uLength,biased=(uniform(0.,1.0)*marge)/360.0))*numpy.array([1,1,0])
            test=True
            newPt = None
            if self.prev_alt is not None :
                if dihedral is not None :
                    newPt = self.pickAlternateHalton(pt1,pt2,None)
                elif nextPoint is not None :                                              
                    newPt= self.prev_alt_pt
                if newPt is None :
                    print ("no sphere points available with prev_alt",dihedral,nextPoint)
                    self.prev_alt = None
                    return None,False
                    attempted+=1
                    continue
            elif alternate :
                if marge_in is not None :
                    newPt = self.pickAlternateHalton(pt1,pt2,length)
                    if newPt is None :
                        print ("no sphere points available with marge_in")
                        return None,False
                    jtrans,rotMatj = self.get_alternate_position(alternate,ia,v,pt2,newPt)
                elif nextPoint is not None :
                    newPts,jtrans,rotMatj = self.place_alternate(alternate,ia,v,pt1,pt2)    
                    #found = self.check_alternate_collision(pt2,newPts,jtrans,rotMatj)
                    newPt = newPts[0]
                    if newPt is None :
                        print ("no  points available place_alternate")
                        return None,False
                else : #no constraint we just place alternate relatively to the given halton new points
                    newPt = self.pickAlternateHalton(pt1,pt2,length)
                    if newPt is None :
                        print ("no sphere points available with marge_in")
                        return None,False
                    jtrans,rotMatj = self.get_alternate_position(alternate,ia,v,pt2,newPt)                        
            elif self.useHalton:
                newPt = self.pickHalton(pt1,pt2)
                if newPt is None :
                    print ("no sphere points available with marge_in")
                    return None,False
            else :
                newPt = self.pickRandomSphere(pt1,pt2,marge,v)                 
            if histoVol.runTimeDisplay :
                self.vi.setTranslation(sp,newPt)
                self.vi.update()
            print ("picked point",newPt)
            if newPt is None :
                print ("no  points available")
                return None,False
            r=[False]
            test =  self.testPoint(newPt)
            if test :
                if not self.constraintMarge :
                    if marge >=175 :
                        attempted +=1
                        continue
                        #print ("no second point not constraintMarge 1 ", marge)
                        #self.prev_alt = None
                        #return None,False
                    if attempted % (self.rejectionThreshold/3) == 0 :
                        marge+=1
                        attempted = 0
                        #need to recompute the mask
                        if not alternate and self.useHalton and self.prev_alt is None:
                            self.sphere_points_mask=numpy.ones(self.sphere_points_nb,'i') 
                            self.mask_sphere_points(v,pt2,marge+2.0,liste_nodes,0)
                            if histoVol.runTimeDisplay:
                                points_mask=numpy.nonzero(self.sphere_points_mask)[0]
                                verts=self.sphere_points[points_mask]*self.uLength+pt2            
                                name = "Hcloud"+self.name
                                pc = self.vi.getObject(name)
                                if pc is None :
                                    pc=self.vi.PointCloudObject(name,vertices=verts)[0]
                                else :
                                    self.vi.updateMesh(pc,vertices=verts)
                                self.vi.update()
                attempted +=1
                print ("rejected boundary")
                continue
            #optionally check for collision
            if checkcollision:  
#                    print ("check collision",alternate)                     
                collision=False
                cutoff=histoVol.largestProteinSize+self.uLength
                if not alternate :
                    prev=None
                    if len(histoVol.rTrans) > 2:
                        prev = histoVol.rTrans[-1]
                    a=numpy.array(newPt)-numpy.array(pt2).flatten()
                    b=numpy.array(pt2).flatten()+a
                    #this s where we use panda
                    rotMatj=[[ 1.,  0.,  0.,  0.],
                       [ 0.,  1.,  0.,  0.],
                       [ 0.,  0.,  1.,  0.],
                       [ 0.,  0.,  0.,  1.]]
                    jtrans = [0.,0.,0.]
#                    rbnode = self.get_rapid_model()
                    rotMatj,jtrans=self.getJtransRot(numpy.array(pt2).flatten(),newPt)                    
                    collision,liste_nodes = self.collision_rapid(jtrans,rotMatj,
                                    cutoff=cutoff,usePP=usePP,point=newPt,
                                    prevpoint=prev)
#                        print "collision",collision,cutoff
#                        raw_input()
                    if not collision :
                        self.alternate_interval+=1
                        if self.alternate_interval >= self.mini_interval:
                            self.prev_was_alternate = False
                        self.prev_alt = None
                        self.prev_vec = None
                        self.update_data_tree(jtrans,rotMatj,pt1=pt2,pt2=newPt)
                        return newPt,True
                else :
                    rotMatj1,jtrans1=self.getJtransRot_r(numpy.array(pt2).flatten(),newPt) 
                    #collision,liste_nodes = self.collision_rapid(jtrans1,rotMatj1,cutoff=cutoff,usePP=usePP,point=newPt)
                    #the collision shouldnt look for previous cylinder                        
                    collision_alternate,liste_nodes = self.partners[alternate].ingr.collision_rapid(jtrans,rotMatj,usePP=usePP)#,liste_nodes=liste_nodes)
                    collision = collision_alternate#(collision or collision_alternate)
#                        print "collision",collision,collision_alternate,len(liste_nodes)
                    if not collision :
                        #what about point in curve and result
                        #self.update_data_tree(jtrans1,rotMatj1,pt1=pt2,pt2=newPt)
                        #self.update_data_tree(jtrans1,rotMatj1,pt1=newPt,pt2=newPts[1])
                        self.partners[alternate].ingr.update_data_tree(jtrans,rotMatj)
                        self.compartment.molecules.append([ jtrans, rotMatj.transpose(), self.partners[alternate].ingr, 0 ])     
                        newv,d1 = self.vi.measure_distance(pt2,newPt,vec=True)
                        #v,d2 = self.vi.measure_distance(newPt,newPts[1],vec=True)
                        #self.currentLength += d1
                        if dihedral is not None :
                            self.prev_alt=alternate
                        self.prev_vec = v
                        if nextPoint is not None and dihedral is None:                                              
                            self.prev_alt_pt = newPts[1]
                        #treat special case of starting other ingr
                        start_ingr_name = self.partners[alternate].getProperties("st_ingr")
                        if start_ingr_name is not None :
                            #add new starting positions
                            start_ingr = self.histoVol.getIngrFromName(start_ingr_name)
                            matrice = numpy.array(rotMatj)#.transpose()
                            matrice[3,:3] = jtrans
                            snewPts = self.get_alternate_starting_point(matrice,alternate)
                            start_ingr.start_positions.append([snewPts[0],snewPts[1]])
                            start_ingr.nbMol+=1
                            #add a mol
                        #we need to store 
                        self.alternate_interval = 0    
                        return newPt,True
                    else : 
                        self.prev_alt_pt = None
                if collision : #increment the range
                    if not self.constraintMarge :
                        if marge >=180 : #pi
                            attempted +=1
                            continue
                            #print ("no second point not constraintMarge 2 ", marge)
                            #return None,False
#                            print ("upate marge because collision ", marge)
                        if attempted %  (self.rejectionThreshold/3) == 0 :
                            marge+=1
                            attempted = 0 
                            #need to recompute the mask
                            if not alternate and self.useHalton and self.prev_alt is None:
                                self.sphere_points_mask=numpy.ones(self.sphere_points_nb,'i') 
                                self.mask_sphere_points(v,pt2,marge+2.0,liste_nodes,0)
                                if self.runTimeDisplay:
                                    points_mask=numpy.nonzero(self.sphere_points_mask)[0]
                                    v=self.sphere_points[points_mask]*self.uLength+pt2            
                                    name = "Hcloud"+self.name
                                    sp = self.vi.getObject(name)
                                    if sp is None :
                                        pc=self.vi.PointCloudObject("bbpoint",vertices=v)[0]
                                    else :
                                        self.vi.updateMesh(pc,vertices=v)
                                    self.vi.update()
                        else :
                            attempted +=1
                    else :
                        attempted +=1
                    print ("rejected collision")
                    continue
            else :
                found = True
                return numpy.array(pt2).flatten()+numpy.array(pt),True
            print ("end loop add attempt ",attempted)
            attempted += 1
        return numpy.array(pt2).flatten()+numpy.array(pt),True
        
    def resetLastPoint(self,listePtCurve):
        self.histoVol.nb_ingredient-=1
        self.histoVol.rTrans.pop(len(self.histoVol.rTrans)-1)
        self.histoVol.rRot.pop(len(self.histoVol.rRot)-1)#rotMatj 
        self.histoVol.rIngr.pop(len(self.histoVol.rIngr)-1)
        self.histoVol.result.pop(len(self.histoVol.result)-1)
        if self.histoVol.treemode == "bhtree":# "cKDTree"
            if len(self.histoVol.rTrans) > 1 : bhtreelib.freeBHtree(self.histoVol.close_ingr_bhtree)
            if len(self.histoVol.rTrans) : self.histoVol.close_ingr_bhtree=bhtreelib.BHtree( self.histoVol.rTrans, None, 10)
        else :
            #rebuild kdtree
            if len(self.histoVol.rTrans) > 1 :self.histoVol.close_ingr_bhtree= spatial.cKDTree(self.histoVol.rTrans, leafsize=10)

        #also remove from the result ?
        self.results.pop(len(self.results)-1)
        self.currentLength -=self.uLength
        #not enought the information is still here
        listePtCurve.pop(len(listePtCurve)-1)

    def grow(self,previousPoint,startingPoint,secondPoint,listePtCurve,
              listePtLinear,histoVol, ptInd, 
              freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,r=False,usePP=False):
        #r is for reverse growing
        Done = False
        runTimeDisplay = histoVol.runTimeDisplay
        gridPointsCoords = histoVol.grid.masterGridPositions
        if runTimeDisplay:
            parent = histoVol.afviewer.orgaToMasterGeom[self]
        k=0
        success = False
        safetycutoff = self.safetycutoff
#        if self.constraintMarge:
#            safetycutoff = 50
        counter = 0
        mask=None
        if self.walkingMode == "lattice" and self.compNum > 0:
            o = self.histoVol.compartments[abs(self.compNum)-1]
            v=o.surfacePointsCoords
            mask=numpy.ones(len(v),int)
        alternate=False
        previousPoint_store = None
        while not Done:
            #rest the mask 
            self.sphere_points_mask=numpy.ones(10000,'i') 
            alternate=False
            print ("attempt K ",k)
            if k > safetycutoff :
                print ("break safetycutoff",k)
                return success, nbFreePoints,freePoints
            previousPoint_store = previousPoint
            previousPoint = startingPoint
            startingPoint = secondPoint
            if runTimeDisplay :# or histoVol.afviewer.doSpheres:
                name = str(len(listePtLinear))+"sp"+self.name + str(ptInd)
                if r :
                    name = str(len(listePtLinear)+1)+"sp"+self.name + str(ptInd)
                sp=self.vi.Sphere(name,radius=self.radii[0][0],parent=parent)[0]
                self.vi.setTranslation(sp,pos=startingPoint)
#                sp.SetAbsPos(self.vi.FromVec(startingPoint))
                #sp=self.vi.newInstance(name,histoVol.afviewer.pesph,
#                                       location=startingPoint,parent=parent)
                #self.vi.scaleObj(sp,self.radii[0][0])                
                self.vi.update()
#            print "ok,p,start",previousPoint,startingPoint
            #pick next point and test collision.
            if self.walkingMode == "sphere":
                if self.placeType == "pandaBullet" :
                    secondPoint,success = self.walkSpherePanda(previousPoint,startingPoint,
                                                      distance,histoVol,
                                                      marge = self.marge,
                                                      checkcollision=True,usePP=usePP)  
                    if secondPoint is None :
                        return False, nbFreePoints,freePoints
                elif self.placeType == "RAPID" :
                    #call function
                    t1 = time()
                    secondPoint,success = self.walkSphereRAPID(previousPoint,startingPoint,
                                                      distance,histoVol,
                                                      marge = self.marge,
                                                      checkcollision=True,usePP=usePP) 
                    print ("walk rapid ",time()-t1)
                    if secondPoint is None :
                        return False, nbFreePoints,freePoints
                else :
                    secondPoint,success = self.walkSphere(previousPoint,startingPoint,
                                                      distance,histoVol,
                                                      marge = self.marge,
                                                      checkcollision=True)            
            if self.walkingMode == "lattice" and self.compNum > 0:
                secondPoint,success,mask = self.walkLatticeSurface(startingPoint,distance,histoVol, 
                                                  2,mask,marge=self.marge,
                                            checkcollision=False,saw=True)            
            elif self.walkingMode == "lattice":
                secondPoint,success = self.walkLattice(startingPoint,distance,histoVol, 
                                                  2,marge=self.marge,
                                            checkcollision=False,saw=True)
            if secondPoint is None or not success : #no available point? try again ?
                secondPoint = numpy.array(previousPoint)
                startingPoint = previousPoint_store
                k+=1
                continue

            if len(secondPoint) == 2 :
                alternate=True
                startingPoint = secondPoint[0]
                secondPoint = secondPoint[1]
#                print ("should add ",startingPoint,secondPoint)
            #print ("accepted ",success)
            v,d = self.vi.measure_distance(startingPoint,secondPoint,vec=True)
            
            rotMatj,jtrans=self.getJtransRot(startingPoint,secondPoint)
            if r :
                #reverse mode
                rotMatj,jtrans=self.getJtransRot(secondPoint,startingPoint)
            cent1T = self.transformPoints(jtrans, rotMatj, self.positions[-1])
            cent2T = self.transformPoints(jtrans, rotMatj, self.positions2[-1])
            print ("here is output of walk",secondPoint,startingPoint,success,alternate,len(secondPoint))             
            print  (cent1T,cent2T,jtrans, rotMatj)           
            if success:
                print ("success grow")
                #do not append if alternate was used
                if self.prev_alt is None :
                    self.results.append([jtrans,rotMatj])
                if r :
                    if alternate :
                        listePtLinear.insert(0,startingPoint)
                    listePtLinear.insert(0,secondPoint)
                else :
                    if alternate :
                        listePtLinear.append(startingPoint)
                    listePtLinear.append(secondPoint)
                self.currentLength +=d
#                self.completion = float(self.currentLength)/self.length
#                print "compl",self.completion
#                cent1T = startingPoint#self.transformPoints(jtrans, rotMatj, self.positions[-1])
#                cent2T = secondPoint#self.transformPoints(jtrans, rotMatj, self.positions2[-1])
#                if len(cent1T) == 1 :
#                    cent1T=cent1T[0]
#                if len(cent2T) == 1 :
#                    cent2T=cent2T[0]
#                print cent1T,cent2T
                if runTimeDisplay:
                    print ("cylinder with",cent1T,cent2T) 
                    name = str(len(listePtLinear))+self.name + str(ptInd)+"cyl"
                    if r :
                        name = str(len(listePtLinear)+1)+self.name + str(ptInd)+"cyl"                    
                    cyl=self.vi.oneCylinder(name,cent1T,cent2T,parent=parent,
                                            instance=histoVol.afviewer.becyl,
                                            radius=self.radii[0][0])
                                            
                    #self.vi.updateTubeMesh(cyl,cradius=self.radii[0][0])
                    self.vi.update()
                if r :
                    listePtCurve.insert(0,jtrans)
                else :   
                    listePtCurve.append(jtrans)
#        if success:
#            for jtrans,rotMatj in self.results:
                #every two we update distance from the previous
                if len(self.results) >= 1 :
                    #jtrans, rotMatj = self.results[-1]
#                    print "trasnfor",jtrans,rotMatj
                    #cent1T=self.transformPoints(jtrans, rotMatj, self.positions[-1])
                    insidePoints = {}
                    newDistPoints = {}
#                    rotMatj=[[ 1.,  0.,  0.,  0.],
#                       [ 0.,  1.,  0.,  0.],
#                       [ 0.,  0.,  1.,  0.],
#                       [ 0.,  0.,  0.,  1.]]
#                    jtrans = [0.,0.,0.]
#                    #move it or generate it inplace
#                    oldpos1=self.positions
#                    oldpos2=self.positions2
#                    if len(cent1T) == 1 :
#                        cent1T=cent1T[0]
#                    if len(cent2T) == 1 :
#                        cent2T=cent2T[0]                    
#                    self.positions=[[cent1T],]
#                    self.positions2=[[cent2T],]
                    #rbnode = histoVol.callFunction(histoVol.addRB,(self, numpy.array(jtrans), numpy.array(rotMatj),),{"rtype":self.Type},)#cylinder
                    #histoVol.callFunction(histoVol.moveRBnode,(rbnode, jtrans, rotMatj,))
#                    print ("before getInsidePoints ok")
                    insidePoints,newDistPoints = self.getInsidePoints(histoVol.grid,
                                        gridPointsCoords,dpad,distance,
                                        centT=cent1T,
                                        jtrans=jtrans, 
                                        rotMatj=rotMatj)
#                    print ("getInsidePoints ok",len(insidePoints))
#                    
#                    self.positions=oldpos1
#                    self.positions2=oldpos2
                    # update free points
                    #print "update free points",len(insidePoints)
                    nbFreePoints = self.updateDistances(histoVol,insidePoints, newDistPoints, 
                                freePoints, nbFreePoints, distance, 
                                gridPointsCoords, verbose)
#                    print ("updateDistances ok",nbFreePoints)
    #                print "distance",d
#                    print "free",nbFreePoints
            #Start Graham on 5/16/12 This progress bar doesn't work properly... compare with my version in HistoVol
                if histoVol.afviewer is not None and hasattr(histoVol.afviewer,"vi"):
                    histoVol.afviewer.vi.progressBar(progress=int((self.currentLength / self.length)*100),
                                                         label=self.name+str(self.currentLength / self.length)+" "+str(self.nbCurve)+"/"+str(self.nbMol))
                else :
                    autopack.helper.progressBar(progress=int((self.currentLength / self.length)*100),
                                                         label=self.name+str(self.currentLength / self.length)+" "+str(self.nbCurve)+"/"+str(self.nbMol))
  
            #Start Graham on 5/16/12 This progress bar doesn't work properly... compare with my version in HistoVol
                if self.currentLength >= self.length:
                    Done = True
                    self.counter = counter +1
#                    return success, nbFreePoints
#                print ("end while loop ok",self.currentLength)
            else :
                secondPoint = startingPoint
                break
        return success, nbFreePoints,freePoints

    def updateGrid(self,rg,histoVol,dpad,freePoints, nbFreePoints, distance, 
                        gridPointsCoords, verbose):
        for i in range(rg):#len(self.results)):
            jtrans, rotMatj = self.results[-i]
            cent1T=self.transformPoints(jtrans, rotMatj, self.positions[-1])
            insidePoints = {}
            newDistPoints = {}
            insidePoints,newDistPoints = self.getInsidePoints(histoVol.grid,
                                gridPointsCoords,dpad,distance,
                                centT=cent1T,
                                jtrans=jtrans, 
                                rotMatj=rotMatj)
            # update free points
            nbFreePoints = self.updateDistances(histoVol,insidePoints, newDistPoints, 
                        freePoints, nbFreePoints, distance, 
                        gridPointsCoords, verbose)  
            return nbFreePoints,freePoints

    def getFirstPoint(self,ptInd,seed=0):
        if self.compNum>0 : #surfacegrowing: first point is aling to the normal:
            v2 = self.histoVol.compartments[abs(self.compNum)-1].surfacePointsNormals[ptInd]
            secondPoint = numpy.array(self.startingpoint)+numpy.array(v2)*self.uLength
        else :               
            #randomize the orientation in the hemisphere following the direction.
            v = self.vi.rotate_about_axis(numpy.array(self.orientation),
                                          random()*math.radians(self.marge),#or marge ?
                                          axis=list(self.orientation).index(0))
            self.vector = numpy.array(v).flatten()*self.uLength*self.jitterMax# = (1,0,0)self.vector.flatten()
            secondPoint = self.startingpoint+self.vector
            print ("newPoints",self.startingpoint,secondPoint)
            #seed="F"
            if seed :
                seed="R"
                secondPoint = self.startingpoint-self.vector
            else :
                seed="F"
            inside = self.histoVol.grid.checkPointInside(secondPoint,dist=self.cutoff_boundary,jitter=self.jitterMax)
            closeS=False
            if inside and self.compNum<=0: 
                #only if not surface ingredient
                closeS = self.checkPointSurface(secondPoint,cutoff=self.cutoff_surface)
            success = False
            if not inside or closeS:
                safety = 30
                k=0
                while not inside or closeS:
                    if k > safety :
                        print ("cant find first inside point",inside,closeS)
                        return None
                    else :
                        k+=1
                    p = self.vi.advance_randpoint_onsphere(self.uLength,
                                                           marge=math.radians(self.marge),
                                                            vector=self.vector)
                    print ("p is ",p,k,safety)
                    print (self.uLength,self.marge,self.vector)
                    pt = numpy.array(p)*numpy.array(self.jitterMax)
                    secondPoint = self.startingpoint+numpy.array(pt)
                    print ("secdpoint ",secondPoint)
                    inside = self.histoVol.grid.checkPointInside(secondPoint,dist=self.cutoff_boundary,jitter=self.jitterMax)
                    if self.compNum<=0 : closeS = self.checkPointSurface(secondPoint,cutoff=self.cutoff_surface)
                    print (inside,closeS)
            if self.runTimeDisplay:
                parent = self.histoVol.afviewer.orgaToMasterGeom[self]
                name = "Startsp"+self.name+seed
                #sp=self.vi.Sphere(name,radius=self.radii[0][0],parent=parent)[0]
                if seed=="F" :
                    #sp=self.vi.newInstance(name,self.histoVol.afviewer.pesph,
                    #                   location=self.startingpoint,parent=parent)
                    #self.vi.scaleObj(sp,self.radii[0][0])
                    sp=self.vi.Sphere(name,radius=self.radii[0][0],parent=parent)[0]
                    self.vi.setTranslation(sp,pos=self.startingpoint)
    #            sp.SetAbsPos(self.vi.FromVec(startingPoint))
                cyl = self.vi.oneCylinder(name+"cyl",self.startingpoint,secondPoint,
                                          instance=self.histoVol.afviewer.becyl,
                                          parent=parent,radius=self.radii[0][0])
    #            self.vi.updateTubeMesh(cyl,cradius=self.radii[0][0])
                #laenge,mx=self.getTubeProperties(head,tail)
                self.vi.update()        
        return secondPoint

    #isit the jitter place ? I guess  Why are there two jitter_place functions?  What is this one?
    def jitter_place(self, histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,
              stepByStep=False, verbose=False,
              sphGeom=None, labDistGeom=None, debugFunc=None,
              sphCenters=None,  sphRadii=None, sphColors=None,usePP=False):
        #ptInd is the starting point
        #loop over until length reach or cant grow anymore
#        self.nbMol = 1               
#        print("JitterPlace ingr Grow....................................................................")
        if type(self.compMask) is str :
            self.compMask = eval(self.compMask) 
            self.prepare_alternates()
        success = True
        self.vi = autopack.helper
#        if histoVol.afviewer != None :
#            self.vi = histoVol.afviewer.vi
        afvi = histoVol.afviewer
        self.histoVol=histoVol
        windowsSize = histoVol.windowsSize
        simulationTimes = histoVol.simulationTimes
        gridPointsCoords = histoVol.grid.masterGridPositions
        self.runTimeDisplay=runTimeDisplay = histoVol.runTimeDisplay
        normal = None
        #self.startingpoint = previousPoint = startingPoint = numpy.array(histoVol.grid.masterGridPositions[ptInd])
        #jitter the first point
        if self.compNum > 0 :
            normal = histoVol.compartments[abs(self.compNum)-1].surfacePointsNormals[ptInd]
        self.startingpoint = previousPoint = startingPoint =self.jitterPosition(numpy.array(histoVol.grid.masterGridPositions[ptInd]), histoVol.smallestProteinSize, normal = normal)

#        print ("PTID ", numpy.array(histoVol.grid.masterGridPositions[ptInd]), "first ",startingPoint)
        v,u = self.vi.measure_distance(self.positions,self.positions2,vec=True)
#        self.uLength = u
        self.vector = numpy.array(self.orientation) * self.uLength
        
        if u != self.uLength :
            self.positions2 = [[self.vector]]      
        if self.compNum == 0 :
            compartment = self.histoVol
        else :
            compartment = self.histoVol.compartments[abs(self.compNum)-1]
        self.compartment = compartment
#        print "get",self.name, compartment.name,len(compartment.molecules)
        secondPoint = self.getFirstPoint(ptInd)
        #check collision ?
        #if we have starting position available use it
        #print ("start pos?",len(self.start_positions),self.nbCurve)        
        if self.nbCurve < len(self.start_positions) :
            self.startingpoint = previousPoint = startingPoint = self.start_positions[self.nbCurve][0]
            secondPoint = self.start_positions[self.nbCurve][1]
#            print ("use Starting pos",self.start_positions[self.nbCurve])
        if secondPoint is None :
#            self.completion = float(self.nbCurve)/self.nbMol
#            print ("no second point ", self.completion,success)
            return success, nbFreePoints
        rotMatj,jtrans=self.getJtransRot(startingPoint,secondPoint)
        #test for collision
        #return success, nbFreePoints
        self.results.append([jtrans,rotMatj])
        if self.placeType =="pandaBullet":
           self.histoVol.nb_ingredient+=1
           self.histoVol.rTrans.append(numpy.array(startingPoint).flatten())
           self.histoVol.rRot.append(numpy.array(numpy.identity(4)))#rotMatj 
           self.histoVol.rIngr.append(self)
           self.histoVol.result.append([ [numpy.array(startingPoint).flatten(),secondPoint], rotMatj, self, 0 ])
#           if len(self.histoVol.rTrans) > 1 : bhtreelib.freeBHtree(self.histoVol.close_ingr_bhtree)
#           if len(self.histoVol.rTrans) : self.histoVol.close_ingr_bhtree=bhtreelib.BHtree( self.histoVol.rTrans, None, 10)
        elif self.placeType =="RAPID":
           self.histoVol.nb_ingredient+=1
           self.histoVol.rTrans.append(numpy.array(jtrans).flatten())
           self.histoVol.rRot.append(numpy.array(rotMatj))#rotMatj 
           self.histoVol.rIngr.append(self)
           self.histoVol.result.append([ [numpy.array(startingPoint).flatten(),secondPoint], rotMatj, self, 0 ])

        if self.histoVol.treemode == "bhtree":# "cKDTree"
            if len(self.histoVol.rTrans) > 1 : bhtreelib.freeBHtree(self.histoVol.close_ingr_bhtree)
            if len(self.histoVol.rTrans) : self.histoVol.close_ingr_bhtree=bhtreelib.BHtree( self.histoVol.rTrans, None, 10)
        else :
            #rebuild kdtree
            if len(self.histoVol.rTrans) > 1 :self.histoVol.close_ingr_bhtree= spatial.cKDTree(self.histoVol.rTrans, leafsize=10)

        self.currentLength = 0.
        counter = self.counter
#        self.Ptis=[ptInd,histoVol.grid.getPointFrom3D(secondPoint)]
        dist,pid = histoVol.grid.getClosestGridPoint(secondPoint)
        self.Ptis=[ptInd,pid]
        listePtCurve=[jtrans]
        listePtLinear=[startingPoint,secondPoint]
        Done = False
        k=0
        v=[0.,0.,0.]
        d=0.
#        print ("so far :",listePtLinear)
        #grow until reach self.currentLength >= self.length
        #or attempt > safety
        success, nbFreePoints,freePoints = self.grow(previousPoint,startingPoint,secondPoint,
                                          listePtCurve,listePtLinear,histoVol, 
                                          ptInd, freePoints, nbFreePoints, distance, 
                                          dpad,stepByStep=False, verbose=False,usePP=usePP)
        nbFreePoints,freePoints = self.updateGrid(2,histoVol,dpad,freePoints, nbFreePoints, distance, 
                        gridPointsCoords, verbose)
        if self.seedOnMinus :
#            print "reverse"
            #secondPoint = self.getFirstPoint(ptInd,seed=1)
            success, nbFreePoints,freePoints = self.grow(previousPoint,listePtLinear[1],
                                          listePtLinear[0],
                                          listePtCurve,listePtLinear,histoVol, 
                                          ptInd, freePoints, nbFreePoints, distance, 
                                          dpad,stepByStep=False, verbose=False,r=True)
            nbFreePoints,freePoints = self.updateGrid(2,histoVol,dpad,freePoints, nbFreePoints, distance, 
                               gridPointsCoords, verbose)
        #store result in molecule
        if verbose:
            print("res",len(self.results))
        for i in range(len(self.results)):
            jtrans, rotMatj = self.results[-i]
            dist,ptInd = histoVol.grid.getClosestGridPoint(jtrans)
            compartment.molecules.append([ jtrans, rotMatj, self, ptInd ])
            #reset the result ?
        self.results=[]
#        print ("After :",listePtLinear)
        self.listePtCurve.append(listePtCurve)
        self.listePtLinear.append(listePtLinear)
        self.nbCurve+=1
        self.completion = float(self.nbCurve)/float(self.nbMol)
        if verbose:
            print("completion",self.completion,self.nbCurve,self.nbMol)
        return success, nbFreePoints            
        
    def prepare_alternates(self,):
        if len(self.partners) :
            self.alternates_names = self.partners.keys()#self.partners_name#[p.name for p in self.partners.values()]
            #self.alternates_weight = [self.partners[name].weight for name in self.partners]
            self.alternates_weight = [self.partners[name].ingr.weight for name in self.partners]
            self.alternates_proba = [self.partners[name].ingr.proba_binding for name in self.partners]

    def prepare_alternates_proba(self,):
        thw=[]
        tw=0.0
        weights = self.alternates_proba#python3?#dict.copy().keys()
        for i, w in enumerate(weights):
            tw+=w
            thw.append(tw) 
        self.alternates_proba = thw

    def pick_random_alternate(self,):
        if not len(self.alternates_names):
            return None,0
        r =  uniform(0,1.0)
        ar =  uniform(0,1.0)
        #weights = self.alternates_weight[:]
        proba = self.alternates_proba[:]
        alti = int((r*(len(self.alternates_names))))#round?
        #print (alti,proba[alti],ar,self.alternates_names[alti])
        if ar < proba[alti]:
            return self.alternates_names[alti],alti
        return None,0
        
    def pick_alternate(self,):
        #whats he current length ie number of point so far
        #whats are he number of alternate and theyre proba
        #pick an alternate according length and proba
        #liste_proba
        #liste_alternate
        #dice = uniform(0.0,1.0)
        #int(uniform(0.0,1.0)*len(self.sphere_points_mask))
        alt_name = None
        #if it is the first two segment dont do it
        if self.currentLength <= self.uLength * 2.0 :
            return None,0
#        return self.pick_random_alternate()
        weights = self.alternates_weight#python3?#dict.copy().keys()
        rnd = uniform(0,1.0) * sum(weights)# * (self.currentLength / self.length)
#        print ("alter",weights,rnd) 
        i=0        
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                r =  uniform(0,1.0)
#                print (r,self.alternates_proba[i])
                if r < self.alternates_proba[i]:
                    alt_name = self.alternates_names[i]
                break
        #alternates_names point to an ingredients id?
        return alt_name,i



    def get_alternate_position(self,alternate,alti,v,pt1,pt2):
        length = self.partners[alternate].getProperties("length")
        #rotation that align snake orientation to current segment        
        rotMatj,jtrans=self.getJtransRot_r(numpy.array(pt1).flatten(),
                                           numpy.array(pt2).flatten(),
                                           length = length)
        #jtrans is the position between pt1 and pt2
        prevMat = numpy.array(rotMatj)
        #jtrans=autopack.helper.ApplyMatrix([jtrans],prevMat.transpose())[0]
        prevMat[3,:3] = numpy.array(pt1)#jtrans   
        rotMatj = numpy.identity(4)
        #oldv is v we can ether align to v or newv
        newv = numpy.array(pt2) - numpy.array(pt1)
        #use v ? for additional point ?
        pta=self.partners[alternate].getProperties("pt1")        
        ptb=self.partners[alternate].getProperties("pt2")
        ptc=self.partners[alternate].getProperties("pt3")
        toalign = numpy.array(ptc) -numpy.array(ptb)
        m = numpy.array(rotVectToVect(toalign,newv)).transpose()
        m[3,:3] = numpy.array(pt1)#jtrans
        pts=autopack.helper.ApplyMatrix([ptb],m.transpose())#transpose ?
        v = numpy.array(pt1)-pts[0]
        m[3,:3] = numpy.array(pt1)+v# - (newpt1-pts[0])
        
        #rotMatj,jt=self.getJtransRot_r(numpy.array(ptb).flatten(),
        #                                   numpy.array(ptc).flatten(),
        #                                   length = length)
        #rotMatj[3,:3] = -numpy.array(ptb)
        #globalM1 = numpy.array(matrix(rotMatj)*matrix(prevMat))
        #
        #offset = numpy.array(ptb)+toalign/2.0
        #npts=numpy.array([pta,ptb,offset,ptc])-numpy.array([ptb])
        #pts=autopack.helper.ApplyMatrix(npts,globalM1.transpose())#transpose ?        
        #trans = numpy.array(jtrans)-1.5*pts[1]-pts[2]
        #now apply matrix and get the offset
        #prevMat = numpy.array(globalM1)
        #jtrans=autopack.helper.ApplyMatrix([jtrans],prevMat.transpose())[0]
        #prevMat[3,:3] = jtrans   
        #npt2=autopack.helper.ApplyMatrix([ptb],prevMat.transpose())[0]
        #offset = numpy.array(npt2) -numpy.array(pt1)
        #jtrans=numpy.array(jtrans)-offset
        #toalign = numpy.array(ptb) -numpy.array(pta)
        #globalM2 = numpy.array(rotVectToVect(toalign,v))
        #compare to superimposition_matrix
        #print globalM1,quaternion_from_matrix(globalM1).tolist()
        #print globalM2,quaternion_from_matrix(globalM2).tolist()
        #center them
        #c1=autopack.helper.getCenter([ptb,ptc])
        #c2=autopack.helper.getCenter([pt1,pt2])
        #globalM = superimposition_matrix(numpy.array([ptb,ptc])-c1,numpy.array([pt1,pt2])-c2)
        #print globalM,quaternion_from_matrix(globalM).tolist()
        rotMatj = numpy.identity(4)
        rotMatj[:3,:3] = m[:3,:3].transpose()
        jtrans = m[3,:3]
        #print ("will try to place alterate at ",jtrans) 
        return jtrans,rotMatj
        
    def get_alternate_position_p(self,alternate,alti,v,pt1,pt2):
        length = self.partners[alternate].getProperties("length")
        rotMatj,jtrans=self.getJtransRot_r(numpy.array(pt1).flatten(),
                                           numpy.array(pt2).flatten(),
                                           length = length)
        prevMat = numpy.array(rotMatj)
        #jtrans=autopack.helper.ApplyMatrix([jtrans],prevMat.transpose())[0]
        prevMat[3,:3] = jtrans   
        rotMatj = numpy.identity(4)
        
        localMR = self.partners_position[alti]
        #instead use rotVectToVect from current -> to local ->
        #align  p2->p3 vector to pt1->pt2
        globalM = numpy.array(matrix(localMR)*matrix(prevMat))
        jtrans = globalM[3,:3]
        rotMatj[:3,:3] = globalM[:3,:3]
        #print ("will try to place alterate at ",jtrans) 
        return jtrans,rotMatj

    def get_alternate_starting_point(self,pt1,pt2,alternate):
        spt1 = self.partners[alternate].getProperties("st_pt1")
        spt2 = self.partners[alternate].getProperties("st_pt2") 

        length = self.partners[alternate].getProperties("length")
        rotMatj,jtrans=self.getJtransRot_r(numpy.array(pt1).flatten(),
                                           numpy.array(pt2).flatten(),
                                           length = length)
        #jtrans is the position between pt1 and pt2
        prevMat = numpy.array(rotMatj)
        prevMat[3,:3] = numpy.array(pt1)#jtrans   
        newv = numpy.array(pt2) - numpy.array(pt1)
        ptb=self.partners[alternate].getProperties("pt2")
        ptc=self.partners[alternate].getProperties("pt3")
        toalign = numpy.array(ptc) -numpy.array(ptb)
        m = numpy.array(rotVectToVect(toalign,newv)).transpose()
        m[3,:3] = numpy.array(pt1)#jtrans
        pts=autopack.helper.ApplyMatrix([ptb],m.transpose())#transpose ?
        v = numpy.array(pt1)-pts[0]
        m[3,:3] = numpy.array(pt1)+v
        newPts=autopack.helper.ApplyMatrix([spt1,spt2],m.transpose())#transpose ?
        return newPts

    def place_alternate(self,alternate,alti,v,pt1,pt2):
        pta=self.partners[alternate].getProperties("pt1")        
        ptb=self.partners[alternate].getProperties("pt2")
        ptc=self.partners[alternate].getProperties("pt3")
        ptd=self.partners[alternate].getProperties("pt4")
        prevMat = numpy.identity(4)        
        if ptb is not None :
            rotMatj,jtrans=self.getJtransRot_r(numpy.array(pt1).flatten(),
                                               numpy.array(pt2).flatten())
            toalign = numpy.array(ptb) -numpy.array(pta)
            newv = numpy.array(pt2) - numpy.array(pt1)        
            prevMat = numpy.array(rotVectToVect(toalign,newv))
            newPts = autopack.helper.ApplyMatrix([ptc,ptd],prevMat.transpose()) #partner positions ?
            prevMat[3,:3] = jtrans   
        else :
            newPt = self.pickHalton(pt1,pt2)
            newPts=[newPt]
            rotMatj,jtrans=self.getJtransRot_r(numpy.array(pt2).flatten(),
                                               numpy.array(newPt).flatten())
            prevMat = numpy.array(rotMatj)
            #jtrans=autopack.helper.ApplyMatrix([jtrans],prevMat.transpose())[0]
            prevMat[3,:3] = jtrans   
        rotMatj= numpy.identity(4)
        return newPts,jtrans,rotMatj
        
    def place_alternate_p(self,alternate,alti,v,pt1,pt2):
        #previou transformation
#        distance,mat = autopack.helper.getTubePropertiesMatrix(pt1,pt2)
#        prevMat = numpy.array(mat)
        #rotMatj,jtrans=self.getJtransRot(numpy.array(pt2).flatten(),numpy.array(pt1).flatten())
        #should apply first to get the new list of point, and length
        p_alternate = self.partners[alternate]#self.histoVol.getIngrFromNameInRecipe(alternate,self.recipe )
        #if p_alternate.getProperties("bend"):
        out1=p_alternate.getProperties("pt1")
        out2=p_alternate.getProperties("pt2")   
        if out1 is not None :         
            rotMatj,jtrans=self.getJtransRot_r(numpy.array(pt1).flatten(),
                                               numpy.array(pt2).flatten())
            prevMat = numpy.array(rotMatj)
            #jtrans=autopack.helper.ApplyMatrix([jtrans],prevMat.transpose())[0]
            prevMat[3,:3] = jtrans   
            newPts = autopack.helper.ApplyMatrix([out1,out2],prevMat.transpose()) #partner positions ?
        else :
            newPt = self.pickHalton(pt1,pt2)
            newPts=[newPt]
            rotMatj,jtrans=self.getJtransRot_r(numpy.array(pt2).flatten(),
                                               numpy.array(newPt).flatten())
            prevMat = numpy.array(rotMatj)
            #jtrans=autopack.helper.ApplyMatrix([jtrans],prevMat.transpose())[0]
            prevMat[3,:3] = jtrans   
        rotMatj = numpy.identity(4)
        #print ("som math",out1,out2,newPts)
        #need also to get the alternate_ingredint new position and add it.
        localMR = self.partners_position[alti]
        globalM = numpy.array(matrix(localMR)*matrix(prevMat))
        jtrans = globalM[3,:3]
        rotMatj[:3,:3] = globalM[:3,:3]
#        print ("will try to place alterate at ",jtrans) 
        return newPts,jtrans,rotMatj
        #we need to add this guy a new mol and tak in account his molarity ?
        #should actually the partner system and the step/step
#        return newPts,None,None
        
class ActinIngrediant(GrowIngrediant):
    def __init__(self, molarity, radii=[[50.],], positions=None, positions2=None,
                 sphereFile=None, packingPriority=0, name=None,
                 pdb=None, color=None, nbJitter=5, jitterMax=(1,1,1),
                 perturbAxisAmplitude = 0.1,length = 10.,closed = False,
                 modelType="Cylinders",biased=1.0,Type="Actine",
                 principalVector=(1,0,0), meshFile=None, packingMode='random',
                 placeType="jitter",marge=35.0,influenceRad =100.0,meshObject=None,
                 orientation = (1,0,0),nbMol=0,**kw):
#        Ingredient.__init__(self, molarity, radii, positions, positions2,
#                            sphereFile,
#                            packingPriority, name, pdb, color, nbJitter,
#                            jitterMax, perturbAxisAmplitude, principalVector,
#                            meshFile, packingMode,placeType,meshObject)

        GrowIngrediant.__init__(self, molarity, radii, positions, positions2,
                            sphereFile,
                            packingPriority, name, pdb, color, nbJitter,
                            jitterMax, perturbAxisAmplitude, 
                            length,closed,
                            modelType,biased,
                            principalVector,
                            meshFile, packingMode,placeType,marge,meshObject,
                            orientation ,nbMol,Type,**kw)
        if name == None:
            name = "Actine_%s_%f"% (str(radii),molarity)
        self.isAttractor = True
        self.constraintMarge = True
        self.seedOnMinus = True
        self.influenceRad = influenceRad
        self.oneSuperTurn = 825.545#cm from c4d graham file
        self.oneDimerSize = 100.0#200 =2 
        self.cutoff_surface = 50.
        self.cutoff_boundary = 1.0
        #how to setup the display from the 3dmesh?should I update radius, position from it ?
        #also the mesh should preoriented correctly...
#        self.positions=[[[0,0,0]]], 
#        self.positions2=[[[0.,100.,0]]],
        

#    def build
        
    def updateFromBB(self,grid):
        return
        r = grid.getRadius()
        self.positions = [0.,0.,0.]
        self.positions2 = [r,0.,0.]
        self.principalVector = [1.,0.,0.]
        self.uLength = r
        self.length = 2*r

#import autopack
#from autopack.Recipe import Recipe
import os
#from DejaVu  import colors
class IngredientDictionary:
    """
    list all available ingrediant, and permit to add user ingrediant
    the listing is based on a scanning of the autopack package directory
    under recipe/membrane and recipe/cyto
    """
    def __init__(self):
        """Set up from the default directory"""
        #get the directory and scan all the name
        self.rdir = autopack.__path__[0]+'/recipes/'
        self.update()
        self.knownIngr = {}
        self.knownIngr['membrane'] = {
        '1h6i':{'name':'AUQAPORINE','priority':-10,'pVector':(0,0,1),'jitterMax':(0.3,0.3,0.3),'mol':.04,'type':'singleS','mode':'random'},
        '1zll':{'name':'PHOSPHOLAMBAN','priority':-5,'pVector':(0,0,-1),'jitterMax':(1,1,1),'mol':0.04,'type':'singleS','mode':'random'},
        '2afl':{'name':'PROTON TRANSPORT','priority':-10,'pVector':(0,0,-1),'jitterMax':(1,1,1),'mol':0.04,'type':'singleS','mode':'random'},
        '2uuh':{'name':'C4 SYNTHASE','priority':-5,'pVector':(0,0,-1),'jitterMax':(1,1,1),'mol':0.04,'type':'multiS','mode':'random'},
        '1yg1':{'name':'FACILITATIVE GLUCOSE','priority':-5,'pVector':(0,0,-1),'jitterMax':(1,1,1),'mol':.04,'type':'multiS','mode':'random'},
        '1ojc':{'name':'OXIDOREDUCTASE','priority':-10,'pVector':(0,0,-1),'jitterMax':(1,1,1),'mol':0.04,'type':'multiS','mode':'random'},
        'ves34':{'name':'VES34','priority':-100,'pVector':(0,0,-1),'jitterMax':(1,1,1),'mol':0.04,'type':'multiS','mode':'random'},
        '1qo1':{'name':'ATP SYNTHASE','priority':-150,'pVector':(0,0,-1),'jitterMax':(1,1,1),'mol':0.01,'type':'multiS','mode':'random'},
        '2abm':{'name':'AQUAPORIN TETRAMER','priority':-2,'pVector':(0,0,-1),'jitterMax':(1,1,1),'mol':0.04,'type':'multiS','mode':'random'},
        '3g61':{'name':'P-GLYCOPROTEIN','priority':-2,'pVector':(0,0,1),'jitterMax':(1,1,1),'mol':0.04,'type':'multiS','mode':'random'},
        '2bg9':{'name':'ION CHANNEL/RECEPTOR','priority':-2,'pVector':(0,0,-1),'jitterMax':(1,1,1),'mol':0.04,'type':'multiS','mode':'random'},
        '2a79':{'name':'POTASSIUM CHANNEL','priority':-1,'pVector':(0,0,1),'jitterMax':(0.3,0.3,0.3),'mol':0.04,'type':'multiS','mode':'random'}
        }
        self.knownIngr['cyto'] = {
        '1AON_centered':{'name':'GROEL','priority':1,'pVector':(1,0,0),'jitterMax':(1,1,1),'mol':.0004,'type':'multiS','mode':'random'},
        '2CPK_centered':{'name':'PHOSPHOTRANSFERASE','priority':1,'pVector':(1,0,0),'jitterMax':(1,1,1),'mol':.0004,'type':'multiS','mode':'close'},
        '1CZA_centered':{'name':'TRANSFERASE','priority':1,'pVector':(1,0,0),'jitterMax':(1,1,1),'mol':.0002,'type':'multiS','mode':'random'},
        '2OT8_centered':{'name':'TRANSPORT PROTEIN','priority':1,'pVector':(1,0,0),'jitterMax':(1,1,1),'mol':.0002,'type':'multiS','mode':'random'},
        '1TWT_centered':{'name':'30S RIBOSOME','priority':2,'pVector':(1,0,0),'jitterMax':(1,1,1),'mol':.0001,'type':'multiS','mode':'random'},
        '1ABL_centered':{'name':'KINASE','priority':0,'pVector':(1,0,0),'jitterMax':(1,1,1),'mol':.01,'type':'singleS','rad':16.,'mode':'random'}
        }
        self.knownIngr['matrix'] = {
        #'1ABL_centered':{'name':'GLUTAMATE','priority':.0001,'pVector':(1,0,0),'jitterMax':(1,1,1),'mol':.150,'type':'singleS','rad':3.61,'mode':'random'}
        'Glutamate_centered':{'name':'GLUTAMATE','priority':.0001,'pVector':(1,0,0),'jitterMax':(1,1,1),'mol':.150,'type':'singleS','rad':3.61,'mode':'random'}
        }
        self.MSca = 1.0
    
    def filterList(self,dir):
        #if no .sph file -> single Sphere ingr
        #else ->multiSphIngr
        listDir=os.listdir(dir)
        listeIngr={}
        for item in listDir:
            if item == 'CVS' :
                continue
            sp = item.split('.')
            if sp[0] not in listeIngr :
                if sp[0]+'.sph' in listDir :
                    nbSph = self.getNumberOfSphere(dir+'/'+sp[0]+'.sph')
                    if nbSph > 1 :
                        listeIngr[sp[0]]=['Multi',True]
                    else : 
                        listeIngr[sp[0]]=['Single',True]
                else :
                    listeIngr[sp[0]]=['Single',False]
        return listeIngr
    
    def changeDir(self,newDir):
        self.rdir = newDir
        
    def update(self,newDir=None):
        if newDir is not None :
            self.changeDir(newDir)
        self.listename={}
        self.listename["membrane"] = self.filterList(self.rdir+'membrane')
        self.listename["cyto"] = self.filterList(self.rdir+'cyto')

    def getNumberOfSphere(self,sphFile):
        f = open(sphFile,'r')
        lines = f.readlines()
        f.close()
        for i in range(len(lines)):
            if lines[i] == '# number of spheres in level 1\n':
                nbSphere = eval(lines[i+1])
                return nbSphere

    def makeBiLayerIngrediant(self,molarity,positions=[-.5,0,0], 
                            positions2=[15,0,0],radii=[3],priority = .001,
                            jitterMax=(0.3,0.3,0.3)):
        ############## cylinder Bilayer here##########
        #  jitterMax CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
        cyl1IngrU = MultiCylindersIngr(molarity,  pdb='LipOutPdb',
                                        name='LipOut', radii=[radii],
                                        positions=[[positions]], positions2=[[positions2]],
                                        packingPriority=priority,
                                        jitterMax=jitterMax, 
                                        principalVector=(1,0,0)#vcolor=colors.bisque,
                                      )        
        cyl1IngrD = MultiCylindersIngr(molarity,  pdb='LipInPdb',
                                        name='LipIn', radii=[radii],
                                        positions=[[positions]], positions2=[[positions2]],
                                        packingPriority=priority,
                                        jitterMax=jitterMax,
                                        principalVector=(-1,0,0)
                                      )#color=colors.burlywood,
        return cyl1IngrU,cyl1IngrD

    def makeKnownIngrediant(self,MSca,listeCol):
        rSurf1 = Recipe()
        rCyto = Recipe()
        rMatrix = Recipe()
        rRec = {'membrane':rSurf1,
                'cyto':rCyto,
                'matrix':rMatrix}
        wrkDir = self.rdir
        i=0
        for k in list(self.knownIngr.keys()):
            for ing_name in list(self.knownIngr[k].keys()):
                ingDic = self.knownIngr[k][ing_name]
                if ingDic['type'] == 'singleS':
                    if 'rad' in ingDic:
                         ingr = SingleSphereIngr( MSca*ingDic['mol'], color=listeCol[i], pdb=ing_name,
                                         name=ingDic['name'],
                                         radius=ingDic['rad'],
                                         meshFile=wrkDir+'/'+k+'/'+ing_name,
                                         packingPriority=ingDic['priority'],
                                         jitterMax=ingDic['jitterMax'],
                                         principalVector=ingDic['pVector'],
                                         packingMode=ingDic['mode'])                       
                    else :
                        ingr = SingleSphereIngr( MSca*ingDic['mol'], color=listeCol[i], pdb=ing_name,
                                         name=ingDic['name'],
                                         sphereFile=wrkDir+'/'+k+'/'+ing_name+'.sph',
                                         meshFile=wrkDir+'/'+k+'/'+ing_name,
                                         packingPriority=ingDic['priority'],
                                         jitterMax=ingDic['jitterMax'],
                                         principalVector=ingDic['pVector'],
                                         packingMode=ingDic['mode'])
                elif ingDic['type'] == 'multiS':
                    ingr = MultiSphereIngr( MSca*ingDic['mol'], color=listeCol[i], pdb=ing_name,
                                         name=ingDic['name'],
                                         sphereFile=wrkDir+'/'+k+'/'+ing_name+'.sph',
                                         meshFile=wrkDir+'/'+k+'/'+ing_name,
                                         packingPriority=ingDic['priority'],
                                         jitterMax=ingDic['jitterMax'],
                                         principalVector=ingDic['pVector'],
                                         packingMode=ingDic['mode'])
                rRec[k].addIngredient(ingr)
                i = i +1
        return rRec


    def getIngrediants(self,name,color):
        wrkDir = self.rdir
        ingr=None
        for k in list(self.knownIngr.keys()):
            for ing_name in list(self.knownIngr[k].keys()):
                ingDic = self.knownIngr[k][ing_name]
                if name == ing_name or name == ingDic['name']:
                    if ingDic['type'] == 'singleS':
                        if 'rad' in ingDic:
                             ingr = SingleSphereIngr( self.MSca*ingDic['mol'], color=color, pdb=ing_name,
                                             name=ingDic['name'],
                                             radius=ingDic['rad'],
                                             meshFile=wrkDir+'/'+k+'/'+ing_name,
                                             packingPriority=ingDic['priority'],
                                             jitterMax=ingDic['jitterMax'],
                                             principalVector=ingDic['pVector'],
                                             packingMode=ingDic['mode'])                       
                        else :
                            ingr = SingleSphereIngr( self.MSca*ingDic['mol'], color=color, pdb=ing_name,
                                             name=ingDic['name'],
                                             sphereFile=wrkDir+'/'+k+'/'+ing_name+'.sph',
                                             meshFile=wrkDir+'/'+k+'/'+ing_name,
                                             packingPriority=ingDic['priority'],
                                             jitterMax=ingDic['jitterMax'],
                                             principalVector=ingDic['pVector'],
                                             packingMode=ingDic['mode'])
                    elif ingDic['type'] == 'multiS':
                        ingr = MultiSphereIngr( self.MSca*ingDic['mol'], color=color, pdb=ing_name,
                                             name=ingDic['name'],
                                             sphereFile=wrkDir+'/'+k+'/'+ing_name+'.sph',
                                             meshFile=wrkDir+'/'+k+'/'+ing_name,
                                             packingPriority=ingDic['priority'],
                                             jitterMax=ingDic['jitterMax'],
                                             principalVector=ingDic['pVector'],
                                             packingMode=ingDic['mode'])
                    break
        return ingr

    def makeIngrediants(self,MSca,listeCol,wanted):
        rSurf1 = Recipe()
        rCyto = Recipe()
        rRec = {'membrane':rSurf1,
                'cyto':rCyto}
        wrkDir = self.rdir
        #excluded = ['ves34']
        i=0
        for k in list(self.listename.keys()):
            for ing_name in list(self.listename[k].keys()):
                #print ing_name
                if ing_name in wanted:
                    if self.listename[k][ing_name][0] == 'Single':
                        sphereFile = None
                        if self.listename[k][ing_name][1] :
                            sphereFile = wrkDir+'/'+k+'/'+ing_name+'.sph'
                        ingr = SingleSphereIngr( MSca*.04, color=listeCol[i], pdb=ing_name,
                                         name=ing_name,
                                         sphereFile=sphereFile,
                                         meshFile=wrkDir+'/'+k+'/'+ing_name,
                                         packingPriority=-2,
                                         jitterMax=(0.3,0.3,0.3),
                                         principalVector=(0,0,1))
                    else :
                        ingr = MultiSphereIngr( MSca*.04, color=listeCol[i], name=ing_name,
                                        sphereFile=wrkDir+'/'+k+'/'+ing_name+'.sph',
                                        meshFile=wrkDir+'/'+k+'/'+ing_name, pdb=ing_name,
                                        packingPriority=-1,
                                        #jitterMax=(1,1,.2),
                                        principalVector=(0,0,-1))
                    rRec[k].addIngredient(ingr)
                    i = i +1
        return rRec

from autopack import IOutils as io 
from xml.dom.minidom import getDOMImplementation 
import json   
class IOingredientTool:
    #xml parser that can return an ingredient
    def __init__(self):
        pass

    def read(self,filename):
        pass

    def write(self,ingr,filename,ingr_format="xml"):
        if ingr_format == "json" :
            ingdic = self.ingrJsonNode(ingr)
            with open(filename+".json", 'w') as fp :#doesnt work with symbol link ?
                json.dump(ingdic,fp,indent=4, separators=(',', ': '))#,indent=4, separators=(',', ': ')            
        elif ingr_format == "xml" :
            ingrnode,xmldoc = self.ingrXmlNode(ingr)
            f = open(filename+".xml","w")        
            xmldoc.writexml(f, indent="\t", addindent="", newl="\n")
            f.close()
        elif ingr_format == "python" :
            ingrnode = self.ingrPythonNode(ingr)
            f = open(filename+".py","w")        
            f.write(ingrnode)
            f.close()            
        elif ingr_format == "all" :
            ingdic = self.ingrJsonNode(ingr)
            with open(filename+".json", 'w') as fp :#doesnt work with symbol link ?
                json.dump(ingdic,fp,indent=4, separators=(',', ': '))#,indent=4, separators=(',', ': ')            
            ingrnode,xmldoc = self.ingrXmlNode(ingr)
            f = open(filename+".xml","w")        
            xmldoc.writexml(f, indent="\t", addindent="", newl="\n")
            f.close()
            ingrnode = self.ingrPythonNode(ingr)
            f = open(filename+".py","w")        
            f.write(ingrnode)
            f.close()

    def makeIngredientFromXml(self,inode=None,filename=None, recipe="Generic"):
        if filename is None and inode is not None :
            f=str(inode.getAttribute("include"))
            if f != '':
                filename = str(f)
        if filename is not None :
            filename=autopack.retrieveFile(filename,
                            destination = recipe+os.sep+"recipe"+os.sep+"ingredients"+os.sep,
                            cache="recipes")
            from xml.dom.minidom import parse
            xmlingr = parse(filename) # parse an XML file by name
            ingrnode = xmlingr.documentElement
        elif inode is not None:
            ingrnode = inode
        else :
            print ("filename is None")
            return None
        kw=self.parseIngrXmlNode(ingrnode)
        #check for overwritten parameter
        overwrite_node = inode.getElementsByTagName("overwrite")
        if overwrite_node :
            kwo=self.parseIngrXmlNode(overwrite_node[0])
            kw.update(kwo)
        ingre = self.makeIngredient(**kw)                    
        return ingre

    def parseIngrXmlNode(self,ingrnode):
        name = str(ingrnode.getAttribute("name"))
        kw = {}
        for k in KWDS:
            v=io.getValueToXMLNode(KWDS[k]["type"],ingrnode,k)
#example of debugging...
#            if k == "rejectionThreshold" :
#                print "rejectionThreshold",KWDS[k]["type"],v,v is not None
#                print "rejectionThreshold",ingrnode.getAttribute(k)
            if v is not None :
                kw[k]=v                   
        #create the ingredient according the type
#        ingre = self1.makeIngredient(**kw)                    
        return kw 
                
    def ingrXmlNode(self,ingr,xmldoc=None):
        rxmldoc=False
        if xmldoc is None : 
            rxmldoc = True
            impl = getDOMImplementation()
            #what about afviewer
            xmldoc = impl.createDocument(None, "ingredient", None)
            ingrnode = xmldoc.documentElement
            ingrnode.setAttribute("name",str(ingr.name))
        else :
            ingrnode = xmldoc.createElement("ingredient")
            ingrnode.setAttribute("name",str(ingr.name))
        for k in ingr.KWDS:
            v = getattr(ingr,k)
            io.setValueToXMLNode(v,ingrnode,k)
        if rxmldoc :
            return ingrnode,xmldoc
        else :
            return ingrnode

    def ingrJsonNode(self,ingr):
        ingdic={}
        for k in ingr.KWDS:
            v = getattr(ingr,k)
            if hasattr(v,"tolist"):
                v=v.tolist()
            ingdic[k] = v
        ingdic["results"]=[] 
        for r in ingr.results:  
            if hasattr(r[0],"tolist"):
                r[0]=r[0].tolist()
            if hasattr(r[1],"tolist"):
                r[1]=r[1].tolist()
            ingdic["results"].append([r[0],r[1]])
        if isinstance(ingr, GrowIngrediant) or isinstance(ingr, ActinIngrediant):
            ingdic["nbCurve"]=ingr.nbCurve
            for i in range(ingr.nbCurve):
                lp = numpy.array(ingr.listePtLinear[i])
                ingr.listePtLinear[i]=lp.tolist()                 
                ingdic["curve"+str(i)] = ingr.listePtLinear[i]
        return ingdic

    def ingrPythonNode(self,ingr,recipe="recipe"):
        inrStr="#include as follow : execfile('pathto/"+ingr.name+".py',globals(),{'recipe':recipe_variable_name})\n"
        if ingr.Type == "MultiSphere":
            inrStr+="from autopack.Ingredient import SingleSphereIngr, MultiSphereIngr\n"
            inrStr+=ingr.name+"= MultiSphereIngr( \n"  
        if ingr.Type == "MultiCylinder":
            inrStr+="from autopack.Ingredient import MultiCylindersIngr\n"
            inrStr+=ingr.name+"= MultiCylindersIngr( \n"                
        for k in ingr.KWDS:
            v = getattr(ingr,k)
            aStr = io.setValueToPythonStr(v,k)
            if aStr is not None :
                inrStr+=aStr+",\n" 
        inrStr+=")\n"
        inrStr+=recipe+".addIngredient("+ingr.name+")\n"        
        return inrStr

    def makeIngredient(self,**kw):
#        from autopack.Ingredient import SingleSphereIngr, MultiSphereIngr,SingleCubeIngr
#        from autopack.Ingredient import MultiCylindersIngr, GrowIngrediant
        ingr = None
        if kw["Type"]=="SingleSphere":
            kw["position"] = kw["positions"][0][0]
            kw["radius"]=kw["radii"][0][0]
            del kw["positions"]
            del kw["radii"]
            ingr = SingleSphereIngr(**kw)
        elif kw["Type"]=="MultiSphere":
            ingr = MultiSphereIngr(**kw)                    
        elif kw["Type"]=="MultiCylinder":
            ingr = MultiCylindersIngr(**kw)                    
        elif kw["Type"]=="SingleCube":
            kw["positions"]=[[[0,0,0],[0,0,0],[0,0,0],]]
            kw["positions2"]=None
            ingr = SingleCubeIngr(**kw)                    
        elif kw["Type"]=="Grow":
            ingr = GrowIngrediant(**kw)                    
        elif kw["Type"]=="Actine":
            ingr = ActinIngrediant(**kw)       
        if "gradient" in kw and kw["gradient"] != "" and kw["gradient"]!= "None":
            ingr.gradient = kw["gradient"]           
        return ingr    
    
