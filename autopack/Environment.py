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
@author: Graham Johnson, Ludovic Autin, & Michel Sanner

# Hybrid version merged from Graham's Sept 2011 and Ludo's April 2012 version on May 16, 2012
# Updated with final thesis HistoVol.py file from Sept 25, 2012 on July 5, 2012 with correct analysis tools

# TODO: fix the save/restore grid
"""
print ('histovol is on***********************************')

import os

from math import floor, ceil, fabs, sqrt,exp,cos
import time
from time import time

from random import randint, random, uniform, gauss, seed
import pdb
import bisect

import numpy, pickle, weakref

print ('AF import')
from autopack.Compartment import CompartmentList
from autopack.Recipe import Recipe
from autopack.Ingredient import GrowIngrediant,ActinIngrediant
from autopack.ray import vlen, vdiff, vcross
from autopack.IOutils import GrabResult

import math
import sys

if sys.version > "3.0.0":
    xrange = range
    import urllib.request as urllib# , urllib.parse, urllib.error
else :
    import urllib

SEED=14

from operator import itemgetter, attrgetter
from .randomRot import RandomRot

import autopack
try :
    helper = autopack.helper
except :
    helper = None
print ("HistoVol helper is "+str(helper))

LOG = False
verbose = 0

#PANDA3D Physics engine ODE and Bullet
try :
    import panda3d
    print ("got Panda3D raw")
except :    
    #MAC PATH
    p="/Developer/Panda3D/lib"#sys.path.append("/Developer/Panda3D/lib/")
    sys.path.append(p)
    print ("Trying Panda3D Except")
try :
    import panda3d
    
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
    print ("Failed to get Panda")

#coul replace by a faster json python library
import json

LISTPLACEMETHOD = autopack.LISTPLACEMETHOD
    
def linearDecayProb():
    """ return a number from 0 (higest probability) to 1 (lowest probability)
    with a linear fall off of probability
    """
    r1 = uniform(-0.25,0.25)
    r2 = uniform(-0.25,0.25)
    return abs(r1+r2)
    #r3 = uniform(-0.25,0.25)
    #r4 = uniform(-0.25,0.25)
    #return abs(r1+r2+r3+r4)  # Gaussian decay

def ingredient_compare1(x, y):
    """
    sort ingredients using decreasing priority and decreasing radii for
    priority ties and decreasing completion for radii ties
    for priority > 0
    """
    p1 = x.packingPriority
    p2 = y.packingPriority
    if p1 < p2: # p1 > p2
        return 1
    elif p1==p2: # p1 == p1
       r1 = x.minRadius
       r2 = y.minRadius
       if r1 > r2: # r1 < r2
           return 1
       elif r1==r2: # r1 == r2
           c1 = x.completion
           c2 = y.completion
           if c1 > c2: # c1 > c2
               return 1
           elif c1 == c2:
               return 0
           else:
               return -1
       else:
           return -1
    else:
       return -1

def ingredient_compare0(x, y):
    """
    sort ingredients using decreasing priority and decreasing radii for
    priority ties and decreasing completion for radii ties
    for priority < 0
    """
    p1 = x.packingPriority
    p2 = y.packingPriority
    if p1 > p2: # p1 > p2
        return 1
    elif p1==p2: # p1 == p1
       r1 = x.minRadius
       r2 = y.minRadius
       if r1 > r2: # r1 < r2
           return 1
       elif r1==r2: # r1 == r2
           c1 = x.completion
           c2 = y.completion
           if c1 > c2: # c1 > c2
               return 1
           elif c1 == c2:
               return 0
           else:
               return -1
       else:
           return -1
    else:
       return -1


def ingredient_compare2(x, y):
    """
    sort ingredients using decreasing radii and decresing completion
    for radii matches:
    priority = 0
    """
    c1 = x.minRadius
    c2 = y.minRadius
    if c1 < c2:
        return 1
    elif c1==c2:
       r1 = x.completion
       r2 = y.completion
       if r1 > r2:
          return 1
       elif r1 == r2:
          return 0
       else:
          return -1
    else:  #x < y
       return -1
       
def cmp2key(mycmp):
    """Converts a cmp= function into a key= function"""
    class K:
        def __init__(self, obj, *args):
            self.obj = obj
        def __cmp__(self, other):
            return mycmp(self.obj, other.obj)
    return K       


def vector_norm(data, axis=None, out=None):
    """get the norm of a vector data"""
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
        
def angleVector(v0,v1, directed=True, axis=0):
    """get the angle between two vector v0 and v1"""
    v0 = numpy.array(v0, dtype=numpy.float64, copy=False)
    v1 = numpy.array(v1, dtype=numpy.float64, copy=False)
    dot = numpy.sum(v0 * v1, axis=axis)
    dot /= vector_norm(v0, axis=axis) * vector_norm(v1, axis=axis)
    return numpy.arccos(dot if directed else numpy.fabs(dot))

def vdistance(c0,c1):
    """get the distance between two points c0 and c1"""
    d = numpy.array(c1) - numpy.array(c0)
    s = numpy.sum(d*d)
    return math.sqrt(s)

class Gradient:
    """
    The Gradient class
    ==========================
    This class handle the use of gradient to control the packing.
    The class define different and setup type of gradient, as well as the sampling function
    """
    
    def __init__(self,name,mode="X",description="",direction=None,bb=None,**kw):
        self.name=name
        self.description=description
        self.start=[]
        self.end=[]
        self.bb=[[],[]]
        if bb is not None:
            self.computeStartEnd()
        self.function=self.defaultFunction #lambda ? 
        self.weight = None
        self.liste_mode = ["X","Y","Z","-X","-Y","-Z","direction","radial"]
        self.mode = mode #can X,Y,Z,-X,-Y,-Z,"direction" custom vector 
        self.weight_mode = "gauss"#"linear" #linear mode for weight generation linearpos linearneg
        if "weight_mode" in kw :
            self.weight_mode = kw["weight_mode"]
        self.pick_mode = "rnd"
        if "pick_mode" in kw :
            self.pick_mode = kw["pick_mode"]
        self.axes = {"X":0,"-X":0,"Y":1,"-Y":1,"Z":2,"-Z":2}
        self.directions = {"X":[1,0,0],"-X":[-1,0,0],"Y":[0,1,0],"-Y":[0,-1,0],"Z":[0,0,1],"-Z":[0,0,-1]}
        self.radius=10.0
        if "radius" in kw :
            self.radius=kw["radius"]
        self.weight_threshold = 0.0
        if direction is None :
            self.direction = self.directions[self.mode]
        else :
#            self.mode = "direction"
            self.direction=direction#from direction should get start and end point of the gradient
        self.distance=0.0
        self.gblob = 4.0
        #Note : theses functions could also be used to pick an ingredient
        self.pick_functions = {"max":self.getMaxWeight,
                               "min":self.getMinWeight,
                               "rnd":self.getRndWeighted,
                               "linear":self.getLinearWeighted,
                               "binary":self.getBinaryWeighted,
                               "sub":self.getSubWeighted,
                               "reg":self.getForwWeight,
                               }    
        self.liste_weigth_mode = self.pick_functions.keys()
        self.liste_options = ["mode","weight_mode","pick_mode","direction","radius","gblob"]
        self.OPTIONS = {
                    "mode":{"name":"mode","values":self.liste_mode,"default":"X",
                                           "type":"liste","description":"gradient direction",
                                           "min":0,"max":0},
                    "weight_mode":{"name":"weight_mode","values":["linear","square","cube","gauss","half-gauss"],"default":"linear","type":"liste",
                                           "description":"calcul of the weight method","min":0,"max":0},
                    "pick_mode":{"name":"weight_mode","values":self.liste_weigth_mode,"default":"linear","type":"liste",
                                           "description":"picking random weighted method","min":0,"max":0},
                    "direction":{"name":"direction","value":[0.5,0.5,0.5],"default":[0.5,0.5,0.5],
                                 "type":"vector","description":"gradient custom direction","min":-2000.0,"max":2000.0},
                    "description":{"name":"description","value":self.description,"default":"a gradient",
                                 "type":"label","description":None,"min":0,"max":0}, 
                    "radius":{"name":"radius","value":self.radius,"default":100.0,
                                 "type":"float","description":"radius for the radial mode","min":0,"max":2000.0}, 
                    "gblob":{"name":"gblob","value":self.gblob,"default":4.0,
                                 "type":"float","description":"bobliness the gaussian mode","min":0.1,"max":2000.0}, 
                                 
                    }

    def getCenter(self):
        """get the center of the gradient grid"""
        center=[0.,0.,0.]
        for i in range(3):
            center[i]=(self.bb[0][i]+self.bb[1][i])/2.
        return center
                
    def computeStartEnd(self):
        """get the overal direction of the gradient"""
        #using bb and direction
        self.start=numpy.array(self.bb[0])
        self.end = numpy.array(self.bb[1])*numpy.array(self.direction)
        self.vgradient = self.end - self.start
        #self.distance = math.sqrt(numpy.sum(d*d))
        
    def defaultFunction(self,xyz):
        """
        #linear function 0->0.1
        #project xyz on direction
        """
        x = numpy.dot(xyz,self.direction)
        v = (x * 1.0) / (self.distance)
        return v

    def pickPoint(self,listPts):
        """
        pick next random point according to the choosen function
        """
        return self.pick_functions[self.pick_mode](listPts)
                
    def buildWeigthMap(self,bb,MasterPosition):
        """
        build the actual gradient value according the gradint mode
        """
        print ("gradient ",self.name,self.mode)
        if self.mode in self.axes :
            self.buildWeigthMapAxe(bb,MasterPosition)
        elif self.mode == "direction":
            self.buildWeigthMapDirection(bb,MasterPosition)
        elif  self.mode == "radial":
            self.buildWeigthMapRadial(bb,MasterPosition)
    
    def get_gauss_weights(self,N,degree=5):
        """
        given a number of point compute the gaussian weight for each
        """
        degree = N/2
        window=N#degree*2#-1  
        weight=numpy.array([1.0]*window)          
        weightGauss=[]  
        for i in range(window):  
            i=i-degree+1  
            frac=i/float(window)  
            gauss=1/(numpy.exp((self.gblob*(frac))**2))  
            weightGauss.append(gauss) 
        return numpy.array(weightGauss)*weight  

    def get_gauss_weights1(self,N):
        """
        given a number of point compute the gaussian weight for each (alternative function)
        """
        support_points = [(float(3 * i)/float(N))**2.0 for i in range(-N,N + 1)]
        gii_factors = [exp(-(i/2.0)) for i in support_points]
        ki = float(sum(gii_factors))
        return [giin/ki for giin in gii_factors]

        
    def getDirectionLength(self,bb=None,direction=None):
        if direction is None :
            direction = self.direction
        if bb is None :
            bb = self.bb
        print (bb)    
        #assume grid orthogonal
        maxinmini=[]
        a=[]
        axes = ["X","Y","Z"]
        for i,ax in enumerate( axes ):
            angle=angleVector(self.directions[ax],direction)
            a.append(angle)
            maxi=max(bb[1][i],bb[0][i])
            mini=min(bb[1][i],bb[0][i])
            maxinmini.append([mini,maxi])
        m = min(a)
        axi = a.index(m)
        L = maxinmini[axi][1]-maxinmini[axi][0]
        vdot = numpy.dot(numpy.array(self.directions[axes[axi]]),numpy.array(direction))#cos a * |A|*|B|
        Ld = (1.0/vdot)*(cos(m)*L)
        return Ld,maxinmini
        
    def get_gauss_weights1(self,N):
        support_points = [(float(3 * i)/float(N))**2.0 for i in range(-N,N + 1)]
        gii_factors = [exp(-(i/2.0)) for i in support_points]
        ki = float(sum(gii_factors))
        return [giin/ki for giin in gii_factors]
 
    def buildWeigthMapRadial(self,bb,MasterPosition):
        """
        from a given point (self.direction) build a radial weigth according the choosen mode
        (linear, gauss, etc...)
        """
        N = len(MasterPosition)
        self.bb= bb
        radial_point = self.direction
        NW=N/3
        self.weight =[]
        center = self.getCenter()
        xl,yl,zl = bb[0]
        xr,yr,zr = bb[1]

        if self.weight_mode == "gauss" :#0-1-0
            d = self.get_gauss_weights(NW)#numpy.random.normal(0.5, 0.1, NW) #one dimension 
        elif self.weight_mode == "half-gauss" :#0-1
            d = self.get_gauss_weights(NW*2)[NW:]
        print (self.name,self.radius)
        for ptid in range(N) :
            dist=vdistance(MasterPosition[ptid],radial_point)
            if self.weight_mode == "linear" :
                w = (1.0-(abs(dist)/self.radius)) if abs(dist) < self.radius else 0.0
                self.weight.append( w )#
            elif self.weight_mode == "square" :
                w = math.pow((1.0-(abs(dist)/self.radius)),2) if abs(dist) < self.radius else 0.0  
                self.weight.append( w )#
            elif self.weight_mode == "cube" :
                w = math.pow((1.0-(abs(dist)/self.radius)),3) if abs(dist) < self.radius else 0.0  
                self.weight.append( w )#
            elif self.weight_mode == "gauss" :
                w = abs(dist)/self.radius if abs(dist) < self.radius else 1.0
                i = int(w*N/3) if int(w*N/3) < len(d) else len(d)-1
                self.weight.append( d[i] )
            elif self.weight_mode == "half-gauss":
                w = abs(dist)/self.radius if abs(dist) < self.radius else 1.0
#                w = (1.0-(abs(dist)/self.radius)) if abs(dist) < self.radius else 0.0
                i = int(w*N/3) if int(w*N/3) < len(d) else len(d)-1
                self.weight.append( d[i] )
       
    def buildWeigthMapDirection(self,bb,MasterPosition):
        """
        from a given direction build a linear weigth according the choosen mode
        (linear, gauss, etc...)
        """
        N = len(MasterPosition)
        self.bb= bb
        axe = self.direction
        NW=N/3
        self.weight =[]
        center = self.getCenter()
        L,maxinmini = self.getDirectionLength(bb)
        if self.weight_mode == "gauss" :
            d = self.get_gauss_weights(NW)#numpy.random.normal(0.5, 0.1, NW) #one dimension 
        elif self.weight_mode == "half-gauss" :#0-1
            d = self.get_gauss_weights(NW*2)[NW:]
        for ptid in range(N) :
            pt = numpy.array(MasterPosition[ptid])-numpy.array(center)#[maxinmini[0][0],maxinmini[1][0],maxinmini[2][0]])
            vdot = numpy.dot(pt,numpy.array(axe))
            p = ((L/2.0)+vdot)/L
            if self.weight_mode == "linear" :
                self.weight.append( p )#-0.5->0.5 axe value normalized?
            elif self.weight_mode == "square" :
                self.weight.append( math.pow(p,2) )#
            elif self.weight_mode == "cube" :
                self.weight.append( math.pow(p,3) )#
            elif self.weight_mode == "gauss" :
#                p goes from 0.0 to 1.0
                if p < 0.1 : p = 0.0
                i = int(p*NW) if int(p*NW) < len(d) else len(d)-1
                #w = d[i] if d[i] > 0.9 else 0.0
                self.weight.append( d[i] )

    def buildWeigthMapAxe(self,bb,MasterPosition,Axe="X"):
        """
        from a given axe (X,Y,Z) build a linear weigth according the choosen mode
        (linear, gauss, etc...)
        """
        N = len(MasterPosition)
        self.bb= bb
        ind = self.axes[self.mode]
        maxi=max(bb[1][ind],bb[0][ind])
        mini=min(bb[1][ind],bb[0][ind])
        maxix=max(bb[1][0],bb[0][0])
        minix=min(bb[1][0],bb[0][0])
        self.weight =[]
        if self.weight_mode == "gauss" :
            d = self.get_gauss_weights(N/3)#d = numpy.random.normal(0.5, 0.1, N/3) #one dimension 
        elif self.weight_mode == "half-gauss" :#0-1
            d = self.get_gauss_weights(NW*2)[NW:]
        for ptid in range(N) :
            p=(MasterPosition[ptid][ind]-mini)/(maxi-mini)
            if self.weight_mode == "linear" :
                self.weight.append( p )#-0.5->0.5 axe value normalized?
            elif self.weight_mode == "square" :
                self.weight.append( math.pow(p,2) )#
            elif self.weight_mode == "cube" :
                self.weight.append( math.pow(p,3) )#
            elif self.weight_mode == "gauss" :
                vax=p#(MasterPosition[ptid][ind]-mini)/(maxi-mini) #0-1 on the axes
                i = int(vax*N/3) if int(vax*N/3) < len(d) else len(d)-1
                self.weight.append( d[i] )
            elif self.weight_mode == "half-gauss":
                #w = abs(dist)/self.radius if abs(dist) < self.radius else 1.0
#                w = (1.0-(abs(dist)/self.radius)) if abs(dist) < self.radius else 0.0
                i = int(p*N/3) if int(p*N/3) < len(d) else len(d)-1
                self.weight.append( d[i] )
                
    def getMaxWeight(self,listPts):
        """
        from the given list of grid point indice, get the point with the maximum weight
        """
        ptInd = listPts[0]
        m=0.0                
        for pi in listPts :
            if self.weight[pi] > m :
                m=self.weight[pi]
                ptInd = pi
        if self.weight[ptInd] < self.weight_threshold:
            ptInd= None
        #print "picked",self.weight[ptInd]
        return ptInd
        
    def getMinWeight(self,listPts):   
        """
        from the given list of grid point indice, get the point with the minimum weight
        """
        ptInd = listPts[0]
        m=1.1         
        for pi in  listPts :
            if self.weight[pi] < m and self.weight[pi] != 0:
                m=self.weight[pi]
                ptInd = pi
        if self.weight[ptInd] < self.weight_threshold:
            return None
        return ptInd
                
    def getRndWeighted(self,listPts):
        """
        From http://glowingpython.blogspot.com/2012/09/weighted-random-choice.html
        Weighted random selection
        returns n_picks random indexes.
        the chance to pick the index i 
        is give by the weight weights[i].
        """
        weight = numpy.take(self.weight,listPts)
        t = numpy.cumsum(weight)
        s = numpy.sum(weight)
        i = numpy.searchsorted(t,numpy.random.rand(1)*s)[0]
        return listPts[i]

    def getLinearWeighted(self,listPts):
        """
        From http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/
        The following is a simple function to implement weighted random selection in Python. 
        Given a list of weights, it returns an index randomly, according to these weights [2].
        For example, given [2, 3, 5] it returns 0 (the index of the first element) with probability 0.2, 
        1 with probability 0.3 and 2 with probability 0.5. 
        The weights need not sum up to anything in particular, 
        and can actually be arbitrary Python floating point numbers.
        """
        totals = []
        running_total = 0
        weights = numpy.take(self.weight,listPts)
        for w in weights:
            running_total += w
            totals.append(running_total)
    
        rnd = random() * running_total
        for i, total in enumerate(totals):
            if rnd < total:
                return listPts[i]

    def getBinaryWeighted(self,listPts):
        """
        From http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/
        Note that the loop in the end of the function is simply looking 
        for a place to insert rnd in a sorted list. Therefore, it can be 
        speed up by employing binary search. Python comes with one built-in, 
        just use the bisect module.
        """
        totals = []
        running_total = 0
        
        weights = numpy.take(self.weight,listPts)
        for w in weights:
            running_total += w
            totals.append(running_total)
    
        rnd = random() * running_total
        i = bisect.bisect_right(totals, rnd)
        return listPts[i]

    def getForwWeight(self,listPts):
        dice = random()
        #sorted ?
        for i in listPts :
            if self.weight[i] > dice and self.weight[i] != 0 :
                return i

    def getSubWeighted(self,listPts):
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
        
        weights = numpy.take(self.weight,listPts)
        rnd = random() * sum(weights)
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                return listPts[i]
                
class Grid:
    """
    The Grid class
    ==========================
    This class handle the use of grid to control the packing. The grid keep information of
    3d positions, distances, freePoints and inside/surface points from organelles.
    NOTE : thi class could be completly replace if openvdb is wrapped to python.
    """
    def __init__(self):
        #a grid is attached to an environement
        self.boundingBox=([0,0,0], [.1,.1,.1])
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
        self.gridSpacing = None
        self.nbGridPoints = None
        self.nbSurfacePoints = 0
        self.gridVolume = 0 # will be the toatl number of grid points
        # list of (x,y,z) for each grid point (x index moving fastest)
        self.masterGridPositions = []
        
        #this are specific for each compartment
        self.aInteriorGrids = []
        self.aSurfaceGrids = []
        #bhtree
        self.surfPtsBht=None
        self.ijkPtIndice = []
        self.filename=None          #used for storing before fill so no need rebuild
        self.result_filename=None   #used after fill to store result

        self.encapsulatingGrid = 1
        
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

    def create3DPointLookup(self, boundingBox=None):
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
        self.ijkPtIndice = numpy.zeros( (nx*ny*nz, 3), 'i')
        space = self.gridSpacing
        # Vector for lower left broken into real of only the z coord.
        i = 0
        for zi in range(nz):
            for yi in range(ny):
                for xi in range(nx):
                    pointArrayRaw[i] = (xl+xi*space, yl+yi*space, zl+zi*space)
                    self.ijkPtIndice[i] = (xi,yi,zi)
                    i+=1
        self.masterGridPositions = pointArrayRaw

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
            return False
        else :
            if dist is not None:
                #distance to closest wall
                d1 = P - O
                s1=min(x for x in (d1*jitter) if x != 0)
                #s1 = numpy.sum(d1*d1)
                d2 = E - P
                s2=min(x for x in (d2*jitter) if x != 0)
                #s2 = numpy.sum(d2*d2)
                if s1 <= dist or s2 <=dist:
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
        
    def getPointsInCube(self, bb, pt, radius,addSP=True,info=False):
        """
        Return all grid points indicesinside the given bouding box.
        """        
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

        encapsulatingGrid = self.encapsulatingGrid  #Graham Added on Oct17 to allow for truly 2D grid for test fills... may break everything!
        
        from math import ceil
        nx = int(ceil((xr-xl)/space))+encapsulatingGrid
        ny = int(ceil((yr-yl)/space))+encapsulatingGrid
        nz = int(ceil((zr-zl)/space))+encapsulatingGrid
        return nx*ny*nz,(nx,ny,nz)

#==============================================================================
# TO DO File IO 
#==============================================================================
    def save(self):
        pass
    def restore(self):
        pass
        
class Environment(CompartmentList):
    """
    The Environment class
    ==========================
    This class is main class in autopack. The class handle all the setup, initialisation
    and process of the packing. We use xml or python for the setup.
    A environments is made of :
        a grid and the gradients if any
        a list of compartment and their recipes (surface and interior)
        a exterior recipe
        each recipe are made of a list of ingredients
    """

    def __init__(self,name="H"):
        CompartmentList.__init__(self)
#        self.compartments = []
        self.verbose = verbose  #Graham added to try to make universal "global variable Verbose" on Aug 28
        self.timeUpDistLoopTotal = 0 #Graham added to try to make universal "global variable Verbose" on Aug 28
        self.name = name        
        self.exteriorRecipe = None
        self.hgrid=[]
        self.world = None   #panda world for collision
        self.octree = None  #ongoing octree test, no need if openvdb wrapp to python
        self.grid = Grid()  # the main grid
        self.encapsulatingGrid = 1  # Only override this with 0 for 2D fills- otherwise its very unsafe!
        self.nbCompartments = 1 # 0 is the exterior, 1 is compartment 1 surface, -1 is compartment 1 interior, etc.
        self.name = "out"

        self.order={}#give the order of drop ingredient by ptInd from molecules
        self.lastrank=0

        # smallest and largest protein radii acroos all recipes
        self.smallestProteinSize = 99999999
        self.largestProteinSize = 0
        self.scaleER = 20.0 # hack in case problem with encapsulating radius
        self.computeGridParams = True
        
        self.EnviroOnly = False
        self.EnviroOnlyCompartiment  =  -1
        # bounding box of the HistoVol

        self.boundingBox = [[0,0,0], [.1,.1,.1]]
        self.fbox_bb = None #used for estimating the volume
        
        self.fbox = None #Oct 20, 2012 Graham wonders if this is part of the problem
        self.fillBB = None # bounding box for a given fill
        self.fillbb_insidepoint =None #Oct 20, 2012 Graham wonders if this is part of the problem
        self.freePointMask =None
        self.molecules = [] # list of ( (x,y,z), rotation, ingredient) triplet generated by packing 
        self.ingr_result = {}
        self.ingr_added = {}
        
        self.randomRot = RandomRot()# the class used to generate random rotation
        self.activeIngr= [] 
        self.activeIngre_saved = []
        
        #optionally can provide a host and a viewer
        self.host = None
        self.afviewer = None
        
        #version of setup used
        self.version = "1.0"
        
        #option for filling using host dynamics capability
        self.windowsSize = 100
        self.windowsSize_overwrite = False

        
        self.runTimeDisplay = False
        self.placeMethod="jitter"
        self.innerGridMethod = "bhtree" #or sdf
        self.orthogonalBoxType = 0
        self.overwritePlaceMethod = False

        #if use C4D RB dynamics, should be genralized
        self.springOptions={}
        self.dynamicOptions={}
        self.setupRBOptions()
        self.simulationTimes = 2.0
        
        #saving/pickle option
        self.saveResult = False
        self.resultfile = "fillResult"
        self.setupfile = ""
        self.grid_filename =None#
        self.grid_result_filename = None#str(gridn.getAttribute("grid_result"))

        #cancel dialog-> need to be develop more
        self.cancelDialog = False

        self.grab_cb = None 
        self.pp_server = None
        
        #
        self.nFill = 0
        self.cFill = 0
        self.FillName=[]
        
        ##do we sort the ingrediant or not see  getSortedActiveIngredients
        self.pickWeightedIngr = True 
        self.pickRandPt = True       ##point pick randomly or one after the other?
        self.currtId = 0
        
        #gradient
        self.gradients={}
        
        self.use_gradient=False #gradient control is also per ingredient
        self.use_halton=False #use halton for grid point distribution
        
        self.ingrLookForNeighbours = True #Old Features to be test
        
        #debug with timer function
        self._timer = False
        self._hackFreepts = False   #hack for speed up ?
        self.freePtsUpdateThrehod = 0.0
        self.nb_ingredient=0
        self.totalNbIngr = 0
        self.treemode = "cKDTree"# "bhtree"
        self.close_ingr_bhtree=None#RBHTree(a.tolist(),range(len(a)),10,10,10,10,10,9999)
        self.rTrans=[]  
        self.result=[]
        self.rIngr=[] 
        self.rRot=[] 
        self.listPlaceMethod = LISTPLACEMETHOD
        #should be part of an independant module 
        self.panda_solver="bullet" #or bullet
        #could be a problem here for pp
        #can't pickle this dictionary
        self.rb_func_dic={}
        
        self.OPTIONS = {
                    "smallestProteinSize":{"name":"smallestProteinSize","value":15,"default":15,
                                           "type":"int","description":"Smallest ingredient packing radius override (low=accurate | high=fast)",
                                           "mini":1.0,"maxi":100.0,
                                           "width":30},
                    "largestProteinSize":{"name":"largestProteinSize","value":0,"default":0,"type":"int","description":"largest Protein Size","width":30},
                    "computeGridParams":{"name":"computeGridParams","value":True,"default":True,"type":"bool","description":"compute Grid Params","width":100},
                    "EnviroOnly":{"name":"EnviroOnly","value":False,"default":False,"type":"bool","description":"Histo volume Only","width":30},
#                    "EnviroOnlyCompartiment":{"name":"EnviroOnlyCompartiment","value":-1,"default":-1,"type":"int","description":"Histo volume Only compartiment"},
                    "windowsSize":{"name":"windowsSize","value":100,"default":100,"type":"int","description":"windows Size","width":30},
#                    "simulationTimes":{"name":"simulationTimes","value":90,"default":90,"type":"int","description":"simulation Times"},
                    "runTimeDisplay":{"name":"runTimeDisplay","value":False,"default":False,"type":"bool","description":"Display packing in realtime (slow)","width":150},
                    "placeMethod": {"name":"placeMethod","value":"jitter","values":self.listPlaceMethod,"default":"placeMethod","type":"liste","description":"     Overriding Packing Method = ","width":30},
                    "use_gradient":{"name":"use_gradient","value":False,"default":False,"type":"bool","description":"Use gradients if defined","width":150},
                    "gradients":{"name":"gradients","value":"","values":[],"default":"","type":"liste","description":"Gradients available","width":150},
                    "innerGridMethod": {"name":"innerGridMethod","value":"bhtree","values":["bhtree","sdf","jordan","jordan3","pyray"],"default":"innerGridMethod","type":"liste","description":"     Method to calculate the inner grid:","width":30},
                    "overwritePlaceMethod":{"name":"overwritePlaceMethod","value":False,"default":False,"type":"bool","description":"Overwrite per-ingredient packing method with Overriding Packing Method:","width":300},
                    "saveResult": {"name":"saveResult","value":False,"default":False,"type":"bool","description":"Save packing result to .apr file (enter full path below):","width":200},
                    "resultfile": {"name":"resultfile","value":"fillResult","default":"fillResult","type":"filename","description":"result filename","width":200},
                    #cancel dialog
                    "cancelDialog": {"name":"cancelDialog","value":False,"default":False,"type":"bool","description":"compute Grid Params","width":30},
                    ##do we sort the ingrediant or not see  getSortedActiveIngredients
                    "pickWeightedIngr":{"name":"pickWeightedIngr","value":True,"default":True,"type":"bool","description":"Prioritize ingredient selection by packingWeight","width":200},
                    "pickRandPt":{"name":"pickRandPt","value":True,"default":True,"type":"bool","description":"Pick drop position point randomly","width":200},
                    #gradient
                    "ingrLookForNeighbours":{"name":"ingrLookForNeighbours","value":False,"default":False,"type":"bool","description":"Look for ingredients attractor and partner","width":30},
                    #debug with timer function
                    "_timer": {"name":"_timer","value":False,"default":False,"type":"bool","description":"evaluate time per function","width":30},
                    "_hackFreepts": {"name":"_hackFreepts","value":False,"default":False,"type":"bool","description":"no free point update","width":30},
                    "freePtsUpdateThrehod":{"name":"freePtsUpdateThrehod","value":0.15,"default":0.15,"type":"float","description":"Mask grid while packing (0=always | 1=never)","mini":0.0,"maxi":1.0,"width":30},
                        }

    def Setup(self,setupfile):
        #parse the given fill for
        #1-fillin option
        #2-recipe
        #use XML with tag description of the setup:
        #filling name root
        #histoVol option
        #cytoplasme recipe if any and its ingredient
        #compartment name= mesh ?
        #orga surfaceingr#file or direct
        #orga interioringr#file or direct
        #etc...
        pass

    def reportprogress(self,label=None,progress=None):
        if self.afviewer is not None and hasattr(self.afviewer,"vi"):
            self.afviewer.vi.progressBar(progress=progress,label=label)
        
    def makeIngredient(self,**kw):
        """
        Helper function to make an ingredient, pass all arguments as keywords.
        """
        from autopack.Ingredient import SingleSphereIngr, MultiSphereIngr,SingleCubeIngr
        from autopack.Ingredient import MultiCylindersIngr, GrowIngrediant
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
            ingr = ActineIngrediant(**kw)       
        if "gradient" in kw and kw["gradient"] != "" and kw["gradient"]!= "None":
            ingr.gradient = kw["gradient"]           
        return ingr

    def set_partners_ingredient(self,ingr):
        if ingr.partners_name :
            for i,iname in enumerate(ingr.partners_name) :
                print (iname)
                ingr_partner = self.getIngrFromName(iname)
                print (ingr_partner) 
                partner = ingr.addPartner(ingr_partner,weight=ingr_partner.weight,
                                properties={"position":ingr.partners_position[i]})
                for p in ingr_partner.properties:
                    partner.addProperties(p,ingr_partner.properties[p])
            if ingr.Type == "Grow" : ingr.prepare_alternates()
        if ingr.excluded_partners_name :
            for iname in ingr.excluded_partners_name :
                ingr.addExcludedPartner(iname)
        ingr.histoVol = self
        
    def set_recipe_ingredient(self,xmlnode,recipe,io_ingr):
        #get the defined ingredient
        ingrnodes = xmlnode.getElementsByTagName("ingredient")
        for ingrnode in ingrnodes:
            ingre = io_ingr.makeIngredientFromXml(inode = ingrnode , recipe=self.name)
            if ingre : recipe.addIngredient(ingre) 
            else : print ("PROBLEM creating ingredient from ",ingrnode,) 
        #check for includes 
        ingrnodes_include = xmlnode.getElementsByTagName("include")
        for inclnode in ingrnodes_include:
            xmlfile = str(inclnode.getAttribute("filename"))
            ingre = io_ingr.makeIngredientFromXml(filename = xmlfile, recipe=self.name)
            if ingre : recipe.addIngredient(ingre)
            else : print ("PROBLEM creating ingredient from ",ingrnode)
            #look for overwritten attribute
        
    def load_XML(self,setupfile):
        """
        Setup the environment according the given xml file. 
        """
        self.setupfile = setupfile
        from autopack import Ingredient as ingr
        from autopack.Ingredient import IOingredientTool
        io_ingr = IOingredientTool()
        from xml.dom.minidom import parse
        self.xmldoc = parse(setupfile) # parse an XML file by name
        root = self.xmldoc.documentElement
        self.name = str(root.getAttribute("name"))
        options=root.getElementsByTagName("options")
        if len(options) :
            options=options[0]
            for k in self.OPTIONS:
                if k == "gradients": 
                    continue
                v=self.getValueToXMLNode(self.OPTIONS[k]["type"],options,k)
                if v is not None :
                    setattr(self,k,v)
            v=self.getValueToXMLNode("vector",options,"boundingBox")
            self.boundingBox = v
            v=self.getValueToXMLNode("string",options,"version")
            self.version = v
            
        gradientsnode=root.getElementsByTagName("gradients")
        if len(gradientsnode) :
            gradientnode=gradientsnode[0]
            grnodes = gradientnode.getElementsByTagName("gradient")
            for grnode in grnodes:
                name = str(grnode.getAttribute("name"))
                mode = str(grnode.getAttribute("mode"))
                weight_mode = str(grnode.getAttribute("weight_mode"))
                pick_mode = str(grnode.getAttribute("pick_mode"))
                direction = str(grnode.getAttribute("direction"))#vector
                description=str(grnode.getAttribute("description"))
                radius=float(str(grnode.getAttribute("radius")))
#                print "weight_mode",weight_mode
                self.setGradient(name=name,mode=mode, direction=eval(direction),
                            weight_mode=weight_mode,description=description,
                            pick_mode=pick_mode,radius=radius)

        gridnode=root.getElementsByTagName("grid")
        if len(gridnode) :
            gridn=gridnode[0]
            self.grid_filename = str(gridn.getAttribute("grid_storage"))
            self.grid_result_filename = str(gridn.getAttribute("grid_result"))

        rnode=root.getElementsByTagName("cytoplasme")
        if len(rnode) :
            rCyto = Recipe()
            rnode=rnode[0]
            #check for include list of ingredients
            ingredients_xmlfile = str(rnode.getAttribute("include"))
            if ingredients_xmlfile :#open the file and parse the ingredient:
                #check if multiple include filename, aumngo',' in the path
                liste_xmlfile = ingredients_xmlfile.split(",")
                for xmlf in liste_xmlfile :                    
                    xmlfile = autopack.retrieveFile(xmlf,
                            destination = self.name+os.sep+"recipe"+os.sep,
                            cache="recipes")
                    if xmlfile :
                        xmlinclude = parse(xmlfile).documentElement
                        self.set_recipe_ingredient(xmlinclude,rCyto,io_ingr)
                
            self.set_recipe_ingredient(rnode,rCyto,io_ingr)                
            #setup recipe
            self.setExteriorRecipe(rCyto)
            
        onodes = root.getElementsByTagName("compartment")#Change to Compartment
        if not len(onodes) :
            #backward compatibility
            onodes = root.getElementsByTagName("organelle")#Change to Compartment            
        from autopack.Compartment import Compartment
        for onode in onodes:
            name = str(onode.getAttribute("name"))
            geom = str(onode.getAttribute("geom"))
            rep =  str(onode.getAttribute("rep"))
            rep_file=str(onode.getAttribute("rep_file"))
            if rep != "None" :
                rname =  rep_file.split("/")[-1]
                fileName, fileExtension = os.path.splitext(rname)
                if fileExtension == "" :
                    fileExtension = autopack.helper.hext
                    if fileExtension == "" :
                        rep_file = rep_file+fileExtension
                    else :
                        rep_file = rep_file+"."+fileExtension   
            else :
                rep=None
                rep_file=None
            o = Compartment(name,None, None, None,filename=geom,object_name=rep,object_filename=rep_file)
            self.addCompartment(o)
            rsnodes = onode.getElementsByTagName("surface")
            if len(rsnodes) :
                rSurf = Recipe(name=o.name+"_surf")
                rsnodes=rsnodes[0]
                ingredients_xmlfile = str(rsnodes.getAttribute("include"))
                if ingredients_xmlfile :#open the file and parse the ingredient:
                    #check if multiple include filename, aumngo',' in the path
                    liste_xmlfile = ingredients_xmlfile.split(",")
                    for xmlf in liste_xmlfile :                    
                        xmlfile = autopack.retrieveFile(xmlf,
                                destination = self.name+os.sep+"recipe"+os.sep,
                            cache="recipes")
                        if xmlfile :
                            xmlinclude = parse(xmlfile).documentElement
                            self.set_recipe_ingredient(xmlinclude,rSurf,io_ingr)
                self.set_recipe_ingredient(rsnodes,rSurf,io_ingr)                
                o.setSurfaceRecipe(rSurf)                
            rinodes = onode.getElementsByTagName("interior")
            if len(rinodes) :
                rMatrix = Recipe(name=o.name+"_int")
                rinodes=rinodes[0]
                ingredients_xmlfile = str(rinodes.getAttribute("include"))
                if ingredients_xmlfile :#open the file and parse the ingredient:
                    #check if multiple include filename, aumngo',' in the path
                    liste_xmlfile = ingredients_xmlfile.split(",")
                    for xmlf in liste_xmlfile :
                        xmlfile = autopack.retrieveFile(xmlf,
                                destination = self.name+os.sep+"recipe"+os.sep,
                            cache="recipes")
                        if xmlfile :
                            xmlinclude = parse(xmlfile).documentElement
                            self.set_recipe_ingredient(xmlinclude,rMatrix,io_ingr)
                self.set_recipe_ingredient(rinodes,rMatrix,io_ingr)                
                o.setInnerRecipe(rMatrix)
        #Go through all ingredient and setup the partner  
        self.loopThroughIngr(self.set_partners_ingredient)
        if self.placeMethod.find("panda") != -1 :
            self.setupPanda()
            
    def getValueToXMLNode(self,vtype,node,attrname):
        """
        Helper function to get the value from a given xml node attribute of a given type 
        """
        value = node.getAttribute(attrname)
        value = str(value)
        if not len(value):
            return None
        if vtype not in ["liste","filename","string"] :
            value=eval(value)
        else :
            value=str(value)
        return value

    def getStringValueOptions(self,value,attrname):
        """
        Helper function to return the given environment option as a string to
        be write in the xml file.
        """
        if value is None:
            return "None"
        if attrname == "color" :
            if type(value) != list and type(value) != tuple :
                if autopack.helper is not None : 
                    value=helper.getMaterialProperty(value,["color"])[0]
                else :
                    value = [1.,0.,0.]
        if type (value) == numpy.ndarray :
            value = value.tolist()
        elif type(value) == list :
            for i,v in enumerate(value) :
                if type(v) == numpy.ndarray :
                    value[i] = v.tolist()
                elif type(v) == list :
                    for j,va in enumerate(v) :
                        if type(va) == numpy.ndarray :
                            v[j] = va.tolist() 
        if type(value) == str :
            value = '"'+value+'"'
        return str(value)
        
    def setValueToXMLNode(self,value,node,attrname):
        """
        Helper function to apply the given environment option value to the 
        given xml node.
        """
        if value is None:
            return
        if attrname == "color" :
            if type(value) != list and type(value) != tuple :
                if autopack.helper is not None : 
                    value=helper.getMaterialProperty(value,["color"])[0]
                else :
                    value = [1.,0.,0.]
        if type (value) == numpy.ndarray :
            value = value.tolist()
        elif type(value) == list :
            for i,v in enumerate(value) :
                if type(v) == numpy.ndarray :
                    value[i] = v.tolist()
                elif type(v) == list :
                    for j,va in enumerate(v) :
                        if type(va) == numpy.ndarray :
                            v[j] = va.tolist()                        
#        print ("setValueToXMLNode ",attrname,value,str(value))  
        node.setAttribute(attrname,str(value))
            
    def save_asXML(self,setupfile,useXref=True):
        """
        Save the current environment setup as an xml file.
        """
        from autopack.Ingredient import IOingredientTool
        io_ingr = IOingredientTool()
        self.setupfile = setupfile
        pathout=os.path.dirname(os.path.abspath(self.setupfile))
        #export all information as xml
        #histovol is a tag, option are attribute of the tag
        from xml.dom.minidom import getDOMImplementation
        impl = getDOMImplementation()
        #what about afviewer
        self.xmldoc = impl.createDocument(None, "autopackSetup", None)
        root = self.xmldoc.documentElement
        root.setAttribute("name",str(self.name))
        options=self.xmldoc.createElement("options")
        for k in self.OPTIONS:
            v = getattr(self,k)
            if k == "gradients" :
                v = self.gradients.keys()
#            elif k == "runTimeDisplay"
            self.setValueToXMLNode(v,options,k)
        #add the boundin box
        self.setValueToXMLNode(self.boundingBox,options,"boundingBox")
        self.setValueToXMLNode(self.version,options,"version")#version?
        root.appendChild(options)
        
        if len(self.gradients):
            gradientsnode=self.xmldoc.createElement("gradients")
            root.appendChild(gradientsnode)
            for gname in self.gradients:
                g = self.gradients[gname]
                grnode = self.xmldoc.createElement("gradient")
                gradientsnode.appendChild(grnode)
                grnode.setAttribute("name",str(g.name))
                for k in g.OPTIONS:
                    v = getattr(g,k)
                    self.setValueToXMLNode(v,grnode,k)      

        #grid path information
        if self.grid.filename is not None or self.grid.result_filename is not None:
            gridnode=self.xmldoc.createElement("grid")
            root.appendChild(gridnode)
            gridnode.setAttribute("grid_storage",str(self.grid.filename))
            gridnode.setAttribute("grid_result",str(self.grid.result_filename))
        
        r =  self.exteriorRecipe
        if r :
            rnode=self.xmldoc.createElement("cytoplasme")
            root.appendChild(rnode)
            for ingr in r.ingredients:                
                if useXref :
                    io_ingr.write(ingr,pathout+os.sep+ingr.name,ingr_format="xml")
                    ingrnode = self.xmldoc.createElement("ingredient")
                    rnode.appendChild(ingrnode)
                    ingrnode.setAttribute("include",str(pathout+os.sep+ingr.name+".xml"))                    
                else :
                    ingrnode = self.xmldoc.createElement("ingredient")
                    rnode.appendChild(ingrnode)
                    ingrnode.setAttribute("name",str(ingr.name))
                    for k in ingr.KWDS:
                        v = getattr(ingr,k)
    #                    print ingr.name+" keyword ",k,v
                        self.setValueToXMLNode(v,ingrnode,k)
        for o in self.compartments:
            onode=self.xmldoc.createElement("compartment")
            root.appendChild(onode)
            onode.setAttribute("name",str(o.name))
            onode.setAttribute("geom",str(o.filename))#should point to the used filename
            onode.setAttribute("rep",str(o.representation))#None
            if o.representation is not None :
                fileName, fileExtension = os.path.splitext(o.representation_file)
            else :
                fileName = None
            onode.setAttribute("rep_file",str(fileName))#None
            rs = o.surfaceRecipe
            if rs :
                onodesurface=self.xmldoc.createElement("surface")
                onode.appendChild(onodesurface)
                for ingr in rs.ingredients: 
                    if useXref :
                        io_ingr.write(ingr,pathout+os.sep+ingr.name,ingr_format="xml")
                        ingrnode = self.xmldoc.createElement("ingredient")
                        onodesurface.appendChild(ingrnode)
                        ingrnode.setAttribute("include",str(pathout+os.sep+ingr.name+".xml")) 
                    else :
                        ingrnode = self.xmldoc.createElement("ingredient")
                        onodesurface.appendChild(ingrnode)
                        ingrnode.setAttribute("name",str(ingr.name))                       
                        for k in ingr.KWDS:
                            v = getattr(ingr,k)
                            self.setValueToXMLNode(v,ingrnode,k)
            ri = o.innerRecipe
            if ri :
                onodeinterior=self.xmldoc.createElement("interior")
                onode.appendChild(onodeinterior)             
                for ingr in ri.ingredients: 
                    if useXref :
                        io_ingr.write(ingr,pathout+os.sep+ingr.name,ingr_format="xml")
                        ingrnode = self.xmldoc.createElement("ingredient")
                        onodeinterior.appendChild(ingrnode)
                        ingrnode.setAttribute("include",str(pathout+os.sep+ingr.name+".xml")) 
                    else :
                        ingrnode = self.xmldoc.createElement("ingredient")
                        onodeinterior.appendChild(ingrnode)
                        ingrnode.setAttribute("name",str(ingr.name))                       
                        for k in ingr.KWDS:
                            v = getattr(ingr,k)
                            self.setValueToXMLNode(v,ingrnode,k)
        f = open(setupfile,"w")        
        self.xmldoc.writexml(f, indent="\t", addindent="", newl="\n")
        f.close()

    def save_asPython(self,setupfile,useXref=True):
        """
        Save the current environment setup as a python script file.
        """
        from autopack.Ingredient import IOingredientTool
        io_ingr = IOingredientTool()
        self.setupfile = setupfile
        pathout=os.path.dirname(os.path.abspath(self.setupfile))
        #add the import statement
        setupStr="""
import sys
import os
#autopack
import autopack
localdir = wrkDir = autopack.__path__[0]
from autopack.Ingredient import SingleSphereIngr, MultiSphereIngr
from autopack.Ingredient import MultiCylindersIngr,GrowIngrediant,ActinIngrediant
from autopack.Compartment import Compartment
from autopack.Recipe import Recipe
from autopack.Environment import Environment
from autopack.Graphics import AutopackViewer as AFViewer
#access the helper
helper = autopack.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper =helperClass()
#create the viewer
ViewerType=autopack.helper.host    
afviewer = AFViewer(ViewerType=helper.host,helper=helper)#long ?
#make some option here     
afviewer.doPoints = False
afviewer.doSpheres = False
afviewer.quality = 1 #lowest quality for sphere and cylinder
afviewer.visibleMesh = True #mesh default visibility 
#create the env
h1 = Environment()
"""
        setupStr+="h1.name='"+self.name+"'\n"    
        for k in self.OPTIONS:
            v = getattr(self,k)
            if k == "gradients" :
                v = self.gradients.keys()
            vstr=self.getStringValueOptions(v,k)#self.setValueToXMLNode(v,options,k)
            setupStr+="h1.%s=%s\n" % (k,vstr)
        #add the boundin box
        vstr=self.getStringValueOptions(self.boundingBox,"boundingBox")#self.setValueToXMLNode(v,options,k)
        setupStr+="h1.%s=%s\n" % ("boundingBox",vstr)
        vstr=self.getStringValueOptions(self.version,k)#self.setValueToXMLNode(v,options,k)
        setupStr+="h1.%s=%s\n" % ("version",vstr)
        
#TODO : GRADIENT
#        if len(self.gradients):
#            gradientsnode=self.xmldoc.createElement("gradients")
#            root.appendChild(gradientsnode)
#            for gname in self.gradients:
#                g = self.gradients[gname]
#                grnode = self.xmldoc.createElement("gradient")
#                gradientsnode.appendChild(grnode)
#                grnode.setAttribute("name",str(g.name))
#                for k in g.OPTIONS:
#                    v = getattr(g,k)
#                    self.setValueToXMLNode(v,grnode,k)      
#
#        grid path information
#        if self.grid.filename is not None or self.grid.result_filename is not None:
#            gridnode=self.xmldoc.createElement("grid")
#            root.appendChild(gridnode)
#            gridnode.setAttribute("grid_storage",str(self.grid.filename))
#            gridnode.setAttribute("grid_result",str(self.grid.result_filename))
#        
        r =  self.exteriorRecipe
        if r :
            setupStr+="cytoplasme = Recipe()\n"
            for ingr in r.ingredients:                
                if useXref :
                    io_ingr.write(ingr,pathout+os.sep+ingr.name,ingr_format="python")
                    setupStr+="execfile('"+pathout+os.sep+ingr.name+".py',globals(),{'recipe':cytoplasme})\n"
                else :
                    ingrnode = io_ingr.ingrPythonNode(ingr,recipe="cytoplasme")
                    setupStr+=ingrnode
            setupStr+="h1.setExteriorRecipe(cytoplasme)\n"                    
        for o in self.compartments:
            setupStr+=o.name+" = Compartment('"+o.name+"',None, None, None,\n"
            setupStr+="         filename='"+o.filename+"',\n"
            if o.representation is not None:
                setupStr+="         object_name ='"+o.representation+"',\n"
                setupStr+="         object_filename ='"+o.representation_file+"'\n"
            setupStr+="         )\n"
            setupStr+="h1.addCompartment("+o.name+")\n"
            rs = o.surfaceRecipe
            if rs :
                setupStr+=o.name+"_surface = Recipe(name='"+o.name+"_surf')\n"
                for ingr in rs.ingredients:                
                    if useXref :
                        io_ingr.write(ingr,pathout+os.sep+ingr.name,ingr_format="python")
                        setupStr+="execfile('"+pathout+os.sep+ingr.name+".py',globals(),{'recipe':"+o.name+"_surface})\n"
                    else :
                        ingrnode = io_ingr.ingrPythonNode(ingr,recipe=o.name+"_surface")
                        setupStr+=ingrnode 
                setupStr+=o.name+".setSurfaceRecipe("+o.name+"_surface)\n"
            ri = o.innerRecipe
            if ri :
                setupStr+=o.name+"_inner = Recipe(name='"+o.name+"_int')\n"
                for ingr in rs.ingredients:                
                    if useXref :
                        io_ingr.write(ingr,pathout+os.sep+ingr.name,ingr_format="python")
                        setupStr+="execfile('"+pathout+os.sep+ingr.name+".py',globals(),{'recipe':"+o.name+"_inner})\n"
                    else :
                        ingrnode = io_ingr.ingrPythonNode(ingr,recipe=o.name+"_inner")
                        setupStr+=ingrnode   
                setupStr+=o.name+".setInnerRecipe("+o.name+"_inner)\n"
        setupStr+="afviewer.SetHistoVol(h1,0,display=False)\n"
        setupStr+="afviewer.displayPreFill()\n"
        setupStr+="bbox = afviewer.helper.getObject('histvolBB')\n"
        setupStr+="if bbox is None : bbox = afviewer.helper.box('histvolBB',cornerPoints=h1.boundingBox)\n"
        setupStr+="helper = afviewer.helper\n"
        setupStr+="noGUI = False\n"
        setupStr+="try :\n"
        setupStr+="    print ('try')\n"
        setupStr+="    AFGui.Set('"+self.name+"',helper=afviewer.helper,afviewer=afviewer,histoVol=h1,bbox=bbox)\n"
        setupStr+="except:\n"
        setupStr+="    print ('no GUI')\n"
        setupStr+="    noGUI = True\n"
        f = open(setupfile,"w")        
        f.write(setupStr)
        f.close()
            
    def includeIngrRecipes(self,ingrname, include):
        """
        Include or Exclude the given ingredient from the recipe. 
        (similar to an active state toggle)
        """
        r =  self.exteriorRecipe
        if (self.includeIngrRecipe(ingrname, include,r)) :
            return
        for o in self.compartments:
            rs =  o.surfaceRecipe
            if (self.includeIngrRecipe(ingrname, include,rs)) :
                return
            ri =  o.innerRecipe
            if (self.includeIngrRecipe(ingrname, include,ri)) :
                return
                
    def includeIngrRecipe(self,ingrname, include,rs):
        """
        Include or Exclude the given ingredient from the given recipe. 
        (similar to an active state toggle)
        """
        for ingr in rs.exclude :
            if ingr.name==ingrname:
                if not include :
                    return True
                else :
                    rs.addIngredient(ingr)
                    return True
        for ingr in rs.ingredients:
            if ingrname == ingr.name:
                if not include :
                    rs.delIngredients(ingr)
                    return True
                else :
                    return True

    def includeIngredientRecipe(self,ingr, include):
        """
        Include or Exclude the given ingredient from the recipe. 
        (similar to an active state toggle)
        """
        r=ingr.recipe#()
        if include :
            r.addIngredient(ingr)
        else :
            r.delIngredient(ingr)

    def sortIngredient(self,reset = False):
        # make sure all recipes are sorted from large to small radius
        if self.exteriorRecipe:
            self.exteriorRecipe.sort()
        for o in self.compartments:
            o.molecules = []
            if reset :
                o.reset()
            if o.innerRecipe:
                o.innerRecipe.sort()
            if o.surfaceRecipe:
                o.surfaceRecipe.sort()

    def setGradient(self,**kw):
        """
        create a grdaient
        assign weight to point
        listorganelle influenced
        listingredient influenced
        """
        if "name" not in kw :
            print ("name kw is required")
            return
        gradient = Gradient(**kw) 
        #default gradient 1-linear Decoy X
        self.gradients[kw["name"]] = gradient
        
    def setDefaultOptions(self):
        """reset all the options to their default values"""
        for options in self.OPTIONS:
            setattr(self,options,self.OPTIONS[options]["default"])
        
    def callFunction(self,function,args=[],kw={}):
        """
        helper function to callback another function with the give arguments and keywords. 
        Optionally time stamp it.
        """
        if self._timer:
            res = self.timeFunction(function,args,kw)
        else :
            if len(kw):
                res = function(*args,**kw)
            else :
                res = function(*args)
        return res
        
    def timeFunction(self,function,args,kw):
        """
        Mesure the time for performing the provided function. 
    
        @type  function: function
        @param function: the function to execute
        @type  args: liste
        @param args: the liste of arguments for the function
        
    
        @rtype:   list/array
        @return:  the center of mass of the coordinates
        """
         
        t1=time()
        if len(kw):
            res = function(*args,**kw)
        else :
            res = function(*args)
        print(("time "+function.__name__, time()-t1))
        return res

    def SetRBOptions(self,obj="moving",**kw):
        """
        Change the rigid body options
        """
        key = ["shape","child",
                    "dynamicsBody", "dynamicsLinearDamp", 
                    "dynamicsAngularDamp", 
                    "massClamp","rotMassClamp"]
        for k in key :
            val = kw.pop( k, None)
            if val is not None:  
                self.dynamicOptions[obj][k]=val
                
    def SetSpringOptions(self,**kw):
        """
        Change the spring options, mainly used by C4D.
        """        
        key = ["stifness","rlength","damping"]
        for k in key :
            val = kw.pop( k, None)
            if val is not None:
                self.springOptions[k] = val
        
    def setupRBOptions(self):
        """
        Set default value for rigid body options
        """
        self.springOptions["stifness"] = 1.
        self.springOptions["rlength"] = 0.
        self.springOptions["damping"] = 1.
        self.dynamicOptions["spring"]={}
        self.dynamicOptions["spring"]["child"] = True
        self.dynamicOptions["spring"]["shape"] = "auto"
        self.dynamicOptions["spring"]["dynamicsBody"] = "on"
        self.dynamicOptions["spring"]["dynamicsLinearDamp"] = 0.0
        self.dynamicOptions["spring"]["dynamicsAngularDamp"] = 0.0
        self.dynamicOptions["spring"]["massClamp"] = 1.
        self.dynamicOptions["spring"]["rotMassClamp"] = 1.        
        self.dynamicOptions["moving"]={}
        self.dynamicOptions["moving"]["child"] = True
        self.dynamicOptions["moving"]["shape"] = "auto"
        self.dynamicOptions["moving"]["dynamicsBody"] = "on"
        self.dynamicOptions["moving"]["dynamicsLinearDamp"] = 1.0
        self.dynamicOptions["moving"]["dynamicsAngularDamp"] = 1.0
        self.dynamicOptions["moving"]["massClamp"] = .001
        self.dynamicOptions["moving"]["rotMassClamp"] = .1        
        self.dynamicOptions["static"]={}
        self.dynamicOptions["static"]["child"] = True
        self.dynamicOptions["static"]["shape"] = "auto"
        self.dynamicOptions["static"]["dynamicsBody"] = "off"
        self.dynamicOptions["static"]["dynamicsLinearDamp"] = 0.0
        self.dynamicOptions["static"]["dynamicsAngularDamp"] = 0.0
        self.dynamicOptions["static"]["massClamp"] = 100.
        self.dynamicOptions["static"]["rotMassClamp"] = 1
        self.dynamicOptions["surface"]={}
        self.dynamicOptions["surface"]["child"] = True
        self.dynamicOptions["surface"]["shape"] = "auto"
        self.dynamicOptions["surface"]["dynamicsBody"] = "off"
        self.dynamicOptions["surface"]["dynamicsLinearDamp"] = 0.0
        self.dynamicOptions["surface"]["dynamicsAngularDamp"] = 0.0
        self.dynamicOptions["surface"]["massClamp"] = 100.
        self.dynamicOptions["surface"]["rotMassClamp"] = 1
                    
    def writeArraysToFile(self, f):
        """write self.gridPtId and self.distToClosestSurf to file. (pickle) """
        pickle.dump(self.grid.masterGridPositions, f)
        pickle.dump(self.grid.gridPtId, f)
        pickle.dump(self.grid.distToClosestSurf, f)

    def readArraysFromFile(self, f):
        """write self.gridPtId and self.distToClosestSurf to file. (pickle) """
        pos = pickle.load(f)
        self.grid.masterGridPositions   = pos     
        
        id = pickle.load(f)
        #assert len(id)==len(self.gridPtId)
        self.grid.gridPtId = id

        dist = pickle.load(f)
        #assert len(dist)==len(self.distToClosestSurf)
        self.grid.distToClosestSurf_store = self.grid.distToClosestSurf[:]
        if len(dist):
            self.grid.distToClosestSurf = dist#grid+organelle+surf       
        self.grid.freePoints = list(range(len(id)))
        
    def saveGridToFile(self,gridFileOut):
        """
        Save the current grid and the compartment grid information in a file. (pickle) 
        """
        d = os.path.dirname(gridFileOut)
        if not os.path.exists(d):
            print ("gridfilename path problem",gridFileOut)
            return
        f = open(gridFileOut, 'wb')#'w'
        self.writeArraysToFile(f) #save self.gridPtId and self.distToClosestSurf
        
        for compartment in self.compartments:
            compartment.saveGridToFile(f)
        f.close()

    def restoreGridFromFile(self, gridFileName):
        """
        Read and setup the grid from the given filename. (pickle) 
        """
        from bhtree import bhtreelib
        aInteriorGrids = []
        aSurfaceGrids = []
        f = open(gridFileName,'rb')
        self.readArraysFromFile(f) #read gridPtId and distToClosestSurf
        #self.BuildGrids()        
        for compartment in self.compartments:
            compartment.readGridFromFile(f) 
            aInteriorGrids.append(compartment.insidePoints)
            aSurfaceGrids.append(compartment.surfacePoints)
            compartment.OGsrfPtsBht = bhtreelib.BHtree(tuple(compartment.vertices), None, 10)
        f.close()
        self.grid.aInteriorGrids = aInteriorGrids
        self.grid.aSurfaceGrids = aSurfaceGrids
 
    
    
    def setMinMaxProteinSize(self):
        """
        Retrieve and store mini and maxi ingredient size
        """
        for compartment in self.compartments:
            mini, maxi = compartment.getMinMaxProteinSize()
            if mini < self.smallestProteinSize:
                self.computeGridParams = True
                self.smallestProteinSize = mini

            if maxi > self.largestProteinSize:
                self.computeGridParams = True
                self.largestProteinSize = maxi

        if self.exteriorRecipe:
            smallest, largest = self.exteriorRecipe.getMinMaxProteinSize()

            if smallest < self.smallestProteinSize:
                self.smallestProteinSize = smallest

            if largest > self.largestProteinSize:
                self.largestProteinSize = largest

    def extractMeshComponent(self, obj):
        """
        Require host helper. Return the v,f,n of the given object
        """
        print ("extractMeshComponent",helper.getType(obj))
        if helper is None : 
            print ("no Helper found")            
            return None,None,None
        if helper.getType(obj) == helper.EMPTY: #compartment master parent?
            childs = helper.getChilds(obj)
            for ch in childs:
                name = helper.getName(ch)
#                print ("childs ",name)
                if helper.getType(ch) == helper.EMPTY:
                    c = helper.getChilds(ch)
                    #should be all polygon
                    faces=[]
                    vertices=[]
                    vnormals=[]
                    for pc in c :
                        f,v,vn = helper.DecomposeMesh(pc,edit=False,
                                                copy=False,tri=True,transform=True)
                        faces.extend(f)
                        vertices.extend(v)
                        vnormals.extend(vn)
                    return vertices, faces, vnormals
                elif helper.getType(ch) == helper.POLYGON :
                    faces,vertices,vnormals = helper.DecomposeMesh(ch,
                                    edit=False,copy=False,tri=True,transform=True)
                    return vertices, faces, vnormals
                else :
                    continue
        elif helper.getType(obj) == helper.POLYGON :
            name = helper.getName(obj)
            #helper.triangulate(c4dorganlle)
#            print ("polygon ",name)
            faces,vertices,vnormals = helper.DecomposeMesh(obj,
                                    edit=False,copy=False,tri=True,transform=True)
#            print ("returning v,f,n",len(vertices))
            return vertices, faces, vnormals
        else :
            print ("extractMeshComponent",helper.getType(obj),helper.POLYGON,
                   helper.getType(obj)==helper.POLYGON)
            return None,None,None
            
    def setCompartmentMesh(self, compartment, ref_obj):
        """
        Require host helper. Change the mesh of the given compartment and recompute
        inside and surface point.
        """
        if compartment.ref_obj == ref_obj : return
        if os.path.isfile(ref_obj):
            fileName, fileExtension = os.path.splitext(ref_obj)            
            if helper is not None:#neeed the helper
#                print ("read withHelper")
                helper.read(ref_obj)
                geom = helper.getObject(fileName)
                #reparent to the fill parent
                #rotate ?
                if helper.host != "c4d" and geom is not None:
                    #need to rotate the transform that carry the shape
                    helper.rotateObj(geom,[0.,-math.pi/2.0,0.0])
        else :
            geom = helper.getObject(ref_obj)
        if geom is not None :
            vertices, faces, vnormals = self.extractMeshComponent(geom)
            compartment.setMesh(filename=ref_obj,vertices=vertices, faces=faces, vnormals=vnormals )
                
    def addCompartment(self, compartment):
        """
        Add the given compartment to the environment. Extend the main bounding box if need
        """
        compartment.setNumber(self.nbCompartments)
        self.nbCompartments += 1

        fits, bb = compartment.inBox(self.boundingBox)
        
        if not fits:
            self.boundingBox = bb
        CompartmentList.addCompartment(self, compartment)

    def getPointCompartmentId(self,point,ray=3):
        #check if point inside  of the compartments
        #closest grid point is 
        d,pid=self.grid.getClosestGridPoint(point)
        cid = self.grid.gridPtId[pid]
        return cid
        ncomp = len(self.compartments)
        if ncomp:
            comp = ncomp
            for i in range(ncomp):                
                inside = self.compartments[comp-1].checkPointInside_rapid(point,self.grid.diag,ray=ray)
                if inside :
                    return -(comp)
                comp=comp-1
            #the point is not inside , is it on the surface ? ie distance to surface < X?
            for i in range(ncomp):                
                distance,nb = self.compartments[i].OGsrfPtsBht.query(point)
                if distance < 1.0 :
                    return i+1
        return 0

    def longestIngrdientName(self):
        """
        Helper function for gui. Return the size of the longest ingredient name
        """        
        M=20        
        r = self.exteriorRecipe
        if r :
            for ingr in r.ingredients:
                if len(ingr.name) > M :
                    M = len(ingr.name)
        for o in self.compartments:
            rs = o.surfaceRecipe
            if rs :
                for ingr in rs.ingredients:
                    if len(ingr.name) > M :
                        M = len(ingr.name)  
            ri = o.innerRecipe
            if ri :
                for ingr in ri.ingredients:
                    if len(ingr.name) > M :
                        M = len(ingr.name)
        return M

    def loopThroughIngr(self,cb_function):
        """
        Helper function that loop through all ingredients of all recipe and apply the give 
        callback functio on each ingredients.
        """
        r = self.exteriorRecipe
        if r :
            for ingr in r.ingredients:
                cb_function(ingr)
        for o in self.compartments:
            rs = o.surfaceRecipe
            if rs :
                for ingr in rs.ingredients:
                    cb_function(ingr)  
            ri = o.innerRecipe
            if ri :
                for ingr in ri.ingredients:
                    cb_function(ingr)
                    
    def getIngrFromNameInRecipe(self,name,r):
        """
        Given an ingredient name and a recipe, retrieve the ingredient object instance
        """
        if r :
            #check if name start with comp name
            #backward compatibility
            #print name,r.name
            if name.find(r.name) == -1 :
                name = r.name+"__"+name
            for ingr in r.ingredients:
                #print ingr.name
                if name == ingr.name :
                    return ingr
            for ingr in r.exclude:
                if name == ingr.name :
                    return ingr        
        return None
        
    def getIngrFromName(self,name,compNum=None):
        """
        Given an ingredient name and optionally the compartment number, retrieve the ingredient object instance
        """
        if compNum == None :
            r = self.exteriorRecipe
            ingr = self.getIngrFromNameInRecipe(name,r)
            if ingr is not None : return ingr
            for o in self.compartments :            
                rs = o.surfaceRecipe
                ingr = self.getIngrFromNameInRecipe(name,rs)
                if ingr is not None : return ingr
                ri = o.innerRecipe
                ingr = self.getIngrFromNameInRecipe(name,ri)
                if ingr is not None : return ingr            
        elif compNum == 0 :
            r = self.exteriorRecipe
            ingr = self.getIngrFromNameInRecipe(name,r)
            if ingr is not None : return ingr
            else : return None
        elif compNum > 0 :
            o=self.compartments[compNum-1]
            rs = o.surfaceRecipe
            ingr = self.getIngrFromNameInRecipe(name,rs)
            if ingr is not None : return ingr
            else : return None
        else : #<0
            o=self.compartments[(compNum*-1)-1]
            ri = o.innerRecipe
            ingr = self.getIngrFromNameInRecipe(name,ri)
            if ingr is not None : return ingr
            else : return None
            
    def setExteriorRecipe(self, recipe):
        """
        Set the exterior recipe with the given one. Create the weakref.
        """
        assert isinstance(recipe, Recipe)
        self.exteriorRecipe = recipe
        recipe.compartment = self#weakref.ref(self)
        for ingr in recipe.ingredients:
            ingr.compNum = 0

    def BuildCompartmentsGrids(self):  
        """
        Build the comparmtents grid (intrior and surface points) to be merged with the main grid
        """        
        aInteriorGrids = []
        aSurfaceGrids = []

        for compartment in self.compartments:
            if autopack.verbose :
                print("in HistoVol, compartment.isOrthogonalBoudingBox =", 
                      compartment.isOrthogonalBoudingBox)
            a,b = compartment.BuildGrid(self)
            aInteriorGrids.append(a)
            aSurfaceGrids.append(b)
    
        self.grid.aInteriorGrids = aInteriorGrids
        self.grid.aSurfaceGrids = aSurfaceGrids
        if autopack.verbose:
            print("I'm out of the loop and have build my grid with inside points")
            print ("build Grids",self.innerGridMethod,
                   len(self.grid.aSurfaceGrids))
 

    def buildGrid(self, boundingBox=None, gridFileIn=None, rebuild=True,
                  gridFileOut=None, previousFill=False,previousfreePoint=None):
        """
        The main build grid function. Setup the main grid and merge the 
        compartment grid. The setup is de novo or using previously builded grid 
        or restored using given file.
        """ 
        if self.use_halton:
            from autopack.Grid import HaltonGrid as Grid
        else :
            from autopack.Grid import Grid            
        #check viewer, and setup the progress bar               
        self.reportprogress(label="Building the Master Grid")
        #get and test the bounding box 
        if boundingBox is None:
            boundingBox = self.boundingBox
        else:
            assert len(boundingBox)==2
            assert len(boundingBox[0])==3
            assert len(boundingBox[1])==3
        self.sortIngredient(reset=rebuild)
        self.reportprogress(label="Computing the number of grid points")
        if rebuild or gridFileIn is not None or self.grid is None:
            # save bb for current fill
            self.fillBB = boundingBox
            self.grid = Grid(boundingBox=boundingBox,
                               space=self.smallestProteinSize)
            nbPoints = self.grid.gridVolume
            if autopack.verbose : 
                print ("new Grid with  ",boundingBox,self.grid.gridVolume)       
            if rebuild :
                verts = []            
                for orga in self.compartments:
                    if len(orga.vertices):
                        for pt3d in orga.vertices:
                            verts.append( pt3d )
                self.grid.set_surfPtsBht(verts)
                self.grid.distToClosestSurf_store = self.grid.distToClosestSurf[:] 
                nbPoints = self.grid.gridVolume
        else :
            self.grid.reset()
            nbPoints = len(self.grid.freePoints)
        if autopack.verbose : 
            print("$$$$$$$$  gridVolume = nbPoints = ", nbPoints, 
              " grid.nbGridPoints = ", self.grid.nbGridPoints)
        #self.grid.distToClosestSurf_store = self.grid.distToClosestSurf[:] 

        if gridFileIn is not None :#and not rebuild:
            if autopack.verbose : 
                print ("file in for building grid but it doesnt work well")
            self.grid.filename = gridFileIn
            if self.nFill == 0 :#first fill, after we can just reset
                print ("restore from file")
                self.restoreGridFromFile(gridFileIn)
        elif gridFileIn is None and rebuild:
            # assign ids to grid points
            if autopack.verbose : 
                print ("file is None thus re/building grid distance")
            self.BuildCompartmentsGrids()
        else :
            print ("file is not rebuild nor restore from file")
        if gridFileOut is not None:
            self.saveGridToFile(gridFileOut)
            self.grid.filename = gridFileOut
        self.exteriorVolume = self.grid.computeExteriorVolume(compartments=self.compartments,space=self.smallestProteinSize,fbox_bb=self.fbox_bb)
        r = self.exteriorRecipe
        if r:
            r.setCount(self.exteriorVolume)#should actually use the fillBB
            
        if not rebuild :
            self.grid.distToClosestSurf = self.grid.distToClosestSurf_store[:]   
        else :
            self.grid.distToClosestSurf_store = self.grid.distToClosestSurf[:]   

        if self.use_gradient and len(self.gradients) and rebuild :
            for g in self.gradients:
                self.gradients[g].buildWeigthMap(boundingBox,self.grid.masterGridPositions)
        if previousFill:
            distance = self.grid.distToClosestSurf#[:]
            nbFreePoints = nbPoints#-1              #Graham turned this off on 5/16/12 to match August Repair for May Hybrid
            for i,mingrs in enumerate(self.molecules) :#( jtrans, rotMatj, self, ptInd )
                nbFreePoints=self.onePrevIngredient(i,mingrs,distance,nbFreePoints,self.molecules)
            for organelle in self.compartments:
                for i,mingrs in enumerate(organelle.molecules) :#( jtrans, rotMatj, self, ptInd )
                    nbFreePoints=self.onePrevIngredient(i,mingrs,distance,nbFreePoints,organelle.molecules)
            self.grid.nbFreePoints = nbFreePoints
        self.setCompatibility()

    def BuildGrids(self):  
        """
        Build the comparmtents grid (intrior and surface points) to be merged with the main grid
        Note :
        #New version allows for orthogonal box to be used as an organelle requireing no expensive InsidePoints test
        # FIXME make recursive?
        """        
        aInteriorGrids = []
        aSurfaceGrids = []

        for compartment in self.compartments:
            print("in HistoVol, compartment.isOrthogonalBoudingBox =", compartment.isOrthogonalBoudingBox)
            b = []
            if compartment.isOrthogonalBoudingBox==1:
                self.EnviroOnly = True
                print(">>>>>>>>>>>>>>>>>>>>>>>>> Not building a grid because I'm an Orthogonal Bounding Box")
                a =self.grid.getPointsInCube(compartment.bb, None, None) #This is the highspeed shortcut for inside points! and no surface! that gets used if the fillSelection is an orthogonal box and there are no other compartments.
                self.grid.gridPtId[a] = -compartment.number
                compartment.surfacePointsCoords = None
                bb0x, bb0y, bb0z = compartment.bb[0]
                bb1x, bb1y, bb1z = compartment.bb[1]
                AreaXplane = (bb1y-bb0y)*(bb1z-bb0z)
                AreaYplane = (bb1x-bb0x)*(bb1z-bb0z)
                AreaZplane = (bb1y-bb0y)*(bb1x-bb0x)
                vSurfaceArea = abs(AreaXplane)*2+abs(AreaYplane)*2+abs(AreaZplane)*2
                print("vSurfaceArea = ", vSurfaceArea)
                compartment.insidePoints = a
                compartment.surfacePoints = b
                compartment.surfacePointsCoords = []
                compartment.surfacePointsNormals = []
                print(' %d inside pts, %d tot grid pts, %d master grid'%( len(a),len(a), len(self.grid.masterGridPositions)))
                compartment.computeVolumeAndSetNbMol(self, b, a,areas=vSurfaceArea)               
                #print("I've built a grid in the compartment test with no surface", a)
                print("The size of the grid I build = ", len(a))

            if self.innerGridMethod =="sdf" and compartment.isOrthogonalBoudingBox!=1: # A fillSelection can now be a mesh too... it can use either of these methods
                a, b = compartment.BuildGrid_utsdf(self) # to make the outer most selection from the master and then the compartment
            elif self.innerGridMethod == "bhtree" and compartment.isOrthogonalBoudingBox!=1:  # surfaces and interiors will be subtracted from it as normal!
                a, b = compartment.BuildGrid(self)
            elif self.innerGridMethod == "jordan" and compartment.isOrthogonalBoudingBox!=1:  # surfaces and interiors will be subtracted from it as normal!
                a, b = compartment.BuildGrid_jordan(self)
            elif self.innerGridMethod == "jordan3" and compartment.isOrthogonalBoudingBox!=1:  # surfaces and interiors will be subtracted from it as normal!
                a, b = compartment.BuildGrid_jordan(self,ray=3)
            elif self.innerGridMethod == "pyray" and compartment.isOrthogonalBoudingBox!=1:  # surfaces and interiors will be subtracted from it as normal!
                a, b = compartment.BuildGrid_pyray(self)    
            aInteriorGrids.append(a)
            print("I'm ruther in the loop")
            aSurfaceGrids.append(b)
    
        self.grid.aInteriorGrids = aInteriorGrids
        print("I'm out of the loop and have build my grid with inside points")
        self.grid.aSurfaceGrids = aSurfaceGrids
        print ("build Grids",self.innerGridMethod,len(self.grid.aSurfaceGrids))
       
    def buildGridOld(self, boundingBox=None, gridFileIn=None, rebuild=True,
                  gridFileOut=None, previousFill=False,previousfreePoint=None):
        """
        The main build grid function. Setup the main grid and merge the 
        compartment grid. The setup is de novo or using previously builded grid 
        or restored using given file. This funcion should be
        split in smaller function for clarity.
        """                
        if self.afviewer is not None and hasattr(self.afviewer,"vi"):
            self.afviewer.vi.progressBar(label="Building the Master Grid")
        if boundingBox is None:
            boundingBox = self.boundingBox
        else:
            assert len(boundingBox)==2
            assert len(boundingBox[0])==3
            assert len(boundingBox[1])==3
        # make sure all recipes are sorted from large to small radius
        if self.exteriorRecipe:
            self.exteriorRecipe.sort()

        for o in self.compartments:
            o.molecules = []
            if rebuild :
                o.reset()
            if o.innerRecipe:
                o.innerRecipe.sort()
            if o.surfaceRecipe:
                o.surfaceRecipe.sort()

        if self.afviewer is not None and hasattr(self.afviewer,"vi"):
            self.afviewer.vi.progressBar(label="Computing the number of grid points")
        if rebuild or gridFileIn is not None:
            # save bb for current fill
            self.fillBB = boundingBox
            grid = Grid()
            self.grid = grid
            grid.boundingBox = boundingBox
            # compute grid spacing
            grid.gridSpacing = space = self.smallestProteinSize*1.1547  # 2/sqrt(3)
            print ("$$$$$$$$  ",boundingBox,space,self.smallestProteinSize)
            grid.gridVolume,grid.nbGridPoints = self.callFunction(grid.computeGridNumberOfPoint,(boundingBox,space))
        grid =self.grid
        nbPoints = self.grid.gridVolume
        print("$$$$$$$$  gridVolume = nbPoints = ", nbPoints, " grid.nbGridPoints = ", self.grid.nbGridPoints)
        # compute 3D point coordiantes for all grid points
        if rebuild or gridFileIn is not None:
            self.callFunction(grid.create3DPointLookup) #generate grid.masterGridPositions 
#            print('grid size', grid.nbGridPoints)
            grid.nbSurfacePoints = 0
            #self.isFree = numpy.ones( (nbPoints,), 'i') # Will never shrink
#            print('nb freePoints', nbPoints)#-1)
            # Id is set set to None initially
            grid.gridPtId = numpy.zeros(nbPoints)#[0]*nbPoints

        xl,yl,zl = boundingBox[0]
        xr,yr,zr = boundingBox[1]
        # distToClosestSurf is set to self.diag initially
        self.grid.diag = diag = vlen( vdiff((xr,yr,zr), (xl,yl,zl) ) )
        if rebuild or gridFileIn is not None:
            self.grid.distToClosestSurf = [diag]*nbPoints#surface point too?
            self.grid.distToClosestSurf = numpy.array(self.grid.distToClosestSurf)
            self.grid.freePoints = list(range(nbPoints))
        else :
            #just reset
            self.grid.distToClosestSurf = [diag]*len(self.grid.distToClosestSurf)#surface point too?
            self.grid.distToClosestSurf = numpy.array(self.grid.distToClosestSurf)
            self.grid.freePoints = list(range(len(self.grid.freePoints)))
            nbPoints = len(self.grid.freePoints)
#        print 'DIAG', diag,self.grid.distToClosestSurf
        self.grid.distToClosestSurf_store = self.grid.distToClosestSurf[:] 
#        if gridFileIn is None :
#            gridFileIn = self.grid_filename
#        if gridFileOut is None :
#            gridFileOut= self.grid_filename
#        if self.grid_filename is not None and not os.path.isfile(self.grid_filename):
#            gridFileIn = None

#        if rebuild :
            #this restore/store the grid information of the organelle.
        if gridFileIn is not None :#and not rebuild:
            print ("file in for building grid but it doesnt work well")
            self.grid.filename = gridFileIn
            if self.nFill == 0 :#?:
                print ("restore from file")
                self.restoreGridFromFile(gridFileIn)
        elif gridFileIn is None and rebuild:
            # assign ids to grid points
            print ("file is None for building grid")
            self.BuildGrids()
        else :
            print ("file is not rebuild")
        if gridFileOut is not None:
            self.saveGridToFile(gridFileOut)
            self.grid.filename = gridFileOut
        # get new set of freePoints which includes surface points
 #       nbPoints = nbPoints-1            #Graham Turned off this redundant nbPoints-1 call on 8/27/11
 #       nbPoints = nbPoints-1          #Graham Turned this one off on 5/16/12 to match August repair in Hybrid
        grid.nbFreePoints = nbPoints#-1
        grdPts = grid.masterGridPositions
        grid.nbFreePoints = len(grdPts)
        # build BHTree for surface points (off grid)
        if rebuild :
            verts = []            
            for orga in self.compartments:
                if orga.surfacePointsCoords:
                    for pt3d in orga.surfacePointsCoords:
                        verts.append( pt3d )
    
            from bhtree import bhtreelib
            grid.surfPtsBht = None
            if verts :
               grid.surfPtsBht = bhtreelib.BHtree( verts, None, 10)
           
        # build list of compartments without a recipe#????
        noRecipe = []
        if self.exteriorRecipe is None:
            noRecipe.append( 0 )
        for o in self.compartments:
            if o.surfaceRecipe is None:
                noRecipe.append( o.number )
            if o.innerRecipe is None:
                noRecipe.append( -o.number )

        # compute exterior volume 
        unitVol = grid.gridSpacing**3
        totalVolume = grid.gridVolume*unitVol
        if self.fbox_bb is not None :
                V,nbG = self.callFunction(grid.computeGridNumberOfPoint,(self.fbox_bb,space))
                totalVolume = V*unitVol
        for o in self.compartments:
            #totalVolume -= o.surfaceVolume
            totalVolume -= o.interiorVolume
        self.exteriorVolume = totalVolume

        r = self.exteriorRecipe
        if r:
            r.setCount(totalVolume)#should actually use the fillBB
            
        if self.use_gradient and len(self.gradients) and rebuild :
            for g in self.gradients:
                self.gradients[g].buildWeigthMap(boundingBox,grid.masterGridPositions)
        if not rebuild :
            self.grid.distToClosestSurf = self.grid.distToClosestSurf_store[:]   
        else :
            self.grid.distToClosestSurf_store = self.grid.distToClosestSurf[:]   
            
        #we should be able here to update the number of free point using a previous grid 
        #overlap
#        print("previousFill",previousFill)
        if previousFill:#actually if there is a previous fill
            #get the intersecting point and update freePoints from this one if they are not free
            #previousfreePoint
            #compute the intersection bounding box and get ptindice for both grid 
            #by getPointsInCube
            #check which one are in freePoints from previous, and update the current one
            #update the curentpass
#            #how to update the distance for each prest ingr ?
            distance = self.grid.distToClosestSurf#[:]
#            nbFreePoints = nbPoints-1                #This already comes from the Point Volume- no subtraction needed (Graham turned off on 8/27/11)
            nbFreePoints = nbPoints#-1              #Graham turned this off on 5/16/12 to match August Repair for May Hybrid
#            backupmol = self.molecules[:]          
#            molecules=self.molecules
##            print("molecules",molecules)
#            for organelle in self.organelles:
#                molecules.extend(organelle.molecules)
            for i,mingrs in enumerate(self.molecules) :#( jtrans, rotMatj, self, ptInd )
                nbFreePoints=self.onePrevIngredient(i,mingrs,distance,nbFreePoints,self.molecules)
            for organelle in self.compartments:
                for i,mingrs in enumerate(organelle.molecules) :#( jtrans, rotMatj, self, ptInd )
                    nbFreePoints=self.onePrevIngredient(i,mingrs,distance,nbFreePoints,organelle.molecules)

#                jtrans, rotMatj, ingr, ptInd = mingrs
##                print ("OK",jtrans, rotMatj, ingr, ptInd)
#                centT = ingr.transformPoints(jtrans, rotMatj, ingr.positions[-1])
#                insidePoints = {}
#                newDistPoints = {}
#                mr = self.get_dpad(ingr.compNum)
#                spacing = self.smallestProteinSize
#                jitter = ingr.getMaxJitter(spacing)
#                dpad = ingr.minRadius + mr + jitter
#                insidePoints,newDistPoints = ingr.getInsidePoints(self.grid,
#                                    self.grid.masterGridPositions,dpad,distance,
#                                    centT=centT,
#                                    jtrans=jtrans, 
#                                    rotMatj=rotMatj)
#                # update free points
#                if len(insidePoints) and self.placeMethod.find("panda") != -1:
#                       print (ingr.name," is inside")
#                       self.checkPtIndIngr(ingr,insidePoints,i,ptInd)
#                       #ingr.inside_current_grid = True
#                else  :
#                    #not in the grid
#                    print (ingr.name," is outside")
#                    #rbnode = ingr.rbnode[ptInd]
#                    #ingr.rbnode.pop(ptInd)
#                    molecules[i][3]=-1
#                    #ingr.rbnode[-1] = rbnode
#                #(self, histoVol,insidePoints, newDistPoints, freePoints,
#                #        nbFreePoints, distance, masterGridPositions, verbose)
#                nbFreePoints = ingr.updateDistances(self,insidePoints, newDistPoints, 
#                            self.grid.freePoints, nbFreePoints, distance, 
#                            self.grid.masterGridPositions,0)
            self.grid.nbFreePoints = nbFreePoints
        #self.hgrid.append(self.grid)
        self.setCompatibility()

    def onePrevIngredient(self,i,mingrs,distance,nbFreePoints,marray):
        """
        Unused
        """
        jtrans, rotMatj, ingr, ptInd = mingrs
        centT = ingr.transformPoints(jtrans, rotMatj, ingr.positions[-1])
        insidePoints = {}
        newDistPoints = {}
        mr = self.get_dpad(ingr.compNum)
        spacing = self.smallestProteinSize
        jitter = ingr.getMaxJitter(spacing)
        dpad = ingr.minRadius + mr + jitter
        insidePoints,newDistPoints = ingr.getInsidePoints(self.grid,
                            self.grid.masterGridPositions,dpad,distance,
                            centT=centT,
                            jtrans=jtrans, 
                            rotMatj=rotMatj)
        # update free points
        if len(insidePoints) and self.placeMethod.find("panda") != -1:
               print (ingr.name," is inside")
               self.checkPtIndIngr(ingr,insidePoints,i,ptInd,marray)
               #ingr.inside_current_grid = True
        else  :
            #not in the grid
            print (ingr.name," is outside")
            #rbnode = ingr.rbnode[ptInd]
            #ingr.rbnode.pop(ptInd)
            marray[i][3]=-ptInd#uniq Id ?
            #ingr.rbnode[-1] = rbnode
        #(self, histoVol,insidePoints, newDistPoints, freePoints,
        #        nbFreePoints, distance, masterGridPositions, verbose)
        #doesnt seem to work properly...
        nbFreePoints = ingr.updateDistances(self,insidePoints, newDistPoints, 
                    self.grid.freePoints, nbFreePoints, distance, 
                    self.grid.masterGridPositions,0)
        #should we reset the ingredient ? completion ?
        if not ingr.is_previous:
            ingr.firstTimeUpdate = True
            ingr.counter = 0
            ingr.rejectionCounter = 0
            ingr.completion= 0.0            #should actually count it
            if hasattr(ingr,"allIngrPts"):  #Graham here on 5/16/12 are these two lines safe?
                del ingr.allIngrPts         #Graham here on 5/16/12 are these two lines safe?        
        return nbFreePoints
        

    def checkPtIndIngr(self,ingr,insidePoints,i,ptInd,marray):
        """
        We need to check if the point indice is correct in the case of panda packing.
        as the pt indice in the result array have a different meaning.
        """
        #change key for rbnode too
        rbnode = None
        if ptInd in ingr.rbnode:
          rbnode = ingr.rbnode[ptInd]
          ingr.rbnode.pop(ptInd)
        elif -ptInd in ingr.rbnode:
          rbnode = ingr.rbnode[-ptInd]
          ingr.rbnode.pop(-ptInd)
        else :
            print ("ptInd "+str(ptInd)+" not in ingr.rbnode")
        if i < len(marray):
            marray[i][3]=insidePoints.keys()[0]
            ingr.rbnode[insidePoints.keys()[0]] = rbnode
#        else :
#            nmol = len(self.molecules)
#            for j,organelle in enumerate(self.organelles):
#                print (i,nmol+len(organelle.molecules))
#                if i < nmol+len(organelle.molecules):
#                    organelle.molecules[i-nmol][3]=insidePoints.keys()[0]
#                    ingr.rbnode[insidePoints.keys()[0]] = rbnode
#                else :
#                    nmol+=len(organelle.molecules)
            
    def setCompatibility(self):
        """
        in earlier version the grid was part of the environment class. 
        Since we split the grid in her own class, to avoid some error during the transition 
        we alias all the function and attribute.
        """
        self.getPointsInCube = self.grid.getPointsInCube
        self.boundingBox=self.grid.boundingBox
        self.gridPtId = self.grid.gridPtId
        self.freePoints = self.grid.freePoints
        self.diag=self.grid.diag
        self.gridSpacing = self.grid.gridSpacing
        self.nbGridPoints = self.grid.nbGridPoints
        self.nbSurfacePoints = self.grid.nbSurfacePoints
        self.gridVolume = self.grid.gridVolume # will be the toatl number of grid points
        self.masterGridPositions = self.grid.masterGridPositions
        self.aInteriorGrids = self.grid.aInteriorGrids
        self.aSurfaceGrids = self.grid.aSurfaceGrids
        self.surfPtsBht=self.grid.surfPtsBht
        self.gridPtId = self.grid.gridPtId = numpy.array(self.grid.gridPtId,int)
        
    def getSortedActiveIngredients(self, allIngredients, verbose=0):
        """
        Sort the active ingredient according their pirority and radius.
        # first get the ones with a packing priority
        # Graham- This now works in concert with ingredient picking
        
        # Graham here- In the new setup, priority is infinite with abs[priority] increasing (+)
        # An ingredients with (-) priority will pack from greatest abs[-priority] one at a time
        #     to lease abs[-priority]... each ingredient will attempt to reach its molarity
        #     before moving on to the next ingredient, and all (-) ingredients will try to
        #     deposit before other ingredients are tested.
        # An ingredient with (+) priority will recieve a weighted value based on its abs[priority]
        #     e.g. an ingredient with a priority=10 will be 10x more likely to be picked than
        #     an ingredient with a priority=1.
        # An ingredient with the default priority=0 will recieve a weighted value based on its
        #     complexity. (currently complexity = minRadius), thus a more 'complex' ingredient
        #     will more likely try to pack before a less 'complex' ingredient.
        #     IMPORTANT: the +priority list does not fully mix with the priority=0 list, but this
        #     should be an option... currently, the priority=0 list is normalized against a range
        #     up to the smallest +priority ingredient and appended to the (+) list
        # TODO: Add an option to allow + ingredients to be weighted by assigned priority AND complexity
        #     Add an option to allow priority=0 ingredients to fit into the (+) ingredient list
        #       rather than appending to the end.
        #     Even better, add an option to set the max priority for the 0list and then plug the results
        #       into the (+) ingredient list.
        #     Get rid of the (-), 0, (+) system and recreate this as a new flag and a class function
        #        so we can add multiple styles of sorting and weighting systems.
        #     Make normalizedPriorities and thresholdPriorities members of Ingredient class to avoid
        #        building these arrays.
        """   
        ingr1 = []  # given priorities
        priorities1 = []
        ingr2 = []  # priority = 0 or none and will be assigned based on complexity
        priorities2 = []
        ingr0 = []  # negative values will pack first in order of abs[packingPriority]
        priorities0 = []
        for ing in allIngredients:
            if ing.completion >= 1.0: continue # ignore completed ingredients
            if ing.packingPriority is None or ing.packingPriority == 0 :
                ingr2.append(ing)
                priorities2.append(ing.packingPriority)
            elif ing.packingPriority > 0 :
                ingr1.append(ing)
                priorities1.append(ing.packingPriority)
            else:
                #ing.packingPriority    = -ing.packingPriority    
                ingr0.append(ing)
                priorities0.append(ing.packingPriority)

        if self.pickWeightedIngr:            
    #Graham here on 5/16/12. Double check that this new version is correct- it uses a very different function than the working September version from 2011
            #sorted(ingr1, key=attrgetter('packingPriority', 'minRadius','completion'), reverse=True)
            #ingr1.sort(key=ingredient_compare1)  #Fails 5/21/12
            if sys.version < "3.0.0" :
                ingr1.sort(ingredient_compare1)
                # sort ingredients with no priority based on radius and completion
                #sorted(ingr2, key=attrgetter('packingPriority', 'minRadius','completion'), reverse=True)
                #ingr2.sort(key=ingredient_compare2)  #Fails 5/21/12
                ingr2.sort(ingredient_compare2)
                #sorted(ingr0, key=attrgetter('packingPriority', 'minRadius','completion'), reverse=True)
                #ingr0.sort(key=ingredient_compare0)
                ingr0.sort(ingredient_compare0)  #Fails 5/21/12
    #            ingr0.sort()
            else :
                try :
                    ingr1.sort(key=ingredient_compare1)
                    ingr2.sort(key=ingredient_compare2)
                    ingr0.sort(key=ingredient_compare0)
                except :
                    print ("ATTENTION INGR NOT SORTED")
        #for ing in ingr3 : ing.packingPriority    = -ing.packingPriority
        #GrahamAdded this stuff in summer 2011, beware!
        if len(ingr1) != 0:
            lowestIng = ingr1[len(ingr1)-1]
            self.lowestPriority = lowestIng.packingPriority
        else :
            self.lowestPriority = 1.
        if verbose:
            print('self.lowestPriority for Ing1 = ', self.lowestPriority)
        self.totalRadii = 0
        for radii in ingr2:
            if radii.modelType=='Cylinders':
                r = max(radii.length/2.,radii.minRadius)
            elif radii.modelType=='Spheres':
                r = radii.minRadius
            elif radii.modelType=='Cube':
                r = radii.minRadius
            self.totalRadii = self.totalRadii + r
            if verbose:
                print('self.totalRadii += ', r, "=",self.totalRadii)
            if r==0 : 
                print (radii,radii.name)
                #safety 
                self.totalRadii = self.totalRadii + 1.0
            
        self.normalizedPriorities0 = []
        for priors2 in ingr2:
            if priors2.modelType=='Cylinders':
                r = max(priors2.length/2.,priors2.minRadius)
            elif priors2.modelType=='Spheres':
                r = priors2.minRadius            
            np = float(r)/float(self.totalRadii) * self.lowestPriority
            self.normalizedPriorities0.append(np)
            priors2.packingPriority = np
            if verbose:
                print('self.normalizedPriorities0 = ', self.normalizedPriorities0)
        activeIngr0 = ingr0#+ingr1+ingr2  #cropped to 0 on 7/20/10
        
        if verbose:
            print('len(activeIngr0)', len(activeIngr0))
        activeIngr12 = ingr1+ingr2
        if verbose:
            print('len(activeIngr12)', len(activeIngr12))
        packingPriorities = priorities0+priorities1+priorities2
        if verbose:
            print('priorities0 is ', priorities0)
            print('priorities1 is ', priorities1)
            print('priorities2 is ', priorities2)
            print('packingPriorities', packingPriorities)

#        if verbose>0:
#            print 'Ingredients:'
#            for i, ingr in enumerate(activeIngr):
#                packPri = ingr.packingPriority
#                if packPri is None:
#                    packPri = -1
#                print '  comp:%2d #:%3d pri:%3d compl:%.2f mRad:%5.1f t:%4d n:%4d '%(
#                    ingr.compNum, i, packPri, ingr.completion, ingr.minRadius,
#                    ingr.nbMol, ingr.counter)
#            raw_input("hit enter")

        return activeIngr0, activeIngr12

#    import fill3isolated # Graham cut the outdated fill3 from this document and put it in a separate file. turn on here if you want to use it.
            
    def updateIngr(self,ingr,completion=0.0,nbMol=0,counter=0):
        """helper function for updating the ingredient completion, nbmol and counter """
        ingr.counter = counter
        ingr.nbMol = nbMol
        ingr.completion = completion

    def reset(self):
        """Reset everything to empty and not done"""
        self.fbox_bb = None
        self.totnbJitter = 0
        self.jitterLength = 0.0
        r =  self.exteriorRecipe
        self.resetIngrRecip(r)
        self.molecules=[]
        for orga in self.compartments:
            #orga.reset()
            rs =  orga.surfaceRecipe
            self.resetIngrRecip(rs)
            ri =  orga.innerRecipe
            self.resetIngrRecip(ri)

        self.ingr_result = {}
        if self.world is not None :
            #need to clear all node
            nodes = self.rb_panda[:]
            for node in nodes:
                self.delRB(node)
            self.static = []
            self.moving = None
        if self.octree is not None :
            del self.octree
            self.octree = None
        #the reset doesnt touch the grid...
#        del self.close_ingr_bhtree
        from bhtree import bhtreelib    
#        if len(self.rTrans) :
#            try :
#                bhtreelib.freeBHtree(self.close_ingr_bhtree)
#            except :
#                print ("problem freeBHtree")
        self.rTrans=[]   
        self.rIngr=[] 
        self.rRot=[] 
        self.result = []
        #rapid node ?
        
    def resetIngrRecip(self,recip):
        """Reset all ingredient of the given recipe"""
        if recip:            
            for ingr in recip.ingredients:
                ingr.results=[]
                ingr.firstTimeUpdate = True
                ingr.counter = 0
                ingr.rejectionCounter = 0
                ingr.completion= 0.0
                ingr.prev_alt = None
                ingr.start_positions=[]
                if hasattr(ingr,"allIngrPts"):  #Graham here on 5/16/12 are these two lines safe?
                    del ingr.allIngrPts         #Graham here on 5/16/12 are these two lines safe?
                if hasattr(ingr,"isph"):  
                    ingr.isph = None
                if hasattr(ingr,"icyl"):  
                    ingr.icyl = None
                if hasattr(ingr,"allIngrPts"):
                    delattr(ingr, "allIngrPts")
            for ingr in recip.exclude:
                ingr.start_positions=[]
                ingr.prev_alt = None
                ingr.results=[]
                ingr.firstTimeUpdate = True
                ingr.counter = 0
                ingr.rejectionCounter = 0
                ingr.completion= 0.0
                if hasattr(ingr,"allIngrPts"):  #Graham here on 5/16/12 are these two lines safe?
                    del ingr.allIngrPts           
                if hasattr(ingr,"isph"):  
                    ingr.isph = None    
                if hasattr(ingr,"icyl"):  
                    ingr.icyl = None
                if hasattr(ingr,"allIngrPts"):
                    delattr(ingr, "allIngrPts")
                    
    def resetIngr(self,ingr):
        """Reset the given ingredient (count, completion, nmol)"""
        ingr.counter = 0
        ingr.nbMol = 0
        ingr.completion = 0.0

    def getActiveIng(self):
        """Return all remaining active ingredients"""
        allIngredients = []
        r = self.exteriorRecipe
        if r is not None:
            if not hasattr(r,"molecules") :
                r.molecules = []
        if r:
            for ingr in r.ingredients:
                ingr.counter = 0 # counter of placed molecules
                if  ingr.nbMol > 0:
                    ingr.completion = 0.0
                    allIngredients.append(ingr)
                else:
                    ingr.completion = 1.0
            
        for o in self.compartments:
            if not hasattr(o,"molecules") :
                o.molecules = []
            r = o.surfaceRecipe
            if r:
                for ingr in r.ingredients:
                    ingr.counter = 0 # counter of placed molecules
                    if  ingr.nbMol > 0:
                        ingr.completion = 0.0
                        allIngredients.append(ingr)
                    else:
                        ingr.completion = 1.0

            r = o.innerRecipe
            if r:
                for ingr in r.ingredients:
                    ingr.counter = 0 # counter of placed molecules
#                    print "nbMol",ingr.nbMol
                    if  ingr.nbMol > 0:
                        ingr.completion = 0.0
                        allIngredients.append(ingr)
                    else:
                        ingr.completion = 1.0
        return allIngredients

    def pickIngredient(self,vThreshStart, verbose=0):
        """
        Main function that decide the next ingredient the packing will try to 
        drop. The picking is weighted or random
        """
        if self.pickWeightedIngr : 
            if self.thresholdPriorities[0] == 2 :
                # Graham here: Walk through -priorities first
                ingr = self.activeIngr[0]
            else:
                #prob = uniform(vRangeStart,1.0)  #Graham 9/21/11 This is wrong...vRangeStart is the point index, need active list i.e. thresholdPriority to be limited
                prob = uniform(0,1.0)
                ingrInd = 0
                for threshProb in self.thresholdPriorities:
                    if prob <= threshProb:
                        break
                    ingrInd = ingrInd + 1
                if ingrInd <  len(self.activeIngr):
                    ingr = self.activeIngr[ingrInd]
                else :
                    print("error in histoVol pick Ingredient",ingrInd)
                    ingr = self.activeIngr[0]
                if verbose:
                    print ('weighted',prob, vThreshStart, ingrInd,ingr.name)
        else :
            #if verbose:
            #    print "random in activeIngr"
            r = random()#randint(0, len(self.activeIngr)-1)#random()
            #n=int(r*(len(self.activeIngr)-1))
            n=int(r*len(self.activeIngr))
            ingr = self.activeIngr[n]
#            print (r,n,ingr.name,len(self.activeIngr)) #Graham turned back on 5/16/12, but may be costly
        return ingr

    def get_dpad(self,compNum):
        """Return the largest encapsulatingRadius and use it for padding"""
        mr = 0.0
        if compNum==0: # cytoplasm -> use cyto and all surfaces
            for ingr1 in self.activeIngr:
                if ingr1.compNum>=0:
                    r = ingr1.encapsulatingRadius
                    if r>mr:
                        mr = r
        else:
            for ingr1 in self.activeIngr:
                if ingr1.compNum==compNum or ingr1.compNum==-compNum:
                    r = ingr1.encapsulatingRadius
                    if r>mr:
                        mr = r
        return mr

    def checkIfUpdate(self,ingr,nbFreePoints,verbose=False):
        """Check if we need to update the distance array. Part of the hack free points"""
        if hasattr(ingr,"nbPts"):
            if hasattr(ingr,"firstTimeUpdate") and not ingr.firstTimeUpdate:
                ratio = float(ingr.nbPts)/float(nbFreePoints)
#                print('freePtsUpdateThrehod = ', self.freePtsUpdateThrehod)
                if verbose:
                    print("checkIfUpdate: ratio = ",ratio,"nbFreePoints = ", nbFreePoints, "ingr.nbPts = ",ingr.nbPts)
                if ratio > self.freePtsUpdateThrehod :
                    return True
                else :
                    if ingr.haveBeenRejected and ingr.rejectionCounter > 5:
                        ingr.haveBeenRejected = False
                        return True
                    #do we check to total freepts? or crowded state ?
                    else :
                        return False
            else :
                ingr.firstTimeUpdate = False
                return True
        else :
            return True


    def getPointToDrop(self,ingr,radius,jitter,freePoints,nbFreePoints,distance,
                       compId,compNum,vRangeStart,vThreshStart, verbose=False):
        """
        Decide next point to use for dropping a given ingredent. The picking can be 
        random, based on closest distance, based on gradients, ordered.
        This function also update the available free point except when hack is on.
        """
        verbose = True
        if ingr.packingMode=='close':
            t1 = time()
            allIngrPts = []
            allIngrDist = []
            if ingr.modelType=='Cylinders' and ingr.useLength :
                cut = ingr.length-jitter
#            if ingr.modelType=='Cube' : #radius iactually the size
#                cut = min(self.radii[0]/2.)-jitter
            elif ingr.cutoff_boundary is not None :
                cut  = radius+ingr.cutoff_boundary-jitter
            else :
                cut  = radius-jitter
            for pt in freePoints:#[:nbFreePoints]:
                d = distance[pt]#look up the distance
                if compId[pt]==compNum and d>=cut:
                    allIngrPts.append(pt)
                    allIngrDist.append(d)
            if verbose:
                print("time to filter using for loop ", time()-t1)
        else:
            t1 = time()
            allIngrPts = []
#            print("allIngrPts = ", allIngrPts)
#            print("len (allIngrPts) = ", len(allIngrPts))
            if ingr.modelType=='Cylinders' and ingr.useLength :
                cut = ingr.length-jitter
            elif ingr.cutoff_boundary is not None :
                cut  = radius+ingr.cutoff_boundary-jitter                
            else :
                cut  = radius-jitter
            #for pt in freePoints[:nbFreePoints]:
            if hasattr(ingr,"allIngrPts") and self._hackFreepts:
                allIngrPts = ingr.allIngrPts
#                print("hasattr(ingr,allIngrPts)")
#                print ("Running nofreepoint HACK")
            else :
                #use periodic update according size ration grid
                update = self.checkIfUpdate(ingr,nbFreePoints)
#                print("in update ",update)
#                print "update ", update,nbFreePoints,hasattr(ingr,"allIngrPts"),cut
                if update :
#                    print("in update loop")
                    for i in range(nbFreePoints):
#                        print("in i range of update loop",i,freePoints[i],distance[i])
                        pt = freePoints[i]
                        d = distance[pt]
#                        print("in update for/if")
#                        print pt,compId[pt], d,cut,compNum,(compId[pt]==compNum and d>=cut)
                        if compId[pt]==compNum and d>=cut:
                            allIngrPts.append(pt)
                    #allIngrDist.append(d)
                    ingr.allIngrPts = allIngrPts
                    ingr.cut = cut
#                    if verbose:
#                    print("getPointToDrop len(allIngrPts) = ", len(allIngrPts))
                else :
#                    print ("else")
                    if hasattr(ingr,"allIngrPts"):
#                        print("allIngrPts = ingr.allIngrPts two elses deep")
                        allIngrPts = ingr.allIngrPts
                    else :    #Graham Here on 5/16/12, double check that this is safe as its not in the September version
                        allIngrPts = freePoints[:nbFreePoints]
                        ingr.allIngrPts = allIngrPts
#                        print("in the last else")
#                        print('freepoint routine here may be unsafe, not in old version, so doublecheck')
                        #compltly unsafe due to surface points !!!
#        print(("time to filter ",nbFreePoints," using lambda ", time()-t1))
        # no point left capable of accomodating this ingredient
#        print("allIngrPts = ", allIngrPts)
        print("len (allIngrPts) = ", len(allIngrPts))
        if len(allIngrPts)==0:
            t=time()
            ingr.completion = 1.0
            ind = self.activeIngr.index(ingr)
            #if ind == 0:
            vRangeStart = vRangeStart + self.normalizedPriorities[0]
            if ind > 0:
                #j = 0
                for j in range(ind):                
                    self.thresholdPriorities[j] = self.thresholdPriorities[j] + self.normalizedPriorities[ind]
            self.activeIngr.pop(ind)
# Start of massive overruling section from corrected thesis file of Sept. 25, 2012
            #this function also depend on the ingr.completiion that can be restored ?
            self.activeIngr0, self.activeIngr12 = self.callFunction(self.getSortedActiveIngredients, (self.activeIngr,False))
            if verbose:
                print('No point left for ingredient %s %f minRad %.2f jitter %.3f in component %d'%(
                ingr.name, ingr.molarity, radius, jitter, compNum))
                print ('len(allIngredients', len(self.activeIngr))
                print ('len(self.activeIngr0)', len(self.activeIngr0))
                print ('len(self.activeIngr12)', len(self.activeIngr12))
            self.activeIngre_saved = self.activeIngr[:]

            self.totalPriorities = 0 # 0.00001
            for priors in self.activeIngr12:
                pp = priors.packingPriority
                self.totalPriorities = self.totalPriorities + pp
                if verbose :
                    print ('totalPriorities = ', self.totalPriorities)
            previousThresh = 0
            self.normalizedPriorities = []
            self.thresholdPriorities = [] 
            # Graham- Once negatives are used, if picked random# 
            # is below a number in this list, that item becomes 
            #the active ingredient in the while loop below
            for priors in self.activeIngr0:
                self.normalizedPriorities.append(0)
                if self.pickWeightedIngr :
                    self.thresholdPriorities.append(2)
            for priors in self.activeIngr12:
                #pp1 = 0
                pp = priors.packingPriority
                if self.totalPriorities != 0:
                    np = float(pp)/float(self.totalPriorities)
                else:
                    np=0.
                self.normalizedPriorities.append(np)
                if verbose :
                    print ('np is ', np, ' pp is ', pp, ' tp is ', np + previousThresh)
                self.thresholdPriorities.append(np + previousThresh)
                previousThresh = np + float(previousThresh)
            self.activeIngr = self.activeIngr0 + self.activeIngr12
            
#            nls=0
#            totalNumMols = 0
#            for threshProb in self.thresholdPriorities:
#                nameMe = self.activeIngr[nls]
#                if verbose:
#                    print ('threshprop Get Point is %f for ingredient: %s %s %d'%(threshProb, nameMe,nameMe.name,nameMe.nbMol))
#                totalNumMols += nameMe.nbMol
#                if verbose:
#                    print ('totalNumMols Get Point= ', totalNumMols)
#                nls+=1

            #print 'vThreshStart before = ', vThreshStart
            #vThreshStart = self.thresholdPriorities[0]
            #print 'vThreshStart after = ', vThreshStart
            #print 'because vself.thresholdPriorities[0] = ', self.thresholdPriorities[0]

            #self.thresholdPriorities.pop(ind)
            #self.normalizedPriorities.pop(ind)
            if verbose:
                print ("time to reject the picking", time()-t)
# End of massive overruling section from corrected thesis file of Sept. 25, 2011
# this chunk overwrites the next three lines from July version. July 5, 2012
#            self.thresholdPriorities.pop(ind)                    
#            self.normalizedPriorities.pop(ind)
#            print(("time to reject the picking", time()-t))
            return False,vRangeStart
#        ptInd = allIngrPts[0]       #turned this off when I imported the large overrulling 
# code from Sept 25 2011 thesis version on July 5, 2012
        if self.pickRandPt:
            t2=time()
            if ingr.packingMode=='close':
                order = numpy.argsort(allIngrDist)
                # pick point with closest distance
                ptInd = allIngrPts[order[0]]
                if (ingr.rejectionCounter % 300 == 0):
                    ptIndr = int(random()*len(allIngrPts))
                    ptInd = allIngrPts[ptIndr]                            
            elif ingr.packingMode=='gradient' and self.use_gradient:  
                #get the most probable point using the gradient                
                #use the gradient weighted map and get mot probabl point
                print ("pick point from gradients",(len(allIngrPts)))
                ptInd = self.gradients[ingr.gradient].pickPoint(allIngrPts) 
            else:
                # pick a point randomly among free points
                ptIndr = int(random()*len(allIngrPts))
                ptInd = allIngrPts[ptIndr]            
            if ptInd is None :
                t=time()
                if verbose:
                    print('No point left for ingredient %s %f minRad %.2f jitter %.3f in component %d'%(
                    ingr.name, ingr.molarity, radius, jitter, compNum))
                ingr.completion = 1.0
                ind = self.activeIngr.index(ingr)
                #if ind == 0:
                vRangeStart = vRangeStart + self.normalizedPriorities[0]
                if ind > 0:
                    #j = 0
                    for j in range(ind):                
                        self.thresholdPriorities[j] = self.thresholdPriorities[j] + self.normalizedPriorities[ind]
                self.activeIngr.pop(ind)
                if verbose:
                    print('popping this gradient ingredient array must be redone using Sept 25, 2011 thesis version as above for nongraient ingredients, TODO: July 5, 2012')
                self.thresholdPriorities.pop(ind)
                self.normalizedPriorities.pop(ind)
                if verbose:
                        print(("time to reject the picking", time()-t))
                print(("vRangeStart",vRangeStart))
                return False,vRangeStart                    

#            print(("time to random pick a point", time()-t2))
        else :
            t3=time()
            allIngrPts.sort()
            ptInd = allIngrPts[0]
#            print(("time to sort and pick a point", time()-t3))
        return True,ptInd

#    import fill4isolated # Graham cut the outdated fill4 from this document and put it in a separate file. turn on here if you want to use it.
    def removeOnePoint(self, pt,freePoints,nbFreePoints):
        try :
                # New system replaced by Graham on Aug 18, 2012
                nbFreePoints -= 1  
                vKill = freePoints[pt]
                vLastFree = freePoints[nbFreePoints]
                freePoints[vKill] = vLastFree
                freePoints[vLastFree] = vKill
                # End New replaced by Graham on Aug 18, 2012
        except:
                pass
        return nbFreePoints

    def getTotalNbObject(self,allIngredients):
        totalNbIngr = 0
        for ingr in allIngredients:
            if ingr.Type == "Grow" :
                totalNbIngr += int (ingr.nbMol*(ingr.length/ingr.uLength))
            else :    
                totalNbIngr += ingr.nbMol
        return totalNbIngr
    
    def fill5(self, seedNum=14, stepByStep=False, verbose=False, sphGeom=None,
              labDistGeom=None, debugFunc=None,name = None, vTestid = 3,vAnalysis = 0,**kw):
        """
        the latest packing loop 
        ## Fill the grid by picking an ingredient first and then
        ## this packing should be able to continue from a previous one
        ## find a suitable point using the ingredient's placer object
        """
        import time
        t1=time.time()
        self.timeUpDistLoopTotal = 0 #Graham added to try to make universal "global variable Verbose" on Aug 28
        self.static=[]
        if self.grid is None:
            print("no grid setup")
            return
        # create a list of active ingredients indices in all recipes to allow
        # removing inactive ingredients when molarity is reached
        allIngredients = self.callFunction(self.getActiveIng)
        usePP = False
        if "usePP" in kw :
            usePP = kw["usePP"]
        nbIngredients = len(allIngredients)
        self.cFill = self.nFill
        if name == None :
            name = "F"+str(self.nFill)
        self.FillName.append(name)
        self.nFill+=1
        # seed random number generator
        SEED=seedNum
        numpy.random.seed(SEED)#for gradient
        seed(seedNum)
        self.randomRot.setSeed(seed=seedNum)
        # create copies of the distance array as they change when molecules
        # are added, theses array can be restored/saved before feeling
        freePoints = self.grid.freePoints[:]
        nbFreePoints = len(freePoints)#-1
#        self.freePointMask = numpy.ones(nbFreePoints,dtype="int32")
        if "fbox" in kw :  # Oct 20, 2012  This is part of the code that is breaking the grids for all meshless compartment fills
            self.fbox = kw["fbox"]
        if self.fbox is not None and not self.EnviroOnly :
            self.freePointMask = numpy.ones(nbFreePoints,dtype="int32")
            bb_insidepoint = self.grid.getPointsInCube(self.fbox, [0,0,0], 1.0)[:]#center and radius ?3,runTime=self.runTimeDisplay
            self.freePointMask[bb_insidepoint]=0
            bb_outside = numpy.nonzero(self.freePointMask)
            self.grid.gridPtId[bb_outside] = 99999
        compId = self.grid.gridPtId
        #why a copy? --> can we split ?
        distance = self.grid.distToClosestSurf[:]

        spacing = self.smallestProteinSize

        # DEBUG stuff, should be removed later
        self.jitterVectors = []
        self.jitterLength = 0.0
        self.totnbJitter = 0
        self.maxColl = 0.0
        self.successfullJitter = []
        self.failedJitter = []
        
        #this function also depend on the ingr.completiion that can be restored ?
        self.activeIngr0, self.activeIngr12 = self.callFunction(self.getSortedActiveIngredients, (allIngredients,verbose))

        print('len(allIngredients', len(allIngredients))
        print('len(self.activeIngr0)', len(self.activeIngr0))
        print('len(self.activeIngr12)', len(self.activeIngr12))
        self.activeIngre_saved = self.activeIngr[:]

        self.totalPriorities = 0 # 0.00001
        for priors in self.activeIngr12:
            pp = priors.packingPriority
            self.totalPriorities = self.totalPriorities + pp
            print('totalPriorities = ', self.totalPriorities)
        previousThresh = 0
        self.normalizedPriorities = []
        self.thresholdPriorities = [] 
        # Graham- Once negatives are used, if picked random# 
        # is below a number in this list, that item becomes 
        # the active ingredient in the while loop below
        for priors in self.activeIngr0:
            self.normalizedPriorities.append(0)
            if self.pickWeightedIngr :#why ?
                self.thresholdPriorities.append(2)
        for priors in self.activeIngr12:
            #pp1 = 0
            pp = priors.packingPriority
            if self.totalPriorities != 0:
                np = float(pp)/float(self.totalPriorities)
            else:
                np=0.
            self.normalizedPriorities.append(np)
            print('np is ', np, ' pp is ', pp, ' tp is ', np + previousThresh)
            self.thresholdPriorities.append(np + previousThresh)
            previousThresh = np + float(previousThresh)
        self.activeIngr = self.activeIngr0 + self.activeIngr12

        nls=0
        totalNumMols = 0
        self.totalNbIngr = self.getTotalNbObject(allIngredients)
        if len(self.thresholdPriorities ) == 0:
            for ingr in allIngredients:
                totalNumMols += ingr.nbMol
            print('totalNumMols Fill5if = ', totalNumMols)
        else :                
            for threshProb in self.thresholdPriorities:
                nameMe = self.activeIngr[nls]
                print('threshprop Fill5else is %f for ingredient: %s %s %d'%(threshProb, nameMe,nameMe.name,nameMe.nbMol))
                totalNumMols += nameMe.nbMol
                print('totalNumMols Fill5else = ', totalNumMols)
                nls+=1
            
        vRangeStart = 0.0
        tCancelPrev=time.time()
        test = True
        kk=0
        ptInd = 0

        PlacedMols = 0
        vThreshStart = 0.0   # Added back by Graham on July 5, 2012 from Sept 25, 2011 thesis version
        
        #if bullet build the organel rbnode
        if self.placeMethod == "pandaBullet":
            self.setupPanda()
            for o in self.compartments:
                if o.rbnode is None :
                    o.rbnode = self.addMeshRBOrganelle(o)
        if usePP :
            import pp
            self.grab_cb = GrabResult() 
            self.pp_server = pp.Server(ncpus=autopack.ncpus)
#==============================================================================
#         #the big loop
#==============================================================================
        while nbFreePoints:
            print (".........At start of while loop, with vRangeStart = ", vRangeStart)
#            for o in self.compartments:
#                print ("compartments = ", o.name)
#            print("freePoints = ", freePoints, "nbFreePoints = ", nbFreePoints)
            if verbose > 1:
                print('Points Remaining', nbFreePoints, len(freePoints))
                print('len(self.activeIngr)', len(self.activeIngr))                
            
            #breakin test
            if len(self.activeIngr)==0:
                print('broken by len****')
                if hasattr(self,"afviewer"):
                    if self.afviewer is not None and hasattr(self.afviewer,"vi"):
                        self.afviewer.vi.resetProgressBar()
                        self.afviewer.vi.progressBar(label="Filling Complete")       
                break
            if vRangeStart>1:
                print('broken by vRange and hence Done!!!****')
                break   
            if self.cancelDialog :
                tCancel = time.time()
                if tCancel-tCancelPrev > 10.:
                    cancel=self.displayCancelDialog()
                    if cancel:
                        print("canceled by user: we'll fill with current objects up to time", tCancel)
                        break
                    #if OK, do nothing, i.e., continue loop (but not the function continue)
                    tCancelPrev= time.time()
            ## pick an ingredient
            
            ingr =  self.callFunction(self.pickIngredient,(vThreshStart,))
            print("picked Ingr ",ingr.name)
            #if ingr.completion >= 1.0 or ingr.is_previous:
            #    continue
#            ingr =  self.callFunction(self.pickIngredient,(vRangeStart,))   # Replaced this with previous line from Sept 25, 2011 thesis version on July 5, 2012
            if hasattr(self,"afviewer"):
                # C4D safety check added by Graham on July 10, 2012 until we can fix the uPy status bar for C4D
                #if self.host == 'c4d':
                    try :
                        import c4d
        #               Start working chunk pasted in by Graham from Sept 2011 on 5/16/12- Resorting to this because it fixes the status bar- must be a uPy problem?
        #               -Ludo: we need to use the helpr for that
        #               -Graham: The helper progressBar is broken at least for C4D
        #                   The self.afviewer.vi.progressBar(progress=int(p),label=ingr.name) functions shows zero progress until everything is ended, so using C4D for now
                        p = float(PlacedMols)/float(totalNumMols)*100.
                        c4d.StatusSetBar(int(p))
                        c4d.StatusSetText(ingr.name+" "+str(ingr.completion))
        #                print("PlacedMols= ", PlacedMols, ", while totalNumMols= ", totalNumMols, " so % = ", p)
        #               End working chunk pasted in by Graham from Sept 2011 on 5/16/12.
#                        print ("using c4d override for Status Bar")
                    except :
        #               p = ((float(t)-float(len(self.activeIngr)))/float(t))*100.
                        p = ((float(PlacedMols))/float(totalNumMols))#*100.    #This code shows 100% of ingredients all the time
                        if self.afviewer is not None and hasattr(self.afviewer,"vi"):
                            self.afviewer.vi.progressBar(progress=int(p),label=ingr.name+" "+str(ingr.completion))
                            if self.afviewer.renderDistance:
                                self.afviewer.vi.displayParticleVolumeDistance(distance,self)
                        #pass
                # End C4D safety check for Status bar added July 10, 2012
            compNum = ingr.compNum
            radius = ingr.minRadius
            jitter = self.callFunction(ingr.getMaxJitter,(spacing,))

            # compute dpad which is the distance at which we need to update
            # distances after the drop is successfull
            mr = self.get_dpad(compNum)
            dpad = ingr.minRadius + mr + jitter

            if verbose > 2:
                print('picked Ingr radius compNum dpad',radius,compNum,dpad)
            
            ## find the points that can be used for this ingredients
            ##
            res=self.callFunction(self.getPointToDrop,(ingr,radius,jitter,
                                        freePoints,nbFreePoints,
                                        distance,compId,compNum,vRangeStart,vThreshStart))
#                                        distance,compId,compNum,vRangeStart))   # Replaced this with Sept 25, 2011 thesis version on July 5, 2012
            if autopack.verbose : 
                print ("get drop point res",res)
            if res[0] :
                ptInd = res[1]
                if ptInd > len(distance):
                    print ("problem ",ptInd)
                    continue
            else :
                print ("vRangeStart coninue ",res)
                vRangeStart = res[1]
                continue
            print ("picked ",ptInd)
            #place the ingrediant
            if self.overwritePlaceMethod :
                ingr.placeType = self.placeMethod
            #check the largestProteinSize
            if ingr.encapsulatingRadius > self.largestProteinSize : 
                self.largestProteinSize = ingr.encapsulatingRadius 
            #histoVol, ptInd, freePoints, nbFreePoints, distance, dpad,usePP,
            #  stepByStep=False, verbose=False,
            success, nbFreePoints = self.callFunction(ingr.place,(self, ptInd, 
                                freePoints, nbFreePoints, distance, dpad,usePP,
                                stepByStep, verbose),
                                {"debugFunc":debugFunc})
#            print("nbFreePoints after PLACE ",nbFreePoints)
            if success:
                print ("success",ingr.completion)
                #update largest protein size
                #problem when the encapsulatingRadius is actually wrong
                if ingr.encapsulatingRadius > self.largestProteinSize : 
                    self.largestProteinSize = ingr.encapsulatingRadius
                PlacedMols+=1
#                nbFreePoints=self.removeOnePoint(ptInd,freePoints,nbFreePoints)  #Hidden by Graham on March 1, 2013 until we can test.
            else :
                print ("rejected",ingr.rejectionCounter)
            if ingr.completion >= 1.0 :
                print('completed***************', ingr.name)
                print('PlacedMols = ', PlacedMols)
                ind = self.activeIngr.index(ingr)
                print('activeIngr index of ',ingr.name,ind)
                print('threshold p len ',len(self.thresholdPriorities),
                                  len(self.normalizedPriorities)) 
#                vRangeStart = vRangeStart + self.normalizedPriorities[0]
                if ind > 0:
                    #j = 0
                    for j in range(ind):   
                        if j >= len(self.thresholdPriorities) or j >= len(self.normalizedPriorities):
                            continue
                        self.thresholdPriorities[j] = self.thresholdPriorities[j] + self.normalizedPriorities[ind]
                self.activeIngr.pop(ind)
#                self.thresholdPriorities.pop(ind)  # Replaced these from July SVN version with large chunk of code from thesis Sept 25, 2011 version on July 5, 2012
#                self.normalizedPriorities.pop(ind) # Replaced these from July SVN version with large chunk of code on next lines from thesis Sept 25, 2011 version on July 5, 2012                
# BEGIN large chunk of code from proper thesis Sept 25, 2011 version to correctly replace simple pop above on July 5, 2012
                #this function also depend on the ingr.completiion that can be restored ?
                self.activeIngr0, self.activeIngr12 = self.callFunction(self.getSortedActiveIngredients, (self.activeIngr,verbose))
                if verbose > 2:
                    print ('len(self.activeIngr', len(self.activeIngr))
                    print ('len(self.activeIngr0)', len(self.activeIngr0))
                    print ('len(self.activeIngr12)', len(self.activeIngr12))
                self.activeIngre_saved = self.activeIngr[:]

                self.totalPriorities = 0 # 0.00001
                for priors in self.activeIngr12:
                    pp = priors.packingPriority
                    self.totalPriorities = self.totalPriorities + pp
#                    print ('totalPriorities = ', self.totalPriorities)
                previousThresh = 0
                self.normalizedPriorities = []
                self.thresholdPriorities = [] 
                # Graham- Once negatives are used, if picked random# 
                # is below a number in this list, that item becomes 
                #the active ingredient in the while loop below
                for priors in self.activeIngr0:
                    self.normalizedPriorities.append(0)
                    if self.pickWeightedIngr :
                        self.thresholdPriorities.append(2)
                for priors in self.activeIngr12:
                    #pp1 = 0
                    pp = priors.packingPriority
                    if self.totalPriorities != 0:
                        np = float(pp)/float(self.totalPriorities)
                    else:
                        np=0.
                    self.normalizedPriorities.append(np)
#                    print ('np is ', np, ' pp is ', pp, ' tp is ', np + previousThresh)
                    self.thresholdPriorities.append(np + previousThresh)
                    previousThresh = np + float(previousThresh)
                self.activeIngr = self.activeIngr0 + self.activeIngr12
                
#                nls=0
#                totalNumMols = 0
#                for threshProb in self.thresholdPriorities:
#                    nameMe = self.activeIngr[nls]
#                    print ('threshprop Success is %f for ingredient: %s %s %d'%(threshProb, nameMe,nameMe.name,nameMe.nbMol))
#                    totalNumMols += nameMe.nbMol
#                    print ('totalNumMols Success= ', totalNumMols)
#                    nls+=1
            
                #print 'vThreshStart before = ', vThreshStart
                #vThreshStart = self.thresholdPriorities[0]
                #print 'vThreshStart after = ', vThreshStart
                #print 'because vself.thresholdPriorities[0] = ', self.thresholdPriorities[0]
                #self.thresholdPriorities.pop(ind)
                #self.normalizedPriorities.pop(ind)
# END large chunk of code from proper thesis Sept 25, 2011 version to correctly replace simple pop above on July 5, 2012
#            if nbFreePoints == 0 :
#                break
        #0.8938
        # for debugging purposes
        #whats the difference with distancetosurface which is stored
        self.distancesAfterFill = distance
        self.freePointsAfterFill = freePoints
        self.nbFreePointsAfterFill = nbFreePoints
        self.distanceAfterFill = distance
        #self.rejectionCount = rejectionCount
#        c4d.documents.RunAnimation(doc, True)
        t2 = time.time()
        print('time to fill', t2-t1)
            
        if self.saveResult:
            self.grid.freePoints = freePoints[:]
            self.grid.distToClosestSurf = distance[:]
            #shoul check extension filename for type of saved file
            self.saveGridToFile(self.resultfile+"grid")
            self.grid.result_filename = self.resultfile+"grid"
            self.store()
            self.store_asTxt()
            self.store_asJson()            
            #self.saveGridToFile_asTxt(self.resultfile+"grid")freePointsAfterFill
            #should we save to text as well
            print('time to save in fil5', time.time()-t2)
#            vAnalysis = 0
            if vAnalysis == 1 :
    #    START Analysis Tools: Graham added back this big chunk of code for analysis tools and graphic on 5/16/12 Needs to be cleaned up into a function and proper uPy code            
                unitVol = self.grid.gridSpacing**3
                #totalVolume = self.grid.gridVolume*unitVol
                wrkDirRes= self.resultfile+"_analyze_"
                print('TODO: overwrite wrkDirRes with specific user directory for each run or each script or put in a cache and offer a chance to save it')
                print("self.compartments = ", self.compartments)
                for o in self.compartments: #only for compartment ?
                    #totalVolume -= o.surfaceVolume
                    #totalVolume -= o.interiorVolume
                    innerPointNum = len(o.insidePoints)-1
                    print ('  .  .  .  . ')
                    print ('for compartment o = ', o.name)
                    print ('inner Point Count = ', innerPointNum)
                    print ('inner Volume = ', o.interiorVolume)
                    print ('innerVolume temp Confirm = ', innerPointNum*unitVol)
                    usedPts = 0
                    unUsedPts = 0
                    #fpts = self.freePointsAfterFill
                    vDistanceString = ""
                    insidepointindce = numpy.nonzero(numpy.equal(self.grid.gridPtId,-o.number))[0]
                    for i in insidepointindce:#xrange(innerPointNum):
#                        pt = o.insidePoints[i] #fpts[i]
#                        print (pt,type(pt))
                        #for pt in self.histo.freePointsAfterFill:#[:self.histo.nbFreePointsAfterFill]:
                        d = self.distancesAfterFill[i]
                        vDistanceString += str(d)+"\n"
                        if d <= 0 :  #>self.smallestProteinSize-0.001:
                            usedPts += 1
                        else:
                            unUsedPts +=1
                    filename = wrkDirRes+"vResultMatrix1" + o.name + "_Testid" + str(vTestid) + "_Seed" + str(seedNum) + "_dists.txt" # Used this from thesis to overwrite less informative SVN version on next line on July 5, 2012
        #            filename = wrkDirRes+"/vDistances1.txt"
                    f = open(filename,"w")
                    vMyString = "I am on" + "\nThis is a new line."
                    f.write(vDistanceString)
                    f.close()
                    
                    #result is [pos,rot,ingr.name,ingr.compNum,ptInd]
                    #if resultfilename == None:
                    #resultfilename = self.resultfile
                    resultfilenameT = wrkDirRes+"vResultMatrix1" + o.name + "_Testid" + str(vTestid) + "_Seed" + str(seedNum) + "_Trans.txt" # Used this from thesis to overwrite less informative SVN version on next line on July 5, 2012
                    resultfilenameR = wrkDirRes+"vResultMatrix1" + o.name + "_Testid" + str(vTestid) + "_Seed" + str(seedNum) + "_Rot.txt" # Used this from thesis to overwrite less informative SVN version on next line on July 5, 2012
        #            resultfilenameT = wrkDirRes+"/vResultMatrix1" + o.name + "_Trans.txt"
        #            resultfilenameR = wrkDirRes+"/vResultMatrix1" + o.name + "_Rot.txt"
                    #pickle.dump(self.molecules, rfile)
                    #OR 
                    vTranslationString = ""
                    vRotationString = ""
                    result=[]
                    matCount = 0
                    # Add safety check for C4D until we can get uPy working for this matrix to hbp rotation function?
        #            from c4d import utils   # Removed by Graham on July 10, 2012 because replaced with more recent Thesis code on July 5, 2012 below
                    #what do you save everthing inleft hand ? and you actually dont use it ??
                    # Note July 4, 2012: the results are saved as right handed (see 2, 1, 0 for h, p, b) and used for analysis tools
                    # Note July 5, 2012: I found the better version we made and added it below to override the C4D version!
                    for pos, rot, ingr, ptInd in o.molecules:
                        #vMatrixString += str(result([pos]))+"\n"
        # BEGIN: newer code from Theis version added July 5, 2012
                        if hasattr(self,"afviewer"):        
                            mat = rot.copy()
                            mat[:3, 3] = pos
                            import math
                            from ePMV import comput_util as c
                            r  = c.matrixToEuler(mat)
                            h1 = math.degrees(math.pi + r[0])
                            p1 = math.degrees(r[1])
                            b1 = math.degrees(-math.pi + r[2])
                            #angles[0] = 180.0+angles[0]
                            #angles[2] = 180.0-angles[2] 
                            #hmat = self.afviewer.vi.FromMat(mat,transpose=True)
                            #rot = utils.MatrixToHPB(hmat)
                            print ('rot from matrix = ', r,h1,p1,b1)
        # END: newer code from Theis version added July 5, 2012
                        result.append([pos,rot])
                        pt3d = result[matCount][0]
                        x, y, z = pt3d #  ADDDED this line back from newer code from Theis version added July 5, 2012
        # BEGIN: retired SVN version, retired July 5, 2012
        #                x, y, z = pt3d
        #                rot3d = result[matCount][1][2]
        #                h1 = rot3d[2]
        #                p1 = rot3d[1]
        #                b1 = rot3d[0]
        #                rot3d = result[matCount][1][1]
        #                h2 = rot3d[2]
        #                p2 = rot3d[1]
        #                b2 = rot3d[0]
        #                rot3d = result[matCount][1][0]
        #                h3 = rot3d[2]
        #                p3 = rot3d[1]
        #                b3 = rot3d[0]
        # can we test for C4D for these last 6 lines until we can get same functionality from uPy?
        #                off = c4d.Vector(0)
        #                vec = c4d.Matrix(off, c4d.Vector(h1, p1, b1), c4d.Vector(h2,p2,b2), c4d.Vector(h3,p3,b3) )
        #                print vec  
        #                #m = rot3d #obj.GetMg()
        #                rot = utils.MatrixToHPB(vec)
        #                print 'rot from matrix = ', rot
        # END: retired SVN version, retired July 5, 2012                
                        vTranslationString += str(x)+ ",\t" + str(y) + ",\t" + str(z) + "\n"
                        #vRotationString += str(rot3d) #str(h)+ ",\t" + str(p) + ",\t" + str(b) + "\n"
                        vRotationString += str(h1)+ ",\t" + str(p1) + ",\t" + str(b1) + ",\t" + ingr.name +"\n"  #  ADDDED this line back from newer code from Theis version added July 5, 2012 to replace next line from SVN
        #                vRotationString += str(h1)+ ",\t" + str(p1) + ",\t" + str(b1) + ingr.name +"\n"
                        #vRotationString += str( (result[matCount][1]).x )+"\n"
                        matCount += 1
                    
                    
                    #result.append([pos,rot,ingr.name,ingr.compNum,ptInd])
                    #d = self.distancesAfterFill[pt]
                    #vDistanceString += str(d)+"\n"
                    #pickle.dump(result, rfile)
                    rfile = open(resultfilenameT, 'w')
                    rfile.write( vTranslationString )
                    rfile.close()
                    
                    rfile = open(resultfilenameR, 'w')
                    rfile.write( vRotationString )
                    rfile.close()
                    print ('len(result) = ', len(result))
                    print ('len(self.molecules) = ', len(self.molecules))
                    ### Graham Note:  There is overused disk space- the rotation matrix is 4x4 with an offset of 0,0,0 and we have a separate translation vector in the results and molecules arrays.  Get rid of the translation vector and move it to the rotation matrix to save space... will that slow the time it takes to extract the vector from the matrix when we need to call it?       
                    print ('*************************************************** vDistance String Should be on')
                    print ('unitVolume2 = ', unitVol)
                    print ('Number of Points Unused = ', unUsedPts)
                    print ('Number of Points Used   = ', usedPts)
                    print ('Volume Used   = ', usedPts*unitVol)
                    print ('Volume Unused = ', unUsedPts*unitVol)
                    print ('vTestid = ', vTestid)
                    print ('self.nbGridPoints = ', self.nbGridPoints)
                    print ('self.gridVolume = ', self.gridVolume)
    #        self.exteriorVolume = totalVolume
                        
            print("self.compartments In HistoVol = ", len(self.compartments))
            if self.compartments == [] :
                #o = self.histoVol
#                o = self.exteriorRecipe
                unitVol = self.grid.gridSpacing**3
                innerPointNum = len(freePoints)
                print ('  .  .  .  . ')
#                print ('for compartment o = ', o.name)
                print ('inner Point Count = ', innerPointNum)
#                print ('inner Volume = ', o.interiorVolume)
                print ('innerVolume temp Confirm = ', innerPointNum*unitVol)
                usedPts = 0
                unUsedPts = 0
                #fpts = self.freePointsAfterFill
                vDistanceString = ""
                for i in xrange(innerPointNum):
                    pt = freePoints[i] #fpts[i]
                    #for pt in self.histo.freePointsAfterFill:#[:self.histo.nbFreePointsAfterFill]:
                    d = self.distancesAfterFill[pt]
                    vDistanceString += str(d)+"\n"
                    if d <= 0 :  #>self.smallestProteinSize-0.001:
                        usedPts += 1
                    else:
                        unUsedPts +=1
#                filename = wrkDirRes+"/vResultMatrix1" + o.name + "_Testid" + str(vTestid) + "_Seed" + str(seedNum) + "_dists.txt" # Used this from thesis to overwrite less informative SVN version on next line on July 5, 2012
#                #            filename = wrkDirRes+"/vDistances1.txt"
#                f = open(filename,"w")
#                vMyString = "I am on" + "\nThis is a new line."
#                f.write(vDistanceString)
#                f.close()
#                resultfilenameT = wrkDirRes+"/vResultMatrix1" + o.name + "_Testid" + str(vTestid) + "_Seed" + str(seedNum) + "_Trans.txt" # Used this from thesis to overwrite less informative SVN version on next line on July 5, 2012
#                resultfilenameR = wrkDirRes+"/vResultMatrix1" + o.name + "_Testid" + str(vTestid) + "_Seed" + str(seedNum) + "_Rot.txt" # Used this from thesis to overwrite less informative SVN version on next line on July 5, 2012
#                vTranslationString = ""
#                vRotationString = ""
#                result=[]
#                matCount = 0
#                # Add safety check for C4D until we can get uPy working for this matrix to hbp rotation function?
#                #            from c4d import utils   # Removed by Graham on July 10, 2012 because replaced with more recent Thesis code on July 5, 2012 below
#                #what do you save everthing inleft hand ? and you actually dont use it ??
#                # Note July 4, 2012: the results are saved as right handed (see 2, 1, 0 for h, p, b) and used for analysis tools
#                # Note July 5, 2012: I found the better version we made and added it below to override the C4D version!
#                for pos, rot, ingr, ptInd in o.molecules:
#                    #vMatrixString += str(result([pos]))+"\n"
#                    # BEGIN: newer code from Theis version added July 5, 2012
#                    if hasattr(self,"afviewer"):        
#                        mat = rot.copy()
#                        mat[:3, 3] = pos
#                        import math
#                        from ePMV import comput_util as c
#                        r  = c.matrixToEuler(mat)
#                        h1 = math.degrees(math.pi + r[0])
#                        p1 = math.degrees(r[1])
#                        b1 = math.degrees(-math.pi + r[2])
#                        #angles[0] = 180.0+angles[0]
#                        #angles[2] = 180.0-angles[2] 
#                        #hmat = self.afviewer.vi.FromMat(mat,transpose=True)
#                        #rot = utils.MatrixToHPB(hmat)
#                        print 'rot from matrix = ', r,h1,p1,b1
#                    # END: newer code from Theis version added July 5, 2012
#                    result.append([pos,rot])
#                    pt3d = result[matCount][0]
#                    x, y, z = pt3d #  ADDDED this line back from newer code from Theis version added July 5, 2012             
#                    vTranslationString += str(x)+ ",\t" + str(y) + ",\t" + str(z) + "\n"
#                    vRotationString += str(h1)+ ",\t" + str(p1) + ",\t" + str(b1) + ",\t" + ingr.name +"\n"  #  ADDDED this line back from newer code from Theis version 
#                    matCount += 1
#                rfile = open(resultfilenameT, 'w')
#                rfile.write( vTranslationString )
#                rfile.close()
#                
#                rfile = open(resultfilenameR, 'w')
#                rfile.write( vRotationString )
#                rfile.close()
#                print ('len(result) = ', len(result))
#                print ('len(self.molecules) = ', len(self.molecules))
                ### Graham Note:  There is overused disk space- the rotation matrix is 4x4 with an offset of 0,0,0 and we have a separate translation vector in the results and molecules arrays.  Get rid of the translation vector and move it to the rotation matrix to save space... will that slow the time it takes to extract the vector from the matrix when we need to call it?       
                print ('*************************************************** vDistance String Should be on')
                print ('unitVolume2 = ', unitVol)
                print ('Number of Points Unused = ', unUsedPts)
                print ('Number of Points Used   = ', usedPts)
                print ('Volume Used   = ', usedPts*unitVol)
                print ('Volume Unused = ', unUsedPts*unitVol)
                print ('vTestid = ', vTestid)
                print ('self.nbGridPoints = ', self.nbGridPoints)
                print ('self.gridVolume = ', self.gridVolume)    
                print ('histoVol.timeUpDistLoopTotal = ', self.timeUpDistLoopTotal)

                            
            #totalVolume = self.grid.gridVolume*unitVol
            #fpts = self.nbFreePointsAfterFill
    #        print 'self.freePointsAfterFill = ', self.freePointsAfterFill
            #print 'nnbFreePointsAfterFill = ', self.nbFreePointsAfterFill
            #print 'Total Points = ', self.grid.gridVolume
            #print 'Total Volume = ', totalVolume
    #    END Analysis Tools: Graham added back this big chunk of code for analysis tools and graphic on 5/16/12 Needs to be cleaned up into a function and proper uPy code   
        print('time to save end', time.time()-t2)            
        if self.afviewer is not None and hasattr(self.afviewer,"vi"):
            self.afviewer.vi.progressBar(label="Filling Complete")
            self.afviewer.vi.resetProgressBar()
        ingredients ={}
        for pos, rot, ingr, ptInd in self.molecules:
            if ingr.name not  in ingredients :
                ingredients[ingr.name]=[ingr,[],[],[]]
            mat = rot.copy()
            mat[:3, 3] = pos
            ingredients[ingr.name][1].append(pos)
            ingredients[ingr.name][2].append(rot)
            ingredients[ingr.name][3].append(numpy.array(mat))
        for o in self.compartments:
            for pos, rot, ingr, ptInd in o.molecules:
                if ingr.name not  in ingredients :
                    ingredients[ingr.name]=[ingr,[],[],[]]
                mat = rot.copy()
                mat[:3, 3] = pos
                ingredients[ingr.name][1].append(pos)
                ingredients[ingr.name][2].append(rot)
                ingredients[ingr.name][3].append(numpy.array(mat)) 
        self.ingr_result = ingredients
        if self.treemode == "bhtree" :
            from bhtree import bhtreelib
            bhtreelib.freeBHtree(self.close_ingr_bhtree)
#        bhtreelib.FreeRBHTree(self.close_ingr_bhtree)
#        del self.close_ingr_bhtree
                    
    def displayCancelDialog(self):
        print('Popup CancelBox: if Cancel Box is up for more than 10 sec, close box and continue loop from here')
#        from pyubic.cinema4d.c4dUI import TimerDialog
#        dialog = TimerDialog()
#        dialog.init()
#        dialog.Open(async=True, pluginid=25555589, width=120, height=100)
#        tt=time.time()
        #while dialog.IsOpen():
        #    if time.time()-tt > 5.:
        #        print "time.time()-tt = ", time.time()-tt
        #        dialog.Close()
#        cancel = dialog._cancel
#        cancel=c4d.gui.QuestionDialog('WannaCancel?') # Removed by Graham on July 10, 2012 because it may no longer be needed, but test it TODO
#        return cancel

    def restore(self,result,orgaresult,freePoint):
        #should we used the grid ? the freePoint can be computed
        #result is [pos,rot,ingr.name,ingr.compNum,ptInd]
        #orgaresult is [[pos,rot,ingr.name,ingr.compNum,ptInd],[pos,rot,ingr.name,ingr.compNum,ptInd]...]
        #after restore we can build the grid and fill!
        #ingredient based dictionary
        ingredients={}
        molecules=[]
        for elem in result :
            pos,rot,name,compNum,ptInd = elem
            #needto check the name if it got the comp rule
            ingr = self.getIngrFromName(name,compNum)
            if ingr is not None:
                molecules.append([pos, rot, ingr, ptInd])
                if name not  in ingredients :
                    ingredients[name]=[ingr,[],[],[]]
                mat = rot.copy()
                mat[:3, 3] = pos
                ingredients[name][1].append(pos)
                ingredients[name][2].append(rot)
                ingredients[name][3].append(numpy.array(mat))
        self.molecules = molecules
        if self.exteriorRecipe:
            self.exteriorRecipe.molecules = molecules
        if len(orgaresult) == len(self.compartments):
            for i,o in enumerate(self.compartments):
                molecules=[]
                for elem in orgaresult[i] :
                    pos,rot,name,compNum,ptInd = elem
                    ingr = self.getIngrFromName(name,compNum)
                    if ingr is not None:
                        molecules.append([pos, rot, ingr, ptInd])
                        if name not in ingredients :
                            ingredients[name]=[ingr,[],[],[]]
                        mat = rot.copy()
                        mat[:3, 3] = pos                            
                        ingredients[name][1].append(pos)
                        ingredients[name][2].append(rot)
                        ingredients[name][3].append(numpy.array(mat))
                o.molecules = molecules
        #consider that one filling have occured
        self.cFill = self.nFill
        #if name == None :
        name = "F"+str(self.nFill)
        self.FillName.append(name)
        self.nFill+=1
        self.ingr_result = ingredients
        self.restoreFreePoints(freePoint)
        return ingredients

    def restoreFreePoints(self,freePoint):
        self.freePoints = self.freePointsAfterFill = freePoint
        self.nbFreePointsAfterFill = len(freePoint)   
        self.distanceAfterFill = self.grid.distToClosestSurf
        self.distancesAfterFill= self.grid.distToClosestSurf
        
            
    def load(self,resultfilename=None,restore_grid=True):
        if resultfilename == None:
            resultfilename = self.resultfile
        #check the extension of the filename none, txt or json
        fileName, fileExtension = os.path.splitext(resultfilename)
        if fileExtension == '':
            try :
                result= pickle.load( open(resultfilename,'rb'))
            except :
                print  ("can't read "+resultfilename)
                return [],[],[]
        elif fileExtension == '.apr':     
            try :
                result= pickle.load( open(resultfilename,'rb'))
            except :
                 return self.load_asTxt(resultfilename=resultfilename)
        elif fileExtension == '.txt':     
            return self.load_asTxt(resultfilename=resultfilename)
        elif fileExtension == '.json':
            return self.load_asJson(resultfilename=resultfilename)  
        else :
            print  ("can't read or recognize "+resultfilename)
            return [],[],[]
        #OR 
        #pos, rot, ingr, ptInd = self.molecules
        #pos,rot,ingr.name,ingr.compNum,ptInd
#        orgaresult=[]
#        freePoint=[]
#        for i, orga in enumerate(self.compartments):
#            orfile = open(resultfilename+"ogra"+str(i),'rb')
#            orgaresult.append(pickle.load(orfile))
#            orfile.close()
#        if restore_grid :
#            freePoint = self.loadFreePoint(resultfilename)
#            self.restoreGridFromFile(resultfilename+"grid")#restore grid distance and ptId
#        return result,orgaresult,freePoint
    
    def loadFreePoint(self,resultfilename):
        rfile = open(resultfilename+"freePoints",'rb')
        freePoint = pickle.load(rfile)
        rfile.close()
        return freePoint       
    
    def store(self,resultfilename=None):
        if resultfilename == None:
            resultfilename = self.resultfile
        resultfilename=autopack.fixOnePath(resultfilename)
        rfile = open(resultfilename, 'wb')
        #pickle.dump(self.molecules, rfile)
        #OR 
        result=[]
        for pos, rot, ingr, ptInd in self.molecules:
            result.append([pos,rot,ingr.name,ingr.compNum,ptInd])
        pickle.dump(result, rfile)
        rfile.close()
        for i, orga in enumerate(self.compartments):
            orfile = open(resultfilename+"ogra"+str(i), 'wb')
            result=[]
            for pos, rot, ingr, ptInd in orga.molecules:
                result.append([pos,rot,ingr.name,ingr.compNum,ptInd])
            pickle.dump(result, orfile)
#            pickle.dump(orga.molecules, orfile)
            orfile.close()
        rfile = open(resultfilename+"freePoints", 'wb')
        pickle.dump(self.freePoints, rfile)
        rfile.close()

    @classmethod
    def dropOneIngr(self,pos,rot, ingrname,ingrcompNum,ptInd,rad=1.0):
        line=""
        line+=("<%f,%f,%f>,")% (pos[0],pos[1],pos[2])
        r=rot.reshape(16,)   
        line+=("<")
        for i in range(15):            
            line+=("%f,")% (r[i])
        line+=("%f>,")% (r[15])
        line+="<%f>,<%s>,<%d>,<%d>\n" % (rad,ingrname,ingrcompNum,ptInd)
        return line
   
    @classmethod
    def getOneIngr(self,line):
        elem = line.split("<")
        pos = eval(elem[1][:-2])
        rot = eval(elem[2][:-2])
        rad = eval(elem[3][:-2])
        ingrname = elem[4][:-2]
        ingrcompNum = eval(elem[5][:-2])
        ptInd = eval(elem[6].split(">")[0])
        return pos,rot, ingrname,ingrcompNum,ptInd,rad

#    @classmethod
    def getOneIngrJson(self,ingr,ingrdic):
        for r in ingr.results:  
            ingrdic[ingr.name]["results"].append([r[0]],r[1],)
#        print ("growingr?",ingr,ingr.name,isinstance(ingr, GrowIngrediant))
        if isinstance(ingr, GrowIngrediant) or isinstance(ingr, ActinIngrediant):
            ingr.nbCurve = ingrdic["nbCurve"]
            ingr.listePtLinear = []
            for i in range(ingr.nbCurve):
                ingr.listePtLinear.append( ingrdic["curve"+str(i)] )
#            print ("nbCurve?",ingr.nbCurve,ingrdic["nbCurve"])
        return ingrdic["results"], ingr.name,ingrdic["compNum"],1,ingrdic["encapsulatingRadius"]

    def load_asTxt(self,resultfilename=None):
#        from upy.hostHelper import Helper as helper 
        if resultfilename == None:
            resultfilename = self.resultfile
        rfile = open(resultfilename,'r')
        #needto parse
        result=[]
        orgaresult=[[],]*len(self.compartments)
#        mry90 = helper.rotation_matrix(-math.pi/2.0, [0.0,1.0,0.0])
#        numpy.array([[0.0, 1.0, 0.0, 0.0], 
#                 [-1., 0.0, 0.0, 0.0], 
#                 [0.0, 0.0, 1.0, 0.0], 
#                 [0.0, 0.0, 0.0, 1.0]])
        lines = rfile.readlines()        
        for l in lines :
            if not len(l) or len(l) < 6 : continue
            pos,rot, ingrname,ingrcompNum,ptInd,rad = self.getOneIngr(l)
            #should I multiply here
            r = numpy.array(rot).reshape(4,4)#numpy.matrix(mry90)*numpy.matrix(numpy.array(rot).reshape(4,4))
            if ingrcompNum == 0 :
                result.append([numpy.array(pos),numpy.array(r),ingrname,ingrcompNum,ptInd])
            else :
                orgaresult[abs(ingrcompNum)-1].append([numpy.array(pos),numpy.array(r),ingrname,ingrcompNum,ptInd])
#        for i, orga in enumerate(self.compartments):
#            orfile = open(resultfilename+"ogra"+str(i),'rb')
#            orgaresult.append(pickle.load(orfile))
#            orfile.close()
#        rfile.close()
#        rfile = open(resultfilename+"freePoints",'rb')
        freePoint = []# pickle.load(rfile)
        try :
            rfile = open(resultfilename+"freePoints",'rb')
            freePoint = pickle.load(rfile)
            rfile.close()
        except :
            pass
        return result,orgaresult,freePoint

    def collectResultPerIngredient(self):
        for pos, rot, ingr, ptInd in self.molecules:
            if isinstance(ingr, GrowIngrediant) or isinstance(ingr, ActinIngrediant):
                pass#already store
            else :
                ingr.results.append([pos,rot])
        for i, orga in enumerate(self.compartments):
            for pos, rot, ingr, ptInd in orga.molecules:
                if isinstance(ingr, GrowIngrediant) or isinstance(ingr, ActinIngrediant):
                    pass#already store
                else :
                    ingr.results.append([pos,rot])                                

    def load_asJson(self,resultfilename=None):
#        from upy.hostHelper import Helper as helper 
        if resultfilename == None:
            resultfilename = self.resultfile
        with open(resultfilename, 'r') as fp :#doesnt work with symbol link ?
            self.result_json=json.load(fp)#,indent=4, separators=(',', ': ')
        #needto parse
        result=[]
        orgaresult=[[],]*len(self.compartments)
        r =  self.exteriorRecipe
        if r :
            if "exteriorRecipe" in self.result_json:
                for ingr in r.ingredients:        
                    if ingr.name not in self.result_json["exteriorRecipe"] : continue
                    iresults, ingrname,ingrcompNum,ptInd,rad = self.getOneIngrJson(ingr,self.result_json["exteriorRecipe"][ingr.name])
                    for r in iresults:
                        rot = numpy.array(r[1]).reshape(4,4)#numpy.matrix(mry90)*numpy.matrix(numpy.array(rot).reshape(4,4))
                        result.append([numpy.array(r[0]),rot,ingrname,ingrcompNum,1])
        #organelle ingr
        for orga in self.compartments:
            #organelle surface ingr
            rs =  orga.surfaceRecipe
            if rs :
                if orga.name+"_surfaceRecipe" in self.result_json:
                    for ingr in rs.ingredients:
                        if ingr.name not in self.result_json[orga.name+"_surfaceRecipe"] : continue
                        iresults, ingrname,ingrcompNum,ptInd,rad = self.getOneIngrJson(ingr,self.result_json[orga.name+"_surfaceRecipe"][ingr.name])
                        for r in iresults:
                            rot = numpy.array(r[1]).reshape(4,4)#numpy.matrix(mry90)*numpy.matrix(numpy.array(rot).reshape(4,4))
                            orgaresult[abs(ingrcompNum)-1].append([numpy.array(r[0]),rot,ingrname,ingrcompNum,1])
            #organelle matrix ingr
            ri =  orga.innerRecipe
            if ri :
                if orga.name+"_innerRecipe" in self.result_json:
                    for ingr in ri.ingredients:                    
                        if ingr.name not in self.result_json[orga.name+"_innerRecipe"] : continue
                        iresults, ingrname,ingrcompNum,ptInd,rad = self.getOneIngrJson(ingr,self.result_json[orga.name+"_innerRecipe"][ingr.name])
                        for r in iresults:
                            rot = numpy.array(r[1]).reshape(4,4)#numpy.matrix(mry90)*numpy.matrix(numpy.array(rot).reshape(4,4))
                            orgaresult[abs(ingrcompNum)-1].append([numpy.array(r[0]),rot,ingrname,ingrcompNum,1])
        freePoint = []# pickle.load(rfile)
        try :
            rfile = open(resultfilename+"freePoints",'rb')
            freePoint = pickle.load(rfile)
            rfile.close()
        except :
            pass
        return result,orgaresult,freePoint
        
    def dropOneIngrJson(self,ingr,rdic):
        rdic[ingr.name]={}
        rdic[ingr.name]["compNum"]= ingr.compNum
        rdic[ingr.name]["encapsulatingRadius"]= float(ingr.encapsulatingRadius)
        rdic[ingr.name]["results"]=[] 
        for r in ingr.results:  
            if hasattr(r[0],"tolist"):
                r[0]=r[0].tolist()
            if hasattr(r[1],"tolist"):
                r[1]=r[1].tolist()
            rdic[ingr.name]["results"].append([r[0],r[1]])
        if isinstance(ingr, GrowIngrediant) or isinstance(ingr, ActinIngrediant):
            rdic[ingr.name]["nbCurve"]=ingr.nbCurve
            for i in range(ingr.nbCurve):
                lp = numpy.array(ingr.listePtLinear[i])
                ingr.listePtLinear[i]=lp.tolist()                 
                rdic[ingr.name]["curve"+str(i)] = ingr.listePtLinear[i]
       
    def store_asJson(self,resultfilename=None):
        if resultfilename == None:
            resultfilename = self.resultfile
        resultfilename=autopack.fixOnePath(resultfilename)
        self.collectResultPerIngredient()
        self.result_json={}
        r =  self.exteriorRecipe
        if r :
            self.result_json["exteriorRecipe"]={}
            for ingr in r.ingredients:
                self.dropOneIngrJson(ingr,self.result_json["exteriorRecipe"])

        #compartment ingr
        for orga in self.compartments:
            #compartment surface ingr
            rs =  orga.surfaceRecipe
            if rs :
                self.result_json[orga.name+"_surfaceRecipe"]={}
                for ingr in rs.ingredients:
                    self.dropOneIngrJson(ingr,self.result_json[orga.name+"_surfaceRecipe"])
            #compartment matrix ingr
            ri =  orga.innerRecipe
            if ri :
                self.result_json[orga.name+"_innerRecipe"]={}
                for ingr in ri.ingredients:
                    self.dropOneIngrJson(ingr,self.result_json[orga.name+"_innerRecipe"])
        with open(resultfilename+".json", 'w') as fp :#doesnt work with symbol link ?
            json.dump(self.result_json,fp,indent=4, separators=(',', ': '))#,indent=4, separators=(',', ': ')
        
    def store_asTxt(self,resultfilename=None):
        if resultfilename == None:
            resultfilename = self.resultfile
        resultfilename=autopack.fixOnePath(resultfilename)
        rfile = open(resultfilename+".txt", 'w')#doesnt work with symbol link ?
        #pickle.dump(self.molecules, rfile)
        #OR 
        result=[]
        line=""
        for pos, rot, ingr, ptInd in self.molecules:
            line+=self.dropOneIngr(pos,rot,ingr.name,ingr.compNum,ptInd,rad=ingr.encapsulatingRadius)
            #result.append([pos,rot,ingr.name,ingr.compNum,ptInd])
        rfile.write(line)
        #write the curve point 

        rfile.close()
        for i, orga in enumerate(self.compartments):
            orfile = open(resultfilename+"ogra"+str(i)+".txt", 'w')
            result=[]
            line=""
            for pos, rot, ingr, ptInd in orga.molecules:
                line+=self.dropOneIngr(pos,rot,ingr.name,ingr.compNum,ptInd,rad=ingr.encapsulatingRadius)
            orfile.write(line)
#            pickle.dump(orga.molecules, orfile)
            orfile.close()
#        rfile = open(resultfilename+"freePoints", 'w')
#        pickle.dump(self.freePoints, rfile)
#        rfile.close()

    @classmethod
    def convertPickleToText(self,resultfilename=None,norga=0):
        if resultfilename == None:
            resultfilename = self.resultfile
        rfile = open(resultfilename)
        result= pickle.load( rfile)
        orgaresult=[]
        for i in range(norga):
            orfile = open(resultfilename+"ogra"+str(i))
            orgaresult.append(pickle.load(orfile))
            orfile.close()
        rfile.close()
        rfile = open(resultfilename+"freePoints")
        freePoint = pickle.load(rfile)
        rfile.close()
        rfile = open(resultfilename+".txt", 'w')
        line=""
        for pos, rot, ingrName,compNum, ptInd in result:
            line+=self.dropOneIngr(pos,rot,ingrName,compNum,ptInd)
            #result.append([pos,rot,ingr.name,ingr.compNum,ptInd])
        rfile.write(line)
        rfile.close()
        for i in range(norga):
            orfile = open(resultfilename+"ogra"+str(i)+".txt", 'w')
            result=[]
            line=""
            for pos, rot, ingrName,compNum, ptInd in orgaresult[i]:
                line+=self.dropOneIngr(pos,rot,ingrName,compNum,ptInd)
            orfile.write(line)
#            pickle.dump(orga.molecules, orfile)
            orfile.close()
         #freepoint
         
    def printFillInfo(self):
        r = self.exteriorRecipe
        if r is not None:
            print('    histoVol exterior recipe:')
            r.printFillInfo('        ')
            
        for o in self.compartments:
            o.printFillInfo()

    def finishWithWater(self,freePoints=None,nbFreePoints=None):
        #self.freePointsAfterFill[:self.nbFreePointsAfterFill]
        water = [ ( 0.000, 0.000, 0.0), #0
                  ( 0.757, 0.586, 0.0), #H
                  (-0.757, 0.586, 0.0)] #H
        #object?
        #sphere sphere of 2.9A
        waterDiam = 2.9
        if freePoints is None :
            freePoints = self.freePointsAfterFill
        if nbFreePoints is None :
            nbFreePoints = self.nbFreePointsAfterFill
        #a freepoint is a voxel, how many water in the voxel
        voxelsize = self.grid.gridSpacing
        nbWaterPerVoxel = voxelsize / waterDiam
        #coords masterGridPositions

    def estimateVolume(self,boundingBox, spacing):
        #need to box N point and coordinaePoint
#        xl,yl,zl = boundingBox[0]
#        xr,yr,zr = boundingBox[1]
#        realTotalVol = (xr-xl)*(yr-yl)*(zr-zl)
        grid = Grid()
        grid.boundingBox = boundingBox
        grid.gridSpacing = spacing# = self.smallestProteinSize*1.1547  # 2/sqrt(3)????
        grid.gridVolume,grid.nbGridPoints = self.callFunction(grid.computeGridNumberOfPoint,(boundingBox,spacing))
        unitVol = spacing**3
        realTotalVol = grid.gridVolume*unitVol
        
        r = self.exteriorRecipe
        if r :
            r.setCount(realTotalVol,reset=False)
        for o in self.compartments:
            o.estimateVolume(hBB=grid.boundingBox)
            rs = o.surfaceRecipe
            if rs :
                realTotalVol = o.surfaceVolume
                rs.setCount(realTotalVol,reset=False)
            ri = o.innerRecipe
            if ri :
                realTotalVol = o.interiorVolume
                ri.setCount(realTotalVol,reset=False)

    def estimateVolume_old(self,boundingBox, spacing):
        #need to box N point and coordinaePoint
        pad = 10.0
        grid = Grid()
        grid.boundingBox = boundingBox
        grid.gridSpacing = spacing# = self.smallestProteinSize*1.1547  # 2/sqrt(3)????
        grid.gridVolume,grid.nbGridPoints = self.callFunction(grid.computeGridNumberOfPoint,(boundingBox,spacing))
        nbPoints = grid.gridVolume            
        # compute 3D point coordiantes for all grid points
        self.callFunction(grid.create3DPointLookup) 
        grid.gridPtId = [0]*nbPoints
        xl,yl,zl = boundingBox[0]
        xr,yr,zr = boundingBox[1]
        realTotalVol = (xr-xl)*(yr-yl)*(zr-zl)
        print ("totalVolume %f for %d points" % (realTotalVol,nbPoints))
        # distToClosestSurf is set to self.diag initially
        grid.diag = diag = vlen( vdiff((xr,yr,zr), (xl,yl,zl) ) )
        distance  = grid.distToClosestSurf = [diag]*nbPoints
        #foreach ingredient get estimation of insidepoint and report the percantage of total point in Volume
        r = self.exteriorRecipe
        if r :
            for ingr in r.ingredients:
                insidePoints,newDistPoints=ingr.getInsidePoints(grid,grid.masterGridPositions,pad,distance,
                       centT=ingr.positions[-1],jtrans=[0.,0.,0.], rotMatj=numpy.identity(4))
                ingr.nbPts = len(insidePoints)
                onemol = (realTotalVol* float(ingr.nbPts) )/ float(nbPoints)
                ingr.vol_nbmol = int(ingr.molarity * onemol)
                print ("ingr %s has %d points representing %f for one mol thus %d mol" %(ingr.name, ingr.nbPts, onemol, ingr.vol_nbmol))
#                ingr.vol_nbmol = ?
        for o in self.compartments:
            rs = o.surfaceRecipe
            if rs :
                for ingr in rs.ingredients:
                    insidePoints,newDistPoints=ingr.getInsidePoints(grid,grid.masterGridPositions,pad,distance,
                           centT=ingr.positions[-1],jtrans=[0.,0.,0.], rotMatj=numpy.identity(4))
                    ingr.nbPts = len(insidePoints)
                    onemol = (realTotalVol* float(ingr.nbPts) )/ float(nbPoints)
                    ingr.vol_nbmol = int(ingr.molarity * onemol)
            ri = o.innerRecipe
            if ri :
                for ingr in ri.ingredients:
                    insidePoints,newDistPoints=ingr.getInsidePoints(grid,grid.masterGridPositions,pad,distance,
                           centT=ingr.positions[-1],jtrans=[0.,0.,0.], rotMatj=numpy.identity(4))
                    ingr.nbPts = len(insidePoints)
                    onemol = (realTotalVol* float(ingr.nbPts) )/ float(nbPoints)
                    ingr.vol_nbmol = int(ingr.molarity * onemol)



#==============================================================================
# AFter this point, features development around physics engine and algo
# octree
# panda bullet
# panda ode
#==============================================================================

    def setupOctree(self,):
        if self.octree is None :
#            from autopack.octree import Octree
            from autopack import octree_exteneded as octree
            from autopack.octree_exteneded import Octree
            octree.MINIMUM_SIZE=self.smallestProteinSize
            octree.MAX_OBJECTS_PER_NODE=10
            self.octree = Octree(self.grid.getRadius(),helper=helper)#Octree((0,0,0),self.grid.getRadius())   #0,0,0 or center of grid?         
            
    def setupPanda(self,):
        self.rb_func_dic = {"bullet":{
        "SingleSphere":self.addSingleSphereRB,
        "SingleCube":self.addSingleCubeRB,
        "MultiSphere":self.addMultiSphereRB,
        "MultiCylinder":self.addMultiCylinderRB,
        "Grow":self.addMultiCylinderRB,        
        "Mesh":self.addMeshRB,
        },"ode":
            {"SingleSphere":self.addSingleSphereRBODE,}
            }
        if self.world is None :
            if panda3d is None :
                return
            from panda3d.core import loadPrcFileData
            
            loadPrcFileData("", "window-type none" ) 
            # Make sure we don't need a graphics engine 
            #(Will also prevent X errors / Display errors when starting on linux without X server)
            loadPrcFileData("", "audio-library-name null" ) # Prevent ALSA errors 
#            loadPrcFileData('', 'bullet-enable-contact-events true')
            loadPrcFileData('', 'bullet-max-objects 10240')#10240
            loadPrcFileData('', 'bullet-broadphase-algorithm sap')#aabb
            loadPrcFileData('', 'bullet-sap-extents '+str(self.grid.diag))#
            
            import direct.directbase.DirectStart 
            from panda3d.core import Vec3
#            bullet.bullet-max-objects = 1024 * 10#sum of all predicted n Ingredient ?
            if self.panda_solver == "bullet":
                from panda3d.bullet import BulletWorld               
                self.worldNP = render.attachNewNode('World')            
                self.world = BulletWorld()
                self.BitMask32 = BitMask32
            elif self.panda_solver == "ode":
                from panda3d.ode import OdeWorld, OdeQuadTreeSpace, OdeHashSpace
                self.world = OdeWorld()    
                #space ?OdeQuadTreeSpace 
                #or hashspace ?
                self.ode_space =OdeHashSpace()#OdeQuadTreeSpace(center,extends,depth)
                self.ode_space.set_levels(-2,6)
                self.ode_space.setAutoCollideWorld(self.world)
                
            self.world.setGravity(Vec3(0, 0, 0))
            self.static=[]
            self.moving = None
            self.rb_panda = []

    def delRB(self, node):
        if panda3d is None :
                return
        if self.panda_solver == "bullet":
            self.world.removeRigidBody(node)
            np = NodePath(node)
            if np is not None : 
                np.removeNode()
        elif self.panda_solver == "ode":
            node.destroy()
            
        if node in self.rb_panda: self.rb_panda.pop(self.rb_panda.index(node))
        if node in self.static: self.static.pop(self.static.index(node))
        if node == self.moving: self.moving = None

    def addSingleSphereRBODE(self,ingr,pMat,jtrans,rotMat):        
        body = OdeBody(self.world)
        M = OdeMass()
        M.setSphereTotal(1.0, ingr.encapsulatingRadius)
        body.setMass(M)
        body.setPosition(Vec3(jtrans[0],jtrans[1],jtrans[2]))
        body.setRotation(pMat)
        #the geometry for the collision ?
        geom = OdeSphereGeom(self.ode_space, ingr.encapsulatingRadius)
        geom.setBody(body)
        return geom

    def addSingleSphereRB(self,ingr,pMat,jtrans,rotMat):
        shape = BulletSphereShape(ingr.encapsulatingRadius)
        inodenp = self.worldNP.attachNewNode(BulletRigidBodyNode(ingr.name))
        inodenp.node().setMass(1.0)
#        inodenp.node().addShape(shape)
        inodenp.node().addShape(shape,TransformState.makePos(Point3(0, 0, 0)))#rotation ?
#        spherenp.setPos(-2, 0, 4)
        return inodenp
        
    def addMultiSphereRB(self,ingr,pMat,jtrans,rotMat):
        inodenp = self.worldNP.attachNewNode(BulletRigidBodyNode(ingr.name))
        inodenp.node().setMass(1.0)
        centT = ingr.positions[0]#ingr.transformPoints(jtrans, rotMat, ingr.positions[0])
        for radc, posc in zip(ingr.radii[0], centT):
            shape = BulletSphereShape(radc)
            inodenp.node().addShape(shape, TransformState.makePos(Point3(posc[0],posc[1],posc[2])))#
        return inodenp

    def multiSphereRB(self,name,pos,rad):
        inodenp = self.worldNP.attachNewNode(BulletRigidBodyNode(name))
        inodenp.node().setMass(1.0)
        #centT = ingr.positions[0]#ingr.transformPoints(jtrans, rotMat, ingr.positions[0])
#        for i in range(len(pos)):#
#            posc = pos[i]
#            radc = rad[i]
        for radc, posc in zip(rad,pos):
            shape = BulletSphereShape(radc)
            inodenp.node().addShape(shape, TransformState.makePos(Point3(posc[0],posc[1],posc[2])))#
        return inodenp
        
    def addSingleCubeRB(self,ingr,pMat,jtrans,rotMat):
        rt = "Box"
        halfextents= ingr.bb[1]
        print (halfextents)
        shape = BulletBoxShape(Vec3(halfextents[0], halfextents[1], halfextents[2]))#halfExtents
        inodenp = self.worldNP.attachNewNode(BulletRigidBodyNode(ingr.name))
        inodenp.node().setMass(1.0)
#        inodenp.node().addShape(shape)
        inodenp.node().addShape(shape,TransformState.makePos(Point3(0, 0, 0)))#, pMat)#TransformState.makePos(Point3(jtrans[0],jtrans[1],jtrans[2])))#rotation ?
#        spherenp.setPos(-2, 0, 4)
        return inodenp
        
    def addMultiCylinderRB(self,ingr,pMat,jtrans,rotMat):
        helper = autopack.helper
        inodenp = self.worldNP.attachNewNode(BulletRigidBodyNode(ingr.name))
        inodenp.node().setMass(1.0)
        centT1 = ingr.positions[0]#ingr.transformPoints(jtrans, rotMat, ingr.positions[0])
        centT2 = ingr.positions2[0]#ingr.transformPoints(jtrans, rotMat, ingr.positions2[0])
        for radc, p1, p2 in zip(ingr.radii[0], centT1, centT2):
            length, mat = helper.getTubePropertiesMatrix(p1,p2)
            pMat = self.pandaMatrice(mat)
#            d = numpy.array(p1) - numpy.array(p2)
#            s = numpy.sum(d*d)
            shape = BulletCylinderShape(radc, length,1)#math.sqrt(s), 1)# { XUp = 0, YUp = 1, ZUp = 2 } or LVector3f const half_extents
            inodenp.node().addShape(shape, TransformState.makeMat(pMat))#
        return inodenp

    def addMeshRBOld(self,ingr,pMat,jtrans,rotMat):
        helper = autopack.helper
        if ingr.mesh is None:
            return
        faces,vertices,vnormals = helper.DecomposeMesh(ingr.mesh,
                                edit=False,copy=False,tri=True,transform=True)        
        from panda3d.bullet import BulletTriangleMesh,BulletTriangleMeshShape
        p0 = Point3(-10, -10, 0)
        p1 = Point3(-10, 10, 0)
        p2 = Point3(10, -10, 0)
        p3 = Point3(10, 10, 0)
        mesh = BulletTriangleMesh()
        points3d = [Point3(v[0],v[1],v[2]) for v in vertices]
        for f in faces:
            mesh.addTriangle(points3d[f[0]],points3d[f[1]],points3d[f[2]])
            
        shape = BulletTriangleMeshShape(mesh, dynamic=False)
        inodenp = self.worldNP.attachNewNode(BulletRigidBodyNode(ingr.name))
        inodenp.node().setMass(1.0)
        inodenp.node().addShape(shape,TransformState.makePos(Point3(0, 0, 0)))#, pMat)#TransformState.makePos(Point3(jtrans[0],jtrans[1],jtrans[2])))#rotation ?
        return inodenp

    def setGeomFaces(self,tris,face):                
        #have to add vertices one by one since they are not in order
        if len(face) == 2 :
            face = numpy.array([face[0],face[1],face[1],face[1]],dtype='int')
        for i in face :        
            tris.addVertex(i)
        tris.closePrimitive()


    def addMeshRB(self,ingr,pMat,jtrans,rotMat):
        inodenp = self.worldNP.attachNewNode(BulletRigidBodyNode(ingr.name))
        inodenp.node().setMass(1.0)
        helper = autopack.helper
        if ingr.mesh is None:
            return
        faces,vertices,vnormals = helper.DecomposeMesh(ingr.mesh,
                               edit=False,copy=False,tri=True,transform=True)        
        from panda3d.core import GeomVertexFormat,GeomVertexWriter,GeomVertexData,Geom,GeomTriangles
        from panda3d.core import GeomVertexReader
        from panda3d.bullet import BulletTriangleMesh,BulletTriangleMeshShape,BulletConvexHullShape
        #step 1) create GeomVertexData and add vertex information
        format=GeomVertexFormat.getV3()
        vdata=GeomVertexData("vertices", format, Geom.UHStatic)        
        vertexWriter=GeomVertexWriter(vdata, "vertex")
        [vertexWriter.addData3f(v[0],v[1],v[2]) for v in vertices]

        #step 2) make primitives and assign vertices to them
        tris=GeomTriangles(Geom.UHStatic)
        [self.setGeomFaces(tris,face) for face in faces]

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

    def addMeshRBOrganelle(self,o):
        helper = autopack.helper
        geom =   helper.getObject(o.gname)      
        if geom is None :
            o.gname = '%s_Mesh'%o.name            
            geom = helper.getObject(o.gname)

        inodenp = self.worldNP.attachNewNode(BulletRigidBodyNode(o.name))
        inodenp.node().setMass(1.0)

        faces,vertices,vnormals = helper.DecomposeMesh(geom,
                               edit=False,copy=False,tri=True,transform=True)        
        from panda3d.core import GeomVertexFormat,GeomVertexWriter,GeomVertexData,Geom,GeomTriangles
        from panda3d.core import GeomVertexReader
        from panda3d.bullet import BulletTriangleMesh,BulletTriangleMeshShape,BulletConvexHullShape
        #step 1) create GeomVertexData and add vertex information
        format=GeomVertexFormat.getV3()
        vdata=GeomVertexData("vertices", format, Geom.UHStatic)        
        vertexWriter=GeomVertexWriter(vdata, "vertex")
        [vertexWriter.addData3f(v[0],v[1],v[2]) for v in vertices]

        #step 2) make primitives and assign vertices to them
        tris=GeomTriangles(Geom.UHStatic)
        [self.setGeomFaces(tris,face) for face in faces]

        #step 3) make a Geom object to hold the primitives
        geom=Geom(vdata)
        geom.addPrimitive(tris)
        #step 4) create the bullet mesh and node
        mesh = BulletTriangleMesh()
        mesh.addGeom(geom)    
        shape = BulletTriangleMeshShape(mesh, dynamic=False)#BulletConvexHullShape
        #or
        #shape = BulletConvexHullShape()
        #shape.add_geom(geom)
        print ("shape ok",shape)
        #inodenp = self.worldNP.attachNewNode(BulletRigidBodyNode(ingr.name))
        #inodenp.node().setMass(1.0)
        inodenp.node().addShape(shape)#,TransformState.makePos(Point3(0, 0, 0)))#, pMat)#TransformState.makePos(Point3(jtrans[0],jtrans[1],jtrans[2])))#rotation ?

        if self.panda_solver == "bullet":
            inodenp.setCollideMask(BitMask32.allOn())
            inodenp.node().setAngularDamping(1.0)
            inodenp.node().setLinearDamping(1.0)
#            inodenp.setMat(pmat)
            self.world.attachRigidBody(inodenp.node())
            inodenp = inodenp.node()
        return inodenp

    def pandaMatrice(self,mat):
        mat = mat.transpose().reshape((16,))
#        print mat,len(mat),mat.shape
        pMat = Mat4(mat[0],mat[1],mat[2],mat[3],
                   mat[4],mat[5],mat[6],mat[7],
                   mat[8],mat[9],mat[10],mat[11],
                   mat[12],mat[13],mat[14],mat[15],)
        return pMat
        
    
    def addRB(self,ingr, trans, rotMat, rtype="SingleSphere",static=False):
        # Sphere
        if panda3d is None :
            return None
        mat = rotMat.copy()
#        mat[:3, 3] = trans
#        mat = mat.transpose()
        mat = mat.transpose().reshape((16,))
#        print mat,len(mat),mat.shape
        mat3x3 =Mat3(mat[0],mat[1],mat[2],
                     mat[4],mat[5],mat[6],
                     mat[8],mat[9],mat[10])
        pmat   = Mat4(mat[0],mat[1],mat[2],mat[3],
                                           mat[4],mat[5],mat[6],mat[7],
                                           mat[8],mat[9],mat[10],mat[11],
                                           trans[0],trans[1],trans[2],mat[15],)
        pMat = TransformState.makeMat(pmat)
        if self.panda_solver == "ode":
            pMat = mat3x3
        shape = None
        inodenp = None
#        print (pMat)         
        if ingr.use_mesh_rb:
            rtype = "Mesh"
            #print ("#######RBNode Mesh ####", ingr.name, ingr.rbnode,self.rb_func_dic[rtype])
        inodenp = self.rb_func_dic[self.panda_solver][rtype](ingr,pMat,trans,rotMat)
        if self.panda_solver == "bullet":
            inodenp.setCollideMask(BitMask32.allOn())
            inodenp.node().setAngularDamping(1.0)
            inodenp.node().setLinearDamping(1.0)
            inodenp.setMat(pmat)
            self.world.attachRigidBody(inodenp.node())
            inodenp = inodenp.node()
        elif self.panda_solver == "ode":
            inodenp.setCollideBits(BitMask32(0x00000002))
            inodenp.setCategoryBits(BitMask32(0x00000001))
            #boxGeom.setBody(boxBody)
        self.rb_panda.append(inodenp)
        #self.moveRBnode(inodenp.node(), trans, rotMat)
        return inodenp
    
    def moveRBnode(self,node, trans, rotMat):
        mat = rotMat.copy()
#        mat[:3, 3] = trans
#        mat = mat.transpose()
        mat = mat.transpose().reshape((16,))
#        print mat,len(mat),mat.shape
        if self.panda_solver == "bullet":
            pMat = Mat4(mat[0],mat[1],mat[2],mat[3],
                       mat[4],mat[5],mat[6],mat[7],
                       mat[8],mat[9],mat[10],mat[11],
                       trans[0],trans[1],trans[2],mat[15],)
            pTrans = TransformState.makeMat(pMat)
            nodenp = NodePath(node)
            nodenp.setMat(pMat)
        elif self.panda_solver == "ode":
            mat3x3 =Mat3(mat[0],mat[1],mat[2],
                     mat[4],mat[5],mat[6],
                     mat[8],mat[9],mat[10])
            body = node.get_body()
            body.setPosition(Vec3(trans[0],trans[1],trans[2]))
            body.setRotation(mat3x3)
            
#        nodenp.setPos(trans[0],trans[1],trans[2])
#        print nodenp.getPos()
    
    def getRotTransRB(self,node ):
        nodenp = NodePath(node)
        m=nodenp.getMat()
        M = numpy.array(m)
        rRot = numpy.identity(4)
        rRot[:3,:3] = M[:3,:3]
        rTrans = M[3,:3]
        return rTrans,rRot
        
    def runBullet(self,ingr,simulationTimes, runTimeDisplay):
        done = False
        t1=time()
        simulationTimes = 5.0
        while not done:
            #should do it after a jitter run
#        for i in xrange(10):
            dt = globalClock.getDt()
            self.world.doPhysics(dt, 100,1.0/500.0)#world.doPhysics(dt, 10, 1.0/180.0)100, 1./500.#2, 1.0/120.0
            #check number of contact betwee currrent and rest ?
            r=[ (self.world.contactTestPair(self.moving, n).getNumContacts() > 0 ) for n in self.static]
            done=not ( True in r)
            print (done,dt,time()-t1)
            if runTimeDisplay :
                #move self.moving and update
                nodenp = NodePath(self.moving)
                ma=nodenp.getMat()
                self.afviewer.vi.setObjectMatrix(ingr.moving_geom,numpy.array(ma),transpose=False)#transpose ?
                self.afviewer.vi.update()
            if (time()-t1) > simulationTimes :
                done = True
                break
            
        
