# -*- coding: utf-8 -*-
"""
    Created on Wed Jun 26 10:59:01 2015
    @author: ludo
    Hoe to run autopack from xml/json file
    """
import os
import sys

#pass in argument the output data folder
#you could actually save a copy of the recipe in the output data folder

if len(sys.argv) < 2 :
    outputFile = os.getcwd() + "/dataERROR"
else :
    outputFile = str(sys.argv[1])

import json
import time

settingsFile = open(outputFile + '/settings.json')
#parameter values applied later
parVals = json.load(settingsFile)
settingsFile.close()

mglPath =parVals['running_settings']['mglPath']
sys.path.append(mglPath)
sys.path.append(mglPath+"/PIL")

print mglPath
try :
    from collections import OrderedDict
    #this work with python2.7
    #otherwise need the one from mgltools
except :
    from ordereddict import OrderedDict

import numpy as np
from scipy import spatial
#for plotting results
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
from matplotlib import image as mpimg
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d

import autopack# as autopack
#where is autopack
localdir = wrkDir = autopack.__path__[0]


from autopack.shapeGenerator import *

#autopack setup
from autopack.Environment import Environment
from autopack.Ingredient import KWDS as ingredients_parameters
from autopack.Analysis import AnalyseAP

from math import sqrt
import matplotlib.pyplot as plt

import scipy.cluster.hierarchy as sch

#from numbapro import jit, int32, float32, complex64
#from numbapro import jit,cuda, int16, float32, from_dtype,void
#from timeit import default_timer as timer


#get the helper
#deal with geometry and such
helper = autopack.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper = helperClass(vi="nogui")
autopack.helper = helper
#==============================================================================
# some usefule function
#==============================================================================

from direct.showbase.ShowBase import ShowBase 

#from direct.filter.CommonFilters import CommonFilters
#from direct.filter.FilterManager import FilterManager

from panda3d.core import DirectionalLight,Vec4,Vec3
from panda3d.core import loadPrcFileData 
from panda3d.core import ColorAttrib, Material
from panda3d.core import WindowProperties

#loadPrcFileData("",
#"""
#   load-display p3tinydisplay # to force CPU only rendering (to make it available as an option if everything else fail, use aux-display p3tinydisplay)
#   window-type offscreen # Spawn an offscreen buffer (use window-type none if you don't need any rendering)
#   audio-library-name null # Prevent ALSA errors
#   show-frame-rate-meter 0
#   sync-video 0
#""")
# add origin ? or border
def distanceMatrix(h,filename,avg=True):
    h.collectResultPerIngredient()
    #order result per ingredient id
    #Y = cdist(XA, XB, 'euclidean')
    ing_pos = np.array([np.array(m[0]).tolist() for m in h.molecules])
    ing_names = [m[2].name for m in h.molecules]
    #Returns a condensed distance matrix Y. For each i and j (where i<j<n), the metric dist(u=X[i], v=X[j]) is computed and stored in entry ij.
    #i,j=i+1
    sorted_indices = np.argsort(ing_names)
    s_pos = ing_pos[sorted_indices]
    #md = spatial.distance.pdist(s_pos)
    actualmatrice=spatial.distance.cdist(s_pos,s_pos)
    #sort bvy ingredient type
    plt.matshow(actualmatrice,cmap="Reds")
    plt.savefig(filename+".png")
    plt.close('all')     # closes the cu
    np.savetxt(filename+".csv", actualmatrice, delimiter=",")
    #numpy.unique(ar, return_index=False, return_inverse=False, return_counts=False)[source]
    res=np.unique(np.array(ing_names)[sorted_indices],return_index=True,return_counts=True)
    #res[2] is the count of instance
    newmatrics = np.zeros((len(res[0]),len(res[0])))
    si=0
    sj=0    
    ei=0
    ej=0
    for i in range(len(res[0])):#0-20
        ni = res[2][i]
        si = res[2][:i].sum()
        ei = (si+ni)
        for j in range(len(res[0])):
            #average the distance i,j
            nj = res[2][j]
            sj = res[2][:j].sum()
            ej = sj + nj
            avg = actualmatrice[si:ei,sj:ej].mean()
            newmatrics[i,j]=avg
    plt.close('all')     
    plt.matshow(newmatrics,cmap="Reds")
    plt.savefig(filename+"_avg.png")
    ingredientNames = ','.join(str(a) for a in res[0])
    np.savetxt(filename+"_avg.csv", newmatrics, delimiter=",", header="#," + str(ingredientNames), comments='')
    #write cvs

from matplotlib.colors import Normalize

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)
    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))
        
def distanceMatrixOverlap(h,filename,avg=True):
    h.collectResultPerIngredient()
    #order result per ingredient id
    #Y = cdist(XA, XB, 'euclidean')
    #numpy.matrix.sum
    ing_pos = np.array([np.array(m[0]).tolist() for m in h.molecules])
    ing_rad = np.array([m[2].encapsulatingRadius for m in h.molecules])
    ing_names = [m[2].name for m in h.molecules]
    #Returns a condensed distance matrix Y. For each i and j (where i<j<n), the metric dist(u=X[i], v=X[j]) is computed and stored in entry ij.
    #i,j=i+1
    sorted_indices = np.argsort(ing_names)
    s_pos = ing_pos[sorted_indices]
    s_rad = ing_rad[sorted_indices]
    s_names= np.array(ing_names)[sorted_indices]
    #md = spatial.distance.pdist(s_pos)
    actualmatrice=spatial.distance.cdist(s_pos,s_pos)
    radiusmatrice=np.dstack(np.meshgrid(s_rad, s_rad)).sum(axis=2)
    finalmatrice = actualmatrice - radiusmatrice
    #sort bvy ingredient type
    norm = MidpointNormalize(midpoint=0)
    
#if True:
#    # Compute and plot first dendrogram.
#    fig = plt.figure()
##    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])#left, bottom, width, height
##    Y = sch.linkage(finalmatrice, method='centroid')
##    Z1 = sch.dendrogram(Y, orientation='right')
##    ax1.set_xticks([])
##    ax1.set_yticks([])    
##    # Compute and plot second dendrogram.
##    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])#left, bottom, width, height
##    Y = sch.linkage(finalmatrice, method='single')
##    Z2 = sch.dendrogram(Y)
##    ax2.set_xticks([])
##    ax2.set_yticks([])      
##    # Plot distance matrix.
#    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])#left, bottom, width, height
##    idx1 = Z1['leaves']
##    idx2 = Z2['leaves']
##    D = finalmatrice[idx1,:]
##    D = D[:,idx2]
#    im = axmatrix.matshow(finalmatrice, aspect='auto', origin='lower', cmap=plt.cm.YlGnBu)
#    plt.xticks(np.arange(0,len(s_names)), [])
#    plt.yticks(np.arange(0,len(s_names)), [])
#    axmatrix.set_xticklabels(s_names[idx1],rotation=90)
#    axmatrix.set_yticklabels(s_names[idx2])
#    # Plot colorbar.
##    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
##    plt.colorbar()
#    fig.show()
    fig, axe = plt.subplots()
    fig.set_size_inches(20, 20)
    axe.matshow(finalmatrice,cmap="hot",norm=norm)
    plt.xticks(np.arange(0,len(sorted_indices)), [])
    plt.yticks(np.arange(0,len(sorted_indices)), [])
    axe.set_xticklabels(s_names,rotation=90)
    axe.set_yticklabels(s_names)
#    plt.colorbar(axe)

    plt.savefig(filename+".png")
    plt.close('all')     # closes the cu
    np.savetxt(filename+".csv", actualmatrice, delimiter=",")
    #numpy.unique(ar, return_index=False, return_inverse=False, return_counts=False)[source]
    res=np.unique(np.array(ing_names)[sorted_indices],return_index=True,return_counts=True)
    #res[2] is the count of instance
    newmatrics = np.zeros((len(res[0]),len(res[0])))
    si=0
    sj=0    
    ei=0
    ej=0
    for i in range(len(res[0])):#0-20
        ni = res[2][i]
        si = res[2][:i].sum()
        ei = (si+ni)
        for j in range(len(res[0])):
            #average the distance i,j
            nj = res[2][j]
            sj = res[2][:j].sum()
            ej = sj + nj
            avg = actualmatrice[si:ei,sj:ej].mean()
            newmatrics[i,j]=avg
    plt.close('all')     
    fig, axe = plt.subplots()
    fig.set_size_inches(20, 20)
    axe.matshow(newmatrics,cmap="hot",norm=norm)
    names = res[0].tolist()
    plt.xticks(np.arange(0,len(names)), [])
    plt.yticks(np.arange(0,len(names)), [])
    axe.set_xticklabels(names,rotation=90)
    axe.set_yticklabels(names)
    plt.savefig(filename+"_avg.png")
    ingredientNames = ','.join(str(a) for a in res[0])
    np.savetxt(filename+"_avg.csv", newmatrics, delimiter=",", header="#," + str(ingredientNames), comments='')
    
def distanceMatrixCenter(h,filename,avg=True):
    h.collectResultPerIngredient()
    #order result per ingredient id
    #Y = cdist(XA, XB, 'euclidean')
    ing_pos = np.array([np.array(m[0]).tolist() for m in h.molecules])
    ing_names = [m[2].name for m in h.molecules]
    #Returns a condensed distance matrix Y. For each i and j (where i<j<n), the metric dist(u=X[i], v=X[j]) is computed and stored in entry ij.
    #i,j=i+1
    sorted_indices = np.argsort(ing_names)
    s_pos = ing_pos[sorted_indices]
    #md = spatial.distance.pdist(s_pos)
    bbox = np.array(h.boundingBox)
    center = bbox[1]/2.0   
    s_pos = np.vstack([center,s_pos])
    actualmatrice=spatial.distance.cdist(s_pos,s_pos)
    #sort bvy ingredient type
    plt.matshow(actualmatrice,cmap='Reds')
#    plt.xticks(np.arange(0,len(s_pos)), [])
#    plt.yticks(np.arange(0,len(s_pos)), [])
    plt.colorbar()
#    plt.grid()
    plt.savefig(filename+".png")
    plt.close('all')     # closes the cu
    np.savetxt(filename+".csv", actualmatrice, delimiter=",")
    #numpy.unique(ar, return_index=False, return_inverse=False, return_counts=False)[source]
    res=np.unique(np.array(ing_names)[sorted_indices],return_index=True,return_counts=True)
    #res[2] is the count of instance
    newmatrics = np.zeros((len(res[0])+1,len(res[0])+1))
    si=0
    sj=0    
    ei=0
    ej=0
    ni=1
    for i in range(len(res[0])+1):#0-20
        if i!=0 :
            ni = res[2][i-1]
            si = res[2][:i-1].sum()
            ei = (si+ni)
        for j in range(len(res[0])+1):
            #average the distance i,j
            if j!=0:
                nj = res[2][j-1]
                sj = res[2][:j-1].sum()
                ej = sj + nj
            if i==0:
                if j == 0 :
                    avg = actualmatrice[0,0]#.mean()                
                else :
                    avg = actualmatrice[0,sj:ej].mean()
            else :
                if j == 0 :
                    avg = actualmatrice[si:ej,0].mean()                
                else :
                    avg = actualmatrice[si:ei,sj:ej].mean()
                #avg = actualmatrice[si:ei,sj:ej].mean()
            newmatrics[i,j]=avg
    plt.close('all')   
    fig, ax = plt.subplots()
    fig.set_size_inches(16.5, 16.5)
    m=ax.matshow(newmatrics,cmap='Reds')
    fig.colorbar(m)
    plt.xticks(np.arange(0,len(res[0])+1), [])
    plt.yticks(np.arange(0,len(res[0])+1), [])
    names = res[0].tolist()
    names.insert(0,"center")
    ax.set_xticklabels(names,rotation=90)
    ax.set_yticklabels(names)
#    plt.grid()
    plt.savefig(filename+"_avg.png")
    ingredientNames = ','.join(str(a) for a in names)
    np.savetxt(filename+"_avg.csv", newmatrics, delimiter=",", header="#, " + str(ingredientNames), comments='')
    return newmatrics
    #write cvs    
#==============================================================================
# this creates the result plot
#==============================================================================
def plotOneResult2D(h,filename):
    width = h.boundingBox[1][0]#should be the boundary here ?
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bbox = h.boundingBox
    pos_rad = [[np.array(m[0]).tolist(),m[2].encapsulatingRadius,m[2].color] for m in h.molecules]#jtrans, rotMatj, self, ptInd
    print pos_rad[0]
    for ingr in pos_rad:
        p=ingr[0]
        r=ingr[1]
        #            for i,p in enumerate(ingrpos[ingr]):
        #                print (p,radius[ingr])
        ax.add_patch(Circle((p[0], p[1]),r,edgecolor="black", facecolor=ingr[2]))
        if h.use_periodicity : # autopack.testPeriodicity
            if p[0] < r:
                ax.add_patch(Circle((p[0] + width,p[1]), r, facecolor=ingr[2]))
            elif p[0] > (width-r):
                ax.add_patch(Circle((p[0] - width,p[1]), r, facecolor=ingr[2]))
            if p[1] < r:
                ax.add_patch(Circle((p[0],p[1]+ width), r, facecolor=ingr[2]))
            elif p[1] > (width-r):
                ax.add_patch(Circle((p[0],p[1] - width), r, facecolor=ingr[2]))
        ax.set_aspect(1.0)
    plt.axhline(y=bbox[0][1], color='k')
    plt.axhline(y=bbox[1][1], color='k')
    plt.axvline(x=bbox[0][0], color='k')
    plt.axvline(x=bbox[1][0], color='k')
    plt.axis([bbox[0][0], bbox[1][0],bbox[0][1], bbox[1][1]])
    plt.savefig(filename)
    plt.close('all')     # closes the cu

#==============================================================================
# this creates the result plot 3D case wit XY,XZ,YZ plane
#==============================================================================

def toWhite(color, factor):
    if (color is None) :
        color = [0.0,0.0,1.0]
    newcolor = np.array(color)
    newcolor = newcolor + factor*(np.array([1.0,1.0,1.0])-newcolor)
    return newcolor.tolist()

def toBlack(color, factor):
    if (color is None) :
        color = [0.0,0.0,1.0]
    newcolor = np.array(color)
    newcolor = newcolor - factor*(newcolor)
    return newcolor.tolist()

def getFactor(bb,d,axe):
    L=bb[1][axe]-bb[0][axe]
    D=d-bb[0][axe]
    if (d<bb[0][axe]):
        D=bb[0][axe]
    if (d>bb[1][axe]):
        D=bb[1][axe]
    return D/L
    
def plotOneResult3D(h,filename):
    planes = [[0,1,2],[0,2,1],[1,2,0]]
    p_names = ["XY","XZ","YZ"]
    for i,plane in enumerate(planes) :
        width = h.boundingBox[1][plane[0]]#should be the boundary here ?
        fig = plt.figure()
        ax = fig.add_subplot(111)
        bbox = h.boundingBox
        pos_rad = [[np.array(m[0]).tolist(),m[2].encapsulatingRadius,m[2].color] for m in h.molecules]#jtrans, rotMatj, self, ptInd
        print pos_rad[0]
        for ingr in pos_rad:
            p=ingr[0]
            r=ingr[1]
            p0=p[plane[0]]
            p1=p[plane[1]]
            d=p[plane[2]]
            if (ingr[2] is None):
                ingr[2] = [np.random.uniform(),np.random.uniform(),np.random.uniform()]
            factor = getFactor(bbox,d,plane[2])
            ingr[2]=toBlack(ingr[2], factor)
            ax.add_patch(Circle((p0, p1),r,edgecolor="black", facecolor=ingr[2]))
            if h.use_periodicity : # autopack.testPeriodicity
                if p0 < r:
                    ax.add_patch(Circle((p0 + width,p1), r, facecolor=ingr[2]))
                elif p[0] > (width-r):
                    ax.add_patch(Circle((p0 - width,p1), r, facecolor=ingr[2]))
                if p1 < r:
                    ax.add_patch(Circle((p0,p1+ width), r, facecolor=ingr[2]))
                elif p1 > (width-r):
                    ax.add_patch(Circle((p0,p1 - width), r, facecolor=ingr[2]))
            ax.set_aspect(1.0)
        plt.axhline(y=bbox[0][plane[1]], color='k')
        plt.axhline(y=bbox[1][plane[1]], color='k')
        plt.axvline(x=bbox[0][plane[0]], color='k')
        plt.axvline(x=bbox[1][plane[0]], color='k')
        plt.axis([bbox[0][plane[0]], bbox[1][plane[0]],bbox[0][plane[1]], bbox[1][plane[1]]])
        plt.savefig(filename+"_"+p_names[i]+".png")
        plt.close('all')     # closes the cu


def add3Dcircle(x,y,z,r,color,axe,direction):
    #edgecolor="black"
    p = Circle((x,y), r,facecolor=color)
    axe.add_patch(p)
    art3d.pathpatch_2d_to_3d(p, z=z, zdir=direction)

#could do http://danielrothenberg.com/simple-ray-tracer-in-numbapro-cuda/
#also check out http://vispy.org/index.html    
def plot3DOneResult3D(h,filename):
    azelev=np.array([[0.0,0.0],[90.0,-10.0],[0.0,90.0]])+10.0
    planes = [[0,1,2],[0,2,1],[1,2,0]]
    p_names = ["XY","XZ","YZ"]
    direction = ["x","z","y"]
    for i,plane in enumerate(planes) :
        width = h.boundingBox[1][plane[0]]#should be the boundary here ?
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        bbox = h.boundingBox
        pos_rad = [[np.array(m[0]).tolist(),m[2].encapsulatingRadius,m[2].color] for m in h.molecules]#jtrans, rotMatj, self, ptInd
        print pos_rad[0]
        for ingr in pos_rad:
            p=ingr[0]
            r=ingr[1]
            p0=p[plane[0]]
            p1=p[plane[1]]
            p2=p[plane[2]]
            if (ingr[2] is None):
                ingr[2] = [np.random.uniform(),np.random.uniform(),np.random.uniform()]
            factor = getFactor(bbox,p2,plane[2])
            ingr[2]=toBlack(ingr[2], 1.0-factor)
            add3Dcircle(p0,p1,p2,r,ingr[2],ax,direction[i])
#            ax.add_patch(Circle((p0, p1),r,edgecolor="black", facecolor=ingr[2]))
            if h.use_periodicity : # autopack.testPeriodicity
                if p0 < r:
                    add3Dcircle(p0 + width,p1,p2,r,ingr[2],ax,direction[i])
#                    ax.add_patch(Circle((p0 + width,p1), r, facecolor=ingr[2]))
                elif p[0] > (width-r):
                    add3Dcircle(p0 - width,p1,p2,r,ingr[2],ax,direction[i])
#                    ax.add_patch(Circle((p0 - width,p1), r, facecolor=ingr[2]))
                if p1 < r:
                    add3Dcircle(p0,p1+ width,p2,r,ingr[2],ax,direction[i])
#                    ax.add_patch(Circle((p0,p1+ width), r, facecolor=ingr[2]))
                elif p1 > (width-r):
                    add3Dcircle(p0,p1 - width,p2,r,ingr[2],ax,direction[i])
#                    ax.add_patch(Circle((p0,p1 - width), r, facecolor=ingr[2]))
#            ax.set_aspect(1.0)
#        plt.axhline(y=bbox[0][plane[1]], color='k')
#        plt.axhline(y=bbox[1][plane[1]], color='k')
#        plt.axvline(x=bbox[0][plane[0]], color='k')
#        plt.axvline(x=bbox[1][plane[0]], color='k')
#        plt.axis([bbox[0][plane[0]], bbox[1][plane[0]],bbox[0][plane[1]], bbox[1][plane[1]]])
        ax.legend()
        ax.set_xlim3d(bbox[0][0],bbox[1][0])
        ax.set_ylim3d(bbox[0][1],bbox[1][1])
        ax.set_zlim3d(bbox[0][2],bbox[1][2])
        ax.view_init(elev=azelev[i][0], azim=azelev[i][1])
        plt.savefig(filename+"_"+p_names[i]+".png")
        plt.close('all')     # closes the cu

rnd = lambda x: x*np.random.rand()

def pandaPlot(h,filename,sphereTree=False):
#    import sys
#    sys.path.append("/Developer/Panda3D/lib/")
    h.setupPanda()
    if not hasattr(h.base,"light"):
        dlight = DirectionalLight('dlight')
        dlight.setColor(Vec4(0.8, 0.8, 0.5, 1))
        dlnp = h.base.render.attachNewNode(dlight)
        h.base.render.setLight(dlnp)
        h.base.light = dlnp
    h.base.camLens.setNear(0.1)
    h.base.camLens.setFar(10000000)
    myMaterial = Material()
    myMaterial.setShininess(5.0) #Make this material shiny
    myMaterial.setAmbient((0, 0, 1, 1)) #Make this material blue
    h.base.disableMouse()
    azelev=np.array([[0.0,0.0],[90.0,-10.0],[0.0,90.0]])+10.0
    planes = [[0,1,2],[0,2,1],[1,2,0]]
    p_names = ["XY","XZ","YZ"]
    direction = ["x","z","y"]
    bbox = np.array(h.boundingBox)
    center = bbox[1]/2.0
    d=max(bbox[1])*3.0#what distance so we see everything ?
    cam_pos =[[center[0],center[1],center[2]+bbox[1][2]+d],#XY 
              [center[0],center[1]-bbox[1][1]-d,center[2]],#XZ
              [center[0]+bbox[1][0]+d,center[1],center[2]]]#YZ
    #panda is Z up X right Y look
    pos_rad=[]
    if sphereTree :       
        for pos, rot, ingr, ptInd in h.molecules:
            level = ingr.maxLevel
            px = ingr.transformPoints(pos, rot, ingr.positions[level])
            if ingr.modelType=='Spheres':
                for ii in range(len(ingr.radii[level])):
                    pos_rad.append([px[ii],ingr.radii[level][ii],ingr])
            elif ingr.modelType=='Cylinders':
                px2 = ingr.transformPoints(pos, rot, ingr.positions2[level])
                for ii in range(len(ingr.radii[level])):
                    pos_rad.append([px[ii],ingr.radii[level][ii],ingr])
                    pos_rad.append([px2[ii],ingr.radii[level][ii],ingr])
    else :
        pos_rad = [[np.array(m[0]).tolist(),m[2].encapsulatingRadius,m[2]] for m in h.molecules]#jtrans, rotMatj, self, ptInd           
    spheres = h.base.render.getChildren()        
    for k,ingr in enumerate(pos_rad):
        p=ingr[0]
        r=ingr[1]
        if (ingr[2].color is None):
            ingr[2].color = [np.random.uniform(),np.random.uniform(),np.random.uniform()]            
        c=ingr[2].color
        if k < len(spheres)-3 : 
            sphere = spheres[k+3]# k=0 is the light    
            sphere.show()             
        else :
            sphere = Sphere(1.0)
        if sphere.get_name().find("Geometry") == -1:
            print k,k+3,ingr
        sphere.setPos(p[0],p[1],p[2])
        sphere.setScale(r)
        sphere.setColor(c[0],c[1],c[2],1.0)#alpha ?
        sphere.setMaterial(myMaterial)
        sphere.reparentTo(h.base.render)
    #clean the unusde one
    if (k < len(spheres)-3 ):
        for j in range(k,len(spheres)-3):
            spheres[j].hide()

    for i,plane in enumerate(planes) :
        print "camera is at ", cam_pos[i], center
        h.base.camera.setPos (cam_pos[i][0],cam_pos[i][1],cam_pos[i][2] ) 
        h.base.graphicsEngine.renderFrame() 
        h.base.camera.lookAt(0,0,0)
        h.base.graphicsEngine.renderFrame() 
        h.base.camera.lookAt(center[0],center[1],center[2])
        h.base.light.setHpr(h.base.camera.getHpr())
        h.base.graphicsEngine.renderFrame() 
        h.base.screenshot(namePrefix=filename+"_"+p_names[i]+".png", defaultFilename=0, source=None, imageComment="")

#        print render.getChildren()
#        for m in h.base.render.getChildren():
#            m.removeNode()
#==============================================================================
# this creates the per ingredient plots 3D
#==============================================================================
def plotIngredientResults3D(h,filename):
    planes = [[0,1,2],[0,2,1],[1,3,0]]
    p_names = ["XY","XZ","YZ"]
    r = h.exteriorRecipe
    if r :
        for ingrTYPE in r.ingredients:
            for i,plane in enumerate(planes) :
                width = h.boundingBox[1][plane[0]]#should be the boundary here ?
                fig = plt.figure()
                ax = fig.add_subplot(111)
                bbox = h.boundingBox
                pos_rad = [[np.array(m[0]).tolist(),m[2].encapsulatingRadius,m[2].color] for m in h.molecules if m[2] == ingrTYPE]#jtrans, rotMatj, self, ptInd
                print pos_rad[0]
                for ingr in pos_rad:
                    p=ingr[0]
                    r=ingr[1]
                    p0=p[plane[0]]
                    p1=p[plane[1]]
                    ax.add_patch(Circle((p0, p1),r,edgecolor="black", facecolor=ingr[2]))
                    factor = getFactor(bbox,d,plane[2])
                    ingr[2]=toBlack(ingr[2], factor)
                    if h.use_periodicity : # autopack.testPeriodicity
                        if p0 < r:
                            ax.add_patch(Circle((p0 + width,p1), r, facecolor=ingr[2]))
                        elif p[0] > (width-r):
                            ax.add_patch(Circle((p0 - width,p1), r, facecolor=ingr[2]))
                        if p1 < r:
                            ax.add_patch(Circle((p0,p1+ width), r, facecolor=ingr[2]))
                        elif p1 > (width-r):
                            ax.add_patch(Circle((p0,p1 - width), r, facecolor=ingr[2]))
                    ax.set_aspect(1.0)
                plt.axhline(y=bbox[0][plane[1]], color='k')
                plt.axhline(y=bbox[1][plane[1]], color='k')
                plt.axvline(x=bbox[0][plane[0]], color='k')
                plt.axvline(x=bbox[1][plane[0]], color='k')
                plt.axis([bbox[0][plane[0]], bbox[1][plane[0]],bbox[0][plane[1]], bbox[1][plane[1]]])
                plt.savefig(filename+ ingrTYPE.name +"_"+p_names[i]+".png")
                plt.close('all')     # closes the cu

#==============================================================================
# this creates the per ingredient plots 2D
#==============================================================================
def plotIngredientResults2D(h,filename):
    r = h.exteriorRecipe
    if r :
        for ingrTYPE in r.ingredients:
            width = h.boundingBox[1][0]#should be the boundary here ?
            fig = plt.figure()
            ax = fig.add_subplot(111)
            bbox = h.boundingBox
            pos_rad = [[np.array(m[0]).tolist(),m[2].encapsulatingRadius,m[2].color] for m in h.molecules if m[2] == ingrTYPE]#jtrans, rotMatj, self, ptInd
            for ingr in pos_rad:
                p=ingr[0]
                r=ingr[1]
                #            for i,p in enumerate(ingrpos[ingr]):
                #                print (p,radius[ingr])
                ax.add_patch(Circle((p[0], p[1]),r,edgecolor="black", facecolor=ingr[2]))
                if autopack.testPeriodicity :
                    if p[0] < r:
                        ax.add_patch(Circle((p[0] + width,p[1]), r, facecolor=ingr[2]))
                    elif p[0] > (width-r):
                        ax.add_patch(Circle((p[0] - width,p[1]), r, facecolor=ingr[2]))
                    if p[1] < r:
                        ax.add_patch(Circle((p[0],p[1]+ width), r, facecolor=ingr[2]))
                    elif p[1] > (width-r):
                        ax.add_patch(Circle((p[0],p[1] - width), r, facecolor=ingr[2]))
                ax.set_aspect(1.0)
            plt.axhline(y=bbox[0][1], color='k')
            plt.axhline(y=bbox[1][1], color='k')
            plt.axvline(x=bbox[0][0], color='k')
            plt.axvline(x=bbox[1][0], color='k')
            plt.axis([bbox[0][0], bbox[1][0],bbox[0][1], bbox[1][1]])
            plt.savefig(filename + ingrTYPE.name)
            plt.close('all')     # closes the cu

#==============================================================================
# this creates the density plots per result
#==============================================================================
def plotDensitiesPerRun(h,setn,numberofRun, output,filename):
    r = h.exteriorRecipe
    if r :
        for ingrTYPE in r.ingredients:
            fig = plt.figure(frameon=False)
            for seed in range(0,numberofRun):
                img = mpimg.imread(output+os.sep+"run_" + str(setn) + "_" + str(seed) + "_density_" + ingrTYPE.name + ".png")
                plt.imshow(img, alpha=(1./(seed+1)))
                plt.hold(True)
            plt.axis("off")
            plt.savefig(filename + ingrTYPE.name + ".png", bbox_inches='tight')
            plt.close('all')     # closes the cu


def clear(h,n=0,torestore=None):
    h.reset()
    gfile = None
    h.buildGrid(boundingBox=h.boundingBox,gridFileIn=gfile,rebuild=True ,
                gridFileOut=None,previousFill=False)

def activenmolecules(h,n):
    for i,ingr in enumerate(h.exteriorRecipe.ingredients):
        if i > n-1 :
            ingr.completion=1.0
            ingr.nbMol = 0
        else :
            ingr.completion=0.0

def pack(h,seed,filename,eid):
#    h.boundingBox = [[0.,0.,0.],[200.,200.,200.]]
    h.fill5(seedNum=seed,verbose = 1,usePP=False)
    h.collectResultPerIngredient()
    todump=[[],[]]
    #we can save only the object positions
    todump[0] = [np.array(m[0]).tolist() for m in h.molecules]# [[ingr.nbMol,len(res)],[]]  POS
    with open(filename, 'w') as fp :
        json.dump(todump,fp)
    #    h.store_asJson(resultfilename=filename+"_results.json")
    #we can specify the information we want in the result file
    #h.saveRecipe(filename+"_results.json",useXref=False,mixed=True,result=True,
    #             kwds=["radii"],grid=False,packing_options=False,indent=False,quaternion=True)
    h.saveRecipe(filename+"_results.json",useXref=False,mixed=True,
                             kwds=["source","name","positions","radii"],result=True,
                           grid=False,packing_options=False,indent=False,quaternion=True) 
                           
def one_exp(h,seed,eid=0,setn=1,periodicity=True,output=None):
    clear(h)
    if h.use_periodicity:
        if output is None :
            output = os.getcwd() + "/dataERROR"
    else :
        if output is None :
            output = os.getcwd() + "/dataERROR"
    if not os.path.exists(output):
        os.makedirs(output)
    setn=str(setn).replace("(","").replace(")","").replace(", ","_")
    pack(h,seed,output+os.sep+"run_"+str(setn)+"_"+str(eid),eid)
    #plotOneResult2D(h,output+os.sep+"run_"+str(setn)+"_"+str(eid)+"_figure.png")
#    plotOneResult3D(h,output+os.sep+"run_"+str(setn)+"_"+str(eid)+"_figure")
    #SHOULD CHECK FOR 2D PACK...IF ONE AXIS OF THE BOUNDING BOXE IS 1
#    pandaPlot(h,output+os.sep+"run_"+str(setn)+"_"+str(eid)+"_figure")
    #plotIngredientResults2D(h,output+os.sep+"run_"+str(setn)+"_"+str(eid)+"_density_")
#    distanceMatrixOverlap(h,output+os.sep+"run_"+str(setn)+"_"+str(eid)+"_distances",avg=True)

#eid = seed-id
#setn = parameterSet

##================
## change the defaults of the recipe
##================
def applyGeneralOptionsChangeDefaults(env,parameter_set):
    for k in parameter_set :
        setattr(env,k,parameter_set[k])

def applyGeneralIngredientsOptionsChangeDefaults(env,parameter_set):
    #apply the key options to all ingredients
    r = env.exteriorRecipe
    if r :
        for ingr in r.ingredients:
            for k in parameter_set :
                setattr(ingr,k,parameter_set[k])

def applyIndividualIngredientsOptionsChangeDefaults(env,parameter_liste):
    #apply the key options to specified ingredients
    r = env.exteriorRecipe
    if r :
        for i,params in enumerate(parameter_liste) :
            for k in params :
                setattr(r.ingredients[i],k,params[k])

def applyIndividualIngredientsOptionsDictChangeDefaults(env,parameter_dict):
    #apply the key options to specified ingredients
    r = env.exteriorRecipe
    if r :
        for iname in parameter_dict :
            ingr = env.getIngrFromName(iname)
            print ingr
            if ingr is None :
                continue
            for k in parameter_dict[iname] :
                setattr(ingr,k,parameter_dict[iname][k])
                if k=="nbMol" :
                    ingr.overwrite_nbMol_value = ingr.nbMol


##================
## change the values for sampling
##================
def applyGeneralOptions(env,parameter_set,nset):
    for k in parameter_set :
        setattr(env,k,parameter_set[k][nset])


def applyGeneralIngredientsOptions(env,parameter_set,nset):
    #apply the key options to all ingredients
    r = env.exteriorRecipe
    if r :
        for ingr in r.ingredients:
            for k in parameter_set :
                setattr(ingr,k,parameter_set[k][nset])

def applyIndividualIngredientsOptions(env,parameter_liste,nset):
    #apply the key options to specified ingredients
    r = env.exteriorRecipe
    if r :
        for i,params in enumerate(parameter_liste) :
            for k in params :
                setattr(r.ingredients[i],k,params[k][nset])


def applyIndividualIngredientsOptionsDict(env,parameter_dict,nset):
    #apply the key options to specified ingredients
    r = env.exteriorRecipe
    if r :
        for iname in parameter_dict :
            ingr = env.getIngrFromName(iname)
            if ingr is None :
                continue
            for k in parameter_dict[iname] :
                setattr(ingr,k,parameter_dict[iname][k][nset])
                if k=="nbMol" :
                    ingr.overwrite_nbMol_value = ingr.nbMol

#==============================================================================
# autoPACK setup, and recipe loading
#==============================================================================
#define wher is out recipe setup file
#setupfile = "/Users/ludo/DEV/autopack_git/autoPACK_database_1.0.0/recipes/NM_Analysis_FigureB1.1.xml"
#force downloading the latest recipe file
autopack.forceFetch=True
#Or we can also use the server for gathering the file
#you can pass directly the recipe as an argument
#for isntance python -i analyze_B2.py /opt/data/dev/cellPACK/cellPACK_data/cellPACK_database_1.1.0/recipes/NM_Analysis_FigureB1.1.json
#you can also pass directly the recipe name as well as the version of the recipe
#for instance python -i analyze_B2.py NM_Analysis_FigureB 1.0 or for instance python -i analyze_B2.py NM_Analysis_FigureB 1.1
#you can also ahrd code the recipe file name or name + version if no argument are passed

# recipe = 'https://raw.githubusercontent.com/mesoscope/cellPACK_data/master/cellPACK_database_1.1.0/recipes/NM_Analysis_FigureB1.1.json'
recipename =parVals['running_settings']['recipeName']
version =parVals['running_settings']['recipeVersion']
recipe =parVals['running_settings']['recipePath']


if recipe != '':
    #remove 1.1.json from recipe
    setupfile = recipe
    filename = setupfile
else : # don't do it as else
    #Or we can also use the server for gathering the file
    recipe = recipename # put recipe here
    version = version # put version here
    filename = autopack.RECIPES[recipe][version]["setupfile"]
    setupfile = autopack.retrieveFile(filename,cache="recipes")

fileName, fileExtension = os.path.splitext(setupfile)
n=os.path.basename(fileName)
h = Environment(name=n)
h.loadRecipe(setupfile)
h.dump = False
afviewer=None
#h.placeMethod="pandaBullet"
#h.encapsulatingGrid=0
# default for periodicity
# h.use_periodicity = False
# autopack.testPeriodicity = True
# autopack.biasedPeriodicity = [1,1,0]
analyse = AnalyseAP(env=h, viewer=afviewer, result_file=None)
h.analyse = analyse
analyse.g.Resolution = 1.0
#h.smallestProteinSize=30.0#get it faster? same result ?
h.boundingBox=np.array(h.boundingBox)

#==============================================================================
# the parameter for a run
#==============================================================================
#all option can be access through a Dictionary that hold :
#options name
#options default value
#options value range
#options type e.g. Vector,Float, Int,String
#to get the value for the current value we can use the python function getattr
#to change the value of the option we can use the python function setattr
#all the packing options (General Parameters) are found in h.OPTIONS
#for k in h.OPTIONS.keys():
#    print h.OPTIONS[k]
#    print "\n"
#the main options to explore are probably:
#smallestProteinSize float
#pickWeightedIngr    boolean
#pickRandPt          boolean
#use_periodicity     boolean

#all the ingredients/Object options (Individual Parameters) are found in Ingredients.KWDS
#for k in ingredients_parameters.keys():
#    print k,ingredients_parameters[k]
#    print "\n"
#for instance packingPriority:
# - an optional packing priority. If omited the priority will be based
#    on the radius with larger radii first
#    ham here: (-)packingPriority object will pack from high to low one at a time
#    (+)packingPriority will be weighted by assigned priority value
#    (0)packignPriority will be weighted by complexity and appended to what is left
#    of the (+) values

#You should here prepare your parameter set I guess
#one set should be a dictionary for the packing
#one or several dictionary for the ingredients
#Comment uncomment here
packing_parameter_set=OrderedDict({
                                  #"smallestProteinSize":[5.0,10.0,15.0],
                                  "pickWeightedIngr":[1,0],
                                  "pickRandPt":[1,0],
                                  "use_periodicity":[1,0]
                                  })
#ingredients options you want to change for all objects
ingredients_parameter_set=OrderedDict({
                                      'rejectionThreshold':[10,20,100],
                                      #'packingPriority':[-10,0,10],
                                      #'cutoff_boundary':[0,10,20],
                                      'jitterMax':[[0.1,0.1,0.0],[1.0,1.0,0.0],[1.,1.,0.]],
                                      'nbJitter':[10,20,30],
                                      })


#individuals options. This recipe has 5objetcs, thus you need 5 set of values
#use a dictionary or a list
individuals_ingredients_parameter_set=OrderedDict({
                                                  "ext__Sphere_radius_100":OrderedDict({
                                                                                       'rejectionThreshold':[10,20,30],
                                                                                       'packingPriority':[-10,0,10],
                                                                                       'cutoff_boundary':[0,10,20],
                                                                                       'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
                                                                                       'nbJitter':[5,10,15],
                                                                                       'nbMol':[2,4,6,8],
                                                                                       }),
                                                  "ext__Sphere_radius_200":OrderedDict({'rejectionThreshold':[10,20,30],
                                                                                       'packingPriority':[-10,0,10],
                                                                                       'cutoff_boundary':[0,10,20],
                                                                                       'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
                                                                                       'nbJitter':[5,10,15],
                                                                                       'nbMol':[2,4,6,8],}),
                                                  "ext__Sphere_radius_50":OrderedDict({'rejectionThreshold':[10,20,30],
                                                                                      'packingPriority':[-10,0,10],
                                                                                      'cutoff_boundary':[0,10,20],
                                                                                      'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
                                                                                      'nbJitter':[5,10,15],
                                                                                      'nbMol':[10,14,16,18],}),
                                                  "ext__Sphere_radius_25p":OrderedDict({'rejectionThreshold':[10,20,30],
                                                                                       'packingPriority':[-10,0,10],
                                                                                       'cutoff_boundary':[0,10,20],
                                                                                       'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
                                                                                       'nbJitter':[5,10,15],
                                                                                       'nbMol':[20,40,60,80],}),
                                                  "ext__Sphere_radius_25r":OrderedDict({'rejectionThreshold':[10,20,30],
                                                                                       'packingPriority':[-10,0,10],
                                                                                       'cutoff_boundary':[0,10,20],
                                                                                       'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
                                                                                       'nbJitter':[5,10,15],
                                                                                       'nbMol':[20,40,60,80],}),
                                                  })
individuals_ingredients_parameter_liste=[
                                         OrderedDict({
                                                     'rejectionThreshold':[10,20,30],
                                                     'packingPriority':[-10,0,10],
                                                     'cutoff_boundary':[0,10,20],
                                                     'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
                                                     'nbJitter':[5,10,15],
                                                     'nbMol':[2,4,6,8],
                                                     }),
                                         OrderedDict({'rejectionThreshold':[10,20,30],
                                                     'packingPriority':[-10,0,10],
                                                     'cutoff_boundary':[0,10,20],
                                                     'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
                                                     'nbJitter':[5,10,15],
                                                     'nbMol':[2,4,6,8],}),
                                         OrderedDict({'rejectionThreshold':[10,20,30],
                                                     'packingPriority':[-10,0,10],
                                                     'cutoff_boundary':[0,10,20],
                                                     'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
                                                     'nbJitter':[5,10,15],
                                                     'nbMol':[5,10,15],}),
                                         OrderedDict({'rejectionThreshold':[10,20,30],
                                                     'packingPriority':[-10,0,10],
                                                     'cutoff_boundary':[0,10,20],
                                                     'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
                                                     'nbJitter':[5,10,15],
                                                     'nbMol':[20,40,60,80],}),
                                         OrderedDict({'rejectionThreshold':[10,20,30],
                                                     'packingPriority':[-10,0,10],
                                                     'cutoff_boundary':[0,10,20],
                                                     'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
                                                     'nbJitter':[5,10,15],
                                                     'nbMol':[20,40,60,80],}),
                                         ]

#there is a lot of different to manipulate the options, I just show you the dictionary way.
#you can also comment/uncomment the option you don't want
#
#==============================================================================
# experiment running
#==============================================================================

time_benchmark_file=outputFile + "/time_check.csv"
timer=0.0
time_data="run,seedID,seedValue,time\n";

#define an ouput directory
output=outputFile
#numberofRun per numberofSet containing numberofParameter
numberofParameter=len(packing_parameter_set)+len(ingredients_parameter_set)
# numberofSet=2 #because of the binary option, its correspondent to the possible number for the options,  e.g. True/False or 1,5 etc..
# numberofRun=1 #the number of seeds to use for one set, one set is a combination of options


ingredientNames = []
for ingr in h.exteriorRecipe.ingredients:
    ingredientNames.append(ingr.name)

if parVals:
    packing_parameter_set = {}#parVals['settings']['packing_parameter_set']
    packing_parameter_set_active = parVals['packing_parameter_set_active']
    
    ingredients_parameter_set = {}#parVals['settings']['ingredients_parameter_set']
    ingredients_parameter_set_active = {}#ssdfsdfsdf parVals['settings']['ingredients_parameter_set_active']

    individuals_ingredients_parameter_set = {}#parVals['settings']['individuals_ingredients_parameter_set']
    individuals_ingredients_parameter_set_active = parVals['individuals_ingredients_parameter_set_active']
    
    numberofSet=parVals['running_settings']['numberofSet']
    numberofRun=parVals['running_settings']['numberofRun']

#prepapre an array of 200 or more seed
nseed = 200
if nseed < numberofRun :
    nseed = numberofRun
seeds_i = analyse.getHaltonUnique(nseed)

# applyGeneralOptionsChangeDefaults(h,packing_parameter_set)
# applyGeneralIngredientsOptionsChangeDefaults(h,ingredients_parameter_set)
# applyIndividualIngredientsOptionsChangeDefaults(h,individuals_ingredients_parameter_liste)
# applyIndividualIngredientsOptionsDictChangeDefaults(h,individuals_ingredients_parameter_set)

#how can we restart at a particular run
restart=True
restart_pset = 45

for pset in range(numberofSet):
    if restart :
        if pset < restart_pset :
            continue
    count=0
    # apply the options for the given set
    applyGeneralOptions(h,packing_parameter_set_active,pset)
    applyGeneralIngredientsOptions(h,ingredients_parameter_set_active,pset)
    #applyIndividualIngredientsOptions(h,individuals_ingredients_parameter_liste_active,pset)
    applyIndividualIngredientsOptionsDict(h,individuals_ingredients_parameter_set_active,pset)
    
    #do the X run
    for seed in seeds_i[:numberofRun] :
        # print "seed",seed
        timer=time.time()
        one_exp(h,seed,eid=count,setn=pset,periodicity=False,output=output)
        count+=1
        dT=time.time()-timer
        time_data+=str(pset)+","+str(count)+","+str(seed)+","+str(dT)+"\n"
    #plotDensitiesPerRun(h,pset,numberofRun,output, output+os.sep+"run_" + str(pset) + "_density_")

f=open(time_benchmark_file,"w")
f.write(time_data)
f.close()
#execfile("analyze_B2.py")