# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:59:01 2013

@author: ludo
Hoe to run autopack from xml file
"""
# Show the effect of garbage collection
import sys
import os
import numpy
from numpy import matrix
import json
import math
from scipy import spatial
from random import randrange,seed,random
sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/")
sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/PIL/")

import autopack
#wrkDir = AutoFill.__path__[0]
localdir = wrkDir = autopack.__path__[0]
import numpy


helper = autopack.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper = helperClass(vi="nogui")
autopack.helper = helper

from autopack.Environment import Environment
from autopack import transformation as tr 

import numpy as np

filename = "/home/ludo/hivexp/HIVimature1.0.json"
fileName, fileExtension = os.path.splitext(filename)
n=os.path.basename(fileName)
h = Environment(name=n)
#h.helper = helper
recipe=n
h.loadRecipe(filename)
h.smallestProteinSize = 40

resultfilename = h.resultfile = "/home/ludo/hivexp/hivfull"
previousresult = "/home/ludo/hivexp/HIVimature_results.json"
r=h.loadResult("/home/ludo/hivexp/HIVimature_results.json",transpose=False)#
#ingredients = h.restore(*r)

#env is 0:6
#if need to offset ENV thats here.
#gag is 6:
#every 4th gag try to place gag_pol
#could try to rotate around a point
mrot=[]
angle=0
for i in range(6):
    m=tr.rotation_matrix(math.radians(angle), [0,0,1], point=[5.841,-1.458,133.803])    
    mrot.append(m)
    angle+=60
#tr.rotation_matrix(math.radians(60), [0,0,1], point=[5.841,-1.458,-133.803])
#tr.rotation_matrix(math.radians(120), [0,0,1], point=[5.841,-1.458,-133.803])
#tr.rotation_matrix(math.radians(180), [0,0,1], point=[5.841,-1.458,-133.803])
#tr.rotation_matrix(math.radians(240), [0,0,1], point=[5.841,-1.458,-133.803])
#tr.rotation_matrix(math.radians(300), [0,0,1], point=[5.841,-1.458,-133.803])


def test(d,i,ingr,all_pos):
    jtrans,rotMatj,n,ptid,com = d[i]
    trans = numpy.array([4.282,-9.147,-76.946])
    if not len(all_pos):
        jtrans = jtrans + helper.ApplyMatrix([trans],rotMatj)[0]
        centT = ingr.transformPoints(jtrans, rotMatj, ingr.positions[-1])
        return True,jtrans,rotMatj
#    mat = tr.rotation_matrix(math.radians(-45),[0,0,1]).transpose()
#    trans = numpy.array([0,0,-74])
    #or random?
    for i in range(6):
        rotMatj = numpy.array(matrix(rotMatj)*matrix(mrot[i]))
        jtrans = jtrans + helper.ApplyMatrix([trans],rotMatj)[0]
        centT = ingr.transformPoints(jtrans, rotMatj, ingr.positions[-1])
        tree= spatial.cKDTree(all_pos, leafsize=10)
        distance,nb = tree.query(centT,len(all_pos))#,distance_upper_bound=cutoff)
        del tree
        #see if collide with anything ?
        #look only at the last one 
        res = distance[-1]
    #    for j,res in enumerate(distance):
        v=res.min()
    #    vi=numpy.where(res==v)
    #    print "found",res
        if v < ingr.radii[0][-1]*2:
    #        print "too close reject",v,ingr.radii[0][-1],i
            return False,jtrans,rotMatj
        else :
            return True,jtrans,rotMatj
    return True,jtrans,rotMatj


gag = h.getIngrFromName("GAG")
ingr = gag_pol = h.getIngrFromName("hiv_gag_pol")
env = h.getIngrFromName("iSUTM")

all_pos=[]
ngagpol=134
data = r[1][0][6:]
respol=[]
safety=0
mat = tr.rotation_matrix(math.radians(-45),[0,0,1]).transpose()
trans = numpy.array([6,12,-74])
trans = numpy.array([0,0,-74])
trans = numpy.array([4.282,-9.147,-76.946])
while ngagpol >=0:
    safety+=1
    if safety > 135 :
        break
    for i in range(len(data)):
#        print i,safety
        res,jtrans,rotMatj = test(data,i,gag_pol,all_pos)
        if res :#test(data,i,gag_pol,all_pos) :
            t,rot,n,ptid,com = data[i]
            #modify jtrans and rotMatj
#            rotMatj = numpy.array(matrix(rotMatj)*matrix(mat))
#            jtrans = jtrans + helper.ApplyMatrix([trans],rotMatj)[0]
            centT = ingr.transformPoints(jtrans, rotMatj, ingr.positions[-1])
            all_pos.extend(centT)
            respol.append([jtrans,rotMatj,"hiv_gag_pol",ptid,com])
            data.pop(i)
            ngagpol-=1
            print "place pol",len(respol),i,safety
            i=0
            break
        


envdata = r[1][0][:6]
alldata = []
alldata.extend(envdata)
alldata.extend(data)
alldata.extend(respol)
r[1][0] = alldata
for d in r[1][0]:
    d[1] = d[1].transpose()
ingredients = h.restore(*r)
h.saveRecipe("hivimmature_start_pol.json",useXref=True,mixed=True,
                     kwds=["source","name"],result=True,
                   grid=False,packing_options=False,indent=False,quaternion=True)
#transpose


#
#h.buildGrid(boundingBox=h.boundingBox,gridFileIn=None,rebuild=True ,
#                              gridFileOut="/home/ludo/hivexp/gridOut",previousFill=True)
#
#gag.completion = 1.0
#env.completion = 1.0
##rest should pack around
#gag.nbMol = 0
#env.nbMol = 0
#
#gag.molarity = 0
#env.molarity = 0
#
#h.fill5(verbose = 0,usePP=False)
#h.collectResultPerIngredient()
#rname="/home/ludo/hivexp/hiv_immature_inside.json"
##        rname="C:\Users\ludovic\Documents\CellPackViewer_Cluster\Data\HIV\cellPACK\mycoDNA1.json"
#h.saveRecipe(rname,useXref=True,mixed=True,
#                     kwds=["source","name"],result=True,
#                   grid=False,packing_options=False,indent=False,quaternion=True)
#
#h.saveGridToFile("/home/ludo/hivexp/grid_after_packing")
#transpose ?