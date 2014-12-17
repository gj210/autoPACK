# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 11:47:01 2014

@author: ludo
"""
MAYA=0
import sys
if MAYA:
    sys.path.append("/Users/ludo/Library/Preferences/Autodesk/maya/2015-x64/plug-ins/MGLToolsPckgs")
    sys.path.append("/Users/ludo/Library/Preferences/Autodesk/maya/2015-x64/plug-ins/MGLToolsPckgs/PIL")
    #maya standalone special
    import maya.standalone
    maya.standalone.initialize()
    #load plugin
    import maya
    maya.cmds.loadPlugin("fbxmaya")
    
import upy
import math
from random import random
from numpy import matrix
import numpy
try :
    import simplejson as json
    from simplejson import encoder    
except :
    import json
    from json import encoder

f=open("/Users/ludo/Library/Application Support/autoPACK/cache_results/HIV_16_result.json","r")
#f=open("/home/ludo/.autoPACK/cache_results/HIV_16_result.json","r")
data = json.load(f)
f.close()

def oneMat(tr,rot):
    rot = numpy.array(rot).transpose()
    rot[3][:3]=tr
    return rot
    
Helper = upy.getHelperClass()
helper = Helper()
helper.instance_dupliFace = True

def makeQuad(vector):
    mY=helper.rotation_matrix(-math.pi/2.0,[0.0,1.0,0.0])
    axis=helper.ApplyMatrix([vector,],mY) [0]
    axis = vector#helper.rotatePoint(vector,[0.,0.,0.],[0.0,1.0,0.0,-math.pi/2.0])    
    av=[int(axis[0]),int(axis[1]),int(axis[2])]
#    av=[int(axis[2]),int(axis[1]),int(axis[0])]
    axe=helper.rerieveAxis(av)
    print ("found ",axe,av)
    #axe="+Y"
    i=0
    f=[]
    quad=numpy.array(helper.quad[axe])*10.0
    f.append([i*4+0,i*4+1,i*4+2,i*4+3])
    q=helper.createsNmesh("quad"+str(axe),quad,[],f,color=[1,0,0])
    return q
    
def testOneIngredient(name,oname,principalVector):
    helper.read("/Users/ludo/Library/Application Support/autoPACK/cache_geometries/"+name+".dae")
    obj = helper.getObject(name)
    print (obj,name)
    #blender
    axis=numpy.array(principalVector[:])#helper.ApplyMatrix([principalVector,],mY) [0]
    mY=helper.rotation_matrix(-math.pi/2.0,[0.0,1.0,0.0])   
    if helper.host.find("blender") !=-1:
        mA=helper.rotation_matrix(-math.pi/2.0,[1.0,0.0,0.0])
        mB=helper.rotation_matrix(math.pi/2.0,[0.0,0.0,1.0])    
        print ("principalVector b",principalVector)
        #axis = helper.rotatePoint(principalVector,[0.,0.,0.],[0.0,1.0,0.0,-math.pi/2.0])
        #print ("principalVector b",principalVector,axis)
        m=matrix(mA)*matrix(mB)
        helper.setObjectMatrix(obj,matrice=m)
#    axis=helper.ApplyMatrix([principalVector,],mY) [0]
    elif helper.host=="maya":
        helper.rotateObj(obj,[0.0,-math.pi/2.0,0.0])
    print ("rotated ",helper.ApplyMatrix([principalVector,],mY) [0])
    inst1 = helper.newInstance(oname,obj)
    res2=numpy.array(data['HIV1_envelope_Pack_145_0_2_0_surfaceRecipe']["HIV1_envelope_Pack_145_0_2_0_surf__"+name]["results"])
    listM=[oneMat(d[0],d[1]) for d in res2]
    print ("principalVector",axis)
    helper.instancePolygon("instOf"+oname, matrices=listM, mesh=inst1,
                                   axis=axis,transpose=True)
    

testOneIngredient("HIV1_MA_Hyb_0_1_0","MA",[0,1,0])     #=> Z 90
testOneIngredient("HIV1_ENV_4nco_0_1_1","ENV",[0,0,-1]) #=> X 90
#makeQuad([0,0,-1]) 
testOneIngredient("HIV1_NEF_Hyb_0_2_0","Nef",[0,0,1]) 
#testOneIngredient("HIV1_ENV_4nco_0_1_1","ENV",[0,0,-1]) 
testOneIngredient("HIV1_HLA_1dlh_0_1_0","hla",[0,0,-1]) 

##execfile("/Users/ludo/DEV/autoPACK_github/autopack/scripts/testinstance.py")