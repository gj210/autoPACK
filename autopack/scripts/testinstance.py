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

#f=open("/Users/ludo/Library/Application Support/autoPACK/cache_results/HIV_16_result.json","r")
f=open("/Users/ludo/Library/Application Support/autoPACK/cache_results/HIV-1_0.1.6_3.apr.json","r")
#f=open("/home/ludo/.autoPACK/cache_results/HIV_16_result.json","r")
#f=open("/home/ludo/.autoPACK/cache_results/HIV-1_0.1.6_3.apr.json","r")
data = json.load(f)
f.close()

    
Helper = upy.getHelperClass()
helper = Helper()

def oneMat(tr,rot):
    if helper.instance_dupliFace :
        rot = numpy.array(rot).transpose()
        rot[3][:3]=tr
        return rot
    else :
        rot = numpy.array(rot)
        rot[:3,3]=tr
        return rot
        
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

def allQuad():
    for a in ["+X","-X","+Y","-Y","+Z","-Z"]:
        v=numpy.array(helper.quad[a])*50.0
        f=[[0,1,2,3]]
        me=helper.createsNmesh(a[1]+a[0],v,[],f)
        
def testOneIngredient(name,iname,oname,principalVector):
#    helper.read("/home/ludo/.autoPACK/cache_geometries/"+name+".dae")
    helper.read("/Users/ludo/Library/Application Support/autoPACK/cache_geometries/"+name+".dae")
    obj = helper.getObject(iname)
    print (obj,name,iname)
    #blender
    axis=numpy.array(principalVector[:])#helper.ApplyMatrix([principalVector,],mY) [0]
    mY=helper.rotation_matrix(-math.pi/2.0,[0.0,1.0,0.0])   
    if helper.host.find("blender") !=-1:
#        mA=helper.rotation_matrix(-math.pi/2.0,[1.0,0.0,0.0])
#        mB=helper.rotation_matrix(math.pi/2.0,[0.0,0.0,1.0])    
#        print ("principalVector b",principalVector)
#        #axis = helper.rotatePoint(principalVector,[0.,0.,0.],[0.0,1.0,0.0,-math.pi/2.0])
#        #print ("principalVector b",principalVector,axis)
#        m=matrix(mA)*matrix(mB)
#        helper.setObjectMatrix(obj,matrice=m)
        helper.resetTransformation(obj)
#    axis=helper.ApplyMatrix([principalVector,],mY) [0]
    elif helper.host=="maya":
         helper.resetTransformation(obj)
#        helper.rotateObj(obj,[0.0,-math.pi/2.0,0.0])
#    print ("rotated ",helper.ApplyMatrix([principalVector,],mY) [0])
    inst1 = helper.newInstance(oname,obj)
#    print (iname)
#    print (data['compartments']['HIV1_envelope_Pack_145_0_2_0']['surface']['ingredients'].keys())
    res2=numpy.array(data['compartments']['HIV1_envelope_Pack_145_0_2_0']['surface']['ingredients'][iname]["results"])
    listM=[oneMat(d[0],d[1]) for d in res2]
    print ("principalVector",axis)
    helper.instancePolygon("instOf"+oname, matrices=listM, mesh=inst1,
                                   axis=axis,transpose=True)
    
#allQuad()
#testOneIngredient("HIV1_MA_Hyb_0_1_0","HIV1_MA_Hyb_0_1_0","MA",[0,1,0])     #=> Z 90
helper.instance_dupliFace = True
testOneIngredient("HIV1_ENV_4nco_0_1_2_Left","HIV1_ENV_4nco_0_1_1","ENV",[0,0,-1]) #=> X 90
#helper.instance_dupliFace = False
#testOneIngredient("HIV1_ENV_4nco_0_1_2_Left","HIV1_ENV_4nco_0_1_1","ENV",[0,0,-1]) #=> X 90

#makeQuad([0,0,-1]) 
testOneIngredient("HIV1_NEF_Hyb_0_2_0","HIV1_NEF_Hyb_0_2_0","Nef",[0,0,1]) 
#testOneIngredient("HIV1_ENV_4nco_0_1_1","ENV",[0,0,-1]) 
testOneIngredient("HIV1_HLA_1dlh_0_1_0","HIV1_HLA_1dlh_0_1_0","hla",[0,0,-1]) 

##execfile("/Users/ludo/DEV/autoPACK_github/autopack/scripts/testinstance.py")
##execfile("/opt/data/dev/autoPACK/autopack/scripts/testinstance.py")
#exec(open("/opt/data/dev/autoPACK/autopack/scripts/testinstance.py").read())