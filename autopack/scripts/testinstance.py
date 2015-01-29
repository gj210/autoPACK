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
#f=open("/Users/ludo/Library/Application Support/autoPACK/cache_results/HIV-1_0.1.6_3.apr.json","r")
#f=open("/home/ludo/.autoPACK/cache_results/HIV_16_result.json","r")
f=open("/home/ludo/.autoPACK/cache_results/HIV-1_0.1.6_3.apr.json","r")
data = json.load(f)
f.close()
f=open("/home/ludo/.autoPACK/cache_recipes/HIV-1_0.1.6-6.json","r")
recipe = json.load(f)
f.close()

    
#Helper = upy.getHelperClass()
helper = upy.getHelperClass()()
#r=767.0
#need quarter up : cube 
def up (r):
    pos = numpy.array([r/2,r/2,r/2] )
    size = numpy.array([r,r,r] )      
    bb=numpy.array([[pos-size/2.0],[pos+size/2.0]])
    return pos,size,bb

#need half up    : cube pos = [0,r/2,r/2.] size = [r,r,r*2.0]    bb=[[pos-size],[pos+size]]
def halfup (r):
    pos = numpy.array([0,r/2,r/2.])
    size = numpy.array([r,r,r*2.0] )      
    bb=numpy.array([[pos-size/2.0],[pos+size/2.0]])
    return pos,size,bb

#need half       : cube pos = [0,0,r/2.] size = [r,r*2.0,r*2.0]  bb=[[pos-size],[pos+size]]
def half (r):
    pos = numpy.array([0,0,r/2.])
    size =  numpy.array([r,r*2.0,r*2.0] )     
    bb=numpy.array([[pos-size/2.0],[pos+size/2.0]])
    return pos,size,bb

def testOnePoints(p,bbi):
     res=numpy.logical_and(numpy.greater(p,bbi[0]),numpy.less(p,bbi[1]))
     if False in res :
         return False
     else :
         return True
    
def getPoints(bbi,datai):
    res=numpy.logical_and(numpy.greater(datai,bbi[0]),numpy.less(datai,bbi[1]))
    c=numpy.average(res,1)
    #d=numpy.equal(c,1.)
    #ptinside = numpy.nonzero(d)[0]
    d=numpy.less(c,1.)
    ptoutside = numpy.nonzero(d)[0]
    data_out = numpy.take(datai,ptoutside,0)
    return data_out
    
def oneMat(tr,rot,mask=None):
    if type(mask) != None :
        if testOnePoints(tr,mask) :
            return None
    rot = numpy.array(rot).transpose()
    rot[3][:3]=tr
    return rot
    
def oneMat2(tr,rot,mask=None):
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
        
def testOneIngredient(name,iname,oname,principalVector,comp,cname="surface",mask=None):
#    helper.read("/home/ludo/.autoPACK/cache_geometries/"+name+".dae")
    obj = helper.getObject(iname)
    if obj is None :
#        helper.read("/Users/ludo/Library/Application Support/autoPACK/cache_geometries/"+name+".dae")
        helper.read("/home/ludo/.autoPACK/cache_geometries/"+name+".dae")
        obj = helper.getObject(iname)
    print (obj,name,iname)
    #blender
    axis=numpy.array(principalVector[:],'f')#helper.ApplyMatrix([principalVector,],mY) [0]
#    mY=helper.rotation_matrix(-math.pi/2.0,[0.0,1.0,0.0])   
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
    inst1 = helper.getObject(oname)
    if inst1 is None :
        inst1 = helper.newInstance(oname,obj)
#    print (iname)
#    print (data['compartments']['HIV1_envelope_Pack_145_0_2_0']['surface']['ingredients'].keys())
    res2=numpy.array(data['compartments'][comp][cname]['ingredients'][iname]["results"])
    if helper.instance_dupliFace :
        listM=[]
        for d in res2 :
            m = oneMat(d[0],d[1],mask=mask)
            if type(m) != None and m is not None :
                listM.append(m)
#                print ("ok ",m)
            else :
#                print (" ok NONENONEEN")
                continue
#        print (listM)
#        listM=[oneMat(d[0],d[1],mask=mask) for d in res2]
    else :
        listM=[oneMat2(d[0],d[1],mask=mask) for d in res2]   
    print ("principalVector",axis)
    helper.instancePolygon("instOf"+oname, matrices=listM, mesh=inst1,
                                   axis=axis,transpose=True)
#helper.quad={"+Z" : [ [-1,-1,0],[-1,1,0],[1,1,0],[1,-1,0],],#XY
#               "+Y" :[[-1,0,-1],[-1,0,1],[1,0,1],[1,0,-1] ],#XZ
#               "-X" :[[0,-1,1],[0,1,1],[0,1,-1], [0,-1,-1]],#YZ
#               "-Z" :[[-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0]],#XY
#               "-Y" :[[-1,0,1],[1,0,1],[1,0,-1], [-1,0,-1]],#XZ
#               "+X" :[[0,-1,1],[0,1,1],[0,1,-1], [0,-1,-1]],#YZ
#           } 
#
#helper.eq={"+X":[helper.track_axis_dic["+Y"][1],-math.pi/2.0],
#            "+Y":[helper.track_axis_dic["-X"][1],math.pi/2.0],
#            "+Z":[helper.track_axis_dic["-Z"][1],0.0],
#            "-X":[helper.track_axis_dic["-Y"][1],-math.pi/2.0],
#            "-Y":[helper.track_axis_dic["-X"][1],math.pi/2.0],
#            "-Z":[helper.track_axis_dic["-Z"][1],-math.pi/2.0]}
##allQuad()
#testOneIngredient("HIV1_MA_Hyb_0_1_0","HIV1_MA_Hyb_0_1_0","MA",[0,1,0])     #=> Z 90
helper.instance_dupliFace = True
p,s,bb=up (767.0) #used for lipids
pp,ss,bbsurface = up (700.0)
bbsurface = numpy.array([[p-ss/2.0],[p+ss/2.0]])
#decrease size
#surface mask is up recipe['compartments']['HIV1_envelope_Pack_145_0_2_0']['surface']['ingredients']
comp='HIV1_envelope_Pack_145_0_2_0'
for ingredients in recipe['compartments']['HIV1_envelope_Pack_145_0_2_0']['surface']['ingredients']:
    fname = recipe['compartments']['HIV1_envelope_Pack_145_0_2_0']['surface']['ingredients'][ingredients]['include']
    f=open("/home/ludo/.autoPACK/cache_recipes/"+fname,"r")
    ingrdic = json.load(f)
    f.close()
    fname=ingrdic['meshFile'].split("/")[-1][:-4]
    iname=recipe['compartments']['HIV1_envelope_Pack_145_0_2_0']['surface']['ingredients'][ingredients]['name']
    axe = ingrdic['principalVector']
    if 'meshName' in ingrdic :
        iname=ingrdic['meshName']
    else :
        print ("no meshName for ingredients",fname,iname)
    print ("instance of ",fname,iname)
#    testOneIngredient(fname,iname,iname+"s",axe,comp,cname="surface",mask=bbsurface)

#inside
pp,ss,bbmatrix = up (650.0)
bbmatrix = numpy.array([[p-ss/2.0],[p+ss/2.0]])
for ingredients in recipe['compartments']['HIV1_envelope_Pack_145_0_2_0']['interior']['ingredients']:
    fname = recipe['compartments']['HIV1_envelope_Pack_145_0_2_0']['interior']['ingredients'][ingredients]['include']
    f=open("/home/ludo/.autoPACK/cache_recipes/"+fname,"r")
    ingrdic = json.load(f)
    f.close()
    fname=ingrdic['meshFile'].split("/")[-1][:-4]
    iname=recipe['compartments']['HIV1_envelope_Pack_145_0_2_0']['interior']['ingredients'][ingredients]['name']
    axe = ingrdic['principalVector']
    if 'meshName' in ingrdic :
        iname=ingrdic['meshName']
    else :
        print ("no meshName for ingredients",fname,iname)
    print ("instance of ",fname,iname)
#    testOneIngredient(fname,iname,iname+"i",axe,comp,cname="interior",mask=bbmatrix)
#surface nucleus
pp,ss,bbnucleo = up (600.0)
bbnucleo = numpy.array([[p-ss/2.0],[p+ss/2.0]])
helper.read("/home/ludo/.autoPACK/cache_geometries/HIV1_capside_3j3q_Rep_Med_0_2_1.dae")
geom = helper.getObject("HIV1_capside_3j3q_Rep_Med")
helper.resetTransformation(geom)
#h=helper.getObject("Hexamers")
#p=helper.getObject("Pentamers")
#get all children and test position against bb
#for c in helper.getChilds(h):
#    t=helper.getTranslation(c)
#    if testOnePoints(t,bbnucleo):
#        helper.deleteObject(c)
#for c in helper.getChilds(p):
#    t=helper.getTranslation(c)
#    if testOnePoints(t,bbnucleo):
#        helper.deleteObject(c)
#remove all the instance inside the bb
#inside nucleus
comp='HIV1_capsid_3j3q_PackInner_0_1_0'
for ingredients in recipe['compartments']['HIV1_capsid_3j3q_PackInner_0_1_0']['interior']['ingredients']:
    fname = recipe['compartments']['HIV1_capsid_3j3q_PackInner_0_1_0']['interior']['ingredients'][ingredients]['include']
    f=open("/home/ludo/.autoPACK/cache_recipes/"+fname,"r")
    ingrdic = json.load(f)
    f.close()
    fname=ingrdic['meshFile'].split("/")[-1][:-4]
    iname=recipe['compartments']['HIV1_capsid_3j3q_PackInner_0_1_0']['interior']['ingredients'][ingredients]['name']
    axe = ingrdic['principalVector']
    if 'meshName' in ingrdic :
        iname=ingrdic['meshName']
    else :
        print ("no meshName for ingredients",fname,iname)
    print ("instance of ",fname,iname)
#    testOneIngredient(fname,iname,iname+"in",axe,comp,cname="interior",mask=bbnucleo)

#helper.instance_dupliFace = False
##testOneIngredient("HIV1_ENV_4nco_0_1_2_Left","HIV1_ENV_4nco_0_1_1","ENV",[0,0,-1]) #=> X 90
##testOneIngredient("HIV1_NEF_Hyb_0_2_1","HIV1_NEF_Hyb_0_2_0","Nef",[0,0,1]) 
#testOneIngredient("HIV1_MA_Hyb_0_1_1","HIV1_MA_Hyb_0_1_0","MA",[0,1,0])     #=> Z 90
#testOneIngredient("HIV1_CAhex_0_1_1","HIV1_CAhex_0_1_0","CA",[1,0,0],cname="interior")
#makeQuad([0,0,-1]) 
#testOneIngredient("HIV1_NEF_Hyb_0_2_0","HIV1_NEF_Hyb_0_2_0","Nef",[0,0,1]) 
#testOneIngredient("HIV1_ENV_4nco_0_1_1","ENV",[0,0,-1]) 
#testOneIngredient("HIV1_HLA_1dlh_0_1_0","HIV1_HLA_1dlh_0_1_0","hla",[0,0,-1]) 

##execfile("/Users/ludo/DEV/autoPACK_github/autopack/scripts/testinstance.py")
##execfile("/opt/data/dev/autoPACK/autopack/scripts/testinstance.py")
#exec(open("/opt/data/dev/autoPACK/autopack/scripts/testinstance.py").read())