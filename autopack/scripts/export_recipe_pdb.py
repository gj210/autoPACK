# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 22:08:01 2014

@author: ludo
"""
import sys
import os
import math

##The program will read and write all the packings in the following file format: It should be a binary file which stores sequentially sphere center x, y, z coordinates and diameter for each particle as floating points in double precision in little-endian byte order. If the machine on which the program is being run is big-endian, the program will detect it and will still read and write in little-endian format. 

import sys
#MGLTools import
sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs")
#Maxwell import
sys.path.append("/usr/local/maxwell-3.0//python/pymaxwell/python2.7_UCS4")
#sys.path.append("/usr/local/maxwell-3.0//python/pymaxwell/python2.7_UCS2")#
sys.path.append("/usr/local/maxwell-3.0/")
from pymaxwell import *

import prody
prody.confProDy(verbosity='none')


import numpy
import autopack
#wrkDir = AutoFill.__path__[0]
localdir = wrkDir = autopack.__path__[0]

from autopack.Environment import Environment

TWOD = 1
NOGUI = 1
ANALYSIS = 0
helper = autopack.helper
if helper is None and not NOGUI:
    import upy
    helperClass = upy.getHelperClass()
    helper =helperClass()
else :
    import upy
    helperClass = upy.getHelperClass()
    helper = helperClass(vi="nogui")
    print helper
autopack.helper = helper
autopack.fixpath = True


try :
    import simplejson as json
    from simplejson import encoder    
except :
    import json
    from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.8g')
cellPackserver="https://raw.githubusercontent.com/mesoscope/cellPACK_data/master/cellPACK_database_1.1.0"
try :
    from collections import OrderedDict
except :
    from ordereddict import OrderedDict

import urllib
import struct



CACHE_FOLDER = "/home/ludo/.autoPACK/cache_others/"
PATH="/opt/data/dev/cellPACK/cellPACK_data/cellPACK_database_1.1.0/recipes/"
vdwRadii = { 'N':1.55, 'C':1.70, 'O':1.52,
                  'H':1.20, 
                  'S':1.80}

def getCenter(coords):
    center = sum(coords)/(len(coords)*1.0)
    center = list(center)
    for i in range(3):
        center[i] = round(center[i], 4)
#        print "center =", center
    return numpy.array(center)

def getRadius(coords,center):
    dist = numpy.linalg.norm(coords-center,axis=1)
    return max(dist)


def fetch_pdb(pdbid):
  url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid.upper()
  return urllib.urlopen(url).read()

def fetch_cellPACK(pdbid):
  url = 'https://raw.githubusercontent.com/mesoscope/cellPACK_data/master/cellPACK_database_1.1.0/other/%s' % pdbid
  return urllib.urlopen(url).read()

from upy.transformation import decompose_matrix
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
         

def readOne(fname,gname,scene):
#    iscene = Cmaxwell(mwcallback);
    filename =autopack.retrieveFile(fname,cache="geoms") #geometries
    ok = scene.readMXS(str(filename));
    if not ok :
        print "problem reading",filename
    obj = scene.getObject(str(gname))[0];
#    print obj,gname
#    obj2 = scene.addObject(obj);
#    iscene.freeScene();
    return obj,scene

def create_mxs_from_particle_bins(binfile,scene=None):
    # Create scene.
    if scene is None :
        scene = Cmaxwell(mwcallback);
        camera = scene.addCamera('camera',1,1.0/800.0,0.04,0.035,400,"CIRCULAR",0,0,25,400,300,1);
        camera.setStep(0,Cvector(0.0,0.0,6.0),Cvector(0.0,0.0,0.0),Cvector(0.0,1.0,0.0),0.035,8.0,0);
        camera.setActive();
        # Add physical sky and sun
        scene.getEnvironment().setActiveSky('PHYSICAL')
        sunColor = Crgb()
        sunColor.assign(1,1,1)
        scene.getEnvironment().setSunProperties(SUN_PHYSICAL,5777,1,1,sunColor)
        
    # Create RFRK object.
    mgr = CextensionManager.instance();
    #mgr.setExtensionsDirectory('C:/Program Files/Next Limit/Maxwell 3/extensions/')
    mgr.loadAllExtensions()
    print mgr
    ext = mgr.createDefaultGeometryProceduralExtension('MaxwellParticles')
    print ext
    params = ext.getExtensionData()
    params.setString('FileName',binfile)
    #set the radius
#    params.getFloat("Radius Factor")
    params.setFloat('Radius Factor',120.000)
#    params.setDouble('RadiusMultiplier',120.000)
    pobject = scene.createGeometryProceduralObject(binfile[-3:],params)
    if pobject.isNull():
        print("Error creating particles object")
        return 0;
    return pobject
    


def buildCompartmentsGeom(comp,scene,parent=None):
    if comp.representation_file is None : 
        return
    nr=scene.createNullObject(str(comp.name)+str("rep"))
    if parent is not None :
        nr.setParent(parent)
    filename =autopack.retrieveFile(comp.representation_file,cache="geoms") #geometries    
    gdic=helper.read(filename)
    for nid in gdic :
        matnode=oneMaterial(str(nid),scene,color=gdic[nid]["color"])
        mxmesh = maxwellMesh(str(nid),gdic[nid]["mesh"][0],gdic[nid]["mesh"][1],gdic[nid]["mesh"][2],scene,matnode=matnode)
        mxmesh.setParent(nr)
        if len(gdic[nid]['instances']):
            nri=scene.createNullObject(str(nid)+str("instances"))
            nri.setParent(nr)
            c=0
            for mat in gdic[nid]['instances']:
                instance = scene.createInstancement(str(nid)+"_"+str(c),mxmesh)
                mat = numpy.array(mat,float).transpose()
                scale, shear, euler, translate, perspective=decompose_matrix(mat)
                base,pivot,ok = instance.getBaseAndPivot();
                base.origin = Cvector(float(mat[3][0]),float(mat[3][1]),float(mat[3][2]))
                mat[3,:3] = [0.,0.,0.]
#                mat = numpy.array(mat,float).transpose()
                newaxis=helper.ApplyMatrix([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],mat)
                
                base.xAxis=Cvector(float(newaxis[0][0]),float(newaxis[0][1]),float(newaxis[0][2]))
                base.yAxis=Cvector(float(newaxis[1][0]),float(newaxis[1][1]),float(newaxis[1][2]))
                base.zAxis=Cvector(float(newaxis[2][0]),float(newaxis[2][1]),float(newaxis[2][2]))
#                base.xAxis=Cvector(float(mat[0][0]),float(mat[1][0]),float(mat[2][0]))
#                base.yAxis=Cvector(float(mat[0][1]),float(mat[1][1]),float(mat[2][1]))
#                base.zAxis=Cvector(float(mat[0][2]),float(mat[1][2]),float(mat[2][2]))\
#                base.xAxis=Cvector(float(mat[0][0]),float(mat[0][1]),float(mat[0][2]))
#                base.yAxis=Cvector(float(mat[1][0]),float(mat[1][1]),float(mat[1][2]))
#                base.zAxis=Cvector(float(mat[2][0]),float(mat[2][1]),float(mat[2][2]))
                base.xAxis.normalize();
                base.yAxis.normalize();
                base.zAxis.normalize();

#                base.origin = Cvector(float(mat[0][3]),float(mat[1][3]),float(mat[2][3]))
#                print ("X",mat[0][0],mat[0][1],mat[0][2],base.xAxis.x(),base.xAxis.y(),base.xAxis.z())
#                print ("Y",mat[1][0],mat[1][1],mat[1][2],base.yAxis.x(),base.yAxis.y(),base.yAxis.z())
#                print ("Z",mat[2][0],mat[2][1],mat[2][2],base.zAxis.x(),base.zAxis.y(),base.zAxis.z())
#                print ("base",mat[3][0],mat[3][1],mat[3][2],base.origin.x(),base.origin.y(),base.origin.z())
                instance.setBaseAndPivot(base,base) 
#                pos=Cvector(float(mat[3][0]),float(mat[3][1]),float(mat[3][2]))
#                eul=Cvector(float(euler[0]),float(euler[1]),float(euler[2]))
#                instance.setPosition(pos)
#                instance.setRotation(eul)
#                instance.setParent(nri)
#                wt=instance.getWorldTransform()[0]
                print str(nid)+"_"+str(c),mat
#                print wt.origin.x(),wt.origin.y(),wt.origin.z(),float(mat[3][0]),float(mat[3][1]),float(mat[3][2])
#                print wt.xAxis.x(),wt.xAxis.y(),wt.xAxis.z()
#                print wt.yAxis.x(),wt.yAxis.y(),wt.yAxis.z()
#                print wt.zAxis.x(),wt.zAxis.y(),wt.zAxis.z()
#                print instance.getPosition()[0]
#                print instance.getRotation()[0]
#                base,pivot,ok = instance.getBaseAndPivot();
                if c==5 :
#                    print ("X",mat,base.xAxis.x(),base.xAxis.y(),base.xAxis.z())
                    break
                c+=1

def toBinary(filename,datai):
    f=open(filename,"wb")#"/home/ludo/Tools/lipidwrapper/bins/test_np3.pxy","wb")
    f.write(struct.pack('=f',3.0))                 #version
    f.write(struct.pack('=i',len(datai)))   #nb particle 
    for i in range(len(datai)):#len(a['arr_0'])):#len(a['arr_0'])):
        f.write(struct.pack('=Qffffff',int(i),datai[i][0],datai[i][1],datai[i][2],0.,0.,0.))
    f.close()
    return filename

def loadPDB(fname,scene=None,center=True,transform=None):
    mols=None
    try :
        mols = prody.parsePDB(fname)#display lines 
        print "load "+fname
    except :
        print "problem reading ",fname
        return None
    moveAtoms(mols, to=numpy.zeros(3))
    return mols

def fetchLoadPDB(pdbid,scene=None,transform=None):
    if pdbid.find("SwissModel") != -1 :
        return None
    if len(pdbid) > 4 :
        #actual file name?
        fname=CACHE_FOLDER+str(pdbid)
        if not os.path.isfile(CACHE_FOLDER+str(pdbid)):
            aStr=fetch_cellPACK(pdbid)
            if not os.path.isfile(fname) :
                f=open(fname,"w")
                f.write(aStr)
                f.close()
    else :
        pdbid=pdbid[:4].upper()
        print pdbid
        fname=CACHE_FOLDER+str(pdbid)+".pdb"
        if not os.path.isfile(CACHE_FOLDER+str(pdbid)+".pdb"):
            aStr=fetch_pdb(pdbid)
            if aStr.find("Not Found") != -1 or aStr.find("<html>") != -1 :
                return None,None      
            if not os.path.isfile(fname) :
                f=open(fname,"w")
                f.write(aStr)
                f.close()
    #read it in Pmv?
    print "load ",fname
#    mols=loadPDB(fname,scene=scene,transform=transform)
#        mols = prody.parsePDB(fname)
    return fname

def buildIngredient(ingr):
    #get the source    
    pdbid=ingr.source["pdb"]
    if pdbid is None or pdbid.find("EMDB")!= -1:
        return liste_pdb
#    if pdbid not in liste_pdb:
#         if not sphere :
    mols=fetchLoadPDB(pdbid)
    if mols is None :
        return None     
    return mols

def writeMol(name,mol,biomt):
    fout=CACHE_FOLDER+"/"+name+"_biomt.pdb"
    f=open(fout,"w")
    fin=open(mol,"r")
    f.write(biomt)
    for l in fin.readlines():
        if l.startswith("ATOM"):
            f.write(l)
    f.close()
    fin.close()
    
def buildRecipe(recipe,name):
    if recipe is None : 
        return None    
    for ingr in recipe.ingredients: 
        mols=buildIngredient(ingr,)
        c=0
        g=[]
        biomt="REMARK 350 BIOMOLECULE: 1  \n"
        for pos,rot in ingr.results:#[pos,rot]
#            print c,ingr.o_name+"_"+str(c),master_node
            #test if in bb ?
            #can we do it using a particle cloud
            mat = rot.copy().transpose()
            mat[:3, 3] = pos
            biomt+="REMARK 350   BIOMT1   %d  %.8f %.8f %.8f %.8f\n" % (c,rot[0][0],rot[0][1],rot[0][2],pos[0])
            biomt+="REMARK 350   BIOMT2   %d  %.8f %.8f %.8f %.8f\n" % (c,rot[1][0],rot[1][1],rot[1][2],pos[1])
            biomt+="REMARK 350   BIOMT3   %d  %.8f %.8f %.8f %.8f\n" % (c,rot[2][0],rot[2][1],rot[2][2],pos[2])
            c+=1       
        #addBIOMT
        writeMol(ingr.name,mols,biomt)
    return None

def build_scene(env,mb=None):
    r =  env.exteriorRecipe
    if r : 
        buildRecipe(r,r.name)
    for o in env.compartments:
        rs = o.surfaceRecipe
        if rs : 
            buildRecipe(rs,str(o.name)+"_surface")
        ri = o.innerRecipe
        if ri : 
            buildRecipe(ri,str(o.name)+"_interior")
    
if len(sys.argv) > 1 :#and doit :
    #f=PATH+"HIV-1_0.1.6-7_mixed_pdb.json"?
    filename = sys.argv[1]
    resultfile=None
    if filename in autopack.RECIPES :
        n=filename
        v=sys.argv[2]
        filename = autopack.RECIPES[n][v]["setupfile"]
        resultfile= autopack.RECIPES[n][v]["resultfile"]
    setupfile = autopack.retrieveFile(filename,cache="recipes")
    print ("ok use ",setupfile,filename)
    fileName, fileExtension = os.path.splitext(setupfile)
    n=os.path.basename(fileName)
    h = Environment(name=n)     
    h.loadRecipe(setupfile)
    h.setupfile=filename
    resultfile = PATH+"HIV-1_0.1.6-8_mixed_pdb.json"
    if resultfile is not None :
        h.resultfile=resultfile
    fileName, fileExtension = os.path.splitext(setupfile)
    rfile = h.resultfile
    resultfilename = autopack.retrieveFile(rfile,cache="results")
    if resultfilename is None :
        print ("no result for "+n+" "+h.version+" "+rfile)
        sys.exit()
    print ("get the result file from ",resultfilename)
    result,orgaresult,freePoint=h.loadResult(resultfilename=resultfilename)
#                                             restore_grid=False,backward=True)#load text ?#this will restore the grid  
    ingredients = h.restore(result,orgaresult,freePoint)
    #export the complete recipe as collada. each ingredient -> meshnode. Each instance->node instance
    env=h
    build_scene(env)

#    cxml.freeScene();
#execfile("pathto/export_recipe_collada.py") 
#I usually run this on with pmv,anaconda or mayapy
                

#execfile("/Users/ludo/DEV/git_upy/examples/export_collada.py")
#import upy
#helper = upy.getHelperClass()()
#helper.read("/Users/ludo/DEV/autopack_git/autoPACK_database_1.0.0/geometries/HIV1_capside_3j3q_Rep_Med_0_2_1.dae")
#