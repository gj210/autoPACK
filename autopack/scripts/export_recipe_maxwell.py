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
#sys.path.append("/usr/local/maxwell-3.0//python/pymaxwell/python2.7_UCS4")
sys.path.append("/usr/local/maxwell-3.0//python/pymaxwell/python2.7_UCS2")
sys.path.append("/usr/local/maxwell-3.0/")
from pymaxwell import *


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
         
def oneMaterial(name,scene,color=None):
	# Add a new material to the scene
    if color == None :
         color = numpy.random.random(4)
         color[3] = 1.0
    if len(color) == 3 :
        color = [color[0],color[1],color[2],1.0]
    mat = scene.createMaterial(str(name));
    mat.read("/usr/local/maxwell-3.0/materials database/mxm files/sss_raki.mxm")
    #return mat
    layer=mat.getLayer(0)
    bsdf = layer.getBSDF(0)
    refl = bsdf.getReflectance();
    attr = Cattribute();
    # Coating Reflectance 0
    attr.activeType = MAP_TYPE_RGB;
    attr.rgb.assign(float(color[0]),float(color[1]),float(color[2]));
    bsdf.setAttribute('color',attr);
    refl.setAttribute('color',attr);#0
    refl.setAttribute('scattering',attr);
    attr.rgb.assign(1,1,1);
    bsdf.setAttribute('color.tangential',attr);
    refl.setAttribute('color.tangential',attr);#90
    # BSDF Transmittance
    attr.rgb.assign(1,1,1);
    refl.setAttribute('transmittance.color',attr);
    return mat
 
def oneMaterial2(name,scene,color=None):
    # Add a new material to the scene
    if color == None :
        color = numpy.random.random(4)
        color[3] = 1.0
    if len(color) == 3 :
        color = [color[0],color[1],color[2],1.0]
#    print "what ?"
    mat = scene.createMaterial(name);
#    print "mat ",mat
    mat.read("/usr/local/maxwell-3.0/materials database/mxm files/wall_beige.mxm")
#    print "mat ",mat
    #return mat
    layer=mat.getLayer(0)
    bsdf = layer.getBSDF(0)
    refl = bsdf.getReflectance();
    attr = Cattribute();
    # Coating Reflectance 0
    attr.activeType = MAP_TYPE_RGB;
    attr.rgb.assign(color[0],color[1],color[2]);
    bsdf.setAttribute('color',attr);
    refl.setAttribute('color',attr);#0
    refl.setAttribute('scattering',attr);
    attr.rgb.assign(1,1,1);
    bsdf.setAttribute('color.tangential',attr);
    refl.setAttribute('color.tangential',attr);#90
    # BSDF Transmittance
    attr.rgb.assign(1,1,1);
    return mat 
 
# def oneMaterial(name,collada_xml,color=None):
#    if color == None :
#        color = numpy.random.random(4)
#        color[3] = 1.0
#    if len(color) == 3 :
#        color = [color[0],color[1],color[2],1.0]
#    effect = material.Effect("effect"+name, [], "phong", 
#                                 diffuse=color)
#    mat = material.Material("material"+name, name+"_material", effect)
#    matnode = scene.MaterialNode("material"+name, mat, inputs=[])
#    collada_xml.effects.append(effect)
#    collada_xml.materials.append(mat)
#    return matnode

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

def maxwellMesh(name,v,n,f,scene,matnode=None):
    vertzyx = numpy.array(v,float)# * numpy.array([1,1,-1])
    fi=numpy.array(f,int)
    nVertices = len(vertzyx)
    nNormals  = 0#len(nVertices)
    nTriangles= len(fi)
    nPositionsPerVertex = 1;
    mxobject = scene.createMesh(str("mesh_")+str(name),nVertices,nNormals,nTriangles,nPositionsPerVertex);
#    ok = mxobject.initializeMesh(nVertices,nNormals,nTriangles,nPositionsPerVertex)
    if mxobject.isNull():
        print("Error creating mesh");
        return 0;
    v=vertzyx[0]
    ok = mxobject.setVertex(0,0,Cvector(v[0],v[1],v[2]))
    for i,v in enumerate(vertzyx[1:]):
        ok |= mxobject.setVertex(i+1,0,Cvector(v[0],v[1],v[2])) 
    if not ok:
        print("Error setting vertices");
        return 0;        
    f = fi[0]
    ok = mxobject.setTriangle(0,f[0],f[1],f[2],f[0],f[1],f[2])
    for i,f in enumerate(fi[1:]):
        ok |= mxobject.setTriangle(i+1,f[0],f[1],f[2],f[0],f[1],f[2])
    if not ok:
        print("Error setting faces");
        return 0;            
    	# Add the sphere to the new scene
#	obj2 = new_scene.addObject(obj);	
	# Assign the material to the object of the new scene
    if matnode is not None :
        mxobject.setMaterial(matnode);
    return mxobject        

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
                
def buildIngredientGeom(ingr,scene,matnode):
#    sc = Cmaxwell(mwcallback);
    iname = ingr.o_name
#    if ingr.Type == 'SingleSphere' :
#        return
#    if ingr.Type == 'Grow' :
#        return    
    if not len(ingr.vertices) :
        ingr.getData()
    vertzyx = numpy.array(ingr.vertices,float)# * numpy.array([1,1,-1])
    fi=numpy.array(ingr.faces,int)
    nVertices = len(vertzyx)
#    norzyx=numpy.array(ingr.vnormals,float)
    nNormals  = 0#len(nVertices)
    nTriangles= len(fi)
    nPositionsPerVertex = 1;
    mxobject = scene.createMesh(str("mesh_")+str(iname),nVertices,nNormals,nTriangles,nPositionsPerVertex);
#    ok = mxobject.initializeMesh(nVertices,nNormals,nTriangles,nPositionsPerVertex)
    if mxobject.isNull():
        print("Error creating mesh");
        return 0;
    v=vertzyx[0]
    ok = mxobject.setVertex(0,0,Cvector(v[0],v[1],v[2]))
    for i,v in enumerate(vertzyx[1:]):
        ok |= mxobject.setVertex(i+1,0,Cvector(v[0],v[1],v[2])) 
    if not ok:
        print("Error setting vertices");
        return 0;        
#    if ingr.vnormals :
#    v=norzyx[0]
#    ok = mxobject.setNormal(0,0,Cvector(v[0],v[1],v[2]))
#    for i in range(1,nNormals):
#        v=nNormals[i]
#        ok |= mxobject.setNormal(i+1,0,Cvector(v[0],v[1],v[2])) 
#    if not ok:
#        print("Error setting normals");
#        return 0;     
    f = fi[0]
    ok = mxobject.setTriangle(0,f[0],f[1],f[2],f[0],f[1],f[2])
    for i,f in enumerate(fi[1:]):
        ok |= mxobject.setTriangle(i+1,f[0],f[1],f[2],f[0],f[1],f[2])
    if not ok:
        print("Error setting faces");
        return 0;            
    	# Add the sphere to the new scene
#	obj2 = new_scene.addObject(obj);	
	# Assign the material to the object of the new scene
    if matnode is not None :
        mxobject.setMaterial(matnode);
#    print mxobject,"created",mxobject.isNull()
#    output=""
#    ok = sc.writeMXS("/home/ludo/autopack_test1.mxs");
#    if not ok:
#		print("Error saving MXS ("+output+")");
#		return 0;

    return mxobject    
#    

def buildRecipe(recipe,name,scene,root_node,mask=None):
    if recipe is None : 
        return scene,root_node
    nr=scene.createNullObject(str(name))
    nr.setParent(root_node)
#    n=scene.Node(str(name))

    for ingr in recipe.ingredients: 
        #for each ingredient
        #build the material
        matnode=oneMaterial(ingr.o_name,scene,color=ingr.color)
#        print "material",matnode
        #build a geomedtry node
        master_node=buildIngredientGeom(ingr,scene,matnode)
        if master_node is None :
            continue
        #hide it
        master_node.setHideToCamera(True)
        master_node.setParent(nr)
#        master_node,scene=readOne(ingr.meshFile,ingr.meshName,scene)
#        print "master_node",type(master_node),master_node
#        master_node = scene.addObject(master_node);
#        print "master_node",master_node
        #apply scatter/cloner?
        n=scene.createNullObject(str(ingr.o_name))
        n.setParent(nr)
#        collada_xml.nodes.append(master_node)
        #build the scene instance node
        c=0
        g=[]
        for pos,rot in ingr.results:#[pos,rot]
#            print c,ingr.o_name+"_"+str(c),master_node
            #test if in bb ?
            if mask != None :
                 if testOnePoints(pos,mask) :
                     continue
            instance = scene.createInstancement(str(ingr.o_name)+"_"+str(c),master_node);
            mat = rot.copy()
            mat[:3, 3] = pos
#            if helper.host == 'dejavu':#need to find the way that will work everywhere
#               mry90 = helper.rotation_matrix(-math.pi/2.0, [0.0,1.0,0.0])#?
#               mat = numpy.array(numpy.matrix(mat)*numpy.matrix(mry90)) 
                # mat = numpy.array(mat).transpose()
            scale, shear, euler, translate, perspective=decompose_matrix(mat)
            #instance.setParent(n)
            #instance.setPosition(Cvector(pos[0],pos[1],pos[2]))

            base,pivot,ok = instance.getBaseAndPivot();
            base.origin = Cvector(pos[0],pos[1],pos[2])
                        #pivot =  Cvector(euler[0],euler[1],euler[2])  
            mat = numpy.array(mat).transpose()
#            print ("tr,eiler",translate,euler,mat)
            base.xAxis=Cvector(mat[0][0],mat[0][1],mat[0][2])
            base.yAxis=Cvector(mat[1][0],mat[1][1],mat[1][2])
            base.zAxis=Cvector(mat[2][0],mat[2][1],mat[2][2])
            
            instance.setBaseAndPivot(base,pivot);  
            #instance.setRotation(Cvector(euler[0],euler[1],euler[2]))
            # Scale instance
            #base.xAxis.scale(RADIUS[name[i]]*2*SCALE);
            #base.yAxis.scale(RADIUS[name[i]]*2*SCALE);
            #base.zAxis.scale(RADIUS[name[i]]*2*SCALE);
		
			# Apply translation and scale
    		
#            instance.setBaseAndPivot(base,pivot); 
            instance.setParent(n)
            c+=1
#            print instance
#            break
    return scene,root_node
        
def build_scene(env):
    #architecture is :
    #-env.name
    #--exterior
    #----ingredients
    #--compartments
    #----surface
    #------ingredients
    #----interior
    #------ingredients
    #create the document and a node for rootenv
    scene = Cmaxwell(mwcallback);
    # Create camera
    camera = scene.addCamera('camera',1,1.0/800.0,0.04,0.035,400,"CIRCULAR",0,0,25,10000,10000,1);
    camera.setStep(0,Cvector(0.0,0.0,6.0),Cvector(0.0,0.0,0.0),Cvector(0.0,1.0,0.0),0.035,8.0,0);
    camera.setActive();
    # Add physical sky and sun
    scene.getEnvironment().setActiveSky('physical');
    sunColor = Crgb()
    sunColor.assign(1,1,1)
    scene.getEnvironment().setSunProperties(SUN_PHYSICAL,5777,1,1,sunColor);
    
    root_env=scene.createNullObject(str(env.name))
    r =  env.exteriorRecipe
    if r : scene,root_env = buildRecipe(r,r.name,scene,root_env)
    for o in env.compartments:
        rs = o.surfaceRecipe
        if rs : 
            p,s,bb=up (767.0) #used for lipids
            pp,ss,bbsurface = up (700.0)
            bbsurface = numpy.array([[p-ss/2.0],[p+ss/2.0]])
            scene,root_env = buildRecipe(rs,str(o.name)+"_surface",scene,root_env,mask=bbsurface)
        ri = o.innerRecipe
        if ri : 
            pp,ss,bbmatrix = up (650.0)
            bbmatrix = numpy.array([[p-ss/2.0],[p+ss/2.0]])
            scene,root_env = buildRecipe(ri,str(o.name)+"_interior",scene,root_env,mask=bbsurface)
        #build the compartments geometry
#        buildCompartmentsGeom(o,scene,parent=root_env)
        #['ID']['node', 'color', 'id', 'instances', 'mesh', 'parentmesh']
        
    fname = "/home/ludo/"+env.name+str(env.version)+".mxs"
    ok = scene.writeMXS(str(fname));
    print "write",ok
    if not ok:
		print("Error saving MXS ('/home/ludo/"+env.name+str(env.version)+".mxs')");
		return 0;
    return scene
    
#        ri = o.innerRecipe
#        if ri :
#            n=scene.Node(str(o.name)+"_interior")
#            root_env.children.append(n)
#            
#    node = scene.Node("HIV1_capside_3j3q_Rep_Med")
#    parent_object=helper.getObject("Pentamers")
#    mesh=None
#    if MAYA:
#        mesh=helper.getObject("HIV1_capsid_3j3q_Rep_Med_Pent_0_1_0_1")
#    collada_xml=helper.instancesToCollada(parent_object,collada_xml=None,instance_node=True,parent_node=node,mesh=mesh)
#    parent_object=helper.getObject("Hexamers")
#    mesh=None
#    if MAYA:
#        mesh=helper.getObject("HIV1_capsid_3j3q_Rep_Med_0_1_0")
#    collada_xml=helper.instancesToCollada(parent_object,collada_xml=collada_xml,instance_node=True,parent_node=node,mesh=mesh)
#    #collada_xml.scene.nodes
#    collada_xml.write("/Users/ludo/DEV/autopack_git/autoPACK_database_1.0.0/geometries/HIV1_capside_3j3q_Rep_Med_0_2_1.dae") 
#    


mgr = CextensionManager.instance();
#mgr.setExtensionsDirectory('C:/Program Files/Next Limit/Maxwell 3/extensions/')
mgr.loadAllExtensions()
print mgr
#ext = mgr.createDefaultGeometryProceduralExtension('MaxwellParticles')
#print ext
#params = ext.getExtensionData()
  
doit=True
if len(sys.argv) > 1 :#and doit :
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
    if doit :
        cxml=build_scene(env)
#    cxml.freeScene();
#execfile("pathto/export_recipe_collada.py") 
#I usually run this on with pmv,anaconda or mayapy
                

#execfile("/Users/ludo/DEV/git_upy/examples/export_collada.py")
#import upy
#helper = upy.getHelperClass()()
#helper.read("/Users/ludo/DEV/autopack_git/autoPACK_database_1.0.0/geometries/HIV1_capside_3j3q_Rep_Med_0_2_1.dae")
#