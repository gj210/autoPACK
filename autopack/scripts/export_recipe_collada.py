# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 22:08:01 2014

@author: ludo
"""
MAYA=False
SWAPZ=True
import sys
import os
import math

#windows path
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17_8DE13DAD\plugins\ePMV\mgl64\MGLToolsPckgs")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17_8DE13DAD\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")

if MAYA:
    sys.path.append("/Users/ludo/Library/Preferences/Autodesk/maya/2015-x64/plug-ins/MGLToolsPckgs")
    sys.path.append("/Users/ludo/Library/Preferences/Autodesk/maya/2015-x64/plug-ins/MGLToolsPckgs/PIL")
    #maya standalone special
    import maya.standalone
    maya.standalone.initialize()
    #load plugin
    import maya
    maya.cmds.loadPlugin("fbxmaya")

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
from collada import Collada
from collada import material
from collada import source
from collada import geometry
from collada import scene


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
          

def oneMaterial(name,collada_xml,color=None):
    if color == None :
        color = numpy.random.random(4)
        color[3] = 1.0
    if len(color) == 3 :
        color = [color[0],color[1],color[2],1.0]
    effect = material.Effect("effect"+name, [], "phong", 
                                 diffuse=color)
    mat = material.Material("material"+name, name+"_material", effect)
    matnode = scene.MaterialNode("material"+name, mat, inputs=[])
    collada_xml.effects.append(effect)
    collada_xml.materials.append(mat)
    return matnode

def colladaMesh(name,v,n,f,collada_xml,matnode=None):
    vertzyx = numpy.array(v)# * numpy.array([1,1,-1])
    z,y,x=vertzyx.transpose()
    if SWAPZ :
        vertxyz = numpy.vstack([x,y,z]).transpose()* numpy.array([1,1,-1])
    else :
        vertxyz = numpy.vstack([x,y,z]).transpose()        
    # vertxyz = numpy.vstack([x,y,z]).transpose()* numpy.array([1,1,-1])        
    vert_src = source.FloatSource(name+"_verts-array", vertxyz.flatten(), ('X', 'Y', 'Z'))
    input_list = source.InputList()
    input_list.addInput(0, 'VERTEX', "#"+name+"_verts-array")
    if type(n) != type(None) and len(n)  :
        norzyx=numpy.array(n)
        nz,ny,nx=norzyx.transpose()
        norxyz = numpy.vstack([nx,ny,nz]).transpose()#* numpy.array([1,1,-1])
#        if SWAPZ :
#            norxyz = numpy.vstack([nx,ny,nz]).transpose()* numpy.array([1,1,-1])
#        else :
#            norxyz = numpy.vstack([nx,ny,nz]).transpose()
        normal_src = source.FloatSource(name+"_normals-array", norxyz.flatten(), ('X', 'Y', 'Z'))
        geom = geometry.Geometry(scene, "geometry"+name, name, [vert_src,normal_src])
        input_list.addInput(0, 'NORMAL', "#"+name+"_normals-array")
    else :
        geom = geometry.Geometry(scene, "geometry"+name, name, [vert_src])        
    #invert all the face
    fi=numpy.array(f,int)[:,::-1]
#    if SWAPZ :
#        fi=numpy.array(f,int)[:,::-1]
#    else :
#        fi=numpy.array(f,int)        
    triset = geom.createTriangleSet(fi.flatten(), input_list, name+"materialref")
    geom.primitives.append(triset)
    collada_xml.geometries.append(geom)
    master_geomnode = scene.GeometryNode(geom, [matnode])
    master_node = scene.Node("node_"+name, children=[master_geomnode,])#,transforms=[tr,rz,ry,rx,s])
    return master_node        

def buildIngredientGeom(ingr,collada_xml,matnode):
    iname = ingr.o_name
    vertzyx = numpy.array(ingr.vertices)# * numpy.array([1,1,-1])
    z,y,x=vertzyx.transpose()
    if SWAPZ :
        vertxyz = numpy.vstack([x,y,z]).transpose()* numpy.array([1,1,-1])
    else :
        vertxyz = numpy.vstack([x,y,z]).transpose()    
    vert_src = source.FloatSource(iname+"_verts-array", vertxyz.flatten(), ('X', 'Y', 'Z'))
    if ingr.vnormals :
        norzyx=numpy.array(ingr.vnormals)
        nz,ny,nx=norzyx.transpose()
        if SWAPZ :
            norxyz = numpy.vstack([nx,ny,nz]).transpose()* numpy.array([1,1,-1])
        else :
            norxyz = numpy.vstack([nx,ny,nz]).transpose()
        normal_src = source.FloatSource(iname+"_normals-array", norxyz.flatten(), ('X', 'Y', 'Z'))
    geom = geometry.Geometry(collada_xml, "geometry"+iname, iname, [vert_src])
    input_list = source.InputList()
    input_list.addInput(0, 'VERTEX', "#"+iname+"_verts-array")
#    input_list.addInput(0, 'NORMAL', "#"+iname+"_normals-array")
    #invert all the face 
    fi=numpy.array(ingr.faces,int)[:,::-1]
#    if SWAPZ :
#        fi=numpy.array(ingr.faces,int)[:,::-1]
#    else :
#        fi=numpy.array(ingr.faces,int) 
    triset = geom.createTriangleSet(fi.flatten(), input_list, iname+"materialref")
    geom.primitives.append(triset)
    collada_xml.geometries.append(geom)
    master_geomnode = scene.GeometryNode(geom, [matnode])
    master_node = scene.Node("node_"+iname, children=[master_geomnode,])#,transforms=[tr,rz,ry,rx,s])
    return master_node    
    

def buildRecipe(recipe,name,collada_xml,root_node,mask=None):
    if recipe is None : return collada_xml,root_node 
    n=scene.Node(str(name))
    for ingr in recipe.ingredients: 
        #for each ingredient
        #build the material
        matnode=oneMaterial(ingr.name,collada_xml)
        #build a geomedtry node
        master_node=buildIngredientGeom(ingr,collada_xml,matnode)
        collada_xml.nodes.append(master_node)
        #build the scene instance node
        c=0
        g=[]
        for pos,rot in ingr.results:#[pos,rot]
            if mask != None :
                 if testOnePoints(pos,mask) :
                     continue
            geomnode = scene.NodeNode(master_node)
            mat = rot.copy()
            mat[:3, 3] = pos
            if helper.host == 'dejavu' :#and SWAPZ:#need to find the way that will work everywhere
               mry90 = helper.rotation_matrix(-math.pi/2.0, [0.0,1.0,0.0])#?
               mat = numpy.array(numpy.matrix(mat)*numpy.matrix(mry90)) 
                # mat = numpy.array(mat).transpose()
            scale, shear, euler, translate, perspective=decompose_matrix(mat)
            p=pos#matrix[3,:3]/100.0#unit problem
            tr=scene.TranslateTransform(p[0],p[1],p[2])
            rx=scene.RotateTransform(1,0,0,numpy.degrees(euler[0]))
            ry=scene.RotateTransform(0,1,0,numpy.degrees(euler[1]))
            rz=scene.RotateTransform(0,0,1,numpy.degrees(euler[2]))
            s=scene.ScaleTransform(scale[0],scale[1],scale[2])
            #n = scene.NodeNode(master_node,transforms=[tr,rz,ry,rx,s])
    #            gnode = scene.Node(self.getName(c)+"_inst", children=[geomnode,])
            ne = scene.Node(ingr.o_name+"_"+str(c), children=[geomnode,],transforms=[tr,rz,ry,rx,s]) #scene.MatrixTransform(matrix)                        
            g.append(ne)
            c+=1
        node = scene.Node(ingr.o_name, children=g)
        n.children.append(node)
    root_node.children.append(n)  
    return collada_xml,root_node

def buildCompartmentsGeomRep(comp,collada_xml,root_node,mask=None):
    if comp.representation_file is None : 
        return collada_xml,root_node
    nr=scene.Node(str(comp.name)+str("rep"))
    filename = autopack.retrieveFile(comp.representation_file,cache="geometries") #geometries    
    gdic=helper.read(filename)
    for nid in gdic :
        matnode=oneMaterial(str(nid),collada_xml,color=gdic[nid]["color"])
        master_node = colladaMesh(str(nid),gdic[nid]["mesh"][0],gdic[nid]["mesh"][1],gdic[nid]["mesh"][2],collada_xml,matnode=matnode)
#        master_node = colladaMesh(str(nid),comp.vertices,comp.vnormals,comp.faces,collada_xml,matnode=matnode)
        collada_xml.nodes.append(master_node)
#        mxmesh.setParent(nr)
        if len(gdic[nid]['instances']):
            geomnode = scene.NodeNode(master_node)
#            !n=scene.Node(str(nid)+str("instances"))
#            nri.setParent(nr)
            g=[]
            c=0
            for mat in gdic[nid]['instances']:
                geomnode = scene.NodeNode(master_node)
#                instance = scene.createInstancement(str(nid)+"_"+str(c),mxmesh)
                mat = numpy.array(mat,float)#.transpose()
                if helper.host == 'dejavu':#need to find the way that will work everywhere
                   mry90 = helper.rotation_matrix(-math.pi/2.0, [0.0,1.0,0.0])#?
                   mat = numpy.array(numpy.matrix(mat)*numpy.matrix(mry90)) 
                    # mat = numpy.array(mat).transpose()
                scale, shear, euler, translate, perspective=decompose_matrix(mat)
                p=translate#matrix[3,:3]/100.0#unit problem
                tr=scene.TranslateTransform(p[0],p[1],p[2])
                rx=scene.RotateTransform(1,0,0,numpy.degrees(euler[0]))
                ry=scene.RotateTransform(0,1,0,numpy.degrees(euler[1]))
                rz=scene.RotateTransform(0,0,1,numpy.degrees(euler[2]))
                s=scene.ScaleTransform(scale[0],scale[1],scale[2])
                ne = scene.Node(str(nid)+"_"+str(c), children=[geomnode,],transforms=[tr,rz,ry,rx,s]) #scene.MatrixTransform(matrix)                        
                g.append(ne)
                c+=1
            node = scene.Node(str(nid)+str("instances"), children=g)
            nr.children.append(node)
    root_node.children.append(nr)  
    return collada_xml,root_node

def buildCompartmentsGeom(comp,collada_xml,root_node,mask=None):
    nr=scene.Node(str(comp.name))
    #filename = autopack.retrieveFile(comp.representation_file,cache="geometries") #geometries    
    #gdic=helper.read(filename)
    matnode=oneMaterial(str(comp.name),collada_xml)
    #v,n,f
    master_node = colladaMesh(str(comp.name),comp.vertices*1.0/collada_xml.assetInfo.unitmeter,None,comp.faces,collada_xml,matnode=matnode)
#    collada_xml.nodes.append(master_node)
#    geomnode = scene.NodeNode(master_node)
#    p=[0,0,0]#matrix[3,:3]/100.0#unit problem
#    euler = [0,0,0]        
#    scale = [1,1,1]
#    tr=scene.TranslateTransform(p[0],p[1],p[2])
#    rx=scene.RotateTransform(1,0,0,numpy.degrees(euler[0]))
#    ry=scene.RotateTransform(0,1,0,numpy.degrees(euler[1]))
#    rz=scene.RotateTransform(0,0,1,numpy.degrees(euler[2]))
#    s=scene.ScaleTransform(scale[0],scale[1],scale[2])
#    ne = scene.Node("compartments_"+str(comp.name), children=[geomnode,],transforms=[tr,rz,ry,rx,s]) #scene.MatrixTransform(matrix)                        
    root_node.children.append(master_node)  
    return collada_xml,root_node
                
def build_scene(env,outfilename):
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
    collada_xml = Collada()
    collada_xml.assetInfo.unitname="centimeter"
    collada_xml.assetInfo.unitmeter=0.01
    collada_xml.assetInfo.upaxis="Y_UP"

    root_env=scene.Node(env.name)
    myscene = scene.Scene(env.name+"_Scene", [root_env])
    collada_xml.scenes.append(myscene)
    collada_xml.scene = myscene                
#
    bbsurface = None
    r =  env.exteriorRecipe
    if r : collada_xml,root_env = buildRecipe(r,r.name,collada_xml,root_env)
    for o in env.compartments:
        rs = o.surfaceRecipe
        if rs : 
            #p,s,bb=up (767.0) #used for lipids
            #pp,ss,bbsurface = up (700.0)
            collada_xml,root_env = buildRecipe(rs,str(o.name)+"_surface",collada_xml,root_env,mask=bbsurface)
        ri = o.innerRecipe
        if ri : collada_xml,root_env = buildRecipe(ri,str(o.name)+"_interior",collada_xml,root_env,mask=bbsurface)
    collada_xml,root_env = buildCompartmentsGeom(env.compartments[0],collada_xml,root_env)
    collada_xml.write(outfilename)
    return collada_xml
    
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

resultfile=None    
if len(sys.argv) > 1 :
    filename = sys.argv[1]
    if filename in autopack.RECIPES :
        n=filename
        v=sys.argv[2]
        filename = autopack.RECIPES[n][v]["setupfile"]
        resultfile= autopack.RECIPES[n][v]["resultfile"]
    setupfile = autopack.retrieveFile(filename,cache="recipes")
    print ("ok use ",setupfile,filename)
    fileName, fileExtension = os.path.splitext(setupfile)
    n=os.path.basename(fileName)
else :
    setupfile = filename = "C:\Users\ludov\Google Drive\TSRI HIV Modeling\HIVsphere_simple_r140.json"
#    resultfile = "C:\Users\ludov\Google Drive\TSRI HIV Modeling\cellpack_models\Sphere_140_m330_random_tr.json"
#    resultfile = "C:\Users\ludov\Google Drive\TSRI HIV Modeling\cellpack_models\Sphere_140_10_random_tr.json"
    resultfile = "C:\Users\ludov\Google Drive\TSRI HIV Modeling\cellpack_models\HIV_VLP_140_214_random_tr.json"
#    resultfile = "C:\Users\ludov\OneDrive\Documents\\backupC4DHiv\sphere140_rand_49_tr.json"
#    setupfile = filename = "C:\Users\ludov\Google Drive\TSRI HIV Modeling\HIVsphere_simple_r139.json"
#    resultfile = "C:\Users\ludov\OneDrive\Documents\\backupC4DHiv\sphere139_random_127_tr.json"   
    n="HIV_VLP"
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
h.compartments[0].name = "HIVenvelop_140"
result,orgaresult,freePoint=h.loadResult(resultfilename=resultfilename)
#                                             restore_grid=False,backward=True)#load text ?#this will restore the grid  
ingredients = h.restore(result,orgaresult,freePoint)
#export the complete recipe as collada. each ingredient -> meshnode. Each instance->node instance
env=h
cxml=build_scene(env,"C:\Users\ludov\Google Drive\TSRI HIV Modeling\cellpack_models\HIV_VLP_140_214_random.dae")
#execfile("pathto/export_recipe_collada.py") 
#I usually run this on with pmv,anaconda or mayapy
                    

#execfile("/Users/ludo/DEV/git_upy/examples/export_collada.py")
#import upy
#helper = upy.getHelperClass()()
#helper.read("/Users/ludo/DEV/autopack_git/autoPACK_database_1.0.0/geometries/HIV1_capside_3j3q_Rep_Med_0_2_1.dae")
#