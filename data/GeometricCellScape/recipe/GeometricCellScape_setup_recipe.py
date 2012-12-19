# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 14:25:36 2011

###############################################################################
#
# autoPACK Authors: Graham T. Johnson, Ludovic Autin, Mostafa Al-Alusi, Michel Sanner
#   Based on COFFEE Script developed by Graham Johnson between 2005 and 2010
#   with assistance from Mostafa Al-Alusi in 2009 and periodic input
#   from Arthur Olson's Molecular Graphics Lab
#
# AFGui.py Authors: Ludovic Autin with minor editing/enhancement from Graham Johnson
#
# Copyright: Graham Johnson ©2010
#
# This file is part of autoPACK, cellPACK, and AutoFill.
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
###############################################################################

Name: 'GeometricCellScape_setup_recipe'
@author: Graham Johnson and Ludovic Autin
"""

import sys
#AUTOFILL
import AutoFill
wrkDir=localwrkDir = AutoFill.__path__[0]
#httpwrkDir = "http://grahamj.com/autofill/autoFillData/GeometricCellScape/"
httpwrkDir = "http://autofill.googlecode.com/svn/data/GeometricCellScape/"
#sphdir = "http://grahamj.com/autofill/autoFillData/GeometricCellScape/spheres/"
wrkDirMesh = httpwrkDir+"geoms/"#"http://grahamj.com/autofill/autoFillData/GeometricCellScape/Geometries/"

localdir = AutoFill.__path__[0]+os.sep+"autoFillRecipeScripts"+os.sep+"GeometricCellScape"
#localdir = AutoFill.RECIPES["GeometricCellScape"]["wrkdir"]

from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
from AutoFill.Ingredient import MultiCylindersIngr, GrowIngrediant,ActinIngrediant
from AutoFill.Organelle import Organelle
from AutoFill.Recipe import Recipe
from AutoFill.HistoVol import Environment
from AutoFill.autofill_viewer import AFViewer


helper = AutoFill.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper =helperClass()

#DEJAVU COLORS
from upy.colors import red, aliceblue, antiquewhite, aqua, goldenrod, thistle, violet, skyblue, royalblue, plum, pink, orchid, mediumvioletred, lightpink, \
     aquamarine, azure, beige, bisque, black, blanchedalmond, mediumpurple, mediumorchid, hotpink, lightskyblue, lightblue, lightcyan, dodgerblue, firebrick, \
     blue, blueviolet, brown, burlywood, cadetblue, slateblue, darkslateblue, darkviolet, darkmagenta, powderblue, palevioletred, indigo,\
     chartreuse, chocolate, coral, cornflowerblue, cornsilk, \
     crimson, cyan, darkblue, darkcyan, darkgoldenrod, \
     orange, purple, deeppink, lightcoral, \
     blue, cyan, mediumslateblue, steelblue, darkcyan, lime,\
     limegreen, darkorchid, tomato, khaki, gold, magenta, green, snow, grey


MSca = 1.0
rSurf ={"Ellipse":Recipe(),"Sphere":Recipe()}
rMatrix = {"Ellipse":Recipe(),"Sphere":Recipe()}
rCyto = Recipe()

VesK = 0.56
mode = "jitter"#"rigid-body" #or "rigid-body" except surface that are always jitter.
doCyot = True
dosurface1 = True
dosurface2 = True
doinner1 = True
doinner2 = True
dolipid = True
doPrev = False

inPriority=0
outPriority=0
surfPriority=-10
inNMol=0
outNMol=0
surfNMol=0
dict = {}
listename = ["BigCuboid","TallCuboid","SmallCuboid",
             "LipidOut","LipidIn",
             "Transporter","Docker","Brancher",
             "RegularElipsoid","TallElipsoid","SmallSphere","BigSphere",
             "Pyramid_Small","Pyramid_Tall","Pyramid_Wide","Pyramid_Big"]
             
dict["BigCuboid"] = {"packingPriority":15.,"nMol":15,"color":red}
dict["TallCuboid"] = {"packingPriority":20.,"nMol":30,"color":red}
dict["SmallCuboid"] = {"packingPriority":0.001,"nMol":1000,"color":red}

dict["LipidOut1"] = {"packingPriority":200.,"nMol":0,"color":grey}
dict["LipidIn1"] = {"packingPriority":200.,"nMol":0,"color":grey}

dict["Transporter1"] = {"packingPriority":-8,"nMol":10,"color":grey}
dict["Docker1"] = {"packingPriority":-9,"nMol":10,"color":grey}
dict["Brancher1"] = {"packingPriority":-10,"nMol":10,"color":grey}

dict["LipidOut2"] = {"packingPriority":200.,"nMol":0,"color":grey}
dict["LipidIn2"] = {"packingPriority":200.,"nMol":0,"color":grey}

dict["Transporter2"] = {"packingPriority":-8,"nMol":10,"color":grey}
dict["Docker2"] = {"packingPriority":-9,"nMol":10,"color":grey}
dict["Brancher2"] = {"packingPriority":-10,"nMol":10,"color":grey}

dict["RegularElipsoid"] = {"packingPriority":5,"nMol":20,"color":blue}
dict["TallElipsoid"] = {"packingPriority":20,"nMol":30,"color":blue}
dict["SmallSphere"] = {"packingPriority":0.001,"nMol":0,"color":blue}
dict["BigSphere"] = {"packingPriority":15,"nMol":15,"color":blue}

dict["Pyramid_Small"] = {"packingPriority":0.001,"nMol":0,"color":green}
dict["Pyramid_Tall"] = {"packingPriority":20,"nMol":30,"color":green}
dict["Pyramid_Wide"] = {"packingPriority":5,"nMol":20,"color":green}
dict["Pyramid_Big"] = {"packingPriority":15,"nMol":15,"color":green}


axes = (0.,0.,1.0)
principalVector=(1,0,0)
#h1.boundingBox = [[0.,0.,0.],[300.,300.,300.]]
bbox=helper.getObject("histoVolBB")
if bbox is None : 
    bbox = helper.box("histoVolBB",center=[0.,0.,0.],size=[230.,768.,1024.])#cornerPoints=h1.boundingBox)
if helper.host == "dejavu":
    #create the box
#    bbox = helper.box("HBB",center=[0.,0.,0.],size=[230.,768.,1024.])[0]
    axes = (1.,0.,0.0)	
    principalVector=(0,0,1)	

#this version use meshFiles instead of meshObject
#===============================================================================
# outside exteriorrecipe
#===============================================================================
if doCyot :
    #hactine=helper.getObject("ActineUnit")#None#getActine unit either one turn or on unit.in that case need to use the spline warp.
    #if hactine is None :
       	
    Actineparam={}
    Actineparam["test1"]={"length":600.,"influenceRad":20.,"marge":20.,
                            "weight":0.2,"nMol":5,"color":None}#1000,1000
    dFunc=[lambda x : 1./x,lambda x : 1./(x*x),None]
    
    testName="test1"
    idFunc = 0
    actine = ActinIngrediant(MSca*0.0001,  pdb=None,
                                  name='Actine', 
                                  positions=[[[0,0,0]]], positions2=[[[0.,60.,0]]],#100
                                  principalVector=(0,1,0),
                                  packingPriority=-1,biased=0.5,
                                  modelType="Cylinders",placeType="jitter",
                                  nbJitter=1, jitterMax=(0,0,0),
                                  closed = False,
                                  #meshObject=hactine,
                                  meshFile=wrkDirMesh+'/Actine.dae',	
                                  orientation = (1,0,0),
                                  distExpression = dFunc[idFunc],
                                  **Actineparam[testName]
                                  )
    rCyto.addIngredient( actine ) #UNTIL CODE FIXED
#    actine.isAttractor = False
#    actine.isAttractor = True
#    actine.constraintMarge = True
#    actine.seedOnMinus = True
#    actine.influenceRad = influenceRad
#    actine.oneSuperTurn = 825.545#cm from c4d graham file
#    actine.oneDimerSize = 100.0#200 =2 
#    actine.cutoff_surface = 50.
#    actine.cutoff_boundary = 1.0    
    #c1=helper.getObject("BigCuboid")
    #3 pyramid ingredient - singleIngreSpher rigid-body?
    cube1 = SingleSphereIngr( MSca*.0005,  75.,#50., 
                                name='iBigCuboid', pdb=None,
                                #meshObject=c1,
                                meshFile=wrkDirMesh+'/iBigCuboid',
                                placeType=mode,
                                **dict["BigCuboid"]
                                #packingMode='close'
                                ) #original radius is 3.61
                                
    rCyto.addIngredient( cube1 )
    
    #c2=helper.getObject("TallCuboid")
    #3 pyramid ingredient - singleIngreSpher rigid-body?
    cube2 = MultiCylindersIngr(MSca*.002, pdb='LipOut26',
                                 name='iTallCuboid', radii=[[15.0]],
                                  positions=[[[0,-100.,0]]], positions2=[[[0,100,0]]],
                    				meshFile=wrkDirMesh+'/iTallCuboid',
                                  #meshObject=c2,
                                  jitterMax=(0.3,0.8,0.8), 
                                  principalVector=principalVector,placeType=mode,
                                  **dict["TallCuboid"]
                                  )
    rCyto.addIngredient( cube2 )
    #its too tall
    cube2.useLength = True
    cube2.length = 150.0
    
    c3=helper.getObject("SmallCuboid")
    #3 pyramid ingredient - singleIngreSpher rigid-body?
    cube3 = SingleSphereIngr( MSca*.01,  30.0,#15., 
                                name='iSmallCuboid', pdb=None,
                                meshFile=wrkDirMesh+'/iSmallCuboid',#meshObject=c3,
                                placeType=mode,
                                **dict["SmallCuboid"]
                                #packingMode='close'
                                ) #original radius is 3.61
    rCyto.addIngredient( cube3 )

#===============================================================================
# ellipse pyramid ingredient
#===============================================================================
if doinner1:
    ho1=helper.getObject("Pyramid_Small")
    #3 pyramid ingredient - singleIngreSpher rigid-body?
    pyramid1 = SingleSphereIngr( MSca*.01,  15., 
                                name='iPyramid_Small', pdb=None,
                                meshFile=wrkDirMesh+'/iPyramid_Small',#meshObject=ho1,
                                placeType=mode,
                                **dict["Pyramid_Small"]
                                #packingMode='close'
                                ) #original radius is 3.61
    rMatrix["Ellipse"].addIngredient( pyramid1 )
    
    ho2 = helper.getObject("Pyramid_Tall")
    pyramid2 = MultiCylindersIngr(MSca*.0005, pdb=None,
                                 name='iPyramid_Tall', radii=[[20.]],
                                  positions=[[[0,-100.,0]]], positions2=[[[0,100,0]]],
                                  meshFile=wrkDirMesh+'/iPyramid_Tall',#meshObject=ho2,
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(0,1,0),
                                  placeType=mode,
                                  **dict["Pyramid_Tall"]
                                  )
    rMatrix["Ellipse"].addIngredient( pyramid2 )
    
    ho3 = helper.getObject("Pyramid_Wide")
    pyramid3 = MultiCylindersIngr(MSca*.002,  pdb=None,
                                 name='iPyramid_Wide', radii=[[55.]],
                                  positions=[[[0.,-15.,0]]], positions2=[[[0.,15,0]]],
                                  meshFile=wrkDirMesh+'/iPyramid_Wide',#meshObject=ho3,
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(0,1,0),
                                  placeType=mode,
                                  **dict["Pyramid_Wide"]
                                  )
    rMatrix["Ellipse"].addIngredient( pyramid3 )

#===============================================================================
# spher ellipsoid ingredient
#===============================================================================
if doinner2 :
    sph1 = helper.getObject("RegularElipsoid")
    ingrsph1 = MultiSphereIngr( VesK*.002,  pdb=None,
                                 name='iRegularElipsoid',
                                 radii=[[20.,30.,20.],],
                                 positions=[[[0.,-30.,0.],[0.,0.,0.],[0.,30.,0.]],],
                                 meshFile=wrkDirMesh+'/iRegularElipsoid',#meshObject=sph1,
                                 placeType=mode,
                                 #nMol=10,
                                 **dict["RegularElipsoid"]
                                 #jitterMax=(1,1,0.2),
                                 #principalVector=(0,0,-1)
                                 )
    rMatrix["Sphere"].addIngredient(ingrsph1)
    
    sph2 = helper.getObject("TallElipsoid")
    ingrsph2 = MultiCylindersIngr(MSca*.0005,  pdb=None,
                                 name='iTallElipsoid', radii=[[20.]],
                                  positions=[[[0,-100.,0]]], positions2=[[[0,100,0]]],
                                  meshFile=wrkDirMesh+'/iTallElipsoid',#meshObject=sph2,
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(0,1,0),placeType=mode,
                                  **dict["TallElipsoid"]
                                  )
    rMatrix["Sphere"].addIngredient(ingrsph2)

    sph3 = helper.getObject("SmallSphere")
    ingrsph3 = SingleSphereIngr( MSca*.01,  15., 
                                name='iSmallSphere', pdb=None,
                                meshFile=wrkDirMesh+'/iSmallSphere',#meshObject=sph3,
                                placeType=mode,
                                **dict["SmallSphere"]
                                #packingMode='close'
                                ) #original radius is 3.61
    rMatrix["Sphere"].addIngredient(ingrsph3)

#    sph4 = helper.getObject("BigSphere")
 #   ingrsph4 = SingleSphereIngr( MSca*.01,  50., 
 #                               name='iBigSphere', pdb=None,
 #                               meshFile=wrkDirMesh+'/iBigSphere',#meshObject=sph4,
  #                              placeType=mode,
  #                              **dict["BigSphere"]
  #                              #packingMode='close'
  #                              ) #original radius is 3.61
  #  rMatrix["Sphere"].addIngredient(ingrsph3)

#===============================================================================
# Lipdi from vG10
#===============================================================================
#lipidmeshin = helper.getObject("LipidIn")
#lipidmeshout = helper.getObject("LipidOut")
if dolipid :
    cyl26IngrO = MultiCylindersIngr(MSca*10.,  pdb='LipOut26',
                                 name='LipOut26', radii=[[4.5]],
                                  positions=[[[0,-.5,0]]], positions2=[[[0,25.5,0]]],
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(0,1,0),
                                  placeType='jitter',
                                  meshFile=wrkDirMesh+'/LipOut26',#meshObject=lipidmeshout,
                                **dict["LipidOut1"]
                                  )
    rSurf["Ellipse"].addIngredient(cyl26IngrO)
    
    # Inner Leaflet Lipids as Cyninders
    cyl26IngrI = MultiCylindersIngr(MSca*10.,  pdb='LipIn26',
                                 name='LipIn26', radii=[[4.5]],
                                  positions=[[[0,-25.5,0]]], positions2=[[[0,+5.,0]]],
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(0,1.,0),placeType='jitter',
                                  meshFile=wrkDirMesh+'/LipIn26',#meshObject=lipidmeshin,
                                  **dict["LipidIn1"]
                                  )
    rSurf["Ellipse"].addIngredient(cyl26IngrI)
    
    # Outer Leaflet Lipids as Cylinders
    cyl26IngrO2 = MultiCylindersIngr(MSca*10.,  pdb='LipOut262',
                                 name='LipOut262', radii=[[4.5]],
                                  positions=[[[0,-.5,0]]], positions2=[[[0,25.5,0]]],
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(0,1,0),placeType='jitter',
                                  meshFile=wrkDirMesh+'/LipOut262',#meshObject=lipidmeshout,
                                  **dict["LipidOut2"]
                                  )
    rSurf["Sphere"].addIngredient(cyl26IngrO2)
    
    # Inner Leaflet Lipids as Cyninders
    cyl26IngrI2 = MultiCylindersIngr(MSca*10.,  pdb='LipIn262',
                                 name='LipIn262', radii=[[4.5]],
                                  positions=[[[0,-25.5,0]]], positions2=[[[0,+5.,0]]],
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(0,1,0),placeType='jitter',
                                  meshFile=wrkDirMesh+'/LipIn262',#meshObject=lipidmeshin,
                                  **dict["LipidIn2"]
                                  )
    rSurf["Sphere"].addIngredient(cyl26IngrI2)

#===============================================================================
# surface recipes 
#===============================================================================
cylCoord1 = [ (-35,0,0), (45,0,0),  (45,0,0),    (85, 40, 0),]
cylCoord2 = [ (45,0,0),  (85,40,0), (105, -40,0), (135, 40, 30)]
cylRadii = [25, 20, 16, 16]
#fat1 = helper.getObject("Transporter")
#tall = helper.getObject("Docker")
#br = helper.getObject("Branch")
if dosurface1:
    cyl4Ingr1 = MultiCylindersIngr(MSca*.04,  pdb='1CYL1',
                                  name='Cylinders1_4', radii=[cylRadii],
                                  positions=[cylCoord1], positions2=[cylCoord2],
                                  #meshFile=wrkDirMesh+'/4Cylinders1',#meshObject=br, THE MESH IS WRONG NEED TO BE CHANGED
                                  meshFile=httpwrkDir+"geoms/Cylinders1_4.dae",
                                  principalVector=principalVector,
                                  placeType="jitter",
                                  **dict["Brancher1"]
                                  )
    rSurf["Sphere"].addIngredient(cyl4Ingr1)
    surf2 = MultiCylindersIngr(MSca*.01,  pdb=None,
                                  name='siTransporter', radii=[[50.]],
                                  positions=[[[0.,-30.,0]]], positions2=[[[0.,30.,0]]],
                                  meshFile=wrkDirMesh+'/siTransporter',#meshObject=fat1,
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(0,1,0),
                                  placeType="jitter",
                                  **dict["Transporter1"]
                                  )
    rSurf["Sphere"].addIngredient(surf2)
    
    surf3 = MultiCylindersIngr(MSca*.02,  pdb=None,
                                  name='siDocker', radii=[[50.]],
                                  positions=[[[0,-30.,0]]], positions2=[[[0.,30,0]]],
                                  meshFile=wrkDirMesh+'/siDocker',#meshObject=tall,
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(0,1,0),
                                  placeType="jitter",#jitter?
                                  **dict["Docker1"]
                                  )
    rSurf["Sphere"].addIngredient(surf3)

if dosurface2:
    cyl4Ingr2 = MultiCylindersIngr(MSca*.04,  pdb='1CYL2',
                                 name='Cylinders2_4', radii=[cylRadii],
                                  positions=[cylCoord1], positions2=[cylCoord2],
                                  #meshFile=wrkDirMesh+'/4Cylinders2',#meshObject=br,THE MESH IS WRONG NEED TO BE CHANGED
                                  meshFile=httpwrkDir+"geoms/Cylinders2_4.dae",
       principalVector=principalVector,
                                  **dict["Brancher2"]
                                  )
    rSurf["Ellipse"].addIngredient(cyl4Ingr2)
    
    esurf2 = MultiCylindersIngr(MSca*.015,  pdb=None,
                                  name='iTransporter', radii=[[50.]],
                                  positions=[[[0,-30.,0]]], positions2=[[[0.,30,0]]],
                                  meshFile=wrkDirMesh+'/iTransporter',#meshObject=fat1,
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(0,1,0),
                                  placeType="jitter",
                                  **dict["Transporter2"]
                                  )
    rSurf["Ellipse"].addIngredient(esurf2)
    
    esurf3 = MultiCylindersIngr(MSca*.02,  pdb=None,
                                  name='iDocker', radii=[[50.]],
                                  positions=[[[0,-30.,0]]], positions2=[[[0.,30,0]]],
                                  meshFile=wrkDirMesh+'/iDocker',#meshObject=tall,
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(0,1,0),
                                  placeType="jitter",#jitter?
                                  **dict["Docker2"]
                                  )
    rSurf["Ellipse"].addIngredient(esurf3)


#===============================================================================
# Previous ingredient from Host
#===============================================================================
if doPrev :
    dict["prevIngr"] = {"packingPriority":0,
                        "nMol":0,"color":None}
    prevIngr = helper.getObject("Prev")
    #need to inverse x/z since the object come from c4d
    cylCoord1 = [ (0,0,0), (0,70,-15), ]
    cylCoord2 = [ (0,45,0),  (0,70,60),]
    #center is not at first cylinder
    cylCoord1 = [ (0,-47,-30), (0,20,-38), ]
    cylCoord2 = [ (0,0,-30),  (0,20,40),]    
    cylRadii = [30, 30, ]
#    h1cyl = helper.oneCylinder("h1",cylCoord1[0],cylCoord2[0],radius=cylRadii[0])    
#    h2cyl = helper.oneCylinder("h2",cylCoord1[1],cylCoord2[1],radius=cylRadii[1])    
    previngr = MultiCylindersIngr( VesK*.002,  pdb=None,
                                 name='prevI',
                                 radii=[cylRadii],
                                 positions=[cylCoord1], 
                                 positions2=[cylCoord2],
                                 principalVector=(0,1,0),
                                 meshFile=wrkDirMesh+'/prevI',#meshObject=prevIngr,
                                 placeType=mode,
                                 #nMol=10,
                                 **dict["prevIngr"]
                                 )
    rCyto.addIngredient(previngr)
    #previngr.completion = 2.0
    #previngr.counter  = 0
    #prevI = helper.getObject("PreviousIngredientsOrOrganelles") #the parent of all instance
    #previngrInstance = helper.getChilds(prevI)
    #previngr.mesh_3d = prevIngr
    #this ingredient should be used for updating the distance and freepoint
    #the ingr should be add after the grid have been build ?
    #appendIngrInstance
#    afviewer.displayIngrResults(previngr,doSphere=True,doMesh=True)    
#    res = [afviewer.collectResult(ingr,pos,rot) for pos, rot, ingr, ptInd in h1.molecules if ingr == previngr]
#    afviewer.displayIngrCylinders(previngr,{previngr:verts},{previngr:radii},visible=1)


#afviewer.displayIngrResults(esurf3,doSphere=True,doMesh=False)
# vesicle
from DejaVu.IndexedPolygons import IndexedPolygonsFromFile

# create HistoVol
h1 = Environment()
h1.setExteriorRecipe(rCyto)

#display the organel, the box, and prepare the hierachy...

#===============================================================================
# Organelles Setup
#===============================================================================
oname = ["Ellipse","Sphere"]
for n in oname :
#    geomS = IndexedPolygonsFromFile(localdir+os.sep+"geometries"+os.sep+n, 'vesicle')
#    faces = geomS.getFaces()
#    vertices = geomS.getVertices()
#    vnormals = geomS.getVNormals()
#    o1 = Organelle(n,vertices, faces, vnormals)
    o1 = Organelle(n,None, None, None,
               filename="http://www.grahamj.com/autofill/autoFillData/GeometricCellScape/organelleGeometries/"+n)
    o1.overwriteSurfacePts = True
    h1.addOrganelle(o1)
    if rSurf.has_key(n): 
        if rSurf[n].ingredients:
            r  = rSurf[n]
            o1.setSurfaceRecipe(r)
    if rMatrix.has_key(n):
        if rMatrix[n].ingredients:
            r = rMatrix[n]
            o1.setInnerRecipe(r)

#define the viewer type dejavu,c4d,blender
h1.name="GeometricCellScape"
afviewer = AFViewer(ViewerType=helper.host,helper=helper)#long ?

#make some option here 
afviewer.doPoints = True
afviewer.doSpheres = False
afviewer.quality = 1 #lowest quality for sphere and cylinder
afviewer.visibleMesh = True #mesh default visibility 

h1.setMinMaxProteinSize()
print 'Cyto', rCyto.getMinMaxProteinSize()
#print 'Surf', rSurf1.getMinMaxProteinSize()
#print 'Matrix', rMatrix1.getMinMaxProteinSize()
#print 'o1', o1.getMinMaxProteinSize()
print 'smallest', h1.smallestProteinSize
print 'largest', h1.largestProteinSize
h1.smallestProteinSize = 15#15
print 'smallest via Override', h1.smallestProteinSize
print 'largest via Override', h1.largestProteinSize
#print o1.innerRecipe#

pad = 20.

afviewer.SetHistoVol(h1,pad,display=False)

h1.host=helper.host

h1.pickWeightedIngr = True ##do we sort the ingrediant or not see  getSortedActiveIngredients
h1.pickRandPt = True ##point pick randomly or one after the other?

h1.placeMethod = "jitter"#"spring" #"sphere"#"spring" #or "sphere"#
h1.overwritePlaceMethod = False #do we overwrtite all ingr place method-> Broke the fillling that end in infinite loop
h1.simulationTimes = 30 #number of time setTime(i/fps) is called
h1.windowsSize = 50. #radius + histoV	ol.largestProteinSize + spacing + windowsSize
h1.windowsSize_overwrite = False
if h1.windowsSize_overwrite :
    h1.windowsSize = actine.influenceRad
    #size of the windows to check for neighbours
	
h1.runTimeDisplay = False #interactive display
#Specify here the option for the dynamic if method is spring

#h1.EnviroOnly = False
#h1.EnviroOnlyCompartiment  =  -1
h1.overwriteSurfacePts = True #no need to discretisize the surface mesh as its already really denses
h1.ingrLookForNeighbours = False
h1.hackFreepts = False #SpeedUp?  Must Refer to getPointToDropHack from Fill4

h1._timer = False # Verbose for timing every function
h1._hackFreepts = False#True#bool(args[1])  # Strong hack that will never update the grids... good for very sparse fill!
h1._freePtsUpdateThrehod = 0.01#float(args[3])  # If my object covers more than 1% of 


if h1.placeMethod == "spring":
    h1.SetSpringOptions(rlength=0.1,stifness = 0.1,damping = 1.0)
    # OR
    #h1.springOptions["stifness"] = 1.
    #h1.springOptions["rlength"] = 0.
    #h1.springOptions["damping"] = 1.
    h1.SetRBOptions(obj="moving",shape="auto",child=False,
                    dynamicsBody="on", dynamicsLinearDamp = 1.0, 
                    dynamicsAngularDamp=1.0, 
                    massClamp = .001, rotMassClamp=0.1)
    h1.SetRBOptions(obj="static",shape="auto",child=False,
                    dynamicsBody="off", dynamicsLinearDamp = 0.0, 
                    dynamicsAngularDamp=0.0, 
                    massClamp = 100., rotMassClamp=1.0)
    #OR using the dictionary directly h1.dynamicOptions["moving"] and h1.dynamicOptions["static"]
    h1.SetRBOptions(obj="spring",shape="auto",child=False,
                    dynamicsBody="on", dynamicsLinearDamp = 0.0, 
                    dynamicsAngularDamp=0.0, 
                    massClamp = 10., rotMassClamp=1.0)
#    h1.SetRBOptions(obj="surface",shape="auto",child=False,
#                    dynamicsBody="on", dynamicsLinearDamp = 0.0, 
#                    dynamicsAngularDamp=0.0, 
#                    massClamp = 10., rotMassClamp=1.0)
    
h1.saveResult = False
#resultfilename = h1.resultfile =AutoFill.RECIPES["GeometricCellScape"]["resultfile"]
#resultfilename = h1.resultfile =wrkDir["GeometricCellScape"]["Resultfile"]


afviewer.displayPreFill()
afviewer.printIngrediants()

#execfile(plgDir+'/extension/testAF/c_displayPreFill.py')


try :
    AFGui.Set("GeometricCellScape",helper=afviewer.helper,afviewer=afviewer,histoVol=h1,bbox=bbox)
except:
    print ("no GUI")

doPrev = False

def FILLT(h):
   import thread
   doc =helper.getCurrentScene()
   #box = doc.get_selection()[0]
   box = doc.GetSelection()[0]
   bb=helper.getCornerPointCube(box)
   h.buildGrid(boundingBox=bb)
   #t1 = time()
   t=thread.start_new_thread(h.fill3,(14,))
   #h.fill3(seedNum=14)
   #t2 = time()
   print 'time to fill', t2-t1
   return t

def DISPLAY():
   t2 = time()
   afviewer.displayFill()
   afviewer.vi.toggleDisplay(afviewer.bsph,False)
   print 'time to display', time()-t2


    
def FILL(h,seed=20,forceBuild=True):
    t1 = time()
    if doPrev :
        previngr.completion = 1.1
        previngr.nbMol = 0
        previngr.counter = 0
    h.fill5(seedNum=seed,verbose=4)
    t2 = time()
    print 'time to fill', t2-t1
    afviewer.displayFill()
    print 'time to display', time()-t2
    afviewer.vi.toggleDisplay(afviewer.bsph,False)
    #execfile(plgDir+'/extension/testAF/c_displayFill.py')
    #afviewer.showHide(afviewer.undspMesh)
    afviewer.displayIngrResults(cyl4Ingr1,doSphere=True,doMesh=False)
    afviewer.displayIngrResults(cyl4Ingr2,doSphere=True,doMesh=False)

def GRID(h,forceBuild=True,fill=True):
    doc =helper.getCurrentScene()
    #box = doc.get_selection()[0]
    if bbox is None :
        box=helper.getCurrentSelection()[0]
    else :
        box = bbox
    bb=helper.getCornerPointCube(box)
    gridFileIn=None
    gridFileOut=None
    if forceBuild :
        gridFileOut=wrkDir+"/fig1_grid"
    else :
        gridFileIn=wrkDir+"/fig1_grid"
    if doPrev :
#        d = h.runTimeDisplay
#        h.runTimeDisplay = True
        previngr.histoVol = h
        previngr.vi = afviewer.vi
        #afviewer.appendIngrInstance(previngr,sel = previngrInstance,bb=bb)
        h.buildGrid(boundingBox=bb,gridFileIn=gridFileIn, 
                  gridFileOut=gridFileOut, previousFill=True)
#        h.runTimeDisplay = d
    else :
        h.buildGrid(boundingBox=bb,gridFileIn=gridFileIn, 
                  gridFileOut=gridFileOut)
    afviewer.displayOrganellesPoints()
    #return
    #actine.updateFromBB(h.grid)
    if fill :
        FILL(h)
    
def SecondFill(h):
    h.setMinMaxProteinSize()
#    pgrid = h.grid
    doc =helper.getCurrentScene()
    #box = doc.get_selection()[0]
    box = doc.GetSelection()[0]
    bb=helper.getCornerPointCube(box)
   # if doPrev :
    #    afviewer.appendIngrInstance(previngr,sel = previngrInstance,bb=bb)
    h.buildGrid(boundingBox=bb, previousFill=True)
    t1 = time()
    h.fill4(seedNum=14,verbose=True)
    t2 = time()
    print 'time to fill', t2-t1
    afviewer.displayFill()
    print 'time to display', time()-t2
    afviewer.vi.toggleDisplay(afviewer.bsph,False)

def load(h):
    result,orgaresult,freePoint=h.load()
    h.restore(result,orgaresult,freePoint)
    afviewer.doPoints = False
    afviewer.doSpheres = False
    afviewer.quality = 1 #lowest quality for sphere and cylinder
    afviewer.visibleMesh = True #mesh default visibility 
    afviewer.displayFill()
    
def FILLload(h):
    #load_previous_result(h)
    h.setMinMaxProteinSize()
    pgrid = h.grid
    doc =helper.getCurrentScene()
    #box = doc.get_selection()[0]
    box = doc.GetSelection()[0]
    bb=helper.getCornerPointCube(box)
    h.buildGrid(boundingBox=bb, previousFill=True)
    t1 = time()
    h.fill4(seedNum=14,verbose=True)
    t2 = time()
    print 'time to fill', t2-t1
    afviewer.displayFill()
    print 'time to display', time()-t2
    afviewer.vi.toggleDisplay(afviewer.bsph,False)

#GRID(h1)  
#load(h1)
#for e in [cyl4Ingr1,surf2,surf3,cyl4Ingr2,esurf2,esurf3]
#afviewer.displayIngrGrow(actine,visible=1)#should display the spline
#afviewer.displayIngrResults(cyl4Ingr1,doSphere=True,doMesh=False)
#afviewer.displayIngrResults(cyl4Ingr2,doSphere=True,doMesh=False)
#afviewer.displayIngrResults(esurf2,doSphere=True,doMesh=True)
#afviewer.displayIngrResults(esurf3,doSphere=True,doMesh=True)
