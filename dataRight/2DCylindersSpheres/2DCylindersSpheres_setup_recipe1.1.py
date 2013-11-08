# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 15:39:00 2011

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

Name: '2DCylinders_setup_recipe.py'
@author: Graham Johnson and Ludovic Autin
"""

#execfile("/Library/MGLTools/1.5.6.up/MGLToolsPckgs/AutoFill/figure1.py")

import numpy
from time import time
tStart = time()
import sys

#autopack

from autopack.Ingredient import SingleSphereIngr, MultiSphereIngr
from autopack.Ingredient import MultiCylindersIngr,GrowIngrediant,ActinIngrediant
from autopack.Compartment import Compartment
from autopack.Recipe import Recipe
from autopack.Environment import Environment
from autopack.Gui import AutoPackGui

#Directory

import autopack
wrkDir = autopack.__path__[0]
#from  Pmv import hostappInterface
#plgDir = hostappInterface.__path__[0]
wrkDirGeom = wrkDir+"/autopackRecipeScripts/2DCylinderSphereFill/"
#DEJAVU COLORS
modelFormat = "dae"
from upy.colors import red, aliceblue, antiquewhite, aqua, goldenrod, thistle, violet, skyblue, royalblue, plum, pink, orchid, mediumvioletred, lightpink, \
     aquamarine, azure, beige, bisque, black, blanchedalmond, mediumpurple, mediumorchid, hotpink, lightskyblue, lightblue, lightcyan, dodgerblue, firebrick, \
     blue, blueviolet, brown, burlywood, cadetblue, slateblue, darkslateblue, darkviolet, darkmagenta, powderblue, palevioletred, indigo, snow,\
     chartreuse, chocolate, coral, cornflowerblue, cornsilk, \
     crimson, cyan, darkblue, darkcyan, darkgoldenrod, \
     orange, purple, deeppink, lightcoral, \
     blue, cyan, mediumslateblue, steelblue, darkcyan, lime,\
     limegreen, darkorchid, tomato, khaki, gold, magenta, green,white

helper = autopack.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper =helperClass()

#import pyubic
#helperClass = pyubic.getHelperClass()
#helper =helperClass()

#===============================================================================
# recipes
#===============================================================================

MSca = 1.0
#no Compartment
#rSurf ={"Ellipse":Recipe(),"Sphere":Recipe()}
#rMatrix = {"Ellipse":Recipe(),"Sphere":Recipe()}
rSurf ={"Ellipse":Recipe(),"Sphere":Recipe()}
rMatrix = {"ECM_box":Recipe()}
rECM = Recipe()
rCyto = Recipe()
mode = "jitter"
forcedPriority = 0
#===============================================================================
# ingredient dictionary
#===============================================================================
rgbPurple = helper.addMaterial("gPurpleMat", (0.467, 0.239, 0.972))
rgbOrange = helper.addMaterial("gOrangeMat", (0.858, 0.470, 0.106))
#helper.changeMaterialProperty(gOrangeMat, "specular_width":0.02)
rgbGrey = helper.addMaterial("gGreyMat", (0.827, 0.827, 0.827))
rgbDarkGrey = helper.addMaterial("gDarkGreyMat", (0.498, 0.498, 0.498))
rgbVeryDarkGrey = helper.addMaterial("gVeryDarkGreyMat", (0.298, 0.298, 0.298))

if forcedPriority:
    d={}
    d["10_c25"] = {"rad":25,"height":400,
                    "dic":{"packingPriority":-30,"nbMol":20,
                        "color":rgbOrange,"useLength":False}}#,"cBound":200.
    d["10_s25"] ={"rad":25,"dic":{"packingPriority":-20,"nbMol":150,"color":rgbDarkGrey,}}
    d["10_s15"] ={"rad":15,"dic":{"packingPriority":-10,"nbMol":150,"color":rgbVeryDarkGrey,}}
else:
    d={}
    d["10_c25"] = {"rad":25,"height":400,
                    "dic":{"packingPriority":0,"nbMol":20,
                        "color":rgbOrange,"useLength":False}}#,"cBound":200.

    d["10_s25"] ={"rad":25,"dic":{"packingPriority":0,"nbMol":100 ,"color":rgbDarkGrey,}}
    d["10_s15"] ={"rad":15,"dic":{"packingPriority":0,"nbMol":100,"color":rgbVeryDarkGrey,}}

#===============================================================================
# ingredient setup
#===============================================================================
#parentmesh = helper.getObject("2DCylindersSpheresObjects")
#if parentmesh is None :
#    parentmesh = helper.newEmpty("2DCylindersSpheresObjects",location=[2000.,0,0])
# START New section added by Graham on July 16, 2012 replaces section below
# This version MAY NOT be safe outside of Cinema 4D  Can we test it ???
vBaseGeometryHider = helper.getObject("BaseGeometryHider") #g
if vBaseGeometryHider is None : #g
    vBaseGeometryHider=helper.newEmpty("BaseGeometryHider") #g

# END New section added by Graham on July 16, 2012    



cyl = helper.getObject("c10_c25")
h = d["10_c25"]["height"]
#what about nbJitter
#lingr={}
#for k in d.keys():
#    m = helper.getObject(k)
#    ingr = SingleSphereIngr( MSca*.01,  d[k]["rad"],name='i'+k, pdb=None,
#                                meshObject=m,placeType=mode,**d[k]["dic"])
#    lingr[k]=ingr
#    rMatrix["ECM_box"].addIngredient( ingr )

#that is a problem here the cylinder use the scale for gettin the irght size, 
#thus instance of it dont have the correct size...
if cyl is None :
    #create a cylinder
    cyl,mesh = helper.Cylinder('c10_c25',radius=d["10_c25"]["rad"],
                               length=h,parent=vBaseGeometryHider,axis=[0.0,1.0,0.0])#axi is Y-up
#    cyl = helper.oneCylinder('10_c25', 
#                            [0,-h/2.,0],[0,h/2.,0],
#                            radius=d["10_c25"]["rad"])
ingcyl = MultiCylindersIngr(MSca*.001,  pdb=None,
                              name='c10_c25', radii=[[d["10_c25"]["rad"]]],
                              positions=[[[0,-h/2.,0]]], positions2=[[[0,h/2.,0]]],
                              meshObject=cyl,
                              #meshFile=wrkDirGeom+"/cylinder."+modelFormat,
                              jitterMax=(1.,1.,0.),
                            
                            #   "useRotAxis":{"name":"useRotAxis","value":False,"default":False,"type":"bool","min":0.,"max":0.,"description":"useRotAxis"},
                            #"rotAxis":{"name":"rotAxis","value":[0.,0.,0.],"default":[0.,0.,0.],"min":0,"max":1,"type":"vector","description":"rotAxis"},
                            #   rotAxis=(0.,3.,1.0),
                              rotAxis = (0.,2.,1.),
                              useRotAxis = 1,
                              nbJitter = 10,
                              #  CRITICAL !!! IF jitter is greater than radius of object, 
                              #e.g. 5x(1,1,1) the point may not be consumed!!!
                              principalVector=(0,1,0),
                              placeType=mode,
                              **d["10_c25"]["dic"]
                              )
rMatrix["ECM_box"].addIngredient( ingcyl )
#rCyto.addIngredient(ingcyl)
ingcyl.rotAxis = (0.,0.,1.)
ingcyl.useRotAxis = 1
ingcyl.rejectionThreshold = 20000
#ingcyl.length = 2*h
ingcyl.rotRange = 6.2831#1.5707
#3.14159#0.1745#10-1.5707#90-6.2831-360
sph = helper.getObject("s10_s25")
if sph is None :
    sph,mesh = helper.Sphere("s10_s25",radius=d["10_s25"]["rad"],res=12,parent=vBaseGeometryHider)
ingsph = SingleSphereIngr( MSca*.015,  d["10_s25"]["rad"],
                            name='s10_s25', pdb=None,
                            meshObject=sph,
                            placeType=mode,
                            jitterMax=(1.,1.,0.),
                            **d["10_s25"]["dic"]
                            #packingMode='close'
                            ) #original radius is 3.61
#rCyto.addIngredient(ingsph)
rMatrix["ECM_box"].addIngredient( ingsph )
ingsph.rejectionThreshold = 2500000

sph3 = helper.getObject("s10_s15")
if sph3 is None :
    sph3,mesh = helper.Sphere("s10_s15",radius=d["10_s15"]["rad"],res=12,parent=vBaseGeometryHider)
ingsph3 = SingleSphereIngr( MSca*.015,  d["10_s15"]["rad"],
                            name='s10_s15', pdb=None,
                            meshObject=sph3,
                            jitterMax=(1.,1.,0.),
                            placeType=mode,
                            **d["10_s15"]["dic"]
                            #packingMode='close'
                            ) #original radius is 3.61
#rCyto.addIngredient(ingsph)
rMatrix["ECM_box"].addIngredient( ingsph3 )
ingsph.rejectionThreshold = 2500000

#===============================================================================
# viewer setup
#===============================================================================
#ViewerType='c4d'
ViewerType=helper.host
afviewer = AutopackViewer(ViewerType=ViewerType,helper=helper)#long ?
#make some option here 
afviewer.doPoints = True
afviewer.doSpheres = False
afviewer.quality = 1 #lowest quality for sphere and cylinder
afviewer.visibleMesh = True #mesh default visibility 
afviewer.name = "CylindersSpheres2D"
#===============================================================================
# Compartments setup
#===============================================================================

#from DejaVu.IndexedPolygons import IndexedPolygonsFromFile

# create Environment
h1 = Environment()
#h1.setExteriorRecipe(rMatrix["ECM_box"])



#afviewer.setupCompartment(h1,rSurf,rMatrix)

#===============================================================================
# User selected bounding box setup
#===============================================================================
#create the box
#bbox = afviewer.helper.getObject("EnvironmentBB")
#if bbox is None :
#    if helper.host == "c4d" :
#    import c4d
#    arr1 = c4d.BaseObject(c4d.Oatomarray)
#    arr1.SetName('BoundingBox')
#    arr1[1000] = 1.0 #radius cylinder
#    arr1[1001] = 1.0 #radius sphere
#    arr1[1002] = 3 #subdivision
#    helper.AddObject(arr1)
#    mat2 = helper.getMaterial("white")
#    helper.assignMaterial(arr1, mat2, texture=True) # I figured out how to assign a material!!!!
#    bbox = afviewer.helper.box("EnvironmentBB",center=[0.,0.,0.],size=[1000.,1000.,6.], parent=arr1)[0] #remember x and z will swap for left-handed hosts or files
#    else :
#   bbox = afviewer.helper.box("EnvironmentBB",center=[0.,0.,0.],size=[1000.,1000.,6.])[0] #remember x and z will swap for left-handed hosts or files

pad = [50,50,0]

#h1.boundingBox = [[-550.,-550,-0.5],[550,550,0.5]]
#bbox = afviewer.helper.getObject("fillBB")
#if bbox is None : bbox = afviewer.helper.box("fillBB", cornerPoints=[[-500.,-500,-0.5],[500,500,0.5]])#cornerPoints=h1.boundingBox)
#print "exteriorCompartment"
exteriorCompartment = afviewer.helper.getObject("fillBB")
if exteriorCompartment is None : exteriorCompartment = afviewer.helper.box("fillBB", cornerPoints=[[-400.,-400,-0.5],[400,400,0.5]])#cornerPoints=h1.boundingBox)

#pad = 0.
#
#h1.boundingBox = [[-550.,-550,-0.5],[550,550,0.5]]
#bbox = afviewer.helper.getObject("fillBB")
#if bbox is None : bbox = afviewer.helper.box("fillBB", cornerPoints=[[-500.,-500,-0.5],[500,500,0.5]])#cornerPoints=h1.boundingBox)
Compartment = False
if Compartment:
#    h1.setExteriorRecipe(rECM)
    o1 = Compartment("exteriorCompartment1",None, None, None,isOrthogonalBoudingBox=1)
    o1.bb = helper.getCornerPointCube(afviewer.helper.getObject("fillBB"))
    h1.addCompartment(o1)
    o1.setInnerRecipe(rMatrix["ECM_box"])
else :
    h1.setExteriorRecipe(rMatrix["ECM_box"])
#===============================================================================
# Environmentume options
#===============================================================================
h1.setMinMaxProteinSize()
#print 'Cyto', rCyto.getMinMaxProteinSize()
#print 'Surf', rSurf1.getMinMaxProteinSize()
#print 'Matrix', rMatrix1.getMinMaxProteinSize()
#print 'o1', o1.getMinMaxProteinSize()
#print 'smallest', h1.smallestProteinSize
#print 'largest', h1.largestProteinSize
h1.smallestProteinSize = 10#15
#print 'smallest via Override', h1.smallestProteinSize
#print 'largest via Override', h1.largestProteinSize
#print o1.innerRecipe#

#pad = [0,0,0]
#h1.boundingBox = [[0.,0.,0.],[1000.,1000.,1.]]
h1.name="Test_CylindersSpheres2D"
h1.encapsulatingGrid = 0  #This value can only sapely be set to 0 instead of 1 if doing a 2D restricted fill.
afviewer.SetHistoVol(h1,pad,display=False)



#
bbox = afviewer.helper.getObject("histoVolBB")
##print ("bbox",bbox)
if bbox is None : bbox = afviewer.helper.box("histoVolBB",cornerPoints=[[-500.,-500.,-.5],[500.,500.,.5]])#cornerPoints=h1.boundingBox)
print ("bbox",bbox)

#h1.host='c4d'
#test = [[False,True],[True,True],[True,False],[False,False]]
test = [[False,False],[True,False],[False,True],[True,True]]

testid = 3

h1.pickWeightedIngr = test[testid][0] ##do we sort the ingrediant or not see  getSortedActiveIngredients
h1.pickRandPt = test[testid][1] ##point pick randomly or one after the other?


h1.placeMethod = "jitter"#"spring" #"sphere"#"spring" #or "sphere"#
h1.overwritePlaceMethod = True #do we overwrtite all ingr place method
h1.simulationTimes = 30 #number of time setTime(i/fps) is called
h1.windowsSize = 10. #radius + Environment.largestProteinSize + spacing + windowsSize
h1.windowsSize_overwrite = False
if h1.windowsSize_overwrite :
    h1.windowsSize = actine.influenceRad
    #size of the windows to check for neighbours
h1.runTimeDisplay = 0 #interactive display
#Specify here the option for the dynamic if method is spring

#h1.EnviroOnly = False
#h1.EnviroOnlyCompartiment  =  -1
h1.overwriteSurfacePts = False
h1.ingrLookForNeighbours = False

h1.overwriteSurfacePts = False
h1.ingrLookForNeighbours = False

h1.hackFreepts = 0 #SpeedUp?  Must Refer to getPointToDropHack from Fill4
h1._timer = False # Verbose for timing every function
h1._hackFreepts = 0 #True#bool(args[1])  # Strong hack that will never update the grids... good for very sparse fill!
h1._freePtsUpdateThrehod = 0.0#float(args[3])  # If my object covers more than 1% of the remaining freepoints, update the perIngredientAvailablePoints Lists

#h1.gSeededTable = numpy.random.RandomState(14)  #want to make one random table from a global seed set by user or defaults


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
    
h1.saveResult = True
#resultfilename = h1.resultfile = "/Users/Shared/fillResult.af"
resultfilename = h1.resultfile =wrkDir+os.sep+"cache_results/CylSpherefillResult.apr"
print ("names are", h1.name, afviewer.name)
afviewer.displayPreFill()
afviewer.printIngrediants()

helper.toggleDisplay(vBaseGeometryHider,False)

#===============================================================================
# functions
#===============================================================================
try :
    AFGui.Set("Test_CylindersSpheres2D",helper=afviewer.helper,afviewer=afviewer,histoVol=h1,bbox=bbox)
except:
    print ("no GUI")

#===============================================================================
# functions
#===============================================================================

def GRID(h,forceBuild=True,fill=True):
    t1 = time()
    doc =helper.getCurrentScene()
    #box = doc.get_selection()[0]
    box = helper.getObject("histoVolBB")#doc.GetSelection()[0]
    #box = doc.GetSelection()[0]
    bb=helper.getCornerPointCube(box)
    gridFileIn=None
    gridFileOut=None
    if forceBuild :
        gridFileOut="/Users/Shared/fill_grid"
    else :
        gridFileIn="/Users/Shared/fill_grid"
    h.buildGrid(boundingBox=bb)#,gridFileIn=gridFileIn, 
                  #gridFileOut=gridFileOut)
    t2 = time()
    gridTime = t2-t1
    if fill :
        FILL(h)
#    print 'time to Build Grid', gridTime
    afviewer.displayCompartmentsPoints()
    #return
    #actine.updateFromBB(h.grid)

def FILL(h,seed=0,forceBuild=True):
    t1 = time()
    # seed random number generator
#    seed(seedNum)
#    self.randomRot.setSeed(seed=seedNum)
    h.fill5(seedNum=seed,verbose=0, vTestid = testid)
    t2 = time()
#    print 'time to run Fill5', t2-t1
    afviewer.displayFill()
#    print 'time to display Fill5', time()-t2
    afviewer.vi.toggleDisplay(afviewer.bsph,False)
    #execfile(plgDir+'/extension/testAF/c_displayFill.py')
    #afviewer.showHide(afviewer.undspMesh)

def SecondFill(h):
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
#    print 'time to fill', t2-t1
    afviewer.displayFill()
#    print 'time to display', time()-t2
    afviewer.vi.toggleDisplay(afviewer.bsph,False)

#GRID(h1)
#tEnd = time()
#print 'time to fill', tEnd-tStart
#FILL(h1,seed=25)
