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

Name: '2DSpheres_setup_recipe.py'
@author: Graham Johnson and Ludovic Autin
"""

# To properly run this file, Open ResultFigure2UnwantedGradients3.c4d
# Delete everything in the hierarchy including and above the "Cytoplasm" null object
# Select the parent of ECM_box called "OrganelleSelector"
# Adjust the test ID and seed in the fill() function, etc.
# Set the following line to your correct file path for this .py file, paste it into the C4D python console and hit enter
# execfile("/Users/grahamold/Dev/autoFillSVN/autofillSVNversions/trunk/AutoFillClean/autoFillRecipeScripts/figure2.py")
# to hide the parent meshes, use the magnifying glass in C4D's object manager:
#    type meshsparent into the box, select all the objects, group them, find the group's parent null and make it invisible

#execfile("/Library/MGLTools/1.5.6.up/MGLToolsPckgs/AutoFill/figure1.py")
#execfile("/Users/ludo/DEV/Autofill_svn_test/autofill/trunk/AutoFillClean/autoFillRecipeScripts/2DsphereFill/2DSpheres_setup_recipe.py")
import sys
import os

#import numpy
from time import time
tStart = time()
#import sys
#sys.path.append("/Users/grahamold/Dev/autoFillSVN/autofillSVNversions/trunk/AutoFillClean/")

#AUTOFILL
import AutoFill
localdir = wrkDir = AutoFill.__path__[0]

from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
from AutoFill.Ingredient import MultiCylindersIngr,GrowIngrediant,ActinIngrediant
from AutoFill.Organelle import Organelle
from AutoFill.Recipe import Recipe
from AutoFill.HistoVol import Environment
from AutoFill.autofill_viewer import AFViewer



from upy.colors import red, aliceblue, antiquewhite, aqua, goldenrod, thistle, violet, skyblue, royalblue, plum, pink, orchid, mediumvioletred, lightpink, \
     aquamarine, azure, beige, bisque, black, blanchedalmond, mediumpurple, mediumorchid, hotpink, lightskyblue, lightblue, lightcyan, dodgerblue, firebrick, \
     blue, blueviolet, brown, burlywood, cadetblue, slateblue, darkslateblue, darkviolet, darkmagenta, powderblue, palevioletred, indigo, snow,\
     chartreuse, chocolate, coral, cornflowerblue, cornsilk, \
     crimson, cyan, darkblue, darkcyan, darkgoldenrod, \
     orange, purple, deeppink, lightcoral, \
     blue, cyan, mediumslateblue, steelblue, darkcyan, lime,\
     limegreen, darkorchid, tomato, khaki, gold, magenta, green,grey,white, lightgrey, darkgrey

helper = AutoFill.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper =helperClass()
print (1,time()-tStart)
tStart = time()
#===============================================================================
# recipes
#===============================================================================

MSca = 1.0
#no organelle
#rSurf ={"Ellipse":Recipe(),"Sphere":Recipe()}
#rMatrix = {"Ellipse":Recipe(),"Sphere":Recipe()}
rCyto = Recipe()
mode = "jitter"
forcedPriority = 0
#===============================================================================
# ingredient dictionary
#===============================================================================
print('******* Pre HelperHost check ***********')
vColorRGB = None
#if helper.host == "c":
#    print ('helper host = c')
#    useHostGeom = 0
#    if useHostGeom :
#        pass
rgbRed = helper.addMaterial("gRedMat", (0.784, 0.204, 0.204))
rgbGrey = helper.addMaterial("gGreyMat", (0.827, 0.827, 0.827))
rgbWhite = helper.addMaterial("gWhiteMat", (0.95, 0.95, 0.95))
rgbDarkGrey = helper.addMaterial("gDarkGreyMat", (0.498, 0.498, 0.498))
#helper.changeMaterialProperty(rgbDarkGrey, specular_width = 0.20)
rgbBlue = helper.addMaterial("gBlueMat", (0.306, 0.451, 0.816))
#helper.changeMaterialProperty(rgbBlue, specular_width = 0.20)
rgbPurple = helper.addMaterial("gPurpleMat", (0.467, 0.239, 0.972))
#rgbRed = (0.7843, 0.2039, 0.2039)
print (22,time()-tStart)
tStart = time()

if forcedPriority:
    d={}
    d["5_n200"] =   {"rad":200,"molarity":0.001,"dic":   {"packingPriority":-50,"nbMol":4, "color": rgbWhite}}
    d["10_n100"] =  {"rad":100,"molarity":0.01,"dic":   {"packingPriority":-40,"nbMol":8,"color": rgbGrey}}
    d["200_n50"] =  {"rad":50,"molarity":0.1,"dic":    {"packingPriority":-30,"nbMol":16,"color": rgbBlue}}
    #d["400_n25"] = {"rad":25,"dic":    {"packingPriority":0,"nMol":60,"color":None}}
    d["x400_n25"] = {"rad":25,"molarity":1.0,"dic":    {"packingPriority":-20,"nbMol":32,"color": rgbPurple}}
    d["y400_n25"] = {"rad":25,"molarity":10.,"dic":    {"packingPriority":-20,"nbMol":32,"materials": red}}
#    d["a400_n25"] = {"rad":25,"dic":    {"packingPriority":0,"nMol":16,"color":None}}
#    d["b400_n25"] = {"rad":25,"dic":    {"packingPriority":0,"nMol":16,"color":None}}
else:
    d={}
    #    d["5_n200"] =   {"rad":200,"dic":   {"packingPriority":0,"nMol":4,"color": lightgrey}}
    d["5_n200"] =   {"rad":200,"molarity":0.001,"molarity":0.000,"dic":   {"priority":0,"nbMol":4, "color": rgbGrey}}
    d["10_n100"] =  {"rad":100,"molarity":0.001,"dic":   {"priority":0,"nbMol":8,"color": rgbDarkGrey}}
    d["200_n50"] =  {"rad":50,"molarity":0.01,"dic":    {"packingPriority":0,"nbMol":16,"color": rgbBlue}}
#    d["400_n25"] =  {"rad":25,"dic":    {"packingPriority":0,"nMol":60,"color":None}}
    d["x400_n25"] = {"rad":25,"molarity":0.1,"dic":    {"packingPriority":0,"nbMol":32,"color": red}}
    d["y400_n25"] = {"rad":25,"molarity":1.0,"dic":    {"packingPriority":0,"nbMol":48, "color": rgbPurple}}
#    d["a400_n25"] = {"rad":25,"dic":    {"packingPriority":0,"nMol":16,"color":None}}
#    d["b400_n25"] = {"rad":25,"dic":    {"packingPriority":0,"nMol":16,"color":None}}

    #    d["5_n200"] =   {"rad":200,"dic":   {"packingPriority":0,"nMol":4,"color": lightgrey}}
    #    d["10_n100"] =  {"rad":100,"dic":   {"packingPriority":0,"nMol":8,"color": darkgrey}}
    #    d["200_n50"] =  {"rad":50,"dic":    {"packingPriority":0,"nMol":16,"color": slateblue}}
    #    d["x400_n25"] = {"rad":25,"dic":    {"packingPriority":0,"nMol":32,"color": darkviolet}}
#    d["y400_n25"] = {"rad":25,"dic":    {"packingPriority":0,"nMol":132, "color": firebrick}}

#packingmode = "gradient"#close"
#grname = "gY"
#grname1 = "gY"
#grname2 = "gX"

print (2,time()-tStart)
tStart = time()
#parentmesh = helper.getObject("2DSpheresObjects") 
    #if parentmesh is None :
#    parentmesh = helper.newEmpty("2DSpheresObjects",location=[2000.,0,0])

#===============================================================================
# ingredient setup
#===============================================================================
# START New section added by Graham on July 16, 2012 replaces section below
# This version MAY NOT be safe outside of Cinema 4D  Can we test it ???
vBaseGeometryHider = helper.getObject("BaseGeometryHider") #g
if vBaseGeometryHider is None : #g
    vBaseGeometryHider=helper.newEmpty("BaseGeometryHider") #g

showFillGeoms=1
# END New section added by Graham on July 16, 2012   


lingr={}
for k in list(d.keys()):
#    if helper.host != "c":
    m = helper.getObject(k)
    if m is None :
       m,mesh = helper.Sphere(k,radius=d[k]["rad"],res=12,parent=vBaseGeometryHider)
       helper.toggleDisplay(m,showFillGeoms)	
#    else :
#       m = helper.getObject(k)
    ingr = SingleSphereIngr(MSca*0.0001, d[k]["rad"],name='i'+k, pdb=None,#d[k]["molarity"],
                          meshObject=m,placeType=mode,**d[k]["dic"])
    
    lingr[k]=ingr
    rCyto.addIngredient(ingr)
print (3,time()-tStart)
tStart = time()
#if grname1 is not None :
#    lingr["x400_n25"].gradient = grname2
#    lingr["x400_n25"].packingMode=packingmode
#    lingr["y400_n25"].gradient = grname1
#    lingr["y400_n25"].packingMode=packingmode

helper.toggleDisplay(vBaseGeometryHider,False)

#===============================================================================
# viewer setup
#===============================================================================
afviewer = AFViewer(ViewerType=helper.host,helper=helper)#long ?
#make some option here 	
afviewer.doPoints = True
afviewer.doSpheres = False
if helper.host == "c": afviewer.doSpheres = False
afviewer.quality = 1 #lowest quality for sphere and cylinder
afviewer.visibleMesh = True #mesh default visibility 
print (4,time()-tStart)
tStart = time()

#===============================================================================
# Organelles setup
#===============================================================================

#from DejaVu.IndexedPolygons import IndexedPolygonsFromFile

# create HistoVol
h1 = Environment()
print (554,time()-tStart)
tStart = time()

h1.setExteriorRecipe(rCyto)
h1.name="Test_Spheres2D"
#afviewer.setupOrganelle(h1,rSurf,rMatrix)
#there is no organelle here
print (555,time()-tStart)
tStart = time()

#===============================================================================
# histoVolume options
#===============================================================================
h1.setMinMaxProteinSize()
print('Cyto', rCyto.getMinMaxProteinSize())
#print 'Surf', rSurf1.getMinMaxProteinSize()
#print 'Matrix', rMatrix1.getMinMaxProteinSize()
#print 'o1', o1.getMinMaxProteinSize()
print('smallest', h1.smallestProteinSize)
print('largest', h1.largestProteinSize)
h1.smallestProteinSize = 15#15
print('smallest via Override', h1.smallestProteinSize)
print('largest via Override', h1.largestProteinSize)
#print o1.innerRecipe#

pad = 0.

#h1.boundingBox = [[-550.,-550,-0.5],[550,550,0.5]]
bbox = afviewer.helper.getObject("histoVolBB")
if bbox is None : bbox = afviewer.helper.box("histoVolBB", cornerPoints=[[-500.,-500,-0.5],[500,500,0.5]])[0]#cornerPoints=h1.boundingBox)
#if bbox is None : bbox = afviewer.helper.box("histoVolBB",center=[0.,0.,0.],size=[1000.,1000.,6.])
#bbox.cornerPoints=helper.getCornerPointCube(bbox)

h1.encapsulatingGrid = 0  #This value can only safely be set to 0 instead of 1 if doing a 2D restricted fill.

#if grname is not None :
#    h1.setGradient(name=grname1,mode="direction",direction=[0.5,0.5,0.0])
#    h1.setGradient(name=grname2,mode="X")
#    h1.use_gradient=True
#    print h1.gradients
#in that case the packing mode need to be "gradient"
#each ingredient with packingmode as gradient need ingr.gradient = nameofthegradient

#h1.name="Spheres2D"
afviewer.SetHistoVol(h1,pad,display=False)
print (5,time()-tStart)
tStart = time()

h1.host=helper.host

#test = [[False,True],[True,True],[True,False],[False,False]]
test = [[False,False],[True,False],[False,True],[True,True]]

#test 0 is 3.88
#test 1 is 3.83 filling time (no building time as there is no organelle)
#0 and 1 look exactly he same
#test 2 is longer 94.54 
#test 3 is longer 97.25
#2 and 3 exactly the same
#the difference came from the ramdom picking..need to solve the sorting problem

#longer to pick random point is really slow compare to pick next
testid = 3
h1.pickWeightedIngr = test[testid][0] ##do we sort the ingrediant or not see  getSortedActiveIngredients
h1.pickRandPt = test[testid][1] ##point pick randomly or one after the other?

h1.placeMethod = "jitter"#"spring" #"sphere"#"spring" #or "sphere"# gradient ?
h1.overwritePlaceMethod = True #do we overwrtite all ingr place method
h1.simulationTimes = 30 #number of time setTime(i/fps) is called
h1.windowsSize = 10. #radius + histoVol.largestProteinSize + spacing + windowsSize
h1.windowsSize_overwrite = False
if h1.windowsSize_overwrite :
    h1.windowsSize = actine.influenceRad
    #size of the windows to check for neighbours
h1.runTimeDisplay = 0 #interactive display
#Specify here the option for the dynamic if method is spring

#h1.EnviroOnly = False
#h1.EnviroOnlyCompartiment  =  -1

h1.overwriteSurfacePts = True
h1.ingrLookForNeighbours = False
h1._freePtsUpdateThrehod = 0.0 # If my object covers more than 1% of the remaining freepoints, update the perIngredientAvailablePoints Lists

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

print (6,time()-tStart)
tStart = time()

h1.saveResult = True
#resultfilename = h1.resultfile = wrkDir+os.sep+"autoFillRecipeScripts"+os.sep+"2DsphereFill"+os.sep+"results"+os.sep+"SpherefillResult.afr"
resultfilename = h1.resultfile =wrkDir+os.sep+"autoFillRecipeScripts/2DsphereFill/results/SpherefillResult.afr"
afviewer.displayPreFill()
afviewer.printIngrediants()
#

print (7,time()-tStart)
#===============================================================================
# User selected bounding box setup
#===============================================================================
#bbox = None
#if helper.host != "c":
    #create the box
#   bbox = helper.getObject("fillBox")
#   if bbox is None : 
#       if helper.host == "c4d" :
    #            import c4d
    #        arr1 = c4d.BaseObject(c4d.Oatomarray)
    #        arr1.SetName('BoundingBox')
    #        arr1[1000] = 1.0 #radius cylinder
    #        arr1[1001] = 1.0 #radius sphere
    #        arr1[1002] = 3 #subdivision
    #        helper.AddObject(arr1)    
    #        mat2 = helper.getMaterial("white")
    #        helper.assignMaterial(arr1, mat2, texture=True) # I figured out how to assign a material!!!!
    #        bbox = helper.box("fillBox",center=[0.,0.,0.],size=[1000.,1000.,6.], parent=arr1)[0] #remember x and z will swap for left-handed hosts or files
    #    else :
#        bbox = helper.box("fillBox",center=[0.,0.,0.],size=[1000.,1000.,6.])[0] #remember x and z will swap for left-handed hosts or files

#box = bbox
#helper.toggleRender(bbox,False)

noGui=False
#===============================================================================
# functions
#===============================================================================
try :
    AFGui.Set("Test_Spheres2D",helper=afviewer.helper,afviewer=afviewer,histoVol=h1,bbox=bbox)
except:
    print ("no GUI")
    noGui=True
print (7,time()-tStart)
tStart = time()

#try :#Deprecated
#    #call the gui
#    if mygui is not None :
#    #from AutoFill.af_script_ui import AFScriptUI
#    #mygui = AFScriptUI(title="AFUI")
#        mygui.reset()
#        mygui.Set(helper=helper,afviewer=afviewer,histoVol=h1,bbox=bbox)
#        mygui.updateWidget()
#except:
#    print ("no GUI")
    
def GRID(h,forceBuild=True, fill=0,seed=20,vTestid = 3,vAnalysis = 0):
    t1 = time()
    doc =helper.getCurrentScene()
    #ox = doc.GetSelection()[0]
    if bbox is None :
        box=helper.getCurrentSelection()[0]
    else :
        box = bbox
    bb=helper.getCornerPointCube(box)
    gridFileIn=None
    gridFileOut=None
    if forceBuild :
        gridFileOut=wrkDir+"fill_grid"
    else :
        gridFileIn=wrkDir+"fill_grid"
    h.buildGrid(boundingBox=bb,gridFileIn=None,rebuild=True ,
                      gridFileOut=None,previousFill=False)
#    h.buildGrid(gridFileIn=gridFileIn, 
#                  gridFileOut=gridFileOut)
    t2 = time()
    gridTime = t2-t1
    if fill :
        FILL(h,seed=seed,vTestid = vTestid,vAnalysis = vAnalysis)
    print ('time to Build Grid', gridTime)
#    afviewer.displayOrganellesPoints()
    #return
    #actine.updateFromBB(h.grid)

def FILL(h,seed=20,forceBuild=True,vTestid = 3,vAnalysis = 0):
    t1 = time()
    h.fill5(seedNum=seed,verbose=4, vTestid = vTestid,vAnalysis = vAnalysis)
    t2 = time()
    print('time to run Fill5', t2-t1)
    afviewer.displayFill()
    print('time to display Fill5', time()-t2)
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
    print('time to fill', t2-t1)
    afviewer.displayFill()
    print('time to display', time()-t2)
    afviewer.vi.toggleDisplay(afviewer.bsph,False)

#GRID(h1)
#FILL(h1,seed=25)
#tEnd = time()
#print ('time to fill', tEnd-tStart)
#loop of packing
def doloop(n):
    # doLoop automatically produces result files, images, and documents from the recipe while adjusting parameters
    # To run doLoop, 1) in your host's python console type:
    # execfile(pathothis recipe) # for example, on my computer:
    # " execfile("/Users/grahamold/Dev/autoFillSVN/autofillSVNversions/trunk/AutoFillClean/autoFillRecipeScripts/2DsphereFill/2DSpheres_setup_recipe.py")
    # 2) prepare your scene for the rendering->render output should be 640,480 but you can change the size in the script at then end.  Set up textures lights, and effects as you wish
    #    Results will appear in the result folder of your recipe path
    # where n is the number of loop, seed = i
    #doloop(n) 
    rangeseed=range(n)
    for i in rangeseed:
        basename = localdir+os.sep+"autoFillRecipeScripts/2DsphereFill/results"+os.sep+"results_seed_"+str(i)
        h1.saveResult = True
        resultfilename = h1.resultfile = basename  
        GRID(h1,seed=i, fill=1,vTestid = i,vAnalysis = 1)#build grid and fill
        #render/save scene
        helper.render(basename+".jpg",640,480)
        helper.write(basename+".c4d",[])
        #clean
    #    afviewer.clearAll("Test_Spheres2D")
        afviewer.clearFill("Test_Spheres2D")