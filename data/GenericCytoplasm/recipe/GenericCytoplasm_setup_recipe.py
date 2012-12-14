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

Name: 'GenericCytoplasm_setup_recipe'
@author: Graham Johnson and Ludovic Autin
"""

import sys
#AUTOFILL
import AutoFill
wrkDir=localwrkDir = AutoFill.__path__[0]
httpwrkDir = "http://grahamj.com/autofill/autoFillData/GenericCytoplasm/"
sphdir = "http://grahamj.com/autofill/autoFillData/GenericCytoplasm/CytoSpheres/"

from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
from AutoFill.Ingredient import MultiCylindersIngr, GrowIngrediant
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
     limegreen, darkorchid, tomato, khaki, gold, magenta, green, snow


MSca = 1.0

mode = "jitter"

## Surface:
rSurf1 = Recipe()

##
## Matrix
##
rMatrix1 = Recipe()

##
## Cytoplasm:
##
rCyto = Recipe()

####### BEGIN New generic Cytoplasm
#should cal the generic scrip here
setupfile=localwrkDir+os.sep+"autoFillRecipeScripts"+os.sep+"GenericCytoplasm"+os.sep+"GenericCytoplasmRecipeFile1.py"
exec(open(setupfile,"r").read(),globals(),{"wrkDir2":sphdir,"tempRecipe":rCyto,"httpwrkDir":httpwrkDir})
####### END New generic Cytoplasm


##

## Cytoplasm:
# 1ABL
dnaIngr = GrowIngrediant(MSca*1., color=coral, pdb=None,
                              name='dnaGrow', radii=[[3.],],
                              positions=[[[0,0,0]]], positions2=[[[0.,100.,0]]],
                              packingPriority=-1,biased=0.5,
                              principalVector=(0,1,0),
                              #marge = 180.0,
                              modelType="Cylinders",placeType="jitter",
                              nbJitter=1, jitterMax=(0,0,0),
                              length = 500000,closed = False,
                              nMol=2,
                              )

dnaIngr.cutoff_boundary=3.0
dnaIngr.cutoff_surface=3.0
dnaIngr.seedOnMinus = True
#rCyto.addIngredient( dnaIngr )

#
#kinase1 = SingleSphereIngr( MSca*.009,  12., color=snow, name='kinase1',
#                           meshFile=wrkDir+'/recipes/cyto/1ABL_centered', pdb='1ABL',
#                            packingPriority=1.1
#                            )
##rCyto.addIngredient( kinase1 )
#
#
#
#
#ingr1ABL = MultiSphereIngr( MSca*.009, name='1ABL', pdb='1ABL',
#                            sphereFile=wrkDir+'/recipes/cyto/1ABL_centered_1.sph',#_1
#                            meshFile=wrkDir+'/recipes/cyto/1ABL_centered',
#                            color=snow, packingPriority=1.1)#, packingMode='close')
##rCyto.addIngredient( ingr1ABL )
#
#
#
#ingr2ABL = MultiSphereIngr( MSca*.0045, name='2ABL', pdb='2ABL',
#                            sphereFile=wrkDir+'/recipes/cyto/1ABL_centered_1.sph',#_1
#                            meshFile=wrkDir+'/recipes/cyto/1ABL_centered',
#                            color=gainsboro, packingPriority=1.1)#, packingMode='close')
##rCyto.addIngredient( ingr2ABL )
#
#
#ingr1IRK = MultiSphereIngr( MSca*.004, name='1IRK', pdb='1IRK',
#                            sphereFile=wrkDir+'/recipes/cyto/1IRK_centered_7.sph',
#                            meshFile=wrkDir+'/recipes/cyto/1IRK_centered',
#                            color=ivory, packingPriority=1.3)#, packingMode='close')
##rCyto.addIngredient( ingr1IRK )
#
#
#
#ingr6TNA = MultiSphereIngr( MSca*.004, name='tRNA', pdb='6TNA',
#                            sphereFile=wrkDir+'/recipes/cyto/6TNA_centered_5.sph',
#                            meshFile=wrkDir+'/recipes/cyto/6TNA_centered',
#                            color=darkgray, packingPriority=1.3)#, packingMode='close')
#
#rCyto.addIngredient( ingr6TNA )
#
#
#
#GroelIngr1 = MultiSphereIngr( MSca*.00001, name='Groel', pdb='1AON',
#                              sphereFile=wrkDir+'/recipes/cyto/1AON_centered.sph',
#                              meshFile=wrkDir+'/recipes/cyto/1AON_centered',
#                              color=silver, packingPriority= -.08)
##rCyto.addIngredient( GroelIngr1 )
#
#
#
## 2CPK
#ingr2CPK = MultiSphereIngr( MSca*.005, name='PHOSPHOTRANSFERASE', pdb='2CPK',
#                            sphereFile=wrkDir+'/recipes/cyto/2CPK_centered_7.sph',
#                            meshFile=wrkDir+'/recipes/cyto/2CPK_centered',
#                            color=linen, packingPriority=1.3)#, packingMode='close')
#
##rCyto.addIngredient( ingr2CPK )
#
#
#
## 1CZA
#ingr1CZA = MultiSphereIngr( MSca*.0003, name='TRANSFERASE', pdb='1CZA',
#                            sphereFile=wrkDir+'/recipes/cyto/1CZA_centered_7.sph',
#                            meshFile=wrkDir+'/recipes/cyto/1CZA_centered',
#                            color=grey, packingPriority=1.5)
##rCyto.addIngredient( ingr1CZA )
#
#
#
## 2OT8
#ingr2OT8 = MultiSphereIngr( MSca*.0002, name='TRANSPORT PROTEIN', pdb='2OT8',
#                            sphereFile=wrkDir+'/recipes/cyto/2OT8_centered_11.sph',
#                            meshFile=wrkDir+'/recipes/cyto/2OT8_centered',
#                            color=dimgrey, packingPriority=1.5)
##rCyto.addIngredient( ingr2OT8 )
#
#
#
## 1TWT
#ingr1TWT = MultiSphereIngr( MSca*.000025, name='30S Ribosome', pdb='1TWT',
#                            sphereFile=wrkDir+'/recipes/cyto/1TWT_centered_18.sph',
#                            meshFile=wrkDir+'/recipes/cyto/1TWT_centered',
#                            color=antiquewhite, packingPriority= -.1)
##rCyto.addIngredient( ingr1TWT )
#
#
#
#ingr1TWV = MultiSphereIngr( MSca*.00001, name='50S RIBOSOME', pdb='1TWV',
#                            sphereFile=wrkDir+'/recipes/cyto/1TWV_centered_20.sph',
#                            meshFile=wrkDir+'/recipes/cyto/1TWV_centered',
#                            color=oldlace, packingPriority= -.3)
#
##rCyto.addIngredient( ingr1TWV )
#ingr3B63_14 = MultiSphereIngr( MSca*.00001, name='ActinFilament_14mer', pdb='3B63',
#                            sphereFile=wrkDir+'/recipes/cyto/3B63_centered_36.sph',
#                            meshFile=wrkDir+'/recipes/cyto/3B63_centered',
#                            color=darkseagreen, packingPriority= -8)
##rCyto.addIngredient( ingr3B63_14 )
#
#
#
#ingr3B63_7pre = MultiSphereIngr( MSca*.00001, name='ActinFilament_7merPre', pdb='3B63AtoG',
#                                sphereFile=wrkDir+'/recipes/cyto/3B63AtoG_centered_28.sph',
#                            meshFile=wrkDir+'/recipes/cyto/3B63AtoG_centered',
#                            color=darkseagreen, packingPriority= -8)
##rCyto.addIngredient( ingr3B63_7pre )
#
#
#
#ingr3B63_7 = MultiSphereIngr( MSca*.00003, name='ActinFilament_7mer', pdb='3B63AtoG',
#                            sphereFile=wrkDir+'/recipes/cyto/3B63AtoG_centered_28.sph',
#                            meshFile=wrkDir+'/recipes/cyto/3B63AtoG_centered',
#                            color=darkseagreen, packingPriority= -.2)
#
##rCyto.addIngredient( ingr3B63_7 )
#
#ingr3B63_3 = MultiSphereIngr( MSca*.0002, name='ActinFilament_3mer', pdb='3B63AtoC',
#                            sphereFile=wrkDir+'/recipes/cyto/3B63AtoC_centered_12.sph',
#                            meshFile=wrkDir+'/recipes/cyto/3B63AtoC_centered',
#                            color=darkseagreen, packingPriority= 1.8)
#
##rCyto.addIngredient( ingr3B63_3 )
#
#
#
#ingr3B63_2 = MultiSphereIngr( MSca*.0004, name='ActinFilament_2mer', pdb='3B63AtoB',
#                            sphereFile=wrkDir+'/recipes/cyto/3B63AtoB_centered_8.sph',
#                            meshFile=wrkDir+'/recipes/cyto/3B63AtoB_centered',
#                            color=darkseagreen, packingPriority= 1.4)
##rCyto.addIngredient( ingr3B63_2 )
#
#
#
#ingr3B63_1 = MultiSphereIngr( MSca*.002, name='ActinFilament_1mer', pdb='3B63A',
#                            sphereFile=wrkDir+'/recipes/cyto/3B63A_centered_8.sph',
#                            meshFile=wrkDir+'/recipes/cyto/3B63A_centered',
#                            color=darkseagreen, packingPriority= 1.3)
##rCyto.addIngredient( ingr3B63_1 )

# vesicle
#from DejaVu.IndexedPolygons import IndexedPolygonsFromFile
# create HistoVol

h1 = Environment()

#geomS = IndexedPolygonsFromFile(localwrkDir+'/cache_organelles/vesicle', 'vesicle')
#faces = geomS.getFaces()
#vertices = geomS.getVertices()
#vnormals = geomS.getVNormals()

#helper.triangulate(c4dorganlle)
#faces,vertices,vnormals = helper.DecomposeMesh(c4dorganlle,edit=False,
#                                               copy=False,tri=True,transform=False)

#o1 = Organelle("vesicle",vertices, faces, vnormals)
#h1.addOrganelle(o1)
#o1.setSurfaceRecipe(rSurf1)
#o1.setInnerRecipe(rMatrix1)
h1.setExteriorRecipe(rCyto)
#o1.overwriteSurfacePts=False
h1.overwriteSurfacePts=False

#define the viewer type dejavu,c4d,blender
h1.name="GenericCytoplasm"
afviewer = AFViewer(ViewerType=helper.host,helper=helper)#long ?

#make some option here
afviewer.doPoints = False
afviewer.doSpheres = False
afviewer.quality = 1 #lowest quality for sphere and cylinder
afviewer.visibleMesh = True #mesh default visibility

h1.setMinMaxProteinSize()
print 'Cyto', rCyto.getMinMaxProteinSize()
print 'Surf', rSurf1.getMinMaxProteinSize()
print 'Matrix', rMatrix1.getMinMaxProteinSize()
#print 'o1', o1.getMinMaxProteinSize()
print 'smallest', h1.smallestProteinSize
print 'largest', h1.largestProteinSize
h1.smallestProteinSize = 30

print 'smallest via Override', h1.smallestProteinSize
print 'largest via Override', h1.largestProteinSize

h1.boundingBox = [[0.,0.,0.],[300.,300.,300.]]
bbox=helper.getObject("histoVolBB")
if bbox is None : bbox = helper.box("histoVolBB",cornerPoints=[[-150.,-150.,-150.],[150.,150.,150.]])#cornerPoints=h1.boundingBox)

pad = 10.
afviewer.SetHistoVol(h1,pad,display=False)

h1.pickWeightedIngr = True ##do we sort the ingrediant or not see  getSortedActiveIngredients
h1.pickRandPt = True ##point pick randomly or one after the other?
h1.placeMethod = "jitter"#"spring" #"sphere"#"spring" #or "sphere"#
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
#h1.overwriteSurfacePts = False
#h1.ingrLookForNeighbours = False


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


h1.saveResult = True
#resultfilename = h1.resultfile =AutoFill.RECIPES["GenericCytoplasm"]["resultfile"]
resultfilename = h1.resultfile =localwrkDir+os.sep+"autoFillRecipeScripts"+os.sep+"GenericCytoplasm"+os.sep+"results"+os.sep+"GenericCytoplasmfillResult.afr"
afviewer.displayPreFill()
afviewer.printIngrediants()

#execfile(plgDir+'/extension/testAF/c_displayPreFill.py')

try :
    AFGui.Set("GenericCytoplasm",helper=afviewer.helper,afviewer=afviewer,histoVol=h1,bbox=bbox)
except:
    print ("no GUI")


