# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 15:02:57 2012

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

Name: 'SynapticVesicle_setup_recipe'
@author: Graham Johnson and Ludovic Autin based on file provided by Takamori et al, Cell 2006
"""

#execfile("/Library/MGLTools/mgltools_i86Darwin9_1.5.6/MGLToolsPckgs/Pmv/hostappInterface/extension/testAF/vesicleG10.py") # on Newton
#execfile("/Library/MGLTools/1.5.6.up/MGLToolsPckgs/AutoFill/scripts/vesicleG10.py")
import numpy
from time import time

import sys

#AUTOFILL
import AutoFill
wrkDir=localwrkDir = AutoFill.__path__[0]
httpwrkDir = "http://grahamj.com/autofill/autoFillData/SynapticVesicle/"
#AUTOFILL
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
from AutoFill.Ingredient import MultiCylindersIngr, GrowIngrediant
from AutoFill.Organelle import Organelle
from AutoFill.Recipe import Recipe
from AutoFill.HistoVol import Environment
from AutoFill.autofill_viewer import AFViewer


#DEJAVU COLORS
from upy.colors import red, aliceblue, antiquewhite, aqua, goldenrod, thistle, violet, skyblue, royalblue, plum, pink, orchid, mediumvioletred, lightpink, \
     aquamarine, azure, beige, bisque, black, blanchedalmond, mediumpurple, mediumorchid, hotpink, lightskyblue, lightblue, lightcyan, dodgerblue, firebrick, \
     blue, blueviolet, brown, burlywood, cadetblue, slateblue, darkslateblue, darkviolet, darkmagenta, powderblue, palevioletred, indigo,\
     chartreuse, chocolate, coral, cornflowerblue, cornsilk, \
     crimson, cyan, darkblue, darkcyan, darkgoldenrod, \
     orange, purple, deeppink, lightcoral, \
     blue, cyan, mediumslateblue, steelblue, darkcyan, lime,\
     limegreen, darkorchid, tomato, khaki, gold, magenta, green

helper = AutoFill.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper =helperClass()

modelFormat = 'dae'
modeltype = 'dae'
sphdir=httpwrkDir+"SynapticVesicleSpheres/"
meshdir=httpwrkDir+"SynapticVesicle_0_0_1_IngrGeom_MedRes/"#
#meshdir=httpwrkDir+"SynapticVesicleFaces/"
#define the viewer type dejavu,c4d,blender
dolipid = False
docyto = False
dosurf = True

MSca = 1.0
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

## cyl test 1
cyl1Ingr = MultiCylindersIngr(MSca*.04, color=coral, pdb='1CYL',
                             name='SingleCylinder', radii=[[23]],
                              positions=[[[-35,0,0]]], positions2=[[[45,0,0]]],
                              packingPriority=-1,
                              principalVector=(1,0,0)
                              )
#rSurf1.addIngredient(cyl1Ingr)

cylCoord1 = [ (-35,0,0), (45,0,0),  (45,0,0),    (85, 40, 0),]
cylCoord2 = [ (45,0,0),  (85,40,0), (105, -40,0), (135, 40, 30)]
cylRadii = [25, 20, 16, 16]

cyl4Ingr = MultiCylindersIngr(MSca*.04, color=aquamarine, pdb='1CYL',
                             name='4Cylinders', radii=[cylRadii],
                              positions=[cylCoord1], positions2=[cylCoord2],
                              packingPriority=-1,
                              principalVector=(1,0,0)
                              )

#rSurf1.addIngredient(cyl4Ingr)


#############  CELL paper vesicle recipe  ################

VesK = 0.56
if dosurf:
    # ves34
    ingrSynVes34 = MultiSphereIngr( VesK*.01480, color=mediumvioletred, pdb='Vesicle_3A_4BAlignedDimer2',
                                 name='Ves34MeshsParent',
                                 sphereFile=sphdir+'Vesicle_3A_4BAlignedDimer2_16.sph',
                                 meshFile=meshdir+'Vesicle_3A_4BAlignedDimer2.'+modelFormat,
                                 packingPriority=-10,
                                 #jitterMax=(1,1,0.2),
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes34)
    #is the principal actor depends on the host ?
    # ves34
    ingrSynVes5 = MultiSphereIngr( VesK*.02220, color=deeppink, pdb='Vesicle_3A_4BAlignedDimer2',
                                 name='Ves5MeshsParent',
                                 sphereFile=sphdir+'Vesicle_5AAlignedMonomer3_21.sph',
                                 meshFile=meshdir+'Vesicle_5AAlignedMonomer3.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes5)
    
    ingrSynVes8 = MultiSphereIngr( VesK*.02220, color=purple, pdb='Vesicle_3A_4BAlignedDimer2',
                                 name='Ves8MeshsParent',
                                 sphereFile=sphdir+'Vesicle_8Ato10CAlignedTrimer3_7.sph',
                                 meshFile=meshdir+'Vesicle_8Ato10CAlignedTrimer3.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes8)
    
    ingrSynVes17 = MultiSphereIngr( VesK*.02220, color=mediumslateblue, pdb='ves17',
                                 name='Ves17MeshsParent',
                                 sphereFile=sphdir+'Vesicle_17AAligned3_5.sph',
                                 meshFile=meshdir+'Vesicle_17AAligned3.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes17)
    
    ingrSynVes20 = MultiSphereIngr( VesK*.01480, color=mediumpurple, pdb='ves20',
                                 name='Ves20MeshsParent',
                                 sphereFile=sphdir+'Vesicle_20AAligned2_3.sph',
                                 meshFile=meshdir+'Vesicle_20AAligned2.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes20)
    
    ingrSynVes22 = MultiSphereIngr( VesK*.07401, color=slateblue, pdb='ves22',
                                 name='Ves22MeshsParent',
                                 sphereFile=sphdir+'Vesicle_22AAligned10_7.sph',
                                 meshFile=meshdir+'Vesicle_22AAligned10.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes22)
    
    ingrSynVes32 = MultiSphereIngr( VesK*.11102, color=violet, pdb='ves32',
                                 name='Ves32MeshsParent',
                                 sphereFile=sphdir+'Vesicle_32AAligned15_9.sph',
                                 meshFile=meshdir+'Vesicle_32AAligned15.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes32)
    
    ingrSynVes47 = MultiSphereIngr( VesK*.02220, color=plum, pdb='ves47',
                                 name='Ves47MeshsParent',
                                 sphereFile=sphdir+'Vesicle_47AAligned3_6.sph',
                                 meshFile=meshdir+'Vesicle_47AAligned3.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes47)
    
    ingrSynVes50 = MultiSphereIngr( VesK*.01480, color=mediumorchid, pdb='ves50',
                                 name='Ves50MeshsParent',
                                 sphereFile=sphdir+'Vesicle_50AAligned2_9.sph',
                                 meshFile=meshdir+'Vesicle_50AAligned2.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes50)
    
    ingrSynVes52 = MultiSphereIngr( VesK*.01480, color=magenta, pdb='ves52',
                                 name='Ves52MeshsParent',
                                 sphereFile=sphdir+'Vesicle_52AAligned2_4.sph',
                                 meshFile=meshdir+'Vesicle_52AAligned2.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes52)
    
    ingrSynVes54 = MultiSphereIngr( VesK*.01480, color=darkslateblue, pdb='ves54',
                                 name='Ves54MeshsParent',
                                 sphereFile=sphdir+'Vesicle_54AAligned2_11.sph',
                                 meshFile=meshdir+'Vesicle_54AAligned2.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes54)
    
    ingrSynVes56 = MultiSphereIngr( VesK*.05921, color=darkorchid, pdb='ves56',
                                 name='Ves56MeshsParent',
                                 sphereFile=sphdir+'Vesicle_56AAligned8_8.sph',
                                 meshFile=meshdir+'Vesicle_56AAligned8.'+modelFormat,
                                 packingPriority=-10,
                                 #jitterMax=(1,1,0.2),
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes56)
    
    ingrSynVes64 = MultiSphereIngr( VesK*.11842, color=darkmagenta, pdb='ves64',
                                 name='Ves64MeshsParent',
                                 sphereFile=sphdir+'Vesicle_64AAAligned16_9.sph',
                                 meshFile=meshdir+'Vesicle_64AAAligned16.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes64)
    
    ingrSynVes80 = MultiSphereIngr( VesK*.05921, color=darkviolet, pdb='ves80',
                                 name='Ves80MeshsParent',
                                 sphereFile=sphdir+'Vesicle_80AAligned8_9.sph',
                                 meshFile=meshdir+'Vesicle_80AAligned8.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes80)
    
    ingrSynVes88 = MultiSphereIngr( VesK*.08141, color=blueviolet, pdb='ves88',
                                 name='Ves88MeshsParent',
                                 sphereFile=sphdir+'Vesicle_88AAligned11_8.sph',
                                 meshFile=meshdir+'Vesicle_88AAligned11.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-100,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes88)
    
    ingrSynVes99 = MultiSphereIngr( VesK*.11842, color=orchid, pdb='ves99',
                                 name='Ves99MeshsParent',
                                 sphereFile=sphdir+'Vesicle_99AAligned16_6.sph',
                                 meshFile=meshdir+'Vesicle_99AAligned16.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes99)
    
    ingrSynVes115 = MultiSphereIngr( VesK*.09622, color=hotpink, pdb='ves115',
                                 name='Ves115MeshsParent',
                                 sphereFile=sphdir+'Vesicle_115AAligned13_7.sph',
                                 meshFile=meshdir+'Vesicle_115AAligned13.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes115)
    
    ingrSynVes128 = MultiSphereIngr( VesK*.00740, color=skyblue, pdb='ves128',
                                 name='Ves128MeshsParent',
                                 sphereFile=sphdir+'Vesicle_128to132hexamer1Aligned_35.sph',
                                 meshFile=meshdir+'Vesicle_128to132hexamer1Aligned.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-010,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes128)
    
    ingrSynVes133 = MultiSphereIngr( VesK*.01480, color=royalblue, pdb='ves133',
                                 name='Ves133MeshsParent',
                                 sphereFile=sphdir+'Vesicle_133A2Aligned_4.sph',
                                 meshFile=meshdir+'Vesicle_133A2Aligned.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes133)
    
    ingrSynVes135 = MultiSphereIngr( VesK*.1924, color=dodgerblue, pdb='ves135',
                                 name='Ves135MeshsParent',
                                 sphereFile=sphdir+'Vesicle_135A26monomerAligned_8.sph',
                                 meshFile=meshdir+'Vesicle_135A26monomerAligned.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes135)
    
    ingrSynVes161 = MultiSphereIngr( VesK*.00740, color=cornflowerblue, pdb='ves161',
                                 name='Ves161MeshsParent',
                                 sphereFile=sphdir+'Vesicle_161to166hexamer1Aligned_30.sph',
                                 meshFile=meshdir+'Vesicle_161to166hexamer1Aligned.'+modelFormat,
                                 packingPriority=-100,
                                 #jitterMax=(1,1,0.2),
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes161)
    
    ingrSynVes167 = MultiSphereIngr( VesK*.07401, color=lightskyblue, pdb='ves167',
                                 name='Ves167MeshsParent',
                                 sphereFile=sphdir+'Vesicle_167and174dimers10needHighPriorityAligned_20.sph',
                                 meshFile=meshdir+'Vesicle_167and174dimers10needHighPriorityAligned.'+modelFormat,
                                 packingPriority=-20,
                                 #jitterMax=(1,1,0.2),
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes167)
    
    ingrSynVes187 = MultiSphereIngr( VesK*.04441, color=blue, pdb='ves187',
                                 name='Ves187MeshsParent',
                                 sphereFile=sphdir+'Vesicle_187A6monomerAligned_8.sph',
                                 meshFile=meshdir+'Vesicle_187A6monomerAligned.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes187)
    
    ingrSynVes188 = MultiSphereIngr( VesK*.04441, color=firebrick, pdb='ves188',
                                 name='Ves188MeshsParent',
                                 sphereFile=sphdir+'Vesicle_188A6monomerAligned_11.sph',
                                 meshFile=meshdir+'Vesicle_188A6monomerAligned.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes188)
    
    ingrSynVes199 = MultiSphereIngr( VesK*.02220, color=palevioletred, pdb='ves199',
                                 name='Ves199MeshsParent',
                                 sphereFile=sphdir+'Vesicle_199A3monomerAligned_4.sph',
                                 meshFile=meshdir+'Vesicle_199A3monomerAligned.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes199)
    
    ingrSynVes203 = MultiSphereIngr( VesK*.01, color=red, pdb='ves203',
                                 name='Ves203MeshsParent',
                                 sphereFile=sphdir+'Vesicle_203to251ATPaseAligned_43.sph',
                                 meshFile=meshdir+'Vesicle_203to251ATPaseAligned.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-20,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes203)
    
    ingrSynVes252 = MultiSphereIngr( VesK*.07401, color=steelblue, pdb='ves252',
                                 name='Ves252MeshsParent',
                                 sphereFile=sphdir+'Vesicle_252_VGLUT1_10Aligned_8.sph',
                                 meshFile=meshdir+'Vesicle_252_VGLUT1_10Aligned.'+modelFormat,
                                 packingPriority=-10,
                                 #jitterMax=(1,1,0.2),
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes252)
    
    ingrSynVes262 = MultiSphereIngr( VesK*.01480, color=lightpink, pdb='ves262',
                                 name='Ves262MeshsParent',
                                 sphereFile=sphdir+'Vesicle_262A2monomerAligned_6.sph',
                                 meshFile=meshdir+'Vesicle_262A2monomerAligned.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes262)
    
    ingrSynVes264 = MultiSphereIngr( VesK*.00740, color=indigo, pdb='ves264',
                                 name='Ves264MeshsParent',
                                 sphereFile=sphdir+'Vesicle_264A1monomerAligned_9.sph',
                                 meshFile=meshdir+'Vesicle_264A1monomerAligned.'+modelFormat,
                                 #jitterMax=(1,1,0.2),
                                 packingPriority=-10,
                                 principalVector=(0,0,-1))
    rSurf1.addIngredient(ingrSynVes264)

############## cylinder Bilayer here##########
# Outer Leaflet Lipids as Cylinders
if dolipid:
    cyl26IngrO = MultiCylindersIngr(MSca*7.0, color=bisque, pdb='LipOut26',
                             name='LipOut26', radii=[[4.5]],
                              positions=[[[-.5,0,0]]], positions2=[[[25.5,0,0]]],
                              packingPriority=1,
                              jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                              principalVector=(1,0,0)
                              )
    rSurf1.addIngredient(cyl26IngrO)
    
    cyl22IngrO = MultiCylindersIngr(MSca*3.0, color=bisque, pdb='LipOut22',
                                 name='LipOut22', radii=[[4]],
                                  positions=[[[-.5,0,0]]], positions2=[[[22.5,0,0]]],
                                  packingPriority=1,
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(1,0,0)
                                  )
    rSurf1.addIngredient(cyl22IngrO)
    
    cyl31IngrO = MultiCylindersIngr(MSca*2.0, color=bisque, pdb='LipOut31',
                                 name='LipOut31', radii=[[4.5]],
                                  positions=[[[-.5,0,0]]], positions2=[[[30.5,0,0]]],
                                  packingPriority=1,
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(1,0,0)
                                  )
    rSurf1.addIngredient(cyl31IngrO)
    
    # Inner Leaflet Lipids as Cyninders
    cyl26IngrI = MultiCylindersIngr(MSca*2.0, color=burlywood, pdb='LipIn26',
                                 name='LipIn26', radii=[[4.5]],
                                  positions=[[[-.5,0,0]]], positions2=[[[25.5,0,0]]],
                                  packingPriority=1,
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(-1,0,0)
                                  )
    rSurf1.addIngredient(cyl26IngrI)
    
    cyl22IngrI = MultiCylindersIngr(MSca*5.0, color=burlywood, pdb='LipIn22',
                                 name='LipIn22', radii=[[4]],
                                  positions=[[[-.5,0,0]]], positions2=[[[22.5,0,0]]],
                                  packingPriority=1,
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(-1,0,0)
                                  )
    rSurf1.addIngredient(cyl22IngrI)
    
    cyl31IngrI = MultiCylindersIngr(MSca*1.0, color=burlywood, pdb='LipIn31',
                                 name='LipIn31', radii=[[4.5]],
                                  positions=[[[-.5,0,0]]], positions2=[[[30.5,0,0]]],
                                  packingPriority=1,
                                  jitterMax=(0.3,0.8,0.8), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(-1,0,0)
                                  )
    rSurf1.addIngredient(cyl31IngrI)
    
    # extra because HistoVol is broken?
    cyl3IngrT = MultiCylindersIngr(MSca*.003, color=burlywood, pdb='LipInPdb3',
                                 name='Lip3', radii=[[3]],
                                  positions=[[[-.5,0,0]]], positions2=[[[15,0,0]]],
                                  packingPriority=.01,
                                  jitterMax=(0.3,0.3,0.3), #  CRITICAL !!! IF jitter is greater than radius of object, e.g. 5x(1,1,1) the point may not be consumed!!!
                                  principalVector=(-1.,0,0)
                                  )
    rSurf1.addIngredient(cyl3IngrT)
 
##
## Cytoplasm:
# 1ABL
dnaIngr = GrowIngrediant(MSca*1., color=coral, pdb=None,
                              name='dnaGrow', radii=[[3.],],
                              positions=[[[0,0,0]]], positions2=[[[0.,50.,0]]],
                              packingPriority=-2,
                              biased=0.5,
                              principalVector=(0,1,0),
                              #marge = 180.0,
                              modelType="Cylinders",placeType="jitter",
                              nbJitter=1, jitterMax=(0,0,0),
                              length = 50000,closed = False,
                              nMol=1,
                              )
dnaIngr.cutoff_boundary=50.0
dnaIngr.cutoff_surface=50.0
dnaIngr.seedOnMinus = True
#rCyto.addIngredient( dnaIngr )
nmol = 10
if docyto:
    from upy.colors import brown, peru, saddlebrown, darkred, snow, seashell, darkgray, silver, linen, dimgrey, grey, oldlace, darkseagreen, ivory, gainsboro
    #HERE WE ARE ACTUALLY USINGTHE GENERIC CYTO INGREDIENT  
    #WHERE ARE THESES INGREDIENTS                          
    ingr1ABL = MultiSphereIngr( MSca*.009, name='1ABL', pdb='1ABL',
                                sphereFile=wrkDir+'/recipes/cyto/1ABL_centered_1.sph',#_1
                                meshFile=wrkDir+'/recipes/cyto/1ABL_centered',
                                color=snow, packingPriority=1.1)#, packingMode='close')
    #rCyto.addIngredient( ingr1ABL )
    
    ingr2ABL = MultiSphereIngr( MSca*.0045, name='2ABL', pdb='2ABL',
                                sphereFile=wrkDir+'/recipes/cyto/1ABL_centered_1.sph',#_1
                                meshFile=wrkDir+'/recipes/cyto/1ABL_centered',
                                color=gainsboro, packingPriority=1.1)#, packingMode='close')
    #rCyto.addIngredient( ingr2ABL )
    
#    ingr1IRK = MultiSphereIngr( MSca*.004, name='1IRK', pdb='1IRK',
#                                sphereFile=wrkDir+'/recipes/cyto/1IRK_centered_7.sph',
#                                meshFile=wrkDir+'/recipes/cyto/1IRK_centered',
#                                color=ivory, packingPriority=1.3,
#                                nMol=nmol)#, packingMode='close')
#    rCyto.addIngredient( ingr1IRK )#=>problem ?
#    
    ingr6TNA = MultiSphereIngr( MSca*.004, name='tRNA', pdb='6TNA',
                                sphereFile=wrkDir+'/recipes/cyto/6TNA_centered_5.sph',
                                meshFile=wrkDir+'/recipes/cyto/6TNA_centered',
                                color=darkgray, packingPriority=0.,
                                nMol=200)#, packingMode='close')
#    rCyto.addIngredient( ingr6TNA )
    
    GroelIngr1 = MultiSphereIngr( MSca*.00001, name='Groel', pdb='1AON',
                                  sphereFile=wrkDir+'/recipes/cyto/1AON_centered.sph',
                                  meshFile=wrkDir+'/recipes/cyto/1AON_centered',
                                  color=silver, packingPriority= -1)
#    rCyto.addIngredient( GroelIngr1 )
    
    # 2CPK
    ingr2CPK = MultiSphereIngr( MSca*.005, name='PHOSPHOTRANSFERASE', pdb='2CPK',
                                sphereFile=wrkDir+'/recipes/cyto/2CPK_centered_7.sph',
                                meshFile=wrkDir+'/recipes/cyto/2CPK_centered',nMol=150,
                                color=linen, packingPriority=0)#, packingMode='close')
#    rCyto.addIngredient( ingr2CPK )
#    
    # 1CZA
    ingr1CZA = MultiSphereIngr( MSca*.0003, name='TRANSFERASE', pdb='1CZA',
                                sphereFile=wrkDir+'/recipes/cyto/1CZA_centered_7.sph',
                                meshFile=wrkDir+'/recipes/cyto/1CZA_centered',nMol=nmol,
                                color=grey, packingPriority=1.0)
#    rCyto.addIngredient( ingr1CZA )
#    
#    # 2OT8
#    ingr2OT8 = MultiSphereIngr( MSca*.0002, name='TRANSPORT PROTEIN', pdb='2OT8',
#                                sphereFile=wrkDir+'/recipes/cyto/2OT8_centered_11.sph',
#                                meshFile=wrkDir+'/recipes/cyto/2OT8_centered',nMol=nmol,
#                                color=dimgrey, packingPriority=1.5)
#    rCyto.addIngredient( ingr2OT8 )
    
    # 1TWT
    ingr1TWT = MultiSphereIngr( MSca*.000025, name='30S Ribosome', pdb='1TWT',
                                sphereFile=wrkDir+'/recipes/cyto/1TWT_centered_18.sph',
                                meshFile=wrkDir+'/recipes/cyto/1TWT_centered',
                                color=antiquewhite, packingPriority= -3)
#    rCyto.addIngredient( ingr1TWT )
    
    ingr1TWV = MultiSphereIngr( MSca*.00001, name='50S RIBOSOME', pdb='1TWV',
                                sphereFile=wrkDir+'/recipes/cyto/1TWV_centered_20.sph',
                                meshFile=wrkDir+'/recipes/cyto/1TWV_centered',
                                color=oldlace, packingPriority= -4)
#    rCyto.addIngredient( ingr1TWV )
#    
#    ingr3B63_14 = MultiSphereIngr( MSca*.00001, name='ActinFilament_14mer', pdb='3B63',
#                                sphereFile=wrkDir+'/recipes/cyto/3B63_centered_36.sph',
#                                meshFile=wrkDir+'/recipes/cyto/3B63_centered',
#                                color=darkseagreen, packingPriority= -8)
#    #rCyto.addIngredient( ingr3B63_14 )
#    
#    ingr3B63_7pre = MultiSphereIngr( MSca*.00001, name='ActinFilament_7merPre', pdb='3B63AtoG',
#                                sphereFile=wrkDir+'/recipes/cyto/3B63AtoG_centered_28.sph',
#                                meshFile=wrkDir+'/recipes/cyto/3B63AtoG_centered',
#                                color=darkseagreen, packingPriority= -8)
#    #rCyto.addIngredient( ingr3B63_7pre )
#    
#    ingr3B63_7 = MultiSphereIngr( MSca*.00003, name='ActinFilament_7mer', pdb='3B63AtoG',
#                                sphereFile=wrkDir+'/recipes/cyto/3B63AtoG_centered_28.sph',
#                                meshFile=wrkDir+'/recipes/cyto/3B63AtoG_centered',
#                                color=darkseagreen, packingPriority= -.2)
#    #rCyto.addIngredient( ingr3B63_7 )
#    
#    ingr3B63_3 = MultiSphereIngr( MSca*.0002, name='ActinFilament_3mer', pdb='3B63AtoC',
#                                sphereFile=wrkDir+'/recipes/cyto/3B63AtoC_centered_12.sph',
#                                meshFile=wrkDir+'/recipes/cyto/3B63AtoC_centered',
#                                color=darkseagreen, packingPriority= 1.8)
#    #rCyto.addIngredient( ingr3B63_3 )
#    
#    ingr3B63_2 = MultiSphereIngr( MSca*.0004, name='ActinFilament_2mer', pdb='3B63AtoB',
#                                sphereFile=wrkDir+'/recipes/cyto/3B63AtoB_centered_8.sph',
#                                meshFile=wrkDir+'/recipes/cyto/3B63AtoB_centered',
#                                color=darkseagreen, packingPriority= 1.4)
#    #rCyto.addIngredient( ingr3B63_2 )
#    
#    ingr3B63_1 = MultiSphereIngr( MSca*.002, name='ActinFilament_1mer', pdb='3B63A',
#                                sphereFile=wrkDir+'/recipes/cyto/3B63A_centered_8.sph',
#                                meshFile=wrkDir+'/recipes/cyto/3B63A_centered',
#                                color=darkseagreen, packingPriority= 1.3)
#    #rCyto.addIngredient( ingr3B63_1 )



# vesicle
from DejaVu.IndexedPolygons import IndexedPolygonsFromFile

# create HistoVol
h1 = Environment()

#the organelle name is the geometry name
o1 = Organelle("vesicle_Mesh",None, None, None,
               filename="http://grahamj.com/autofill/autoFillData/SynapticVesicle/SynapticVesicleOrganelles/SynapticVesicleSurfaceMesh1."+modeltype)

h1.addOrganelle(o1)
o1.setSurfaceRecipe(rSurf1)
#o1.setInnerRecipe(rMatrix1)
if docyto:h1.setExteriorRecipe(rCyto)
o1.overwriteSurfacePts=True

#define the viewer type dejavu,c4d,blender
h1.name="SynapticVesicle" 
afviewer = AFViewer(ViewerType=helper.host,helper=helper)#long ?
#make some option here 
afviewer.doPoints = True
afviewer.doSpheres = True
afviewer.quality = 1 #lowest quality for sphere and cylinder
afviewer.visibleMesh = True #mesh default visibility 

## add padding
#bb = h1.boundingBox
#pad = 100.
#x,y,z = bb[0]
#bb[0] = [x-pad+pad, y-pad, z-pad+pad]
#x,y,z = bb[1]
#bb[1] = [x+pad-pad, y+pad, z+pad-pad]
#print 'Bounding box x with padding', h1.boundingBox

h1.setMinMaxProteinSize()
#print 'Cyto', rCyto.getMinMaxProteinSize()
#print 'Surf', rSurf1.getMinMaxProteinSize()
#print 'Matrix', rMatrix1.getMinMaxProteinSize()
#print 'o1', o1.getMinMaxProteinSize()
#print 'smallest', h1.smallestProteinSize
#print 'largest', h1.largestProteinSize
h1.smallestProteinSize = 10

#print 'smallest via Override', h1.smallestProteinSize
#print 'largest via Override', h1.largestProteinSize

#print o1.innerRecipe

pad = [100,100,100] #100.
h1.name = "SynapticVesicle"
afviewer.SetHistoVol(h1,pad,display=False)
h1.pickWeightedIngr = True ##do we sort the ingrediant or not see  getSortedActiveIngredients
h1.pickRandPt = True ##point pick randomly or one after the other?
h1.placeMethod = "jitter"#"spring"#"spring" #"sphere"#"spring" #or "sphere"#
h1.overwritePlaceMethod = 1#False #do we overwrtite all ingr place method
h1.simulationTimes = 30 #number of time setTime(i/fps) is called
h1.windowsSize = 500. #radius + histoVol.largestProteinSize + spacing + windowsSize
h1.windowsSize_overwrite = False
if h1.windowsSize_overwrite :
    h1.windowsSize = actine.influenceRad
    #size of the windows to check for neighbours
h1.runTimeDisplay = False #interactive display
#Specify here the option for the dynamic if method is spring

#h1.EnviroOnly = False
#h1.EnviroOnlyCompartiment  =  -1
#h1.overwriteSurfacePts = False
h1.ingrLookForNeighbours = False

h1.hackFreepts = True #SpeedUp?  Must Refer to getPointToDropHack from Fill4


h1._timer = False # Verbose for timing every function
h1._hackFreepts = False#True#bool(args[1])  # Strong hack that will never update the grids... good for very sparse fill!
h1._freePtsUpdateThrehod = 0.01#float(args[3])  # If my object covers more than 1% of the remaining freepoints, update the perIngredientAvailablePoints Lists


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
resultfilename = h1.resultfile = localwrkDir+os.sep+"autoFillRecipeScripts/SynapticVesicle/results/SynapticVesiclefillResult.afr"

afviewer.displayPreFill()
afviewer.printIngrediants()
#execfile(plgDir+'/extension/testAF/c_displayPreFill.py')
#from Pmv.hostappInterface.cinema4d_dev import helperC4D as helper

bbox=helper.getObject("histoVolBB")
try :
    AFGui.Set("SynapticVesicle",helper=afviewer.helper,afviewer=afviewer,histoVol=h1,bbox=bbox)
except:
    print ("no GUI")


# t=FILLT(h1) #follow this with 'DISPLAY()'
def FILLT(h):
   import thread
   doc =helper.getCurrentScene()
   #box = doc.get_selection()[0]
   box = doc.GetSelection()[0]
   bb=helper.getCornerPointCube(box)
   h.buildGrid(boundingBox=bb)
   #t1 = time()
   t=thread.start_new_thread(h.fill4,(14,))
   #h.fill3(seedNum=14)
   #t2 = time()
#   print 'time to fill', t2-t1
   return t
   
def DISPLAY():
   t2 = time()
   afviewer.displayFill()
   afviewer.vi.toggleDisplay(afviewer.bsph,False)
#   print 'time to display', time()-t2

def FILL(h):
    doc =helper.getCurrentScene()
    #box = doc.get_selection()[0]
    box = doc.GetSelection()[0]
    bb=helper.getCornerPointCube(box)
    h.buildGrid(boundingBox=bb)
    t1 = time()
    h.fill5(seedNum=14)
    t2 = time()
#    print 'time to fill', t2-t1
    afviewer.displayFill()
#    print 'time to display', time()-t2
    #afviewer.vi.toggleDisplay(afviewer.bsph,False)
    #execfile(plgDir+'/extension/testAF/c_displayFill.py')
    #afviewer.showHide(afviewer.undspMesh)

#FILL(h1)
