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

Name: 'GenericCytoplasmRecipeFile1'
@author: Graham Johnson and Ludovic Autin
"""

from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
from AutoFill.Ingredient import MultiCylindersIngr, GrowIngrediant
#DEJAVU COLORS
from upy.colors import red, aliceblue, antiquewhite, aqua, goldenrod, thistle, violet, skyblue, royalblue, plum, pink, orchid, mediumvioletred, lightpink, \
     aquamarine, azure, beige, bisque, black, blanchedalmond, mediumpurple, mediumorchid, hotpink, lightskyblue, lightblue, lightcyan, dodgerblue, firebrick, \
     blue, blueviolet, brown, burlywood, cadetblue, slateblue, darkslateblue, darkviolet, darkmagenta, powderblue, palevioletred, indigo,\
     chartreuse, chocolate, coral, cornflowerblue, cornsilk, \
     crimson, cyan, darkblue, darkcyan, darkgoldenrod, \
     orange, purple, deeppink, lightcoral, \
     blue, cyan, mediumslateblue, steelblue, darkcyan, lime,\
     limegreen, darkorchid, tomato, khaki, gold, magenta, green,\
     snow

modelFormat = "dae"
mode  = "jitter"
MSca=1.0

#before execfile define tempRecipe
#hoo3=helper.getObject("1PFK_4mer1.c4d")
pr = MultiSphereIngr( MSca*.00006,
                     name='tetramer1PFK', pdb="1pfk",
                     sphereFile=wrkDir2 + '/1PFK_4mer1.sph',
#                     meshObject=hoo3,
                     meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/tetramer1PFK."+modelFormat,#+modeltype,#wrkDir+'/iTat',
                     placeType=mode,
                     packingPriority=0,
                     #nMol= 5,
                     color = red
                     #packingMode='close'
                     ) #original radius is 3.61
tempRecipe.addIngredient( pr )

#hoo3=helper.getObject("7TIM_1mer1.c4d")
pr2 = MultiSphereIngr( MSca*.00006,
                      name='dimer7TIM', pdb="7tim",
                      sphereFile=wrkDir2 + '/7TIM_1mer1.sph',
                      meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/dimer7TIM."+modelFormat,
#                      meshObject=hoo3,
                      placeType=mode,
                      packingPriority=0,
                      #nMol= 5,
                      color = red
                      #packingMode='close'
                      ) #original radius is 3.61
tempRecipe.addIngredient( pr2 )

#hoo3=helper.getObject("1A49_4mer1.c4d")
pr3 = MultiSphereIngr( MSca*.00006,
                      name='tetramer1A49', pdb="1a49",
                      sphereFile=wrkDir2 + '/1A49_4mer1.sph',
                      meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/tetramer1A49."+modelFormat,
#                      meshObject=hoo3,
                      placeType=mode,
                      packingPriority=0,
                      #nMol= 5,
                      color = red
                      #packingMode='close'
                      ) #original radius is 3.61
tempRecipe.addIngredient( pr3 )

#hoo3=helper.getObject("1GAX_1mer1.c4d")
pr4 = MultiSphereIngr( MSca*.00006,
                      name='onemer1GAX', pdb="1GAX",
                      sphereFile=wrkDir2 + '/1GAX_1mer1.sph',
                      meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/onemer1GAX."+modelFormat,
#                      meshObject=hoo3,
                      placeType=mode,
                      packingPriority=0,
                      #nMol= 5,
                      color = red
                      #packingMode='close'
                      ) #original radius is 3.61
tempRecipe.addIngredient( pr4 )

#hoo3=helper.getObject("1GAX_1mer1trnaOnly.c4d")
pr5 = MultiSphereIngr( MSca*.00006,
                      name='onemer1GAX_rna', pdb="1GAX",
                      sphereFile=wrkDir2 + '/1GAX_1mer1trnaOnly.sph',
                      meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/onemer1GAX_rna."+modelFormat,
#                      meshObject=hoo3,
                      placeType=mode,
                      packingPriority=0,
                      #nMol= 5,
                      color = red
                      #packingMode='close'
                      ) #original radius is 3.61
tempRecipe.addIngredient( pr5 )

#hoo3=helper.getObject("1GAX_1mer1protOnly.c4d")
pr6 = MultiSphereIngr( MSca*.00006,
                      name='onemer1GAX_prot', pdb="1GAX",
                      sphereFile=wrkDir2 + '/1GAX_1mer1protOnly.sph',
                      meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/onemer1GAX_prot."+modelFormat,
#                      meshObject=hoo3,
                      placeType=mode,
                      packingPriority=0,
                      #nMol= 5,
                      color = red
                      #packingMode='close'
                      ) #original radius is 3.61
tempRecipe.addIngredient( pr6 )

#hoo3=helper.getObject("1TWT_7mer1.c4d")
m1TWT_7mer1 = MultiSphereIngr( MSca*.00006,
                              name='hexamer1TWT', pdb="1TWT",
                              sphereFile=wrkDir2 + '/1TWT_7mer1.sph',
                              meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/hexamer1TWT."+modelFormat,
#                              meshObject=hoo3,
                              placeType=mode,
                              packingPriority=0,
                              #nMol= 5,
                              color = red
                              #packingMode='close'
                              ) #original radius is 3.61
tempRecipe.addIngredient( m1TWT_7mer1 )

#hoo3=helper.getObject("1TTT_1mer1protOnly.c4d")
m1TTT_1mer1protOnly = MultiSphereIngr( MSca*.00006,
                                      name='onemer1TTT_prot', pdb="1TTT",
                                      sphereFile=wrkDir2 + '/1TTT_1mer1protOnly.sph',
                                      meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/onemer1TTT_prot."+modelFormat,
#                                      meshObject=hoo3,
                                      placeType=mode,
                                      packingPriority=0,
                                      #nMol= 5,
                                      color = red
                                      #packingMode='close'
                                      ) #original radius is 3.61
tempRecipe.addIngredient( m1TTT_1mer1protOnly)

#hoo3=helper.getObject("1TTT_1mer1trnaOnly.c4d")
m1TTT_1mer1trnaOnly = MultiSphereIngr( MSca*.00006,
                                      name='onemer1TTT_rna', pdb="1TTT",
                                      sphereFile=wrkDir2 + '/1TTT_1mer1trnaOnly.sph',
                                      meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/onemer1TTT_rna."+modelFormat,
#                                      meshObject=hoo3,
                                      placeType=mode,
                                      packingPriority=0,
                                      #nMol= 5,
                                      color = red
                                      #packingMode='close'
                                      ) #original radius is 3.61
tempRecipe.addIngredient( m1TTT_1mer1trnaOnly)

#hoo3=helper.getObject("1TTT_1mer1.c4d")
m1TTT_1mer1 = MultiSphereIngr( MSca*.00006,
                              name='onemer1TTT', pdb="1TTT",
                              sphereFile=wrkDir2 + '/1TTT_1mer1.sph',
                              meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/onemer1TTT."+modelFormat,
#                              meshObject=hoo3,
                              placeType=mode,
                              packingPriority=0,
                              #nMol= 5,
                              color = red
                              #packingMode='close'
                              ) #original radius is 3.61
tempRecipe.addIngredient( m1TTT_1mer1)

#hoo3=helper.getObject("1EIY_6mer1.c4d")
m1EIY_6mer1 = MultiSphereIngr( MSca*.00006,
                              name='hexamer1EIY', pdb="1EIY",
                              sphereFile=wrkDir2 + '/1EIY_6mer1.sph',
                              meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/hexamer1EIY."+modelFormat,
#                              meshObject=hoo3,
                              placeType=mode,
                              packingPriority=0,
                              #nMol= 5,
                              color = red
                              #packingMode='close'
                              ) #original radius is 3.61
tempRecipe.addIngredient( m1EIY_6mer1)

#hoo3=helper.getObject("3GPD_4mer1.c4d")
m3GPD_4mer1 = MultiSphereIngr( MSca*.00006,
                              name='tetramer3GPD', pdb="3GPD",
                              sphereFile=wrkDir2 + '/3GPD_4mer1.sph',
                              meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/tetramer3GPD."+modelFormat,
#                              meshObject=hoo3,
                              placeType=mode,
                              packingPriority=0,
                              #nMol= 5,
                              color = red
                              #packingMode='close'
                              ) #original radius is 3.61
tempRecipe.addIngredient( m3GPD_4mer1)

#hoo3=helper.getObject("1CZA_1mer1.c4d")
m1CZA_1mer1 = MultiSphereIngr( MSca*.00006,
                              name='monomer1CZA', pdb="1CZA",
                              sphereFile=wrkDir2 + '/1CZA_1mer1.sph',
                              meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/monomer1CZA."+modelFormat,
#                              meshObject=hoo3,
                              placeType=mode,
                              packingPriority=0,
                              #nMol= 5,
                              color = red
                              #packingMode='close'
                              ) #original radius is 3.61
tempRecipe.addIngredient( m1CZA_1mer1)

#hoo3=helper.getObject("3PGK_1mer1.c4d")
m3PGK_1mer1 = MultiSphereIngr( MSca*.00006,
                              name='monomer3PGK', pdb="3PGK_1mer1",
                              sphereFile=wrkDir2 + '/3PGK_1mer1.sph',
                              meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/monomer3PGK."+modelFormat,
#                              meshObject=hoo3,
                              placeType=mode,
                              packingPriority=0,
                              #nMol= 5,
                              color = red
                              #packingMode='close'
                              ) #original radius is 3.61
tempRecipe.addIngredient( m3PGK_1mer1)

#hoo3=helper.getObject("2XSM_16mer.c4d")
m2XSM_16mer = MultiSphereIngr( MSca*.00006,
                              name='ico2XSM', pdb="2XSM",
                              sphereFile=wrkDir2 + '/2XSM_16mer.sph',
                              meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/ico2XSM."+modelFormat,
#                              meshObject=hoo3,
                              placeType=mode,
                              packingPriority=0,
                              #nMol= 5,
                              color = red
                              #packingMode='close'
                              ) #original radius is 3.61
tempRecipe.addIngredient( m2XSM_16mer)

#hoo3=helper.getObject("1ASY_4mer1.c4d")
m1ASY_4mer1 = MultiSphereIngr( MSca*.00006,
                              name='tetramer1ASY', pdb="1ASY_4mer1",
                              sphereFile=wrkDir2 + '/1ASY_4mer1.sph',
                              meshFile=httpwrkDir+"Cytoplasm_IngrGeom_Med/tetramer1ASY."+modelFormat,
#                              meshObject=hoo3,
                              placeType=mode,
                              packingPriority=0,
                              #nMol= 5,
                              color = red
                              #packingMode='close'
                              ) #original radius is 3.61
tempRecipe.addIngredient( m1ASY_4mer1)


####### END New generic Cytoplasm