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

Name: 'BloodSerum_setup_recipe'
@author: Graham Johnson and Ludovic Autin
"""
VERSION = "1.0"

#AUTOFILL
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr, SingleCubeIngr
from AutoFill.Ingredient import MultiCylindersIngr,GrowIngrediant,ActinIngrediant
from AutoFill.Organelle import Organelle
from AutoFill.Recipe import Recipe
from AutoFill.HistoVol import Environment
from AutoFill.autofill_viewer import AFViewer

#Directory

import AutoFill
wrkDir = AutoFill.__path__[0]
#from  Pmv import hostappInterface
#plgDir = hostappInterface.__path__[0]

#DEJAVU COLORS

from upy.colors import red, aliceblue, antiquewhite, aqua, goldenrod, thistle, violet, skyblue, royalblue, plum, pink, orchid, mediumvioletred, lightpink, \
     aquamarine, azure, beige, bisque, black, blanchedalmond, mediumpurple, mediumorchid, hotpink, lightskyblue, lightblue, lightcyan, dodgerblue, firebrick, \
     blue, blueviolet, brown, burlywood, cadetblue, slateblue, darkslateblue, darkviolet, darkmagenta, powderblue, palevioletred, indigo, snow,\
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
    
#wrkDir = "/Users/ludo/Desktop/ResultsFigure8BPollard/ePMVmoleculesBlood"#AutoFill.__path__[0]
#wrkDir='/Users/grahamc4d/Desktop/AutoFillPaper/BloodModel/ePMVmoleculesBlood/'#graham dir
#wrkDir=AutoFill.RECIPES["BloodSerum"]["wrkdir"]
wrkDir = AutoFill.__path__[0]+os.sep+"autoFillRecipeScripts"+os.sep+"BloodSerum"

sphdir = "http://grahamj.com/autofill/autoFillData/BloodSerum/BloodSpheres"#wrkDir+"/BloodSpheres/"
#sphdir = wrkDir+"/OldSphereTreesIncorrectMinMaxRadiiForFasterFill"
httpwrkDir = "http://grahamj.com/autofill/autoFillData/BloodSerum/"
#sphdir = httpwrkDir+"http://www.grahamj.com/autofill/autoFillData/BloodSerum/BloodSpheresTooSmallForSpeed/"

modelFormat = "dae"

MSca = 1.0
rSurf ={"Ellipse":Recipe(),"Sphere":Recipe()}
rMatrix = {"ECM_box":Recipe()}#is it cyto ?
rECM = Recipe() #exterior recipe ie MAtrix

VesK = 1. #0.56
mode = "jitter" #"rigid-body" #or "rigid-body" except surface that are always jitter.
doECM = True
doMatrix = True

dict = {}
 
dict["dHIV"] =                  {"packingPriority":-200.,"nbMol":1,"color":None}#2
dict["d2plv_Polio"] =           {"packingPriority":-11.,"nbMol":1,"color":None}#2
dict["d3ghg_Fibrinogen"] =      {"packingPriority":-2,"nbMol":0,"color":None}
dict["d3gau_FactorH1"] =        {"packingPriority":-2,"nbMol":0,"color":None}
dict["d3irl_Heparin"] =         {"packingPriority":-2,"nbMol":2,"color":None}
dict["d1e7i_SerumAlbumin"] =    {"packingPriority":0.,"nbMol":0,"color":None}
dict["d1atu_AntiTrypsin"] =     {"packingPriority":0.,"nbMol":0,"color":None}
dict["d1lfg_Transferrin"] =     {"packingPriority":0.,"nbMol":0,"color":None}
dict["dIgG_Antibody_1mer"] =    {"packingPriority":2.,"nbMol":0,"color":None}
dict["dIgA_Antibody_2mer"] =    {"packingPriority":1.,"nbMol":0,"color":None}
dict["dIgM_Antibody_5mer"] =    {"packingPriority":-10.,"nbMol":0,"color":None}
dict["d2hui_Insulin"] =         {"packingPriority":0.,"nbMol":2,"color":None}
dict["dLDL_EMDB_5421"] =        {"packingPriority":-3.,"nbMol":0,"color":None}
dict["dLDL_EMDB_5239"] =        {"packingPriority":-3.,"nbMol":0,"color":None}
#dict["ddHDL_Simulated"] =       {"packingPriority":-3.,"nbMol":0,"color":None}
#dict["dsHDL_Simulated"] =       {"packingPriority":-3.,"nbMol":0,"color":None}
#Generic Proteins Suggested in David's Recipe:
dict["d1ysx_23kD"] =            {"packingPriority":0.,"nbMol":0,"color":None}
dict["d1smd_56kD"] =            {"packingPriority":0.,"nbMol":0,"color":None}
dict["d7aat_91kD"] =            {"packingPriority":0.,"nbMol":0,"color":None}
dict["d2tsc_63kD"] =            {"packingPriority":0.,"nbMol":0,"color":None}
dict["d2hhb_hemoglobin"] =      {"packingPriority":0.,"nbMol":0,"color":None}
dict["d1kcw_ceruloplasmin"] =   {"packingPriority":0.,"nbMol":0,"color":None}

#Buffer:
dict["d1e7i_SerumAlbumin2"] =    {"packingPriority":0.,"nbMol":1,"color":None}

#Retired:
#dict["dEM_Fibers"] =            {"packingPriority":-2,"nbMol":0,"color":None}
dict["iLDL2"] = {"packingPriority":-3.,"nbMol":2,"color":None}#2 LDL_EMDB_5421
#dict["iHDL"] = {"packingPriority":-2.,"nbMol":0,"color":None}#2
#dict["i1ITG"] = {"packingPriority":0.,"nbMol":0,"color":None}#22
#dict["CoarseMS_7aat"] = {"packingPriority":0.,"nbMol":0,"color":None}#8
#dict["CoarseMS_2tsc"] = {"packingPriority":0.,"nbMol":0,"color":None}#8
#dict["CoarseMS_1smd"] = {"packingPriority":0.,"nbMol":0,"color":None}#14
#dict["CoarseMS_2hiv"] = {"packingPriority":0.,"nbMol":0,"color":None}#1 


#===============================================================================
# fiber
#===============================================================================
#actually not the fiber mesh
cyl,snakem = helper.Cylinder('snakeCyl',radius=25.,
                               length=150.0,axis=[1.0,0.0,0.0])
snakeIngr = GrowIngrediant(.08*.4*1.18085538E-5, color=coral, pdb=None,
      name='snake', radii=[[25.],],
      positions=[[[0,0,0]]], positions2=[[[0., 150., 0.]]], #Ludo had 370 for y
      packingPriority=-10,
      biased=0.5,
      principalVector=(0,1,0),
      marge = 90.0,
      meshObject=cyl,
      modelType="Cylinders",placeType="jitter",
      nbJitter=1, jitterMax=(0,0,0),
      length = 8000,closed = False,
      #walkingMode = "lattice",
      nbMol=1,#2-3
      )
snakeIngr.cutoff_boundary=25.0
snakeIngr.cutoff_surface=25.0
snakeIngr.seedOnMinus = False#True
snakeIngr.marge = 40.#25.
snakeIngr.constraintMarge = True
rMatrix["ECM_box"].addIngredient(snakeIngr)
#afviewer.displayIngrGrow(snakeIngr)
#===============================================================================
# blood recipe
#0.00000393618461*David's 75nmCube number gives the molarity
#===============================================================================
BloodConc = 1.0
if doMatrix :
    #ho1 = helper.getObject("3ghg_Fibrinogen")
    m,mesh = helper.Sphere("HIV",radius=725.,res=12)
    mHIV = SingleSphereIngr( BloodConc*800.0E-12,  725.,
                            name='HIVsphere', pdb=None,
                            meshObject=m,
                            placeType=mode,  
                            #sphereFile=wrkDir+'/1e7i.sph',
                            **dict["dHIV"]
                            #packingMode='close'
                            )                               
    rMatrix["ECM_box"].addIngredient( mHIV )
    
#    m,mesh = helper.box("testB",center=[0.,0.,0.],size=[500,1000,2000])
#    ingr = SingleCubeIngr( BloodConc*800.0E-12,  radii = [[500,1000,2000],],name="test", pdb=None,
#                              positions=[[[0,0,0],[0,0,0],[0,0,0],]],
#                                meshObject=m,placeType=mode,
#                                rotAxis = (0.,0.,0.),
#                                  useRotAxis = 1,
#                                      ) 
#    rMatrix["ECM_box"].addIngredient( ingr )
#    hoo2=helper.getObject("2plv_Polio")
    #Polio Virus 2plv = 1/75nm cube therefor molarity = 0.00000393618461
    #reducing to 1/3 of David's molarity to make more realistic
    #Appears that viral load should be closer to 40,000/mL = 6.43E-17
    #I'll crank to 750,000/mL for high viral Load by multiplying by 18.75
    m2plv = MultiSphereIngr( BloodConc*18.75*6.43E-17,  50., 
        name='plv_Polio', 
        meshFile=httpwrkDir+"BloodIngredientGeom/plv_Polio."+modelFormat,#+modeltype,#wrkDir+'/iTat',
        placeType=mode, 
        pdb="2plv",
        sphereFile=sphdir+'/2plv.sph',
        **dict["d2plv_Polio"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( m2plv )   

#    hoo3=helper.getObject("3ghg_Fibrinogen")
    #Fibrinogen 3ghg = ~3/75nm cube therefor molarity = 0.00001180855383
    #MSMS volume = 266168.36
    m3ghg = MultiSphereIngr( BloodConc*1.18085538E-5,  
        name='Fibrinogen', 
        sphereFile=sphdir+'/3ghg.sph',
#                                meshObject=hoo3,
        meshFile=httpwrkDir+"BloodIngredientGeom/Fibrinogen."+modelFormat,
        placeType=mode,  
        pdb="3ghg",
        **dict["d3ghg_Fibrinogen"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( m3ghg )   
    
#    hoo3=helper.getObject("3gau_FactorH1")
    #Fibrinogen 3ghg = ~3/75nm cube therefor molarity = 0.00001180855383
    #MSMS volume = 266168.36
    m3gau = MultiSphereIngr( BloodConc*3.23E-6,  
        name='FactorH1', 
        sphereFile=sphdir+'/3gau_FactorH3.sph',
        meshFile=httpwrkDir+"BloodIngredientGeom/FactorH1."+modelFormat,
#                                meshObject=hoo3,
        placeType=mode,
        pdb="3gau",
        **dict["d3gau_FactorH1"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( m3gau )  

#    hoo3=helper.getObject("3irl_Heparin")
    #Fibrinogen 3ghg = ~3/75nm cube therefor molarity = 0.00001180855383
    #MSMS volume = 266168.36
    m3irl = MultiSphereIngr( BloodConc*3.23E-20,  
        name='Heparin',
        sphereFile=sphdir+'/3irl_heparin2.sph',
#                                meshObject=hoo3,
        pdb="3irl",
        meshFile=httpwrkDir+"BloodIngredientGeom/Heparin."+modelFormat,
        placeType=mode,  
        **dict["d3irl_Heparin"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( m3irl ) 
    
#    hoo1=helper.getObject("1e7i_SerumAlbumin")
    #Serum albumin 1e7i = 162/75nm cube therefor molarity = 0.0006376
    #m1e7i = SingleSphereIngr( BloodConc*0.0006376,  50.,
    #MSMS volume = 79320.58A^3
    m1e7i = MultiSphereIngr( BloodConc*0.0006376,  50.,  
        name='SerumAlbumin', 
#                                meshObject=hoo1,
        meshFile=httpwrkDir+"BloodIngredientGeom/SerumAlbumin."+modelFormat,
        placeType=mode,  
        pdb="1e7i",
        sphereFile=sphdir+'/1e7i.sph',
        **dict["d1e7i_SerumAlbumin"]
        #packingMode='close'
        )                               
    rMatrix["ECM_box"].addIngredient( m1e7i )        

#    hoo0=helper.getObject("1atu_AntiTrypsin")
    #Antitrypsin 1atu = 14/75nm cube therefor molarity = 0.00005510658454
    c1atu = MultiSphereIngr( BloodConc*0.00005510,
        name='AntiTrypsin', 
#                                meshObject=hoo0,
        meshFile=httpwrkDir+"BloodIngredientGeom/AntiTrypsin."+modelFormat,
        placeType=mode,  
        pdb="1atu",
        sphereFile=sphdir+'/1atu.sph',
        **dict["d1atu_AntiTrypsin"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( c1atu )

#    hoo0=helper.getObject("1lfg_Transferrin")
    #Transferrin 1lfg = 8/75nm cube therefor molarity = 0.00003148947688
    c1lfg = MultiSphereIngr( BloodConc*3.14894769E-5,
        name='Transferrin', 
        pdb="1lfg",
#                                meshObject=hoo0,
        meshFile=httpwrkDir+"BloodIngredientGeom/Transferrin."+modelFormat,
        placeType=mode,  
        sphereFile=sphdir+'/1lfg.sph',
        **dict["d1lfg_Transferrin"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( c1lfg )

#    hoo3=helper.getObject("IgG_Antibody_1mer")
    #IgG Antibodies antibody.pdb (file from David G.) = 22/75nm cube therefor molarity = 0.00008659606142
    mantB = MultiSphereIngr( BloodConc*0.000086596,  50., 
        name='IgG_Antibody_1mer', pdb=None, 
#                                meshObject=hoo3,
        placeType=mode,  
        meshFile=httpwrkDir+"BloodIngredientGeom/iIgG_Antibody_1mer."+modelFormat,
        sphereFile=sphdir+'/antB.sph',
        **dict["dIgG_Antibody_1mer"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( mantB )   

#    hoo3=helper.getObject("IgA_Antibody_2mer")
    #IgA Antibodies antibody.pdb*2 (file from David G.) = 2/75nm cube therefor molarity = 0.00000787236922
    igad = MultiSphereIngr( BloodConc*7.87236921E-6,
        #radii=[[40,85,40]],
        #positions=[[[-100,0,0],[0,0,0],[+100,0,0]]],
        name='iIgA_Antibody_2mer', pdb=None,
        sphereFile=sphdir+'/IgA_Antibody_2mer.sph',
#                                meshObject=hoo3,
        meshFile=httpwrkDir+"BloodIngredientGeom/iIgA_Antibody_2mer."+modelFormat,
        placeType=mode,  
#                                sphereFile=wrkDir+'/1atu.sph',
        **dict["dIgA_Antibody_2mer"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( igad )
    
#    hoo3=helper.getObject("IgM_Antibody_5mer")
    #IgA Antibodies antibody.pdb*2 (file from David G.) = 0.5/75nm cube therefor molarity = 0.000001968092305
    igmh = MultiSphereIngr( BloodConc*1.9680923E-6,  
    #igmh = MultiCylindersIngr( BloodConc*1.9680923E-6,  
        name='iIgM_Antibody_5mer', pdb=None,
        #radii=[[175]],
        #positions=[[[-50,0,0]]], positions2=[[[50,0,0]]],
#                                meshObject=hoo3,
        #principalVector=(1,0,0),
        meshFile=httpwrkDir+"BloodIngredientGeom/iIgM_Antibody_5mer."+modelFormat,    
        sphereFile=sphdir+'/IgM_Antibody_5mer.sph',
        placeType=mode,
        **dict["dIgM_Antibody_5mer"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( igmh ) 
    
#    hoo3=helper.getObject("2hui_Insulin")
    #Insulin 2hui.pdb NMR ~ >800pMol/L after eating to ~57 to 79pM/L between meals extremely low picoMolar conc
    ihui = MultiSphereIngr( BloodConc*800.0E-12,
        sphereFile=sphdir+'/2hui.sph',
        name='Insulin', pdb="2hui",
#        meshObject=hoo3,
        meshFile=httpwrkDir+"BloodIngredientGeom/Insulin."+modelFormat,   
        placeType=mode,  
        **dict["d2hui_Insulin"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( ihui )
    
#    hoo3=helper.getObject("LDL_EMDB_5421")
    #2.5mmol/L
    ildl1 = MultiSphereIngr( BloodConc*1.4E-6,  140., 
        name='iLDL', 
#                                meshObject=hoo3,
        meshFile=httpwrkDir+"BloodIngredientGeom/iLDL."+modelFormat,
        pdb = "EMDB_5421",
        placeType=mode, 
        sphereFile=sphdir+'/EMDB_5421.sph',
        **dict["dLDL_EMDB_5421"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( ildl1 )  
    
#    hoo3=helper.getObject("LDL_EMDB_5239")
    ildl2 = MultiSphereIngr( BloodConc*1.4E-6,  140., 
        name='dLDL', 
#        meshObject=hoo3,
        pdb="EMDB_5239",
        meshFile=httpwrkDir+"BloodIngredientGeom/dLDL."+modelFormat,
        placeType=mode,
        sphereFile=sphdir+'/EMDB_5239.sph',  
        **dict["dLDL_EMDB_5239"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( ildl2 )  
    
    #Replace these with low res versions by default:
#    idhdl = MultiSphereIngr( BloodConc*3.93618461E-6,  140., 
#        name='idHDLs', pdb=None,
#        meshFile=httpwrkDir+"BloodIngredientGeom/idHDLs."+modelFormat,
#        placeType=mode,
#        sphereFile=sphdir+'/dHDLminusH_simulated.sph',  
#        **dict["ddHDL_Simulated"]
#        #packingMode='close'
#        )
#    rMatrix["ECM_box"].addIngredient( idhdl ) 
#    
#    ishdl = MultiSphereIngr( BloodConc*3.93618461E-6,  140., 
#        name='isHDLs', pdb=None,
#        meshFile=httpwrkDir+"BloodIngredientGeom/isHDLs."+modelFormat,
#        placeType=mode,
#        sphereFile=sphdir+'/sHDLminusH_simulated.sph',  
#        **dict["dsHDL_Simulated"]
#        #packingMode='close'
#                            )
#    rMatrix["ECM_box"].addIngredient( ishdl ) 
    
#    hoo3=helper.getObject("2hhb_hemoglobin")
    #Hemoglobin 2hhb, extremely low conc of 175 to 600 nM
    i2hhb = MultiSphereIngr( BloodConc*300.0E-9,  140., 
        name='Hemoglobin', pdb="2hhb",
#        meshObject=hoo3,
        meshFile=httpwrkDir+"BloodIngredientGeom/Hemoglobin."+modelFormat,
        placeType=mode,
        sphereFile=sphdir+'/2hhb_hemoglobin.sph',  
        **dict["d2hhb_hemoglobin"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( i2hhb )   
      
#    hoo3=helper.getObject("1kcw_ceruloplasmin")
    #Hemoglobin 2hhb, extremely low conc of 175 to 600 nM
    i1kcw = MultiSphereIngr( BloodConc*2.5E-6,  140., 
        name='Ceruloplasmin', pdb="1kcw",
#        meshObject=hoo3,
        meshFile=httpwrkDir+"BloodIngredientGeom/Ceruloplasmin."+modelFormat,
        placeType=mode,
        sphereFile=sphdir+'/1kcw_ceruloplasmin.sph',  
        **dict["d1kcw_ceruloplasmin"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( i1kcw )
    
#==============================================================================
# #Generic Proteins from David's Blood Recipe
#==============================================================================

#    hoo1=helper.getObject("1ysx_23kD")
    #Serum albumin NMR fragment for generic 20kD = 2/75nm cube therefor molarity = 3.93618461E-6*2
    i1ysx = MultiSphereIngr( BloodConc*2*3.93618461E-6,  50.,  
        name='i1ysx_23kD', pdb="1ysx",#is it SerumAlbumin ?
#                                meshObject=hoo1,
        meshFile=httpwrkDir+"BloodIngredientGeom/i1ysx_23kD."+modelFormat,
        placeType=mode,  
        sphereFile=sphdir+'/1ysx_23kD.sph',
        **dict["d1ysx_23kD"]
        #packingMode='close'
        )                               
    rMatrix["ECM_box"].addIngredient( i1ysx )  
    
#    hoo1=helper.getObject("1smd_56kD")
    #Human Salivary Amylase for generic 40kD = 14/75nm cube therefor molarity = 3.93618461E-6*14
    i1smd = MultiSphereIngr( BloodConc*14*3.93618461E-6,  50.,  
        name='i1smd_56kD', pdb="1smd",
#        meshObject=hoo1,
        meshFile=httpwrkDir+"BloodIngredientGeom/i1smd_56kD."+modelFormat,
        placeType=mode,
        sphereFile=sphdir+'/1smd_56kD.sph',
        **dict["d1smd_56kD"]
        #packingMode='close'
        )                                
    rMatrix["ECM_box"].addIngredient( i1smd )  
    
#    hoo1=helper.getObject("7aat_91kD")
    #91kD AminoTransferase for generic 80kD = 8/75nm cube therefor molarity = 3.93618461E-6*8
    i7aat = MultiSphereIngr( BloodConc*8*3.93618461E-6,  50.,  
        name='i7aat_91kD', pdb="7aat",
#        meshObject=hoo1,
        meshFile=httpwrkDir+"BloodIngredientGeom/i7aat_91kD."+modelFormat,
        placeType=mode,  
        sphereFile=sphdir+'/7aat_91kD.sph',
        **dict["d7aat_91kD"]
        #packingMode='close'
        )                                
    rMatrix["ECM_box"].addIngredient( i7aat )  
    
#    hoo1=helper.getObject("2tsc_63kD")
    #63kD AminoTransferase for generic 60kD = 8/75nm cube therefor molarity = 3.93618461E-6*8
    i2tsc = MultiSphereIngr( BloodConc*8*3.93618461E-6,  50.,  
        name='i2tsc_63kD', pdb="2tsc",
#        meshObject=hoo1,
        meshFile=httpwrkDir+"BloodIngredientGeom/i2tsc_63kD."+modelFormat,
        placeType=mode,  
        sphereFile=sphdir+'/2tsc_63kD.sph',
        **dict["d2tsc_63kD"]
        #packingMode='close'
        )                                
    rMatrix["ECM_box"].addIngredient( i2tsc )  
    


#Making buffer around ECM with fake organelle for bounding box??
#if doECM :    
##    hoo1=helper.getObject("1e7i_SerumAlbumin2")
#    #Serum albumin 1e7i = 162/75nm cube therefor molarity = 0.0006376
#    #m1e7i = SingleSphereIngr( BloodConc*0.0006376,  50.,
#    m1e7i2 = MultiSphereIngr(0.01* BloodConc*0.0006376,  50.,  
#        name='SerumAlbumin2', pdb="1e7i",
##        meshObject=hoo1,
#        meshFile=httpwrkDir+"BloodIngredientGeom/SerumAlbumin2."+modelFormat,
#        placeType=mode,  
#        sphereFile=sphdir+'/1e7i.sph',
#        **dict["d1e7i_SerumAlbumin2"]
#        #packingMode='close'
#        ) 
#    rECM.addIngredient( m1e7i2 )
#    rMatrix1.addIngredient( m1e7i2 ) 



#afviewer.displayIngrResults(esurf3,doSphere=True,doMesh=False)
# vesicle
from DejaVu.IndexedPolygons import IndexedPolygonsFromFile

# create HistoVol
h1 = Environment()
h1.name="BloodSerum"
h1.version = VERSION
#display the organel, the box, and prepare the hierachy...

#define the viewer type dejavu,c4d,blender
afviewer = AFViewer(ViewerType=helper.host,helper=helper)#long ?


#===============================================================================
# Organelles Setup
#===============================================================================
#ok test with no organelle but the fill_box
print (wrkDir+"/Geometries/ECM_box."+modelFormat) 
#the organelle name is the geometry name
organelle = False
if organelle: 
    h1.setExteriorRecipe(rECM)
    o1 = Organelle("ECM_box",None, None, None,
                   filename=httpwrkDir+"BloodIngredientGeom/ECM_box."+modelFormat)
    h1.addOrganelle(o1)
    o1.setInnerRecipe(rMatrix["ECM_box"])
else :
    h1.setExteriorRecipe(rMatrix["ECM_box"])
#    fbbox = afviewer.helper.getObject("fillBB")
#    if fbbox is None : fbbox = afviewer.helper.box("fillBB",cornerPoints=[[-250.,-500.,-1000.],[250,500.,1000.]])#cornerPoints=h1.boundingBox)


#make some option here
afviewer.doPoints = 0
afviewer.doSpheres = 0
afviewer.quality = 1 #lowest quality for sphere and cylinder
afviewer.visibleMesh = 1 #mesh default visibility 

h1.setMinMaxProteinSize()
print 'ECM', rECM.getMinMaxProteinSize()
#print 'Surf', rSurf1.getMinMaxProteinSize()
#print 'Matrix', rMatrix1.getMinMaxProteinSize()
#print 'o1', o1.getMinMaxProteinSize()
print 'smallest', h1.smallestProteinSize
print 'largest', h1.largestProteinSize
h1.smallestProteinSize = 25 #20 takes ~20 sec #10 takes 23minutes   #15
print 'smallest via Override', h1.smallestProteinSize
print 'largest via Override', h1.largestProteinSize
#print o1.innerRecipe#

pad = 50.
h1.boundingBox = [[-300.,-550.,-1050.],[300,550.,1050.]]

afviewer.SetHistoVol(h1,pad,display=False)
h1.name="BloodSerum"
#h1.host='c4d'

h1.pickWeightedIngr = True ##do we sort the ingrediant or not see  getSortedActiveIngredients
h1.pickRandPt = True ##point pick randomly or one after the other?

h1.placeMethod = "jitter"#"spring" #"sphere"#"spring" #or "sphere"#
h1.overwritePlaceMethod = False #do we overwrtite all ingr place method
h1.simulationTimes = 30 #number of time setTime(i/fps) is called
h1.windowsSize = 1. #radius + histoVol.largestProteinSize + spacing + windowsSize
h1.windowsSize_overwrite = False
if h1.windowsSize_overwrite :
    h1.windowsSize = actine.influenceRad
    #size of the windows to check for neighbours

h1.runTimeDisplay = False #interactive display
#Specify here the option for the dynamic if method is spring

h1.hackFreepts = True #SpeedUp?  Must Refer to getPointToDropHack from Fill4

#h1.EnviroOnly = False
#h1.EnviroOnlyCompartiment  =  -1
h1.overwriteSurfacePts = True #no need to discretisize the surface mesh as its already really denses
h1.ingrLookForNeighbours = False


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
    
afviewer.displayPreFill()
#afviewer.printIngrediants()

h1.saveResult = True
#resultfilename = h1.resultfile = wrkDir+os.sep+"autoFillRecipeScripts"+os.sep+"2DsphereFill"+os.sep+"results"+os.sep+"SpherefillResult.afr"
#resultfilename = h1.resultfile = wrkDir+os.sep+"results/BloodSerumfillResult.afr"

bbox = afviewer.helper.getObject("histoVolBB")
if bbox is None : bbox = afviewer.helper.box("histoVolBB",cornerPoints=h1.boundingBox)#cornerPoints=h1.boundingBox)


#hivsphere = afviewer.helper.getObject("HIV")
afviewer.helper.reParent("HIV","AutoFillHider")
afviewer.helper.reParent('snakeCyl',"AutoFillHider")

try :
    AFGui.Set("BloodSerum",helper=afviewer.helper,afviewer=afviewer,histoVol=h1,bbox=bbox)
except:
    print ("no GUI")


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


    
def FILL(h,seed=3,forceBuild=True):
    t1 = time()
    h.fill5(seedNum=seed,verbose=0)
    t2 = time()
    print 'time to run Fill5', t2-t1
    afviewer.displayFill()
    print 'time to display Fill5', time()-t2
    afviewer.vi.toggleDisplay(afviewer.bsph,False)
    #execfile(plgDir+'/extension/testAF/c_displayFill.py')
    #afviewer.showHide(afviewer.undspMesh)
    afviewer.displayIngrGrow(snakeIngr)

def GRID(h,forceBuild=True,fill=True):
    t1 = time()
    doc =helper.getCurrentScene()
    #box = doc.get_selection()[0]
    box = helper.getObject("histoVolBB")#doc.GetSelection()[0]
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
#    afviewer.displayOrganellesPoints()
    #return
    #actine.updateFromBB(h.grid)
    if fill :
        FILL(h)
    print 'time to Build Grid', gridTime
    
def SecondFill(h):
    h.setMinMaxProteinSize()
#    pgrid = h.grid
    doc =helper.getCurrentScene()
    #box = doc.get_selection()[0]
    box = doc.GetSelection()[0]
    bb=helper.getCornerPointCube(box)
    if doPrev :
        afviewer.appendIngrInstance(previngr,sel = previngrInstance,bb=bb)
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
    afviewer.doSpheres = True
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
    print 'time to fill Load', t2-t1
    afviewer.displayFill()
    print 'time to display Load', time()-t2
    afviewer.vi.toggleDisplay(afviewer.bsph,False)
