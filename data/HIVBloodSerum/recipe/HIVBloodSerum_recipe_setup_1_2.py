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

Name: 'figureHIV3.5mesh2'
@author: Graham Johnson and Ludovic Autin
"""

#execfile("/Users/ludo/DEV/AutoFill/HIV/figureHiv3.5mesh.py")
#execfile("/Users/ludo/DEV/Autofill_svn_test/autofill/trunk/AutoFillClean/autoFillRecipeScripts/HIV/figureHiv3.5mesh2.py")
#neeed PDB and not mesh....
import os
import numpy
from time import time
tStart = time()
import sys
#sys.path.append("/Users/grahamold/Dev/autoFillSVN/autofillSVNversions/trunk/")
#sys.path.append("/Users/ludo/DEV/Autofill_svn_test/autofill/trunk/")

#sys.path.append("/Users/ludo/DEV/af_svn/trunk/")
#sys.path.append("/Users/ludo/DEV/upy/trunk/")
#AUTOFILL
from AutoFill import Ingredient
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
from AutoFill.Ingredient import MultiCylindersIngr,GrowIngrediant,ActinIngrediant
from AutoFill.Organelle import Organelle
from AutoFill.Recipe import Recipe
from AutoFill.HistoVol import Environment
from AutoFill.autofill_viewer import AFViewer

#Directory
import upy
import AutoFill
wrkDir = AutoFill.__path__[0]
#from  Pmv import hostappInterface
#plgDir = hostappInterface.__path__[0]

#DEJAVU COLORS

from upy.colors import red, aliceblue, antiquewhite, aqua, \
    goldenrod, thistle, violet, skyblue, royalblue, plum, \
     pink, orchid, mediumvioletred, lightpink, \
     aquamarine, azure, beige, bisque, black, blanchedalmond, \
     mediumpurple, mediumorchid, hotpink, lightskyblue, \
     lightblue, lightcyan, dodgerblue, firebrick, \
     blue, blueviolet, brown, burlywood, cadetblue, slateblue, \
     darkslateblue, darkviolet, darkmagenta, powderblue, palevioletred, indigo, snow,\
     chartreuse, chocolate, coral, cornflowerblue, cornsilk, \
     crimson, cyan, darkblue, darkcyan, darkgoldenrod, \
     orange, purple, deeppink, lightcoral, \
     blue, cyan, mediumslateblue, steelblue, darkcyan, lime,\
     limegreen, darkorchid, tomato, khaki, gold, magenta, green, grey

helper = AutoFill.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper =helperClass()
    AutoFill.helper = helper

#if AutoFill.helper.host == 'c4d' or AutoFill.helper.host == 'maya' :
#    modeltype = "fbx"#"fbx"#"dae" #fbx
#else :
modelFormat = "fbx"
modeltype = "fbx"
hext= AutoFill.helper.hext
if AutoFill.helper.host == 'c4d':
    modelFormat = "dae"
    
modelFormatSpecial = modelFormat
if AutoFill.helper.host == 'c4d':
    modelFormatSpecial = "_c4d.dae"
 #mesfile?
#wrkdir could be html ?
#replace the workingdir by the url
httpdir = "http://autofill.googlecode.com/svn/data/"
httpwrkDirHIV=httpdir+"HIV/"
httpwrkDirBlood=httpdir+"BloodSerum/"
#httpwrkDir = "http://grahamj.com/autofill/autoFillData/HIV/HIV_0_0_1/"
##sphDir = httpwrkDir+"spheres/"
##meshDir = httpwrkDir+"geoms/"
wrkDirOrga = httpwrkDirHIV+"geoms/"

#if modelFormat = "dae":
#    meshDir = meshDir+"fbxIngredients/"
#else:
#    meshDir = meshDir+"daeIngredients/"

#will have sph and dejavu file
wrkDir3 = AutoFill.__path__[0]+os.sep+".."+os.sep+"AutoFill"+os.sep 
#wrkDir3 should point to original path since it is where the file are
localdir = AutoFill.__path__[0]+os.sep+"autoFillRecipeScripts"+os.sep+"HIVBloodSerum"+os.sep

MSca = 1.0
rSurf ={"Capside":Recipe(),"Nucleus":Recipe()}
#rSurfM ={"Nucleus":Recipe()}
rMatrix = {"Capside":Recipe(),"Nucleus":Recipe(),"ECM_box":Recipe()}
rCyto = Recipe()

VesK = 0.56
mode = "jitter" #"jitter" #"rigid-body" #or "rigid-body" except surface that are always jitter.
dosurf = 1
dolipid = False#True
dolipidCyl = False
domatrix = 1
doMatrix = 1
donucl = 1
useBigVol = 0
useMeshFile = False
rRNA = False

if useBigVol:
    hvb = "BighistoVolBB2"
else :
    hvb = "histoVolBB"
    
dict = {}

#Viral Enzymes pink:
#RT 1hys
#IN 1ex4
#PR 1hpv
#Structural Proteins lightblue
#MA 1hiw #inner surface -1> cylinder fat to do OR sphere fileradii=[[41.]],positions=[[[0,-5,0]]], positions2=[[[0,-60,0]]],
#CA 3h47 #capside part of one big ingredient nucleus, but need in the matrix 
# agglomerate in cyto->attractor of itself
#TODO....
#SU-TM 1g9m (SU, top) and 2ezo (TM, bottom)=> how to build it ?

#NC 1a1t #yellow orange
#Accessory Proteins green
#Vpu 1pi7 and 1vpu # how to build
#Vif 3dcg chain E
#Vpr 1esx -> multicylinder -> build it / extend clusterGui for cylinder ?
#P6 ?
#Nef 1avv and 1qa5
#Rev 1etf
#Tat 1biv and 1jfw.

#LipidSurface
dict["d1uph_MA"] = {"packingPriority":-3,"nMol":0,}#"color":lightblue}#29/row
dict["sutm"] = {"packingPriority":-1,"nMol":0}#,"color":plum}#29/row
dict["Nef"] = {"packingPriority":-1,"nMol":20}#,"color":limegreen}#5/row
if dolipidCyl:
    dict["LipidOut1"] = {"packingPriority":200.,"nMol":0,"color":grey}
    dict["LipidIn1"] = {"packingPriority":200.,"nMol":0,"color":grey}

if domatrix:
    dict["1hys"] = {"packingPriority":0.,"nMol":0,}#"color":mediumpurple}#14
    dict["1ex4"] = {"packingPriority":0.,"nMol":20,}#"color":mediumpurple}#162
    dict["1hpv"] = {"packingPriority":0.,"nMol":20,}#"color":mediumpurple }#8
    dict["1vpu"] = {"packingPriority":0.,"nMol":1,}#"color":limegreen}#5
    dict["3dcg"] = {"packingPriority":0.,"nMol":4,}#"color":limegreen}#2
    dict["Vif"] = {"packingPriority":0.,"nMol":20,}#"color":limegreen}#2
    dict["1esx"] = {"packingPriority":0.,"nMol":10,}#"color":limegreen}#22 flexible ?
    dict["3h47"] = {"packingPriority":0.,"nMol":0,}#"color":lightskyblue}#?
    dict["hex001"] = {"packingPriority":0.,"nMol":10}#,"color":lightskyblue}#?

#NucleoCapsid
dict["1a1t"] = {"packingPriority":0.,"nMol":0,"color":plum}#2
dict["Tat"] = {"packingPriority":0,"nMol":0,"color":limegreen}#5/row
dict["1etf"] = {"packingPriority":0.,"nMol":0,"color":limegreen}#8=

#NucleoSurface
dict["fpen"] = {"packingPriority":0.,"nMol":50}#,"color":plum}#2
dict["hex"] = {"packingPriority":0,"nMol":50}#,"color":lightskyblue}#5/row

#ExperimentalJunk
dict["ingrSynVes34"] = {"packingPriority":1.,"nMol":5}#,"color":deeppink }#8
dict["ingrSynVes34a"] = {"packingPriority":1.,"nMol":5}#,"color":blue }#8

#BloodSerum
#dict["dHIV"] =                  {"packingPriority":-200.,"nbMol":1,"color":None}#2
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
#packingmode = "gradient"#close"
grnameMA = "gradRadMA"
grnameSPIKE = "gradRadSPIKE"


#start with HIV recipe
httpwrkDir = httpwrkDirHIV
sphDir = httpwrkDir+"spheres/"
meshDir = httpwrkDir+"geoms/"

HIVmatrixConc = 1.0  

#==============================================================================
# BloodSerum recipe
#==============================================================================
httpwrkDir = httpwrkDirBlood
sphdir = sphDir = httpwrkDir+"spheres/"
meshDir = httpwrkDir+"geoms/"

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
    m,mesh = helper.Sphere("HIV",radius=725.,res=12)
    mHIV = MultiSphereIngr( BloodConc*800.0E-12,  725.,
                            name='HIVsphere', pdb=None,
#                            meshObject=m,
                            meshFile=httpwrkDirHIV+"geoms/HIV.fbx"#+modelFormat,
                            placeType=mode,  
                            #sphereFile=wrkDir+'/1e7i.sph',
                            **dict["dHIV"]
                            #packingMode='close'
                            )                               
    rMatrix["ECM_box"].addIngredient( mHIV )
    should be a previous ingrdient
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
        meshFile=meshDir+"/plv_Polio."+modelFormat,#+modeltype,#wrkDir+'/iTat',
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
        meshFile=meshDir+"/Fibrinogen."+modelFormat,
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
        meshFile=meshDir+"/FactorH1."+modelFormat,
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
        meshFile=meshDir+"/Heparin."+modelFormat,
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
        meshFile=meshDir+"/SerumAlbumin."+modelFormat,
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
        meshFile=meshDir+"/AntiTrypsin."+modelFormat,
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
        meshFile=meshDir+"/Transferrin."+modelFormat,
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
        meshFile=meshDir+"/iIgG_Antibody_1mer."+modelFormat,
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
        meshFile=meshDir+"/iIgA_Antibody_2mer."+modelFormat,
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
        meshFile=meshDir+"/iIgM_Antibody_5mer."+modelFormat,    
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
        meshFile=meshDir+"/Insulin."+modelFormat,   
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
        meshFile=meshDir+"/iLDL."+modelFormat,
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
        meshFile=meshDir+"/dLDL."+modelFormat,
        placeType=mode,
        sphereFile=sphdir+'/EMDB_5239.sph',  
        **dict["dLDL_EMDB_5239"]
        #packingMode='close'
        )
    rMatrix["ECM_box"].addIngredient( ildl2 )  
    
    #Replace these with low res versions by default:
#    idhdl = MultiSphereIngr( BloodConc*3.93618461E-6,  140., 
#        name='idHDLs', pdb=None,
#        meshFile=meshDir+"/idHDLs."+modelFormat,
#        placeType=mode,
#        sphereFile=sphdir+'/dHDLminusH_simulated.sph',  
#        **dict["ddHDL_Simulated"]
#        #packingMode='close'
#        )
#    rMatrix["ECM_box"].addIngredient( idhdl ) 
#    
#    ishdl = MultiSphereIngr( BloodConc*3.93618461E-6,  140., 
#        name='isHDLs', pdb=None,
#        meshFile=meshDir+"/isHDLs."+modelFormat,
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
        meshFile=meshDir+"/Hemoglobin."+modelFormat,
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
        meshFile=meshDir+"/Ceruloplasmin."+modelFormat,
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
        meshFile=meshDir+"/i1ysx_23kD."+modelFormat,
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
        meshFile=meshDir+"/i1smd_56kD."+modelFormat,
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
        meshFile=meshDir+"/i7aat_91kD."+modelFormat,
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
        meshFile=meshDir+"/i2tsc_63kD."+modelFormat,
        placeType=mode,  
        sphereFile=sphdir+'/2tsc_63kD.sph',
        **dict["d2tsc_63kD"]
        #packingMode='close'
        )                                
    rMatrix["ECM_box"].addIngredient( i2tsc )  
    

#    rMatrix["Nucleus"].addIngredient( dnaIngr )
    
#afviewer.displayIngrGrow(dnaIngr)

#afviewer.displayIngrResults(esurf3,doSphere=True,doMesh=False)
# create HistoVol
h1 = Environment()
h1.setExteriorRecipe(rCyto)

#display the organel, the box, and prepare the hierachy...

#===============================================================================
# Organelles Setup
#===============================================================================
#what abounucleo capside 
# vesicle
from DejaVu.IndexedPolygons import IndexedPolygonsFromFile
#    from DejaVu.IndexedPolygons import IndexedPolygonsFromFile
name=["Capside","Nucleus"]
if helper.host == "3dsmax":
    hext = "FBX"
for n in name :
    f=wrkDirOrga+'/'+n
    if n == "Nucleus" :#and helper.host == "c4d":
        o1 = Organelle(n,None, None, None,
               filename=f,object_name = "HIV_Nucleocapside",
#               object_filename = "http://grahamj.com/autofill/autoFillData/HIV/HIV_0_0_1/HIV_OrganelleGeoms/HIV_1_1_NucleocapsidHostMesh.c4d")
                object_filename = httpwrkDirHIV+"geoms/HIV_1_1_NucleocapsidHostMesh."+hext)
    else :
        o1 = Organelle(n,None, None, None,
               filename=f)
    
    o1.overwriteSurfacePts = True
    h1.addOrganelle(o1)
    if n in rSurf: 
        if rSurf[n].ingredients:
            r  = rSurf[n]
            o1.setSurfaceRecipe(r)
    if n in rMatrix:
        if rMatrix[n].ingredients:
            r = rMatrix[n]
            o1.setInnerRecipe(r)
            
h1.setExteriorRecipe(rMatrix["ECM_box"])

h1.setMinMaxProteinSize()
#print 'Cyto', rCyto.getMinMaxProteinSize()
#print 'Surf', rSurf1.getMinMaxProteinSize()
#print 'Matrix', rMatrix1.getMinMaxProteinSize()
#print 'o1', o1.getMinMaxProteinSize()
#print 'smallest', h1.smallestProteinSize
#print 'largest', h1.largestProteinSize
h1.smallestProteinSize = 50 #10 takes 23minutes   #15
#print 'smallest via Override', h1.smallestProteinSize
#print 'largest via Override', h1.largestProteinSize
#print o1.innerRecipe#

pad = 0.
h1.pickWeightedIngr = True ##do we sort the ingrediant or not see  getSortedActiveIngredients
h1.pickRandPt = True ##point pick randomly or one after the other?

h1.placeMethod = "jitter"#"spring" #"sphere"#"spring" #or "sphere"#
h1.overwritePlaceMethod = False #do we overwrtite all ingr place method
h1.simulationTimes = 30 #number of time setTime(i/fps) is called
h1.windowsSize = 50. #radius + histoVol.largestProteinSize + spacing + windowsSize
h1.windowsSize_overwrite = False
if h1.windowsSize_overwrite :
    h1.windowsSize = actine.influenceRad
    #size of the windows to check for neighbours

h1.runTimeDisplay = 0 #interactive display
#Specify here the option for the dynamic if method is spring

#h1.EnviroOnly = False
#h1.EnviroOnlyrtiment  =  -1
h1.overwriteSurfacePts = True #no need to discretisize the surface mesh as its already really denses
h1.ingrLookForNeighbours = False

h1.hackFreepts = True #SpeedUp?  Must Refer to getPointToDropHack from Fill4

h1._timer = False # Verbose for timing every function
h1._hackFreepts = False#True#bool(args[1])  # Strong hack that will never update the grids... good for very sparse fill!
h1._freePtsUpdateThrehod = 0.01#float(args[3])  # If my object covers more than 1% of the remaining freepoints, update the perIngredientAvailablePoints Lists

if grnameMA is not None :
#    h1.setGradient(name=grname1,mode="Y",direction=[0.5,0.5,0.0])
#    h1.setGradient(name=grname2,mode="X")
    h1.setGradient(name=grnameMA, mode="radial", weight_mode="linear", pick_mode="rnd", direction=[480.0, 480.0, 0.0], radius=1380.)
    h1.setGradient(name=grnameSPIKE, mode="radial", weight_mode="linear", pick_mode="rnd", direction=[-480.0, -480.0, 0.0], radius=1380.)


#def SetUPViewer():
#    import upy
#    helperClass = upy.getHelperClass()
ViewerType=AutoFill.helper.host
#    vi = None
#if ViewerType=='dejavu':
#    from DejaVu import Viewer
#    vi = Viewer()    
#    helper = helperClass(master=vi)#master=vi
#else :
#    helper = helperClass()
#    if ViewerType=='dejavu':
#    from DejaVu import Viewer
#    vi = Viewer()    
#        helper = helperClass()#self.GUI.VIEWER)#master=vi
#    else :
#        helper = helperClass()

afviewer = AFViewer(helper=AutoFill.helper,ViewerType=ViewerType)

#make some option here 
afviewer.doPoints = 1#True
afviewer.doSpheres = 0#1
afviewer.quality = 1 #lowest quality for sphere and cylinder
afviewer.visibleMesh = True #mesh default visibility 
h1.boundingBox = [[-1700.,-850.,-1050.],[1700,850.,1050.]]#1579.83 cm1700.248 cm3500 cm
h1.name="HIVBloodSerum"
afviewer.SetHistoVol(h1,pad,display=False)

afviewer.doSpheres = 0    
afviewer.displayPreFill()
#    afviewer.printIngrediants()
#    return afviewer

#h1.host=ViewerType
#afviewer = SetUPViewer()
h1.saveResult = False

resultfilename = h1.resultfile = httpdir+"HIVBloodSerum/HIVBloodSerum_1_0.apr"
#resultfilename = h1.resultfile =wrkDir3+"/HIV/hivRes1/fillResultHIV.af"
#resultfilename =h1.resultfile = wrkDir3+"/HIV/hivafresult_srf/fillResult.af"
#resultfilename =h1.resultfile =wrkDir3+"/HIV/grahamafrogra0.txt"
#print "resultfilename",resultfilename
bbox = None
#create the box
bbox = afviewer.helper.getObject(hvb)
if bbox is None : 
    h1.boundingBox = [[-1700.,-850.,-1050.],[1700,850.,1050.]]
    bbox = afviewer.helper.box(hvb,cornerPoints=h1.boundingBox)
#print "TEST",bbox
helper = afviewer.helper

noGUI = False
try :
    print ("try")
    AFGui.Set("HIVBloodSerum",helper=afviewer.helper,afviewer=afviewer,histoVol=h1,bbox=bbox)
except:
    print ("no GUI")
    noGUI = True


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
#   print 'time to fill', t2-t1
   return t

def DISPLAY():
   t2 = time()
   afviewer.displayFill()
#   afviewer.vi.toggleDisplay(afviewer.bsph,False)
#   print 'time to display', time()-t2
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
    h.buildGrid(boundingBox=bb,gridFileIn=gridFileIn,rebuild=True ,
                      gridFileOut=gridFileOut,previousFill=False)
#    h.buildGrid(gridFileIn=gridFileIn, 
#                  gridFileOut=gridFileOut)
    t2 = time()
    gridTime = t2-t1
    print ('time to Build Grid', gridTime)
    if fill :
        FILL(h,seed=seed,vTestid = vTestid,vAnalysis = vAnalysis)
    
#    afviewer.displayOrganellesPoints()
    #return
    #actine.updateFromBB(h.grid)

def FILL(h,seed=20,forceBuild=True,vTestid = 3,vAnalysis = 0):
    t1 = time()
    h.fill5(seedNum=seed,verbose=4,fbox=None, vTestid = vTestid,vAnalysis = vAnalysis)
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

def load(h,fname):
    if fname.find("http") != -1 or fname.find("ftp") != -1:
        #http://grahamj.com/autofill/autoFillData/HIV/HIVresult_2_afr.afr
        try :
            import urllib.request as urllib# , urllib.parse, urllib.error
        except :
            import urllib
        name =   fname.split("/")[-1]
        tmpFileName = wrkDir3+os.sep+"cache"+os.sep+name
        if not os.path.isfile(tmpFileName):
            urllib.urlretrieve(fname, tmpFileName)
        fname = tmpFileName
    h.resultfile =fname# wrkDir3+"/HIV/hivafresult_srf/fillResult.af"
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
#    print 'time to fill', t2-t1
    afviewer.displayFill()
#    print 'time to display', time()-t2
#    afviewer.vi.toggleDisplay(afviewer.bsph,False)
def doloop(n):
    # doLoop automatically produces result files, images, and documents from the recipe while adjusting parameters
    # To run doLoop, 1) in your host's python console type:
    # execfile(pathothis recipe) # for example, on my computer:
    #  execfile("/Users/grahamold/Dev/autoFillSVN/autofillSVNversions/trunk/AutoFillClean/autoFillRecipeScripts/HIV/figureHiv3.5mesh2.py")
    # 2) prepare your scene for the rendering->render output should be 640,480 but you can change the size in the script at then end.  Set up textures lights, and effects as you wish
    #    Results will appear in the result folder of your recipe path
    rangeseed=range(n)
    for i in rangeseed:
        basename = localdir+os.sep+"results"+os.sep+"results_seed_"+str(i)
        h1.saveResult = True
        resultfilename = h1.resultfile = basename  
        GRID(h1,seed=i, fill=1,vTestid = i,vAnalysis = 1)#build grid and fill
        #render/save scene
        helper.render(basename+".jpg",640,480)
        helper.write(basename+".c4d",[])
        #clean
        #afviewer.clearAll("Test_Spheres2D")
        afviewer.clearFill("HIVBloodSerum")
