# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:59:01 2013

@author: ludo
Hoe to run autopack from xml file
"""
# Show the effect of garbage collection
import sys
import numpy as np
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14 Student_7B992864/plugins/ePMV/mgl64/MGLToolsPckgs/")
#should have dejavu...
import gc
import pprint
for i in range(2):
    print 'Collecting %d ...' % i
    n = gc.collect()
    print 'Unreachable objects:', n
    print 'Remaining Garbage:', 
    pprint.pprint(gc.garbage)
    del gc.garbage[:]
    print
    
import os
import sys
#MAC
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14 Student_7B992864/plugins/ePMV/mgl64/MGLToolsPckgs")
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14 Student_7B992864/plugins/ePMV/mgl64/MGLToolsPckgs/PIL/")
#windows
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs\lib-tk")

import autopack
#wrkDir = AutoFill.__path__[0]
localdir = wrkDir = autopack.__path__[0]

from autopack.Environment import Environment
from autopack.Graphics import AutopackViewer as AFViewer
from autopack.Analysis import AnalyseAP
TWOD = 0
NOGUI = 1
ANALYSIS = 0
usePP=False
helper = autopack.helper
if helper is None and not NOGUI:
    import upy
    helperClass = upy.getHelperClass()
    helper =helperClass()
else :
    import upy
    helperClass = upy.getHelperClass()
    helper = helperClass(vi="nogui")
autopack.helper = helper
#filename = "/Users/ludo/Desktop/cell.xml"
#filename = wrkDir+os.sep+"autoFillRecipeScripts/2DsphereFill/Test_Spheres2D1.1.xml"
#filename = wrkDir+os.sep+"autoFillRecipeScripts/2DsphereFill/Test_Spheres2Dgradients1.0.xml"
filename = "/Users/ludo/DEV/autopack_git/data/Mycoplasma/recipe/Mycoplasma1.3.xml"
#filename = "/Users/ludo/Desktop/cell_hack.xml"
#filename = "/Users/ludo/DEV/autopack_git/autoPACK_database_1.0.0/recipes/NM_Analysis_FigureC1.3.xml"
filename = "D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_1.0.json"

fileName, fileExtension = os.path.splitext(filename)
n=os.path.basename(fileName)
h = Environment(name=n)
#h.helper = helper
recipe=n
h.loadRecipe(filename)

afviewer=None
if not NOGUI :
    print h,helper
    setattr(h,"helper",helper)
    afviewer = AFViewer(ViewerType=h.helper.host,helper=h.helper)
    afviewer.SetHistoVol(h,20.0,display=False)
    h.host=h.helper.host
    afviewer.displayPreFill()
h.saveResult = False
#resultfilename = h1.resultfile = wrkDir+os.sep+"autoFillRecipeScripts"+os.sep+"2DsphereFill"+os.sep+"results"+os.sep+"SpherefillResult.afr"
#resultfilename = h.resultfile = wrkDir+os.sep+"autoFillRecipeScripts/2DsphereFill/results/2DsphereFill_1.1.apr"
#resultfilename = h.resultfile = wrkDir+os.sep+"autoFillRecipeScripts/Mycoplasma/results/MycoplasmaPackResult_3"
resultfilename = h.resultfile = "D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_Results_july"
proxy_radius = 11.85
#h.smallestProteinSize=15
#h.exteriorRecipe.ingredients[0].uLength = 100.0
#overwrite the jiterMax usin jitterMax = [d[k]["rad"]/(15.*1.1547),d[k]["rad"]/(15.*1.1547),0.0]
#loopThroughIngr
def setJitter(ingr):
    pos = np.array(ingr.positions[0])
    if len(ingr.positions)==2 :
        pos = np.array(ingr.positions[1])
    R = np.linalg.norm(pos,axis=1).max()+proxy_radius
    #compare with vertices radius ?
    if  abs(ingr.encapsulatingRadius - R) > 0.1 :
        ingr.encapsulatingRadius = R 
    if ( sum(ingr.jitterMax) == 1.0) :
        ingr.jitterMax =[1.0,1.0,0.0]
    ingr.rejectionThreshold = 300
    #check the radius ?

#    if ingr.packingMode != "gradient":
#        ingr.molarity = 0.0
#        ingr.nbMol = 0
#        print ingr.name
    
def setCompartment(ingr):
#    ingr.checkCompartment=True
#    ingr.compareCompartment=True#slow down a little.
#    ingr.nbMol*=2
    ingr.rejectionThreshold=100#[1,1,0]#
#    ingr.jitterMax =[ingr.encapsulatingRadius/(25.*1.1547),ingr.encapsulatingRadius/(25.*1.1547),0.0]
#    ingr.cutoff_boundary=ingr.encapsulatingRadius
#    ingr.cutoff_boundary=500+ingr.encapsulatingRadius
    ingr.nbJitter = 6
#   
    
def filterOutFiber(env):
    names =["mpn208","mpn191","mpn529", "DNA","mRNA","peptides"]
    for n in names :
        ingr = env.getIngrFromName(n)
        ingr.molarity = 0.0
        ingr.nbMol = 0
        ingr.overwrite_nbMol_value = 0
    
#h.loopThroughIngr(setCompartment)
#setJitter
#raw_input()
if ANALYSIS:
    h.placeMethod="RAPID"
    h.encapsulatingGrid=0
    autopack.testPeriodicity = False
    analyse = AnalyseAP(env=h, viewer=afviewer, result_file=None)
    output="D:\\Data\\cellPACK_data\\Mycoplasma\\"
    analyse.g.Resolution = 1.0
    h.boundingBox=numpy.array(h.boundingBox)
    fbox_bb=numpy.array(h.boundingBox)
#    h.boundingBox[0]-=numpy.array([500.0,500.0,0.0])    
#    h.boundingBox[1]+=numpy.array([500.0,500.0,0.0])
    d=analyse.doloop(5,h.boundingBox,wrkDir,output,rdf=True,
                     render=False,twod=TWOD,use_file=True)#,fbox_bb=fbox_bb)
#    if not NOGUI :
#        afviewer.displayFill() 
else :
    autopack.testPeriodicity = False
    gridfile = "D:\\Data\\cellPACK_data\\Mycoplasma\\grid_store"
    h.placeMethod="RAPID"
    h.saveResult = True
    resultfilename = h.resultfile = "D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_Results_july"
    h.innerGridMethod = "trimesh" 
    h.smallestProteinSize = 80.0
    h.freePtsUpdateThrehod=0.0
    #remove HUE,RIBO,DNAPOLY,PEPTIDE,RNA,DNA,Lipo
    filterOutFiber(h)        
    h.loopThroughIngr(setJitter)
    h.buildGrid(boundingBox=h.boundingBox,gridFileIn=gridfile,rebuild=True ,gridFileOut=gridfile,previousFill=False)
    h.fill5(verbose = 0,usePP=usePP)
    h.collectResultPerIngredient()
    h.saveRecipe("D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_Results_H_1.json",useXref=False,mixed=True, kwds=["source","name"],result=True, transpose = True,  grid=False,packing_options=False,indent=False,quaternion=True) 
    h.saveRecipe("D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_Results_H_2.json",useXref=False,mixed=True, kwds=["source","name"],result=True, transpose = False, grid=False,packing_options=False,indent=False,quaternion=True) 

    from autopack.IOutils import serializedRecipe, saveResultBinary, toBinary
    djson, all_pos, all_rot = serializedRecipe(h,False,True,True,True)#transpose, use_quaternion, result=False, lefthand=False
    with open("D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_Results_serialized.json","w") as f:
        f.write(djson)
    
    from autopack.Serializable import sCompartment,sIngredientGroup,sIngredient,sIngredientFiber
    sCompartment.static_id = 0;sIngredientFiber.static_id = 0;sIngredient.static_id = [0, 0, 0];sIngredientGroup.static_id= 0;djson, all_pos1, all_rot1 = serializedRecipe(h,True,True,True,True);toBinary(all_pos1, all_rot1,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_2.0_serializedTR_L.bin")     
    sCompartment.static_id = 0;sIngredientFiber.static_id = 0;sIngredient.static_id = [0, 0, 0];sIngredientGroup.static_id= 0;djson, all_pos2, all_rot2 = serializedRecipe(h,False,True,True,True);toBinary(all_pos2, all_rot2,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_2.0_serialized_L.bin")     
    sCompartment.static_id = 0;sIngredientFiber.static_id = 0;sIngredient.static_id = [0, 0, 0];sIngredientGroup.static_id= 0;djson, all_pos3, all_rot3 = serializedRecipe(h,False,True,True,False);toBinary(all_pos3, all_rot3,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_2.0_serialized.bin")     
    sCompartment.static_id = 0;sIngredientFiber.static_id = 0;sIngredient.static_id = [0, 0, 0];sIngredientGroup.static_id= 0;djson, all_pos4, all_rot4 = serializedRecipe(h,True,True,True,False);toBinary(all_pos4, all_rot4,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_2.0_serializedTR.bin")     
    


##    if not NOGUI :
#        afviewer.displayFill()  

#analyse.grid_pack(h.boundingBox,wrkDir,seed=0, fill=1,vTestid = 0,vAnalysis = 1)#build grid and fill
#for i in range(2):
#    print 'Collecting %d ...' % i
#    n = gc.collect()
#    print 'Unreachable objects:', n
#    print 'Remaining Garbage:', 
#    pprint.pprint(gc.garbage)
#    del gc.garbage[:]
#    print
#execfile("/Users/ludo/DEV/autoPACK_github/autopack/scripts/testAF_xml.py")       
