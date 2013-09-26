# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:59:01 2013

@author: ludo
Hoe to run autopack from xml file
"""
# Show the effect of garbage collection
import sys
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
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14 Student_7B992864/plugins/ePMV/mgl64/MGLToolsPckgs")
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14 Student_7B992864/plugins/ePMV/mgl64/MGLToolsPckgs/PIL/")
import autopack
#wrkDir = AutoFill.__path__[0]
localdir = wrkDir = autopack.__path__[0]

from autopack.Environment import Environment
from autopack.Graphics import AutopackViewer as AFViewer
from autopack.Analysis import AnalyseAP
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
autopack.helper = helper
#filename = "/Users/ludo/Desktop/cell.xml"
#filename = wrkDir+os.sep+"autoFillRecipeScripts/2DsphereFill/Test_Spheres2D1.1.xml"
filename = "/Users/ludo/DEV/autopack_git/data/Mycoplasma/recipe/Mycoplasma1.3.xml"
#filename = "/Users/ludo/Desktop/cell_hack.xml"
fileName, fileExtension = os.path.splitext(filename)
n=os.path.basename(fileName)
h = Environment(name=n)
#h.helper = helper
recipe=n
h.load_XML(filename)
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
resultfilename = h.resultfile = wrkDir+os.sep+"autoFillRecipeScripts/Mycoplasma/results/Mycoplasma_1.2_mp.apr"
#h.exteriorRecipe.ingredients[0].uLength = 100.0
if ANALYSIS:
    analyse = AnalyseAP(env=h, viewer=afviewer, result_file=None)
    output=localdir+os.sep+"autoFillRecipeScripts/Cell/results/"
    analyse.g.Resolution = 1.0
    d=analyse.doloop(1,h.boundingBox,wrkDir,output,rdf=True,render=False,twod=TWOD)
else :
    gridfile = localdir+os.sep+"autoFillRecipeScripts/Mycoplasma/results/grid_store"
    h.placeMethod="RAPID"
    h.saveResult = True
    h.innerGridMethod = "jordan"#jordan pure python ? sdf ?
    h.boundingBox = [[-250.0, -6500.0/2.0, -250.0], [250.0, 6500.0/2.0, 250.0]]
#    h.buildGrid(boundingBox=h.boundingBox,gridFileIn=None,rebuild=True ,
#                          gridFileOut=gridfile,previousFill=False)
    h.buildGrid(boundingBox=h.boundingBox,gridFileIn=gridfile,rebuild=True ,
                          gridFileOut=None,previousFill=False)

    h.fill5(verbose = 3,usePP=True)
#    if not NOGUI :
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
#execfile("/Users/ludo/DEV/autopack_git/autopack/scripts/testAF_xml.py")       