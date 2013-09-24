# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:59:01 2013

@author: ludo
Hoe to run autopack from xml file
"""
# Show the effect of garbage collection
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
import AutoFill
#wrkDir = AutoFill.__path__[0]
localdir = wrkDir = AutoFill.__path__[0]

from AutoFill.HistoVol import Environment
from AutoFill.autofill_viewer import AFViewer
from AutoFill.analysis import AnalyseAP
TWOD = 1
NOGUI = 1
helper = AutoFill.helper
if helper is None and not NOGUI:
    import upy
    helperClass = upy.getHelperClass()
    helper =helperClass()
else :
    import upy
    helperClass = upy.getHelperClass()
    helper =helperClass(vi="nogui")
AutoFill.helper = helper
#filename = "/Users/ludo/Desktop/cell.xml"
#filename = "/Users/ludo/DEV/autofill_svn/trunk/AutoFillClean/autoFillRecipeScripts/2DsphereFill/Test_Spheres2D1.0.xml"
filename = "/Users/ludo/Desktop/cell_hack.xml"
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
resultfilename = h.resultfile = wrkDir+os.sep+"autoFillRecipeScripts/Cell/results/CellfillResult.afr"
h.exteriorRecipe.ingredients[0].uLength = 100.0
analyse = AnalyseAP(env=h, viewer=afviewer, result_file=None)
output=localdir+os.sep+"autoFillRecipeScripts/Cell/results/"
analyse.g.Resolution = 1.0
d=analyse.doloop(1,h.boundingBox,wrkDir,output,rdf=True,render=False,twod=TWOD)

#analyse.grid_pack(h.boundingBox,wrkDir,seed=0, fill=1,vTestid = 0,vAnalysis = 1)#build grid and fill
#for i in range(2):
#    print 'Collecting %d ...' % i
#    n = gc.collect()
#    print 'Unreachable objects:', n
#    print 'Remaining Garbage:', 
#    pprint.pprint(gc.garbage)
#    del gc.garbage[:]
#    print
#         