# -*- coding: utf-8 -*-
"""
Created on Wed May  6 10:39:08 2015

@author: ludo
"""
#need panda and bhtree
#
import sys
import os
import json

#panda3d is there :
#sys.path.append("/usr/lib/python2.7/dist-packages/")
#LINUX
#sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs")
#Windows
sys.path.append("C:\Program Files\MAXON\CINEMA 4D R16 Demo\plugins\ePMV\mgl64\MGLToolsPckgs")
sys.path.append("C:\Program Files\MAXON\CINEMA 4D R16 Demo\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")
sys.path.append("C:\Program Files\MAXON\CINEMA 4D R16 Demo\plugins\ePMV\mgl64\MGLToolsPckgs\lib-tk")
#run with C:\>C:\Python26\python.exe "C:\Program Files\MAXON\CINEMA 4D R16 Demo\plugins\ePMV\mgl64\MGLToolsPckgs\autopack\scripts\run_autoPACK.py"  "C:\Users\ludovic\Downloads\Mycoplasma1.5_1.json"

try :
    from collections import OrderedDict
except :
    from ordereddict import OrderedDict

from autopack.Environment import Environment
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.8g')

import autopack
import upy
helperClass = upy.getHelperClass()
helper = helperClass(vi="nogui")
print helper
autopack.helper = helper
#autopack.fixpath = True
autopack.forceFetch = False

#example callback to run on all ingredient
def setJitter(ingr):
#    ingr.jitterMax = 50#[ingr.encapsulatingRadius,ingr.encapsulatingRadius,0.0]
    ingr.nbJitter = 50
    ingr.rejectionThreshold=500
    
def includeIngredient_cb(env,name,include,m,n,p,ingr=None):
    if ingr == None :
        ingr = env.getIngrFromName(name)
    if ingr is not None :
        env.includeIngredientRecipe(ingr, include)
#        ingr.Set(molarity= float(m),
#                     nbMol = int(n),
#                    priority = float(p),)
def excludeIngr(ingr):
    if ingr.o_name != "HIV1_ENV_4nco_0_1_1":
        ingr.recipe.delIngredient(ingr)
        print ingr.o_name,"removed"
 
doit=True    
import datetime
i = datetime.datetime.now()   
#update the json dic if wrong center info
if len(sys.argv) > 1 :
    filename = sys.argv[1]
    resultfile=None
    if filename in autopack.RECIPES :
        n=filename
        v=sys.argv[2]
        filename = autopack.RECIPES[n][v]["setupfile"]
        resultfile= autopack.RECIPES[n][v]["resultfile"]
        setupfile = autopack.retrieveFile(filename,cache="recipes")
    else :
        setupfile =filename
    #overwrite result file
    resultfile="output_autopack_test"
    print ("ok use ",setupfile,filename)
    fileName, fileExtension = os.path.splitext(setupfile)
    print fileName
    n=os.path.basename(fileName)
    recipe=n   
    h = Environment(name=n)     
    h.loadRecipe(setupfile)
    h.setupfile=filename
    if resultfile is not None :
        h.resultfile=resultfile
    h.saveResult = False
#    h.overwriteSurfacePts = False
#    h.compartments[0].overwriteSurfacePts = False
    #build the grid
    #change the grid size ?
    h.smallestProteinSize = 100.0
    h.freePtsUpdateThrehod=0.0
    if doit :
        h.boundingBox=[[ -2482, -2389, -300.26],[ 2495, 2466, 300.02]]
        h.buildGrid(boundingBox=h.boundingBox,gridFileIn=None,rebuild=True ,
                              gridFileOut=None,previousFill=False)
#    h.loopThroughIngr(excludeIngr)
#    h.loopThroughIngr(excludeIngr)
#    h.loopThroughIngr(excludeIngr)
#    h.loopThroughIngr(excludeIngr)    
#    h.loopThroughIngr(setJitter)
#    h.loopThroughIngr(setJitter)
#    h.loopThroughIngr(setJitter)
#    h.loopThroughIngr(setJitter)
    
        h.fill5(verbose = 0,usePP=False)
        h.collectResultPerIngredient()
        rname=fileName+"_"+i.isoformat()+".json"
        rname="C:\Users\ludovic\Documents\CellPackViewer_Cluster\Data\HIV\cellPACK\mycoDNA1.json"
        h.saveRecipe(rname,useXref=True,mixed=True,
                             kwds=["source","name"],result=True,
                           grid=False,packing_options=False,indent=False,quaternion=True)  
        #h.saveRecipe(fileName+"_"+i.isoformat()+".json",useXref=True,mixed=True,
        #             kwds=["source","name","positions","radii"],result=True,
         #          grid=False,packing_options=False,indent=False,quaternion=True)
    else :       
        h.saveRecipe(fileName+"_"+i.isoformat()+".json",useXref=True,mixed=True,result=False,
                       grid=True,packing_options=True,indent=True,quaternion=True)
                   
                   
#h.saveRecipe(fileName+"_mixed_pdb.json",useXref=useXref,mixed=True,
#                     kwds=["source"],result=True,
#                   grid=False,packing_options=False,indent=False,quaternion=True)                   
    #add nucleocapside...
    #add lipids ?
    #RNA?
                   
    #positions and radii?
#    h.saveRecipe("result_mixed_spheres.json",useXref=True,mixed=True,
#                     kwds=["source","positions","radii"],result=True,
#                   grid=False,packing_options=False,indent=False,quaternion=True)
    #optionally, but should force save the sphereTreee
                   
#    f="/home/ludo/Dev/testPanda/result_mixed_spheres.json"
#    with open(f, 'r') as fp :#doesnt work with symbol link ?
#        jsondic=json.load(fp,object_pairs_hook=OrderedDict)#,indent=4, separators=(',', ': ')
#    for o in h.compartments:
#        rs = o.surfaceRecipe
#        if rs :
#            adic = jsondic["compartments"][o.name]["surface"]["ingredients"]
#            for ingr in rs.ingredients:
##                print ingr.source
##                print ingr.positions
##                print ingr.radii
##                print adic[ingr.o_name]["source"] 
#                adic[ingr.o_name]["source"].update(ingr.source)
#                adic[ingr.o_name]["positions"] = ingr.positions
#                adic[ingr.o_name]["radii"] = ingr.radii
##                print jsondic["compartments"][o.name]["surface"]["ingredients"][ingr.o_name]["radii"] 
#        ri = o.innerRecipe
#        if ri :
#            adic = jsondic["compartments"][o.name]["interior"]["ingredients"]
#            for ingr in ri.ingredients:
##                print ingr.source
##                print adic[ingr.o_name]["source"]
#                adic[ingr.o_name]["source"].update(ingr.source)
#                adic[ingr.o_name]["positions"] = ingr.positions
#                adic[ingr.o_name]["radii"] = ingr.radii
##                print jsondic["compartments"][o.name]["interior"]["ingredients"][ingr.o_name]["radii"]
#    output = f
#    with open(output, 'w') as fp :#doesnt work with symbol link ?
#       json.dump(jsondic,fp,separators=(',', ':'))#,indent=4, separators=(',', ': ')

#euler_from_matrix