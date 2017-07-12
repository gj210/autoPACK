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
import numpy as np
np.seterr(all='raise')
#panda3d is there :
#sys.path.append("/usr/lib/python2.7/dist-packages/")
#LINUX
sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs")
#Windows
#sys.path.append("C:\Program Files\MAXON\CINEMA 4D R16 Demo\plugins\ePMV\mgl64\MGLToolsPckgs")
#sys.path.append("C:\Program Files\MAXON\CINEMA 4D R16 Demo\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")
#sys.path.append("C:\Program Files\MAXON\CINEMA 4D R16 Demo\plugins\ePMV\mgl64\MGLToolsPckgs\lib-tk")
##run with C:\>C:\Python26\python.exe "C:\Program Files\MAXON\CINEMA 4D R16 Demo\plugins\ePMV\mgl64\MGLToolsPckgs\autopack\scripts\run_autoPACK.py"  "C:\Users\ludovic\Downloads\Mycoplasma1.5_1.json"

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
 
def removeSphereFile(ingr):
    ingr.sphereFile = None


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
    #any previous results ? build around previous ?
    h.setupfile=filename
    if resultfile is not None :
        h.resultfile=resultfile
    h.saveResult = False
#    h.overwriteSurfacePts = False
#    h.compartments[0].overwriteSurfacePts = False
    #build the grid
    #change the grid size ?
    #h.smallestProteinSize = 50.0
    h.freePtsUpdateThrehod=0.0

    #V≈1.13×1011
    #24848.13
    #DNA
    dna_data = np.loadtxt("/opt/data/dev/flex/data/tps_path_all.txt")   
    dna_ingr = h.getIngrFromName("DNA")
    #dna_ingr = h.compartments[1].innerRecipe.ingredients[1]
    dna_ingr.nbCurve = 1
    dna_ingr.listePtLinear = dna_data.tolist()
    dna_ingr.completion = 1.0
    dna_ingr.nbMol = 0
    dna_ingr.is_previous = True
    
    rotMatj,jtrans=dna_ingr.getJtransRot(np.array(dna_data[0]).flatten(),dna_data[1])# for all points ?
    dna_ingr.update_data_tree(jtrans,rotMatj,pt1=dna_data[0],pt2=dna_data[1],updateTree=False)
    #self.update_data_tree(numpy.array(pt2).flatten(),rotMatj,pt1=pt2,pt2=newPt)
    for i in range(1,len(dna_data)-1):
        rotMatj,jtrans=dna_ingr.getJtransRot(np.array(dna_data[i]).flatten(),dna_data[i+1])# for all points ?
        dna_ingr.update_data_tree(jtrans,rotMatj,pt1=dna_data[i],pt2=dna_data[i+1],updateTree=False)
        if (i % 200) == 0 : 
            print i,len(dna_data)
            
    
    doit = True
    i = datetime.datetime.now()               
    if doit :
        h.smallestProteinSize = 20.0
        #h.boundingBox=[[ -2482, -2389, -500.26],[ 2495, 2466, 500.02]]
        h.buildGrid(boundingBox=h.boundingBox,gridFileIn=None,rebuild=True ,
                              gridFileOut=None,previousFill=True)
        h.fill5(verbose = 1,usePP=False)
        h.collectResultPerIngredient()
        rname=fileName+"_"+i.isoformat()+".json"
#        rname="C:\Users\ludovic\Documents\CellPackViewer_Cluster\Data\HIV\cellPACK\mycoDNA1.json"
#        rname="/home/ludo/Dev/testRun/HIVBloodTest.json"
        rname=fileName+"_result.json"
        h.saveRecipe(rname,useXref=True,mixed=True,
                             kwds=["source","name"],result=True,
                           grid=False,packing_options=False,indent=False,quaternion=True)  
        h.saveRecipe(fileName+"_"+i.isoformat()+".json",useXref=True,mixed=True,
                     kwds=["source","name","positions","radii"],result=True,
                  grid=False,packing_options=False,indent=False,quaternion=True)
    else :       
        h.saveRecipe(fileName+"_"+i.isoformat()+".json",useXref=True,mixed=True,result=False,
                       grid=True,packing_options=True,indent=True,quaternion=True)

#import c4d;h=c4d.af.values()[0].histoVol.values()[0]                
#rname="/opt/data/dev/cellPACK/cellPACK_data_git/cellPACK_database_1.1.0/results/NM_Analysis_FigureA1.0.cpr.json"
#h.saveRecipe(rname,useXref=True,mixed=True,kwds=["source","name","positions","radii"],result=True,grid=False,packing_options=False,indent=False,quaternion=True)                     
#rname="/opt/data/dev/cellPACK/cellPACK_data_git/cellPACK_database_1.1.0/results/NM_Analysis_FigureB1.0.cpr.json"
#rname="/opt/data/dev/cellPACK/cellPACK_data_git/cellPACK_database_1.1.0/results/NM_Analysis_FigureC1.4.cpr.json"
#                     

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