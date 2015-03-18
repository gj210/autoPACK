# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 22:08:01 2014

@author: ludo
"""
import sys
import os
import math

##The program will read and write all the packings in the following file format: It should be a binary file which stores sequentially sphere center x, y, z coordinates and diameter for each particle as floating points in double precision in little-endian byte order. If the machine on which the program is being run is big-endian, the program will detect it and will still read and write in little-endian format. 

import numpy
import autopack
#wrkDir = AutoFill.__path__[0]
localdir = wrkDir = autopack.__path__[0]

from autopack.Environment import Environment

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
    print helper
autopack.helper = helper
autopack.fixpath = True

load=False
pack=True
if len(sys.argv) > 1 :#and doit :
    filename = sys.argv[1]
    resultfile=None
    if filename in autopack.RECIPES :
        n=filename
        v=sys.argv[2]
        filename = autopack.RECIPES[n][v]["setupfile"]
        resultfile= autopack.RECIPES[n][v]["resultfile"]
    setupfile = autopack.retrieveFile(filename,cache="recipes")
    print ("ok use ",setupfile,filename)
    fileName, fileExtension = os.path.splitext(setupfile)
    n=os.path.basename(fileName)
    h = Environment(name=n)     
    h.loadRecipe(setupfile)
    h.setupfile=filename
    if resultfile is not None :
        h.resultfile=resultfile
    fileName, fileExtension = os.path.splitext(setupfile)
    if load :
        rfile = h.resultfile
        resultfilename = autopack.retrieveFile(rfile,cache="results")
        if resultfilename is None :
            print ("no result for "+n+" "+h.version+" "+rfile)
            sys.exit()
        print ("get the result file from ",resultfilename)
        result,orgaresult,freePoint=h.loadResult(resultfilename=resultfilename)
    #                                             restore_grid=False,backward=True)#load text ?#this will restore the grid  
        ingredients = h.restore(result,orgaresult,freePoint)
    elif pack :
        #pack
        gridfile = "/home/ludo/Dev/bdbox/grid_store"
        h.placeMethod="RAPID"#does it use all vertices ?
        h.saveResult = True
        h.resultfile = "/home/ludo/Dev/bdbox/"+n+v
        h.innerGridMethod = "bhtree"#jordan pure python ? sdf ?
    #    h.boundingBox = [[-250.0, -6500.0/2.0, -250.0], [250.0, 6500.0/2.0, 250.0]]
#        h.boundingBox =[[-2482, -2389.0, 100.0], [2495, 2466, 2181.0]]
        h.buildGrid(boundingBox=h.boundingBox,gridFileIn=gridfile,rebuild=True ,
                          gridFileOut=None,previousFill=False)
        h.fill5(verbose = 0,usePP=True)
        h.collectResultPerIngredient()
        h.saveRecipe(h.resultfile+"_mixed.json",useXref=True,mixed=True,
                     kwds=["compNum","encapsulatingRadius"],result=True,
                   grid=False,packing_options=False,indent=False)
    #display ?
    #export the complete recipe as collada. each ingredient -> meshnode. Each instance->node instance
    env=h
    from autopack import bd_box
    #bd_box.executable_dmatrix="dmatrix"
    #env.exportToBD_BOX(res_filename="/home/ludo/Dev/bdbox/"+n+v,bd_type="rigid")
    #require no overlapp at all...
    
#execfile("pathto/export_recipe_collada.py") 
#I usually run this on with pmv,anaconda or mayapy
                

#execfile("/Users/ludo/DEV/git_upy/examples/export_collada.py")
#import upy
#helper = upy.getHelperClass()()
#helper.read("/Users/ludo/DEV/autopack_git/autoPACK_database_1.0.0/geometries/HIV1_capside_3j3q_Rep_Med_0_2_1.dae")
#