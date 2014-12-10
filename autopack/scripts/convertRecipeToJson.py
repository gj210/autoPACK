# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:59:01 2013

@author: ludo
convert xml recipe to json
could also use to test  recipe if they correctly load 
"""
import sys
import os

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
check_result = True
export_json = False
useXref = True
if len(sys.argv) > 1 :
    filename = sys.argv[1]
    resultfile=None
    if filename in autopack.RECIPES :
        n=filename
        v=sys.argv[2]
        filename = autopack.RECIPES[n][v]["setupfile"]
        resultfile= autopack.RECIPES[n][v]["resultfile"]
    setupfile = autopack.retrieveFile(filename,cache="recipes")
    print ("ok use ",setupfile)
    fileName, fileExtension = os.path.splitext(setupfile)
    n=os.path.basename(fileName)
    h = Environment(name=n)     
    h.loadRecipe(setupfile)
    h.setupfile=filename
    h.resultfile=resultfile
    fileName, fileExtension = os.path.splitext(setupfile)
    if export_json:
        h.saveRecipe(fileName+".json",useXref=useXref)
    if check_result:
        rfile = h.resultfile
        resultfilename = autopack.retrieveFile(rfile,cache="results")
        if resultfilename is None :
            print ("no result for "+n+" "+h.version+" "+rfile)
            sys.exit()
        result,orgaresult,freePoint=h.load(resultfilename=resultfilename,restore_grid=False)#load text ?#this will restore the grid  
        ingredients = h.restore(result,orgaresult,freePoint)
#        print ("json ?",h.result_json)
        if export_json :
            fileName, fileExtension = os.path.splitext(resultfilename)
            h.store_asJson(n+"_results.json",indent=False)#fileName+"1.json",indent=False)     
            print ("ok ",n+"_results.json",fileName+".json")
else :    
    for recipe in autopack.RECIPES:
        rname = recipe
        recipes = autopack.RECIPES[recipe]
        for version in recipes:
            print rname,version,autopack.RECIPES[recipe][version]["setupfile"]
            setupfile = autopack.RECIPES[recipe][version]["setupfile"]
            fileName, fileExtension = os.path.splitext(setupfile)
            if fileExtension != ".json":
                setupfile = autopack.retrieveFile(setupfile,cache="recipes") 
                if setupfile is None :
                    print ("setup file doesnt exist ",rname,version,setupfile, autopack.RECIPES[recipe][version]["setupfile"])
                    continue
                print (fileName, fileExtension)
                n= rname
                h = Environment(name=n)            
                recipe=n
                h.loadRecipe(setupfile)
                fileName, fileExtension = os.path.splitext(setupfile)
    #            if not os.path.isfile(fileName+".json"):
                if export_json:
                    h.saveRecipe(fileName+".json")
            if check_result:
                rfile = autopack.RECIPES[recipe][version]["resultfile"]
                resultfilename = autopack.retrieveFile(h.resultfile+".json",cache="results")
                if resultfilename is None :
                    print ("no result for "+rname+" "+version+" "+rfile)
                    continue
                result,orgaresult,freePoint=h.load(resultfilename=resultfilename,restore_grid=False)#load text ?#this will restore the grid  
                ingredients = h.restore(result,orgaresult,freePoint)
                if export_json :
                    fileName, fileExtension = os.path.splitext(resultfilename)
                    h.store_asJson(fileName+".json",indent=False)            
                #what about result ?
    #            break
    #    break
#execfile("/Users/ludo/DEV/autoPACK_github/autopack/scripts/testAF_xml.py")       