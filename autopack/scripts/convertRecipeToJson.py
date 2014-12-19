# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:59:01 2013

@author: ludo
convert xml recipe to json
could also use to test  recipe if they correctly load 
"""
import sys
import os

Maya= True
if Maya :
    import sys
    sys.path.append("/Users/ludo/Library/Preferences/Autodesk/maya/2015-x64/plug-ins/MGLToolsPckgs")
    sys.path.append("/Users/ludo/Library/Preferences/Autodesk/maya/2015-x64/plug-ins/MGLToolsPckgs/PIL")    
    #maya standalone special
    import maya.standalone
    maya.standalone.initialize()
    #load plugin
    import maya
    maya.cmds.loadPlugin("fbxmaya")
    
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
export_json = True
useXref = False
mixedJson = True

#def convertOneRecipe()

if len(sys.argv) > 1 :
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
    h.resultfile=resultfile
    fileName, fileExtension = os.path.splitext(setupfile)
    if export_json:
        print ("expot json recipe ",fileName)
        h.saveRecipe(fileName+".json",useXref=useXref)
    if check_result:
        rfile = h.resultfile
        resultfilename = autopack.retrieveFile(rfile,cache="results")
        if resultfilename is None :
            print ("no result for "+n+" "+h.version+" "+rfile)
            sys.exit()
        print ("get the result file from ",resultfilename)
        result,orgaresult,freePoint=h.loadResult(resultfilename=resultfilename,
                                                 restore_grid=False,backward=True)#load text ?#this will restore the grid  
        ingredients = h.restore(result,orgaresult,freePoint)
#        print ("json ?",h.result_json)
        if export_json :
            fileName, fileExtension = os.path.splitext(resultfilename)
            h.store_asJson("/Users/ludo/"+n+"_result.json",indent=False)#fileName+"1.json",indent=False)     
            print ("ok ","/Users/ludo/"+n+"_result.json",fileName+".json")
    if mixedJson and export_json:
        h.setupfile=filename
        h.saveRecipe(fileName+"_mixed.json",useXref=useXref,mixed=True,
                     kwds=["compNum","encapsulatingRadius"],result=True,
                   grid=False,packing_options=False,indent=False)
        print ("mixed file ",filename)
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
#I usually run this on with pmv,anaconda or mayapy