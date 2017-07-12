# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:59:01 2013

@author: ludo
Hoe to run autopack from xml file
"""
# Show the effect of garbage collection
import sys
import os

sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/")
sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/PIL/")

import autopack
#wrkDir = AutoFill.__path__[0]
localdir = wrkDir = autopack.__path__[0]
path="/opt/data/dev/BD_EXP/"

helper = autopack.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper = helperClass(vi="nogui")
autopack.helper = helper

from autopack.Environment import Environment
filename = path+"smallBlood.json"
fileName, fileExtension = os.path.splitext(filename)
n=os.path.basename(fileName)
h = Environment(name=n)
#h.helper = helper
recipe=n
h.loadRecipe(filename)
h.placeMethod="pandaBullet"



rname=path+"smallblood_results.json"

#h.buildGrid(boundingBox=h.boundingBox,gridFileIn=None,rebuild=True ,gridFileOut=path+"gridOut",previousFill=False,lookup=2)
#
#def setNB(ingr):
#    ingr.nbMol = 0
#    ingr.molarity = 0.0
#    ingr.overwrite_nbMol_value = 0.0
#    ingr.overwrite_nbMol = True
#
#h.loopThroughIngr(setNB)
#h.loopThroughIngr(setNB)
#h.loopThroughIngr(setNB)
#h.loopThroughIngr(setNB)
#
#h.exteriorRecipe.ingredients[0].nbMol=2
#h.exteriorRecipe.ingredients[0].overwrite_nbMol_value=1
#
#h.exteriorRecipe.ingredients[1].nbMol=2
#h.exteriorRecipe.ingredients[1].overwrite_nbMol_value=1
#
#h.exteriorRecipe.ingredients[2].nbMol=2
#h.exteriorRecipe.ingredients[2].overwrite_nbMol_value=1

#h.fill5(verbose = 0,usePP=False)

r=h.loadResult(rname,transpose=True)#)# or not ?
ingredients = h.restore(*r)
#
h.collectResultPerIngredient()

h.saveRecipe(path+"smallblood_results_tr.json",useXref=True,mixed=True,kwds=["source","name"],result=True,grid=False,packing_options=False,indent=False,quaternion=True)
#h.exportToBD_BOX(res_filename=path+"test_tr",bd_type="rigid")#transpose ?
from autopack.bd_box import rigid_box as bd_box            
res_filename = path+"test_tr"
bd=bd_box(res_filename,bounding_box=h.boundingBox)
bd.offset = [0,0,0]
#bd.params["xbox"]["value"]=4000.0
#bd.params["ybox"]["value"]=4000.0
#bd.params["zbox"]["value"]=4000.0
bd.makePrmFile()
r =  h.exteriorRecipe
#bd.addAutoPackIngredient(r.ingredients[1])
#bd.addAutoPackIngredient(r.ingredients[2])
#bd.addAutoPackIngredient(r.ingredients[3])
if r :
    for ingr in r.ingredients:
#        if ingr.name.find("Heparin") != -1 :
#            continue
#        if ingr.name.find("i1smd_56kD") != -1 :
#            continue       
#        if ingr.name.find("snake") != -1 :
#            continue       
        bd.addAutoPackIngredient(ingr)
#compartment ingr
#for orga in h.compartments:
#    #compartment surface ingr
#    rs =  orga.surfaceRecipe
#    if rs :
#        for ingr in rs.ingredients:
#            bd.addAutoPackIngredient(ingr)
#    #compartment matrix ingr
#    ri =  orga.innerRecipe
#    if ri :
#        for ingr in ri.ingredients:
#            bd.addAutoPackIngredient(ingr)
bd.write()

#
#from bd_box import rigid_box as bd_box            
#res_filename = path+"test_tr"
#bd=bd_box(res_filename,bounding_box=h.boundingBox)
#r =  h.exteriorRecipe
#if r :
#    for ingr in r.ingredients:
#        bd.addAutoPackIngredient(ingr)
##compartment ingr
#for orga in h.compartments:
#    #compartment surface ingr
#    rs =  orga.surfaceRecipe
#    if rs :
#        for ingr in rs.ingredients:
#            bd.addAutoPackIngredient(ingr)
#    #compartment matrix ingr
#    ri =  orga.innerRecipe
#    if ri :
#        for ingr in ri.ingredients:
#            bd.addAutoPackIngredient(ingr)
#bd.write()


#-36.4229 19.1282 -84.4719
#-15.1086 -38.0522 83.3628

#then run ~/Tools/bd_box-2.2/bin/bd_rigid test.prm --out_filename=err.log
#TODO 
#FIX Heparin
#FIX i1smd_56kD