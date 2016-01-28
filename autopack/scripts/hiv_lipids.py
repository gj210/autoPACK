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

helper = autopack.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper = helperClass(vi="nogui")
autopack.helper = helper

from autopack.Environment import Environment
filename = "/home/ludo/hivexp/BloodHIV1.0.json"
#filename ="/home/ludo/hivexp/HIVimature1.0.json"
fileName, fileExtension = os.path.splitext(filename)
n=os.path.basename(fileName)
h = Environment(name=n)
#h.helper = helper
recipe=n
h.loadRecipe(filename)
#previousresult = "/home/ludo/hivexp/BloodHIV1.0_mixed.json"
#previousresult = "/home/ludo/hivexp/BloodHIV1.0_centered.json"
previousresult ="/home/ludo/hivexp/BloodHIV1.0_mixed_tr.json"
r=h.loadResult(previousresult,transpose=True)#
ingredients = h.restore(*r)

from autopack.lipids import lipidsCG
#generate the membrane
lcg = lipidsCG()
comp = h.compartments[0]
lcg.coverShape(comp.vertices,comp.vnormals)
#remove lipids overlaping
posrot=lcg.arrayOfBeads
#save again with lipids ?


#fix HIV1_CA_mono_0_1_0
#center is [405.609,1000.253,377.306]
# if apply the center on all instance ? meanings rotate the center and apply it

#h.saveRecipe("/home/ludo/hivexp/BloodHIV1.0_noref.json",useXref=False,mixed=True,result=False,grid=False,packing_options=True,indent=True,quaternion=True)

h.saveRecipe("/home/ludo/hivexp/BloodHIV1.0_mixed_tr.json",useXref=True,mixed=True,kwds=["source","name"],result=True,grid=False,packing_options=False,indent=False,quaternion=True)
#h.exportToBD_BOX(res_filename="/opt/data/dev/hivbdbox/test_ntr",bd_type="rigid")
#then run ~/Tools/bd_box-2.2/bin/bd_rigid test.prm --out_filename=err.log
from autopack.bd_box import rigid_box as bd_box            
res_filename = "/opt/data/dev/hivbdbox/test_full_hack_tr"
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
        if ingr.name.find("Heparin") != -1 :
            continue
        if ingr.name.find("i1smd_56kD") != -1 :
            continue       
        if ingr.name.find("snake") != -1 :
            continue       
        bd.addAutoPackIngredient(ingr)
#compartment ingr
for orga in h.compartments:
    #skip nucleocapside fro now as it is not properly done
    #compartment surface ingr
#    if orga.name.find("capsid_3j3q_PackOuter") == -1 :
#        continue
    rs =  orga.surfaceRecipe
    if rs :
        for ingr in rs.ingredients:
            bd.addAutoPackIngredient(ingr)
    #compartment matrix ingr
    ri =  orga.innerRecipe
    if ri :
        for ingr in ri.ingredients:
            bd.addAutoPackIngredient(ingr)
bd.write()

#draw sphere {0 0 0} radius 100.0
#pbc box
#pbc set {4000 4000 4000}
#transpose
#
#h.buildGrid(boundingBox=h.boundingBox,gridFileIn=None,rebuild=True ,
#                              gridFileOut="/home/ludo/hivexp/gridOut",previousFill=True)
#
#gag.completion = 1.0
#env.completion = 1.0
##rest should pack around
#gag.nbMol = 0
#env.nbMol = 0
#
#gag.molarity = 0
#env.molarity = 0
#
#h.fill5(verbose = 0,usePP=False)
#h.collectResultPerIngredient()
#rname="/home/ludo/hivexp/hiv_immature_inside.json"
##        rname="C:\Users\ludovic\Documents\CellPackViewer_Cluster\Data\HIV\cellPACK\mycoDNA1.json"
#h.saveRecipe(rname,useXref=True,mixed=True,
#                     kwds=["source","name"],result=True,
#                   grid=False,packing_options=False,indent=False,quaternion=True)
#
#h.saveGridToFile("/home/ludo/hivexp/grid_after_packing")
#transpose ?