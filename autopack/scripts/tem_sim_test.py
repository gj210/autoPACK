# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 21:58:54 2014

@author: ludo
execfile("C:\\Users\\ludov\\OneDrive\\Documents\\autoPACK\\autopack\\scripts\\tem_sim_test.py")
"""
import sys
import os

from autopack.tem_sim import tem_sim
import autopack
#wrkDir = AutoFill.__path__[0]
localdir = wrkDir = autopack.__path__[0]
from autopack.Environment import Environment
helper = autopack.helper
NOGUI = False
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
filename = "C:\\Dev\\hiv_liposome\\forBrett\\HIVspherejson.json"
resultfile="C:\\Dev\\hiv_liposome\\forBrett\\HIVspherejson_result.json"
fileName, fileExtension = os.path.splitext(filename)
n=os.path.basename(fileName)
h = Environment(name=n)     
h.loadRecipe(filename)
h.setupfile=filename
if resultfile is not None :
    h.resultfile=resultfile
result,orgaresult,freePoint=h.loadResult(resultfilename=resultfile)
ingredients = h.restore(result,orgaresult,freePoint)   
tem = tem_sim(name= "C:\\Dev\\hiv_liposome\\forBrett\\HIVliposome",bounding_box=[[-1500, -1500, -200], [1500.0, 1500.0, 200]])
tem.setup()
def setupFromEnv(self,env):
    r =  env.exteriorRecipe
    if r :
        for ingr in r.ingredients:
            self.addAutoPackIngredient(ingr)

    #compartment ingr
    for orga in env.compartments:
        #compartment surface ingr
        rs =  orga.surfaceRecipe
        if rs :
            for ingr in rs.ingredients:
                self.addAutoPackIngredient(ingr)
        #compartment matrix ingr
        ri =  orga.innerRecipe
        if ri :
            for ingr in ri.ingredients:
                self.addAutoPackIngredient(ingr)
setupFromEnv(tem,h)
tem.write()