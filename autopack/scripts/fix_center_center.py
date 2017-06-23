# -*- coding: utf-8 -*-
"""
Created on Thu Jun 09 12:17:52 2016

@author: ludov
"""
import sys
import os
import numpy as np
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs\lib-tk")

NOGUI=1
import autopack
from autopack.Environment import Environment
from autopack import IOutils as io

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
recipe= "D:\\Data\\cellPACK_data\\cellPACK_database_1.1.0\\recipes\\HIV_VLP.1.1.json" #127
filename = recipe#"/Users/ludo/DEV/autopack_git/autoPACK_database_1.0.0/recipes/NM_Analysis_FigureA1.0.xml"
fileName, fileExtension = os.path.splitext(filename)
n=os.path.basename(fileName)
h = Environment(name=n)
#h.helper = helper
recipe=n
h.loadRecipe(filename)
#recenter ingredient position according center.
#should manually do the mesh ?
for c in h.compartments:
    for ingr in c.surfaceRecipe.ingredients:
        ingr.positions = np.array(ingr.positions)-np.array(ingr.offset)
        ingr.sphereFile = None
#save the recipe  
h.version="1.2"
h.saveRecipe("D:\\Data\\cellPACK_data\\cellPACK_database_1.1.0\\recipes\\HIV_VLP.1.2.json",useXref=False,format_output="json")