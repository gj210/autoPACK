# -*- coding: utf-8 -*-
"""
Created on Wed May  6 10:39:08 2015

@author: ludo
"""
import sys
import os
import json
import numpy as np

np.seterr(all='raise')
# sys.path.append("/usr/lib/python2.7/dist-packages/")
# LINUX
# sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs")
# Windows
sys.path.append("C:\Program Files\MAXON\CINEMA 4D R16 Demo\plugins\ePMV\mgl64\MGLToolsPckgs")
sys.path.append("C:\Program Files\MAXON\CINEMA 4D R16 Demo\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")
sys.path.append("C:\Program Files\MAXON\CINEMA 4D R16 Demo\plugins\ePMV\mgl64\MGLToolsPckgs\lib-tk")
# run with C:\>C:\Python26\python.exe "C:\Program Files\MAXON\CINEMA 4D R16 Demo\plugins\ePMV\mgl64\MGLToolsPckgs\autopack\scripts\run_autoPACK.py"  "C:\Users\ludovic\Downloads\Mycoplasma1.5_1.json"

try:
    from collections import OrderedDict
except:
    from ordereddict import OrderedDict

import autopack
from autopack.Environment import Environment
from autopack.Analysis import AnalyseAP

from json import encoder

encoder.FLOAT_REPR = lambda o: format(o, '.8g')

import upy

helperClass = upy.getHelperClass()
helper = helperClass(vi="nogui")
print helper
autopack.helper = helper
# autopack.fixpath = True
autopack.forceFetch = False
loglevel = 0


def applyIndividualIngredientsOptionsDict(env, parameter_dict, nset):
    for iname in parameter_dict:
        ingr = env.getIngrFromName(iname)
        if ingr is None:
            continue
        for k in parameter_dict[iname]:
            setattr(ingr, k, parameter_dict[iname][k][nset])
            if k == "nbMol":
                ingr.overwrite_nbMol_value = ingr.nbMol


def clear(h, n=0, torestore=None):
    h.reset()
    gfile = None
    # what happen if we dont rebuild the grid
    h.buildGrid(boundingBox=h.boundingBox, gridFileIn=gfile, rebuild=False,
                gridFileOut=None, previousFill=False, lookup=2)


def pack(h, seed, filename):
    clear(h)
    h.fill5(seedNum=seed, verbose=loglevel, usePP=False)
    h.collectResultPerIngredient()
    # we can specify the information we want in the result file
    h.saveRecipe(filename + "_results.json", useXref=False, mixed=True, result=True,
                 kwds=["radii"], grid=False, packing_options=False, indent=False, quaternion=True)
    h.saveRecipe(filename + "_results_tr.json", useXref=False, mixed=True, result=True, transpose=True,
                 kwds=["radii"], grid=False, packing_options=False, indent=False, quaternion=True)


import datetime

i = datetime.datetime.now()

# pass the recipe in command line or force define it
filename = ""
resultfile = ""
setupfile = ""

if len(sys.argv) > 1:
    filename = sys.argv[1]
    resultfile = None
    # we can pass recipe name and version from repo
    if filename in autopack.RECIPES:
        n = filename
        v = sys.argv[2]
        filename = autopack.RECIPES[n][v]["setupfile"]
        resultfile = autopack.RECIPES[n][v]["resultfile"]
        setupfile = autopack.retrieveFile(filename, cache="recipes")
    # or the filename directly
    else:
        setupfile = filename
    # overwrite result file
    resultfile = "output_autopack_test"

print ("ok use ", setupfile, filename)
fileName, fileExtension = os.path.splitext(setupfile)
print fileName
n = os.path.basename(fileName)
recipe = n
h = Environment(name=n)
h.loadRecipe(setupfile)
h.setupfile = filename
if resultfile is not None:
    h.resultfile = resultfile
h.saveResult = False

h.freePtsUpdateThrehod = 0.0
afviewer = None
analyse = AnalyseAP(env=h, viewer=afviewer, result_file=None)
h.analyse = analyse
analyse.g.Resolution = 1.0

# firs time we build the grid
h.buildGrid(boundingBox=h.boundingBox, gridFileIn=None, rebuild=True,
            gridFileOut=None, previousFill=False)

# here we can do experiment one after each other or using a loop after changing some parameters.
# you need to choose which ingredient, and then which variable to sample, for instance here I have 3 set of value:
nset = 3
nseed = 5  # number of seed per set
individuals_ingredients_parameter_set = OrderedDict({
    "ingredient1": OrderedDict({
        'rejectionThreshold': [10, 20, 30],
        'packingPriority': [-10, 0, 10],
        'cutoff_boundary': [0, 10, 20],
        'jitterMax': [[0.1, 0.1, 0.1], [0.5, 0.5, 0.5], [1., 1., 1.]],
        'nbJitter': [5, 10, 15],
        'nbMol': [2, 4, 6, 8],
    }),
    "ingredient2": OrderedDict({'rejectionThreshold': [10, 20, 30],
                                'packingPriority': [-10, 0, 10],
                                'cutoff_boundary': [0, 10, 20],
                                'jitterMax': [[0.1, 0.1, 0.1], [0.5, 0.5, 0.5], [1., 1., 1.]],
                                'nbJitter': [5, 10, 15],
                                'nbMol': [2, 4, 6, 8],}),
    "ingredient3": OrderedDict({'rejectionThreshold': [10, 20, 30],
                                'packingPriority': [-10, 0, 10],
                                'cutoff_boundary': [0, 10, 20],
                                'jitterMax': [[0.1, 0.1, 0.1], [0.5, 0.5, 0.5], [1., 1., 1.]],
                                'nbJitter': [5, 10, 15],
                                'nbMol': [10, 14, 16, 18],})
})

# gather uniq seed
seeds_i = analyse.getHaltonUnique(nseed)

for i in range(nset):
    # apply the your setup
    applyIndividualIngredientsOptionsDict(h, individuals_ingredients_parameter_set, i)
    for seed in seeds_i:
        result_filename = "pathtouniqresultfilename_" + str(seed) + "_" + str(i)
        pack(h, seed, result_filename)

    # # OR change manually some parameter and run
    # h.compartments[0].surfaceRecipe.ingredients[0].nbMol = 0
    # h.compartments[0].surfaceRecipe.ingredients[1].nbMol = 200
    # h.compartments[0].surfaceRecipe.ingredients[2].nbMol = 0
    # seed = 14
    # result_filename = "myrun10_"+str(seed)
    # pack(h, seed, result_filename)
    #
    # # change again some parameters...
