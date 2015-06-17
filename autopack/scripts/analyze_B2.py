# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:59:01 2013

@author: ludo
Hoe to run autopack from xml file
"""
import os
import sys
import json
import numpy as np
import math
import time
from collections import OrderedDict
#we need the modules from MGLToolsPckgs available from the command Lines
#change accordingly your system, on Mac we can use the one distributed with C4D-ePMV
#sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14_95696CEA/plugins/ePMV/mgl64/MGLToolsPckgs")
#on Linux we can use the regular MGLTools
sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs")

#for plotting results
from matplotlib import pyplot as plt
from matplotlib.patches import Circle

import autopack
#where is autopack
localdir = wrkDir = autopack.__path__[0]

#autopack setup
from autopack.Environment import Environment
from autopack.Ingredient import KWDS as ingredients_parameters
from autopack.Graphics import AutopackViewer as AFViewer
from autopack.Analysis import AnalyseAP

#get the helper
#deal with geometry and such
helper = autopack.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper = helperClass(vi="nogui")   
autopack.helper = helper
#==============================================================================
# some usefule function
#==============================================================================

def plotOneResult2D(h,filename):
    width = h.boundingBox[1][0]#should be the boundary here ?
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bbox = h.boundingBox
    radius={}
    ingrrot={}
    ingrpos={}
    pos_rad = [[np.array(m[0]).tolist(),m[2].encapsulatingRadius,m[2].color] for m in h.molecules]#jtrans, rotMatj, self, ptInd
    for ingr in pos_rad: 
        p=ingr[0]
        r=ingr[1]
#            for i,p in enumerate(ingrpos[ingr]): 
#                print (p,radius[ingr])
        ax.add_patch(Circle((p[0], p[1]),r,
                        edgecolor="black", facecolor=ingr[2])) 
        if autopack.testPeriodicity :
            if p[0] < r:
                ax.add_patch(Circle((p[0] + width,p[1]), r, facecolor=ingr[2]))
            elif p[0] > (width-r):
                ax.add_patch(Circle((p[0] - width,p[1]), r, facecolor=ingr[2]))
            if p[1] < r:
                ax.add_patch(Circle((p[0],p[1]+ width), r, facecolor=ingr[2]))
            elif p[1] > (width-r):
                ax.add_patch(Circle((p[0],p[1] - width), r, facecolor=ingr[2]))
        ax.set_aspect(1.0)
    plt.axhline(y=bbox[0][1], color='k')
    plt.axhline(y=bbox[1][1], color='k')
    plt.axvline(x=bbox[0][0], color='k')
    plt.axvline(x=bbox[1][0], color='k')
    plt.axis([bbox[0][0], bbox[1][0],
                 bbox[0][1], bbox[1][1]])
    plt.savefig(filename)
    plt.close('all')     # closes the cu
    
def clear(h,n=0,torestore=None):
    h.reset()
    gfile = None
    h.buildGrid(boundingBox=h.boundingBox,gridFileIn=gfile,rebuild=True ,
                gridFileOut=None,previousFill=False)

def activenmolecules(h,n):
    for i,ingr in enumerate(h.exteriorRecipe.ingredients):
        if i > n-1 :
            ingr.completion=1.0
            ingr.nbMol = 0
        else :
            ingr.completion=0.0

def pack(h,seed,filename,eid):
    h.fill5(seedNum=seed,verbose = 0,usePP=False)
    h.collectResultPerIngredient()
    todump=[[],[]]
    #we can save only the object positions
    todump[0] = [np.array(m[0]).tolist() for m in h.molecules]# [[ingr.nbMol,len(res)],[]]  POS     
    with open(filename, 'w') as fp :
        json.dump(todump,fp)      
#    h.store_asJson(resultfilename=filename+"_results.json")
    #we can specify the information we want in the result file    
    h.saveRecipe(filename+"_results.json",useXref=False,mixed=True,result=True,
                   kwds=["radii"],grid=False,packing_options=False,indent=False,quaternion=True)
                   
def one_exp(h,seed,eid=0,setn=1,periodicity=True,output=None):
    clear(h)
    if h.use_periodicity:
        if output is None :
            output="/Users/ludo/DEV/autoPACKresults/NM_Analysis_3B_P_"
    else :
        if output is None :
            output="/Users/ludo/DEV/autoPACKresults/NM_Analysis_3B_nP_"
    if not os.path.exists(output):
        os.makedirs(output)
    setn=str(setn).replace("(","").replace(")","").replace(", ","_")
    pack(h,seed,output+os.sep+"run_"+str(setn)+"_"+str(eid),eid)
    plotOneResult2D(h,output+os.sep+"run_"+str(setn)+"_"+str(eid)+"figure.png")

def applyGeneralOptions(env,parameter_set,nset):
    for k in parameter_set :
        setattr(env,k,parameter_set[k][nset])
        
def applyGeneralIngredientsOptions(env,paremeter_set,nset):
    #apply the key options to all ingredients
    r = env.exteriorRecipe
    if r :
        for ingr in r.ingredients:
            for k in paremeter_set :
                setattr(ingr,k,paremeter_set[k][nset])
            
def applyIndividualIngredientsOptions(env,paremeter_liste,nset):
    #apply the key options to specified ingredients
    r = env.exteriorRecipe
    if r :
        for i,params in enumerate(paremeter_liste) :
            for k in params :
                setattr(r.ingredients[i],k,params[k][nset])
                
def applyGeneralOptions_product(env,parameter_set,nset):
    for i,k in enumerate(parameter_set) :
        setattr(env,k,parameter_set[k][nset[i]])
        
def applyGeneralIngredientsOptions_product(env,paremeter_set,nset,offset=0):
    #apply the key options to all ingredients
    r = env.exteriorRecipe
    if r :
        for ingr in r.ingredients:
            for i,k in enumerate(paremeter_set) :
                setattr(ingr,k,paremeter_set[k][nset[i+offset]])   
                             
#==============================================================================
# autoPACK setup, and recipe loading
#==============================================================================
#define wher is out recipe setup file
#setupfile = "/Users/ludo/DEV/autopack_git/autoPACK_database_1.0.0/recipes/NM_Analysis_FigureB1.1.xml"
#force downloading the latest recipe file
autopack.forceFetch=True
#Or we can also use the server for gathering the file
recipe = "NM_Analysis_FigureB"
version = "1.0"
filename = autopack.RECIPES[recipe][version]["setupfile"]
resultfile= autopack.RECIPES[recipe][version]["resultfile"]

setupfile = autopack.retrieveFile(filename,cache="recipes")

fileName, fileExtension = os.path.splitext(setupfile)
n=os.path.basename(fileName)
h = Environment(name=n)
h.loadRecipe(setupfile)
afviewer=None
h.placeMethod="pandaBullet"
h.encapsulatingGrid=0
#default for periodicity
h.use_periodicity = False
autopack.testPeriodicity = True
autopack.biasedPeriodicity = [1,1,0]
analyse = AnalyseAP(env=h, viewer=afviewer, result_file=None)
h.analyse = analyse
analyse.g.Resolution = 1.0
#h.smallestProteinSize=30.0#get it faster? same result ?
h.boundingBox=np.array(h.boundingBox)


#==============================================================================
# the parameter for a run
#==============================================================================
#all option can be access through a Dictionary that hold :
#options name
#options default value
#options value range
#options type e.g. Vector,Float, Int,String
#to get the value for the current value we can use the python function getattr
#to change the value of the option we can use the python function setattr
#all the packing options (General Parameters) are found in h.OPTIONS
for k in h.OPTIONS.keys():
    print h.OPTIONS[k]
#the main options to explore are probably:
#smallestProteinSize float
#pickWeightedIngr    boolean
#pickRandPt          boolean
#use_periodicity     boolean

#all the ingredients/Object options (Individual Parameters) are found in Ingredients.KWDS
for k in ingredients_parameters.keys():
    print k,ingredients_parameters[k]
#for instance packingPriority:
# - an optional packing priority. If omited the priority will be based
#    on the radius with larger radii first
#    ham here: (-)packingPriority object will pack from high to low one at a time
#    (+)packingPriority will be weighted by assigned priority value
#    (0)packignPriority will be weighted by complexity and appended to what is left
#    of the (+) values
    


#You should here prepare your parameter set I guess
#one set should be a dictionary for the packing
#one or several dictionary for the ingredients
#Comment uncomment here
packing_parameter_set=OrderedDict({
#"smallestProteinSize":[5.0,10.0,15.0],
"pickWeightedIngr":[1,0],
"pickRandPt":[1,0],
"use_periodicity":[1,0]
})
#ingredients options you want to change for all objects
ingredients_paremeter_set=OrderedDict({
'rejectionThreshold':[10,20,100],
#'packingPriority':[-10,0,10],
#'cutoff_boundary':[0,10,20],
'jitterMax':[[0.1,0.1,0.0],[1.0,1.0,0.0],[1.,1.,0.]],
'nbJitter':[10,20,30],
})










#individuals options. This recipe has 5objetcs, thus you need 5 set of values
#use a dictionary or a list
individuals_ingredients_paremeter_set=OrderedDict({
"ext__Sphere_radius_100":OrderedDict({
'rejectionThreshold':[10,20,30],
'packingPriority':[-10,0,10],
'cutoff_boundary':[0,10,20],
'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
'nbJitter':[5,10,15],
'nbMol':[2,4,6,8],
}),
"ext__Sphere_radius_200":OrderedDict({'rejectionThreshold':[10,20,30],
'packingPriority':[-10,0,10],
'cutoff_boundary':[0,10,20],
'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
'nbJitter':[5,10,15],
'nbMol':[2,4,6,8],}),
"ext__Sphere_radius_50":OrderedDict({'rejectionThreshold':[10,20,30],
'packingPriority':[-10,0,10],
'cutoff_boundary':[0,10,20],
'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
'nbJitter':[5,10,15],
'nbMol':[10,14,16,18],}),
"ext__Sphere_radius_25p":OrderedDict({'rejectionThreshold':[10,20,30],
'packingPriority':[-10,0,10],
'cutoff_boundary':[0,10,20],
'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
'nbJitter':[5,10,15],
'nbMol':[20,40,60,80],}),
"ext__Sphere_radius_25r":OrderedDict({'rejectionThreshold':[10,20,30],
'packingPriority':[-10,0,10],
'cutoff_boundary':[0,10,20],
'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
'nbJitter':[5,10,15],
'nbMol':[20,40,60,80],}),
})
individuals_ingredients_paremeter_liste=[
OrderedDict({
'rejectionThreshold':[10,20,30],
'packingPriority':[-10,0,10],
'cutoff_boundary':[0,10,20],
'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
'nbJitter':[5,10,15],
'nbMol':[2,4,6,8],
}),
OrderedDict({'rejectionThreshold':[10,20,30],
'packingPriority':[-10,0,10],
'cutoff_boundary':[0,10,20],
'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
'nbJitter':[5,10,15],
'nbMol':[2,4,6,8],}),
OrderedDict({'rejectionThreshold':[10,20,30],
'packingPriority':[-10,0,10],
'cutoff_boundary':[0,10,20],
'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
'nbJitter':[5,10,15],
'nbMol':[5,10,15],}),
OrderedDict({'rejectionThreshold':[10,20,30],
'packingPriority':[-10,0,10],
'cutoff_boundary':[0,10,20],
'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
'nbJitter':[5,10,15],
'nbMol':[20,40,60,80],}),
OrderedDict({'rejectionThreshold':[10,20,30],
'packingPriority':[-10,0,10],
'cutoff_boundary':[0,10,20],
'jitterMax':[[0.1,0.1,0.1],[0.5,0.5,0.5],[1.,1.,1.]],
'nbJitter':[5,10,15],
'nbMol':[20,40,60,80],}),
]




#a variable that define if you want to use the individual option
use_individual_options = False 

#there is a lot of different to manipulate the options, I just show you the dictionary way.
#you can also comment/uncomment the option you don't want
#
#==============================================================================
# experiment running
#==============================================================================
#prepapre an array of 200 seed
seeds_i = analyse.getHaltonUnique(200)

#define an ouput directory
output="/home/ludo/Dev/testRun/"
#numberofRun per numberofSet containing numberofParameter
numberofParameter=len(packing_parameter_set)+len(ingredients_paremeter_set)
numberofSet=2 #because of the binary option, its correspondent to the possible number for the options,  e.g. True/False or 1,5 etc..
numberofRun=1 #the number of seeds to use for one set, one set is a combination of options

usecombinatorial=True

time_benchmark_file=output+"time_check.txt"
timer=0.0
time_data=""

if usecombinatorial:
    #TODO:
    #ideally we will have the combinatorial sort per possible range for option
    #start with binary
    #then smaller to bigger range
    from itertools import product
    #we can generate the number of set by  combinatorial
    sets=product(range(numberofSet),repeat=numberofParameter)
    for s in sets :
        count=0
        print s,str(s).replace("(","").replace(")","").replace(", ","_")
        applyGeneralOptions_product(h,packing_parameter_set,s)
        applyGeneralIngredientsOptions_product(h,ingredients_paremeter_set,s,offset=len(packing_parameter_set))
        for seed in seeds_i[:numberofRun] :
            print "seed",seed
            timer=time.time()
            one_exp(h,seed,eid=count,setn=s,periodicity=False,output=output)
            count+=1  
            dT=time.time()-timer
            time_data+=str(count)+"\t"+str(seed)+"\t"+str(dT)+"\n"
else :
    for pset in range(numberofSet):
        count=0
        #apply the options for the given set
        applyGeneralOptions(h,packing_parameter_set,pset)
        applyGeneralIngredientsOptions(h,ingredients_paremeter_set,pset)
        if use_individual_options:
            applyIndividualIngredientsOptions(h,individuals_ingredients_paremeter_liste,pset)    
        #do the X run
        for seed in seeds_i[:numberofRun] :
            print "seed",seed
            timer=time.time()
            one_exp(h,seed,eid=count,setn=pset,periodicity=False,output=output)
            count+=1
            dT=time.time()-timer
            time_data+=str(count)+"\t"+str(seed)+"\t"+str(dT)+"\n"
f=open(time_benchmark_file,"w")
f.write(time_data)
f.close()
#execfile("analyze_B2.py")