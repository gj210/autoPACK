# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:59:01 2013

@author: ludo
Hoe to run autopack from xml file
"""
# Show the effect of garbage collection
import sys
import json
import numpy
import numpy as np
import math
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14_95696CEA/plugins/ePMV/mgl64/MGLToolsPckgs")
#should have dejavu...
import gc
import pprint
for i in range(2):
    print 'Collecting %d ...' % i
    n = gc.collect()
    print 'Unreachable objects:', n
    print 'Remaining Garbage:', 
    pprint.pprint(gc.garbage)
    del gc.garbage[:]
    print
from matplotlib import pyplot as plt
from matplotlib.patches import Circle

import os
#import sys
#sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14 Student_7B992864/plugins/ePMV/mgl64/MGLToolsPckgs")
#sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14 Student_7B992864/plugins/ePMV/mgl64/MGLToolsPckgs/PIL/")
import autopack
#wrkDir = AutoFill.__path__[0]
localdir = wrkDir = autopack.__path__[0]

from autopack.Environment import Environment
from autopack.Graphics import AutopackViewer as AFViewer
from autopack.Analysis import AnalyseAP
TWOD = 1
NOGUI = 0
ANALYSIS = 1
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
#filename = "/Users/ludo/Desktop/cell.xml"
#filename = wrkDir+os.sep+"autoFillRecipeScripts/2DsphereFill/Test_Spheres2D1.1.xml"
#filename = wrkDir+os.sep+"autoFillRecipeScripts/2DsphereFill/Test_Spheres2Dgradients1.0.xml"
#filename = "/Users/ludo/DEV/autopack_git/data/Mycoplasma/recipe/Mycoplasma1.3.xml"
#filename = "/Users/ludo/Desktop/cell_hack.xml"
filename = "/Users/ludo/DEV/autopack_git/autoPACK_database_1.0.0/recipes/NM_Analysis_FigureB1.0.xml"
filename = "/Users/ludo/DEV/autopack_git/autoPACK_database_1.0.0/recipes/NM_Analysis_FigureB1.1.xml"

fileName, fileExtension = os.path.splitext(filename)
n=os.path.basename(fileName)
h = Environment(name=n)
#h.helper = helper
recipe=n
h.load_XML(filename)
afviewer=None
if not NOGUI :
    print h,helper
    setattr(h,"helper",helper)
    afviewer = AFViewer(ViewerType=h.helper.host,helper=h.helper)
    afviewer.SetHistoVol(h,20.0,display=False)
    h.host=h.helper.host
    afviewer.displayPreFill()
h.saveResult = False

#resultfilename = h1.resultfile = wrkDir+os.sep+"autoFillRecipeScripts"+os.sep+"2DsphereFill"+os.sep+"results"+os.sep+"SpherefillResult.afr"
#resultfilename = h.resultfile = wrkDir+os.sep+"autoFillRecipeScripts/2DsphereFill/results/2DsphereFill_1.1.apr"
#resultfilename = h.resultfile = wrkDir+os.sep+"autoFillRecipeScripts/Mycoplasma/results/MycoplasmaPackResult_3"
resultfilename = h.resultfile = "/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_rev/NM_Analysis_B"

#h.smallestProteinSize=15
#h.exteriorRecipe.ingredients[0].uLength = 100.0
#overwrite the jiterMax usin jitterMax = [d[k]["rad"]/(15.*1.1547),d[k]["rad"]/(15.*1.1547),0.0]
#loopThroughIngr
def setJitter(ingr):
    ingr.jitterMax =[ingr.encapsulatingRadius,ingr.encapsulatingRadius,0.0]
#    if ingr.packingMode != "gradient":
#        ingr.molarity = 0.0
#        ingr.nbMol = 0
#        print ingr.name
global threed
threed=False


def setCompartment(ingr):
#    ingr.checkCompartment=True
#    ingr.compareCompartment=True#slow down a little.
#    ingr.nbMol*=2
    ingr.rejectionThreshold=100#[1,1,0]#
    if threed : ingr.jitterMax = [1,1,1]
    else : ingr.jitterMax = [1,1,0]
    ingr.cutoff_boundary=0#ingr.encapsulatingRadius/2.0
#    ingr.cutoff_boundary=500+ingr.encapsulatingRadius
    ingr.nbJitter = 6
    #set the nbMol == 30% area total
    print ingr,ingr.jitterMax
    
#   
h.loopThroughIngr(setCompartment)
#setJitter
#raw_input()
if ANALYSIS:
    def clear(h,n=0,torestore=None):
        if not NOGUI :
            afviewer.clearFill("")
        else :
            h.reset()
        gfile = "/Users/ludo/Downloads/autopack_ouput/hiv_experiment/hiv_exp_grid"
        h.buildGrid(boundingBox=h.boundingBox,gridFileIn=gfile,rebuild=True ,
                    gridFileOut=None,previousFill=False)

    def activenmolecules(h,n):
        for i,ingr in enumerate(h.exteriorRecipe.ingredients):
            if i > n-1 :
                ingr.completion=1.0
                ingr.nbMol = 0
            else :
                ingr.completion=0.0

    def setNBMol(h,nmol):
#        if threed use volume and not area
        nt=0
        nb=[]
        na=[]
        weight={1:[1.,],
           2:[0.5,0.5],
            3:[0.37037037037037035,0.37037037037037035,0.18518518518518517+0.07407407407407407],
            4:[0.37037037037037035,0.37037037037037035,0.18518518518518517,0.07407407407407407],
            5:[0.36144578,  0.36144578,  0.18072289,  0.07228915,  0.02409638],
            }
        for ingr in h.exteriorRecipe.ingredients:
            if ingr.completion == 1.0: continue
            if threed :
                ingrarea = (3.0/4.0)*math.pi*ingr.encapsulatingRadius**3.0
            else :
                ingrarea = math.pi*ingr.encapsulatingRadius**2.0
#            nt+=ingr.nbMol
#            nb.append(ingr.nbMol)
            na.append(ingrarea)   
        if threed :
            totalarea = (h.boundingBox[1][0]-h.boundingBox[0][0])**3#(bb[0][1]-bb[0][0])**2            
        else :
            totalarea = (h.boundingBox[1][0]-h.boundingBox[0][0])**2#(bb[0][1]-bb[0][0])**2
        w=weight[nmol]#np.array(nb,float)/nt #weight for each
        na=np.array(na,float)
        tarea=totalarea*0.6#40%
        nb=tarea/sum(na*w)
        c=0
        for i,ingr in enumerate(h.exteriorRecipe.ingredients):
            if ingr.completion == 1.0: continue
            ingr.nbMol = int(w[c]*nb)
            ingr.molarity = 0.0
            print ("WEIGHT ",ingr.name,ingr.nbMol,w[c],nb,tarea,na)
            c+=1
            
    def pack(h,seed,filename,eid):
        h.fill5(seedNum=seed,verbose = 0,usePP=False)
        h.collectResultPerIngredient()
#        r=ingr.results#[pos,rot]
#        ingrpos = numpy.array([ numpy.array(rl) for rl in  numpy.array(r).transpose()[0]])
        #save it 
#        filename="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/gauss_%d_%d_%d.ccp4" % (i,j,k)
#        computeGauss(filename,ingrpos)
#        filename="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/lorentz_%d_%d_%d.ccp4" % (i,j,k)
#        ingrpos=grabSpherePosition(h,ingr)
#        res=ingr.results#[pos,rot]
        #sph_pos = h.compartments[0].surfaceRecipe.ingredients[0].positions[-1]
        todump=[[],[]]
        todump[0] = [np.array(m[0]).tolist() for m in h.molecules]# [[ingr.nbMol,len(res)],[]]  POS     
#        todump[1] = [[tuple(r[0]),r[1].tolist(),"iSUTM",1,0] for r in res ]
#        todump[1]=h.analyse.rdf(np.array(todump[0]),dr=20,rMax=None).tolist()
        with open(filename, 'w') as fp :#doesnt work with symbol link ?
            json.dump(todump,fp)      
#        if threed :
#            plotOneResult3D(h,filename+".png",width=h.boundingBox[1][0])
#        else :
#            plotOneResult2D(h,filename+".png")
        h.store_asJson(resultfilename=filename+"_results.json")
        #rdf ?
#        numpy.savetxt(filename+".txt",h.molecules,delimiter=",")

    def plotOneResult3D(h,filename,width = 1000.0):
        plt.close('all')     # closes the current figure  
        import mpl_toolkits.mplot3d.axes3d as p3
        pos=[]
        s=[]
        c=[]
        for i in range(len(h.molecules)):
            m=h.molecules[i]
            pos.append(np.array(m[0]).tolist())
            s.append(m[2].encapsulatingRadius**2)
            c.append(m[2].color)
        #should be the boundary here ?
        fig = plt.figure()
        bbox = h.boundingBox
        ax = fig.gca(projection='3d')
        x,y,z = np.array(pos).transpose()
        ax.scatter(x, y, z,s=s,c=c)       
        ax.legend()
        ax.set_xlim3d(0, width)
        ax.set_ylim3d(0, width)
        ax.set_zlim3d(0, width)       
        plt.savefig(filename)
        #plt.savefig("test.svg")          
        return x, y, z,s,c
        
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
        plt.close('all')     # closes the current figure        
        
    def one_exp(h,seed,eid=0,nmol=1,periodicity=True):
        clear(h)
#        activenmolecules(h,nmol)
#        setNBMol(h,nmol)
        if periodicity:
            output="/Users/ludo/DEV/autoPACKresults/NM_Analysis_3B_P_"+str(nmol)
            if threed :
                output="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_3D_P_"+str(nmol)            
            h.use_periodicity = True	
            autopack.testPeriodicity = True
            if threed :
                autopack.biasedPeriodicity = [1,1,1]
            else :                
                autopack.biasedPeriodicity = [1,1,0]
        else :
            output="/Users/ludo/DEV/autoPACKresults/NM_Analysis_3B_nP_"+str(nmol)
            if threed :
                output="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_3D_nP_"+str(nmol)            
            h.use_periodicity = False	
            autopack.testPeriodicity = False
            if threed :
                autopack.biasedPeriodicity = [1,1,1]
            else :                
                autopack.biasedPeriodicity = [1,1,0]
#        d = os.path.dirname(output)
        if not os.path.exists(output):
            os.makedirs(output)
        pack(h,seed,output+os.sep+"pos_"+str(eid),eid)
        #plot as 2D figure this pack and save it

    def testPeriodicity(o):
#        h.buildGrid(boundingBox=h.boundingBox,gridFileIn=None,rebuild=True , gridFileOut=None,previousFill=False)
        pos=helper.ToVec(helper.getTranslation(o)) 
        ppos=h.grid.getPositionPeridocity(pos,[1,1,1],100.0)  
        a=helper.Spheres("testP",vertices=ppos,radii=[100.,]*len(ppos))                

    #h.placeMethod="jitter"
    h.placeMethod="pandaBullet"
    h.encapsulatingGrid=0
    h.use_periodicity = True	
    autopack.testPeriodicity = True
    autopack.biasedPeriodicity = [1,1,1]
    analyse = AnalyseAP(env=h, viewer=afviewer, result_file=None)
    h.analyse = analyse
#    output="/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_rev"
    analyse.g.Resolution = 1.0
    h.smallestProteinSize=30.0#get it faster? same result ?
    h.boundingBox=np.array(h.boundingBox)
    seeds_i = analyse.getHaltonUnique(200)
    count=0
#    h.buildGrid(boundingBox=h.boundingBox,gridFileIn=None,rebuild=True , gridFileOut=None,previousFill=False)
    seed = seeds_i[0]
#    count=0
#    one_exp(h,seed,eid=count,nmol=5,periodicity=True)
#    one_exp(h,seed,eid=count,nmol=5,periodicity=False)
#    if not NOGUI :
#        afviewer.displayPreFill();afviewer.displayFill()
#        afviewer.displayFill()

#    for seed in seeds_i :
#        print "seed",seed
#        one_exp(h,seed,eid=count,nmol=1,periodicity=True)
#        one_exp(h,seed,eid=count,nmol=2,periodicity=True)
#        one_exp(h,seed,eid=count,nmol=3,periodicity=True)
#        one_exp(h,seed,eid=count,nmol=4,periodicity=True)
#        one_exp(h,seed,eid=count,nmol=5,periodicity=True)
#        one_exp(h,seed,eid=count,nmol=1,periodicity=False)
#        one_exp(h,seed,eid=count,nmol=2,periodicity=False)
#        one_exp(h,seed,eid=count,nmol=3,periodicity=False)
#        one_exp(h,seed,eid=count,nmol=4,periodicity=False)
#        one_exp(h,seed,eid=count,nmol=5,periodicity=False)
#        count+=1
#    count=0
#    for seed in seeds_i :
#        print "seed",seed
#        one_exp(h,seed,eid=count,nmol=2,periodicity=True)
#        count+=1
#    count=0
#    for seed in seeds_i :
#        print "seed",seed
#        one_exp(h,seed,eid=count,nmol=3,periodicity=True)
#        count+=1
#    count=0
#    for seed in seeds_i[579:] :
#        print "seed",seed
#        one_exp(h,seed,eid=count,nmol=4,periodicity=True)
#        count+=1
#    count=0
#    for seed in seeds_i :
#        print "seed",seed
    seed = seeds_i[0]
    count=0
#    one_exp(h,seed,eid=count,nmol=1,periodicity=False)
#    one_exp(h,seed,eid=count,nmol=2,periodicity=False)
#    one_exp(h,seed,eid=count,nmol=3,periodicity=False)
#    one_exp(h,seed,eid=count,nmol=4,periodicity=False)
#    one_exp(h,seed,eid=count,nmol=5,periodicity=False)
#    one_exp(h,seed,eid=count,nmol=1,periodicity=True)
#    one_exp(h,seed,eid=count,nmol=2,periodicity=True)
#    one_exp(h,seed,eid=count,nmol=3,periodicity=True)
#    one_exp(h,seed,eid=count,nmol=4,periodicity=True)
#    one_exp(h,seed,eid=count,nmol=5,periodicity=True)
        #increase BB 
#        h.boundingBox =[[0., 0., 0.], [2000., 2000., 1.]]
#        one_exp(h,seed,eid=count,nmol=5,periodicity=True)
#        one_exp(h,seed,eid=count,nmol=5,periodicity=False)
#        count+=1
#    count=0
#    for seed in seeds_i :
#        print "seed",seed
#        one_exp(h,seed,eid=count,nmol=2,periodicity=False)
#        count+=1
#    count=0
#    for seed in seeds_i :
#        print "seed",seed
#        one_exp(h,seed,eid=count,nmol=3,periodicity=False)
#        count+=1
#    count=0
#    for seed in seeds_i :
#        print "seed",seed
#        one_exp(h,seed,eid=count,nmol=4,periodicity=False)
#        count+=1

#    fbox_bb=numpy.array(h.boundingBox)
#    h.boundingBox[0]-=numpy.array([100.0,100.0,0.0])    
#    h.boundingBox[1]+=numpy.array([100.0,100.0,0.0])#?
#    nt=0
#    nb=[]
#    na=[]
#    for ingr in h.exteriorRecipe.ingredients:
#        ingrarea = math.pi*ingr.encapsulatingRadius**2.0
#        nt+=ingr.nbMol
#        nb.append(ingr.nbMol)
#        na.append(ingrarea)   
#    totalarea = 1200.0**2#(bb[0][1]-bb[0][0])**2
#    w=np.array(nb,float)/nt #weight for each
#    na=np.array(na,float)
#    tarea=totalarea*0.4#40%
#    nb=tarea/sum(na*w)
#    for i,ingr in enumerate(h.exteriorRecipe.ingredients):
#        ingr.nbMol = int(w[i]*nb)
#        ingr.molarity = 0.0
#    d=analyse.doloop(1,h.boundingBox,wrkDir,output,rdf=True,plot=True,
#                     render=False,twod=TWOD,use_file=True)#,fbox_bb=fbox_bb)
#    if not NOGUI :
#        afviewer.displayFill() 
else :
    gridfile = localdir+os.sep+"autoFillRecipeScripts/Mycoplasma/results/grid_store"
    h.placeMethod="RAPID"
    h.saveResult = True
    h.innerGridMethod = "bhtree"#jordan pure python ? sdf ?
#    h.boundingBox = [[-250.0, -6500.0/2.0, -250.0], [250.0, 6500.0/2.0, 250.0]]
    h.boundingBox =[[-2482, -2389.0, 100.0], [2495, 2466, 2181.0]]
#    h.buildGrid(boundingBox=h.boundingBox,gridFileIn=None,rebuild=True ,
#                          gridFileOut=gridfile,previousFill=False)
    h.buildGrid(boundingBox=h.boundingBox,gridFileIn=gridfile,rebuild=True ,
                          gridFileOut=None,previousFill=False)

    h.fill5(verbose = 0,usePP=True)    
#execfile("/Users/ludo/DEV/autoPACKresults/analyze_B1.py")