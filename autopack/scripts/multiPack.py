# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 13:04:24 2013

@author: ludo
"""
import os
import sys
import pickle
import weakref
from random import randint, random, uniform, gauss, seed
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14 Student_7B992864/plugins/ePMV/mgl64/MGLToolsPckgs")
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14 Student_7B992864/plugins/ePMV/mgl64/MGLToolsPckgs/PIL/")
import AutoFill
AutoFill.usePP = True
#wrkDir = AutoFill.__path__[0]
localdir = wrkDir = AutoFill.__path__[0]

from AutoFill.HistoVol import Environment
from AutoFill.autofill_viewer import AFViewer
from AutoFill.analysis import AnalyseAP

from time import time
import numpy
import pp

import thread
#change AutoPACK if we want multiprocess:
#TypeError: can't pickle weakref objects - remove weakref
#TypeError: can't pickle instancemethod objects - remove / change attribute pointing to function

#self is env
class GrabResult(object):
    def __init__(self,nbFreePoints,histoVol):
        self.queue = []
        self.totalnbFreePoints=nbFreePoints
        self.nbFreePoints=nbFreePoints
        self.distance =[]
        self.freePoints=[]
        self.molecules=[]#result
        self.histoVol=weakref.ref(histoVol)
        self.lock = thread.allocate_lock()
       
    def reset(self,nbFreePoints):
        self.totalnbFreePoints=nbFreePoints
        self.nbFreePoints=nbFreePoints
        self.queue = []
        
    def add(self,value):
        #called at completion of the job
        print ("get ",len(value))
        h=self.histoVol()
#        self.lock.aquire()# and update the grid ?        
        #value is success, ingr,insidePts, newDistancePts , histoVol,organelle.molecules
        ingr = value[1]  
        print ("ingr",ingr.name,ingr.completion,ingr.rejectionCounter)
        memingr = h.getIngrFromName(ingr.name)
#        print ("ingr ",ingr,memingr,value[0])
        if value[0] : 
            if value[1].compNum == 0 :
                organelle = h
            else :
                organelle = h.organelles[abs(value[1].compNum)-1]
#            print ("organelle",organelle)
#            print ("len",len(value[5]))
            #molecules have changed
            if not len(organelle.molecules):
                organelle.molecules = value[5]
            else :
                mol = [x for x in value[5] if x[3] not in numpy.array(organelle.molecules)[:,3]]                       
#               mol = [x for x in value[4] if x not in self.molecules]
                organelle.molecules.extend(mol)
#            print ("omolecules",len(organelle.molecules))
#            print organelle.molecules
            #order[ptInd] changed
#            print ("do we have an hitoVol",value[4])
#            print ("do we have an hitoVol",len(value[4].order))
#            print (value[4].order,h.order)
            h.order = value[4].order
#            print ("order",len(h.order))
            #lastrank changed
            h.lastrank =value[4].lastrank
#            print ("lastrank",h.lastrank)
            #nb_ingredient changed
            h.nb_ingredient =value[4].nb_ingredient
#            print ("nb_ingredient",h.nb_ingredient)
            #update grid  
#            print ("update freePoints and Distance",self.nbFreePoints)
            self.nbFreePoints = h.callFunction(value[1].updateDistances,(h,value[2],
                                                value[3],self.freePoints,self.nbFreePoints, self.distance,
                                                h.masterGridPositions, 0))
#            print ("should be updated",self.nbFreePoints)
#        else :
            #we should update the ingredient in both case
        memingr.counter = ingr.counter
        memingr.completion = ingr.completion
        h.successfullJitter = value[4].successfullJitter
        memingr.rejectionCounter = ingr.rejectionCounter
        memingr.haveBeenRejected = ingr.haveBeenRejected
#            diff = self.totalnbFreePoints - value[1]
#            self.nbFreePoints-=diff
#            self.distance = value[3]
#            self.freePoints= value[2]
            
        self.queue.append([value[0]])#value is success,nbFreepoint,freePoints, distance
        del value
#        self.lock.release()
       # print ("queue extended")

def pack_multi(env, ncpus=1,seedNum=14, stepByStep=False, verbose=False, sphGeom=None,
              labDistGeom=None, debugFunc=None,name = None, vTestid = 3,vAnalysis = 0,**kw):
        ## Fill the grid by picking an ingredient first and then
        ## this filling should be able to continue from a previous one
        ## find a suitable point suing hte ingredient's placer object
        self=env
        import time
        import pp
        nparts = 2
        job_server = pp.Server(ncpus=ncpus)
        t1=time.time()
        self.timeUpDistLoopTotal = 0 #Graham added to try to make universal "global variable Verbose" on Aug 28
        self.static=[]
        if self.grid is None:
            print("no grid setup")
            return
        # create a list of active ingredients indices in all recipes to allow
        # removing inactive ingredients when molarity is reached
        allIngredients = self.callFunction(self.getActiveIng)

        nbIngredients = len(allIngredients)
        self.cFill = self.nFill
        if name == None :
            name = "F"+str(self.nFill)
        self.FillName.append(name)
        self.nFill+=1
        # seed random number generator
        SEED=seedNum
        numpy.random.seed(SEED)#for gradient
        seed(seedNum)
        self.randomRot.setSeed(seed=seedNum)
        # create copies of the distance array as they change when molecules
        # are added, theses array can be restored/saved before feeling
        freePoints = self.grid.freePoints[:]
        nbFreePoints = len(freePoints)#-1
        grab_callback = GrabResult(nbFreePoints,self)
        grab_callback.freePoints=freePoints
        
#        self.freePointMask = numpy.ones(nbFreePoints,dtype="int32")
        if "fbox" in kw :  # Oct 20, 2012  This is part of the code that is breaking the grids for all meshless organelle fills
            self.fbox = kw["fbox"]
        if self.fbox is not None and not self.EnviroOnly :
            self.freePointMask = numpy.ones(nbFreePoints,dtype="int32")
            bb_insidepoint = self.grid.getPointsInCube(self.fbox, [0,0,0], 1.0)[:]#center and radius ?3,runTime=self.runTimeDisplay
            self.freePointMask[bb_insidepoint]=0
            bb_outside = numpy.nonzero(self.freePointMask)
            self.grid.gridPtId[bb_outside] = 99999
        compId = self.grid.gridPtId
        #why a copy? --> can we split ?
        distance = self.grid.distToClosestSurf[:]
        grab_callback.distance = distance
        spacing = self.smallestProteinSize

        # DEBUG stuff, should be removed later
        self.jitterVectors = []
        self.jitterLength = 0.0
        self.totnbJitter = 0
        self.maxColl = 0.0
        self.successfullJitter = []
        self.failedJitter = []
        
        #this function also depend on the ingr.completiion that can be restored ?
        self.activeIngr0, self.activeIngr12 = self.callFunction(self.getSortedActiveIngredients, (allIngredients,verbose))

        print('len(allIngredients', len(allIngredients))
        print('len(self.activeIngr0)', len(self.activeIngr0))
        print('len(self.activeIngr12)', len(self.activeIngr12))
        self.activeIngre_saved = self.activeIngr[:]

        self.totalPriorities = 0 # 0.00001
        for priors in self.activeIngr12:
            pp = priors.packingPriority
            self.totalPriorities = self.totalPriorities + pp
            print('totalPriorities = ', self.totalPriorities)
        previousThresh = 0
        self.normalizedPriorities = []
        self.thresholdPriorities = [] 
        # Graham- Once negatives are used, if picked random# 
        # is below a number in this list, that item becomes 
        # the active ingredient in the while loop below
        for priors in self.activeIngr0:
            self.normalizedPriorities.append(0)
            if self.pickWeightedIngr :#why ?
                self.thresholdPriorities.append(2)
        for priors in self.activeIngr12:
            #pp1 = 0
            pp = priors.packingPriority
            if self.totalPriorities != 0:
                np = float(pp)/float(self.totalPriorities)
            else:
                np=0.
            self.normalizedPriorities.append(np)
            print('np is ', np, ' pp is ', pp, ' tp is ', np + previousThresh)
            self.thresholdPriorities.append(np + previousThresh)
            previousThresh = np + float(previousThresh)
        self.activeIngr = self.activeIngr0 + self.activeIngr12

        nls=0
        totalNumMols = 0
        self.totalNbIngr = self.getTotalNbObject(allIngredients)
        if len(self.thresholdPriorities ) == 0:
            for ingr in allIngredients:
                totalNumMols += ingr.nbMol
            print('totalNumMols Fill5if = ', totalNumMols)
        else :                
            for threshProb in self.thresholdPriorities:
                nameMe = self.activeIngr[nls]
                print('threshprop Fill5else is %f for ingredient: %s %s %d'%(threshProb, nameMe,nameMe.name,nameMe.nbMol))
                totalNumMols += nameMe.nbMol
                print('totalNumMols Fill5else = ', totalNumMols)
                nls+=1
        print ("tobj = ",self.totalNbIngr) 
        a=numpy.ones((self.totalNbIngr,3))*999999999.9#*100.0#max ingredient excepted.

        vRangeStart = 0.0
        tCancelPrev=time.time()
        test = True
        kk=0
        ptInd = 0

        PlacedMols = 0
        vThreshStart = 0.0   # Added back by Graham on July 5, 2012 from Sept 25, 2011 thesis version
        
        #if bullet build the organel rbnode
        if self.placeMethod == "pandaBullet":
            self.setupPanda()
            for o in self.organelles:
                if o.rbnode is None :
                    o.rbnode = self.addMeshRBOrganelle(o)
#==============================================================================
#         #the big loop
#==============================================================================
        #self.largestProteinSize = 0 #before starting reset largest size
        self.grabedvalue=[]
        while nbFreePoints:
            #breakin test
            print ("nbFreePoints",nbFreePoints,len(self.activeIngr),vRangeStart)
            if len(self.activeIngr)==0:
                print('broken by len****')
                break
            if vRangeStart>1:
                print('broken by vRange and hence Done!!!****')
                break   
            #we do two pass by dividing the grid on Y by 2*ncpus
            for p in range(nparts):  
                print ("pass n ",p)
                grab_callback.reset(grab_callback.nbFreePoints)
                picked_ingredients=[]
                for n in range(ncpus):
                    print ("prepare job",n)
                    ingr =  self.callFunction(self.pickIngredient,(vThreshStart,))
                    picked_ingredients.append(ingr)
                    compNum = ingr.compNum
                    radius = ingr.minRadius
                    jitter = self.callFunction(ingr.getMaxJitter,(spacing,))
    
                    # compute dpad which is the distance at which we need to update
                    # distances after the drop is successfull
                    mr = self.get_dpad(compNum)
                    dpad = ingr.minRadius + mr + jitter
                
                    ## find the points that can be used for this ingredients
                    ## in the slice ? 
                    res=self.callFunction(self.getPointToDrop,(ingr,radius,jitter,
                                                freePoints,nbFreePoints,
                                                distance,compId,compNum,vRangeStart,vThreshStart))
        #                                        distance,compId,compNum,vRangeStart))   # Replaced this with Sept 25, 2011 thesis version on July 5, 2012
                    print ("pick",ingr,res)
                    if res[0] :
                        ptInd = res[1]
                        if ptInd > len(distance):
                            print ("problem ",ptInd)
                            continue
                    else :
                        print ("vRangeStart coninue ",res)
                        vRangeStart = res[1]
                        continue
    #                continue
        #            print ("picked ",ptInd)
                    #place the ingrediant
                    if self.overwritePlaceMethod :
                        ingr.placeType = self.placeMethod
                    #check the largestProteinSize
                    #worker ?
                    if ingr.encapsulatingRadius > self.largestProteinSize : 
                        self.largestProteinSize = ingr.encapsulatingRadius
                    print ("submit job",n)
                    p=job_server.submit(ingr.place_mp,(self, ptInd, 
                                        freePoints, nbFreePoints, distance, dpad,True,
                                        stepByStep, verbose), modules=("AutoFill",),
                                        callback=grab_callback.add)
    #                success, nbFreePoints = self.callFunction(ingr.place,(self, ptInd, 
    #                                    freePoints, nbFreePoints, distance, dpad,
    #                                    stepByStep, verbose),
    #                                    {"debugFunc":debugFunc})
                print ("after place nbFreePoints",nbFreePoints)
                job_server.wait() 
                #grab result
                results=grab_callback.queue[:]#success, nbFreePoints
#                print ("results",results)               
                nbFreePoints=grab_callback.nbFreePoints
                print ("cumul freePts",grab_callback.nbFreePoints)
                #need to cumul freepoints and distance
                distance = grab_callback.distance[:]
                freePoints = grab_callback.freePoints[:] 
#                self.molecules = grab_callback.molecules[:]
#                self.grabedvalue.append(grab_callback.queue[:]) 
                for n in range(ncpus):                
        #            print("nbFreePoints after PLACE ",nbFreePoints)
                    ingr = picked_ingredients[n]
                    success = results[n][0]
#                    nbFreePoints = results[n][1]
                    if success:
                        if ingr.encapsulatingRadius > self.largestProteinSize : 
                            self.largestProteinSize = ingr.encapsulatingRadius
                        PlacedMols+=1
                    print ("ingr",ingr.completion,ingr.name)
                    if ingr.completion >= 1.0 :
                        ind = self.activeIngr.index(ingr)
                        if ind > 0:
                            #j = 0
                            for j in range(ind):   
                                if j >= len(self.thresholdPriorities) or j >= len(self.normalizedPriorities):
                                    continue
                                self.thresholdPriorities[j] = self.thresholdPriorities[j] + self.normalizedPriorities[ind]
                        self.activeIngr.pop(ind)
                        self.activeIngr0, self.activeIngr12 = self.callFunction(self.getSortedActiveIngredients, (self.activeIngr,verbose))
                        self.activeIngre_saved = self.activeIngr[:]
        
                        self.totalPriorities = 0 # 0.00001
                        for priors in self.activeIngr12:
                            pp = priors.packingPriority
                            self.totalPriorities = self.totalPriorities + pp
        
                        previousThresh = 0
                        self.normalizedPriorities = []
                        self.thresholdPriorities = [] 
                        # Graham- Once negatives are used, if picked random# 
                        # is below a number in this list, that item becomes 
                        #the active ingredient in the while loop below
                        for priors in self.activeIngr0:
                            self.normalizedPriorities.append(0)
                            if self.pickWeightedIngr :
                                self.thresholdPriorities.append(2)
                        for priors in self.activeIngr12:
                            #pp1 = 0
                            pp = priors.packingPriority
                            if self.totalPriorities != 0:
                                np = float(pp)/float(self.totalPriorities)
                            else:
                                np=0.
                            self.normalizedPriorities.append(np)
                            self.thresholdPriorities.append(np + previousThresh)
                            previousThresh = np + float(previousThresh)
                        self.activeIngr = self.activeIngr0 + self.activeIngr12
            #end for part
#            break
        self.distancesAfterFill = distance
        self.freePointsAfterFill = freePoints
        self.nbFreePointsAfterFill = nbFreePoints
        self.distanceAfterFill = distance

        t2 = time.time()
        print('time to fill', t2-t1)
            
        if self.saveResult:
            self.grid.freePoints = freePoints[:]
            self.grid.distToClosestSurf = distance[:]
            #shoul check extension filename for type of saved file
            self.saveGridToFile(self.resultfile+"grid")
            self.grid.result_filename = self.resultfile+"grid"
            self.store()
            self.store_asTxt()
            self.store_asJson()            
        print('time to save end', time.time()-t2)            
        ingredients ={}
        for pos, rot, ingr, ptInd in self.molecules:
            if ingr.name not  in ingredients :
                ingredients[ingr.name]=[ingr,[],[],[]]
            mat = rot.copy()
            mat[:3, 3] = pos
            ingredients[ingr.name][1].append(pos)
            ingredients[ingr.name][2].append(rot)
            ingredients[ingr.name][3].append(numpy.array(mat))
        for o in self.organelles:
            for pos, rot, ingr, ptInd in o.molecules:
                if ingr.name not  in ingredients :
                    ingredients[ingr.name]=[ingr,[],[],[]]
                mat = rot.copy()
                mat[:3, 3] = pos
                ingredients[ingr.name][1].append(pos)
                ingredients[ingr.name][2].append(rot)
                ingredients[ingr.name][3].append(numpy.array(mat)) 
        self.ingr_result = ingredients
        if self.treemode == "bhtree" :
            from bhtree import bhtreelib
            bhtreelib.freeBHtree(self.close_ingr_bhtree)
#        bhtreelib.FreeRBHTree(self.close_ingr_bhtree)
#        del self.close_ingr_bhtree
                    
                    
def pack(env,seed=20,forceBuild=True,vTestid = 3,vAnalysis = 0):
    t1 = time()
    env.fill5(seedNum=seed,verbose=4, vTestid = vTestid,vAnalysis = vAnalysis)
    t2 = time()
    print('time to run Fill5', t2-t1)

def grid_pack(bb,env,wrkDir,forceBuild=True, fill=0,seed=20,vTestid = 3,vAnalysis = 0):
    t1 = time()
#        if bbox is None :
#            box=self.helper.getCurrentSelection()[0]
#        else :
#            box = bbox[0]
#        bb=self.helper.getCornerPointCube(box)
    gridFileIn=None
    gridFileOut=None
    if forceBuild :
        gridFileOut=wrkDir+"fill_grid"
    else :
        gridFileIn=wrkDir+"fill_grid"
    env.buildGrid(boundingBox=bb,gridFileIn=None,rebuild=True ,
                      gridFileOut=None,previousFill=False)
#    h.buildGrid(gridFileIn=gridFileIn, 
#                  gridFileOut=gridFileOut)
    t2 = time()
    gridTime = t2-t1
    if fill :
        pack(seed=seed,vTestid = vTestid,vAnalysis = vAnalysis)
    print ('time to Build Grid', gridTime)
#    afviewer.displayOrganellesPoints()
    #return
    #actine.updateFromBB(h.grid)


TWOD = 1
NOGUI = 1
helper = AutoFill.helper
if helper is None and not NOGUI:
    import upy
    helperClass = upy.getHelperClass()
    helper =helperClass()
else :
    import upy
    helperClass = upy.getHelperClass()
    helper =helperClass(vi="nogui")
AutoFill.helper = helper
filename = "/Users/ludo/DEV/autofill_svn/trunk/AutoFillClean/autoFillRecipeScripts/2DsphereFill/Test_Spheres2D1.0.xml"
fileName, fileExtension = os.path.splitext(filename)
n=os.path.basename(fileName)
h = Environment(name=n)
recipe=n
h.load_XML(filename)
afviewer=None
if not NOGUI :
    setattr(h,"helper",helper)
    afviewer = AFViewer(ViewerType=h.helper.host,helper=h.helper)
    afviewer.SetHistoVol(h,20.0,display=False)
    h.host=h.helper.host
    afviewer.displayPreFill()
h.saveResult = True
ingr=h.exteriorRecipe.ingredients[0]
resultfilename = h.resultfile = wrkDir+os.sep+"autoFillRecipeScripts/2DsphereFill/results/Test.afr"
h.buildGrid(boundingBox=h.boundingBox,gridFileIn=None,rebuild=True ,
                      gridFileOut=None,previousFill=False)
pack_multi(h,ncpus=1)
if not NOGUI :
    afviewer.displayFill()  

#load ?
#result,orgaresult,freePoint=h.load(resultfilename=resultfilename,restore_grid=False)#load text ?#this will restore the grid  
#ingredients = h.restore(result,orgaresult,freePoint)
               
#filename = "/Users/ludo/Desktop/cell.xml"
#filename = "/Users/ludo/DEV/autofill_svn/trunk/AutoFillClean/autoFillRecipeScripts/2DsphereFill/Test_Spheres2D1.0.xml"
#execfile("/Users/ludo/DEV/autofill_svn/trunk/AutoFillClean/multiPack")