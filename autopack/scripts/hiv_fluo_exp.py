# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:59:01 2013

@author: ludo
Hoe to run autopack from xml file
"""
# Show the effect of garbage collection
import sys
import numpy
from numpy import matrix
import json
import math
from random import randrange,seed,random
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14_95696CEA/plugins/ePMV/mgl64/MGLToolsPckgs/")
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
    #meshFile="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/MA_hexagone.c4d" 
#                rep="EnvelopeBilayer" rep_file="http://autofill.googlecode.com/svn/data/HIV/geoms/EnvelopeBilayer3"
import os
import sys
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14_95696CEA/plugins/ePMV/mgl64/MGLToolsPckgs")
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14_95696CEA/plugins/ePMV/mgl64/MGLToolsPckgs/PIL/")
import autopack
#wrkDir = AutoFill.__path__[0]
localdir = wrkDir = autopack.__path__[0]

from Volume.IO.volWriters import WriteCCP4
from Volume.Grid3D import Grid3DUC, Grid3DSI, Grid3DF
import numpy
import scipy
from scipy import ndimage 

from autopack.Environment import Environment
from autopack.Graphics import AutopackViewer as AFViewer
from autopack.Analysis import AnalyseAP
from autopack import transformation as tr 
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
autopack.helper = helper
import numpy as np
from scipy import fftpack
from scipy import stats

def convolve(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft*psf_fft)))

def deconvolve(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft/psf_fft)))

def convoluteGauss(filename,ingpos,boundingBox=[[-693.41299629211426, -651.68601989746094, -735.077], [653.65499687194824, 664.00502014160156, 694.65200000000004]]):
    w=WriteCCP4()
    numpy.savetxt(filename+".txt",ingpos,delimiter=",")
    space=50.0
    boundingBox=2.0*numpy.array(boundingBox)
    #inpos need to be translated along the normal in order to get position of FAB+fluorochrome
    #can use an extrasphere in the ingredient that could represent it ?
#    center = (boundingBox[1]-boundingBox[0])/2.0
    x = numpy.arange(boundingBox[0][0], boundingBox[1][0] + space, space)#*1.1547) gridspacing is already multiplied by 1.1547
    y = numpy.arange(boundingBox[0][1], boundingBox[1][1] + space, space)#*1.1547)
    z = numpy.arange(boundingBox[0][2], boundingBox[1][2] + space, space)#*1.1547)
    nx = len(x) # sizes must be +1 or the right, top, and back edges don't get any points using this numpy.arange method
    ny = len(y)
    nz = len(z)
    zvalues = numpy.zeros(nx*ny*nz)
    values = zvalues.reshape((nx,ny,nz))
    ax,ay,az = ingpos.transpose()
    ix = np.searchsorted(x, ax) - 1
    iy = np.searchsorted(y, ay) - 1
    iz = np.searchsorted(z, az) - 1
    gaussian_blur_sigma =  (400.0/space)/2.35482#16.986436330590024/10.0   0.42lambdan or 1.222lamdan
    cauchy_sigma = (400.0/space)/2.0
    #square or power 
    values[ix,iy,iz]=600.0
#    gauss=ndimage.filters.gaussian_filter(values,gaussian_blur_sigma)
    X, Y , Z= np.ogrid[-nx/2:nx/2, -ny/2:ny/2, -nz/2:nz/2]#-27:27
    psf_gauss = stats.norm.pdf(np.sqrt((X*space)**2 + (Y*space)**2 +(Z*space)**2), 0, gaussian_blur_sigma*space)#norm is gaussian
#    psf_cauchy = stats.cauchy.pdf(np.sqrt((X*space)**2 + (Y*space)**2 +(Z*space)**2), 0, cauchy_sigma*space)#norm is gaussian
    gauss_conv = convolve(values,psf_gauss)
#    cauchy_conv = convolve(values,psf_cauchy)
    h={}
    maskGrid = Grid3DF( numpy.array(np.real(gauss_conv),'f'), [0,0,0], [space,space,space] , h)
    maskGrid.stepSize = [space,space,space]
    maskGrid.origin = (boundingBox[1]-boundingBox[0])/2.0
    #maskGrid.normalize()
    h['amin'], h['amax'],h['amean'],h['arms']= maskGrid.stats()
    w.write(filename,maskGrid)

def computeGauss(filename,ingpos,boundingBox=[[-693.41299629211426, -651.68601989746094, -735.077], [653.65499687194824, 664.00502014160156, 694.65200000000004]]):
    w=WriteCCP4()
    numpy.savetxt(filename+".txt",ingpos,delimiter=",")
    space=200.0
    boundingBox=2.0*numpy.array(boundingBox)
    #inpos need to be translated along the normal in order to get position of FAB+fluorochrome
    #can use an extrasphere in the ingredient that could represent it ?
#    center = (boundingBox[1]-boundingBox[0])/2.0
    x = numpy.arange(boundingBox[0][0], boundingBox[1][0] + space, space)#*1.1547) gridspacing is already multiplied by 1.1547
    y = numpy.arange(boundingBox[0][1], boundingBox[1][1] + space, space)#*1.1547)
    z = numpy.arange(boundingBox[0][2], boundingBox[1][2] + space, space)#*1.1547)
    nx = len(x) # sizes must be +1 or the right, top, and back edges don't get any points using this numpy.arange method
    ny = len(y)
    nz = len(z)
    zvalues = numpy.zeros(nx*ny*nz)
    values = zvalues.reshape((nx,ny,nz))
    ax,ay,az = ingpos.transpose()
    ix = np.searchsorted(x, ax) - 1
    iy = np.searchsorted(y, ay) - 1
    iz = np.searchsorted(z, az) - 1
    gaussian_blur_sigma =  (400.0/space)/2.35482#16.986436330590024/10.0
    values[ix,iy,iz]=1
    #DeltaX =DeltaY= lambda / 2*n*sin(alpha)
    #DeltaZ = lambda / 2*n*sin(alpha)**2
    #lambda is the wavelength
    #n is indice of refraction air = 1, water is 1.333
    #alpha is hald the aperture ie angle microscope
    gauss=ndimage.filters.gaussian_filter(values,gaussian_blur_sigma)
    h = {}
    #span ?
    #gaussian noise
    #noise = numpy.random.normal(0.5,0.5, gauss.shape)
    noise=numpy.random.uniform(0,1,gauss.shape)
    g=gauss#*noise
    maskGrid = Grid3DF( numpy.array(g,'f'), [0,0,0], [space,space,space] , h)
    maskGrid.stepSize = [space,space,space]
    maskGrid.origin = (boundingBox[1]-boundingBox[0])/2.0
    #maskGrid.normalize()
    h['amin'], h['amax'],h['amean'],h['arms']= maskGrid.stats()
    w.write(filename,maskGrid)

#filename = "/Users/ludo/Desktop/cell.xml"
#filename = wrkDir+os.sep+"autoFillRecipeScripts/2DsphereFill/Test_Spheres2D1.1.xml"
#filename = wrkDir+os.sep+"autoFillRecipeScripts/2DsphereFill/Test_Spheres2Dgradients1.0.xml"
#filename = "/Users/ludo/DEV/autopack_git/data/Mycoplasma/recipe/Mycoplasma1.3.xml"
#filename = "/Users/ludo/Desktop/cell_hack.xml"
filename = "/Users/ludo/Downloads/autopack_ouput/hiv_experiment/HIV1.6.xml"
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
    #afviewer.displayPreFill()
h.saveResult = True

#resultfilename = h1.resultfile = wrkDir+os.sep+"autoFillRecipeScripts"+os.sep+"2DsphereFill"+os.sep+"results"+os.sep+"SpherefillResult.afr"
#resultfilename = h.resultfile = wrkDir+os.sep+"autoFillRecipeScripts/2DsphereFill/results/2DsphereFill_1.1.apr"
#resultfilename = h.resultfile = wrkDir+os.sep+"autoFillRecipeScripts/Mycoplasma/results/MycoplasmaPackResult_3"
resultfilename = h.resultfile = "/Users/ludo/Downloads/autopack_ouput/hiv_experiment/hiv_exp"

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
def setCompartment(ingr):
#    ingr.checkCompartment=True
#    ingr.compareCompartment=True#slow down a little.
#    ingr.nbMol*=2
    ingr.rejectionThreshold=100#[1,1,0]#
#    ingr.jitterMax =[ingr.encapsulatingRadius/(25.*1.1547),ingr.encapsulatingRadius/(25.*1.1547),0.0]
#    ingr.cutoff_boundary=ingr.encapsulatingRadius
#    ingr.cutoff_boundary=500+ingr.encapsulatingRadius
    ingr.nbJitter = 6

    
h.loopThroughIngr(setCompartment)
#setJitter
#raw_input()
if ANALYSIS:
    h.placeMethod="RAPID"
    h.encapsulatingGrid=0
    autopack.testPeriodicity = False
    analyse = AnalyseAP(env=h, viewer=afviewer, result_file=None)
    output="/Users/ludo/DEV/autoPACKresults/NM_Analysis_C1"
    analyse.g.Resolution = 1.0
    h.boundingBox=numpy.array(h.boundingBox)
    fbox_bb=numpy.array(h.boundingBox)
#    h.boundingBox[0]-=numpy.array([500.0,500.0,0.0])    
#    h.boundingBox[1]+=numpy.array([500.0,500.0,0.0])
    d=analyse.doloop(5,h.boundingBox,wrkDir,output,rdf=True,
                     render=False,twod=TWOD,use_file=True)#,fbox_bb=fbox_bb)
#    if not NOGUI :
#        afviewer.displayFill() 
else :
    global fluo_pos
    fluo_pos=numpy.array([55.,190., -30.0])
    #grab previous result
    #organelle.molecules
    def grabSpherePosition(h,ingr):
        h.collectResultPerIngredient()
        res=ingr.results#[pos,rot]
        #sph_pos = h.compartments[0].surfaceRecipe.ingredients[0].positions[-1]
        pos = [ingr.transformPoints(r[0],r[1], [fluo_pos])[0] for r in res ]
        return pos

    def render():
        #render each axes 
        pass

    def displayResult(h,i,j,k):
        if not NOGUI :
            afviewer.clearFill("HIV")
        else :
            h.reset()
        with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/ma_models/model_%d.json"%i, 'r') as fp :#doesnt work with symbol link ?
            ma_results=json.load(fp)
        #treat ma_reslt for the issing 90degree rotation
        mat = tr.rotation_matrix(math.pi/2.0,[0,0,1]).transpose()
        for ma in ma_results :
            ma[1] = numpy.array(matrix(ma[1])*matrix(mat)).tolist()
        with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/env_models/model_%d_%d_%d.json"%(i,j,k), 'r') as fp :#doesnt work with symbol link ?
            env_results=json.load(fp)
#        mat = tr.rotation_matrix(math.pi/2.0,[0,0,1]).transpose()
#        mat = tr.rotation_matrix(-math.pi/2.0,[0,1,0]).transpose()
#        for env in env_results[1] :
#            env[1] = numpy.array(matrix(env[1])*matrix(mat)).tolist()
        #can also use the load_result function
        #filename="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/lorentz_%d_%d_%d.ccp4" % (i,j,k)
        #a=numpy.loadtxt(filename+".txt",delimiter=",")#position of ENV
        #env_results=[ [p,numpy.identity(4),"iSUTM",1,0] for p in a]
        ma_results.extend(env_results[1])
        ingredients = h.restore([],[ma_results],[])
        if not NOGUI :
            afviewer.doPoints = False #self.getVal(self.points_display)
            afviewer.doSpheres = True
            afviewer.quality = 1 #lowest quality for sphere and cylinder
            afviewer.visibleMesh = True #mesh default visibility 
            afviewer.displayPreFill()
            afviewer.displayFill()
        
    def clear(h,n=0,torestore=None):
        if not NOGUI :
            afviewer.clearFill("HIV")
        else :
            h.reset()
        if torestore is not None :
            if len(torestore):
                if n!= 0:
                    ingredients = h.restore([],[torestore[:n]],[])#for testing only use the first three.
                else :
                    ingredients = h.restore([],[torestore],[])#for testing only use the first three. 
            gfile = "/Users/ludo/Downloads/autopack_ouput/hiv_experiment/hiv_exp_grid"
            h.buildGrid(boundingBox=h.boundingBox,gridFileIn=gfile,rebuild=False ,
                        gridFileOut=None,previousFill=True)
            ingr2 = h.getIngrFromName("Envelope_surf__MA_hexagone")
            ingr2.completion=1.0
            ingr2.nbMol = 0
            

    def pack(h,sd,i,j,k,maresults):#model number / ENV mode nb / ENV distribution nb
        #k,i,eid,k        
        clear(h,n=0,torestore=maresults)  
        print "############################################"
        print i,j,k
        print "############################################"
        if not NOGUI:
            autopack.helper.update()
        #choose randomly nb of ENV between 7-12
        seed(int(i*100*6 + j*100 + k))#unique seed
        ingr = h.getIngrFromName("Envelope_surf__iSUTM")
        ingr.nbMol = randrange(5,15)#work ?
        ingr.overwrite_distFunc=False
        h.resultfile = "/Users/ludo/Downloads/autopack_ouput/hiv_experiment/res/hiv_exp_%d_%d_%d"% (i,j,k)
        h.fill5(seedNum=int(k*10*8 + j*10 + i),verbose = 0,usePP=False)
        h.collectResultPerIngredient()
        r=ingr.results#[pos,rot]
#        ingrpos = numpy.array([ numpy.array(rl) for rl in  numpy.array(r).transpose()[0]])
        #save it 
#        filename="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/gauss_%d_%d_%d.ccp4" % (i,j,k)
#        computeGauss(filename,ingrpos)
        filename="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/lorentz_%d_%d_%d.ccp4" % (i,j,k)
        ingrpos=grabSpherePosition(h,ingr)
        res=ingr.results#[pos,rot]
        #sph_pos = h.compartments[0].surfaceRecipe.ingredients[0].positions[-1]
        todump = [[ingr.nbMol,len(res)],[]]       
        todump[1] = [[tuple(r[0]),r[1].tolist(),"iSUTM",1,0] for r in res ]
        with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/env_models/model_%d_%d_%d.json"%(i,j,k), 'w') as fp :#doesnt work with symbol link ?
            json.dump(todump,fp)        
        numpy.savetxt(filename+".txt",ingrpos,delimiter=",")
        #convoluteGauss(filename,ingrpos)
        #save
        
    def modeRandom(h,n,i,maresults,eid=0):
        ingr = h.getIngrFromName("Envelope_surf__iSUTM")
        ingr.partners={}
        ingr.partners_name=[]
        ingr.packingMode="random"
        h.ingrLookForNeighbours = False 
        for k in range(n):
            pack(h,k,i,eid,k,[])
                
    def modeRandomEmpty(h,n,i,maresults,eid=1):
        ingr = h.getIngrFromName("Envelope_surf__iSUTM")
        ingr.partners={}
        ingr.partners_name=[]
        ingr.packingMode="random"
        h.ingrLookForNeighbours = False 
        for k in range(n):
            pack(h,k,i,eid,k,maresults)
            
    def modeCloseMA(h,n,i,maresults,eid=2):
        ingr = h.getIngrFromName("Envelope_surf__iSUTM")
        ingr.partners={}
        ingr.partners_name=["MA_hexagone"]
        ingr.packingMode="closePartner"                
        ingr.proba_not_binding = 0#0.5 ?
        ingr.partners_position=[[],]
        ingr2 = h.getIngrFromName("Envelope_surf__MA_hexagone")
        ingr2.weight=1.0        
        h.set_partners_ingredient(ingr)
        h.ingrLookForNeighbours = True 
        for k in range(n):
            pack(h,k,i,eid,k,maresults)

    def modeCloseMA_ENV(h,n,i,maresults,eid=3,w1=0.5,w2=0.5,w3=0):
        ingr = h.getIngrFromName("Envelope_surf__iSUTM")
        ingr.partners={}
        ingr.partners_name=["MA_hexagone","iSUTM"]
        ingr.packingMode="closePartner"                
        ingr.proba_not_binding = w3#0.5 ?
        ingr.partners_position=[[],[]]
        ingr.weight=w1
        ingr2 = h.getIngrFromName("Envelope_surf__MA_hexagone")
        ingr2.weight=w2 
        h.set_partners_ingredient(ingr)
        h.ingrLookForNeighbours = True 
        for k in range(n):
            pack(h,k,i,eid,k,maresults)

    def modeClosePlane(h,n,i,maresults,eid=3,w1=0.5,w2=0.5,w3=0):
        ingr = h.getIngrFromName("Envelope_surf__iSUTM")
        ingr.partners={}
        ingr.partners_name=["MA_hexagone","iSUTM"]
        ingr.packingMode="closePartner"                
        ingr.proba_not_binding = w3#0.5 ?
        ingr.partners_position=[[],[]]
        ingr.weight=w1
        ingr2 = h.getIngrFromName("Envelope_surf__MA_hexagone")
        ingr2.weight=w2 
        h.set_partners_ingredient(ingr)
        h.ingrLookForNeighbours = True 
        for k in range(n):
            pack(h,k,i,eid,k,maresults)

    def modeCloseENV(h,n,i,maresults,eid=4):
        ingr = h.getIngrFromName("Envelope_surf__iSUTM")
        ingr.partners={}
        ingr.partners_name=["iSUTM"]
        ingr.packingMode="closePartner"                
        ingr.proba_not_binding = 0.0#0.5 ?
        ingr.partners_position=[[],]
        ingr.weight=1.0
        h.set_partners_ingredient(ingr)
        h.ingrLookForNeighbours = True 
        for k in range(n):
            pack(h,k,i,eid,k,maresults)

    def random_subset( iterator, K ):
        #from http://stackoverflow.com/questions/2612648/reservoir-sampling
        result = []
        N = 0   
        for item in iterator:
            N += 1
            if len( result ) < K:
                result.append( item )
            else:
                s = int(random() * N)
                if s < K:
                    result[ s ] = item   
        return result
        
    def pickRandomFromMApos(h,sd,i,j,k,maresults):
        #there is a missing rotation of 90degree around Z
        ingr = h.getIngrFromName("Envelope_surf__iSUTM")
        n=len(maresults)
        mat = tr.rotation_matrix(-math.pi/2.0,[1,0,0]).transpose()
        imat= []
        for ma in maresults :
            imat.append( [ma[0],numpy.array(matrix(ma[1])*matrix(mat)).tolist() ])
        seed(int(i*100*6 + j*100 + k))#seed(sd)
        nbMol = randrange(5,15)
        res=random_subset(imat,nbMol)
        ingrpos = [ingr.transformPoints(r[0],r[1], [fluo_pos])[0] for r in res ]
        #pick nbMol random position from maresults
        filename="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/lorentz_%d_%d_%d.ccp4" % (i,j,k)
        numpy.savetxt(filename+".txt",ingrpos,delimiter=",")
#        res=[]
        todump = [[nbMol,len(res)],[]]
        todump[1] = [[r[0],r[1],"iSUTM",1,0] for r in res ]
        with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/env_models/model_%d_%d_%d.json"%(i,j,k), 'w') as fp :#doesnt work with symbol link ?
            json.dump(todump,fp)        
        #with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/env_models/model_%d.json"%m, 'w') as fp :#doesnt work with symbol link ?
        #    json.dump(fp,res)        
        #convoluteGauss(filename,ingrpos)
                
    def modeMA_hex_center(h,n,i,maresults,eid=5):
        #use a mask for position OR pick randomly n number of hexagon position
#        h.compartments[0].surfaceRecipe.ingredients[0].partners={}
#        h.compartments[0].surfaceRecipe.ingredients[1].packingMode="random"
#        h.ingrLookForNeighbours = False 
        for k in range(n):
            pickRandomFromMApos(h,k,i,eid,k,maresults)
            
packing = False
disp = False
if not NOGUI :
    helper.update()
if packing :    
#    h.use_halton = True          
    autopack.verbose = 0
    h.runTimeDisplay=0
    h.saveResult = False
    for m in range(10):
        with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/ma_models/model_%d.json"%m, 'r') as fp :#doesnt work with symbol link ?
            ma_results=json.load(fp)
        ma_copy = ma_results[:]    
        #treat ma_reslt for the issing 90degree rotation
        mat = tr.rotation_matrix(math.pi/2.0,[0,0,1]).transpose()
        for ma in ma_results :
            ma[1] = numpy.array(matrix(ma[1])*matrix(mat)).tolist()
#        mat = tr.rotation_matrix(math.pi/2.0,[1,0,0])
#        for ma in ma_results :
#            ma[1] = numpy.array(matrix(mat)*matrix(ma[1])).tolist()
#        modeRandom(h,100,m,ma_results,eid=0)
#        modeRandomEmpty(h,100,m,ma_results,eid=1)
        modeMA_hex_center(h,100,m,ma_copy,eid=2)
#        modeCloseMA(h,100,m,ma_results,eid=3) #close to edge but need to increase de speed of packing
#        modeCloseMA_ENV(h,100,m,ma_results,eid=4)#50/50
#        modeCloseENV(h,100,m,ma_results,eid=5)#50/50
#        modeCloseMA_ENV(h,100,m,ma_results,eid=6,w1=0.3333,w2=0.33333,w3=0.44)
#        modeCloseMA_ENV(h,100,m,ma_results,eid=7,w1=0.6,w2=0.2,w3=0.2)
#        a) MA free 0%, MA-MA 100%, ENV free 0%, ENV-ENV 0%, ENV-MA 100%
#        modeCloseMA_ENV(h,1,m,[],eid=7,w1=0.6,w2=0.2,w3=0.2)
#    if not NOGUI :
#        afviewer.displayPreFill()
#        afviewer.displayFill()
#    displayResult(h,0,2,0)
elif disp :
    #execfile("/Users/ludo/DEV/testGaussCluster.py")
    h.FillName=["test",]*(h.cFill+1)
    n=0;displayResult(h,n,n,50);
    #a,fluos,locis,npx,limits=displayOne(helper,n,n,50,20.0,None,space=100,redo=True)
#    execfile("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/hiv_exp.py")
else :
    resultfilename = h.resultfile = "/Users/ludo/Downloads/autopack_ouput/hiv_experiment/hiv_exp"    
    h.buildGrid(boundingBox=h.boundingBox,gridFileIn=None,rebuild=False ,
                        gridFileOut=None,previousFill=True)
    h.fill5(seedNum=10,verbose = 0,usePP=False) 
    if not NOGUI :
        afviewer.displayPreFill()
        afviewer.displayFill()
    