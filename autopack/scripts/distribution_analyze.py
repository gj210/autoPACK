# -*- coding: utf-8 -*-
"""
Created on Fri May  9 14:19:45 2014

@author: ludo
"""
import sys
import json
import numpy as np
import math
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14_95696CEA/plugins/ePMV/mgl64/MGLToolsPckgs")
import os
import scipy
from scipy import spatial
import matplotlib
from matplotlib import pyplot as plt
from autopack.Analysis import AnalyseAP
from autopack.GeometryTools import GeometriTools,Rectangle
analyse = AnalyseAP(env=None, viewer=None, result_file=None)
analyse.bbox=[[0, 0, 0], [1000, 1000, 1]]

from matplotlib import rcParams

#set plot attributes
fig_width = 12  # width in inches
fig_height = 10  # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'osx',
          'axes.labelsize': 8,
          'axes.titlesize': 8,
          'font.size': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'figure.figsize': fig_size,
          'savefig.dpi' : 300,
          'font.family': 'sans-serif',
          'axes.linewidth' : 0.5,
          'xtick.major.size' : 2,
          'ytick.major.size' : 2,
          'font.size' : 8
          }
rcParams.update(params)
import numpy
def vector_norm(data, axis=None, out=None):
    """Return length, i.e. Euclidean norm, of ndarray along axis.
    >>> v = numpy.random.random(3)
    >>> n = vector_norm(v)
    >>> numpy.allclose(n, numpy.linalg.norm(v))
    True
    >>> v = numpy.random.rand(6, 5, 3)
    >>> n = vector_norm(v, axis=-1)
    >>> numpy.allclose(n, numpy.sqrt(numpy.sum(v*v, axis=2)))
    True
    >>> n = vector_norm(v, axis=1)
    >>> numpy.allclose(n, numpy.sqrt(numpy.sum(v*v, axis=1)))
    True
    >>> v = numpy.random.rand(5, 4, 3)
    >>> n = numpy.empty((5, 3))
    >>> vector_norm(v, axis=1, out=n)
    >>> numpy.allclose(n, numpy.sqrt(numpy.sum(v*v, axis=1)))
    True
    >>> vector_norm([])
    0.0
    >>> vector_norm([1])
    1.0
    """
    data = numpy.array(data, dtype=numpy.float64, copy=True)
    if out is None:
        if data.ndim == 1:
            return math.sqrt(numpy.dot(data, data))
        data *= data
        out = numpy.atleast_1d(numpy.sum(data, axis=axis))
        numpy.sqrt(out, out)
        return out
    else:
        data *= data
        numpy.sum(data, axis=axis, out=out)
        numpy.sqrt(out, out) 
def angle_between_vectors(v0, v1, directed=True, axis=0):
    """Return angle between vectors.
    If directed is False, the input vectors are interpreted as undirected axes,
    i.e. the maximum angle is pi/2.
    >>> a = angle_between_vectors([1, -2, 3], [-1, 2, -3])
    >>> numpy.allclose(a, math.pi)
    True
    >>> a = angle_between_vectors([1, -2, 3], [-1, 2, -3], directed=False)
    >>> numpy.allclose(a, 0)
    True
    >>> v0 = [[2, 0, 0, 2], [0, 2, 0, 2], [0, 0, 2, 2]]
    >>> v1 = [[3], [0], [0]]
    >>> a = angle_between_vectors(v0, v1)
    >>> numpy.allclose(a, [0, 1.5708, 1.5708, 0.95532])
    True
    >>> v0 = [[2, 0, 0], [2, 0, 0], [0, 2, 0], [2, 0, 0]]
    >>> v1 = [[0, 3, 0], [0, 0, 3], [0, 0, 3], [3, 3, 3]]
    >>> a = angle_between_vectors(v0, v1, axis=1)
    >>> numpy.allclose(a, [1.5708, 1.5708, 1.5708, 0.95532])
    True
    """
    v0 = numpy.array(v0, dtype=numpy.float64, copy=False)
    v1 = numpy.array(v1, dtype=numpy.float64, copy=False)
    dot = numpy.sum(v0 * v1, axis=axis)
    dot /= vector_norm(v0, axis=axis) * vector_norm(v1, axis=axis)
    return numpy.arccos(dot if directed else numpy.fabs(dot))    

def signed_angle_between_vectors(Vn,v0, v1, directed=True, axis=0):
    Vn = numpy.array(Vn)
    angles = angle_between_vectors(v0, v1, directed=directed, axis=axis)
    cross = numpy.cross(v0,v1)
    dot = numpy.dot(cross,Vn)
#    ind=numpy.nonzero(dot < 0)
    if dot < 0 :
        angles*=-1.0
    return angles
    
def ApplyMatrix(coords,mat):
        """
        Apply the 4x4 transformation matrix to the given list of 3d points.
    
        @type  coords: array
        @param coords: the list of point to transform.
        @type  mat: 4x4array
        @param mat: the matrix to apply to the 3d points
    
        @rtype:   array
        @return:  the transformed list of 3d points
        """
    
        #4x4matrix"
        mat = numpy.array(mat)
        coords = numpy.array(coords)
        one = numpy.ones( (coords.shape[0], 1), coords.dtype.char )
        c = numpy.concatenate( (coords, one), 1 )
        return numpy.dot(c, numpy.transpose(mat))[:, :3]
#output="/Users/ludo/DEV/autoPACKresults/NM_Analysis_B5"
#position_file=output+os.sep+"pos"
#total_positions = numpy.genfromtxt(position_file, delimiter=',')
#[u'ext__Sphere_radius_50', u'ext__Sphere_radius_25p', u'ext__Sphere_radius_200', u'ext__Sphere_radius_25r', u'ext__Sphere_radius_100']
def gatherResultPerIngredient(path,N):
    all_pos_ingr={}
#    for n in [u'ext__Sphere_radius_50', u'ext__Sphere_radius_25p', u'ext__Sphere_radius_200', u'ext__Sphere_radius_25r', u'ext__Sphere_radius_100']:
    for n in [u'ext__Bacteria_Rad25_1_3', u'ext__IngredientC_1_1', u'ext__snake']:
       all_pos_ingr[n] =[] 
    for i in range(N):
        filename=path+os.sep+"_posIngr_"+str(i)+".json"
        with open(filename, 'r') as fp :#doesnt work with symbol link ?
            d=json.load(fp)
        for k in d :
            all_pos_ingr[k].append(np.array(d[k]))
    return all_pos_ingr
    
def gatherAllResultPerIngredient(path,N):
    all_pos_rot_ingr={}
#    for n in [u'ext__Sphere_radius_50', u'ext__Sphere_radius_25p', u'ext__Sphere_radius_200', u'ext__Sphere_radius_25r', u'ext__Sphere_radius_100']:
    for n in [u'ext__Bacteria_Rad25_1_3', u'ext__IngredientC_1_1', u'ext__snake']:
       all_pos_rot_ingr[n] =[] 
    for i in range(N):
        filename=path+os.sep+"results_seed_"+str(i)+".json"
        with open(filename, 'r') as fp :#doesnt work with symbol link ?
            d=json.load(fp)
        R=d["exteriorRecipe"]
        for k in R :
            all_pos_rot_ingr[k].append(np.array(R[k]["results"]))
    return all_pos_rot_ingr

def gatherResult(path,N):
    ALL_POS=[]
    for i in range(N):
        p=[]
        filename=path+os.sep+"_posIngr_"+str(i)+".json"
        with open(filename, 'r') as fp :#doesnt work with symbol link ?
            d=json.load(fp)
        for k in d :
            p.extend(np.array(d[k]))
        ALL_POS.append(p)
    return np.array(ALL_POS)
    
def dumpResultsCSV(M,path,N):
    all_pos=[]
    output=path+str(M)
    for i in range(N):
        filename=output+os.sep+"pos_"+str(i)
        with open(filename, 'r') as fp :#doesnt work with symbol link ?
            pos=json.load(fp)
        all_pos.extend(np.array(pos[0]))
    #    g+=analyse.rdf(pos[0],dr=dr)
    todump=[]
    for p in all_pos :
        todump.append(p.tolist())
    numpy.savetxt(output+os.sep+"pos_total.csv", np.array(all_pos), delimiter=",")

def dumpResults(M,path,N):
    all_pos=[]
    output=path+str(M)
    for i in range(N):
        filename=output+os.sep+"pos_"+str(i)
        with open(filename, 'r') as fp :#doesnt work with symbol link ?
            pos=json.load(fp)
        all_pos.append(np.array(pos[0]))
    #    g+=analyse.rdf(pos[0],dr=dr)
    todump=[]
    for p in all_pos :
        todump.append(p.tolist())
    with open(output+os.sep+"pos_total", 'w') as fp :
        json.dump(todump,fp)
    return todump

def getOccurencePerIngredient(M,path,N):
    ingrOcc={}
    ingr_name=[ u'ext__Sphere_radius_25r',
         u'ext__Sphere_radius_25p',
         u'ext__Sphere_radius_50',
         u'ext__Sphere_radius_200',
         u'ext__Sphere_radius_100']
    for iname in  ingr_name :
        ingrOcc[iname]=[]
    output=path+str(M)
    for i in range(N):
        filename=output+os.sep+"pos_"+str(i)+"_results.json.json"#pos_999_results.json.json
        with open(filename, 'r') as fp :#doesnt work with symbol link ?
            r=json.load(fp)#po,rot,ingr,
        for iname in  ingr_name :     
            ingrOcc[iname].append(len(r['exteriorRecipe'][iname]['results']))
    return ingrOcc,ingr_name
    
def plotOccurencePerIngredient(occurence,name,filename):
    fig, ax = plt.subplots()
    fig.set_size_inches(12, 10)
    #plt.title("Distribution along axis X and Y for 1000 run with %d ingredients" % M)
    B=[]
    V=[]
    S=[]
    for i,n in enumerate(name) :
        v=np.array(occurence[n]).mean()
        s=np.array(occurence[n]).std()
        B.append(i*10)
        V.append(v)
        S.append(s)
    w = 5
    ba=ax.bar(B, V, w,color='r', yerr=S)
    ax.set_xticks(B)
    ax.set_xticklabels(name) 
    autolabel(ax,ba,err=None)
    autolabel(ax,ba,err=S)
    plt.savefig(filename,dpi=300)
    
def getHistData(data,axis,width):
    b=range(0,width+100,100)
    N=np.zeros(len(b)-1)
    Nn=[]   
    for run in data :
        if not len(run) :
            n=np.zeros(len(b)-1)
        else :
            n,bi = np.histogram(np.array(run).transpose()[axis],bins=b)
        N+=n
        Nn.append(n)
        print len(Nn)
    v=np.average(Nn,axis=0)
    s=np.std(Nn,axis=0)#sqrt of N axis 0
#    s=stats.stem
    s=1.96*(s/np.sqrt(1000))
    return v,s,Nn,np.average(np.sum(Nn,axis=1))/10.0

def autolabel(ax,rects,err=None):
    # attach some text labels
    for i,rect in enumerate(rects):
        height = rect.get_height()
        v='%.1f'%height
        y=0.5*height
        if err is not None :
            v='%.2f'%err[i]
            y=1.05*height
        ax.text(rect.get_x()+rect.get_width()/2., y, v,
                ha='center', va='bottom')
                
def plot_all_histogram(M,data,path,width=1000,name="histXY.png",labels=False):
    print "Distribution along axis X for 1000 run with %d ingredients"
    output=path+str(M)
    b=np.arange(0,width+100,100)
    #fig, ax = plt.subplots()
    vx,sx,Nnx,Ax = getHistData(data,0,width)
    vy,sy,Nny,Ay = getHistData(data,1,width)
    w=width/len(b[:-1])
    plt.close('all')
    fig, ax = plt.subplots()
#    fig.set_size_inches(12, 10)
    plt.title("Distribution along axis X and Y for 1000")
    #ax.hist(mydata, weights=np.zeros_like(data) + 1. / data.size)
    bx=ax.bar(b[:-1], vx, w/2.0,color='grey',ecolor='grey', yerr=sx)
    by=ax.bar(np.array(b[:-1])+w/2.0, vy, w/2.0,color='darkgrey',ecolor='darkgrey', yerr=sy)
    plt.hlines(Ax,1,999)
    plt.hlines(Ay,1,999)
    #ax.legend( (bx[0], by[0]), ('X', 'Y') )
    if labels :
        autolabel(ax,bx,err=None)
        autolabel(ax,bx,err=sx)
        autolabel(ax,by,err=None)
        autolabel(ax,by,err=sy)
    plt.savefig(output+os.sep+name,dpi=300)
    return vx,sx,vy,sy,Nnx,Nny

def getHistDataAngle(data,b=np.arange(-180,210,30)):
    N=np.zeros(len(b)-1)
    Nn=[]   
    for run in data :
        if not len(run) :
            n=np.zeros(len(b)-1)
        else :
            n,bi = np.histogram(np.array(run),bins=b)
        N+=n
        Nn.append(n)
        print len(Nn)
    v=np.average(Nn,axis=0)
    s=np.std(Nn,axis=0)#sqrt of N axis 0
#    s=stats.stem
    s=1.96*(s/np.sqrt(1000))
    return bi,v,s,Nn,np.average(np.sum(Nn,axis=1))/10.0

def plotAngle(data,output,name,labels=False):
#    s=stats.stem
    b=np.arange(-180,210,30)
    #fig, ax = plt.subplots()
    bi,vx,sx,Nnx,Ax = getHistDataAngle(data[0])
    bi,vy,sy,Nny,Ay = getHistDataAngle(data[1])
    w=360.0/len(b)
    plt.close('all')
    fig, ax = plt.subplots()
#    fig.set_size_inches(12, 10)
    plt.title("Z Angle Distribution for 1000 run")
    #ax.hist(mydata, weights=np.zeros_like(data) + 1. / data.size)
#    bx=ax.bar(b[:-1], vx, w/2.0,color='r', yerr=sx)
#    by=ax.bar(np.array(b[:-1])+w/2.0, vy, w/2.0,color='b', yerr=sy)
    bincenters = 0.5*(b[1:]+b[:-1])
    bx=ax.bar(bincenters-w/2.0, vx, w/2.0,color='grey',ecolor='grey', yerr=sx)
    by=ax.bar(np.array(bincenters), vy, w/2.0,color='darkgrey',ecolor='darkgrey', yerr=sy)
#    plt.hlines(Ax,1,999)
#    plt.hlines(Ay,1,999)
#    ax.legend( (bx[0], by[0]), ('X', 'Y') )
    ax.set_xticks(b )
    ax.set_xticklabels(range(-180,210,30)) 
    if labels:
        autolabel(ax,bx,err=None)
        autolabel(ax,bx,err=sx)
        autolabel(ax,by,err=None)
        autolabel(ax,by,err=sy)
    plt.savefig(output+os.sep+name,dpi=300)  
    return bi,vx,sx,Nnx,Ax,vy,sy,Nny,Ay

def plot_xyz_histogram(M,data,path,width=1000,filename=None):
    print "Distribution along axis X for 1000 run with %d ingredients"
    output=path+str(M)
    b=range(0,width+100,100)
    #fig, ax = plt.subplots()
    plt.close('all')
    fig, ax = plt.subplots()
    fig.set_size_inches(12, 10)
    plt.title("Distribution along axis X Y Z for 1000 run with %d ingredients" % M)
    vx,sx = getHistData(data,0,width)
    vy,sy = getHistData(data,1,width)
    vz,sz = getHistData(data,2,width)
    w=width/len(b[:-1])
    bx=ax.bar(b[:-1], vx, w/3.0,color='r', yerr=sx)
    by=ax.bar(np.array(b[:-1])+w/3.0, vy, w/3.0,color='b', yerr=sy)
    bz=ax.bar(np.array(b[:-1])+2*(w/3.0), vz, w/3.0,color='g', yerr=sz)
    ax.legend( (bx[0], by[0], bz[0]), ('X', 'Y','Z') )
    autolabel(ax,bx,err=None)
    autolabel(ax,bx,err=sx)
    autolabel(ax,by,err=None)
    autolabel(ax,by,err=sy)
    autolabel(ax,bz,err=None)
    autolabel(ax,bz,err=sz)
    if filename is None :
        filename=output+os.sep+"hist.png"
    plt.savefig(filename,dpi=300)
    return vx,sx,vy,sy

def plot_histogram(M,data,path,axis=0,width=1000):
    print "Distribution along axis X for 1000 run with %d ingredients"
    b=range(0,width+100,100)
    N=np.zeros(len(b)-1)
    Nn=[]   
    output=path+str(M)
    for run in data :
        n,bi = np.histogram(np.array(run).transpose()[axis],bins=b)
        N+=n
        Nn.append(n)
        print len(Nn)
    val=N/len(data)
    v=np.average(Nn,axis=0)
    s=np.std(Nn,axis=0)#sqrt of N axis 0
    menStd = np.sqrt(val)
    #fig, ax = plt.subplots()
    width=width/len(b[:-1])
    plt.close('all')
    fig, ax = plt.subplots()
    fig.set_size_inches(12, 10)
    plt.title("Distribution along axis X for 1000 run with %d ingredients" % M)
    plt.bar(b[:-1], v, width,color='r', yerr=s)
    plt.savefig(output+os.sep+"hist"+str(axis)+".png",dpi=300)

def doit(M,n,ALL_POS):
    path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_P_"
    all_pos=dumpResults(M,path,n)
    ALL_POS.append(all_pos)
    path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_nP_"
    all_pos=dumpResults(M,path,n)
    ALL_POS.append(all_pos)
    with open("/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_pos_total.json", 'w') as fp :
        json.dump(ALL_POS,fp)
    path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_P_"
#    plot_histogram(M,ALL_POS[0],path)
    plot_all_histogram(M,ALL_POS[0],path)
    path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_nP_"
#    plot_histogram(M,ALL_POS[1],path)
    plot_all_histogram(M,ALL_POS[1],path)
    return ALL_POS

ALL_POS=[]
ALL_POSB=doit(5,1000,ALL_POS)
path="/Users/ludo/DEV/autoPACKresults"
vx,sx,vy,sy,Nnx,Nny = plot_all_histogram("",ALL_POSB[0],path,width=1000,name="FigureB_histXY.svg",labels=True)
#path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_nP_"
#vx,sx,vy,sy,Nnx,Nny = plot_all_histogram("",ingr_pos[n],path,width=1000,name="FigureB_histXY.png")
#path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_P_"
#ALL_POSA=gatherResult(path,1000)
#dumpResultsCSV(5,path,1000)
#ingrOcc,ingr_name=getOccurencePerIngredient(5,path,1000)
#plotOccurencePerIngredient(ingrOcc,ingr_name,"/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_P_5/occurences.png")
#path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_nP_"
#dumpResultsCSV(5,path,1000)
#ingrOcc,ingr_name=getOccurencePerIngredient(5,path,1000)
#plotOccurencePerIngredient(ingrOcc,ingr_name,"/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_nP_5/occurences.png")
#
path ="/Users/ludo/DEV/autoPACKresults/NM_Analysis_A2"
ALL_POSA=gatherResult(path,1000)
#path="/Users/ludo/DEV/autoPACKresults"
###plot_all_histogram("",ALL_POSA,path)
vx,sx,vy,sy,Nnx,Nny = plot_all_histogram("",ALL_POSA,path,width=1000,name="FigureA_histXY.svg",labels=True)
#
#
path ="/Users/ludo/DEV/autoPACKresults/NM_Analysis_C3"
ingr_pos=gatherResultPerIngredient(path,1000)
ipr = gatherAllResultPerIngredient(path,1000)
#corner is [0,0]
#angle corner-position vs direction = 1,0,0 * rot
#ApplyMatrix([[1,0,0]],R[0][0])
#rot = [ run.transpose()[1] for run in ipr.values()[0] ]
#rot = [ r.tolist() for r in rot ]
#pos = [ run.transpose()[1] for run in ipr.values()[0] ]
#pos = [ r.tolist() for r in rot ]
corner=np.array([0,0,0])
ANGLES=[]
for j in range(2):
    angles=[]
    for run in ipr.values()[j] :
        p,r=run.transpose()
        an=[]
        for i in range(len(p)):
            v1=corner-p[i]
            v2=ApplyMatrix([[0,1,0]],r[i])[0]
#            a=signed_angle_between_vectors([0,0,1],v1,v2)
            a=angle_between_vectors(v1,v2)
            an.append(a)
        angles.append(an)
    ANGLES.append(np.degrees(np.array(angles,'f')))
SANGLES=np.array(ANGLES)
UANGLES=np.array(ANGLES)
from autopack.transformation import euler_from_matrix,matrixToEuler,euler_matrix,superimposition_matrix,rotation_matrix,affine_matrix_from_points
R=np.array(rot,'f')
rad=[]
for run in R :
    a=[euler_from_matrix(ro)[2] for ro in run]
    rad.append(a)
A=np.array(angles,'f')
AR=np.array(rad,'f')
##plot_all_histogram("",ALL_POS,path)
for n in [u'ext__Bacteria_Rad25_1_3', u'ext__IngredientC_1_1', u'ext__snake']:
#    vx,sx,vy,sy,Nnx,Nny=plot_all_histogram("",ingr_pos[n],path,width=1000,name=n+"_histXY_nolabels.ps")
    vx,sx,vy,sy,Nnx,Nny=plot_all_histogram("",ingr_pos[n],path,width=1000,name=n+"_histXY_labels.svg",labels=True)
#    
##Z=1.96
#Z*std/np.sqrt(N)
#z,pval = mstats.normaltest(mx)

#if(pval < 0.055):
#    print "Not normal distribution"
ANGLES=[[],[],[]]
for j,n in enumerate([u'ext__Bacteria_Rad25_1_3', u'ext__IngredientC_1_1', u'ext__snake']):
    f="/Users/ludo/DEV/autoPACKresults/eulerFigure3C/"+n+"_euler_Z.csv"
    angles = numpy.loadtxt(f)
    #divide by number of individus
    v=0
    for i in range(1000):
        nv=len(ingr_pos[n][i])
        ANGLES[j].append(angles[v:v+nv])
        v+=nv
#    ANGLES.append(angles)
output ="/Users/ludo/DEV/autoPACKresults/"
plotAngle(ANGLES,output,"angle_hist.svg")

plotAngle(ANGLES,output,"angle_hist.svg",labels=True)

#make the svg figure
import svgutils.transform as sg
import sys 

def makeSVG(output,sc=0.435,W2=680/2 ,H1=300):
    #create new SVG figure
    fig = sg.SVGFigure("19.05cm", "22.86cm")
    
    # load matpotlib-generated figures
    figA2 = sg.fromfile(output+os.sep+'FigureB_histXY.svg')
    figA1= sg.fromfile(output+os.sep+'FigureA_histXY.svg')
    figB4 = sg.fromfile(output+os.sep+'angle_hist.svg')
    figB=[]
    for n in [u'ext__Bacteria_Rad25_1_3', u'ext__IngredientC_1_1', u'ext__snake']:
        fig2 = sg.fromfile(output+os.sep+n+"_histXY_labels.svg")
        figB.append(fig2) 
        
    figB.append(figB4)    
    
#    W2 = 680/2    
#    H1 = 350
    # get the plot objects
    plot1 = figA1.getroot()
    plot1.moveto(5, 0, scale=sc)
    plot2 = figA2.getroot()
    plot2.moveto(W2, 0, scale=sc)
    
    plot3 = figB[0].getroot()
    plot3.moveto(5, H1, scale=sc)
    plot4 = figB[1].getroot()
    plot4.moveto(W2, H1, scale=sc)
    
    plot5 = figB[2].getroot()
    plot5.moveto(5, H1*2, scale=sc)
    plot6 = figB[3].getroot()
    plot6.moveto(W2, H1*2, scale=sc)
    
    # add text labels
    txt1 = sg.TextElement(25,20, "A.1", size=12, weight="bold")
    txt2 = sg.TextElement(W2+5,20, "A.2", size=12, weight="bold")
    
    txt3 = sg.TextElement(25,20+H1, "B.1", size=12, weight="bold")
    txt4 = sg.TextElement(W2+5,20+H1, "B.2", size=12, weight="bold")
    
    txt5 = sg.TextElement(25,20+H1*2, "B.3", size=12, weight="bold")
    txt6 = sg.TextElement(W2+5,20+H1*2, "B.4", size=12, weight="bold")
    
    # append plots and labels to figure
    fig.append([plot1, plot2])
    fig.append([plot3, plot4])
    fig.append([plot5, plot6])
    fig.append([txt1, txt2])
    fig.append([txt3, txt4])
    fig.append([txt5, txt6])
    
    # save generated SVG files
    fig.save(output+os.sep+"fig_final.svg")

TWOD=False
THREED=False
if TWOD:
    ALL_POSP=[]
    ALL_POSnP=[]
    for M in range(1,6):
        if M == 4 :
            n=579
        else :
            n=1000
        print M,n
        path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_P_"
        dumpResultsCSV(M,path,n)
        all_pos=dumpResults(M,path,n)
        ALL_POSP.append(all_pos)
        path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_nP_"
        all_pos=dumpResults(M,path,n)
        dumpResultsCSV(M,path,n)
        ALL_POSnP.append(all_pos)
#    with open("/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_pos_total.json", 'w') as fp :
#        json.dump(ALL_POS,fp)
    
    #with open("/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_pos_total.json", 'r') as fp :
    #    ALL_POS=json.load(fp)
    #path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_P_"
    #plot_histogram(0,ALL_POS[0],path)
    ##histogram for all on X and Y + error Bar
    for M in range(1,6):
        path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_P_"
#        plot_histogram(M,ALL_POSP[M-1],path)
        plot_all_histogram(M,ALL_POSP[M-1],path)
        path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_nP_"
        plot_all_histogram(M,ALL_POSnP[M-1],path)
#        plot_histogram(M,ALL_POSnP[M-1],path)
if THREED :
    #3D
    ALL_POS=[]
    n=52
    M=5
    path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_3D_P_"
    all_pos=dumpResults(M,path,n)
    ALL_POS.append(all_pos)
    path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_3D_nP_"
    all_pos=dumpResults(M,path,n)
    ALL_POS.append(all_pos)
    with open("/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_pos_total.json", 'w') as fp :
        json.dump(ALL_POS,fp)
    path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_3D_P_"
    plot_xyz_histogram(M,ALL_POS[0],path,width=1000,filename="/Users/ludo/DEV/autoPACKresults/histXYZ3DP.png")
#    plot_histogram(M,ALL_POS[0],path,axis=0,width=1000)
#    plot_histogram(M,ALL_POS[0],path,axis=1,width=1000)
#    plot_histogram(M,ALL_POS[0],path,axis=2,width=1000)
    path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_2B_3D_nP_"
    plot_xyz_histogram(M,ALL_POS[1],path,width=1000,filename="/Users/ludo/DEV/autoPACKresults/histXYZ3DnP")
#    plot_histogram(M,ALL_POS[1],path,axis=0,width=1000)
#    plot_histogram(M,ALL_POS[1],path,axis=1,width=1000)
#    plot_histogram(M,ALL_POS[1],path,axis=2,width=1000)


#output="/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_P_"+str(1)
#filename=output+os.sep+"pos_"+str(0)
#with open(filename, 'r') as fp :#doesnt work with symbol link ?
#    pos=json.load(fp)        
#positions=pos[0]
#rdf=pos[1]
#dr=40.0
#rMax =  np.sqrt(1000**2+1000**2)
#edges = np.arange(0., rMax+1.1*dr, dr)
##G=analyse.rdf(positions,dr=dr)
##plt.plot(edges[:-1],G);plt.show()
#K,L=analyse.ripley(positions,dr=dr)
##plt.plot(edges[:-1],L);plt.show()
#bbox = analyse.bbox
#rect=Rectangle(bbox[1][1],bbox[0][1],bbox[1][0],bbox[0][0])
#p=[503.43258941010185, 872.58598615819722, 0.0]
#m=[p[0],p[1]]
#r=140.0
#area = math.pi*r**2
#chs = analyse.g.check_sphere_inside(rect,m,r)
#ch=analyse.g.check_rectangle_oustide(rect,m,r)
#leftBound,rightBound = analyse.g.getBoundary(rect,m,r)
#area = analyse.g.get_rectangle_cercle_area(rect,m,r,rightBound,leftBound)
#self = analyse.g
##draw the rectangle and the circle
#fig = plt.figure()
#ax = fig.add_subplot(111)
#from matplotlib.patches import Circle
#from matplotlib.patches import Rectangle
#ax.add_patch(Rectangle([0.0,0.0],1000,1000))
#ax.add_patch(Circle((p[0],p[1]), r, facecolor='red'))
#ax.vlines(leftBound,0,1000)
#ax.vlines(rightBound,0,1000)
#ax.set_aspect(1.0)
#plt.axhline(y=bbox[0][1], color='k')
#plt.axhline(y=bbox[1][1], color='k')
#plt.axvline(x=bbox[0][0], color='k')
#plt.axvline(x=bbox[1][0], color='k')
#plt.axis([bbox[0][0], bbox[1][0],
#             bbox[0][1], bbox[1][1]])
#
#
#a=0.0
#i=leftBound + self.Resolution
#while (i <= rightBound):
#    upperBound = min(self.UpperRectangleFunction(rect, i - self.Resolution / 2.0), 
#                     self.UpperCircleFunction(m, r, i - self.Resolution / 2.0));  
#    lowerBound = max(self.LowerRectangleFunction(rect, i - self.Resolution / 2.0), 
#                     self.LowerCircleFunction(m, r, i - self.Resolution / 2.0));  
#    a += (upperBound - lowerBound ) * self.Resolution;
#    print upperBound,lowerBound,(upperBound - lowerBound ) * self.Resolution,a
##    ax.add_patch(Rectangle([0.0,0.0],1000,1000))
#    i+=self.Resolution
#    
#
#
#
#ingrpositions={}
#rangeseed=range(1000)
#def merge(d1, d2, merge=lambda x,y:y):
#        result = dict(d1)
#        for k,v in d2.iteritems():
#            #print k
#            if k in result:
#                result[k].append(v)
#            else:
#                result[k] = [v,]
#        return result
#def rdf(positions,dr=10,rMax=None):
#    N=len(positions)
#    V=1000**2
#    diag = maxdist = np.sqrt(1000**2+1000**2)
##    r=np.array([25,25,50,75])
#    #all_distance = scipy.spatial.distance.pdist(pos)
#    dr=dr#all_distance.min()
#    if rMax is None :
#        rMax=diag
#    edges = np.arange(0., rMax+1.1*dr, dr)
##    n,radii = np.histogram(all_distance,bins=edges)
#    g=np.zeros((N,len(edges)-1))
#    dv=[]
#    density=float(N)/float(V)
#    for i,p in enumerate(positions) :
#        di = scipy.spatial.distance.cdist(positions, [p,], 'euclidean')
#        dN,bins = np.histogram(di,bins=edges)
#        dV = np.array(analyse.getAreaShell(analyse.bbox,edges,p))
#        dv.append(dV)
#        g[i] = dN/(dV*density)
#    return np.average(g,axis=0)#/np.array(dv)
#
#def ripleyK(positions):
#    #K(t) = lambda^-1*SUM(I(dij<t)/n)
#    #K(t) = A*SUM(wij*I(i,j)/n**2)
#    #lambda = n/A A is the area of the region containing all points
#    #I indicator function 1 if its operand is true, 0 otherwise
#    #t is the search radius
#    #if homogenous K(s) = pi*s**2
#    #L(t) = (K(t)/pi)**1/2
#    #A common plot is a graph of t - \hat{L}(t) against t
#    #which will approximately follow the horizontal zero-axis with constant 
#    #dispersion if the data follow a homogeneous Poisson process.
#    N=len(positions)
#    V=1000**2
#    diag = maxdist = np.sqrt(1000**2+1000**2)
##    r=np.array([25,25,50,75])
#    #all_distance = scipy.spatial.distance.pdist(pos)
#    dr=dr#all_distance.min()
#    if rMax is None :
#        rMax=diag
#    edges = np.arange(0., rMax+1.1*dr, dr)
##    n,radii = np.histogram(all_distance,bins=edges)
#    k=np.zeros((N,len(edges)-1))
#    dv=[]
#    density=float(N)/float(V)
#    for i,p in enumerate(positions) :
#        di = scipy.spatial.distance.cdist(positions, [p,], 'euclidean')
#        #dV = np.array(analyse.getAreaShell(analyse.bbox,edges,p))
#        for j,e in enumerate(edges) :
#            area0=math.pi*r**2#complete circle
#            area1=analyse.rectangle_circle_area(analyse.bbox,p,edges)
#            w=area1/area0
#            k[i,j]=w*np.sum(np.nonzero(di < t))/N**2
#    Kt=V*np.sum(k,axis=0)
#    Lt=(Kt/pi)**0.5
#    #K(t) = lambda^-1*SUM(I(dij<t)/n)      
##    lambda1 = (n/A)**-1.0
##    i = d < t
##    I=np.array(i,dtype=int)/n
##    Kt=lambda1*I
##    Kl=(Kt/pi)**0.5
#    return Kt,Lt
#
#
#    
#pos=[]
#for i in rangeseed:
#    dict1= analyse.loadJSON(output+os.sep+"_posIngr_"+str(i)+".json")
#    p = [np.array(p) for p in dict1.values()]
#    pos.append(np.vstack(p))
#    #ingrpositions=dict(merge(ingrpositions,dict1))
#ingrnames=ingrpositions.keys() 
#dr=5.0
#rMax =  np.sqrt(1000**2+1000**2)
#edges = np.arange(0., rMax+1.1*dr, dr)
#g=np.zeros(len(edges)-1)
#for i,p in enumerate(pos[:100]) :
#    print i
#    g+=rdf(p,dr=dr)
#G=g/len(pos)
#plt.plot(edges[1:-1],G[1:]);plt.show()
#
#dr=40.0
#rMax =  np.sqrt(1000**2+1000**2)
#edges = np.arange(0., rMax+1.1*dr, dr)
#g=np.zeros(len(edges)-1)
#all_pos=[]
#M=4
#output="/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_P_"+str(M)
#for i in range(580):
#    filename=output+os.sep+"pos_"+str(i)
#    with open(filename, 'r') as fp :#doesnt work with symbol link ?
#        pos=json.load(fp)
#    all_pos.append(np.array(pos[0]))
##    g+=analyse.rdf(pos[0],dr=dr)
#todump=[]
#for p in all_pos :
#    todump.append(p.tolist())
#with open(output+os.sep+"pos_total", 'w') as fp :
#    json.dump(todump,fp)
#with open(output+os.sep+"pos_total", 'r') as fp :
#    all_pos1=pos=json.load(fp)
#
#
#dr=40
#rMax=800.0
#analyse.g.Resolution = 6.0
#count=0
#edges = np.arange(0., rMax+1.1*dr, dr)
#g=np.zeros(len(edges)-1)
#for pos in all_pos:
#    print ">",count    
#    g+=analyse.rdf(pos,dr=dr,rMax=rMax)
#    count+=1
#G=g/len(pos)
#plt.plot(edges[1:-1],G[1:]);plt.show()
#
#
#
#pos=np.array([np.array(p[0]) for p in h.molecules])
#N=len(pos)
#V=1000**2
#diag = maxdist = np.sqrt(1000**2+1000**2)
#r=np.array([25,25,50,75])
#all_distance = scipy.spatial.distance.pdist(pos)
#dr=20.0#all_distance.min()
#rMax=600.0
#edges = np.arange(0., rMax+1.1*dr, dr)
#n,radii = np.histogram(all_distance,bins=edges)
#g=np.zeros((N,len(n)))
#dv=[]
#density=float(N)/float(V)
#for i,p in enumerate(pos) :
#    di = scipy.spatial.distance.cdist(pos, [p,], 'euclidean')
#    dN,bins = np.histogram(di,bins=radii)
#    dV = np.array(analyse.getAreaShell(analyse.bbox,radii,p))
#    dv.append(dV)
#    g[i] = (dN/density)/dV
##    g[i]= ((dN/float(N))/dV)/density
#    #whatthe number density ?
##g_average=np.zeros(len(n))
##for i in range(len(n)): 
##    g_average[i] = np.average(g[:,i])
#G=np.average(g,axis=0)#/np.array(dv)
##G=np.sum(g,axis=0)#/ len(g)#(dV/V)
#plt.plot(radii[1:-1],G[1:]);plt.show()
##doesnt tend to 0 ?
#bins = []
#for i in range(len(r)-1) :
#    for j in range(i,len(r)) :
#        print i,j,r[i]+r[j]
#        bins.append(r[i]+r[j])
#
#ALL_POS=[]
#for M in range(5):
#    path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_P_"
#    all_pos=dumpResults(M,path,1000)
#    ALL_POS.append(all_pos)
#    path="/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_nP_"
#    all_pos=dumpResults(M,path,1000)
#    ALL_POS.append(all_pos)
#with open("/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_pos_total.json", 'w') as fp :
#    json.dump(ALL_POS,fp)
#
#
#all_pos=[]
#M=3
#output="/Users/ludo/DEV/autoPACKresults/NM_Analysis_B_nP_"+str(M)
#for i in range(1000):
#    filename=output+os.sep+"pos_"+str(i)
#    with open(filename, 'r') as fp :#doesnt work with symbol link ?
#        pos=json.load(fp)
#    all_pos.append(np.array(pos[0]))
##    g+=analyse.rdf(pos[0],dr=dr)
#todump=[]
#for p in all_pos :
#    todump.append(p.tolist())
#with open(output+os.sep+"pos_total", 'w') as fp :
#    json.dump(todump,fp)
#
#
#N=np.zeros(10)   
#b=range(0,1100,100)
#Nn=[]   
#for run in all_pos :
#    n,bi = np.histogram(run.transpose()[0],bins=b)
#    N+=n
#    Nn.append(n)
#val=N/len(all_pos)
#v=np.average(Nn,axis=0)
#s=np.std(Nn,axis=0)
#menStd = np.sqrt(val)
##fig, ax = plt.subplots()
#width=1000.0/len(b[:-1])
#loci1 = plt.bar(b[:-1], v, width,color='r', yerr=s)
##error bar ?
#yerr=
# [50, 50, 75, 100, 50, 75, 100, 100, 125]
    