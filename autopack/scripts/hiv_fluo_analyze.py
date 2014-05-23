# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 09:54:06 2014

@author: ludo
"""
import sys
sys.path.append("/Users/ludo/anaconda/envs/p2.6/lib/python2.6/site-packages/")
import json
import numpy
import numpy as np 
import scipy
from scipy import ndimage
from scipy import signal
from scipy import spatial
from scipy import misc

import numpy as np
from scipy.stats import cauchy
import matplotlib

matplotlib.use('Agg')
#matplotlib.use('tkAgg')

from matplotlib import pyplot as plt
import matplotlib.cm as cm
from scipy import fftpack
from scipy import stats
#from sklearn.cluster import DBSCAN
#from sklearn import metrics
#from sklearn.datasets.samples_generator import make_blobs
#from sklearn.preprocessing import StandardScaler
import numpy as np
#from sklearn.cluster import MeanShift, estimate_bandwidth
#from sklearn.datasets.samples_generator import make_blobs
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib.offsetbox import AnchoredOffsetbox
#from skimage.segmentation import random_walker
from skimage import measure
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14_95696CEA/plugins/ePMV/mgl64/MGLToolsPckgs")

from Volume.IO.volWriters import WriteCCP4
from Volume.Grid3D import Grid3DUC, Grid3DSI, Grid3DF
# -*- coding: utf-8 -*-
# -*- mode: python -*-
# Adapted from mpl_toolkits.axes_grid2
# LICENSE: Python Software Foundation (http://docs.python.org/license.html)

#Once in MatLab, Ripley's K function was calculated and subsequently normalized to L(r), 
#where r is the radius, using the equation L(r) = sqrt[K(r)/π]. 
#Plotting L(r) − r [called H(r)] versus r allowed the peak of maximum clustering 
#to be determined and represents a number between the radius and diameter of the clusters

#Statistical analysis of H(r) maximums was carried out in GraphPad Prism 5 using a one-way analysis 
#of variance (ANOVA) with a Bonferroni post hoc test. Values were considered significantly 
#different with a P value of ≤0.05.
W=4195 
def generate_hex_circle_packing(a, width):
    """Generate a domain of a given width filled with hexagonally packed
    circles.  The height will be determined so that the vertical 
    boundary condition is periodic.

    Arguments:
    a       particle radius
    width   domain width, in terms of particle radius

    Returns:
    x_list  list of x coordinates
    y_list  list of y coordinates
    x_size  width of domain (equal to argument width)
    y_size  height of domain
    """
    numParticles = 0

    x_list = []
    y_list = []
    y = a
    x = a
    rowNumber = 0
    # Create a row
    while y <= width*1.01:
        # Create circles in a row
        while x < width:
            x_list.append(x)
            x = x + 2*a
            y_list.append(y)
            numParticles = numParticles + 1
        y = y + a*sqrt(3.)
        rowNumber = rowNumber + 1
        if rowNumber%2 == 0:
            x = a
        else:
            x = 0
    x_size = width
    y_size = rowNumber*a*sqrt(3)

    return array(x_list), array(y_list), x_size, y_size

def plot_adsorbed_circles(adsorbed_x, adsorbed_y, radius, width, reference_indices=[]):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    # Plot each run
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for p in range(len(adsorbed_x)):
        if len(np.where(reference_indices == p)[0]) > 0:
            ax.add_patch(Circle((adsorbed_x[p], adsorbed_y[p]), radius, 
                edgecolor='black', facecolor='black'))
        else:
            ax.add_patch(Circle((adsorbed_x[p], adsorbed_y[p]), radius,
                edgecolor='black', facecolor='white'))

    ax.set_aspect(1.0)
    plt.axhline(y=0, color='k')
    plt.axhline(y=width, color='k')
    plt.axvline(x=0, color='k')
    plt.axvline(x=width, color='k')
    plt.axis([-0.1*width, width*1.1, -0.1*width, width*1.1])
    plt.xlabel("non-dimensional x")
    plt.ylabel("non-dimensional y")

    return ax
#xy = numpy.meshgrid(self.x,self.y,copy=False)
#locs = numpy.vstack(xyz).reshape(2,-1).T 
#dist = scipy.spatial.distance.pdist(np.sum(gauss_conv,axis=0))
def RipleysK(locs,dist,box,method=0,DIST=None):
    #http://www.colorado.edu/geography/class_homepages/geog_4023_s07/labs/lab10/RipleysK.m
    #% RipleysK: Calculate K statistic
    #% 
    #% K = RipleysK(locs,dist, box,method) calculates G, the K statistic at each 
    #% distance for the data with x-y coordinates given by locs, and the
    #% bounding rectangle given by box=[minX maxX minY maxY].
    #% If method=0, no edge correction is applied.
    #% If method=1, points are only used if they are at least h units away from
    #% the edge.
    #%
    #% Note: The L statistic may be calculated from the K statistic as follows: 
    #%   L = sqrt(K/pi)-h;
    #%   
    
#    if nargin<4 and method=1: return
    N,k = locs.shape;
#    if k~=2, error('locs must have two columns'); end
#    % rbox is distance to box
#    bigx = tile(x,(2,1))
    if DIST is None :
        DX = np.tile(locs[:,0],(0,N))-np.tile(locs[:,0].conj().transpose(),(N,0));
        DY = np.tile(locs[:,1],(0,N))-np.tile(locs[:,1].conj().transpose(),(N,0));
        DIST = np.sqrt(DX**2+DY**2);
        DIST = np.sort(DIST);
    
    if method==1:
        rbox = np.min(np.array([locs[:,0].conj().transpose()-box[0],box[1]-locs[:,0].conj().transpose(),locs[:,1].conj().transpose()-box[2], box[3]-locs[:,1].conj().transpose()]) );
        K = np.zeros(len(dist));
        for k in range(len(K)):
            I = find(rbox>dist[k]);
            if len(I)>0:
                K[k] = np.sum(np.sum(DIST[1:,I]<dist[k]))/len(I);
    elif method==0:
        K = np.zeros(len(dist),1);
        for k in range(len(K)):#for k inrange(_=1:length(K)
            K[k] = np.sum(np.sum(DIST[1:,:]<dist[k]))/N;
    
    lambda1 = N/((box[1]-box[0])*(box[3]-box[2]));
    K = K/lambda1;
    return K,(K/pi)**0.5

#d=scipy.spatial.distance.pdist(np.sum(gauss_conv,axis=0))
def ripleyK(d,A,n,t):
    #K(t) = lambda^-1*SUM(I(dij<t)/n)
    #K(t) = A*SUM(wij*I(i,j)/n**2)
    #lambda = n/A A is the area of the region containing all points
    #I indicator function 1 if its operand is true, 0 otherwise
    #t is the search radius
    #if homogenous K(s) = pi*s**2
    #L(t) = (K(t)/pi)**1/2
    #A common plot is a graph of t - \hat{L}(t) against t
    #which will approximately follow the horizontal zero-axis with constant 
    #dispersion if the data follow a homogeneous Poisson process.
    lambda1 = (n/A)**-1.0
    i = d < t
    I=np.array(i,dtype=int)/n
    Kt=lambda1*I
    Kl=(Kt/pi)**0.5
    return Kt,Kl
    
def PairCorrelationFunction_2D(x,y,S,rMax,dr):
    """Compute the two-dimensional pair correlation function, also known 
    as the radial distribution function, for a set of circular particles 
    contained in a square region of a plane.  This simple function finds 
    reference particles such that a circle of radius rMax drawn around the 
    particle will fit entirely within the square, eliminating the need to 
    compensate for edge effects.  If no such particles exist, an error is
    returned. Try a smaller rMax...or write some code to handle edge effects! ;) 
    
    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        S               length of each side of the square region of the plane
        rMax            outer diameter of largest annulus
        dr              increment for increasing radius of annulus

    Returns a tuple: (g, radii, interior_x, interior_y)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        annuli used to compute g(r)
        interior_x      x coordinates of reference particles
        interior_y      y coordinates of reference particles
    """
    from numpy import zeros, sqrt, where, pi, average, arange, histogram
    # Number of particles in ring/area of ring/number of reference particles/number density
    # area of ring = pi*(r_outer**2 - r_inner**2)

    # Find particles which are close enough to the box center that a circle of radius
    # rMax will not cross any edge of the box
    bools1 = x>1.1*rMax
    bools2 = x<(S-1.1*rMax)
    bools3 = y>rMax*1.1
    bools4 = y<(S-rMax*1.1)
    interior_indices, = where(bools1*bools2*bools3*bools4)
    num_interior_particles = len(interior_indices)

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a circle of radius rMax\
                will lie entirely within a square of side length S.  Decrease rMax\
                or increase the size of the square.")

    edges = arange(0., rMax+1.1*dr, dr)
    num_increments = len(edges)-1
    g = zeros([num_interior_particles, num_increments])
    radii = zeros(num_increments)
    numberDensity = len(x)/S**2

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = sqrt((x[index]-x)**2 + (y[index]-y)**2)
        d[index] = 2*rMax

        (result,bins) = histogram(d, bins=edges, normed=False)
        g[p,:] = result/numberDensity
        
    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1])/2.        
        rOuter = edges[i+1]
        rInner = edges[i]
        g_average[i] = average(g[:,i])/(pi*(rOuter**2 - rInner**2))

    return (g_average, radii, interior_indices)
####

def PairCorrelationFunction_3D(x,y,z,S,rMax,dr):
    """Compute the three-dimensional pair correlation function for a set of
    spherical particles contained in a cube with side length S.  This simple 
    function finds reference particles such that a sphere of radius rMax drawn
    around the particle will fit entirely within the cube, eliminating the need
    to compensate for edge effects.  If no such particles exist, an error is 
    returned.  Try a smaller rMax...or write some code to handle edge effects! ;) 
    
    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        z               an array of z positions of centers of particles
        S               length of each side of the cube in space
        rMax            outer diameter of largest spherical shell
        dr              increment for increasing radius of spherical shell

    Returns a tuple: (g, radii, interior_x, interior_y, interior_z)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        spherical shells used to compute g(r)
        interior_x      x coordinates of reference particles
        interior_y      y coordinates of reference particles
        interior_z      z coordinates of reference particles
    """
    from numpy import zeros, sqrt, where, pi, average, arange, histogram

    # Find particles which are close enough to the cube center that a sphere of radius
    # rMax will not cross any face of the cube
    bools1 = x>rMax
    bools2 = x<(S-rMax)
    bools3 = y>rMax
    bools4 = y<(S-rMax)
    bools5 = z>rMax
    bools6 = z<(S-rMax)

    interior_indices, = where(bools1*bools2*bools3*bools4*bools5*bools6)
    num_interior_particles = len(interior_indices)

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a sphere of radius rMax\
                will lie entirely within a cube of side length S.  Decrease rMax\
                or increase the size of the cube.")

    edges = arange(0., rMax+1.1*dr, dr)
    num_increments = len(edges)-1
    g = zeros([num_interior_particles, num_increments])
    radii = zeros(num_increments)
    numberDensity = len(x)/S**3

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = sqrt((x[index]-x)**2 + (y[index]-y)**2 + (z[index]-z)**2)
        d[index] = 2*rMax

        (result,bins) = histogram(d, bins=edges, normed=False)
        g[p,:] = result/numberDensity
        
    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1])/2.        
        rOuter = edges[i+1]
        rInner = edges[i]
        g_average[i] = average(g[:,i])/(4./3.*pi*(rOuter**3 - rInner**3))

    return (g_average, radii, x[interior_indices], y[interior_indices], z[interior_indices])
    # Number of particles in shell/total number of particles/volume of shell/number density
    # shell volume = 4/3*pi(r_outer**3-r_inner**3)
####


class fluoSim:
    def __init__(self,values=None,psf_type="gauss",boundingBox=4195,res=400.0,space=50.0 ):
        self.psf = None
        self.res = res
        self.values = values
        self.data = None
        self.makeGrid(values=values,boundingBox=boundingBox,res=res,space=space)
        self.psf_type =  psf_type       
        if psf_type == "gauss":
            self.setPSF_gauss()
        elif psf_type == "lorentz":
            self.setPSF_cauchy()
        if self.values is not None :
            self.setValues()
        self.shifted_cmap = shiftedColorMap(cm.get_cmap('hot', boundingBox), start=0, midpoint=0.5, stop=0.55,  name='shifted')
        self.extE=[-self.W,self.W,-self.W,self.W]
#        self.shifted_cmap = shiftedColorMap(cm.get_cmap('hot', 4000), start=-1.0, midpoint=0.5, stop=0.65,  name='shifted')
                    
    def makeGrid(self,values=None,boundingBox=2000,res=400.0,space=200.0 ):
        self.W=boundingBox
        self.boundingBox=np.array([[-self.W, -self.W, -self.W], [self.W, self.W, self.W]])
        self.x = np.linspace(self.boundingBox[0][0], self.boundingBox[1][0] , int((self.W*2.)/space))#*1.1547) gridspacing is already multiplied by 1.1547
        self.y = np.linspace(self.boundingBox[0][1], self.boundingBox[1][1] , int((self.W*2.)/space))#*1.1547)
        self.z = np.linspace(self.boundingBox[0][2], self.boundingBox[1][2] , int((self.W*2.)/space))#*1.1547)
#        self.x = np.arange(self.boundingBox[0][0], self.boundingBox[1][0] + space, space)#*1.1547) gridspacing is already multiplied by 1.1547
#        self.y = np.arange(self.boundingBox[0][1], self.boundingBox[1][1] + space, space)#*1.1547)
#        self.z = np.arange(self.boundingBox[0][2], self.boundingBox[1][2] + space, space)#*1.1547)
        self.nx = len(self.x) # sizes must be +1 or the right, top, and back edges don't get any points using this numpy.arange method
        self.ny = len(self.y)
        self.nz = len(self.z)
        self.space = space

    def setValues(self,values=None):
        if values is None :
            values = self.values
        if values is not None :
            zvalues=np.zeros(self.nx*self.ny*self.nz)
            self.data = zvalues.reshape((self.nx,self.ny,self.nz))
            ax,ay,az = values.transpose()
            ix = np.searchsorted(self.x, ax) - 1
            iy = np.searchsorted(self.y, ay) - 1
            iz = np.searchsorted(self.z, az) - 1
            self.data[ix,iy,iz]=1.0

    def setValuesAccurate(self,values=None):
        if values is None :
            values = self.values
        if values is not None :
            xyz = numpy.meshgrid(self.x,self.y,self.z,copy=False)
            grid = np.array(xyz).T
            tr=spatial.cKDTree(grid.reshape((self.nx*self.ny*self.nz,3)), leafsize=10)
            d,i=tr.query(values)
            zvalues=np.zeros(self.nx*self.ny*self.nz)
            zvalues[i]=1.0
            self.data = zvalues.reshape((self.nx,self.ny,self.nz))
#            ax,ay,az = values.transpose()
#            ix = np.searchsorted(self.x, ax) - 1
#            iy = np.searchsorted(self.y, ay) - 1
#            iz = np.searchsorted(self.z, az) - 1
#            self.data[ix,iy,iz]=1.0
#map_coords = np.broadcast_arrays(np.arange(a.shape[0])[:, None], y, x)
#   .....: ndimage.map_coordinates(a, map_coords, order=1)    
            
    def setPSF_gauss(self,res=None):
        if res == None :
            res = self.res
        gaussian_blur_sigma =  (res/self.space)/2.35482
        X, Y , Z= np.ogrid[-self.nx/2:self.nx/2, -self.ny/2:self.ny/2, -self.nz/2:self.nz/2]#-27:27
        self.psf = stats.norm.pdf(np.sqrt((X*self.space)**2 + (Y*self.space)**2 +(Z*self.space)**2), 0, gaussian_blur_sigma*self.space)#norm is gaussian
        return self.psf        

    def lorentzian(self,x,*p) :
        # A lorentzian peak with:
        #   Constant Background          : p[0] #0
        #   Peak height above background : p[1] #1 amplitude is 1/(pi*gamma)
        #   Central value                : p[2] #x0 = 0
        #   Full Width at Half Maximum   : p[3] #gamma/2.?
        return p[0]+(p[1]/(numpy.pi*p[3]))/(1.0+((x-p[2])/p[3])**2)
        
    def setPSF_cauchy(self,res=None):
        if res == None :
            res = self.res       
        cauchy_sigma =  (res/self.space)/2.0
        X, Y , Z= np.ogrid[-self.nx/2:self.nx/2, -self.ny/2:self.ny/2, -self.nz/2:self.nz/2]#-27:27
#        self.psf = stats.cauchy.pdf(np.sqrt((X*self.space)**2 + (Y*self.space)**2 +(Z*self.space)**2), 0, cauchy_sigma*self.space)#norm is gaussian
        self.psf = lorentzian(np.sqrt((X*self.space)**2 + (Y*self.space)**2 +(Z*self.space)**2), 0, 600,0,cauchy_sigma*self.space)#norm is gaussian
        return self.psf

    def setPSF_standardcauchy(self,res=None):
        self.psf = np.random.standard_cauchy((self.nx,self.ny,self.nz))
#        self.psf = self.psf[(self.psf>-25) & (self.psf<25)]  # truncate distribution so it plots well
        return self.psf
        
    def convolve(self,star, psf):
        star_fft = fftpack.fftshift(fftpack.fftn(star))
        psf_fft = fftpack.fftshift(fftpack.fftn(psf))
        return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft*psf_fft)))
    
    def deconvolve(self,star, psf):
        star_fft = fftpack.fftshift(fftpack.fftn(star))
        psf_fft = fftpack.fftshift(fftpack.fftn(psf))
        return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft/psf_fft)))
    
    def getConvolution(self):
#        self.data_conv = np.real(self.convolve(self.data,self.psf))
        self.data_conv=ndimage.filters.gaussian_filter(self.data,(self.res/self.space)/2.35482)
        return self.data_conv

    def getCluster1DXYZ(self,limit=None,ecc=0.5,uplimit=30,saveimgeFile=None):
        locis=[0,0,0]
        nps = [0,0,0]
        limits=[0,0,0]
        for i in range(3):
            loci,npx,limit = self.getCluster1D(limit=limit,ecc=ecc,uplimit=uplimit,ax=i)
            locis[i] = loci
            nps[i] = npx
            limits[i] = limit
            if saveimgeFile != None:
                self.saveFluoImage(np.sum(self.data_conv,axis=i),saveimgeFile+"_"+str(i)+".jpg")
        return locis,nps,limits

    def getCluster1D_contour(self,data1d=None,axe=None,outputfile=None,ax=0,level=1):
        if data1d is None :
            data1d=np.sum(self.data_conv,axis=ax)
        data1d[data1d < data1d.mean()]=0
        if axe is not None :
            ctr=axe.contour(data1d,level,hold='on',colors='g',linewidths=2)
        else :
            ctr=plt.contour(data1d,level,hold='on',colors='g',linewidths=2)
        n = len(ctr.collections[-1].properties()['paths'])
        if n > 3 :
            n=3
        if axe is not None :
            axe.text(2, 3, "loci nb is "+str(n),color='white')
        return n

    def analyzeGauss(self,data1d,limit=None,ax=0):
        if data1d is None :
            data1d=np.sum(self.data_conv,axis=ax)
        if limit is None:
            limit=self.data_conv.max()+data1d.std()
        data1d[data1d < limit]=0#removebg
        labeled, nr_objects = ndimage.label(data1d,structure=np.ones((3,3), dtype="bool8"))#,structure=s)# > data.max()/2.0);
        m=measure.regionprops(labeled, properties=['Area', 'Perimeter','EquivDiameter','FilledArea','Eccentricity','Extent','Solidity'])
        return  m      
        
    def getCluster1D(self,data1d=None,limit=None,ecc=0.5,uplimit=30,axe=None,
                     outputfile=None,ax=0,hide=False):
        #data=data/data.max()
        if data1d is None :
            data1d=np.sum(self.data_conv,axis=ax)
        if limit is None:
            limit=self.data_conv.max()+data1d.std()
#            limit=self.data_conv.mean()+self.data_conv.std()
#        print limit
        data1d[data1d < limit]=0#removebg
    #    data[data > 0.55*data.max() ] = 0.55
        npx=len(np.nonzero(data1d > 0)[0])
        #erosion
        #a=ndimage.binary_opening(data1d).astype(np.int)
        labeled, nr_objects = ndimage.label(data1d,structure=np.ones((3,3), dtype="bool8"))#,structure=s)# > data.max()/2.0);
        m=measure.regionprops(labeled, properties=['Area', 'Perimeter','EquivDiameter','FilledArea','Eccentricity'])
        #1 loci circular with npix < 32 -> foci1
        #2 loci circular each withnpix < 20 -> foci2
        #else loci 3
        loci=0
        if nr_objects == 1: 
            score = m[0]['Eccentricity'] #0 circle 1 extend ellipse
            if score < 0.5 and m[0]['FilledArea'] < uplimit:
                loci=1
            elif m[0]['FilledArea'] < uplimit:#up to 30?
                loci=1
            else :
                loci=3
        elif nr_objects == 2:
            if m[0]['Eccentricity'] < ecc :
                #first blob spherical
                if m[1]['FilledArea'] < uplimit :
                    if m[1]['FilledArea'] > 2 and  m[0]['FilledArea'] > 2:
                        loci = 2
                    else :
                        loci = 1
                elif m[1]['Eccentricity'] < ecc:
                    loci = 2
                else :
                    loci=3 
            else :
                if m[0]['FilledArea'] < uplimit :
                    if m[1]['FilledArea'] < uplimit :
                        if m[1]['FilledArea'] > 2 and  m[0]['FilledArea'] > 2:
                            loci = 2
                        else :
                            loci = 1
                    elif m[1]['Eccentricity'] < ecc:
                        loci = 2
                    else :
                        loci=3 
                else :
                    loci = 3               
        elif nr_objects >= 3:
    #        print "nr_objects 3"
            loci=3
        else :
            print "?",nr_objects,loci
    #    print "LOCI ",loci
        if type(axe) == int and axe == 1 :
            plt.close('all')            
            f, axe = plt.subplots()
#            f.set_size_inches(12, 10)
        if axe is not None :
            axe.imshow(labeled);
            #axe.text(2, 3, "loci nb is "+str(loci)+" "+str(nr_objects)+" features "+str(npx)+" px ",color='white')
            if hide : 
                axe.xaxis.set_visible(False)
                axe.yaxis.set_visible(False)
            print "loci nb is "+str(loci)+" "+str(nr_objects)+" features "+str(npx)+" px "
        return loci,npx,limit 
    
    def writeCCP4(self,filename, data=None):  
        if data is None :
            data = self.data_conv
        w=WriteCCP4()       
        h={}
        maskGrid = Grid3DF( np.array(data,'f'), [0,0,0], [self.space,self.space,self.space] , h)
        maskGrid.stepSize = [self.space,self.space,self.space]
        maskGrid.origin = (self.boundingBox[1]-self.boundingBox[0])/2.0
        #maskGrid.normalize()
        h['amin'], h['amax'],h['amean'],h['arms']= maskGrid.stats()
        w.write(filename,maskGrid)

    def writeDATA(self,filename, data=None):
        if data is None :
            data = self.data_conv
        np.save(filename,data)

    def readDATA(self,filename):
        self.data_conv=np.load(filename)

    def addSBar(self,ax,pos=3,scaleb=1000,hideax=False):
        bar = AnchoredSizeBar(ax.transData, scaleb, '',loc=pos,pad=1, sep=1, borderpad=0.01, frameon=False)#,prop=None)
        rect = bar.size_bar._children[0]
        rect.set_ec("white")
        rect.set_height(200)
        rect.set_fc('white')
        ax.add_artist(bar)
        if hideax:
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
    
    def plotPSF(self,psf=None,ax=0):
        plt.close('all')
        extE=[-self.W,self.W,-self.W,self.W]
        f, axe = plt.subplots()
        #shifted_cmap = shiftedColorMap(cm.get_cmap('hot', 30), start=0, midpoint=0.5, stop=0.55,  name='shifted')
        if psf is None :
            psf = self.psf
        axe.imshow(numpy.sum(psf,axis=ax),extent=extE)#cmap=shifted_cmap,interpolation='none',
        self.addSBar(axe)#,pos=3,scaleb=1000)

    def complexePlot(self,data,ax):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        plt.close('all')
        extE=[-self.W,self.W,-self.W,self.W]
        axScatter = plt.subplot(111)
        #axScatter.scatter(x, y)
        axScatter.imshow(numpy.sum(data,axis=ax),extent=extE)
        axScatter.set_aspect(1.)        
        # create new axes on the right and on the top of the current axes.
        divider = make_axes_locatable(axScatter)
        axHistx = divider.append_axes("top", size=1.2, pad=0.1, sharex=axScatter)
        axHisty = divider.append_axes("right", size=1.2, pad=0.1, sharey=axScatter)
        # the scatter plot:
        # histograms
        size =self.nx/2.0
        sizey=self.ny/2.0
        x, y = mgrid[-size:size+1, -sizey:sizey+1]
        g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
        #bins = np.arange(-lim, lim + binwidth, binwidth)
        axHistx.plot(g)      
        axHisty.plot(g)      
#        axHistx.hist(x, bins=bins)
#        axHisty.hist(y, bins=bins, orientation='horizontal')

    def plotDistrib(self,distrib,values=None):
        plt.close('all')
        x = np.linspace(-self.W, self.W, len(distrib))
        f, axe = plt.subplots()
        axe.plot(x,distrib, ls='-', color='red',label=r'$\mu=%i,\ \gamma=%.1f$' % (0, 400/2))
        plt.show()

    def saveFluoImage(self,data,filename):
        plt.close('all')            
        f, axe = plt.subplots()
#        ticks=np.linspace(-self.W,self.W,5)
        axe.imshow(data,cmap=self.shifted_cmap,interpolation='none',extent=self.extE)
#        axe.set_xticks(ticks)
#        axe.set_yticks(ticks)
        self.addSBar(axe,hideax=True)  
        plt.savefig(filename)
                
    def plotData(self,values=None,ax=0,axe=None,hide=False):
        if values is None :
            return
        if axe is None :
            plt.close('all')            
            f, axe = plt.subplots()
            f.set_size_inches(12, 10)
        data=np.sum(values,axis=ax)
        ticks=np.linspace(-self.W,self.W,5)
        axe.imshow(data,cmap=self.shifted_cmap,interpolation='none',extent=self.extE)
        axe.set_xticks(ticks)
        axe.set_yticks(ticks)
        if hide : 
            axe.xaxis.set_visible(False)
            axe.yaxis.set_visible(False)
        self.addSBar(axe)
        
    def setAtScale(self,space):
        sc=( ((self.W*2.)/space)/((self.W*2.)/self.space) )
        r=scipy.ndimage.interpolation.zoom(self.data_conv,sc)#50->200
        return r

    def setDataAtScale(self,space):
        #is the interpolation translate ?
        sc=( ((self.W*2.)/space)/((self.W*2.)/self.space) )
        self.data=scipy.ndimage.interpolation.zoom(self.data,sc)#50->200
#        self.makeGrid(boundingBox=self.W,res=400.0,space=space )
#        self.setPSF_gauss()
        self.data [self.data < 0.01] = 0.0                     
        
    def plotRotXYZ(self,values=None,values_conv=None,pos=3,limit=None,ecc=0.5,uplimit=25,
                   filename=None):
        plt.close('all')
        if values is None :
            values = self.data
        if values_conv is None :
            values_conv = self.data_conv
        f, axes = plt.subplots(3,3)
        locis=[]
        for i in range(3):
            self.plotData(values=values,ax=i,axe=axes[0,i])
            self.plotData(values=values_conv,ax=i,axe=axes[1,i])
            loci,npx,limits=self.getCluster1D(limit=limit,ecc=ecc,uplimit=uplimit,axe=axes[2,i],ax=i)
#            n=self.getCluster1D_contour(axe=axes[2,i],ax=i)
        if filename is not None :
            f.set_size_inches(12, 10);
            plt.savefig(filename,dpi=300)
        return f

    def parseAPresults(self,filename,recipe_name,ingr_name):
        with open(filename, 'r') as fp :#doesnt work with symbol link ?
            ap_results=json.load(fp)
        #ap_results['HIV1_envelope_Pack_106_0_2_0_surfaceRecipe']['HIV1_envelope_Pack_106_0_2_0_surf__HIV1_ENV_4nco_0_1_1']['results']
        ingr_res = ap_results[recipe_name][ingr_name]['results']
        ingr_pos = np.array(ingr_res).transpose()[0]
        return np.array(ingr_pos.tolist(),dtype='f')

#from https://gist.github.com/dmeliza/3251476 
class AnchoredScaleBar(AnchoredOffsetbox):
    def __init__(self, transform, sizex=0, sizey=0, labelx=None, labely=None, loc=4,
                 pad=0.1, borderpad=0.1, sep=2, prop=None, **kwargs):
        """
        Draw a horizontal and/or vertical  bar with the size in data coordinate
        of the give axes. A label will be drawn underneath (center-aligned).
 
        - transform : the coordinate frame (typically axes.transData)
        - sizex,sizey : width of x,y bar, in data units. 0 to omit
        - labelx,labely : labels for x,y bars; None to omit
        - loc : position in containing axes
        - pad, borderpad : padding, in fraction of the legend font size (or prop)
        - sep : separation between labels and bars in points.
        - **kwargs : additional arguments passed to base class constructor
        """
        from matplotlib.patches import Rectangle
        from matplotlib.offsetbox import AuxTransformBox, VPacker, HPacker, TextArea, DrawingArea
        bars = AuxTransformBox(transform)
        if sizex:
            bars.add_artist(Rectangle((0,0), sizex, 0, fc="none"))
        if sizey:
            bars.add_artist(Rectangle((0,0), 0, sizey, fc="none"))
 
        if sizex and labelx:
            bars = VPacker(children=[bars, TextArea(labelx, minimumdescent=False)],
                           align="center", pad=0, sep=sep)
        if sizey and labely:
            bars = HPacker(children=[TextArea(labely), bars],
                            align="center", pad=0, sep=sep)
 
        AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                   child=bars, prop=prop, frameon=False, **kwargs)


def gaussian(x,*p) :
    # A gaussian peak with:
    #   Constant Background          : p[0]
    #   Peak height above background : p[1]
    #   Central value                : p[2]
    #   Standard deviation           : p[3]
    return p[0]+p[1]*numpy.exp(-1*(x-p[2])**2/(2*p[3]**2))
#self.nx=100,self.ny=100,self.nz=100
#X, Y , Z= np.ogrid[-self.nx/2:self.nx/2, -self.ny/2:self.ny/2, -self.nz/2:self.nz/2]#-27:27
#C=np.sqrt((X*self.space)**2 + (Y*self.space)**2 +(Z*self.space)**2)
#
def lorentzianPSF(C,A,B,FWHM) :
    # A lorentzian peak with:
    #   Constant Background          : B #0
    #   Amplitude Peak height above background : A #1
    #   Central value                : p[2] #x0 = 0, y0 =0
    #   Full Width at Half Maximum   : p[3] #gamma/2.?
    (1/(numpy.pi*FWHM/2.0))*(1+(sqrt(2)-1)*(C/FWHM))+B
    (1.0/(numpy.pi*FWHM/2.0))/(1.0+(C/FWHM/2.0)**2)
    
    return p[0]+(p[1]/numpy.pi)/(1.0+((x-p[2])/p[3])**2)
        
def lorentzian(x,*p) :
    # A lorentzian peak with:
    #   Constant Background          : p[0] #0
    #   Peak height above background : p[1] #1 amplitude is 1/(pi*gamma)
    #   Central value                : p[2] #x0 = 0
    #   Full Width at Half Maximum   : p[3] #gamma/2.?
    return p[0]+(p[1]/(numpy.pi*p[3]))/(1.0+((x-p[2])/p[3])**2)
    
def linear(x,*p) :
    # A linear fit with:
    #   Intercept                    : p[0]
    #   Slope                        : p[1]
    return p[0]+p[1]*x
def power(x,*p) :
    # A power law fit with:
    #   Normalization                : p[0]
    #   Offset                       : p[1]
    #   Constant                     : p[3]
    return p[0]*(x-p[1])**p[2]+p[3]

    
def add_scalebar(ax, matchx=True, matchy=True, hidex=True, hidey=True, **kwargs):
    """ Add scalebars to axes
 
    Adds a set of scale bars to *ax*, matching the size to the ticks of the plot
    and optionally hiding the x and y axes
 
    - ax : the axis to attach ticks to
    - matchx,matchy : if True, set size of scale bars to spacing between ticks
                    if False, size should be set using sizex and sizey params
    - hidex,hidey : if True, hide x-axis and y-axis of parent
    - **kwargs : additional arguments passed to AnchoredScaleBars
 
    Returns created scalebar object
    """
    def f(axis):
        l = axis.get_majorticklocs()
        return len(l)>1 and (l[1] - l[0])
    
    if matchx:
        kwargs['sizex'] = f(ax.xaxis)
        kwargs['labelx'] = str(kwargs['sizex'])
    if matchy:
        kwargs['sizey'] = f(ax.yaxis)
        kwargs['labely'] = str(kwargs['sizey'])
        
    sb = AnchoredScaleBar(ax.transData, **kwargs)
    ax.add_artist(sb)
 
    if hidex : ax.xaxis.set_visible(False)
    if hidey : ax.yaxis.set_visible(False)
 
    return sb
    
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap
    
def convolve(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft*psf_fft)))

def deconvolve(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft/psf_fft)))

#from scipy.cluster.vq import kmeans,vq
#centroids,_ = kmeans(data,2)
#centroids,_ = kmeans(data,3)
#idx,_ = vq(data,centroids)
#d = sch.distance.pdist(X) distance between value
#distance between position?
#extract position from data > 0.0
#linkage_matrix = linkage(D, method='complete')
#clusters = fcluster(linkage_matrix, threshold)
#len(numpy.unique(scipy.cluster.hierarchy.fcluster(linkage_matrix, threshold, criterion='distance')))
#b = ((a[:,:,0] == 255.0) * (a[:,:,1] == 0) * (a[:,:,2] == 0))*1
#labeled_array, num_features = scipy.ndimage.measurements.label(b.astype('Int8'))
def cluster2(data,axe=None):
    markers = np.zeros(data.shape, dtype=np.uint)
    markers[data < 0.1] = 1
    markers[data > 1.3] = 2
    struct=np.ones((3,3), dtype="bool8")
    labeled, nr_objects = ndimage.label(markers,structure=struct)#,structure=s)
    print nr_objects,
    X=data
    bandwidth = estimate_bandwidth(data, quantile=0.2)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(X)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_
    labels_unique = np.unique(labels)
    n_clusters_ = len(np.nonzero(labels_unique))
    labels = random_walker(data, markers, beta=10, mode='bf')
    n,bins=np.histogram(labels.flatten(), bins=[1,2,3], density=False)
    #a circular blob up to 25px is a 1 foci (16.266666666666666+-4)
    if n[1] < 18 :
#        print n,n[1],1
        loci = 1
    elif n[1] < 25 :
#        print n,n[1],2
        loci = 2
    elif n[1] >=25:
#        print n,n[1],3
        loci = 3
    if axe is not None :
        axe.imshow(labels);
        axe.text(2, 3, "loci nb is "+str(loci)+" "+str(n[1])+" px"+" "+str(nr_objects)+" "+str(n_clusters_),color='white')
    return loci,n[1]
    
def cluster (data,limit,nbup,nbdown,axe=None,outputfile=None):
    #1 loci is around 20px
    d={1:"loci1",2:"loci2",3:"loci3M"}
    #data=data/data.max()
    data[data < limit]=0#removebg
#    data[data > 0.55*data.max() ] = 0.55
    npx=len(numpy.nonzero(data > 0)[0])
#    s = [[1,1,0],
#         [1,1,1],
#         [1,1,0]]
#    struct=np.ones((3,3), dtype="bool8")
    labeled, nr_objects = ndimage.label(data,structure=np.ones((3,3), dtype="bool8"))#,structure=s)# > data.max()/2.0);
    m=measure.regionprops(labeled, properties=['Area', 'Perimeter','EquivDiameter','FilledArea','Eccentricity'])
    #1 loci circular with npix < 32 -> foci1
    #2 loci circular each withnpix < 20 -> foci2
    #else loci 3
    loci=0
    if nr_objects == 1: 
        score = m[0]['Eccentricity'] #0 circle 1 extend ellipse
#        score=m[0]['Area']/m[0]['Perimeter']
        #compare to perfect circle?
#        pscore=m[0]['EquivDiameter']/4.0
#        delta = abs(pscore - score)
#        print "nr_objects 1",score,m[0]['FilledArea']
        if score < 0.5 and m[0]['FilledArea'] < 25:
            loci=1
        elif m[0]['FilledArea'] < 25:#up to 30?
            loci=1
        else :
            loci=3
    elif nr_objects == 2:
#        score1=m[0]['Area']/m[0]['Perimeter']
#        pscore1=m[0]['EquivDiameter']/4.0
#        delta1 = abs(pscore1 - score1)
#        score2=m[1]['Area']/m[1]['Perimeter']
#        pscore2=m[1]['EquivDiameter']/4.0
#        delta2 = abs(pscore2 - score2)
#        print "nr_objects 2",  m[0]['Eccentricity'],m[0]['Area'],  m[1]['Eccentricity'],m[1]['Area']
        if m[0]['Eccentricity'] < 0.5 :
            #first blob spherical
            if m[1]['FilledArea'] < 20 :
                loci = 2
            elif m[1]['Eccentricity'] < 0.5:
                loci = 2
            else :
                loci=3 
        else :
            if m[0]['FilledArea'] < 20 :
                if m[1]['FilledArea'] < 20 :
                    loci = 2
                elif m[1]['Eccentricity'] < 0.5:
                    loci = 2
                else :
                    loci=3 
            else :
                loci = 3               
#        if m[0]['FilledArea'] < 30 and m[1]['FilledArea'] < 30 :
#            loci =2 
#        else :
#            loci =3
#        if delta1 < 0.1 and delta2 < 0.1 and npx < 30:
#            loci = 2
#        else :
#            loci = 3
    elif nr_objects >= 3:
#        print "nr_objects 3"
        loci=3
    else :
        print "?",nr_objects
#    print "LOCI ",loci
    if axe is not None :
        axe.imshow(labeled);
        axe.text(2, 3, "loci nb is "+str(loci)+" "+str(nr_objects)+" px",color='white')
    return loci,npx 
    maxp=labeled.shape[0]*labeled.shape[1]
    if nr_objects > 3 : nr_objects=3
    #return labeled,nr_objects
    if nr_objects == 1 :
        b = range(1,nr_objects+2)
    else :
        b = range(1,nr_objects+1)
    b = range(1,nr_objects+2)
    n,bins=np.histogram(labeled.flatten(), bins=b, density=False)
    r=np.nonzero(n > nbdown)[0]
    loci = len(r)
    if loci == 1 :
        if n[0] > nbup:# or npx > 20 :
            loci = 3             
    if loci == 3 and npx < 20 :
        loci = len(r)
    l=0
    if npx < 50 and npx >= 30: l=3
    elif npx < 30 and npx >= 20: l=2    
    elif npx < 20 : l=1
    print "loci nb is ",loci,"nb is ", len(r),n,bins,data.max()/2.0,npx,l,nr_objects
    if axe is not None :
        axe.imshow(labeled);
        axe.text(2, 3, "loci nb is "+str(loci)+" "+str(npx)+" px",color='white')
    if outputfile is not None :
        f=open(outputfile,'w')
        f.write("%d\n" % (loci) )
        f.close()#    plt.show()
    
    return labeled,loci,npx

def convoluteGauss(filename,boundingBox=[[-693.41299629211426, -651.68601989746094, -735.077], 
                                         [653.65499687194824, 664.00502014160156, 694.65200000000004]]):
    #w=WriteCCP4()
    a=numpy.loadtxt(filename+".txt",delimiter=",")
    space=200.0
    boundingBox=numpy.array([[-W, -W, -W], [W, W, W]])
    center = (boundingBox[1]-boundingBox[0])/2.0
    x = numpy.arange(boundingBox[0][0], boundingBox[1][0] + space, space)#*1.1547) gridspacing is already multiplied by 1.1547
    y = numpy.arange(boundingBox[0][1], boundingBox[1][1] + space, space)#*1.1547)
    z = numpy.arange(boundingBox[0][2], boundingBox[1][2] + space, space)#*1.1547)
    nx = len(x) # sizes must be +1 or the right, top, and back edges don't get any points using this numpy.arange method
    ny = len(y)
    nz = len(z)
    zvalues=numpy.zeros(nx*ny*nz)
    values = zvalues.reshape((nx,ny,nz))
    ax,ay,az = a.transpose()
    ix = np.searchsorted(x, ax) - 1
    iy = np.searchsorted(y, ay) - 1
    iz = np.searchsorted(z, az) - 1
    #square or power 
    gaussian_blur_sigma =  (400.0/space)/2.35482#16.986436330590024/10.0   0.42lambdan or 1.222lamdan
    cauchy_sigma = (400.0/space)/2.0
    #square or power 
    values[ix,iy,iz]=1.0
#    gauss_conv=ndimage.filters.gaussian_filter(values,gaussian_blur_sigma)
    X, Y , Z= np.ogrid[-nx/2:nx/2, -ny/2:ny/2, -nz/2:nz/2]#-27:27
    psf_gauss = stats.norm.pdf(np.sqrt((X*space)**2 + (Y*space)**2 +(Z*space)**2), 0, gaussian_blur_sigma*space)#norm is gaussian
    psf_cauchy = []#stats.cauchy.pdf(np.sqrt((X*space)**2 + (Y*space)**2 +(Z*space)**2), 0, (cauchy_sigma*space)/2.0)#norm is gaussian
    gauss_conv = np.real(convolve(values,psf_gauss/psf_gauss.max()))
    cauchy_conv = []#np.real(convolve(values,psf_cauchy/psf_cauchy.max()))
#    h={}
#    maskGrid = Grid3DF( numpy.array(np.real(cauchy_conv),'f'), [0,0,0], [space,space,space] , h)
#    maskGrid.stepSize = [space,space,space]
#    maskGrid.origin = (boundingBox[1]-boundingBox[0])/2.0
#    #maskGrid.normalize()
#    h['amin'], h['amax'],h['amean'],h['arms']= maskGrid.stats()
#    w.write(filename,maskGrid)
    return values,psf_gauss,psf_cauchy,gauss_conv,cauchy_conv

def addSBar(ax,pos=3):
    bar = AnchoredSizeBar(ax.transData, 1000, '',loc=pos,pad=1, sep=1, borderpad=0.01, frameon=False)#,prop=None)
    rect = bar.size_bar._children[0]
    rect.set_ec("white")
    rect.set_height(200)
    rect.set_fc('white')
    ax.add_artist(bar)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

def plotRot(values,gauss_conv,limit=0.05,nbup=40,nbdown=0,pos=3):
    shifted_cmap = shiftedColorMap(cm.get_cmap('hot', 30), start=0, midpoint=0.5, stop=0.55,  name='shifted')
    extE=[-W,W,-W,W]
    plt.close('all')
    f, axes = plt.subplots(3,3)
    for i in range(3):
        vdata=numpy.sum(values,axis=i)
        axes[0,i].imshow(vdata,cmap=shifted_cmap,interpolation='none',extent=extE)
        addSBar(axes[0,i],pos=pos)  
        data=numpy.sum(np.real(gauss_conv),axis=i)
        axes[1,i].imshow(data,cmap=shifted_cmap,interpolation='none',extent=extE)
        addSBar(axes[0,i],pos=pos)  
        fmax=np.real(gauss_conv).max()+data.std()
        loci,npx=cluster (data,fmax,30,2,axe=axes[2,i],outputfile=None)
#        loci,npx=cluster2(data,axe=axes[1,i])
#        labeled,loci=cluster (data,limit,nbup,nbdown,axe=axes[1,i],outputfile=None)
#    return axes
       
def plot (values,psf_gauss,psf_cauchy,gauss_conv,cauchy_conv,ax=1):
    shifted_cmap = cm.get_cmap('hot', 30)#shiftedColorMap(cm.get_cmap('hot', 30), start=0, midpoint=0.5, stop=0.55,  name='shifted')
    plt.close('all')
    extE=[-W,W,-W,W]
    f, axes = plt.subplots(4,2)
    axes[0,0].imshow(numpy.sum(values,axis=ax),cmap=shifted_cmap,interpolation='none',extent=extE)
    axes[1,0].imshow(numpy.sum(psf_gauss,axis=ax),cmap=shifted_cmap,interpolation='none',extent=extE)
    axes[1,1].imshow(numpy.sum(psf_cauchy,axis=ax),cmap=shifted_cmap,interpolation='none', extent=extE)
    #axes[1,1].imshow(psf_cauchy_1 )
    axes[2,0].imshow(numpy.sum(np.real(gauss_conv),axis=ax),cmap=shifted_cmap,interpolation='none',extent=extE)
    addSBar(axes[2,0])  
    axes[2,1].imshow(numpy.sum(np.real(cauchy_conv),axis=ax),cmap=shifted_cmap,interpolation='none',extent=extE)
    addSBar(axes[2,1]) 
    #scatter
    x = np.linspace(-W, W, values.shape[0])
    dist1=stats.cauchy(0,400/2)
    dist2=stats.norm(0,400/2.35482)
    d=dist1.pdf(x)
    axes[0,1].plot(x, d/d.max(), ls='-', color='red',
                 label=r'$\mu=%i,\ \gamma=%.1f$' % (0, 400/2))
    d=dist2.pdf(x)
    axes[0,1].plot(x, d/d.max(), ls=':', color='blue',
                 label=r'$\mu=%i,\ \gamma=%.1f$' % (0, 400/2.35482))
#    data=numpy.sum(np.real(gauss_conv),axis=1)
#    labeled,loci=doit (data,axes[3,0],0.3,20,5)
#    data=numpy.sum(np.real(cauchy_conv),axis=1)
#    labeled,loci=doit (data,axes[3,1],0.3,20,5)
    return axes

def countLoci(x):
    #stacked bar
    #error is sqrt(n[2])/float(n.sum()) Rare events - Poisson distribution
    #x = np.loadtxt("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/loci_output/loci_res.txt")
    n, bins = np.histogram(x, 3, normed=0)
    plt.close("all")
    ind = np.arange(len(n))
#    rects1 = plt.bar(ind, n/float(n.sum()), 1, color='r', yerr=0)
    v=n/float(n.sum())
    loci3 = plt.bar([0], [v[2]], 0.35, color='grey', yerr=sqrt(n[2])/float(n.sum()))
    loci2 = plt.bar([0], [v[1]], 0.35, color='white',bottom=v[2], yerr=sqrt(n[1])/float(n.sum()))
    loci1 = plt.bar([0], [v[0]], 0.35, color='black', yerr=sqrt(n[0])/float(n.sum()),bottom=v[2]+v[1])
    plt.ylabel('% foci ENV')
##    plt.title('Scores by group and gender')
    plt.xticks(ind+0.35/2., ('random','') )
#    plt.yticks(np.arange(0,81,10))
    plt.legend( (loci1[0], loci2[0], loci3[0]), ('foci1', 'foci2', 'foci3') )
    plt.show()
#    ax.xaxis.set_visible(False)
#    ax.yaxis.set_visible(False)
    return n,bins,n/float(n.sum())

def autolabel(rects,ax):
    #from http://matplotlib.org/examples/api/barchart_demo.html
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., height/2.0, '%d'%int(height),
                ha='center', va='bottom')

def autolabelyerr(ax,rects,err=None):
    # attach some text labels
    for i,rect in enumerate(rects):
        height = rect.get_height()
        v='%.2f'%height
        y=0.5*height
        if err is not None :
            v='%.2f'%err[i]
            y=1.05*height
        ax.text(rect.get_x()+rect.get_width()/2., y, v,
                ha='center', va='bottom')
                
def autolabels(loci1,loci2,loci3,ax,yerr1,yerr2,yerr3):
    #from http://matplotlib.org/examples/api/barchart_demo.html
    # attach some text labels
    for i in range(len(loci1)):#rects:
        rect1 = loci1[i]
        rect2 = loci2[i]
        rect3 = loci3[i]
        height1 = rect1.get_height()
        height2 = rect2.get_height()
        height3 = rect3.get_height()
        ax.text(rect1.get_x()+rect1.get_width()/2., height1/2.0, '%2.1f'% (height1*100.0),
                ha='center', va='bottom',color='black')
        ax.text(rect2.get_x()+rect2.get_width()/2., height2/2.0+height1, '%2.1f'% (height2*100.0),
                ha='center', va='bottom',color='black')
        ax.text(rect3.get_x()+rect2.get_width()/2., height3/2.0+height1+height2, '%2.1f'% (height3*100.0),
                ha='center', va='bottom',color='white')
        ax.text(rect1.get_x()+rect1.get_width()/2., 1.01*height1, '%2.1f'% (yerr1[i]*100.0),
                ha='center', va='bottom',color='black')
        ax.text(rect2.get_x()+rect2.get_width()/2., 1.01*(height2+height1), '%2.1f'% (yerr2[i]*100.0),
                ha='center', va='bottom',color='white')
        ax.text(rect3.get_x()+rect2.get_width()/2., 1.01*(height3+height1+height2), '%2.1f'% (yerr3[i]*100.0),
                ha='center', va='bottom',color='black')

def chisquare(o,e):
    return (o-e)**2/e
#The corresponding probability is 0.05<P<0.02. This is smaller than the conventionally accepted significance level of 0.05 or 5%, so the null hypothesis that the two distributions are the same is rejected.
#Df	0.5	    0.10	0.05   0.02	0.01	0.001
#1	0.455	2.706	3.841	5.412	6.635	10.827
#2	1.386	4.605	5.991	7.824	9.210	13.815
#3	2.366	6.251	7.815	9.837	11.345	16.268
#4	3.357	7.779	9.488	11.668	13.277	18.465
#5	4.351	9.236	11.070	13.388	15.086	20.517
def countLocis(x,nb=4,xlabels=('RAND-noMA', 'RAND-offMA', 'RAND-inMA', 'closeMA', 'closeMA-ENV','close-ENV','PR-')):
    #stacked bar
    #x = np.loadtxt("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/loci_output/loci_res.txt")
    data1=[]
    data2=[]
    data3=[]
    yerr1=[]
    yerr2=[]
    yerr3=[]
    nn=[]
    ind=np.arange(nb+1)
    for k in range(nb) :
#        if not len(run) :
#            n=np.zeros(len(b)-1)
#        else :
#        
        n, bins = np.histogram(x[k], 3, normed=0)
        if np.sum(n) == 0 :
            v=[0,0,0]
        else :
            v=n/float(n.sum())
        data3.append(v[0])
        data2.append(v[1])
        data1.append(v[2])
        yerr1.append(sqrt(n[0])/float(n.sum()))
        yerr2.append(sqrt(n[1])/float(n.sum()))
        yerr3.append(sqrt(n[2])/float(n.sum()))        
        print n,bins,v
    data1.append(48.0/100.0)
    data2.append(24.0/100.0)
    data3.append(28.0/100.0)
    yerr1.append(0.3/100.0)
    yerr2.append(1.6/100.0)
    yerr3.append(1.9/100.0)    
#    obs,n,v =sampleObservation(3000)   
#    data3.append(v[0])
#    data2.append(v[1])
#    data1.append(v[2])
#    yerr1.append(sqrt(n[0])/float(n.sum()))
#    yerr2.append(sqrt(n[1])/float(n.sum()))
#    yerr3.append(sqrt(n[2])/float(n.sum()))            
    plt.close("all")
    fig, ax = plt.subplots()
    fig.set_size_inches(12, 10)
    width=0.35
    loci1 = ax.bar(ind, data1, width, color='grey', yerr=yerr1)
    loci2 = ax.bar(ind, data2, width, color='white', bottom=numpy.array(data1), yerr=yerr2)
    loci3 = ax.bar(ind, data3, width, color='black', yerr=yerr3,bottom=numpy.array(data1)+numpy.array(data2))
    ax.set_ylabel('% loci ENV')
#    ax.title('Scores by group and gender')
    ax.set_xticks(ind+width/2. )
    ax.set_xticklabels(xlabels)    
    #plt.yticks(np.arange(0,1,10))
    ax.legend( (loci1[0], loci2[0], loci3[0]), ('3 foci', '2 focis', '1 focis') )
#    autolabelyerr(ax,loci1,err=yerr1)
#    autolabelyerr(ax,loci2,err=yerr2)
#    autolabelyerr(ax,loci1,err=yerr1)
#    autolabelyerr(ax,loci2,err=None)
#    autolabelyerr(ax,loci3,err=None)
#    autolabelyerr(ax,loci3,err=None)
    autolabels(loci1,loci2,loci3,ax,yerr1,yerr2,yerr3)    
    plt.show()
    return  data1,data2,data3,yerr1,yerr2,yerr3

def plotFociEnv(foci,nbENV):
    env=[]
    env.extend(nbENV)
    env.extend(nbENV)
    env.extend(nbENV)
    datas=np.zeros((15-5,3))
    yerrors=np.zeros((15-5,3))
    #foreach env number (5,15) count foci nb
    for i,f in enumerate(foci) :
        datas[env[i]-5][f-1]+=1
    width=0.35
    ind=np.arange(5,15)
    s=np.sum(datas,axis=1)
    plt.close("all")
    data1,data2,data3=datas.transpose()/s
    yerr1,yerr2,yerr3=np.sqrt(datas.transpose())/s
    fig, ax = plt.subplots()
    loci1 = ax.bar(ind, data3, width, color='grey', )#yerr=yerr3)
    loci2 = ax.bar(ind, data2, width, color='white', bottom=numpy.array(data3))#, yerr=yerr2)
    loci3 = ax.bar(ind, data1, width, color='black', bottom=numpy.array(data3)+numpy.array(data2))#yerr=yerr1,
    ax.set_ylabel('% loci ENV')
#    ax.title('Scores by group and gender')
    ax.set_xticks(ind+width/2. )
    ax.set_xticklabels(ind)    
    #plt.yticks(np.arange(0,1,10))
#    autolabels(loci1,loci2,loci3,ax,yerr1,yerr2,yerr3)
#    ax.legend( (loci1[0], loci2[0], loci3[0]), ('3 foci', '2 focis', '1 focis') )
    plt.xticks(np.arange(5,15))
    
def getRndWeighted(listPts,weight,yerr):
        w=[yerr[i]*np.random.random()+weight[i] for i in range(len(weight)) ]
        t = numpy.cumsum(w)
        s = numpy.sum(w)
        i = numpy.searchsorted(t,numpy.random.rand(1)*s)[0]
        return listPts[i]
        
def getForwWeight(listPts,weight):
        dice = random()
        #sorted ?
        for i in range(len(listPts)) :
            if weight[i] > dice and weight[i] != 0 :
                return listPts[i]
                
def getSubWeight(data,weights):
#    weights = numpy.take(self.weight,listPts)
    rnd = random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return data[i]
            
def sampleObservation(N):
    values=[1,2,3]
    weight = [0.48,0.24,0.28]
    yerr=[0.3/100.0,1.6/100.0,1.9/100.0]
#    data = [getSubWeight(values,weight,yerr) for i in range(100)]
    data = [getRndWeighted(values,weight,yerr) for i in range(N)]
    n, bins = np.histogram(data, 3, normed=0)
    v=n/float(n.sum())
    return data,n,v         
    
    
def countLocisAx(ax,x,nb=4,xlabels=('RAND-noMA', 'RAND-offMA', 'RAND-inMA', 'closeMA', 'closeMA-ENV','close-ENV','PR-')):
    #stacked bar
    #x = np.loadtxt("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/loci_output/loci_res.txt")
    data1=[]
    data2=[]
    data3=[]
    yerr1=[]
    yerr2=[]
    yerr3=[]
    nn=[]
    ind=np.arange(nb+1)
    for k in range(nb) :
        n, bins = np.histogram(x[k], 3, normed=0)
        v=n/float(n.sum())
        data3.append(v[0])
        data2.append(v[1])
        data1.append(v[2])
        yerr1.append(sqrt(n[0])/float(n.sum()))
        yerr2.append(sqrt(n[1])/float(n.sum()))
        yerr3.append(sqrt(n[2])/float(n.sum()))        
        print n,bins,v
    data1.append(48.0/100.0)
    data2.append(24.0/100.0)
    data3.append(28.0/100.0)
    yerr1.append(0.3/100.0)
    yerr2.append(1.6/100.0)
    yerr3.append(1.9/100.0)

    #plt.close("all")
    #fig, ax = plt.subplots()
    width=0.35
    loci1 = ax.bar(ind, data1, width, color='grey', yerr=yerr1)
    loci2 = ax.bar(ind, data2, width, color='white', bottom=numpy.array(data1), yerr=yerr2)
    loci3 = ax.bar(ind, data3, width, color='black', yerr=yerr3,bottom=numpy.array(data1)+numpy.array(data2))
    ax.set_ylabel('% loci ENV')
#    ax.title('Scores by group and gender')
    ax.set_xticks(ind+width/2. )
    ax.set_xticklabels(xlabels)    
    #plt.yticks(np.arange(0,1,10))
    ax.legend( (loci1[0], loci2[0], loci3[0]), ('3 foci', '2 focis', '1 focis') )
    autolabels(loci1,loci2,loci3,ax,yerr1,yerr2,yerr3)
#    autolabel(loci1,ax)
#    autolabel(loci2,ax)
#    autolabel(loci3,ax)
#    autolabels(loci1,loci2,loci3,ax)    
#    plt.show()

def doOne(filename,rname,iname):
    fluos = fluoSim(boundingBox=4195,res=400.0,space=100.0) 
    data = fluos.parseAPresults(filename,recipe_name = rname,ingr_name=iname)
    fluos.setValues(values=data)
    gauss_conv=fluos.getConvolution()
    g=fluos.setAtScale(200.)
    fluos.data_conv = g
    locis,npx,limits = fluos.getCluster1DXYZ()
    fluos.plotRotXYZ()
    return data,gauss_conv

def drawAsColoredInstance3D(helper,data):#or instance  
    shifted_cmap = cm.get_cmap('hot', 30)#shiftedColorMap(cm.get_cmap('hot', 30), start=0, midpoint=0.25/10., stop=0.55/10.0,  name='shifted') 
    colors = shifted_cmap(data/data.max()) 
    b=helper.box("base",size=[200,200,200])[0]    
    #X, Y , Z= np.ogrid[nx/2:nx/2, -ny/2:ny/2, -nz/2:nz/2]
    count=0    
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            for k in range(data.shape[2]):
                a=helper.newInstance("map"+str(i)+"_"+str(j)+"_"+str(k),b,location=[i*200.0,j*200.0,k*200.0])
                helper.changeObjColorMat(a,colors[i,j,k][:3])
                p=(count/float(data.shape[0]*data.shape[1]*data.shape[2]))*100.0
                helper.progressBar(progress=p,label=str(p))

def drawAsColoredInstance2D(helper,data):#or instance  
    shifted_cmap = cm.get_cmap('hot', 30)#shiftedColorMap(cm.get_cmap('hot', 30), start=0, midpoint=0.25/10., stop=0.55/10.0,  name='shifted') 
    colors = shifted_cmap(data/data.max()) 
    b=helper.box("base",size=[200,200,200])[0]    
    #X, Y , Z= np.ogrid[nx/2:nx/2, -ny/2:ny/2, -nz/2:nz/2]
    count=0    
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            a=helper.newInstance("map"+str(i)+"_"+str(j),b,location=[i*200.0,j*200.0,0.0])
            helper.changeObjColorMat(a,colors[i,j][:3])
            p=(count/float(data.shape[0]*data.shape[1]))*100.0
            helper.progressBar(progress=p,label=str(p))
#            print shifted_cmap(data[i,j])[:3]
#                    parent = None,material=None,**kw)

def makeImage(helper,name,data,w,imfilename):
    shifted_cmap = cm.get_cmap('hot', 30)#shiftedColorMap(cm.get_cmap('hot', 4195), start=0, midpoint=0.5, stop=0.55,  name='shifted')#cm.get_cmap('hot', 30)#shiftedColorMap(cm.get_cmap('hot', 30), start=0, midpoint=0.25/10., stop=0.55/10.0,  name='shifted') 
    #apply colorsmap to data
    colors = shifted_cmap(data/data.max()) 
    print colors.shape
    import Image
    im = Image.fromarray(np.uint8(colors*255),mode = None)
    print im.size
    im.save(imfilename)
    return im
    
def drawAsColoredPlane(helper,name,data,w,imfilename):#or instance  
    p = helper.getObject(name)
    if p is None :
        p,mpl = helper.plane(name,center = [0,0,0],size=[w,w],parent=None)    
    im=makeImage(helper,name,data,w,imfilename)
    mat = helper.createTexturedMaterial(name+"_Mat",imfilename)
    #assign the material to the plane
    helper.assignMaterial(p,mat,texture=True)
    return p

def updateP(helper,p,data,w,imfilename):
    p=helper.getObject(p)
    im=makeImage(helper,name,data,w,imfilename)
    mat = helper.createTexturedMaterial(name+"_Mat",imfilename)
    helper.assignMaterial(p,mat,texture=True)
    
#def finalResult(count,sample)::
#    plt.close('all')
#    extE=[-1200,1200,-1200,1200]
#    f, axes = plt.subplots(4,2)
#    #plot    
    
i=0
j=0
k=0
locis_g={}
locis_c={}
npx_g={}
limits_g={}
nbENV={}
rnbENV={}
gauss_conv_sample={}
for j in range(9):
    locis_g[j]=[]
    locis_c[j]=[]
    gauss_conv_sample[j]=None
    npx_g[j]=[]
    limits_g[j]=[]
    nbENV[j]=[]
    rnbENV[j]=[]
doplot1 = False
doplot = False
a=None
shifted_cmap = shiftedColorMap(cm.get_cmap('hot', 30), start=0, midpoint=0.25, stop=0.55,  name='shifted')
#plt.close('all')
extE=[-W,W,-W,W]

#execfile("/Users/ludo/DEV/func.py")
N=3000
d=np.random.random_integers(1, 3, N)
locis_g[0]=d              

def fig2data ( fig , ax,nx=41,ny=41):
    """
    @brief Convert a Matplotlib figure to a 4D numpy array with RGBA channels and return it
    @param fig a matplotlib figure
    @return a numpy 3D array of RGBA values
    from http://www.icare.univ-lille1.fr/wiki/index.php/How_to_convert_a_matplotlib_figure_to_a_numpy_array_or_a_PIL_image
    """
    matplotlib.use('tkAgg')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    fig.canvas.setMaximumHeight(ny)
    fig.canvas.setMaximumWidth(nx)
    # draw the renderer
    fig.canvas.draw ( ) 
    # Get the RGBA buffer from the figure
    w,h = fig.canvas.get_width_height()
    buf = numpy.fromstring ( fig.canvas.tostring_argb(), dtype=numpy.uint8 )
    buf.shape = ( w, h,4 ) 
    # canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
    buf = numpy.roll ( buf, 3, axis = 2 )
#    matplotlib.use('osx')
    return buf
    
def drawEmittor(data,bbox):
    matplotlib.use('tkAgg')
    from matplotlib.patches import Circle,Rectangle
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for p in data :
#        ax.add_patch(Circle((p[0], p[1]), 20.0,
#                        edgec,olor="red", facecolor="red"))
        ax.add_patch(Rectangle((p[0], p[1]),100.0, 100.0,
                        edgecolor="red", facecolor="red"))
    ax.set_aspect(1.0)
    plt.axhline(y=bbox[0][1], color='k')
    plt.axhline(y=bbox[1][1], color='k')
    plt.axvline(x=bbox[0][0], color='k')
    plt.axvline(x=bbox[1][0], color='k')
    plt.axis([bbox[0][0], bbox[1][0],
                 bbox[0][1], bbox[1][1]])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)
    ax.spines['bottom'].set_linewidth(0)
    ax.spines['left'].set_linewidth(0)
    ax.spines['bottom'].set_visible(False)
    return fig,ax             

def crop(image, x1, x2, y1, y2):
    """
    Return the cropped image at the x1, x2, y1, y2 coordinates
    """
    if x2 == -1:
        x2=image.shape[1]-1
    if y2 == -1:
        y2=image.shape[0]-1

    mask = np.zeros(image.shape)
    mask[y1:y2+1, x1:x2+1]=1
    m = mask>0
    return image[m].reshape((y2+1-y1, x2+1-x1))

def FigureM6():
    labels = ('RAND-noMA', 'RAND-offMA', 'RAND-inMA', 'closeMA', 'closeMA-ENV','close-ENV')
    locis_g,npx_g,limits_g,nbENV = analyseAll(threshold=None, limit=20)
    for i in range(1,7):
        plotFociEnv(locis_g[i],nbENV[i])
        plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/results_foci_env_%s.ps" % labels[i-1],dpi=300)
    return locis_g,npx_g,limits_g,nbENV
    
def FigureM7a(name,i,j,k):
    fluos,a,gauss_conv,env_results=analyse_one(i,j,k,0,space=100)
    res=[]
    binspx=[10,15,20,25,30,35,40,45,50,60,70,80]
    for i in range(3):
#        f, axes = plt.subplots(6,1)
        fluos.plotData(values=fluos.data,ax=i,axe=None,hide=True)
        plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/figure7/"+name+"_data_"+str(i)+".svg",dpi=300)
        fluos.plotData(values=fluos.data_conv,ax=i,axe=None,hide=True)
        plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/figure7/"+name+"_gauss_"+str(i)+".svg",dpi=300)
        #plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/figure7/input_2_2_10_"+str(i)+".ps",dpi=300)
        f, axes = plt.subplots(4,2)    
#        for j,limit in enumerate([0,0.005,0.010,0.015,0.020,0.025,0.030,0.035]):
        for j,limit in enumerate([0,0.010,0.020,0.030]):
            if limit == 0 :
                limit=None
            fo=[]
            for k,l in enumerate(binspx):
                loci,npx,limits=fluos.getCluster1D(limit=limit,uplimit=l,ax=i)#axes[k][j]
#                plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/figure7/feature_2_2_10_%d_%s_%d.ps" % (i,str(limit),l),dpi=100)
                res.append(loci)
                fo.append(loci)
            print f
            loci,npx,limits=fluos.getCluster1D(limit=limit,uplimit=l,axe=axes[j][0],ax=i,hide=True)#axes[k][j] 
            axes[j][1].bar(binspx,fo)
            axes[j][1].set_xticks(binspx )
            axes[j][1].set_xticklabels(binspx)    
            axes[j][1].set_yticks([1,2,3] )
            axes[j][1].set_yticklabels([1,2,3])    
#            axes[j][1].xaxis.set_visible(False)
#            axes[j][1].yaxis.set_visible(False)
            #the othe imae should be blanck with the 4 focis values 
        plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/figure7/%s_feature_%d.ps" % (name,i),dpi=300)
    return res

def formatResultM7a(res):
    c=0
    for i in range(3):
        for limit in [0,0.010,0.020,0.030]:
            for l in [10,15,20,25]:
                print "Axe "+str(i)+" threshold "+str(limit)+" uplimit "+str(l)+" focis "+str(res[c])
                c+=1

                
def analyse_one(i,j,k,ax,space=100.0):
    fluos = fluoSim(boundingBox=4195,res=400.0,space=space)
    filename="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/lorentz_%d_%d_%d.ccp4" % (i,j,k)
    a=numpy.loadtxt(filename+".txt",delimiter=",")
    fluos.setValuesAccurate(values=a)
#    fluos.setDataAtScale(200.0)
    gauss_conv=fluos.getConvolution()
    g=fluos.setAtScale(200.)
    fluos.data_conv = g#    fluos.plotRotXYZ()
#    data1=np.sum(gauss_conv,axis=ax)   
#    drawAsColoredInstance2D(helper,data1) 
    with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/env_models/model_%d_%d_%d.json"%(i,j,k), 'r') as fp :#doesnt work with symbol link ?
        env_results=json.load(fp)
    
    return fluos,a,gauss_conv,env_results

def getOne(i,j,k,limit,threshold,space=100,redo=False):
    fluos = fluoSim(boundingBox=4195,res=400.0,space=space)
    filename="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/lorentz_%d_%d_%d.ccp4" % (i,j,k)
    a=numpy.loadtxt(filename+".txt",delimiter=",")
    print a
    fluos.setValuesAccurate(values=a)#space 100 for value ie 80pixels
    if redo :
        gauss_conv=fluos.getConvolution()
        g=fluos.setAtScale(200.)
        fluos.data_conv = g#    fluos.plotRotXYZ()
    else :        
        f="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output2/gauss_conv_%d_%d_%d.npy" % (i,j,k)
        fluos.data_conv=numpy.load(f) #this array come from computato earier
        #    g=fluos.setAtScale(200.0)
        fluos.setDataAtScale(200.0)
#    fluos.data_conv=g
    print fluos.data_conv.shape
    fim=None#"/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output/fluo_gauss_conv_%d_%d_%d.jpg" % (i,j,k)
    #0.02
    locis,npx,limits = fluos.getCluster1DXYZ(limit=threshold,uplimit=limit,saveimgeFile=fim)
    return a,fluos,locis,npx,limits

def displayOne(helper,i,j,k,limit,threshold,space=100,redo=False):
    a,fluos,locis,npx,limits = getOne(i,j,k,limit,None,space=space,redo=redo)
    ax=0
    print fluos.W,fluos.data_conv.shape,fluos.data.shape
    drawAsColoredPlane(helper,"gaussX",np.sum(fluos.data_conv,axis=ax),fluos.W*2.,"/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output_%d_%d_%d_%d.jpg" % (i,j,k,ax))
    drawAsColoredPlane(helper,"posX",np.sum(fluos.data,axis=ax),fluos.W*2.,"/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output_p_%d_%d_%d_%d.jpg" % (i,j,k,ax))
    ax=1
    drawAsColoredPlane(helper,"gaussY",np.sum(fluos.data_conv,axis=ax),fluos.W*2.,"/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output_%d_%d_%d_%d.jpg" % (i,j,k,ax))
    drawAsColoredPlane(helper,"posY",np.sum(fluos.data,axis=ax),fluos.W*2.,"/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output_p_%d_%d_%d_%d.jpg" % (i,j,k,ax))
    ax=2
    drawAsColoredPlane(helper,"gaussZ",np.sum(fluos.data_conv,axis=ax),fluos.W*2.,"/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output_%d_%d_%d_%d.jpg" % (i,j,k,ax))
    drawAsColoredPlane(helper,"posZ",np.sum(fluos.data,axis=ax),fluos.W*2.,"/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output_p_%d_%d_%d_%d.jpg" % (i,j,k,ax))
    return a,fluos,locis,npx,limits



def analyseAll(threshold=None, limit=15):
    i=0
    j=0
    k=0
    locis_g={}
    locis_c={}
    npx_g={}
    limits_g={}
    nbENV={}
    rnbENV={}
    gauss_conv_sample={}
    for j in range(9):
        locis_g[j]=[]
        locis_c[j]=[]
        gauss_conv_sample[j]=None
        npx_g[j]=[]
        limits_g[j]=[]
        nbENV[j]=[]
        rnbENV[j]=[]
    fluos = fluoSim(boundingBox=4195,res=400.0,space=100.0)
    n=0
    for i in range(10):#MA distrib
        for j in range(6):#ENV hypothesis
            for k in range(100):#ENV distrib    
#                with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/env_models/model_%d_%d_%d.json"%(i,j,k), 'r') as fp :#doesnt work with symbol link ?
#                    env_results=json.load(fp)
#                rnbENV[j+1].append(env_results[0][0])            
                filename="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/lorentz_%d_%d_%d.ccp4" % (i,j,k)
                a=numpy.loadtxt(filename+".txt",delimiter=",")
                nbENV[j+1].append(len(a))
                f="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output2/gauss_conv_%d_%d_%d.npy" % (i,j,k)
                fluos.data_conv=numpy.load(f)
                fim=None#"/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output/fluo_gauss_conv_%d_%d_%d" % (i,j,k)
                #0.02
                locis,npx,limits = fluos.getCluster1DXYZ(limit=threshold,uplimit=limit,saveimgeFile=fim)
                locis_g[j+1].extend(locis)
                npx_g[j+1].extend(npx) 
                limits_g[j+1].extend(limits)   
    return locis_g,npx_g,limits_g,nbENV

def FigureM7b(results):
    labels = ('R','RAND-noMA', 'RAND-offMA', 'RAND-inMA', 
                        'closeMA', 'closeMA-ENV','close-ENV','PR-')
    for i,thr in enumerate([0.010,0.020,0.030]):
        for j,l in enumerate([10,15,20,25]):
            locis_g = results[0][i][j]
            countLocis(locis_g,nb=len(labels)-1,xlabels=labels)
            plt.title("Focis Count at threshold = %.2f and uplimit = %d" % (thr,l))
            plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/results_final%d_%d.svg" % (i,j),dpi=300)
#    
    
def sampleSize(Z,std,margin):
    #    Necessary Sample Size = (Z-score)² – StdDev*(1-StdDev) / (margin of error)²
    #((1.96)² x .5(.5)) / (.05)²
    return (Z**2*std**2)/margin**2

    
def runAllThr():
    dLOCIS=[]
    dNPX=[]
    dlimits=[]
    dNENV=[]
    for thr in [0.010,0.020,0.030]:
        LOCIS=[]
        NPX=[]
        limits=[]
        NENV=[]
        print "test thr ",thr
        for l in [10,15,20,25]:
            print "test limit l",l
            locis_g,npx_g,limits_g,nbENV = analyseAll(threshold=thr, limit=l)
            LOCIS.append(locis_g)
            NPX.append(npx_g)
            limits.append(limits_g)
        dLOCIS.append(LOCIS)
        dNPX.append(NPX)
        dlimits.append(limits)
    results=[dLOCIS,dNPX,dlimits,dNENV]
    with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/results_final.json", 'w') as fp : 
        json.dump(results,fp)  
    return results
    
#helper = upy.getHelperClass()()
#displayOne(helper,0,0,0,15.0,None,space=100)
#fluos,a,gauss_conv,env_results=analyse_one(2,2,10,0,space=100);
#fluos.plotRotXYZ(uplimit=15)

#uplimit=25)
#plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/plotXYZ_000.png",dpi=300)
#fig,ax=drawEmittor(a[:,:2],fluos.boundingBox)    
#buf=fig2data ( fig , ax,nx=fluos.nx*2,ny=fluos.ny*2)
#plt.close('all')
##plt.imshow(buf)
#RGB=np.sum(buf,axis=2)
#labeled, nr_objects = ndimage.label(buf[:,:,1],structure=np.ones((3,3), dtype="bool8"))
#labeled[labeled<=1]=0
#labeled[labeled>1]=1
#plt.imshow(labeled)
#convolve labeled
#data_conv = np.real(fluos.convolve(labeled,fluos.psf))
#
anal = False
#clean buf from balck axis line, and keep only the red
if anal == True:
    fluos = fluoSim(boundingBox=4195,res=400.0,space=100.0)
    n=0
    for i in range(10):#MA distrib
        for j in range(6):#ENV hypothesis
            for k in range(100):#ENV distrib    
    #            print i,j,k,int(i*100*6 + j*100 + k),n
    #            n+=1
                with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/env_models/model_%d_%d_%d.json"%(i,j,k), 'r') as fp :#doesnt work with symbol link ?
                    env_results=json.load(fp)
    #            nbENV[j+1].append(env_results[0][0])
                #apply a random rotation instead ? and compare ?
                rnbENV[j+1].append(env_results[0][0])            
                filename="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/lorentz_%d_%d_%d.ccp4" % (i,j,k)
                a=numpy.loadtxt(filename+".txt",delimiter=",")
                nbENV[j+1].append(len(a))
                fluos.setValuesAccurate(values=a)
    #            fluos.setDataAtScale(200.0)
                gauss_conv=fluos.getConvolution()
                g=fluos.setAtScale(200.)
                fluos.data_conv = g
                locis,npx,limits = fluos.getCluster1DXYZ(uplimit=20)
                locis_g[j+1].extend(locis)
                npx_g[j+1].extend(npx)
                limits_g[j].extend(limits)
                fluos.writeDATA("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output2/gauss_conv_%d_%d_%d" % (i,j,k))
    #            f="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/gauss_conv_%d_%d_%d.npy" % (i,j,k)
    #            fluos.data_conv=numpy.load(f)
    #            fim="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output/fluo_gauss_conv_%d_%d_%d.jpg" % (i,j,k)
    #            0.02
    #            locis,npx,limits = fluos.getCluster1DXYZ(limit=None,uplimit=15,saveimgeFile=fim)
    #            locis_g[j+1].extend(locis)
    #            npx_g[j+1].extend(npx) 
    #            limits_g[j+1].extend(limits)
LOCIS=[]
NPX=[]
limits=[]
NENV=[]
labels = ('RAND-noMA', 'RAND-offMA', 'RAND-inMA', 
                        'closeMA', 'closeMA-ENV','close-ENV','PR-')

#for l in [10,15,20,25]:
#    print "test limit l",l
#    locis_g,npx_g,limits_g,nbENV = analyseAll(threshold=None, limit=l)
#    LOCIS.append(locis_g)
#    NPX.append(npx_g)
#    limits.append(limits_g)
#    plt.close('all')
#    fig, ax = plt.subplots()
#    fig.set_size_inches(12, 10)
#    plt.title("Focis Count at threshold = auto and uplimit = %d" % l)
#    countLocisAx(ax,[locis_g[1],locis_g[2],locis_g[3],locis_g[4],locis_g[5],locis_g[6]],nb=len(labels)-1,xlabels=labels)
#    plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/results_final%d_%s_%d_%f.png" % (7,"auto",j,l),dpi=300)
#    print "done"
##locis_g,npx_g,limits_g,nbENV = analyseAll(threshold=None, limit=20)
#for j,l in enumerate([10,15,20,25]):
#    np.savetxt("foci_counts_"+str(l)+".csv",[LOCIS[i][1],LOCIS[i][2],LOCIS[i][3],LOCIS[i][4],LOCIS[i][5],LOCIS[i][6]])
#for i in range(1,7):
#    plotFociEnv(locis_g[i],nbENV[i])
#    plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/results_foci_env_%d.ps" % i,dpi=300)
#
#data,n,v=sampleObservation(N)
#locis_g = np.loadtxt("foci_counts_20.csv")
#Z=1.96
#Z*std/np.sqrt(N)
#You can also use this handy formula in finding the confidence interval: x̅ ± Za/2 * σ/√(n). Here, x̅ represents the mean.

#LOCIS=[]
#for j,l in enumerate([10,15,20,25]):
#    l=np.loadtxt("foci_counts_"+str(l)+".csv",[LOCIS[i][1],LOCIS[i][2],LOCIS[i][3],LOCIS[i][4],LOCIS[i][5],LOCIS[i][6]])
#    LOCIS.append()
#save csv
#from scipy import stats
#chisq, pval = stats.chisquare(r, np.array(data))
#likelihood = ('possible' if pval > 0.05 else
#                  'unlikely' if pval > 0.01 else 'HIGHLY UNLIKELY')
##make the figures
#labels = ('R','RAND-noMA', 'RAND-offMA', 'RAND-inMA', 
#                        'closeMA', 'closeMA-ENV','close-ENV','PR-')
#fig, ax = plt.subplots(7,5)
#for i,thr in enumerate([0.005,0.010,0.015,0.020,0.025,0.030,0.035]):
#    #histgram of NPX
#    n, bins, patches = ax[i][0].hist( [results[1][i][0][1],results[1][i][0][2],results[1][i][0][3],
#                                 results[1][i][0][4],results[1][i][0][5],results[1][i][0][6]], 10, 
#                                histtype='bar',label=('RAND-noMA', 'RAND-offMA', 'RAND-inMA', 'closeMA', 
#                                'closeMA-ENV','closeENV'));
#    #ax[i][0].legend()
#    for j,l in enumerate([10,15,20,25]):
#        countLocisAx(ax[i][j+1],results[0][i][j],nb=len(labels)-1,xlabels=labels)

##fig, ax = plt.subplots(7,5)
#for i,thr in enumerate([0.005,0.010,0.015,0.020,0.025,0.030,0.035]):
#    #histgram of NPX
#    plt.close('all')
#    fig, ax = plt.subplots()
#    fig.set_size_inches(12, 10)
#    plt.title("Area Coverage at threshold = %.3f" % thr)
#    n, bins, patches = ax.hist( [results[1][i][0][1],results[1][i][0][2],results[1][i][0][3],
#                                 results[1][i][0][4],results[1][i][0][5],results[1][i][0][6]], 10, 
#                                histtype='bar',label=('RAND-noMA', 'RAND-offMA', 'RAND-inMA', 'closeMA', 
#                                'closeMA-ENV','closeENV'));
#    ax.legend()
#    fig.set_size_inches(12, 10);plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/results_final%d_%f.png" % (i,thr),dpi=300)
#    for j,l in enumerate([10,15,20,25]):
#        plt.close('all')
#        fig, ax = plt.subplots()
#        fig.set_size_inches(12, 10)
#        plt.title("Focis Count at threshold = %.3f and uplimit = %d" % (thr,l))
#        countLocisAx(ax,results[0][i][j],nb=len(labels)-1,xlabels=labels)
#        plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/results_final%d_%f_%d_%f.png" % (i,thr,j,l),dpi=300)
# 
#locis_g,npx_g,limits_g,nbENV = analyseAll(threshold=None, limit=15.0)
#LOCIS=[]
#NPX=[]
#limits=[]
#NENV=[]
##print "test thr ",thr
#for l in [10,15,20,25]:
#    print "test limit l",l
#    locis_g,npx_g,limits_g,nbENV = analyseAll(threshold=None, limit=l)
#    LOCIS.append(locis_g)
#    NPX.append(npx_g)
#    limits.append(limits_g)
#
#dLOCIS.append(LOCIS)
#dNPX.append(NPX)
#dlimits.append(limits)
#results=[dLOCIS,dNPX,dlimits,dNENV]
#with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/results_final_7.json", 'w') as fp : 
#    json.dump(results,fp)  


#plt.close('all')
#fig, ax = plt.subplots()
#fig.set_size_inches(12, 10)
#plt.title("Area Coverage at threshold = auto")
#n, bins, patches = ax.hist( [npx_g[1],npx_g[2],npx_g[3],
#                             npx_g[4],npx_g[5],npx_g[6]], 10, 
#                            histtype='bar',label=('RAND-noMA', 'RAND-offMA', 'RAND-inMA', 'closeMA', 
#                            'closeMA-ENV','closeENV'));
#ax.legend()
#plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/results_final%d_%s.png" % (7,"auto"),dpi=300)
#
##correlatio nbEnv <-> loci
#for i in range(1,7):
#    plotFociEnv(locis_g[i],nbENV[i])
#    plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/results_foci_env_%d.png" % i,dpi=300)
##locis_g,npx_g,limits_g,nbENV = analyseAll(threshold=0.05, limit=15)
##  
#import upy
#helper = upy.getHelperClass()()
#
#for i in range(1,8):
#    Y,X=np.histogram(npx_g[i]);
#    plt.bar(X[:-1]+(i*0.1),Y,width=0.5,color=np.random.rand(3,1))
#filename1="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/farady/HIV_1_7result1.json"
#data1,gauss_conv=doOne(filename1,'HIV1_envelope_Pack_106_0_2_0_surfaceRecipe','HIV1_envelope_Pack_106_0_2_0_surf__HIV1_ENV_4nco_0_1_1')
#filename2="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/farady/HIV_1_8result1.json"
#data2,gauss_conv2=doOne(filename2,'HIV1_envelope_Pack_183_0_2_0_surfaceRecipe','HIV1_envelope_Pack_183_0_2_0_surf__HIV1_ENV_4nco_0_1_1')
###sp=helper.Sphere("sph")[0]
###isp=helper.instancesSphere("sted",data,[200,]*len(data),sp,[[1,0,0]],None)
#fluos = fluoSim(psf_type="gauss",boundingBox=4195,res=400.0,space=200.0) 
#fluos.setValues(values=data1);
#gauss_conv=fluos.getConvolution()
##locis = fluos.getCluster1DXYZ()
#
##fluos.plotRotXYZ()
#points=data1
#values=np.ones(data1.shape[0])
#X, Y , Z= np.ogrid[-fluos.nx/2:fluos.nx/2, -fluos.ny/2:fluos.ny/2, -fluos.nz/2:fluos.nz/2]
##X, Y , Z= np.ogrid[-fluos.W:fluos.W, -fluos.W:fluos.W, -fluos.W:fluos.W]
#grid_z0 = scipy.interpolate.griddata(points, values, (X*fluos.space, Y*fluos.space, Z*fluos.space), method='nearest')
#
##data1=np.sum(gauss_conv,axis=1)   
#drawAsColoredInstance2D(helper,data1)    
#for i in range(10):#MA distrib
#    for j in range(8):#ENV hypothesis
#        for k in range(100):#ENV distrib    
#            print i,j,k        
#            filename="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/lorentz_%d_%d_%d.ccp4" % (i,j,k)
#            values,psf_gauss,psf_cauchy,gauss_conv,cauchy_conv=convoluteGauss(filename)
#            if doplot1 :
#                a=plotRot(values,gauss_conv,limit=0.4,nbup=30)
#                plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output/lorentz_%d_%d_%d.jpg" % (i,j,k))
#            else :
#                for ax in range(3):
#                    if doplot :
##                        a=plotRot(gauss_conv,limit=0.2,nbup=30)
##                        plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output/lorentz_%d_%d_%d_%d.jpg" % (i,j,k,ax))                    
#                        axes=plot (values,psf_gauss,psf_cauchy,gauss_conv,cauchy_conv,ax=ax)
#                    data=numpy.sum(np.real(gauss_conv),axis=ax)
#                    gauss_conv_sample[j]=data
#                    fmax=np.real(gauss_conv).max()+data.std()
#                    if doplot :
#                        a=axes[3,0]
#                    #0.35 - 20 - 1
#                    #0.2 - 27 - 1
#                    loci,npx=cluster (data,fmax,30,2,axe=a,outputfile=None)
##                    loci,npx=cluster2(data)
#                    locis_g[j+1].append(loci)
#                    npx_g[j].append(npx)
#    #                data=numpy.sum(np.real(cauchy_conv),axis=ax)
#    #                if doplot :
#    #                    a=axes[3,1]
#    #                labeled,loci=cluster (data,0.35,27,5,axe=a,outputfile=None)
#    #                locis_c[j].append(loci)
#                    if doplot :
#                        plt.savefig("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output/lorentz_%d_%d_%d_%d.jpg" % (i,j,k,ax))
#                        plt.show()
#countLoci(locis_g[0])
#countLoci(locis_c)
#countLocis(locis_g,nb=3)
#make a random on.
#n, bins, patches = plt.hist( [nbENV[1],nbENV[2],nbENV[3],nbENV[4],nbENV[5],nbENV[6]], 10, histtype='bar',label=('RAND-noMA', 'RAND-offMA', 'RAND-inMA', 'closeMA', 'closeMA-ENV','closeENV'));plt.legend()
#n, bins, patches = plt.hist( [rnbENV[1],rnbENV[2],rnbENV[4],rnbENV[5]], 10, histtype='bar',label=('RAND-noMA', 'RAND-offMA',  'closeMA', 'closeMA-ENV','closeENV'));plt.legend()
#
#n, bins, patches = plt.hist( [npx_g[1],npx_g[2],npx_g[3],npx_g[4],npx_g[5],npx_g[6]], 10, histtype='bar',label=('RAND-noMA', 'RAND-offMA', 'RAND-inMA', 'closeMA', 'closeMA-ENV','closeENV'));plt.legend()
#nbe=[nbENV[1],nbENV[2],nbENV[3],nbENV[4],nbENV[5]]#,nbENV[7]]
#ynbe=[]
#ynbes=[]
#data=[]
#for e in nbe:
#    y,binEdges=np.histogram(data,bins=10)
#    menStd = np.sqrt(y)
#    ynbe.append(np.array(e).std())
#    ynbes.append(menStd)
#    data.append(y)
#    
#rnbe=[rnbENV[1],rnbENV[2],rnbENV[3],rnbENV[4],rnbENV[5]]#,rnbENV[6]]#,rnbENV[7]]
#
labels = ('R','RAND-noMA', 'RAND-offMA', 'RAND-inMA', 
                        'closeMA', 'closeMA-ENV','close-ENV','PR-')
##                        ('R','RAND-noMA', 'RAND-offMA', 'RAND-inMA', 
##                        'closeMA', 'closeMA-ENV','close-ENV','MA-ENV30/30/40','MA-ENV60/20/20','PR-')
#n, bins, patches = plt.hist( nbe, 10, histtype='bar',label=('RAND-noMA', 'RAND-offMA', 'RAND-inMA', 'closeMA', 'closeMA-ENV','closeENV'));plt.legend()
#n, bins, patches = plt.hist( , 10, histtype='bar',label=('RAND-noMA', 'RAND-offMA',  'closeMA', 'closeMA-ENV','closeENV'));plt.legend()
#countLocis(locis_g,nb=len(labels)-1,xlabels=labels)
##with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/loci_output/locis_g.json", 'w') as fp : json.dump(locis_g,fp)        

#single: 28 +-1.9%
#double: 24 +-1.6%
#multiple: 48 +-0.3%
#(taken from Table S1 from the Chojnakci paper)
#0 39px->3
#1 31
#2 31 24+7
#3 1 14 
#4 12  1
#n=10
#k=3
#X = scipy.randn(n,2)
#d = sch.distance.pdist(X)
#Z= sch.linkage(d,method='complete')
#T = sch.fcluster(Z, k, 'maxclust')

#d=numpy.sum(np.real(gauss_conv),axis=1)
#d[d>d.mean()]=d.mean()
#orig_cmap = plt.cm.coolwarm
#orig_cmap = plt.cm.hot
#d=d/d.max()
#shifted_cmap = shiftedColorMap(orig_cmap, start=0, midpoint=0.5, stop=0.55,  name='shifted')
#shrunk_cmap = shiftedColorMap(orig_cmap, start=0.15, midpoint=0.75, stop=0.85, name='shrunk')
#plt.close('all')
#extE=[-1200,1200,-1200,1200]
##CS = plt.contourf(numpy.sum(np.real(gauss_conv),axis=1), 10, cmap=plt.cm.hot)
#plt.imshow(d,cmap=shifted_cmap,interpolation='none',extent=extE)
#plt.colorbar()
#plt.show()
#close('all')
#f, axes = plt.subplots()
#axes.imshow(numpy.sum(np.real(gauss_conv),axis=1),cmap=shifted_cmap,interpolation='none',extent=extE)
#bar = AnchoredSizeBar(axes.transData, 1000, '',loc=3,pad=2, sep=1, borderpad=0.01, frameon=False)#,prop=None)
#rect = bar.size_bar._children[0]
#rect.set_ec("white")
#rect.set_height(1)
#rect.set_fc('white')
#axes.add_artist(bar)
##__init__(self, transform, size, label, loc, pad=0.10000000000000001, borderpad=0.10000000000000001, sep=2, prop=None, frameon=True, **kwargs)
##import matplotlib.font_manager as fm
##fontprops = fm.FontProperties(size=14, family='monospace')
##bar = AnchoredSizeBar(axes.transData, 1000, '3 units', 4, pad=0.5, sep=5, borderpad=0.5, frameon=False, size_vertical=0.5, color='white', fontproperties=fontprops)
##bar.patch.set(alpha=0.5, boxstyle='round')
#axes.add_artist(bar)
#X=data
#bandwidth = estimate_bandwidth(data, quantile=0.2, n_samples=500)
#ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
#ms.fit(X)
#labels = ms.labels_
#cluster_centers = ms.cluster_centers_
#labels_unique = np.unique(labels)
#n_clusters_ = len(np.nonzero(labels_unique))
#
#print("number of estimated clusters : %d" % n_clusters_)
#from itertools import cycle
#
#pl.figure(1)
#pl.clf()
#
#colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
#for k, col in zip(range(n_clusters_), colors):
#    my_members = labels == k
#    cluster_center = cluster_centers[k]
#    pl.plot(X[my_members, 0], X[my_members, 1], col + '.')
#    pl.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
#            markeredgecolor='k', markersize=14)
#pl.title('Estimated number of clusters: %d' % n_clusters_)
#pl.show()
#markers = np.zeros(data.shape, dtype=np.uint)
#markers[data < 0.1] = 1
#markers[data > 1.3] = 2
#struct=np.ones((3,3), dtype="bool8")#starting point is the marker 
#scipy.ndimage.measurements.watershed_ift(input,marker,structure)
#m=1
#with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/ma_models/model_%d.json"%m, 'r') as fp :#doesnt work with symbol link ?
#            ma_results=json.load(fp)
#import scipy
#import pylab
#import scipy.ndimage.measurements as measurements
#data = scipy.ones((20,20), dtype=scipy.uint8)
#for i in range(20):
#    for j in range(20):
#        data[i,j] = (1+scipy.sin(i/5.)*scipy.sin(j/5.))*128
#
#startingpoints = scipy.zeros(scipy.shape(data), dtype=int)
#startingpoints[5,5]=1
#startingpoints[15,15]=2
#
#segmented = measurements.watershed_ift(data, startingpoints)
#
#pylab.pcolor(segmented)
#
#i=0;j=3;k=0
#filename="/Users/ludo/Downloads/autopack_ouput/hiv_experiment/gauss_output/lorentz_%d_%d_%d.ccp4" % (i,j,k)
#values,psf_gauss,psf_cauchy,gauss_conv,cauchy_conv=convoluteGauss(filename)
#data=numpy.sum(np.real(gauss_conv),axis=2)
#d=numpy.sum(np.real(psf_gauss),axis=2)
##imshow(data)
#fmax=np.real(gauss_conv).max()+data.std()
#data[data<fmax]=0
#labeled, nr_objects = ndimage.label(data,structure=np.ones((3,3), dtype="bool8"))
##imshow(labeled)
#imgX = misc.imread("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output/sample_%d_%d_%d_X.png"% (i,j,k),flatten=0);#Back rotate -90
#imgY = misc.imread("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output/sample_%d_%d_%d_Y.png"% (i,j,k),flatten=0);#Bottom rotate 180
#imgZ = misc.imread("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/images_output/sample_%d_%d_%d_Z.png"% (i,j,k),flatten=0);#Right rotated +90
#plotRot(values,gauss_conv,imgs=[imgX,imgY,imgZ])
###m=cm.get_cmap('hot', 30)
##cdata = numpy.array(m(data)*255,numpy.uint8)
##startingpoints=numpy.sum(np.real(values),axis=0)
##startingpoints= numpy.array(startingpoints,dtype=int)
##segmented = measurements.watershed_ift(cdata[:,:,0], startingpoints)
##from skimage import morphology
#data=data/data.max()
#blobs = data > 0.4#data.mean()
#all_labels = morphology.label(blobs)
#blobs_labels = morphology.label(blobs, background=0)
#imshow(blobs_labels)
#struct=np.ones((3,3), dtype="bool8")
#cd=data[:]
#cd[cd < 0.4] = 0
#labeled, nr_objects = ndimage.label(data,structure=struct)#,structure=s)# > data.max()/2.0);
#m=measure.regionprops(labeled, properties=['Area', 'Perimeter','EquivDiameter','FilledArea','Eccentricity'])
#imshow(labeled)
#8000 Federal
#3137 SSTax
#2180 State
#execfile("/Users/ludo/DEV/testGaussCluster.py")