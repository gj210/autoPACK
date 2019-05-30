# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 09:54:06 2014

@author: Ludovic Autin.
"""
import os
import sys
import json
import numpy as np
import scipy
from scipy import ndimage
from scipy import signal
from scipy import spatial
from scipy import misc
from scipy import stats

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
#from sklearn.cluster import MeanShift, estimate_bandwidth
#from sklearn.datasets.samples_generator import make_blobs
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib.offsetbox import AnchoredOffsetbox
#from skimage.segmentation import random_walker
#from skimage import measure

from Volume.IO.volWriters import WriteCCP4
from Volume.Grid3D import Grid3DUC, Grid3DSI, Grid3DF
import autopack
from autopack import Analysis
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

class fluoSim:
    def __init__(self,values=None,psf_type="gauss",
                 boundingBox=4195,res=400.0,space=50.0 ):
        self.mapping=None
        self.values = values
        self.gridxyz=None
        self.overwrite_spacing = None
#        self.setup(values=values,psf_type=psf_type,boundingBox=boundingBox,
#                   res=res,space=space)

    def makeGrid(self,values=None,boundingBox=2000,res=400.0,space=200.0 ):
        self.W=float(boundingBox)
        self.boundingBox=np.array([[-self.W, -self.W, -self.W], [self.W, self.W, self.W]])
        sp=int((self.W*2.)/space)
        if self.overwrite_spacing is not None :
            sp=self.overwrite_spacing
        self.x = np.linspace(self.boundingBox[0][0], self.boundingBox[1][0] , sp)#*1.1547) gridspacing is already multiplied by 1.1547
        self.y = np.linspace(self.boundingBox[0][1], self.boundingBox[1][1] , sp)#*1.1547)
        self.z = np.linspace(self.boundingBox[0][2], self.boundingBox[1][2] , sp)#*1.1547)
#        self.x = np.arange(self.boundingBox[0][0], self.boundingBox[1][0] + space, space)#*1.1547) gridspacing is already multiplied by 1.1547
#        self.y = np.arange(self.boundingBox[0][1], self.boundingBox[1][1] + space, space)#*1.1547)
#        self.z = np.arange(self.boundingBox[0][2], self.boundingBox[1][2] + space, space)#*1.1547)
        self.nx = len(self.x) # sizes must be +1 or the right, top, and back edges don't get any points using this numpy.arange method
        self.ny = len(self.y)
        self.nz = len(self.z)
        self.space = space

    def makeDiscreteGrid(self,boundingBox,space):
        self.W=float(boundingBox)
        self.boundingBox=np.array([[-self.W, -self.W, -self.W], [self.W, self.W, self.W]])
        self.x = np.arange(self.boundingBox[0][0]- space, self.boundingBox[1][0]+ space, space)#*1.1547) gridspacing is already multiplied by 1.1547
        self.y = np.arange(self.boundingBox[0][1]- space, self.boundingBox[1][1] + space, space)#*1.1547)
        self.z = np.arange(self.boundingBox[0][2]- space, self.boundingBox[1][2] + space, space)#*1.1547)
        self.nx = len(self.x) # sizes must be +1 or the right, top, and back edges don't get any points using this numpy.arange method
        self.ny = len(self.y)
        self.nz = len(self.z)
        self.space = space

    def setup(self,values=None,psf_type="gauss",boundingBox=4195,res=400.0,space=50.0):
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
        self.shifted_cmap = shiftedColorMap(cm.get_cmap('hot', 100), start=0, midpoint=0.5, stop=0.55,  name='shifted')
        self.extE=[-self.W,self.W,-self.W,self.W]
#        self.shifted_cmap = shiftedColorMap(cm.get_cmap('hot', 4000), start=-1.0, midpoint=0.5, stop=0.65,  name='shifted')
#        shifted_cmap = cm.get_cmap('hot', 30)#shiftedColorMap(cm.get_cmap('hot', 30), start=0, midpoint=0.25/10., stop=0.55/10.0,  name='shifted')

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
            self.mapping=[ix,iy,iz]

    def setValuesAccurate(self,values=None):
        if values is None :
            values = self.values
        if values is not None :
            #filter values not inside the Box boundingBox
            xyz = np.meshgrid(self.x,self.y,self.z,copy=False)
            grid = np.array(xyz).T
            self.gridxyz= grid.reshape((self.nx*self.ny*self.nz,3))
            tr=spatial.cKDTree(self.gridxyz, leafsize=10)
            d,i=tr.query(values)
            zvalues=np.zeros(self.nx*self.ny*self.nz)
            zvalues[i]=1.0
            self.data = zvalues.reshape((self.nx,self.ny,self.nz))
            self.mapping=i
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

    def getCluster1D(self,data1d=None,limit=None,ecc=0.5,uplimit=30,axe=None,outputfile=None,ax=0):
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
        if axe is not None :
            axe.imshow(labeled);
            axe.text(2, 3, "loci nb is "+str(loci)+" "+str(nr_objects)+" features "+str(npx)+" px ",color='white')
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

    def plotData(self,values=None,ax=0,axe=None):
        if values is None :
            return
        if axe is None :
            plt.close('all')
            f, axe = plt.subplots()
        data=np.sum(values,axis=ax)
        ticks=np.linspace(-self.W,self.W,5)
        axe.imshow(data,cmap=self.shifted_cmap,interpolation='none',extent=self.extE)
        axe.set_xticks(ticks)
        axe.set_yticks(ticks)
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

    def parseAPresults(self,recipe_name,ingr_name,filename=None,results=None):
        if filename is not None:
            with open(filename, 'r') as fp :#doesnt work with symbol link ?
                ap_results=json.load(fp)
        elif results is not None :
            ap_results= results
        #ap_results['HIV1_envelope_Pack_106_0_2_0_surfaceRecipe']['HIV1_envelope_Pack_106_0_2_0_surf__HIV1_ENV_4nco_0_1_1']['results']
        ingr_res = ap_results[recipe_name][ingr_name]['results']
        ingr_pos = np.array(ingr_res).transpose()[0]
        return np.array(ingr_pos.tolist(),dtype='f')

    def simulateFrom(self,rname,iname,bbox,res,space,rspace,filename=None,results=None):
        self.setup(boundingBox=bbox,res=res,space=space)
        data = self.parseAPresults(filename=filename,results=results,
                                    recipe_name = rname,ingr_name=iname)
        self.setValuesAccurate(values=data)
        gauss_conv=self.getConvolution()
        g=self.setAtScale(rspace)
        self.data_conv = g
#        locis,npx,limits = self.getCluster1DXYZ()
#        self.plotRotXYZ()
        return data,gauss_conv


    def drawAsColoredInstance3D(self,helper,data,rspace):#or instance
        if helper is None :
            return
        colors = self.shifted_cmap(data/data.max())
        b=helper.box("base",size=[rspace,rspace,rspace])[0]
        count=0
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                for k in range(data.shape[2]):
                    a=helper.newInstance("map"+str(i)+"_"+str(j)+"_"+str(k),b,location=[i*rspace,j*rspace,k*rspace])
                    helper.changeObjColorMat(a,colors[i,j,k][:3])
                    p=(count/float(data.shape[0]*data.shape[1]*data.shape[2]))*100.0
                    helper.progressBar(progress=p,label=str(p))

    def drawAsColoredInstance2D(self,helper,data,rspace):#or instance
        if helper is None :
            return
        colors = self.shifted_cmap(data/data.max())
        b=helper.box("base",size=[rspace,rspace,rspace])[0]
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

    def makeImage(self,helper,name,data,w,imfilename):
        if helper is None :
            return
        colors = self.shifted_cmap(data/data.max())
        from PIL import Image
        im = Image.fromarray(np.uint8(colors*255),mode = None)
        im.save(imfilename)
        return im

    def drawAsColoredPlane(self,helper,name,ax,data,w,imfilename):#or instance
        if helper is None :
            return
        p = helper.getObject(name)
        r=[0.,0.,0.]
        tr=[0,0,0]
        if p is None :
            if ax == 0 :
                #this may be only valid for C4D
                ax = 1#-X
                tr=[0,0,-w/2.0]
                r=[0.,np.pi,0.]
            elif ax == 1 :
                ax = 5#-Z
                tr=[w/2.,0,0]
                r=[-np.pi/2.0,0.,0.]
            elif ax == 2 :
                ax = 2#+Y
                tr=[0,-w/2.,0]
                r=[0.,np.pi/2.0,0.]
            p,mpl = helper.plane(name,center = [0,0,0],size=[w,w],parent=None,axis=ax)
            helper.rotateObj(p,r)
            helper.setTranslation(p,tr)
        im=self.makeImage(helper,name,data,w,imfilename)
        mat = helper.createTexturedMaterial(name+"_Mat",imfilename)
        #assign the material to the plane
        helper.assignMaterial(p,mat,texture=True)
        #could use c4d viewport viewport[c4d.symbols.BASEDRAW_DATA_PICTURE] = "/tmp/arpmv.jpg"
        return p

    def drawAPlane(self,helper,output,ax):
        self.drawAsColoredPlane(helper,"gauss"+str(ax),ax,np.sum(self.data_conv,axis=ax),self.W*2.,output+"_"+str(ax)+".jpg")
        self.drawAsColoredPlane(helper,"posX"+str(ax),ax,np.sum(self.data,axis=ax),self.W*2.,output+"_p_"+str(ax)+".jpg")


    def drawAsPlaneXYZ(self,helper,output):
        if helper is None :
            return
        for i in range(3):
            self.drawAPlane(helper,output,i)

import upy
uiadaptor = upy.getUIClass()

class fluoSimGui(uiadaptor):
    def CreateLayout(self):
        self._createLayout()
        return 1

    def Command(self,*args):
        self._command(args)
        return 1

    def setup(self,**kw):
#        self.subdialog = True
        self.block = True
#        self.scrolling = True
        self.title = "Fluoresence simulation"
        self.parent = None
        self.bbox = 0
        self.space= 0
        self.rspace=0
        self.res=0
        self.psf_type=""
        self.helper = autopack.helper
        self.env=None
        if "environment" in kw :
            self.env = kw["environment"]
        self.fluos = None
        if "parent" in kw :
            self.parent = kw["parent"]
        #should load the user saved preferences
        self.SetTitle(self.title)
        self.initWidget()
        self.setupLayout()

    def initWidget(self, ):
        #values=None,psf_type="gauss",
        #boundingBox=4195,res=400.0,space=50.0
        self.Widget={}
        self.Widget["label"]={}
        self.Widget["options"]={}

        self.Widget["label"]["psf_type"] = self._addElemt(name="psf_typeLabel",label="choose a PSF",width=120)
        self.Widget["options"]["psf_type"] = self._addElemt(name="psf_type",
                                    value=["gauss","lorentz"],
                                    width=100,height=10,action=None,
                                    variable=self.addVariable("int",0),
                                    type="pullMenu",)

        self.Widget["label"]["boundingBox"] = self._addElemt(name="boundingBoxLabel",
                                            label="box length",width=120)
        self.Widget["options"]["boundingBox"] = self._addElemt(name="boundingBox",
                                value=4000.0, width=100,height=10,action=None,
                                    mini=1,maxi=10000,type="inputFloat")

        self.Widget["label"]["HWFM"] = self._addElemt(name="HWFMLabel",label="Fluo resolution (HWFM) :",width=120)
        self.Widget["options"]["HWFM"] = self._addElemt(name="HWFM",
                                value=400.0, width=100,height=10,action=None,
                                    mini=1,maxi=1000,type="inputFloat")

        self.Widget["label"]["space"] = self._addElemt(name="spaceLabel",label="grid spacing:",width=120)
        self.Widget["options"]["space"] = self._addElemt(name="space",
                                value=100.0, width=100,height=10,action=None,
                                    mini=1,maxi=1000,type="inputFloat")

        self.Widget["label"]["rspace"] = self._addElemt(name="rspaceLabel",label="pixel size (A):",width=120)
        self.Widget["options"]["rspace"] = self._addElemt(name="rscpace",
                                value=200.0, width=100,height=10,action=None,
                                    mini=1,maxi=1000,type="inputFloat")

        liste_input = ["selection","object name"]
        label = "Object name:"
        if self.env is not None :
            liste_input = ["selection","ingredient name","object name"]
            label = "Object/Ingredient name:"
        self.Widget["label"]["fluorophore"] = self._addElemt(name="fluorophoreLsbel",label="Fluorophore positions:",width=120)
        self.Widget["options"]["fluorophore"] = self._addElemt(name="psf_type",
                                    value=liste_input,
                                    width=100,height=10,action=None,
                                    variable=self.addVariable("int",0),
                                    type="pullMenu",)
        self.Widget["label"]["selection"] = self._addElemt(name="selectionLsbel",label=label,width=120)
        self.Widget["options"]["selection"] = self._addElemt(name="selection",
                                    value="",
                                    width=100,height=10,action=None,
                                    variable=self.addVariable("str",""),
                                    type="inputStr",)
        self.Widget["label"]["offset"] = self._addElemt(name="offsetLsbel",label="Fluorophore offset position:",width=120)
        self.Widget["options"]["offset"] = self._addElemt(name="offset",
                                    value="",
                                    width=100,height=10,action=None,
                                    variable=self.addVariable("str",""),
                                    type="inputStr",)

        self.Widget["label"]["fluos"] = self._addElemt(name="fluoLabel",label="",width=120)

        #define the buttons
        self.BTN={}
        self.BTN["close"]=self._addElemt(name="Close",width=50,height=10,
                         action=self.close,type="button",icon=None,
                                     variable=self.addVariable("int",0))
        self.BTN["setup"]=self._addElemt(name="Setup",width=50,height=10,
                         action=self.Setup,type="button",icon=None,
                                     variable=self.addVariable("int",0))
        self.BTN["generate"]=self._addElemt(name="Generate",width=50,height=10,
                         action=self.Generate,type="button",icon=None,
                                     variable=self.addVariable("int",0))

        self.BTN["displayXYZ"]=self._addElemt(name="displayXYZ",width=50,height=10,
                         action=self.displayXYZ,type="button",icon=None,
                                     variable=self.addVariable("int",0))
#        self.BTN["displayCsutom"]=self._addElemt(name="Generate",width=50,height=10,
#                         action=self.Generate,type="button",icon=None,
#                                     variable=self.addVariable("int",0))

    def setupLayout(self, ):
        self._layout = []
        for wname in ["psf_type","boundingBox","space","HWFM","rspace","fluorophore",
                      "selection","offset"]:
            widget =[self.Widget["label"][wname],self.Widget["options"][wname]]
            self._layout.append(widget)
        self._layout.append([self.BTN["setup"],self.BTN["generate"]])
        self._layout.append([self.Widget["label"]["fluos"],])
        self._layout.append([self.BTN["displayXYZ"],])
        self._layout.append([self.BTN["close"]])

    def getAllVal(self,):
        self.bbox = self.getVal(self.Widget["options"]["boundingBox"])
        self.space= self.getVal(self.Widget["options"]["space"])
        self.rspace=self.getVal(self.Widget["options"]["rspace"])
        self.res=self.getVal(self.Widget["options"]["HWFM"])
        self.psf_type=self.getVal(self.Widget["options"]["psf_type"])

    def Setup(self,*args,**kw):
        #get value from widget
        self.getAllVal()
        if self.fluos is None :
            self.fluos = fluoSim(values=None,psf_type=self.psf_type,
                 boundingBox=self.bbox,res=self.res,space=self.space )
            self.setVal(self.Widget["label"]["fluos"],"fluos initialized")
        else :
            self.fluos.setup(values=None,psf_type=self.psf_type,boundingBox=self.bbox,
                   res=self.res,space=self.space)
            self.setVal(self.Widget["label"]["fluos"],"fluos updated")

    def generate(self,data):
        if not len(data) :
            return
        self.fluos.setValuesAccurate(values=data)
        gauss_conv=self.fluos.getConvolution()
        g=self.fluos.setAtScale(self.rspace)
        self.fluos.data_conv = g

    def getData(self,):
        option = self.getVal(self.Widget["options"]["fluorophore"])
        print ("option for fluorophore is ",option)
        #"selection","ingredient name","object name"
        #if selection get the selection position
        pos=[]
        rot=[]
        s_name=self.getVal(self.Widget["options"]["selection"])
        offset=self.getVal(self.Widget["options"]["offset"])
        if option == "ingredient name":
            if self.env == None :
                return
            ingr = self.env.getIngrFromName(s_name)
            self.env.collectResultPerIngredient()
            res=ingr.results#[pos,rot]
            #sph_pos = h.compartments[0].surfaceRecipe.ingredients[0].positions[-1]
            if offset != "" :
                pos = [ingr.transformPoints(r[0],r[1], [eval(offset)])[0] for r in res ]
            else :
                pos = [r[0] for r in res]
            return pos
        else :
            if option == "selection" :
                selections = self.helper.getCurrentSelection()
                if not len(selections):
                    return []
                o=selections[0]
            elif option == "object name" :
                o = self.helper.getObject(s_name)
                if o is None :
                    return []
            if self.helper.getType(o) == self.helper.EMPTY :
                #get child position
                print ("null object")
                pos=[]
                rot=[]
                for ch in self.helper.getChilds(o):
                    print ("object",self.helper.getName(ch))
                    p,r=self.helper.getPropertyObject(ch,pos=True,rotation=True)
                    pos.append(p)
                    rot.append(r)
            elif self.helper.getType(o) == self.helper.POLYGON :
                #rotation should be the normal ?
                f,v,vn = self.helper.DecomposeMesh(o,edit=False,copy=False,tri=True,transform=True)
                vn,fn,area=self.helper.getFaceNormalsArea(v, f)
                pos = v[:]
                #rot align to vn
            if offset != "" :
               pos = [self.helper.ApplyMatrix([eval(offset)],r) for r in rot]
            return pos
        #if selection Null, get children position
        #if ingredient name, get ingredient instance position
        #if objectName get the position, if objec Null get child poition
        #if polygin should we get the vertices as position ?

    def Generate(self,*args,**kw):
        #where does come from the data...
        #can be ingredient name, object
#        data = fluos.parseAPresults(filename,recipe_name = rname,ingr_name=iname)
        if self.fluos is None:
            self.Setup()
        data = self.getData()
        print data
        if not len(data):
            self.setVal(self.Widget["label"]["fluos"],"fluos failed")
            return
        self.generate(np.array(data))
        self.setVal(self.Widget["label"]["fluos"],"fluos generated")

    def displayXYZ(self,*args,**kw):
        if self.fluos is None:
            self.Generate()
        output=autopack.cache_results+os.sep
        self.fluos.drawAsPlaneXYZ(self.helper,output)
