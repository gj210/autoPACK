# -*- coding: utf-8 -*-
"""
Created on Mon May  6 22:58:44 2013

@author: ludo
"""
import os
import math
import numpy
import csv
import json

from time import time
try :
    import matplotlib
    # Force matplotlib to not use any Xwindows backend.
    #matplotlib.use('Agg')   
    from matplotlib import pylab
    from matplotlib import pyplot
    import matplotlib.mlab as mlab
    from matplotlib.patches import Circle
except :
    matplotlib=None
    pylab = None
    pyplot = None
    Circle = None
    mlab = None
    
import autopack
from autopack.GeometryTools import GeometriTools,Rectangle
import Image

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
    ind=numpy.nonzero(dot < 0)
    angles[ind]*=-1.0
    return angles

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
        
        
class AnalyseAP:
    def __init__(self, env=None, viewer=None, result_file=None):
        self.env=None
        self.smallest=99999.0
        self.largest = 0.0
        if env :
            self.env = env
            self.smallest,self.largest = self.getMinMaxProteinSize()
        self.afviewer = viewer
        self.helper = None
        if viewer :
            self.helper = self.afviewer.vi
        self.resutl_file = result_file
        self.center=[0,0,0]
        self.bbox=[[0,0,0],[1,1,1]]
        self.g=GeometriTools()
        self.g.Resolution = 1.0#or grid step?
        self.current_pos=None
        self.current_distance=None
        autopack._colors = None
        
    def getMinMaxProteinSize(self):
        smallest=999999.0
        largest= 0.0
        for organelle in self.env.compartments:
            mini, maxi = organelle.getMinMaxProteinSize()
            if mini < smallest:
                smallest = mini
            if maxi > largest:
                largest = maxi

        if self.env.exteriorRecipe:
            mini, maxi = self.env.exteriorRecipe.getMinMaxProteinSize()
            if mini < smallest:
                smallest = mini
            if maxi > largest:
                largest = maxi
        return smallest,largest
        
    def getPositionsFromResFile(self,):
        #could actually restore file using histoVol.
        #or not
        #need to parse apr file here anyway
        return []
        
    def getPositionsFromObject(self,parents): 
        positions=[]
        for parent in parents:
            obparent = self.helper.getObject(parent)
            childs = self.helper.getChilds(obparent)
            for ch in childs:
                ingr_name = self.helper.getName(ch)
                meshp = self.helper.getObject("Meshs_"+ingr_name.split("_")[0])
                if meshp is None :
                    c = self.helper.getChilds(ch)
                    if not len(c) :
                        continue                        
                    meshpchilds = self.helper.getChilds(c[0])#continue #should get sphere/cylnder parent ?
                else :
                    meshpchilds = self.helper.getChilds(meshp)
                for cc in meshpchilds:
                    pos = self.helper.ToVec(self.helper.getTranslation(cc))
                    positions.append(pos)
        return positions
        
    def getDistanceFrom(self,target,parents=None,**options):
        """
        target : name or host object target or target position
        parent : name of host parent object for the list of object to measre distance from
        objects : list of object or list of points
        """
        #get distance from object to the target.
        #all object are in h.molecules and orga.molecules
        #get options
        targetPos = [0,0,0]
        usePoint = False
        threshold = 99999.
        if "usePoint" in options:
            usePoint = options["usePoint"]        
        if "threshold" in options:
            threshold = options["threshold"]
        if type(target) == list or type(target) == tuple:
            targetPos = target
        elif type(target) == unicode or type(target) == str : 
            o = self.helper.getObject(target)
            if o is not None :
                targetPos = self.helper.ToVec(self.helper.getTranslation(o)) #hostForm
        else :
            o = self.helper.getObject(target)
            if o is not None :
                targetPos = self.helper.ToVec(self.helper.getTranslation(o)) #hostForm
        listeObjs=[]
        listeDistances = []
        listeCenters =[]
        if self.resutl_file is None :
            if parents is None and self.resutl_file is None:
                listeParent = [self.env.name+"_cytoplasm"]
                for o in self.env.compartments :
                    listeParent.append(o.name+"_Matrix")
                    listeParent.append(o.name+"_surface")
            elif parents is not None and self.resutl_file is None:       
                listeParent = parents
            listeCenters =self.getPositionsFromObject(listeParent)
        else :
            #use data from file
            listeCenters =self.getPositionsFromResFile(listeParent) 

        delta = numpy.array(listeCenters)-numpy.array(targetPos)
        delta *= delta
        distA = numpy.sqrt( delta.sum(1) )           
        return distA

    def getClosestDistance(self,parents=None,**options):
        if self.resutl_file is None :
            if parents is None and self.resutl_file is None:
                listeParent = [self.env.name+"_cytoplasm"]
                for o in self.env.compartments :
                    listeParent.append(o.name+"_Matrix")
                    listeParent.append(o.name+"_surface")
            elif parents is not None and self.resutl_file is None:       
                listeParent = parents
            listeCenters =self.getPositionsFromObject(listeParent)
        else :
            #use data from file
            listeCenters =self.getPositionsFromResFile(listeParent) 
        #is the distance in the result array ?
        listeDistance=numpy.zeros(len(listeCenters))+99999
        for i in range(len(listeCenters) ):
            for j in range(i+1,len(listeCenters)):
                #should use point
                d = self.helper.measure_distance(listeCenters[i],listeCenters[j])
                if d < listeDistance[i] :
                    listeDistance[i] = d
        return listeDistance



    def displayDistance(self,ramp_color1=[1,0,0],ramp_color2=[0,0,1],
                        ramp_color3=None,cutoff=60.0):
        distances = numpy.array(self.env.grid.distToClosestSurf[:])
        positions = numpy.array(self.env.grid.masterGridPositions[:])
        #map the color as well ?
        from upy import colors as col
        from DejaVu.colorTool import Map
        ramp = col.getRamp([ramp_color1,ramp_color2],size=255)   #color
        mask = distances > cutoff
        ind=numpy.nonzero(mask)[0]
        distances[ind]=cutoff
        mask = distances < 0#-cutoff
        ind=numpy.nonzero(mask)[0]
        distances[ind]=0#cutoff   
        newd=numpy.append(distances,cutoff)
        colors = Map(newd, ramp)[:-1]#1D array of the grid x,y,1
#        colors = Map(distances, ramp)
#        sphs = self.helper.Spheres("distances",vertices=numpy.array(positions),radii=distances,colors=colors)
        base=self.helper.getObject(self.env.name+"distances_base")
        if base is None :
            base=self.helper.Sphere(self.env.name+"distances_base")[0]
        p=self.helper.getObject(self.env.name+"distances")
        if p is not None :
            self.helper.deleteObject(p)#recursif?
        p = self.helper.newEmpty(self.env.name+"distances")
        sphs=self.helper.instancesSphere(self.env.name+"distances",positions,distances,base,colors,None,parent=p)
        #can use cube also 

#    def displayFillBoxCubeGrid(self,):
#        bb_insidepoint = self.env.grid.getPointsInCubeFillBB(self.env.fbox, [0,0,0], 1.0)[:]#center and radius ?3,runTime=self.runTimeDisplay
#        positions = numpy.array(self.env.grid.masterGridPositions[bb_insidepoint])    
        
    def displayDistanceCube(self,ramp_color1=[1,0,0],ramp_color2=[0,0,1],
                        ramp_color3=None,cutoff=60.0):
#        if self.env.fbox is not None :
#            bb_insidepoint = self.env.grid.getPointsInCubeFillBB(self.env.fbox, [0,0,0], 1.0)[:]#center and radius ?3,runTime=self.runTimeDisplay
#            distances = numpy.array(self.env.grid.distToClosestSurf[bb_insidepoint])
#            positions = numpy.array(self.env.grid.masterGridPositions[bb_insidepoint])    
#        else :
        distances = numpy.array(self.env.grid.distToClosestSurf[:])
        positions = numpy.array(self.env.grid.masterGridPositions[:])
        #map the color as well ?
        from upy import colors as col
        from DejaVu.colorTool import Map
        ramp = col.getRamp([ramp_color1,ramp_color2],size=255)   #color
        mask = distances > cutoff
        ind=numpy.nonzero(mask)[0]
        distances[ind]=cutoff
        mask = distances < 0#-cutoff
        ind=numpy.nonzero(mask)[0]
        distances[ind]=0#cutoff   
        newd=numpy.append(distances,cutoff)
        colors = Map(newd, ramp)[:-1]#1D array of the grid x,y,1
#        colors = Map(distances, ramp)
#        sphs = self.helper.Spheres("distances",vertices=numpy.array(positions),radii=distances,colors=colors)
        base=self.helper.getObject(self.env.name+"distances_base_cube")
        if base is None :
#            base=self.helper.Sphere(self.env.name+"distances_base")[0]
            size = self.env.grid.gridSpacing
            base=self.helper.box(self.env.name+"distances_base_cube",
                                 center=[0.,0.,0.],size=[size,size,size])[0]
        parent_cube=self.helper.getObject(self.env.name+"distances_cubes")
        if parent_cube is not None :
            self.helper.deleteObject(parent_cube)#recursif?
        parent_cube = self.helper.newEmpty(self.env.name+"distances_cubes")
        #sphs=self.helper.instancesSphere(self.env.name+"distances_cubes",positions,distances,base,colors,None,parent=p)
        #can use cube also 
        for i,p in enumerate(positions):
            mat = self.helper.addMaterial("matcube"+str(i),colors[i])
            c = self.helper.newInstance(self.env.name+"distances_cubes_"+str(i),base,
                                      location=p,material=mat,parent=parent_cube)

    def displayDistancePlane(self,ramp_color1=[1,0,0],ramp_color2=[0,0,1],
                        ramp_color3=None,cutoff=60.0):
        #which axis ?
        distances = numpy.array(self.env.grid.distToClosestSurf[:])
        max_distance= max(distances)
        min_distance= min(distances)
        
#        positions = numpy.array(self.env.grid.masterGridPositions[:])
        #positions = self.env.grid.masterGridPositions[:]
        #map the color as well ?
        from upy import colors as col
        from DejaVu.colorTool import Map
        ramp = col.getRamp([ramp_color1,ramp_color2],size=255)   #color
        mask = distances > cutoff
        ind=numpy.nonzero(mask)[0]
        distances[ind]=cutoff
        mask = distances < 0#-cutoff
        ind=numpy.nonzero(mask)[0]
        distances[ind]=0#cutoff 
        newd=numpy.append(distances,cutoff)
        colors = Map(newd, ramp)[:-1]#1D array of the grid x,y,1
        autopack._colors = colors
#        sphs = self.helper.Spheres("distances",vertices=numpy.array(positions),radii=distances,colors=colors)
        p=self.helper.getObject(self.env.name+"distances")
        if p is not None :
            self.helper.deleteObject(p)#recursif?
        p = self.helper.newEmpty(self.env.name+"distances_p")
        #sphs=self.helper.instancesSphere(self.env.name+"distances_cubes",positions,distances,base,colors,None,parent=p)
        #can use cube also 
        print ("grid is ",self.env.grid.nbGridPoints)
        print ("colors shape is ",colors.shape)
        d = numpy.array(self.env.grid.boundingBox[0]) - numpy.array(self.env.grid.boundingBox[1])
        p,mpl = self.helper.plane(self.env.name+"distances_plane",
                                  center = self.env.grid.getCenter(),
                                size=[math.fabs(d[0]),math.fabs(d[1])],
                                    parent=p)
        self.helper.rotateObj(p,[0,0,-math.pi/2.0])
        filename = autopack.cache_results+os.sep+self.env.name+"distances_plane_texture.png"
        c=colors.reshape((self.env.grid.nbGridPoints[0],
                          self.env.grid.nbGridPoints[1],
                          self.env.grid.nbGridPoints[2],3))
        
        im = Image.fromstring("RGB",(c.shape[0],c.shape[1]),numpy.uint8(c*255.0).tostring())
        im.save(str(filename))
        mat = self.helper.createTexturedMaterial(self.env.name+"planeMat",str(filename))
        #assign the material to the plane
        self.helper.assignMaterial(p,mat,texture=True)


    def writeJSON(self,filename,data):
        with open(filename, 'w') as fp :#doesnt work with symbol link ?
            json.dump(data,fp,indent=4, separators=(',', ': '))#,indent=4, separators=(',', ': ')
                
    def loadJSON(self,filename):
        with open(filename) as data_file:    
            data = json.load(data_file)  
        return data

    def grabResultFromJSON(self,n):
        ingrrot={}
        ingrpos={}
        for i in range(n):
            with open("results_seed_"+str(i)+".json") as data_file:  
                data = json.load(data_file)
            for recipe in data:
                for ingrname in data[recipe]:
                    for k in range(len(data[recipe][ingrname]['results'])):
                        if ingrname not in ingrrot:
                            ingrrot[ingrname]=[]
                            ingrpos[ingrname]=[]                            
                        ingrrot[ingrname].append(data[recipe][ingrname]['results'][k][1])
                        ingrpos[ingrname].append(data[recipe][ingrname]['results'][k][0])
            print i
        return ingrpos,ingrrot

    def grabResultFromTXT(self,n,doanalyze=False):
        from autopack import transformation as t
        ingrrot={}
        ingrpos={}
        ingrpos3=[]
        ingrrot3=[]
        for i in range(1000):
            print i    
            files=open("results_seed_"+str(i)+".txt","r")
            lines = files.readlines()
            files.close()
            for l in lines :
                l=l.replace("<"," ").replace(">"," ")
                elem=l.split()
                ingrname=elem[-5]
                if ingrname not in ingrrot:
                    ingrrot[ingrname]=[]
                    ingrpos[ingrname]=[]                            
                ingrrot[ingrname].append(eval(elem[2])) 
                ingrpos[ingrname].append(eval(elem[0]))
        for ingrname in ingrrot:
            ingrrot[ingrname] = [numpy.array(m).reshape((4,4)) for m in ingrrot[ingrname]]
        if doanalyze :
            for ingrname in ingrrot:
                eulers3 = [t.euler_from_matrix(m, 'rxyz') for m in ingrrot[ingrname]]
                e3=numpy.degrees(numpy.array(eulers3)).transpose()
                numpy.savetxt(ingrname+"_euler_X.csv", numpy.array(e3[0]), delimiter=",") 
                numpy.savetxt(ingrname+"_euler_Y.csv", numpy.array(e3[1]), delimiter=",") 
                numpy.savetxt(ingrname+"_euler_Z.csv", numpy.array(e3[2]), delimiter=",") 
                self.histo(e3[0],ingrname+"_euler_X.png",bins=12,size=max(e3[0]))
                self.histo(e3[1],ingrname+"_euler_Y.png",bins=12,size=max(e3[1]))
                self.histo(e3[2],ingrname+"_euler_Z.png",bins=12,size=max(e3[2]))
#                ingrpositions,distA,angles3=self.getDistanceAngle(ingrpos3, ingrrot3)
#                numpy.savetxt(ingrname+"_angle_X.csv", numpy.array(angles3[1]), delimiter=",") 
#                numpy.savetxt(ingrname+"_angle_Y.csv", numpy.array(angles3[2]), delimiter=",") 
#                numpy.savetxt(ingrname+"_angle_Z.csv", numpy.array(angles3[3]), delimiter=",") 
#                self.histo(angles3[1],ingrname+"_angle_X.png",bins=12,size=max(angles3[1]))
#                self.histo(angles3[2],ingrname+"_angle_Y.png",bins=12,size=max(angles3[2]))
#                self.histo(angles3[3],ingrname+"_angle_Z.png",bins=12,size=max(angles3[3]))
        return ingrpos,ingrrot
                            
    #should take any type of list...
    def save_csv(self,data,filename=None):
        if filename is None :
            filename = "output.csv"
        resultFile = open(filename,'wb')
        wr = csv.writer(resultFile, dialect='excel')
        #wr.writerows(data) list of list ?
        #resultFile.close()
        for item in data: wr.writerow([item,]) 
        resultFile.close()

    def rectangle_circle_area(self,bbox,center,radius):
        #http://www.eex-dev.net/index.php?id=100        
        #[[0.,0,0],[1000,1000,1]]
        rect=Rectangle(bbox[0][0],bbox[1][0],bbox[0][1],bbox[1][1])#top,bottom, right, left
#        rect=Rectangle(bbox[1][1],bbox[0][1],bbox[1][0],bbox[0][0])#top,bottom, right, left        
        m=[center[0],center[1]]
        r=radius
        area = 0.0
        chs = self.g.check_sphere_inside(rect,m,r)
        if chs :
#            print "sphere go outside ",r 
            ch=self.g.check_rectangle_oustide(rect,m,r)
            if ch :
                leftBound,rightBound = self.g.getBoundary(rect,m,r)
                area = self.g.get_rectangle_cercle_area(rect,m,r,rightBound,leftBound)
#                print area,leftBound,rightBound
        return area


    def getAxeValue(self,ingrname,axe=0):
        ingrpositions=[self.env.molecules[i][0][axe] for i in xrange(len(self.env.molecules)) if self.env.molecules[i][2].name == ingrname]
        return ingrpositions

    def getAxesValues(self,positions):
        pp = numpy.array(positions).transpose()
        px=pp[0]
        py=pp[1]
        pz=pp[2]
        return px,py,pz

    def getDistance(self,ingrname, center):
        distA=[]
        ingrpositions=[self.env.molecules[i][0] for i in xrange(len(self.env.molecules)) if self.env.molecules[i][2].name == ingrname]
        ingrpositions=numpy.array(ingrpositions)
        if len(ingrpositions) :
            delta = numpy.array(ingrpositions)-numpy.array(center)
            delta *= delta
            distA = numpy.sqrt( delta.sum(1) ).tolist()
        return ingrpositions,distA
        
    def getDistanceAngle(self,ingr, center):
        #need matrix to euler? then access and plot them?
        #also check the measure angle one
        angles=[]
        distA=[]
        ingrpositions=[self.env.molecules[i][0] for i in xrange(len(self.env.molecules)) if self.env.molecules[i][2].name == ingr.name]
        ingrpositions=numpy.array(ingrpositions)
        ingrrotation=[self.env.molecules[i][1] for i in xrange(len(self.env.molecules)) if self.env.molecules[i][2].name == ingr.name]
        ingrrotation=numpy.array(ingrrotation)
        if len(ingrpositions) :
            #thats m
            orientationX = numpy.array([autopack.helper.ApplyMatrix([[1,0,0],],m)[0] for m in ingrrotation])
            orientationY = numpy.array([autopack.helper.ApplyMatrix([[0,1,0],],m)[0] for m in ingrrotation])
            orientationZ = numpy.array([autopack.helper.ApplyMatrix([[0,0,1],],m)[0] for m in ingrrotation])
#            orientationX = numpy.array([m[0][:3] for m in ingrrotation])
#            orientationY = numpy.array([m[1][:3] for m in ingrrotation])
#            orientationZ = numpy.array([m[2][:3] for m in ingrrotation])
            delta = numpy.array(ingrpositions)-numpy.array(center)
            #lets do it on X,Y,Z and also per positions ?            
            anglesX=signed_angle_between_vectors([0,0,1],m[0][:3],-delta,directed=False,axis=1)
            anglesY=signed_angle_between_vectors([0,0,1],m[1][:3],-delta,directed=False,axis=1)
            anglesZ=signed_angle_between_vectors([1,0,0],m[2][:3],-delta,directed=False,axis=1)
            delta *= delta
            distA = numpy.sqrt( delta.sum(1) ).tolist()  
            angles = numpy.array([distA,anglesX,anglesY,anglesZ])
        return ingrpositions,distA,numpy.degrees(angles)

    def getVolumeShell(self,bbox,radii,center):
        #rectangle_circle_area
        volumes=[]
        box_size0 = bbox[1][0] - bbox[0][0] 
        for i in range(len(radii)-1):
            r1=radii[i]
            r2=radii[i+1]
            v1 = self.g.calc_volume(r1, box_size0 / 2.)
            v2 = self.g.calc_volume(r2, box_size0 / 2.)
#            if v1 == 0 or v2 == 0 :
#                volumes.append((4./3.)*numpy.pi*(numpy.power(r2,3)-numpy.power(r1, 3)))
#            else :
            volumes.append(v2-v1)
        return volumes
        
    def rdf_3d(self,ingr):
        #see for intersection volume here http://crowsandcats.blogspot.com/2013/04/cube-sphere-intersection-volume.html
        #and here http://crowsandcats.blogspot.com/2013/05/extending-radial-distributions.html
        #will require scipy...worth it ?
        #should be pairewise distance ? or not ?
        distances = numpy.array(self.env.distances[ingr.name])
        basename = self.env.basename
        numpy.savetxt(basename+ingr.name+"_pos.csv", numpy.array(self.env.ingrpositions[ingr.name]), delimiter=",") 
        self.histo(distances,basename+ingr.name+"_histo.png")
        numpy.savetxt(basename+ingr.name+"_distances.csv", numpy.array(distances), delimiter=",") 
        # the bin should be not less than the biggest ingredient radius
        #b=int(distances.max()/self.largest)
        b=100
        #bin_edges = numpy.arange(0, min(box_size) / 2, bin_width)
        new_rdf, edges = numpy.histogramdd(distances, bins=b, range=[(distances.min(), distances.max())])
        radii = edges[0]
        #from http://isaacs.sourceforge.net/phys/rdfs.html
        dnr=new_rdf
        N=len(distances)
        V=self.env.grid.nbGridPoints[0]*self.env.grid.nbGridPoints[1]*self.env.grid.nbGridPoints[2]*self.env.grid.gridSpacing**3
        density = 1#len(x)/1000.0**2
#        Vshell = (4./3.)*numpy.pi*(numpy.power(radii[1:],3)-numpy.power(radii[:-1], 3))
        Vshell = numpy.array(self.getVolumeShell(self.bbox,radii,self.center))
        gr = (dnr*V)/(N*Vshell)
        numpy.savetxt(basename+ingr.name+"_rdf.csv", numpy.array(gr), delimiter=",")
        self.plot(gr,radii[:-1],basename+ingr.name+"_rdf.png")

    def getAreaShell(self,bbox,radii,center):
        #rectangle_circle_area
        areas=[]
        for i in range(len(radii)-1):
            r1=radii[i]
            r2=radii[i+1]
            area1=self.rectangle_circle_area(bbox,center,r1)
            area2=self.rectangle_circle_area(bbox,center,r2)
            if area1 == 0 or area2 == 0 :
                areas.append(numpy.pi*(numpy.power(r2,2)-numpy.power(r1, 2)))
            else :
                areas.append(area2-area1)
        return areas
        
    def rdf_2d(self,ingr):
        distances = numpy.array(self.env.distances[ingr.name])
        basename = self.env.basename
        numpy.savetxt(basename+ingr.name+"_pos.csv", numpy.array(self.env.ingrpositions[ingr.name]), delimiter=",") 
        self.histo(distances,basename+ingr.name+"_histo.png")
        numpy.savetxt(basename+ingr.name+"_distances.csv", numpy.array(distances), delimiter=",") 
        # the bin should be not less than the biggest ingredient radius
#        b=int(distances.max()/self.largest)
        b=100
        new_rdf, edges = numpy.histogramdd(distances)#, bins=b, range=[(distances.min(), distances.max())],normed=0)
        radii = edges[0]
#        r=radii.tolist()
#        r.insert(0,0.0)
#        radii = numpy.array(r)
#        rdf= new_rdf.tolist()
#        rdf.insert(0,0)
#        new_rdf = numpy.array(rdf)
        #from http://isaacs.sourceforge.net/phys/rdfs.html
        dnr=new_rdf[:]
        N=len(distances)
        V=self.env.grid.nbGridPoints[0]*self.env.grid.nbGridPoints[1]*self.env.grid.gridSpacing**2
        density = 1#len(x)/1000.0**2
        Vshell = numpy.array(self.getAreaShell(self.bbox,radii,self.center))       
#        print Vshell
#        Vshell1 = numpy.pi*density*(numpy.power(radii[1:],2)-numpy.power(radii[:-1], 2))
#        print Vshell1   
#        print radii
        gr = (dnr*V)/(N*Vshell)
        numpy.savetxt(basename+ingr.name+"_rdf.csv", numpy.array(gr), delimiter=",")
        self.plot(gr,radii[:-1],basename+ingr.name+"_rdf.png")
        #simpl approach Ni/Areai
        G=dnr/Vshell
        numpy.savetxt(basename+ingr.name+"_rdf_simple.csv", numpy.array(G), delimiter=",")
        self.plot(numpy.array(G),radii[:-1],basename+ingr.name+"_rdf_simple.png")
        print G
        
    def axis_distribution_total(self,all_positions):
        basename = self.env.basename
        numpy.savetxt(basename+"total_pos.csv", numpy.array(all_positions), delimiter=",") 
        px,py,pz = self.getAxesValues(all_positions)
#        m1=numpy.nonzero( numpy.logical_and(
#               numpy.greater_equal(px, 0.), numpy.less_equal(px, 1000.0)))
#        m2=numpy.nonzero( numpy.logical_and(
#               numpy.greater_equal(py, 0.), numpy.less_equal(py, 1000.0)))
        self.histo(px,basename+"total_histo_X.png",bins=10)
        self.histo(py,basename+"total_histo_Y.png",bins=10)
        self.histo(pz,basename+"total_histo_Z.png",bins=10)

    def axis_distribution(self,ingr):
        basename = self.env.basename
        px,py,pz = self.getAxesValues(self.env.ingrpositions[ingr.name])
        self.histo(px,basename+ingr.name+"_histo_X.png",bins=10)
        self.histo(py,basename+ingr.name+"_histo_Y.png",bins=10)
        self.histo(pz,basename+ingr.name+"_histo_Z.png",bins=10)
        #do it for all ingredient cumulate?

    def occurence_distribution (self,ingr):
        basename = self.env.basename
        occ = self.env.occurences[ingr.name]
        self.simpleplot(range(len(occ)),occ,basename+ingr.name+"_occurence.png")
        
    def correlation(self,ingr):
        basename = self.env.basename
        posxyz = numpy.array(self.env.ingrpositions[ingr.name]).transpose()
        g_average, radii, x,y,z=self.PairCorrelationFunction_3D(posxyz,1000,900,100)
        self.plot(g_average,radii,basename+ingr.name+"_corr.png")
        
    def PairCorrelationFunction_3D(self,data,S,rMax,dr):
        """Compute the three-dimensional pair correlation function for a set of
    spherical particles contained in a cube with side length S. This simple
    function finds reference particles such that a sphere of radius rMax drawn
    around the particle will fit entirely within the cube, eliminating the need
    to compensate for edge effects. If no such particles exist, an error is
    returned. Try a smaller rMax...or write some code to handle edge effects! ;)
    Arguments:
    x an array of x positions of centers of particles
    y an array of y positions of centers of particles
    z an array of z positions of centers of particles
    S length of each side of the cube in space
    rMax outer diameter of largest spherical shell
    dr increment for increasing radius of spherical shell
    
    Returns a tuple: (g, radii, interior_x, interior_y, interior_z)
    g(r) a numpy array containing the correlation function g(r)
    radii a numpy array containing the radii of the
    spherical shells used to compute g(r)
    interior_x x coordinates of reference particles
    interior_y y coordinates of reference particles
    interior_z z coordinates of reference particles
    """
        from numpy import zeros, sqrt, where, pi, average, arange, histogram
        x=data[0]
        y=data[1]
        z=data[2]
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
            raise RuntimeError ("No particles found for which a sphere of radius rMax\
    will lie entirely within a cube of side length S. Decrease rMax\
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


    def histo(self,distances,filename,bins=100,size=1000.0):
        pylab.clf()
        mu, sigma = numpy.mean(distances) , numpy.std(distances)
        ## the histogram of the data
#        b=numpy.arange(distances.min(), distances.max(), 2)
#        n, bins, patches = pyplot.hist(distances, bins=bins, normed=1, facecolor='green')#, alpha=0.75)
        y,binEdges=numpy.histogram(distances,bins=bins)        
        bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        menStd  = numpy.sqrt(y)#or sigma?
        width=bins
        pyplot.bar(bincenters,y,width=width, color='r', yerr=menStd)
        # add a 'best fit' line?
#        y = mlab.normpdf( bins, mu, sigma)#should be the excepted distribution
#        l = pyplot.plot(bins, y, 'r--', linewidth=3)
        pyplot.savefig(filename)     
#        pylab.close()     # closes the current figure
        
    def plot(self,rdf,radii,filname):
        pylab.clf()
        matplotlib.rc('font', size=14)
        matplotlib.rc('figure', figsize=(5, 4))
#        pylab.clf()
        pylab.plot(radii, rdf, linewidth=3)
        pylab.xlabel(r"distance $r$ in $\AA$")
        pylab.ylabel(r"radial distribution function $g(r)$")
        pylab.savefig(filname)
#        pylab.close()     # closes the current figure

    def simpleplot(self,X,Y,filname,w=3):
        pylab.clf()
        pylab.plot(X, Y, linewidth=w)
        pylab.savefig(filname)
#        pylab.close()     # closes the current figure

    def grid_pack(self,bb,wrkDir,forceBuild=True, fill=0,seed=20,vTestid = 3,
                  vAnalysis = 0,fbox_bb=None):
        t1 = time()
#        if bbox is None :
#            box=self.helper.getCurrentSelection()[0]
#        else :
#            box = bbox[0]
#        bb=self.helper.getCornerPointCube(box)
        gridFileIn=None
        gridFileOut=None
#        self.env.grid.reset()
#        self.grid.reset()self.env.grid = None
#        if forceBuild :
#            gridFileOut=wrkDir+os.sep+"fill_grid"
#        else :
#            gridFileIn=wrkDir+os.sep+"fill_grid"
        print (gridFileIn,gridFileOut,forceBuild)
        self.env.buildGrid(boundingBox=bb,gridFileIn=gridFileIn,rebuild=forceBuild ,
                          gridFileOut=gridFileOut,previousFill=False)
    #    h.buildGrid(gridFileIn=gridFileIn, 
    #                  gridFileOut=gridFileOut)
        t2 = time()
        gridTime = t2-t1
        if fill :
            self.pack(seed=seed,vTestid = vTestid,vAnalysis = vAnalysis,fbox_bb=fbox_bb)
        print ('time to Build Grid', gridTime)
    #    afviewer.displayOrganellesPoints()
        #return
        #actine.updateFromBB(h.grid)
    
    def pack(self,seed=20,forceBuild=True,vTestid = 3,vAnalysis = 0,fbox_bb = None):
        t1 = time()
        print ("seed is ",seed, fbox_bb)
        self.env.fill5(seedNum=seed,verbose=4, vTestid = vTestid,vAnalysis = vAnalysis,fbox = fbox_bb)
        t2 = time()
        print('time to run Fill5', t2-t1)

    def calcDistanceMatrixFastEuclidean2(self,nDimPoints):
        nDimPoints = numpy.array(nDimPoints)
        n,m = nDimPoints.shape
        delta = numpy.zeros((n,n),'d')
        for d in xrange(m):
            data = nDimPoints[:,d]
            delta += (data - data[:,numpy.newaxis])**2
        return numpy.sqrt(delta)

    def flush(self):
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
                 
    def merge(self,d1, d2, merge=lambda x,y:y):
        result = dict(d1)
        for k,v in d2.iteritems():
            print k
            if k in result:
                result[k].extend(v)
            else:
                result[k] = v
        return result

    def plotNResult2D(self,n,bbox=[[0.0, 0, 0.], [1000.0, 1000.0, 1000.0]]):
        for i in range(n):
            f="results_seed_"+str(i)+".json"  
            self.plotOneResult2D(filename="f",bbox=bbox)

    def plotOneResult2D(self,data=None,filename=None,
                        bbox=[[0.0, 0, 0.], [1000.0, 1000.0, 1000.0]]):
        if data is None and filename is None :
            return
        elif data is None and filename is not None :
            with open(filename) as data_file:  
                data = json.load(data_file)
        width = 1000.0#should be the boundary here ?
        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        radius={}
        ingrrot={}
        ingrpos={}
        for recipe in data:
            for ingrname in data[recipe]:
                for k in range(len(data[recipe][ingrname]['results'])):
                    if ingrname not in ingrrot:
                        ingrrot[ingrname]=[]
                        ingrpos[ingrname]=[]  
                        radius[ingrname]=data[recipe][ingrname]['encapsulatingRadius']                 
                    ingrrot[ingrname].append(data[recipe][ingrname]['results'][k][1])
                    ingrpos[ingrname].append(data[recipe][ingrname]['results'][k][0])
        for ingr in ingrpos: 
            for i,p in enumerate(ingrpos[ingr]): 
                print (p,radius[ingr])
                ax.add_patch(Circle((p[0], p[1]),radius[ingr],
                            edgecolor="black", facecolor="red")) 
            ax.set_aspect(1.0)
            pyplot.axhline(y=bbox[0][1], color='k')
            pyplot.axhline(y=bbox[1][1], color='k')
            pyplot.axvline(x=bbox[0][0], color='k')
            pyplot.axvline(x=bbox[1][0], color='k')
            pyplot.axis([bbox[0][0], bbox[1][0],
                         bbox[0][1], bbox[1][1]])
            pyplot.savefig("plot"+ingr+".png")
            pylab.close()     # closes the current figure        
        return ingrpos
#        res=plotOneResult(None,filename="results_seed_8.json")
        
    def getHaltonUnique(self,n):
        from autopack.ldSequence import halton
        seeds_f = numpy.array(halton(int(n*1.5)))*int(n*1.5)
        seeds_int = numpy.array(numpy.round(seeds_f),'int')
        sorted_s,indices_u=numpy.unique(seeds_int,return_index=True)
        seeds_i = numpy.array(seeds_int[numpy.sort(indices_u)])[:n]
        return seeds_i
                
    def doloop(self,n,bbox,wrkDir,output,rdf=True, render=False, 
               plot = True,twod=True,fbox_bb=None,use_file = True,
               seeds_i=None):
        # doLoop automatically produces result files, images, and documents from the recipe while adjusting parameters
        # To run doLoop, 1) in your host's python console type:
        # execfile(pathothis recipe) # for example, on my computer:
        # " execfile("/Users/grahamold/Dev/autoFillSVN/autofillSVNversions/trunk/AutoFillClean/autoFillRecipeScripts/2DsphereFill/2DSpheres_setup_recipe.py")
        # 2) prepare your scene for the rendering->render output should be 640,480 but you can change the size in the script at then end.  Set up textures lights, and effects as you wish
        #    Results will appear in the result folder of your recipe path
        # where n is the number of loop, seed = i
        #analyse.doloop(n) 
        if seeds_i is None :
            seeds_i = self.getHaltonUnique(n)
        numpy.savetxt(output+os.sep+"seeds", seeds_i, delimiter=",")
        angle_file= output+os.sep+"angle"
        position_file=output+os.sep+"pos"
        distance_file=output+os.sep+"dist"
        position_files=output+os.sep+"pos.json"
        distance_files=output+os.sep+"dist.json"
        occurences_file=output+os.sep+"occurence"
        rangeseed=range(n)
        distances={}
        ingrpositions={}
        anglesingr={}
        occurences={}
        total_positions=[]
        total_distances=[]
        total_angles=[]
        self.bbox = bbox
        angles=None
        rebuild = True
#        seeds_i=[]
        for si in range(n):           
#            if i > 0 : rebuild = False #bu need to reset ...
            basename = output+os.sep+"results_seed_"+str(si)
#            self.env.saveResult = False
            resultfilename = self.env.resultfile = basename  
            #Clear
            if self.afviewer :
                self.afviewer.clearFill("Test_Spheres2D")
            else :
                self.env.reset()
            #no need to rebuild the grid ?
            self.env.saveResult = True
            self.env.resultfile = basename
            se=seeds_i[si]#int(time())
#            seeds_i.append(se)
            self.grid_pack(bbox,output,seed=se, fill=1,vTestid = si,vAnalysis = 1,
                           forceBuild=rebuild,fbox_bb=fbox_bb)#build grid and fill
            #store the result 
#            self.env.store_asJson(resultfilename=basename)
            self.center = self.env.grid.getCenter()
            if render :
                #render/save scene if hosted otherwise nothing
                self.helper.render(basename+".jpg",640,480)
                self.helper.write(basename+".c4d",[])
            if rdf :
                if plot and twod:
                    width = 1000.0#should be the boundary here ?
                    fig = pyplot.figure()
                    ax = fig.add_subplot(111)
                center = self.env.grid.getCenter()#[500,500,0.5]#center of the grid
                r = self.env.exteriorRecipe
                d={}
                if r :
#                    print ("DONERUNXXXCYTO!!!!!")
                    for ingr in r.ingredients:
                        if ingr.name not in distances :
                            distances[ingr.name]=[]
                            ingrpositions[ingr.name]=[]
                            anglesingr[ingr.name]=[]
                            occurences[ingr.name]=[]
                        if ingr.packingMode=='gradient' and self.env.use_gradient:
                            self.center = center = self.env.gradients[ingr.gradient].direction
                            #also mesure the angle pos>center pcpalVector
                            ingrpos,d,angles=self.getDistanceAngle(ingr, center)
                            if use_file :
                                f_handle = file(angle_file, 'a')
                                numpy.savetxt(f_handle, angles, delimiter=",")
                                f_handle.close()
                            else :
                                anglesingr[ingr.name].extend(angles)
                                total_angles.extend(angles)
                        else :
                            ingrpos,d=self.getDistance(ingr.name, center)
                        occurences[ingr.name].append(len(ingrpos))
                        if use_file :
                            f_handle = file(position_file, 'a')
                            numpy.savetxt(f_handle, ingrpos, delimiter=",")
                            f_handle.close()
                            f_handle = file(distance_file, 'a')
                            numpy.savetxt(f_handle, d, delimiter=",")
                            f_handle.close() 
                            distances[ingr.name]=d
                            ingrpositions[ingr.name]=ingrpos.tolist()
                        else :
                            distances[ingr.name].extend(d)
                            ingrpositions[ingr.name].extend(ingrpos)
                            total_positions.extend(ingrpos)
                            total_distances.extend(d)
                           
                        #print plot,twod
                        if plot and twod:
                            for i,p in enumerate(ingrpos): 
                                ax.add_patch(Circle((p[0], p[1]), ingr.encapsulatingRadius,
                                                edgecolor="black", facecolor=ingr.color))
#                                 Plot "image" particles to verify that periodic boundary conditions are working
                                r=ingr.encapsulatingRadius
                                if autopack.testPeriodicity :
                                    if p[0] < r:
                                        ax.add_patch(Circle((p[0] + width,p[1]), r, facecolor=ingr.color))
                                    elif p[0] > (width-r):
                                        ax.add_patch(Circle((p[0] - width,p[1]), r, facecolor=ingr.color))
                                    if p[1] < r:
                                        ax.add_patch(Circle((p[0],p[1]+ width), r, facecolor=ingr.color))
                                    elif p[1] > (width-r):
                                        ax.add_patch(Circle((p[0],p[1] - width), r, facecolor=ingr.color))
                                if i == 0:#len(ingrpos)-1:
                                    continue
                                if ingr.Type== "Grow":
                                    pyplot.plot([ingrpos[-i][0], ingrpos[-i-1][0]], 
                                            [ingrpos[-i][1], ingrpos[-i-1][1]], 'k-', lw=2)  
                                    #plot the sphere
                                    if ingr.use_rbsphere :
                                        r,pts = ingr.getInterpolatedSphere(ingrpos[-i-1],ingrpos[-i])
                                        for pt in pts :
                                            ax.add_patch(Circle((pt[0], pt[1]), ingr.minRadius,
                                                edgecolor="black", facecolor=ingr.color))
#                print ("DONERUNXXXbefore!!!!!")
                for o in self.env.compartments:
#                    print ("DONERUNXXXComp!!!!!")
                    rs = o.surfaceRecipe
                    if rs :
                        for ingr in rs.ingredients:
                            if ingr.name not in distances :
                                distances[ingr.name]=[]
                                ingrpositions[ingr.name]=[]
                                occurences[ingr.name]=[]
                            if ingr.packingMode=='gradient' and self.env.use_gradient :
                                center = self.env.gradients[ingr.gradient].direction
                            ingrpos,d=self.getDistance(ingr.name, center)
                            occurences[ingr.name].append(len(ingrpos))
                            if use_file :
                                f_handle = file(position_file, 'a')
                                numpy.savetxt(f_handle, ingrpos, delimiter=",")
                                f_handle.close()
                                f_handle = file(distance_file, 'a')
                                numpy.savetxt(f_handle, d, delimiter=",")
                                f_handle.close()
                                distances[ingr.name]=d.tolist()
                                ingrpositions[ingr.name]=ingrpos.tolist()
                            else :
                                distances[ingr.name].extend(d)
                                ingrpositions[ingr.name].extend(ingrpos)
                                total_positions.extend(ingrpos)
                                total_distances.extend(d)
                            if plot and twod:
                                for p in ingrpos: 
                                    ax.add_patch(Circle((p[0], p[1]), ingr.encapsulatingRadius,
                                                    edgecolor="black", facecolor=ingr.color))
#                    print ("DONERUNXXXYYY!!!!!")
                    ri = o.innerRecipe
                    if ri :
                        for ingr in ri.ingredients:
                            if ingr.name not in distances :
                                distances[ingr.name]=[]
                                ingrpositions[ingr.name]=[]
                                occurences[ingr.name]=[]
                            if ingr.packingMode=='gradient' and self.env.use_gradient:
                                center = self.env.gradients[ingr.gradient].direction
                            ingrpos,d=self.getDistance(ingr.name, center)
                            occurences[ingr.name].append(len(ingrpos))
                            if use_file :
                                f_handle = file(position_file, 'a')
                                numpy.savetxt(f_handle, ingrpos, delimiter=",")
                                f_handle.close()
                                f_handle = file(distance_file, 'a')
                                numpy.savetxt(f_handle, d, delimiter=",")
                                f_handle.close()
                                distances[ingr.name]=d.tolist()
                                ingrpositions[ingr.name]=ingrpos.tolist()
                            else :
                                distances[ingr.name].extend(d)
                                ingrpositions[ingr.name].extend(ingrpos)
                                total_positions.extend(ingrpos)
                                total_distances.extend(d)
                            if plot and twod:
                                for p in ingrpos: 
                                    ax.add_patch(Circle((p[0], p[1]), ingr.encapsulatingRadius,
                                        edgecolor="black", facecolor=ingr.color))
                #write
                if use_file :
                    self.writeJSON(output+os.sep+"_posIngr_"+str(si)+".json",ingrpositions)
                    self.writeJSON(output+os.sep+"_dIngr_"+str(si)+".json",distances)
                    self.writeJSON(output+os.sep+"_angleIngr_"+str(si)+".json",anglesingr)
                #print ("############")                 
                #print (output+os.sep+"_dIngr_"+str(si)+".json",distances)
#                print ("DONERUNXXX!!!!!")
                if plot and twod:
                    ax.set_aspect(1.0)
                    pyplot.axhline(y=bbox[0][1], color='k')
                    pyplot.axhline(y=bbox[1][1], color='k')
                    pyplot.axvline(x=bbox[0][0], color='k')
                    pyplot.axvline(x=bbox[1][0], color='k')
                    pyplot.axis([bbox[0][0], bbox[1][0],
                                 bbox[0][1], bbox[1][1]])
                    pyplot.savefig(basename+".png")
                    pylab.close()     # closes the current figure
#            return                
#            print ("DONERUN!!!!!")
#            self.flush()        
        #plot(x)
#        print ("DONE1!!!!")
        numpy.savetxt(output+os.sep+"seeds", seeds_i, delimiter=",")
        if use_file :
            total_positions = numpy.genfromtxt(position_file, delimiter=',')
            try :
                total_angles = numpy.genfromtxt(angle_file, delimiter=',')
            except :
                total_angles=[]
            #gatherall result
            ingrpositions={}
            distances={}
            anglesingr={}
            for i in rangeseed:
                dict1= self.loadJSON(output+os.sep+"_posIngr_"+str(i)+".json")
                ingrpositions=dict(self.merge(ingrpositions,dict1))
                dict1= self.loadJSON(output+os.sep+"_dIngr_"+str(i)+".json")
                distances=dict(self.merge(distances,dict1))
                dict1=self.loadJSON(output+os.sep+"_angleIngr_"+str(i)+".json")
                anglesingr=dict(self.merge(anglesingr,dict1))
#        print ("DONE2!!!!")
        self.writeJSON(occurences_file,occurences)
        self.env.ingrpositions=ingrpositions
        self.env.distances = distances
        self.env.basename = basename
        self.env.occurences = occurences
        self.env.angles =  total_angles
#        if rdf :
#            if twod : self.env.loopThroughIngr(self.rdf_2d)
#            else : self.env.loopThroughIngr(self.rdf_3d)
#            do the X Y histo averge !        
        self.env.loopThroughIngr(self.axis_distribution)        
        self.env.loopThroughIngr(self.occurence_distribution)        
        self.axis_distribution_total(total_positions)
#        self.env.loopThroughIngr(self.correlation)
#        print ("DONE3!!!!")
        #plot the angle
        if len(total_angles) : 
            self.histo(total_angles[1],output+os.sep+"_anglesX.png",bins=12,size=max(total_angles[2]))
            self.histo(total_angles[2],output+os.sep+"_anglesY.png",bins=12,size=max(total_angles[2]))
            self.histo(total_angles[3],output+os.sep+"_anglesZ.png",bins=12,size=max(total_angles[2]))
#        print ("DONE!!!!")
        return distances
#from bhtree import bhtreelib
#bht = bhtreelib.BHtree( verts, None, 10)
#closest = bht.closestPointsArray(tuple(grdPos), diag, returnNullIfFail)
        
