# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 22:34:09 2014

@author: ludo
"""
import json
import numpy
import math
import upy
from math import sqrt,asin,cos,sin,pi
from scipy import spatial
from numpy import matrix
from autopack.transformation import euler_from_matrix,matrixToEuler,euler_matrix,superimposition_matrix,rotation_matrix,affine_matrix_from_points
from random import randrange,random,seed
from upy.colors import getRamp
from time import time

SEED = 12
helper = upy.getHelperClass()()
DEBUG=False
global HEXAGON_FREE#[iedge_w,]
global HEXAGON
global ALL_HEXAGON
global ALL_HEXAGON_POS

def vdistance(c0,c1):
    """get the distance between two points c0 and c1"""
    d = numpy.array(c1) - numpy.array(c0)
    s = numpy.sum(d*d)
    return math.sqrt(s)
def vdiff(p1, p2):
    # returns p1 - p2
    x1,y1,z1 = p1
    x2,y2,z2 = p2
    return (x1-x2, y1-y2, z1-z2)

def vcross(v1,v2):
    x1,y1,z1 = v1
    x2,y2,z2 = v2
    return (y1*z2-y2*z1, z1*x2-z2*x1, x1*y2-x2*y1)

def dot(v1,v2):
    x1,y1,z1 = v1
    x2,y2,z2 = v2
    return ( x1 * x2 ) + ( y1 * y2 ) +  ( z1 * z2 )
    
from math import sqrt
def vnorm(v1):
    x1,y1,z1 = v1
    n1 = 1./sqrt(x1*x1 + y1*y1 + z1*z1)
    return (x1*n1, y1*n1, z1*n1)


def vlen(v1):
    x1,y1,z1 = v1
    return sqrt(x1*x1 + y1*y1 + z1*z1)

def getFaceNormalsArea( vertices, faces):
    """compute the face normal of the compartment mesh"""
    normals = []
    vnormals = numpy.array(vertices[:])
    areas = [] #added by Graham
    face = [[0,0,0],[0,0,0],[0,0,0]]
    v = [[0,0,0],[0,0,0],[0,0,0]]        
    for f in faces:
        for i in range(3) :
            face [i] = vertices[f[i]]
        for i in range(3) :
            v[0][i] = face[1][i]-face[0][i]
            v[1][i] = face[2][i]-face[0][i]                
        normal = vcross(v[0],v[1])
        n = vlen(normal)
        if n == 0. :
            n1=1.
        else :
            n1 = 1./n
        normals.append( (normal[0]*n1, normal[1]*n1, normal[2]*n1) )
#        The area of a triangle is equal to half the magnitude of the cross product of two of its edges
        for i in range(3) :
            vnormals[f[i]] = [normal[0]*n1, normal[1]*n1, normal[2]*n1]
        areas.append(0.5*vlen(normal)) #added by Graham
    return vnormals,normals, areas
    
def getVertexNormals(vertices, faces):
    vnormals = numpy.array(vertices[:])
    face = [[0,0,0],[0,0,0],[0,0,0]]
    v = [[0,0,0],[0,0,0],[0,0,0]]        
    for f in faces:
        for i in range(3) :
            face [i] = vertices[f[i]]
        for i in range(3) :
            v[0][i] = face[1][i]-face[0][i]
            v[1][i] = face[2][i]-face[0][i]                
        normal = vcross(v[0],v[1])
        n = vlen(normal)
        if n == 0. :
            n1=1.
        else :
            n1 = 1./n
        for i in range(3) :
            vnormals[f[i]] = [normal[0]*n1, normal[1]*n1, normal[2]*n1]
    return vnormals #areas added by Graham


class Hexa:
    def __init__(self,name,center,radius,rotation=[]):
        self.name=name
        self.free_pos=numpy.ones(6)
        self.children=[None,None,None,None,None,None]
        self.nvisit=numpy.zeros(6)
        self.pcpalAxis,self.coordinates,self.neighboors_pos=self.getCoordHexagone(radius)
        self.edge_pos = numpy.array(self.neighboors_pos)/2.0
        self.edge_tree=None
        self.tr_c = self.coordinates[:]
        self.tr_n = self.neighboors_pos[:]
        self.nr=[]
        self.nb=[[1,3],[2,4],[2,0],[3,1],[5,1],[5,3]]
        self.inst=None
        self.edge_nb=[[[0,4],[1,3]],
                      [[1,5],[2,4]],
                      [[2,0],[3,5]],
                      [[3,1],[4,0]],
                      [[4,2],[5,1]],
                      [[5,3],[0,2]]]
        self.cb={0:3,1:4,2:5,3:0,4:1,5:2}
#        self.cbi={3:0,4:1,5:2,0:3,1:4,2:5}
        self.colors=getRamp([[0,0.6,0],[1,1,0]],size=7)
        self.threshold = 500.0
        self.pos=center
        self.rejected=False
        self.mat = numpy.identity(4)
        if not len(rotation) :
            rotation = numpy.identity(4)
        if len(rotation):
            self.mat = numpy.array(rotation)#.transpose()
            self.mat[:3, 3] = center            
            self.tr_c=helper.ApplyMatrix(self.coordinates,self.mat)
            self.tr_n=helper.ApplyMatrix(self.neighboors_pos,self.mat)
            self.tr_edge=helper.ApplyMatrix(self.edge_pos,self.mat)
            self.edge_tree=spatial.cKDTree(self.tr_edge, leafsize=10) 
            self.nr=helper.ApplyMatrix(self.neighboors_pos,rotation)
#        self.displaySpheres()
        self.displayObject()
        
    def displaySpheres(self,):
        root = helper.getObject("all_spheres")
        if root is None:
            root = helper.newEmpty("all_spheres")
        parent = helper.newEmpty(self.name+"_iSpheres",parent=root)
        s,ms = helper.Sphere(self.name+"_base_sphere",radius=5.0,res=12,parent=parent)
        hsphs=helper.instancesSphere(self.name+"_spheresH",self.tr_c,[1,]*6,s,[[0,1,0]],None,parent=parent)
        sphs=helper.instancesSphere(self.name+"_sphereHN",self.tr_n,[1,]*6,s,[[0,0,1]],None,parent=parent)
        return hsphs,sphs
        
    def displayObject(self,mobject=None):
        if mobject is None:
            mobject=helper.getObject("hex")
        parent = helper.getObject("all_hex")
        if parent is None:
            parent = helper.newEmpty("all_hex")
        self.inst=helper.newInstance(self.name+"_obj",mobject,matrice=self.mat,parent=parent)
        return self.inst
                    
    def getCoordHexagone(self,radius,center=[0,0,0]):
        D=2*radius
        dY=D/4.0
        dX=math.sqrt(3)*dY
        self.innerR=dX
        pcpalAxis=[0,0,1]
        hexc=coordinates=numpy.array([[0,2.0*dY,0.0],#up
                    [dX,dY,0],#right-up
                    [dX,-dY,0],#right-down
                    [0,-2.0*dY,0],#down
                    [-dX,-dY,0],#left-down
                    [-dX,dY,0]#left-up
                     ])
        #neighboors_pos
        #hexc = hexc + numpy.array(center)
        neighboors_pos=numpy.zeros([6,3])
        for i in range(5):
            neighboors_pos[i]=hexc[i]+hexc[i+1]
        neighboors_pos[5] = hexc[5]+hexc[0]  
        return pcpalAxis,coordinates,neighboors_pos

    def distanceEdge(self, edge_id):
        #{0:3,1:4,2:5,3:0,4:1,5:2}
        #neighbor = self.children[edge_id]
#        print "edge distances ",D1,D2        
        edge1=self.tr_c[self.edge_nb[edge_id][0][0]]-self.tr_c[self.edge_nb[edge_id][1][0]]
        edge2=self.children[edge_id].tr_c[self.edge_nb[edge_id][0][1]]-self.children[edge_id].tr_c[self.edge_nb[edge_id][1][1]]
        a = angle_between_vectors([edge1], [edge2], axis=1)[0]
        D1 = vdistance(self.tr_c[self.edge_nb[edge_id][0][0]],self.children[edge_id].tr_c[self.edge_nb[edge_id][0][1]])
        D2 = vdistance(self.tr_c[self.edge_nb[edge_id][1][0]],self.children[edge_id].tr_c[self.edge_nb[edge_id][1][1]])
#        print (edge_id,"edge distances and angle ",D1,D2,a,math.degrees(a))      
#        rotMat = numpy.array( rotVectToVect(edge2, edge1 ), 'f') 
        return max(D1,D2),math.degrees(a)

    def updateDistance(self,cutoff=None):
        distances=[]
        angles=[]
        for i in range(6):
            if self.free_pos[i] == 1 : #free pos
                print i,"free"
                continue
            if self.children[i] == None : #no children
                print i,"null"
                continue
            d,a = self.distanceEdge(i)
            distances.append(d)
            angles.append(a)
        D=max(distances)
        A=max(angles)
#        print "update",angles,A
        if cutoff == None :
            cutoff = self.threshold
        if A > cutoff :
            self.rejected = True
#            helper.toggleDisplay(self.inst,False)
            helper.deleteObject(self.inst)
            self.inst=None
            
    def updateOneChildrenFreePos(self,i):
        if self.free_pos[i] == 1 : #free pos
            print i,"free"
            return 0
        if self.children[i] == None : #no children
            print i,"null"
            return 0
        self.children[i].free_pos[self.cb[i]]=0#start point
        self.children[i].children[self.cb[i]]=self
#            print self.children[i].free_pos
        indicedown = i - 1 
        if ( indicedown == -1 ) : indicedown = 5
        indicedup = i + 1
        if ( indicedup == 6 ) : indicedup = 0
#            print indicedown,self.children[indicedown].name
#            print indicedup, self.children[indicedup].name
#            print self.cb[indicedown],self.children[self.cb[indicedown]].name
#            print self.cb[indicedup], self.children[self.cb[indicedup]].name
#            
        if self.free_pos[indicedown] == 0:
            if self.children[i].free_pos[self.cb[indicedown]] == 1:
                self.children[i].free_pos[self.cb[indicedown]]=0
                self.children[i].children[self.cb[indicedown]] = self.children[indicedup]
        if self.free_pos[indicedup] == 0:
            if self.children[i].free_pos[self.cb[indicedup]] == 1:
                self.children[i].free_pos[self.cb[indicedup]]=0
                self.children[i].children[self.cb[indicedup]] = self.children[indicedown]
        return self.distanceEdge(i)  
        
    def updateChildrenFreePos(self,idc=None):
        distance = []
        angles=[]
        ids=0
        if idc == None :
            for i in range(6):
                if self.free_pos[i] == 1 : #free pos
                    #print i,"free"
                    continue
                if self.children[i] == None : #no children
                    #print i,"null"
                    continue
                d,a = self.updateOneChildrenFreePos(i)
                distance.append(d)
                angles.append(a)
            A=0#max(angles)
            ids=0#angles.index(A)
        else :
            d,A=self.updateOneChildrenFreePos(idc)
            ids = idc
#        print "A is ",A,self.threshold
        if A > self.threshold :
            if not self.children[ids].rejected:
#                self.rejected = True
#                helper.toggleDisplay(self.inst,False)
#                helper.deleteObject(self.inst)
#                self.inst=None
                print "reject ", self.children[ids].name
                self.free_pos[ids] = 1
                self.children[ids].rejected = True
                helper.deleteObject(self.children[ids].inst)
                self.children[ids].inst = None
#                self.children[ids].children[self.cb[ids]] = None
#                self.children[ids].free_pos[self.cb[ids]] = 1
                #delete it ?                
                try :
                    idh = HEXAGON.index(self)
                except :
                    HEXAGON.append(self)
                try :
                    idah = ALL_HEXAGON.index(self.children[ids])
                    ALL_HEXAGON.pop(idah)
                    ALL_HEXAGON_POS.pop(idah)
                except :
                    print self.children[ids].name+"not in big list"
                self.children[ids].rejectItSelf()
                self.children[ids] = None 
                   
        if not len(numpy.nonzero(self.free_pos)[0]):
            try :
                idh = HEXAGON.index(self)
                HEXAGON.pop(idh)
#                idah = ALL_HEXAGON.index(self.children[ids])
                ALL_HEXAGON.pop(idh)
                ALL_HEXAGON_POS.pop(idh)
            except :
                print self.name+"not in big list"
#            print self.children[i].name,self.children[i].free_pos
#            for c in self.children[i].children : 
#                if c is not None : print c.name
#                else : print "None"
#        self.updateColor()

    def updateFreePos(self,tiles):
        #grab surrounding object hexa, fin eac edge they occupied
        T=tiles.all_pos_tree#spatial.cKDTree(tiles.ALL_HEXAGON_POS, leafsize=10)
        liste_surrounding = T.query_ball_point(self.pos,125.0)
        for ih in liste_surrounding:
            if tiles.ALL_HEXAGON[ih].name == self.name :
                continue
            d,i=self.edge_tree.query(tiles.ALL_HEXAGON[ih].pos)
            self.free_pos[i]=0
            self.children[i]=tiles.ALL_HEXAGON[ih]
            tiles.ALL_HEXAGON[ih].free_pos[self.cb[i]]=0
            tiles.ALL_HEXAGON[ih].children[self.cb[i]]=self   
        
        self.updateColor()
        
    def updateColor(self):
        helper.changeObjColorMat(self.inst,self.colors[len(numpy.nonzero(self.free_pos)[0])])
        for i in range(6):
            if self.free_pos[i] == 1 : #free pos
                continue
            if self.children[i] == None : #no children
                continue
#            print i,numpy.nonzero(self.children[i].free_pos)
            try :
                helper.changeObjColorMat(self.children[i].inst,self.colors[len(numpy.nonzero(self.children[i].free_pos)[0])])
            except :
                print i,numpy.nonzero(self.children[i].free_pos),self.name

    def rejectItSelf(self):
        for i in range(6):
            if self.free_pos[i] == 1 : #free pos
                continue
            if self.children[i] == None : #no children
                continue
            self.children[i].free_pos[self.cb[i]]=1#start point
            self.children[i].children[self.cb[i]]=None
        
    def nextFree(self):
        return numpy.nonzero(self.free_pos)[0][0]

class Tilling:
    def __init__(self,v,f,vn,inputMesh,hradius,init_seed=12,inputConstraint=None):
        self.HEXAGON_FREE=[]#[iedge_w,]
        self.HEXAGON=[]
        self.ALL_HEXAGON=[]
        self.ALL_HEXAGON_POS=[]
        self.all_pos_tree=None
        self.area = 0
        self.surface_covered = 0
        self.v=[]
        self.f=[]
        self.vn=[]
        self.vtree=None
        self.rc=None
        self.inputMesh = inputMesh
        self.constraintMesh = inputConstraint
        if not len (v) and inputMesh!=None:
            f,v,vn = helper.DecomposeMesh(input_mesh,edit=False,copy=False,tri=True,transform=True)
            vn,fn,area=getFaceNormalsArea(v, f)
            self.area = sum(area)
            self.v=v
            self.f=f
            self.vn=vn
            self.vtree= spatial.cKDTree(v, leafsize=10)
            self.rc = utils.GeRayCollider()
            self.rc.Init(input_mesh)
        elif inputConstraint != None :
            f,v,vn = helper.DecomposeMesh(inputConstraint,edit=False,copy=False,tri=True,transform=True)
            vn,fn,area=getFaceNormalsArea(v, f)
#            self.area = sum(area)
            self.v=v
            self.f=f
            self.vn=vn
            self.vtree= spatial.cKDTree(v, leafsize=10)
            self.rc = utils.GeRayCollider()
            self.rc.Init(inputConstraint)            
        self.hexamer_radius=hradius#55
        self.alpha = hradius# (hradius/2.)*math.sqrt(3)#or radius 47
        self.sphs=[]
        self.distance_threshold = [self.hexamer_radius*1.5,self.hexamer_radius*1.85]
        self.display = True
        self.hexa_area = ((3.0*math.sqrt(3.0))/2.0)*(self.alpha*self.alpha)
        D=2*hradius
        dY=D/4.0
        dX=math.sqrt(3)*dY
        self.innerR=dX
        numpy.random.seed(init_seed)#for gradient
        seed(init_seed)

    def getAllPos(self,):
        return [h.pos for h in self.HEXAGON]

    def getFreePos(self,weight=False):
        T=self.all_pos_tree#spatial.cKDTree(self.ALL_HEXAGON_POS, leafsize=10)
        ids=[]
        colors=[]
        pts=[]
        w=[]
        qw = [0.0,0.001,0.05,0.1,0.25,1,1]
        for i in range(len(self.HEXAGON)):
            #self.HEXAGON[i].updateFreePos()
            for j in range(6):
                if self.HEXAGON[i].free_pos[j] == 1 :#free
                    #find how many neighboors are present for this position
                    pts.append(self.HEXAGON[i].tr_n[j])                    
                    ids.append([i,j])
                    if weight :
                        results = T.query_ball_point(self.HEXAGON[i].tr_n[j],115.0)
                        R= len(results)
                        if R > 6:
                            R=6
                        w.append(qw[R] )
                    else :
                        R=0
                    colors.append(self.HEXAGON[i].colors[R])
        return ids,pts,colors,w

    def displayWeights(self,points,colors):
        parent = helper.getObject("W_iSpheres")
        if parent is None :
            parent = helper.newEmpty("W_iSpheres")
        s = helper.getObject("W_base_sphere")
        if s is None :
            s,ms = helper.Sphere("W_base_sphere",radius=5.0,res=12,parent=parent)
        self.sphs=helper.updateInstancesSphere("W_spheresH",self.sphs,points,[1.0,]*len(points),s,
                            colors,None,parent=parent,delete=False)
        
    def weightFreePos(self,startpos=None):
        T=self.all_pos_tree#spatial.cKDTree(self.ALL_HEXAGON_POS, leafsize=10)
        self.HEXAGON_FREE=numpy.zeros((len(self.HEXAGON),6))
        self.HEXAGON_FREE=self.HEXAGON_FREE.reshape(len(self.HEXAGON)*6)
        #qw=numpy.arange(0.0,1.0+1./6.,1.0/6.0)
        maxd=300.0
        qw = [0.0,0.001,0.05,0.1,0.25,1,1]
        pts=[]
        colors=[]
        for i in range(len(self.HEXAGON)):
            self.HEXAGON[i].updateFreePos(self)
    #        hid=ALL_HEXAGON.index(HEXAGON[i])
            for j in range(6):
                if self.HEXAGON[i].free_pos[j] == 1 :#free
                    #find how many neighboors are present for this position
                    pts.append(self.HEXAGON[i].tr_n[j])
                    results = T.query_ball_point(self.HEXAGON[i].tr_n[j],115.0)
    #                print "r ",j,results
    #                if hid in results :
    #                    print "pop ",hid,results.index(hid)
    #                    results.pop(results.index(hid))
                    R= len(results)
                    if R > 6:
                        R=6
                    w=qw[R]
    #                print "t ",i,j,results,R,HEXAGON[i].name
    #                for ii in results :
    #                    print "ii",ii,HEXAGON[i].name,j,ALL_HEXAGON[ii].name
    #                d=vdistance(startpos,HEXAGON[i].tr_n[j])
    #                if d > maxd  : d = maxd
    #                p=1.0-d/maxd
                    self.HEXAGON_FREE[6 * i + j]=w#(w+p)/2.0
                    colors.append(self.HEXAGON[i].colors[R])
    #    displayWeights(points,colors)
        return pts,colors
        
    def getRndWeighted(self,w=None,startpos=None):
        """
        From http://glowingpython.blogspot.com/2012/09/weighted-random-choice.html
        Weighted random selection
        returns n_picks random indexes.
        the chance to pick the index i 
        is give by the weight weights[i].
        """
        weight = w
        if weight is None :
            weight,pts,colors = self.weightFreePos(startpos=startpos)#numpy.take(self.weight,listPts)
        t = numpy.cumsum(weight)
        s = numpy.sum(weight)
        i = numpy.searchsorted(t,numpy.random.rand(1)*s)[0]
        return i   

    def alignToAdjacentPos(self,h):
        #align the given hexa to his closest adjacent hexa ?
        #only if angle > threshold
        pass

    def constrainPosRot(self,next_pos,h=None,edge=0):
        maxd=600.0
        maxangle=25.0#degree
        mat=numpy.identity(4)
        #closest surface point. depending the distance rotate a little ?
        distance,indice = self.vtree.query(next_pos)
        print ("distance to closest surface point",distance)#,self.distance_threshold,self.alpha)
        #X goes all the way to the surface and stop
        #Y should curve
        #check dx,dy 
        reject=False
        f=abs(next_pos[0])/self.innerR
        print "dX pos factor is",f,math.floor(f),math.ceil(f)
#        print "dX ",next_pos[0],next_pos[0]-self.v[indice][0]
#        print "dY ",next_pos[1],next_pos[1]-self.v[indice][1]
#        if abs(next_pos[0]) > abs(next_pos[1]) :#X direction
            #do nothing until touch surface then stop growing that direction
        if distance < 50.0 :
            reject = True
        elif abs(next_pos[0]) > 100.0 :#distance < 300.0 : #Y direction
            #lut depedngin on how may radius we have
            #
            #increase the rotation angle accorind the surface ?
            #calculate angle depedning on the distance
            #angle = 10.0*(next_pos[0]/abs(next_pos[0]))#((abs(next_pos[0])/maxd)*maxangle)*(next_pos[1]/abs(next_pos[1]))#*-1.0#to get the sign            
            angle = 5.0*(next_pos[0]/abs(next_pos[0]))#(((distance/maxd))*maxangle)*(next_pos[0]/abs(next_pos[0]))#*-1.0#to get the sign
            #should use LUT
            print ("angle next",angle)
            #rotMat1 = euler_matrix(0.0,0.0,math.radians(30.0)).transpose()
            rotMat2 = euler_matrix(0.0,math.radians(angle),0.0).transpose()
            m=numpy.identity(4)
            m[:3,:3]=h.mat[:3,:3]
            rotMat = matrix(rotMat2)*matrix(m)
            mat = numpy.array(rotMat)
            #Hypothenus = 47.63
            #opposite = Hypothenus*math.sin(angle)
            #final vector is up (0,0,1)*R
#            m=numpy.identity(4)
#            m[:3,:3]=h.mat[:3,:3]
#            tr_up=helper.ApplyMatrix([numpy.array([0.0,0.0,1.0]),],m)[0]
            up=numpy.array([0.0,0.0,1.0])*47.63*math.sin(math.radians(angle))
            new_pos = helper.ApplyMatrix([up,],mat)[0]*(next_pos[0]/abs(next_pos[0]))
            next_pos = next_pos+new_pos
            mat[:3, 3] = next_pos
#            p = (h.neighboors_pos[edge]-h.pos)/2.0
##            mat[:3, 3] = next_pos#helper.ApplyMatrix([next_pos,],mat)[0]
#            #the position should be vector from neighbors rotated.
#            #transfor the neighboors_pos
#            ptToRot = (h.neighboors_pos[edge]-h.pos)/2.0#tr_edge[edge]
#            print "rotate ptToRot",edge,ptToRot,h.neighboors_pos[edge],h.tr_edge[edge]
#            n_pos = h.tr_edge[edge]+helper.ApplyMatrix([ptToRot,],mat)[0]
#            print "rotate ptToRot to",n_pos
#            #trasnform this one to actual position
#            tr_n=helper.ApplyMatrix([n_pos,],h.mat)[0]
#            print "rotate ptToRot applied mat ",tr_n
        else :
            next_pos[2]=0.0
            print ("do nothing",next_pos)
            pass
        return indice,next_pos,reject,mat
        
    def computePosRot(self,next_pos):
        distance,indice = self.vtree.query(next_pos)
        if self.rc.Intersect(helper.FromVec(numpy.array(next_pos)*100.0), helper.FromVec(next_pos)*-1, 100000):
            ray_result =self. rc.GetNearestIntersection()
            n = ray_result["f_normal"].GetNormalized()
            v2 = helper.ToVec(n)# vn[indice]
        else :
            v2 = self.vn[indice]#next_pos#vn[indice]
        vnax=numpy.zeros((4,4))
        vnax[0][:3] = vnorm(numpy.cross(v2,[0,1,0]))
        vnax[1][:3] = vnorm(numpy.cross(v2,vnax[0][:3]))
        vnax[2][:3] = vnorm(v2)#pcpal axis is [0,0,1]
        euler = euler_from_matrix(vnax)
        rotMat = euler_matrix(euler[0],euler[1],0.0).transpose()
        mat = numpy.array(rotMat)
        mat[:3, 3] = next_pos
        if self.rc.Intersect(helper.FromVec(numpy.array(next_pos)*10.0), helper.FromVec(next_pos)*-1, 100000):
            ray_result = self.rc.GetNearestIntersection()
            newposition = helper.ToVec(ray_result["hitpos"])
        else :
            newposition=self.v[indice]#(next_pos+tr)#-pm.transpose()[3,:3]#+tr#pos+p
        mat = numpy.array(rotMat)#matrix(R)*matrix(rotMat))#
        mat[:3, 3] = newposition#+tr
        return indice,newposition,rotMat,mat


    def clearOccupied(self):
        i=len(self.HEXAGON)-1
        #T=self.all_pos_tree#spatial.cKDTree(self.ALL_HEXAGON_POS, leafsize=10)
        while i > 0:    
#            results = T.query_ball_point(self.HEXAGON[i].pos,115.0)
#            print i,self.HEXAGON[i].name,len(numpy.nonzero(self.HEXAGON[i].free_pos)[0])
            if not len(numpy.nonzero(self.HEXAGON[i].free_pos)[0]):
                try :
#                    print self.HEXAGON[i].name+" pop"
                    self.HEXAGON.pop(i)
    #                ALL_HEXAGON.pop(i)
    #                ALL_HEXAGON_POS.pop(i)
                except :
                    print self.HEXAGON[i].name+"not in big list"
    #        else :
    #           R= len(results)
    #           if R >= 6:
    #                HEXAGON[i].free_pos=numpy.zeros(6)
    #                try :
    #                    HEXAGON.pop(i)
    #    #                ALL_HEXAGON.pop(i)
    #    #                ALL_HEXAGON_POS.pop(i)
    #                except :
    #                    print HEXAGON[i].name+"not in big list"
                
            i=i-1

    def randomStart(self,point):
        start=Hexa("hexa_"+str(len(self.HEXAGON)),point,self.hexamer_radius)
        self.HEXAGON.append(start)
        self.ALL_HEXAGON.append(start)
        self.ALL_HEXAGON_POS.append(start.pos)
        self.all_pos_tree = spatial.cKDTree(self.ALL_HEXAGON_POS, leafsize=10)
        start.updateFreePos(self)
        
    def firstRandomStart(self):
        p = helper.advance_randpoint_onsphere(650,marge=math.radians(90.0),vector=[0,1,0])
        indice,newposition,rotMat,mat=self.computePosRot(p)
        start=Hexa("hexa_"+str(len(self.HEXAGON)),newposition,hexamer_radius,mat)
        self.HEXAGON.append(start)
        self.ALL_HEXAGON.append(start)
        self.ALL_HEXAGON_POS.append(start.pos)
        self.all_pos_tree = spatial.cKDTree(self.ALL_HEXAGON_POS, leafsize=10)
        start.updateFreePos(self)
#        ids,pts,colors,w=self.getFreePos()
#        self.displayWeights(pts,colors)
        return start 

    def placeHexa(self,start,idc,edge_id,pos,mata):
        i=edge_id
        start.nvisit[i]+=1
        T=self.all_pos_tree#spatial.cKDTree(self.ALL_HEXAGON_POS, leafsize=10)
        if self.display:
            sp1= helper.getObject("test1")
            if sp1 is None :
                sp1=helper.Sphere("test1",pos=numpy.array(pos))[0]
            helper.setTranslation(sp1,numpy.array(pos))
            helper.changeObjColorMat(sp1,[start.nvisit[i]/2.0,0,0])
        distance,indice = T.query(pos)    
        if distance > self.distance_threshold[1] or distance < self.distance_threshold[0]:#
            print "rejected too close ",distance,self.distance_threshold
            if start.nvisit[i] >= 2 :
                start.free_pos[i]=0
            return None
        #test angle edge         
        m2=matrix(mata.transpose())
        h=Hexa("hexa_"+str(idc)+"_"+str(i),pos,self.hexamer_radius,m2.T)#(m1*m2).T)
        self.HEXAGON.append(h)
        self.ALL_HEXAGON.append(h)
        self.ALL_HEXAGON_POS.append(pos)
        self.all_pos_tree = spatial.cKDTree(self.ALL_HEXAGON_POS, leafsize=10)
#        print "created hexa_"+str(idc)+"_"+str(i)
#        start.children[i]=h
#        start.free_pos[i]=0
        h.updateFreePos(self)
        return start    
        
    def nextHexa(self):
        for h in range(len(self.HEXAGON)):
            self.HEXAGON[h].updateFreePos(self)
        self.clearOccupied()
        ids,pts,colors,w = self.getFreePos(weight=True)
        if not len(pts):
            return 2
        i=self.getRndWeighted(w=w)#getSubWeighted()#getRndWeighted()
#        i=randrange(len(pts))
#        idh,edgid = ij(i)
        if self.inputMesh != None :
            indice,pos,rotMat,mata=self.computePosRot(pts[i])
        elif self.constraintMesh != None :
            start=self.HEXAGON[ids[i][0]]
            indice,pos,reject,mata=self.constrainPosRot(pts[i],h=start,edge=ids[i][1])
            indicetouse = self.ALL_HEXAGON.index(start)
            if reject :
                i=ids[i][1]
                start.nvisit[i]+=1
                if start.nvisit[i] >= 2 :
                    start.free_pos[i]=0
                return 1
        else :
            pos = pts[i]
            rotMat = numpy.identity(4)
            mata = numpy.identity(4)
        start=self.HEXAGON[ids[i][0]]
        indicetouse = self.ALL_HEXAGON.index(start)
        print "step ",start.name,indicetouse,ids[i][1],start.free_pos
        new_start=self.placeHexa(start,indicetouse,ids[i][1],pos,mata)
        if new_start is not None :
            start = new_start
            return 0
        else :
            return 1

    def tile(self,doarea=False,countThreshold=10,percantage_surface=100):
        #main loop
        t1 = time()
        done = False
        asa=0
        asac=1500
        countH=0
        parea=0
        while not done:
            helper.progressBar(progress=int((float(asa) / asac)*100.0),
                           label="i "+str(asa)+" / "+str(asac)+" _ "+str(countH) +" _ "+str(parea) )
            countH=len(self.ALL_HEXAGON)
            if not len(self.HEXAGON) :
                print "stop len hexagon"
                break                    
            if asa > asac :
                print "stop asac"
                done = True
                break      
            if self.display:
                helper.deleteObject("W_iSpheres");
            res  = self.nextHexa()
            if self.display:
                helper.update()
            if res == 1 :
#                startpos = start.pos
                asa+=1
                continue
            elif res == 2 :
                break
            else :
                if doarea:
                    self.surface_covered = self.hexa_area*len(self.ALL_HEXAGON)
                    parea = (self.surface_covered*100.0)/self.area
                    if parea >= percantage_surface :
                        print str(parea)+" % is covered,stop" , percantage_surface,  parea >= percantage_surface           
                        break
                else :
                    if countH > countThreshold :
                        print "stop countH"
                        break 
                if self.display:
                    tiles.sphs=[];
                    ids,pts,colors,w=self.getFreePos();
                    self.displayWeights(pts,colors)
                asa+=1        
        print "elapsed time is ",time()-t1,len(self.ALL_HEXAGON)
        if doarea:
            self.surface_covered = self.hexa_area*len(self.ALL_HEXAGON)
            print (self.surface_covered*100.0)/self.area," % is covered"
        helper.resetProgressBar()
        helper.progressBar(label="Tilling Complete in "+str(time()-t1)+" secs")

    def save(self,filename):
        #export [[567.30621337890625, -76.661155700683594, 318.21942138671875], [[0.42438143491744995, 3.5741718014068754e-09, 0.9054835279589204, 0.0], [0.23323895037174225, 0.96625566434317833, -0.10931429396496474, 0.0], [-0.87492853403091431, 0.25758499087995512, 0.41006097498201871, 0.0], [0.0, 0.0, 0.0, 1.0]], "MA_hexagone", 1, 0]
        results=[]        
        for h in self.ALL_HEXAGON:
            results.append([h.pos,h.mat.tolist(),"MA_hexagone",1,0])
        with open(filename, 'w') as fp :#doesnt work with symbol link ?
            json.dump(results,fp)

    def clear(self):
        helper.deleteObject("all_spheres") 
        helper.deleteObject("all_hex") 
        
def getFaceNormals(vertices, faces,fillBB=None):
    """compute the face normal of the compartment mesh"""
    normals = []
    areas = [] #added by Graham
    face = [[0,0,0],[0,0,0],[0,0,0]]
    v = [[0,0,0],[0,0,0],[0,0,0]]        
    for f in faces:
        for i in range(3) :
            face [i] = vertices[f[i]]
        for i in range(3) :
            v[0][i] = face[1][i]-face[0][i]
            v[1][i] = face[2][i]-face[0][i]                
        normal = vcross(v[0],v[1])
        n = vlen(normal)
        if n == 0. :
            n1=1.
        else :
            n1 = 1./n
        normals.append( (normal[0]*n1, normal[1]*n1, normal[2]*n1) )
    return normals, areas #areas added by Graham

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
    
def rotax( a, b, tau, transpose=1 ):
    """
    Build 4x4 matrix of clockwise rotation about axis a-->b
    by angle tau (radians).
    a and b are sequences of 3 floats each
    Result is a homogenous 4x4 transformation matrix.
    NOTE: This has been changed by Brian, 8/30/01: rotax now returns
    the rotation matrix, _not_ the transpose. This is to get
    consistency across rotax, mat_to_quat and the classes in
    transformation.py
    when transpose is 1 (default) a C-style rotation matrix is returned
    i.e. to be used is the following way Mx (opposite of OpenGL style which
    is using the FORTRAN style)
    """

    assert len(a) == 3
    assert len(b) == 3
    if tau <= -2*pi or tau >= 2*pi:
        tau = tau%(2*pi)

    ct = cos(tau)
    ct1 = 1.0 - ct
    st = sin(tau)

    # Compute unit vector v in the direction of a-->b. If a-->b has length
    # zero, assume v = (1,1,1)/sqrt(3).

    v = [b[0]-a[0], b[1]-a[1], b[2]-a[2]]
    s = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]
    if s > 0.0:
        s = sqrt(s)
        v = [v[0]/s, v[1]/s, v[2]/s]
    else:
        val = sqrt(1.0/3.0)
        v = (val, val, val)

    rot = numpy.zeros( (4,4), 'f' )
    # Compute 3x3 rotation matrix

    v2 = [v[0]*v[0], v[1]*v[1], v[2]*v[2]]
    v3 = [(1.0-v2[0])*ct, (1.0-v2[1])*ct, (1.0-v2[2])*ct]
    rot[0][0]=v2[0]+v3[0]
    rot[1][1]=v2[1]+v3[1]
    rot[2][2]=v2[2]+v3[2]
    rot[3][3] = 1.0;

    v2 = [v[0]*st, v[1]*st, v[2]*st]
    rot[1][0]=v[0]*v[1] * ct1-v2[2]
    rot[2][1]=v[1]*v[2] * ct1-v2[0]
    rot[0][2]=v[2]*v[0] * ct1-v2[1]
    rot[0][1]=v[0]*v[1] * ct1+v2[2]
    rot[1][2]=v[1]*v[2] * ct1+v2[0]
    rot[2][0]=v[2]*v[0] * ct1+v2[1]

    # add translation
    for i in (0,1,2):
        rot[3][i] = a[i]
    for j in (0,1,2):
        rot[3][i] = rot[3][i]-rot[j][i]*a[j]
    rot[i][3]=0.0

    if transpose:
        return rot
    else:
        return numpy.transpose(rot)

def rotVectToVect(vect1, vect2, i=None):
    """returns a 4x4 transformation that will align vect1 with vect2
vect1 and vect2 can be any vector (non-normalized)
"""
    v1x, v1y, v1z = vect1
    v2x, v2y, v2z = vect2
    
    # normalize input vectors
    norm = 1.0/sqrt(v1x*v1x + v1y*v1y + v1z*v1z )
    v1x *= norm
    v1y *= norm
    v1z *= norm    
    norm = 1.0/sqrt(v2x*v2x + v2y*v2y + v2z*v2z )
    v2x *= norm
    v2y *= norm
    v2z *= norm
    
    # compute cross product and rotation axis
    cx = v1y*v2z - v1z*v2y
    cy = v1z*v2x - v1x*v2z
    cz = v1x*v2y - v1y*v2x

    # normalize
    nc = sqrt(cx*cx + cy*cy + cz*cz)
    if nc==0.0:
        return [ [1., 0., 0., 0.],
                 [0., 1., 0., 0.],
                 [0., 0., 1., 0.],
                 [0., 0., 0., 1.] ]

    cx /= nc
    cy /= nc
    cz /= nc
    
    # compute angle of rotation
    if nc<0.0:
        if i is not None:
            print ('truncating nc on step:', i, nc)
        nc=0.0
    elif nc>1.0:
        if i is not None:
            print ('truncating nc on step:', i, nc)
        nc=1.0
        
    alpha = asin(nc)
    if (v1x*v2x + v1y*v2y + v1z*v2z) < 0.0:
        alpha = pi - alpha

    # rotate about nc by alpha
    # Compute 3x3 rotation matrix

    ct = cos(alpha)
    ct1 = 1.0 - ct
    st = sin(alpha)
    
    rot = [ [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.] ]


    rv2x, rv2y, rv2z = cx*cx, cy*cy, cz*cz
    rv3x, rv3y, rv3z = (1.0-rv2x)*ct, (1.0-rv2y)*ct, (1.0-rv2z)*ct
    rot[0][0] = rv2x + rv3x
    rot[1][1] = rv2y + rv3y
    rot[2][2] = rv2z + rv3z
    rot[3][3] = 1.0;

    rv4x, rv4y, rv4z = cx*st, cy*st, cz*st
    rot[0][1] = cx * cy * ct1 - rv4z
    rot[1][2] = cy * cz * ct1 - rv4x
    rot[2][0] = cz * cx * ct1 - rv4y
    rot[1][0] = cx * cy * ct1 + rv4z
    rot[2][1] = cy * cz * ct1 + rv4x
    rot[0][2] = cz * cx * ct1 + rv4y
    #rot[2][:3]=[0,0,1]
    return rot

def getCoordHexagone(radius,center=[0,0,0]):
    D=2*radius
    dY=D/4.0
    dX=math.sqrt(3)*dY
    pcpalAxis=[0,0,1]
    hexc=coordinates=numpy.array([[0,2.0*dY,0.0],#up
                [dX,dY,0],#right-up
                [dX,-dY,0],#right-down
                [0,-2.0*dY,0],#down
                [-dX,-dY,0],#left-down
                [-dX,dY,0]#left-up
                 ])
    #neighboors_pos
    #hexc = hexc + numpy.array(center)
    neighboors_pos=numpy.zeros([6,3])
    for i in range(5):
        neighboors_pos[i]=hexc[i]+hexc[i+1]
    neighboors_pos[5] =hexc[5]+hexc[0]  
    return pcpalAxis,coordinates,neighboors_pos

def distanceEdge(hexo, edge_id, coordA,coordB):
    #{0:3,1:4,2:5,3:0,4:1,5:2}
    #neighbor = self.children[edge_id]
    edge1=coordA[hexo.edge_nb[edge_id][0][0]]-coordA[hexo.edge_nb[edge_id][1][0]]
    edge2=coordB[hexo.edge_nb[edge_id][0][1]]-coordB[hexo.edge_nb[edge_id][1][1]]
    a = angle_between_vectors([edge1], [edge2], axis=1)[0]
    D1 = vdistance(coordA[hexo.edge_nb[edge_id][0][0]],coordB[hexo.edge_nb[edge_id][0][1]])
    D2 = vdistance(coordA[hexo.edge_nb[edge_id][1][0]],coordB[hexo.edge_nb[edge_id][1][1]])
#    print (edge_id,"edge distances and angle ",D1,D2,a,math.degrees(a))      
    rotMat = numpy.array( rotVectToVect(edge2, edge1 ), 'f') 
    return a,D1,D2,rotMat

def computePosRot(helper,next_pos,tree):
    distance,indice = tree.query(next_pos)
    vx, vy, vz = v1 = [0,0,1]#
    newv1=v1#helper.ApplyMatrix([v1,],pm)[0]#with threanslation ?
    if rc.Intersect(helper.FromVec(numpy.array(next_pos)*100.0), helper.FromVec(next_pos)*-1, 100000):
        ray_result = rc.GetNearestIntersection()
        n = ray_result["f_normal"].GetNormalized()
        v2 = helper.ToVec(n)# vn[indice]
    else :
        if is_sphere:
            v2 = next_pos
        else :   
            v2 = vn[indice]#next_pos#vn[indice]
    vnax=numpy.zeros((4,4))
    vnax[0][:3] = vnorm(numpy.cross(v2,[0,1,0]))
    vnax[1][:3] = vnorm(numpy.cross(v2,vnax[0][:3]))
    vnax[2][:3] = vnorm(v2)
    euler = euler_from_matrix(vnax)
    rotMat = euler_matrix(euler[0],euler[1],0.0).transpose()
    mat = numpy.array(rotMat)
    mat[:3, 3] = next_pos
    if rc.Intersect(helper.FromVec(numpy.array(next_pos)*10.0), helper.FromVec(next_pos)*-1, 100000):
        ray_result = rc.GetNearestIntersection()
        newposition = helper.ToVec(ray_result["hitpos"])
    else :
        newposition=v[indice]#(next_pos+tr)#-pm.transpose()[3,:3]#+tr#pos+p
    mat = numpy.array(rotMat)#matrix(R)*matrix(rotMat))#
    mat[:3, 3] = newposition#+tr
#    helper.setObjectMatrix(ipol[1],mat.transpose())  
    #test the new post 
    return indice,None,None,newposition,rotMat,mat
    
def one_next(helper,next_pos,tree,rc,v,vn,old_c,old_n,oi,ni,pm,is_sphere=False):
#    ax,old_c,old_n=getCoordHexagone(55.0)    
    distance,indice = tree.query(next_pos)
    #align
    vx, vy, vz = v1 = [0,0,1]#
#    R=rotax( [0,0,0], v1, random()*math.radians(3.0), transpose=1 )
    newv1=v1#helper.ApplyMatrix([v1,],pm)[0]#with threanslation ?
    #next_point should be outisde
    if DEBUG :
        sp1= helper.getObject("test1")
        if sp1 is None :
            sp1=helper.Sphere("test1",pos=numpy.array(next_pos)*100.0)[0]
        helper.setTranslation(sp1,numpy.array(next_pos)*10)
        sp2= helper.getObject("test2")
        if sp2 is None :
            sp2=helper.Sphere("test2",pos=numpy.array(next_pos)*-1.0)[0]
        helper.setTranslation(sp2,numpy.array(next_pos)*-1.0)
        l1= helper.getObject("line1")
        if l1 is None :
            l1=helper.spline("line1", [numpy.array(next_pos)*10,numpy.array(next_pos)*-1.0])[0]
        helper.update_spline(l1,[numpy.array(next_pos)*10,numpy.array(next_pos)*-1.0])
        helper.update()
    if rc.Intersect(helper.FromVec(numpy.array(next_pos)*100.0), helper.FromVec(next_pos)*-1, 100000):
        ray_result = rc.GetNearestIntersection()
        n = ray_result["f_normal"].GetNormalized()
        v2 = helper.ToVec(n)# vn[indice]
#        print ray_result
#    next_pos = v[indice]
    else :
        if is_sphere:
            v2 = next_pos
        else :   
            v2 = vn[indice]#next_pos#vn[indice]
#    v2 = vn[indice]
#    rotMat =superimposition_matrix(newv1,v2)
    vnax=numpy.zeros((4,4))
    vnax[0][:3] = vnorm(numpy.cross(v2,[0,1,0]))
    vnax[1][:3] = vnorm(numpy.cross(v2,vnax[0][:3]))
    vnax[2][:3] = vnorm(v2)
    euler = euler_from_matrix(vnax)
    rotMat = euler_matrix(euler[0],euler[1],0.0).transpose()
    
#    m=numpy.identity(4)
#    m[:3,:3]=pm[:3,:3]
#    rotMat = matrix(rotMat)*matrix(m)
#    rotMat = numpy.array( rotVectToVect(newv1, v2 ), 'f')
    #add some random motion around vx, vy, vz
#    rotMat = numpy.identity(4)
#    euler = euler_from_matrix(rotMat)
#    print (math.degrees(euler[0]),math.degrees(euler[1]),math.degrees(euler[2]),)
    #filter the rotatio, keep only X-Y
#    rotMat = numpy.identity(4)
    mat = numpy.array(rotMat)
    mat[:3, 3] = next_pos
##    nr=helper.ApplyMatrix(n,rotMat)
#    ahsphs=helper.instancesSphere("sphereHexagoneA",newci,[1,]*6,s,[[0,1,0]],None,parent=parent)
#    asphs=helper.instancesSphere("sphereHexagoneNA",newni,[1,]*6,s,[[0,0,1]],None,parent=parent)    
    #new position, compareto current ?
#    Cnew = newci[ni] - newci[ni+1]#edge1
#    Cold = old_c[oi] - old_c[oi-1]#edge2
#    rotMat2 = numpy.array( rotVectToVect(Cnew, Cold ), 'f')
#    mat2 = numpy.array(rotMat2).transpose()
#    mat[:3, 3] = next_pos-pm[3,:3]
#    newci2=helper.ApplyMatrix(newci,mat2)
#    newni2=helper.ApplyMatrix(newni,mat2)    
    tr = numpy.array([0,0,0])
    if rc.Intersect(helper.FromVec(numpy.array(next_pos)*10.0), helper.FromVec(next_pos)*-1, 100000):
        ray_result = rc.GetNearestIntersection()
        newposition = helper.ToVec(ray_result["hitpos"])
    else :
        newposition=v[indice]#(next_pos+tr)#-pm.transpose()[3,:3]#+tr#pos+p
#    helper.setTranslation(sphs[4],pos=newposition,absolue=True)
    #need rotation on [0,0,1]
#    newci=helper.ApplyMatrix(old_c,mat)
#    newni=helper.ApplyMatrix(old_n,mat)
#    C1new=newci[ni]
#    C5old=old_c[oi]
#    tr=C5old-C1new
    mat = numpy.array(rotMat)#matrix(R)*matrix(rotMat))#
    mat[:3, 3] = newposition#+tr
#    helper.setObjectMatrix(ipol[1],mat.transpose())  
    #test the new post 
    return indice,None,None,newposition,rotMat,mat

def testHex(distance,indice,i,start,mata):
    if distance <= 75.0 :
        return False
    
def getHex(start,idc,edge_id,pos,mata):
    i=edge_id
    start.nvisit[i]+=1
    T=spatial.cKDTree(ALL_HEXAGON_POS, leafsize=10)
    sp1= helper.getObject("test1")
    if sp1 is None :
        sp1=helper.Sphere("test1",pos=numpy.array(pos))[0]
    helper.setTranslation(sp1,numpy.array(pos))
    helper.changeObjColorMat(sp1,[start.nvisit[i]/2.0,0,0])
    distance,indice = T.query(pos)    
    if distance <= 75.0 :#
        if start.nvisit[i] >= 2 :
            start.free_pos[i]=0
        return None
    #test angle edge         
    m2=matrix(mata.transpose())
    h=Hexa("hexa_"+str(idc)+"_"+str(i),pos,hexamer_radius,m2.T)#(m1*m2).T)
    HEXAGON.append(h)
    ALL_HEXAGON.append(h)
    ALL_HEXAGON_POS.append(pos)
#    helper.update()
    print "created hexa_"+str(idc)+"_"+str(i)
    start.children[i]=h
    start.free_pos[i]=0
    return start    
    

def oneHex(start,edge_id,tree,rc,v,vn,is_sphere,n=99,idc=0):
    if start is None :
        return
#    m1=matrix(start.mat.transpose())
    i=edge_id
    start.nvisit[i]+=1
    T=spatial.cKDTree(ALL_HEXAGON_POS, leafsize=10)
    ind1,newci,newni,pos,rotMat,mata=one_next(helper,start.tr_n[i],
                              tree,rc,v,vn,start.tr_c,start.tr_n,start.nb[i][0],
                              start.nb[i][1],start.mat,is_sphere=is_sphere)
    sp1= helper.getObject("test1")
    if sp1 is None :
        sp1=helper.Sphere("test1",pos=numpy.array(pos))[0]
    helper.setTranslation(sp1,numpy.array(pos))
    helper.changeObjColorMat(sp1,[start.nvisit[i]/2.0,0,0])
    distance,indice = T.query(pos)    
    if distance <= 77.0 :#
        if start.nvisit[i] >= 2 :
            start.free_pos[i]=0
        return None
    #test angle edge         
    m2=matrix(mata.transpose())
    h=Hexa("hexa_"+str(idc)+"_"+str(i),pos,hexamer_radius,m2.T)#(m1*m2).T)
    HEXAGON.append(h)
    ALL_HEXAGON.append(h)
    ALL_HEXAGON_POS.append(pos)
#    helper.update()
    start.children[i]=h
    start.free_pos[i]=0
    return start    

def oneRing(start,tree,rc,v,vn,is_sphere,n=99,idc=0):
    if start is None :
        return
    m1=matrix(start.mat.transpose())
    count=0
    while len(numpy.nonzero(start.free_pos)[0]):
        i=start.nextFree()
        new_start = oneHex(start,i,tree,rc,v,vn,is_sphere,n=n,idc=idc)
        if new_start is None :
            if n == count :
                break
            count+=1
            continue
        if n == count :
            break
        count+=1
        start = new_start
    start.updateChildrenFreePos()
    return start

def ij(u,width=6):
    return [u/width, u%width]
   
def u(i,j,width=6):
    return width * i + j

def displayWeights(points,colors,sphs=[]):
    parent = helper.getObject("W_iSpheres")
    if parent is None :
        parent = helper.newEmpty("W_iSpheres")
    s = helper.getObject("W_base_sphere")
    if s is None :
        s,ms = helper.Sphere("W_base_sphere",radius=5.0,res=12,parent=parent)
    sphs=helper.updateInstancesSphere("W_spheresH",sphs,points,[1.0,]*len(points),s,
                        colors,None,parent=parent,delete=False)
    return sphs


def getSubWeighted(startpos=None):
    """
    From http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/
    This method is about twice as fast as the binary-search technique, 
    although it has the same complexity overall. Building the temporary 
    list of totals turns out to be a major part of the functions runtime.
    This approach has another interesting property. If we manage to sort 
    the weights in descending order before passing them to 
    weighted_choice_sub, it will run even faster since the random 
    call returns a uniformly distributed value and larger chunks of 
    the total weight will be skipped in the beginning.
    """
    
    weights,pts,colors = weightFreePos(startpos=startpos)#numpy.take(self.weight,listPts)
    rnd = random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i
 

def oneStart(vi):
    pos=v[vi]
    #m=helper.getTransformation(hexamer_obj)
    hax = numpy.identity(4)
    v2 = vn[vi]#pos
    vnax=numpy.zeros((4,4))
    vnax[0][:3] = vnorm(numpy.cross(v2,hax[1][:3]))
    vnax[1][:3] = vnorm(numpy.cross(v2,vnax[0][:3]))
    vnax[2][:3] = vnorm(v2)
    euler = euler_from_matrix(vnax)
    rotMat = euler_matrix(euler[0],euler[1],0.0)#.transpose()
    mat = numpy.array(rotMat)
    mat[3, :3] = pos
    start=Hexa("hexa_"+str(vi),pos,hexamer_radius,mat.transpose())
    HEXAGON.append(start)
    ALL_HEXAGON.append(start)
    ALL_HEXAGON_POS.append(start.pos)
    return start 


def pickNext(startpos):
    #weight according distance from starting point
    #grab all center
    centers = [h.pos for h in HEXAGON]
    ctree = spatial.cKDTree(centers, leafsize=10)
    distance,ind = ctree.query(startpos)
    return ind

def oneStep(start,sphs=[]):
    for h in range(len(HEXAGON)):
        HEXAGON[h].updateFreePos()
    clearOccupied()
    ids,pts,w,HEXAGON = getFreePos(HEXAGON)
#    print ids
#    i = randrange(len(pts))
    i=getRndWeighted(w=w)#getSubWeighted()#getRndWeighted()
#    idh,edgid = ij(i)
    indice,a,b,pos,rotMat,mata=computePosRot(helper,pts[i],tree)
    start=HEXAGON[ids[i][0]]
    new_start=getHex(start,ids[i][0],ids[i][1],pos,mata)
    if new_start is not None :
        start = new_start
#        start.updateFreePos()
        weights,pts,colors = weightFreePos()
        sphs= displayWeights(pts,colors,sphs=sphs)
    return start,new_start    
#    start,new_start=oneStep()
#def recursiveRing(start,inc,r):
#    for i in range(6):
#        start.children[i]=oneRing(start.children[i],idc=str(r)+"_"+str(inc)+"_"+str(i))
   
#a,extendedC,extendedN=getCoordHexagone(n[0])
#extendedN=n+numpy.array(n[0])

#parent = helper.newEmpty("iSpheres")
s,ms = helper.Sphere("sphere",radius=5.0,res=12)
#sphs=helper.instancesSphere("sphereHexagone",c,[1,]*6,s,[[0,1,0]],None,parent=parent)
#sphs=helper.instancesSphere("sphereHexagoneN",n,[1,]*6,s,[[0,0,1]],None,parent=parent)
#sphs=helper.instancesSphere("sphereHexagoneNE",extendedN,[1,]*6,s,[[1,0,0]],None,parent=parent)

hexamer_obj = helper.getObject("hex")
hexamer_pos=[]
hexamer_free_nb=[1,1,1,1,1,1]#all free
hexamer_radius=55#inner radus 0f 40, distance h-h ~80
ax,c,n=getCoordHexagone(hexamer_radius)
hexamer_local_pos=[[],]
from c4d import utils
sel = helper.getCurrentSelection()
if not len(sel):
    input_mesh = helper.getObject("inputMesh")
else :
    input_mesh = sel[0]

c_mesh = helper.getObject("envMesh")
#
input_mesh = None
nmodels = 1
for n in range(nmodels):
    tiles=Tilling([],[],[],input_mesh,hexamer_radius,init_seed=n,inputConstraint=c_mesh)   
#    tiles.firstRandomStart()
#    tiles.firstRandomStart()
#    tiles.firstRandomStart()
    tiles.randomStart([0,0,0])
    tiles.display=True
    tiles.tile(countThreshold=500,doarea=False,percantage_surface=60.0)
#    tiles.save("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/ma_models/model_"+str(n)+".json")
#    tiles.clear()
    
#export the model

#600 -> 46 area 60%
#448 -> 55 area 60%

#helper.deleteObject("W_iSpheres");start,new_start  = tiles.nextHexa();tiles.sphs=[];ids,pts,colors,w=tiles.getFreePos();tiles.displayWeights(pts,colors)
#tiles.tile(10)
#print len(tiles.ALL_HEXAGON)
#        tiles.sphs=[];ids,pts,colors,w=tiles.getFreePos();tiles.displayWeights(pts,colors)
#        tiles.displayWeights(pts,colors)

#p = helper.advance_randpoint_onsphere(650,marge=math.radians(45.0),vector=[0,1,0])
#d,ind=tree.query(p)
#start=oneStart(ind)
#startpos=start.pos
#sphs=[]
#start=oneStep(start,sphs=sphs)
#p = helper.advance_randpoint_onsphere(650,marge=math.radians(90.0),vector=[0,1,0])
#d,ind=tree.query(p)
#start=oneStart(ind)
#p = helper.advance_randpoint_onsphere(650,marge=math.radians(90.0),vector=[0,1,0])
#d,ind=tree.query(p)
#start=oneStart(ind)
#p = helper.advance_randpoint_onsphere(650,marge=math.radians(90.0),vector=[0,1,0])
#d,ind=tree.query(p)
#start=oneStart(ind)

#start=oneRing(start,tree,rc,v,vn,is_sphere,idc="S2")#n=1,
##print HEXAGON,len(HEXAGON)
#for i in range(6):
#    start.children[i]=oneRing(start.children[i],tree,rc,v,vn,is_sphere,idc=i)
#for j in range(6):    
#    for k in range(6):
#        start.children[j].children[k]=oneRing(start.children[j].children[k],tree,rc,v,vn,is_sphere,idc=str(j)+"_"+str(k))
#for j in range(6):    
#    for k in range(6):
#        for l in range(6):
#            start.children[j].children[k].children[l]=oneRing(start.children[j].children[k].children[l],tree,rc,v,vn,is_sphere,idc=str(j)+"_"+str(k)+"_"+str(l))

#done=True
#countH=0
##1500 Matrix ->  350
##2500 GAG 
#countThreshold=50#414
#safetyThreshold=300
#safety=0
#asa=0
#asac=1500
#sphs=[]
#while not done:
#    helper.progressBar(progress=int((float(asa) / asac)*100.0),
#                   label="i "+str(asa)+" / "+str(asac)+" _ "+str(countH) )
#    countH=len(ALL_HEXAGON)
#    if countH > countThreshold :
#        print "stop countH"
#        break
#    if not len(HEXAGON) :
#        print "stop len hexagon"
#        break
#    if safety > safetyThreshold :
#        print "stop safety"
#        done = True
#        break
#    if asa > asac :
#        print "stop asac"
#        done = True
#        break       
#    for h in range(len(HEXAGON)):
#        HEXAGON[h].updateFreePos()
#    clearOccupied()
#    ids,pts,w,HEXAGON = getFreePos(HEXAGON)
##    print ids
##    i = randrange(len(pts))
#    i=getRndWeighted(w=w)#getSubWeighted()#getRndWeighted()
##    idh,edgid = ij(i)
#    indice,a,b,pos,rotMat,mata=computePosRot(helper,pts[i],tree)
#    start=HEXAGON[ids[i][0]]
#    new_start=getHex(start,ids[i][0],ids[i][1],pos,mata)
#    #print "sel ",i,idh,edgid
##    idh = randrange(len(HEXAGON))#  pickNext(startpos)#pickNext(startpos)#randrange(len(HEXAGON))
##    edgid = randrange(6)
##    start = HEXAGON[idh]
##    start.updateFreePos()
##    new_start = oneHex(start,edgid,tree,rc,v,vn,is_sphere,n=99,idc=idh)
#    helper.update()
#    if new_start == None :
#        startpos = start.pos
#        asa+=1
#        safety+=1
#        continue
#    print "step ",start.name,ids[i][0],ids[i][1],start.free_pos
#    safety = 0
##    start = new_start
##    start.updateChildrenFreePos()
##    start.updateFreePos()
##    weights,pts,colors = weightFreePos(startpos=startpos)
#    for h in range(len(HEXAGON)):
#        HEXAGON[h].updateFreePos()    
#    clearOccupied()
#    sphs= displayWeights(pts,[[1,0,0],]*len(pts),sphs=sphs)
#    countH=len(ALL_HEXAGON)
#    if countH > countThreshold :
#        print "stop len countH"
#        break
#    if not len(HEXAGON) :
#        print "stop len HEXAGON"
#        break
#    asa+=1

#    start = oneRing(start,tree,rc,v,vn,is_sphere,idc=countH)
#    countH=len(ALL_HEXAGON)
#    if countH > countThreshold :
#        break
#    if not len(HEXAGON) :
#        break
#    if (safety % 10) == 0 : 
#        helper.progressBar(progress=int((float(safety) / safetyThreshold)*100.0),
#                       label="i "+str(safety)+" / "+str(safetyThreshold)+" _ "+str(countH) )
#        startpos = HEXAGON[randrange(len(HEXAGON))].pos
#        safety = 0
#    safety+=1
#    asa+=1
# check if any frePos remain ?
#print len(ALL_HEXAGON)
#[h.updateChildrenFreePos() for h in HEXAGON]

#a,d1,d2,rMat2 = distanceEdge(start, 0, start.tr_c,start.children[0].tr_c)
#rMat1 = numpy.identity(4)
#rMat1[:3,:3] = start.children[0].mat[:3,:3]
#m=matrix(rMat1.transpose())*matrix(rMat2.transpose());mata = numpy.array(m);mata[3, :3] = start.children[0].pos;helper.setObjectMatrix(start.children[0].inst,mata)
#        mata = numpy.array(m).transpose()#matrix(R)*matrix(rotMat))#
#        mata[:3, 3] = pos#+tr        


#    [h.updateChildrenFreePos() for h in HEXAGON]
#[h.updateChildrenFreePos(cutoff=1.0) for h in HEXAGON]

#
#r=recursiveRing(start,0,0)

##start.children[0]=oneRing(start.children[0],idc=0)
##start.children[1]=oneRing(start.children[1],idc=1)
##start.children[0]=oneRing(start.children[0],idc=0)
#for i in range(6):
#    start.children[i]=oneRing(start.children[i],tree,rc,v,vn,is_sphere,idc=i)
#for j in range(6):    
#    for k in range(6):
#        start.children[j].children[k]=oneRing(start.children[j].children[k],tree,rc,v,vn,is_sphere,idc=str(j)+"_"+str(k))
#    
#for j in range(6):    
#    for k in range(6):
#        for l in range(6):
#            start.children[j].children[k].children[l]=oneRing(start.children[j].children[k].children[l],tree,rc,v,vn,is_sphere,idc=str(j)+"_"+str(k)+"_"+str(l))

#for j in range(6):    
#    for k in range(6):
#        for l in range(6):
#            for m in range(6):
#                start.children[j].children[k].children[l].children[m]=oneRing(start.children[j].children[k].children[l].children[m],tree,rc,v,vn,is_sphere,idc=str(j)+"_"+str(k)+"_"+str(l)+"_"+str(m))
#
#start = start.children[j].children[0].children[l].children[0]#.children[0]
#start=oneRing(start,tree,rc,v,vn,is_sphere,idc="AX")
#for i in range(6):
#    start.children[i]=oneRing(start.children[i],tree,rc,v,vn,is_sphere,idc="AX"+str(i))
#for j in range(6):    
#    for k in range(6):
#        start.children[j].children[k]=oneRing(start.children[j].children[k],tree,rc,v,vn,is_sphere,idc="AX"+str(j)+"_"+str(k))
#
#start = start.children[j].children[1].children[l].children[1]#.children[0]
#start=oneRing(start,tree,rc,v,vn,is_sphere,idc="BX")
#for i in range(6):
#    start.children[i]=oneRing(start.children[i],tree,rc,v,vn,is_sphere,idc="BX"+str(i))
#for j in range(6):    
#    for k in range(6):
#        start.children[j].children[k]=oneRing(start.children[j].children[k],tree,rc,v,vn,is_sphere,idc="BX"+str(j)+"_"+str(k))
#
#start = start.children[j].children[2].children[l].children[2]#.children[0]
#start=oneRing(start,tree,rc,v,vn,is_sphere,idc="CX")
#for i in range(6):
#    start.children[i]=oneRing(start.children[i],tree,rc,v,vn,is_sphere,idc="CX"+str(i))
#for j in range(6):    
#    for k in range(6):
#        start.children[j].children[k]=oneRing(start.children[j].children[k],tree,rc,v,vn,is_sphere,idc="CX"+str(j)+"_"+str(k))

#for j in range(6):    
#    for k in range(6):
#        for l in range(6):
#            for m in range(6):
#                for n in range(6):
#                    start.children[j].children[k].children[l].children[m].children[n]=oneRing(start.children[j].children[k].children[l].children[m].children[n],idc=str(j)+"_"+str(k)+"_"+str(l)+"_"+str(m)+"_"+str(n))
#for j in range(6):    
#    for k in range(6):
#        for l in range(6):
#            for m in range(6):
#                for n in range(6):
#                    for o in range(6):
#                        start.children[j].children[k].children[l].children[m].children[n].children[o]=oneRing(start.children[j].children[k].children[l].children[m].children[n].children[o],idc=str(j)+"_"+str(k)+"_"+str(l)+"_"+str(m)+"_"+str(n)+"_"+str(o))
#for j in range(6):
#    #new starting point 
#    start = ahexa.children[j] 
#    m1=matrix(start.mat.transpose())
#    ahexaic=[]
#    print j,start.free_pos,start.name
#    while len(numpy.nonzero(start.free_pos)[0]):
#        i=start.nextFree()
#        ind1,newci,newni,pos,rotMat,mata=one_next(helper,start.tr_n[i],
#                                                  tree,vn,start.tr_c,start.tr_n,start.nb[i][0],
#                                            start.nb[i][1],start.mat,is_sphere=is_sphere)
#        m2=matrix(mata.transpose())
#        h=Hexa("hexa_"+str(j)+"_"+str(i),pos,hexamer_radius,(m1*m2).T)
#        start.children[i]=h
#        start.free_pos[i]=0
#    start.updateChildrenFreePos()
 
#for j in range(6):    
#    for k in range(6):
#        start = ahexa.children[j].children[k]
#        if start is None : continue
#        m1=matrix(start.mat.transpose())
#        print j,k,start.free_pos,start.name
#        while len(numpy.nonzero(start.free_pos)[0]):
#            i=start.nextFree()
#            ind1,newci,newni,pos,rotMat,mata=one_next(helper,start.tr_n[i],
#                                                      tree,vn,start.tr_c,start.tr_n,start.nb[i][0],
#                                                start.nb[i][1],start.mat,is_sphere=is_sphere)
#            m2=matrix(mata.transpose())
#            h=Hexa("hexa_"+str(k)+"_"+str(j)+"_"+str(i),pos,hexamer_radius,(m1*m2).T)
#            start.children[i]=h
#            start.free_pos[i]=0
#        start.updateChildrenFreePos()
#        
##whats next pos+rot
##apply rotMat to c
#mat = numpy.array(rotMat)#.transpose()
#mat[:3, 3] = pos
#
#newc=helper.ApplyMatrix(c,mat)
#newn=helper.ApplyMatrix(n,mat)
#nr=helper.ApplyMatrix(n,rotMat)
#hsphs=helper.instancesSphere("sphereHexagone",newc,[1,]*6,s,[[0,1,0]],None,parent=parent)
#sphs=helper.instancesSphere("sphereHexagoneN",newn,[1,]*6,s,[[0,0,1]],None,parent=parent)
##planar tiling
#mat2 = numpy.array(rotMat)
#mat2[:3, 3] = mat[:3, 3] + nr[4]
#mat3 = numpy.array(rotMat)
#mat3[:3, 3] = mat[:3, 3] + nr[5]
#mat4 = numpy.array(rotMat)
#mat4[:3, 3] = mat[:3, 3] + nr[0]
#mat5 = numpy.array(rotMat)
#mat5[:3, 3] = mat[:3, 3] + nr[1]
#mat6 = numpy.array(rotMat)
#mat6[:3, 3] = mat[:3, 3] + nr[2]
#mat7 = numpy.array(rotMat)
#mat7[:3, 3] = mat[:3, 3] + nr[3]
#ipol=helper.instancePolygon("hexI", matrices=[mat,mat2,mat3,mat4,mat5,mat6,mat7],
#                            hmatrices=None, mesh=hexamer_obj,transpose=True)
##helper.instancePolygon("hexII", matrices=[mat2],hmatrices=None, mesh=hexamer_obj,transpose=True)
#M1=matrix(mat2.transpose())
##helper.instancePolygon("hexII", matrices=[mat],hmatrices=None, mesh=hexamer_obj,transpose=True)
#
##
##next
##adjustements
#m1=matrix(mat.transpose())
#ind1,newci,newni,pos,mata=one_next(helper,newn[4],tree,vn,newc,newn,5,1,mat)
#m2=matrix(mata.transpose())
#helper.setTranslation(sphs[4],pos=pos,absolue=True)
#ahsphs=helper.instancesSphere("sphereHexagoneA",newci,[1,]*6,s,[[0,1,0]],None,parent=parent)
##asphs=helper.instancesSphere("sphereHexagoneNA",newni,[1,]*6,s,[[0,0,1]],None,parent=parent)
#helper.setObjectMatrix(ipol[1],m1*m2)
#
##next
#ind2,newci,newni,pos,matb=one_next(helper,newn[5],tree,vn,newc,newn,5,3,mat)
#m3=matrix(matb.transpose())
#helper.setTranslation(sphs[5],pos=pos,absolue=True)
#ahsphs=helper.instancesSphere("sphereHexagoneB",newci,[1,]*6,s,[[0,1,0]],None,parent=parent)
##asphs=helper.instancesSphere("sphereHexagoneNB",newni,[1,]*6,s,[[0,0,1]],None,parent=parent)
#helper.setObjectMatrix(ipol[2],m1*m3)
#
##next
#ind3,newci,newni,pos,matb=one_next(helper,newn[0],tree,vn,newc,newn,1,3,mat)
#m3=matrix(matb.transpose())
#helper.setTranslation(sphs[0],pos=pos,absolue=True)
#ahsphs=helper.instancesSphere("sphereHexagoneB",newci,[1,]*6,s,[[0,1,0]],None,parent=parent)
##asphs=helper.instancesSphere("sphereHexagoneNB",newni,[1,]*6,s,[[0,0,1]],None,parent=parent)
#helper.setObjectMatrix(ipol[3],m1*m3)
#
#
##next
#ind4,newci,newni,pos,matb=one_next(helper,newn[1],tree,vn,newc,newn,2,4,mat)
#m3=matrix(matb.transpose())
#helper.setTranslation(sphs[1],pos=pos,absolue=True)
#ahsphs=helper.instancesSphere("sphereHexagoneB",newci,[1,]*6,s,[[0,1,0]],None,parent=parent)
##asphs=helper.instancesSphere("sphereHexagoneNB",newni,[1,]*6,s,[[0,0,1]],None,parent=parent)
#helper.setObjectMatrix(ipol[4],m1*m3)
##
##next
#ind5,newci,newni,pos,matb=one_next(helper,newn[2],tree,vn,newc,newn,2,0,mat)
#m3=matrix(matb.transpose())
#helper.setTranslation(sphs[2],pos=pos,absolue=True)
#ahsphs=helper.instancesSphere("sphereHexagoneC",newci,[1,]*6,s,[[0,1,0]],None,parent=parent)
##asphs=helper.instancesSphere("sphereHexagoneNB",newni,[1,]*6,s,[[0,0,1]],None,parent=parent)
#helper.setObjectMatrix(ipol[5],m1*m3)
#
##next
#ind6,newci,newni,pos,matb=one_next(helper,newn[3],tree,vn,newc,newn,3,1,mat)
#m3=matrix(matb.transpose())
#helper.setTranslation(sphs[3],pos=pos,absolue=True)
#ahsphs=helper.instancesSphere("sphereHexagoneD",newci,[1,]*6,s,[[0,1,0]],None,parent=parent)
##asphs=helper.instancesSphere("sphereHexagoneNB",newni,[1,]*6,s,[[0,0,1]],None,parent=parent)
#helper.setObjectMatrix(ipol[6],m1*m3)
#
#nor=[]
#nor.extend([v[len(v)/2],v[len(v)/2]+vn[len(v)/2]*20.0])
#nor.extend([v[ind1],v[ind1]+vn[ind1]*20.0])
#nor.extend([v[ind2],v[ind2]+vn[ind2]*20.0])
#nor.extend([v[ind3],v[ind3]+vn[ind3]*20.0])
#nor.extend([v[ind4],v[ind4]+vn[ind4]*20.0])
#nor.extend([v[ind5],v[ind5]+vn[ind5]*20.0])
#nor.extend([v[ind6],v[ind6]+vn[ind6]*20.0])
#parent = helper.newEmpty("normals")
#ahsphs=helper.instancesSphere("NORMALS",nor,[1,]*6,s,[[0,1,0]],None,parent=parent)
#

#pick a hexamer to drop
#pick where to drop ie pos 0-1-2-3-4-5 of closest hexagone position if free 

#
##next
#next_p=newn[5]
##find closest point + normal
#distance,indice = tree.query(next_p)
#
#va=next_p-pos
#vb=v[indice]-pos
#rotab = numpy.array( rotVectToVect(va, vb ), 'f')
#p=helper.ApplyMatrix([va,],rotab)[0]
#newposition=pos+p
#helper.setTranslation(sphs[5],pos=newposition,absolue=True)
#m = numpy.array(rotab)
#m[:3, 3] = newposition
##helper.setObjectMatrix(ipol[1],m.transpose())
#helper.setTranslation(ipol[2],pos=newposition,absolue=True)
#
#vx, vy, vz = v1 = [0,0,1]#
#v2 = vn[indice]
#rotMat = numpy.array( rotVectToVect(v1, v2 ), 'f')
#p=helper.ApplyMatrix([va,],rotMat)[0]
#newposition=pos+p
#mat = numpy.array(rotMat)
#mat[:3, 3] = newposition
#helper.setObjectMatrix(ipol[2],mat.transpose())
#
##rotMat = numpy.array( rotVectToVect(v1, v2 ), 'f')
##a=angle_between_vectors(v1,v2)#sign ?
##print a
##MR=rotax( newc[4], newc[5], a, transpose=1 )
##mr=matrix(MR)
##M2=matrix(numpy.identity(4))
##M2[:3,:3]=MR[:3,:3]
##M2[3, :3] = mat2[:3, 3]
##helper.setObjectMatrix(ipol[1],M2)
##p=helper.ApplyMatrix([next_p,],MR)[0]

#                execfile("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/testTiling.py")
#execfile("/opt/data/dev/autoPACK/autopack/scripts/testTiling.py")