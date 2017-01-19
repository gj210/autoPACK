# -*- coding: utf-8 -*-
"""
Created on Tue May 17 09:34:54 2016

@author: ludov
"""
import h5py
import json
import sys
import os
import math
import urllib2
import requests
from scipy import spatial

sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs\lib-tk")
import numpy as np
#import autopack
#import upy

#from upy import hostHelper
#autopack.helper = hostHelper.Helper()
#autopack.helper.host = "none"
R=35000.0
#helper = upy.getHelperClass()(vi="nogui")
#autopack.helper = helper
#pick objectrbc.dae
spacing = 60.0 #encapsulating radius of hemoglobin
bbsize = [70000,20000,70000]
bb=[[-R,-10000,-R],[R,10000,R]]

R1=[ 0, 7626.9, 10855.487, 14218.598, 17178.135,19599.575]
R2=[35000.0, 34262.739, 33186.543, 31706.774, 29688.908,16640.038]
H = 2000.0
TotalH = 20000
from autopack.Compartment import Compartment
##name, vertices, faces, vnormals
o = Compartment("rbc", None, None, None,filename="D:\\Data\\cellPACK_data\\cellPACK_database_1.1.0\\geometries\\rbc.dae")
o.OGsrfPtsBht = spatial.cKDTree(tuple(o.vertices), leafsize=10)               
from autopack.Grid import HaltonGrid, Grid
#bb=[[-R,-35,-R],[R,35,R]]
hg = HaltonGrid(boundingBox=bb, space=spacing, setup=False)
size = hg.getNBgridPoints()
hg.create3DPointLookup()
N=np.linalg.norm(hg.masterGridPositions,axis=1)
mask = N < R
inside_points = np.nonzero(mask)[0]
##print size
##print size1[0]*size[1]*size[2]
#hg.create3DPointLookup()
##hg.computeGridNumberOfPoint(bb, spacing)
#print hg.nbGridPoints#,len(hg.masterGridPositions)
##hg = HaltonGrid(boundingBox=bb, space=spacing, setup=False)
##hg.create3DPointLookup()
##keep only the one in circle
##which norm < R
#N=np.linalg.norm(hg.masterGridPositions,axis=1)
#mask = N < R
#inside_points = np.nonzero(mask)[0]
##hg.masterGridPositions[inside_points]
##8000000points
## D = o.OGsrfPtsBht.query_ball_point(hg.masterGridPositions, spacing*1.2)
##d = o.OGsrfPtsBht.query(hg.masterGridPositions)
##count = 0
##inside_points=[]
##prev_inside=False
##inside=False
##for p in hg.masterGridPositions:
##    inside = o.checkPointInside_rapid(p, hg.diag, ray=3)
##    if inside :
##        if d[0][count] > 30.0 :
##            inside_points.append(count)
###    if d[0][count] < spacing*2.0:#len(D[count]):
###        #check if inside
###        inside = o.checkPointInside_rapid(p, hg.diag, ray=3)
###        if inside :
###            if d[0][count] > 30.0 :
###                inside_points.append(count)
###        prev_inside=inside
###    else :
###        if prev_inside :
###            inside_points.append(count)            
##    count+=1
##
#N = len(inside_points)
#ones = np.ones(N)
#axes = np.random.random((N,3))
#quat = np.column_stack([axes,ones])
#pos = np.column_stack([hg.masterGridPositions[inside_points],ones])
#allpos = np.vstack([pos,quat])
###np.save("D:\\Data\\cellPACK_data\\cellPACK_database_1.1.0\\geometries\\rbc_points")
#fptr = open("D:\\Data\\cellPACK_data\\cellPACK_database_1.1.0\\geometries\\rbc_points", "wb")
#np.array(pos, 'f').flatten().tofile(fptr)  # 4float position
#np.array(quat, 'f').flatten().tofile(fptr)  # 4flaot quaternion
#fptr.close()
#print N
z
from autopack.Grid import HaltonGrid, Grid
bb=[[-850,-850,-850],[850,850,850]]
spacing = 50
hg = Grid(boundingBox=bb, space=spacing, setup=False)
#size = hg.getNBgridPoints()
#print size
#print size1[0]*size[1]*size[2]
hg.create3DPointLookup()
#hg.computeGridNumberOfPoint(bb, spacing)
print hg.nbGridPoints#,len(hg.masterGridPositions)

#check if C:\Users\ludov\AppData\Roaming\autoPACK\cache_geometries\HIV_VLP.binvox exist
#f='C:\Users\ludov\Downloads\HIV_VLP.binvox'
ff="C:\\Users\\ludov\\AppData\\Roaming\\autoPACK\\cache_geometries\\HIV_VLP.binvox"
from autopack import binvox_rw
with open(ff, 'rb') as f:
     model = binvox_rw.read_as_coord_array(f)
with open(ff, 'rb') as f:     
     m,r = binvox_rw.read(f)   

#indice of inside BB
#ijk = model.xyzToijk(xyz_Data[0])

m1 = (hg.masterGridPositions < bb[0]).any(axis=1)
m2 = (hg.masterGridPositions > bb[1]).any(axis=1)
m3 = m1 | m2
#outside indice
outsidebb = np.nonzero(m3)[0]
insidebb = np.nonzero(m3 == False)[0]

ijk = m.xyzToijk(hg.masterGridPositions[insidebb]).astype(int)
i = m.ijkToIndex(ijk).astype(int)
#inside
inbb_inside = np.nonzero(m.data[i] == True)[0]

inside_points = insidebb[inbb_inside]

hg.gridPtId = np.zeros(hg.gridVolume, 'i')
#hg.gridPtId

#import c4d;env=c4d.af.values()[0].histoVol.values()[0];self = env.compartments[0];v=env.grid.masterGridPositions[self.insidePoints];s = env.afviewer.vi.Points("test", vertices=v, );s = env.afviewer.vi.Points("test1", vertices=self.binvox_3d, )
#d,nb=env.grid.getClosestGridPoint(self.ogsurfacePoints);v2=env.grid.masterGridPositions[nb];s = env.afviewer.vi.Points("test3", vertices=v, )      
                
#from autopack import binvox_rw
#with open('C:\Users\ludov\Downloads\RBC_forLudo_12.binvox', 'rb') as f:
#     model = binvox_rw.read_as_coord_array(f)
##with open('C:\Users\ludov\Downloads\RBC_forLudo_1.binvox', 'rb') as f:
##     bin_model = binvox_rw.read_as_3d_array(f,fix_coords=False)
##with open('C:\Users\ludov\Downloads\RBC_forLudo_1.binvox', 'rb') as f:     
##     m,r = binvox_rw.read(f)
##model.translate=[0,0,0]
#model.axis_order="xzy"
##model.data = m.ijk.transpose()
#xyz_Data = model.ijkToxyz()
##ones = np.ones(N)
#
#print len(xyz_Data)
#fptr = open("D:\\Data\\cellPACK_data\\cellPACK_database_1.1.0\\geometries\\rbc_points", "wb")
#np.array(xyz_Data, 'f').flatten().tofile(fptr)  # 4float position
#fptr.close()

#np.array(quat, 'f').flatten().tofile(fptr)  # 4flaot quaternion

#halton inside mesh, with flood fill
#get the mesh
#build the grid
#process
#floodfill based on distance from surface
##first point is outside
# for all points in hg.masterGridPosition
#o.checkPointInside_rapid(point, diag, ray=1)