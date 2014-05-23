# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 21:58:54 2014

@author: ludo
"""
import sys
import json
import numpy
import math
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14_95696CEA/plugins/ePMV/mgl64/MGLToolsPckgs/")
# axis sequences for Euler angles
# epsilon for testing whether a number is close to zero
_EPS = numpy.finfo(float).eps * 4.0

_NEXT_AXIS = [1, 2, 0, 1]
# map axes strings to/from tuples of inner axis, parity, repetition, frame
_AXES2TUPLE = {
    'sxyz': (0, 0, 0, 0), 'sxyx': (0, 0, 1, 0), 'sxzy': (0, 1, 0, 0),
    'sxzx': (0, 1, 1, 0), 'syzx': (1, 0, 0, 0), 'syzy': (1, 0, 1, 0),
    'syxz': (1, 1, 0, 0), 'syxy': (1, 1, 1, 0), 'szxy': (2, 0, 0, 0),
    'szxz': (2, 0, 1, 0), 'szyx': (2, 1, 0, 0), 'szyz': (2, 1, 1, 0),
    'rzyx': (0, 0, 0, 1), 'rxyx': (0, 0, 1, 1), 'ryzx': (0, 1, 0, 1),
    'rxzx': (0, 1, 1, 1), 'rxzy': (1, 0, 0, 1), 'ryzy': (1, 0, 1, 1),
    'rzxy': (1, 1, 0, 1), 'ryxy': (1, 1, 1, 1), 'ryxz': (2, 0, 0, 1),
    'rzxz': (2, 0, 1, 1), 'rxyz': (2, 1, 0, 1), 'rzyz': (2, 1, 1, 1)}

_TUPLE2AXES = dict((v, k) for k, v in _AXES2TUPLE.items())

def euler_from_matrix(matrix, axes='sxyz'):
    """Return Euler angles from rotation matrix for specified axis sequence.

    axes : One of 24 axis sequences as string or encoded tuple

    Note that many Euler angle triplets can describe one matrix.

    >>> R0 = euler_matrix(1, 2, 3, 'syxz')
    >>> al, be, ga = euler_from_matrix(R0, 'syxz')
    >>> R1 = euler_matrix(al, be, ga, 'syxz')
    >>> numpy.allclose(R0, R1)
    True
    >>> angles = (4*math.pi) * (numpy.random.random(3) - 0.5)
    >>> for axes in _AXES2TUPLE.keys():
    ...    R0 = euler_matrix(axes=axes, *angles)
    ...    R1 = euler_matrix(axes=axes, *euler_from_matrix(R0, axes))
    ...    if not numpy.allclose(R0, R1): print(axes, "failed")

    """
    try:
        firstaxis, parity, repetition, frame = _AXES2TUPLE[axes.lower()]
    except (AttributeError, KeyError):
        _TUPLE2AXES[axes]  # validation
        firstaxis, parity, repetition, frame = axes

    i = firstaxis
    j = _NEXT_AXIS[i+parity]
    k = _NEXT_AXIS[i-parity+1]

    M = numpy.array(matrix, dtype=numpy.float64, copy=False)[:3, :3]
    if repetition:
        sy = math.sqrt(M[i, j]*M[i, j] + M[i, k]*M[i, k])
        if sy > _EPS:
            ax = math.atan2( M[i, j],  M[i, k])
            ay = math.atan2( sy,       M[i, i])
            az = math.atan2( M[j, i], -M[k, i])
        else:
            ax = math.atan2(-M[j, k],  M[j, j])
            ay = math.atan2( sy,       M[i, i])
            az = 0.0
    else:
        cy = math.sqrt(M[i, i]*M[i, i] + M[j, i]*M[j, i])
        if cy > _EPS:
            ax = math.atan2( M[k, j],  M[k, k])
            ay = math.atan2(-M[k, i],  cy)
            az = math.atan2( M[j, i],  M[i, i])
        else:
            ax = math.atan2(-M[j, k],  M[j, j])
            ay = math.atan2(-M[k, i],  cy)
            az = 0.0

    if parity:
        ax, ay, az = -ax, -ay, -az
    if frame:
        ax, az = az, ax
    return [math.degrees(ax), math.degrees(ay), math.degrees(az)]
    
def eulerToMatrix(euler): 
      # Assuming the angles are in radians.
    heading=euler[0]
    attitude=euler[1]
    bank=euler[2]
    m=[[ 1.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.],
       [ 0.,  0.,  1.,  0.],
       [ 0.,  0.,  0.,  1.]]
    ch = math.cos(heading)
    sh = math.sin(heading)
    ca = math.cos(attitude)
    sa = math.sin(attitude)
    cb = math.cos(bank)
    sb = math.sin(bank)
    m[0][0] = ch * ca
    m[0][1] = sh*sb - ch*sa*cb
    m[0][2] = ch*sa*sb + sh*cb
    m[1][0] = sa
    m[1][1] = ca*cb
    m[1][2] = -ca*sb
    m[2][0] = -sh*ca
    m[2][1] = sh*sa*cb + ch*sb
    m[2][2] = -sh*sa*sb + ch*cb
    return m
  
def getOneIngr(line):
    elem = line.split("<")
    pos = eval(elem[1][:-2])
    rot = eval(elem[2][:-2])
    rad = eval(elem[3][:-2])
    ingrname = elem[4][:-2]
    ingrcompNum = eval(elem[5][:-2])
    ptInd = eval(elem[6].split(">")[0])
    return pos,rot, ingrname,ingrcompNum,ptInd,rad

def toTEM():
    f="/Users/ludo/DEV/autopack_git/data/HIV/results/HIVresult_1_4.apr"
    rfile=open(f,"r")
    result=[]
    orgaresult=[[],]*2
    lines = rfile.readlines()   
    p=[]
    rr=[]     
    for l in lines :
        if not len(l) or len(l) < 6 : continue
        pos,rot, ingrname,ingrcompNum,ptInd,rad = getOneIngr(l)
        #should I multiply here
        r = numpy.array(rot).reshape(4,4)#numpy.matrix(mry90)*numpy.matrix(numpy.array(rot).reshape(4,4))
        if ingrcompNum == 0 :
            result.append([numpy.array(pos),numpy.array(r),ingrname,ingrcompNum,ptInd])
        else :
            orgaresult[abs(ingrcompNum)-1].append([numpy.array(pos),numpy.array(r),ingrname,ingrcompNum,ptInd])
        if ingrname == "iSUTM":
            p.append(numpy.array(pos)*0.05)
            rr.append(euler_from_matrix(r))
    fpdb="iSUTM_gp120_0.0.pdb"
    output="iSutm_coordinate.txt"
    aStr="# File created for TEM-simulator, version 1.3.\n"
    aStr+=str(len(p))+" 6\n"
    aStr+="#            x             y             z           phi         theta           psi\n"
    for i in range(len(p)):
        aStr+="{0:14.4f}{1:14.4f}{2:14.4f}{3:14.4f}{4:14.4f}{5:14.4f}\n".format(p[i][0],p[i][1],p[i][2],rr[i][0],rr[i][1],rr[i][2])
    f=open(output,"w")
    f.write(aStr)

def parseTEM(filename):
    f=open(filename,"r")
    regexp = r"[-+]?[0-9]*\.?[0-9]+"
    res=numpy.fromregex(filename, regexp,[('num', numpy.float), ])
    res2=res['num'][3:].reshape((int(res['num'][1]),6))
    pos,rot = numpy.hsplit(res2, 2)
    M = []#[eulerToMatrix(r) for r in rot]
    for i in range(len(pos)):
        mrot = eulerToMatrix(rot[i])
        mrot[3][:3] = pos[i]
        M.append(mrot)
    return M

def prepareBD_BOX():
    #should start with 2DSpheres experiments
    #then some simple molecules
    #two files:
    #one prm and str file which is the autoPACK ouput
    #only sphere and position. so we need to actually compute the spheres position foreach ingredient
    #this should go in the save method of the env and do the show sphereTree function.
    #sub ATM 491 -245.8562 -65.9199 -287.9397 30.0000 0.0000 60.0000 0.5922 1
    #sub name id x y z  Q 2R LJ m
    #then bond lines - each line de
    #ne a single bond and its parameters (Equation 1)
    #bond id(1) id(2) ro rmax H
    return

#import upy
#helper = upy.getHelperClass()()
#M=parseTEM("/Users/ludo/DEV/autopack_git/data/HIV/PDBsOld/Chimera_Dev/autoPACK_Chimera_HIV_1_3/iSutm_coordinate.txt")
#toinst=helper.getObject("iSUTM")
##left-hand / right hand problem ?
#ipoly = helper.instancePolygon("instOfObj", matrices=M, mesh=toinst)
#toTEM()
#execfile("/Users/ludo/DEV/tem_test.py")
def completeMapping(h,traj):
    r =  h.exteriorRecipe
    if r :
        for ingr in r.ingredients:
            if not len(ingr.results) : continue
            if ingr.Type == "SingleSphere" or ingr.Type == "SingleCube":
                traj.makeIngrMapping(ingr.name,len(ingr.results))
            elif ingr.Type == "MultiSphere":
                level = ingr.maxLevel
                nbPrim = len(ingr.radii[level])
                traj.makeIngrMapping(ingr.name,len(ingr.results)*nbPrim)
    #compartment ingr
    for orga in h.compartments:
        #compartment surface ingr
        rs =  orga.surfaceRecipe
        if rs :
            for ingr in rs.ingredients:
                if not len(inr.results) : continue
                if ingr.Type == "SingleSphere" or ingr.Type == "SingleCube":
                    traj.makeIngrMapping(ingr.name,len(ingr.results))
                elif ingr.Type == "MultiSphere":
                    level = ingr.maxLevel
                    nbPrim = len(ingr.radii[level])
                    traj.makeIngrMapping(ingr.name,len(ingr.results)*nbPrim)
        #compartment matrix ingr
        ri =  orga.innerRecipe
        if ri :
            for ingr in ri.ingredients:
                if not len(inr.results) : continue
                if ingr.Type == "SingleSphere" or ingr.Type == "SingleCube":
                    traj.makeIngrMapping(ingr.name,len(ingr.results))
                elif ingr.Type == "MultiSphere":
                    level = ingr.maxLevel
                    nbPrim = len(ingr.radii[level])
                    traj.makeIngrMapping(ingr.name,len(ingr.results)*nbPrim)

import autopack
def getIngredientInstancePos(traj, name,instance_id,frame):
        indice = traj.mapping[name][1][instance_id]
        pos = traj.getPosAt(indice,frame)
        return pos
        
def applyState(traj, h, frame):
        r =  h.exteriorRecipe
        indice = 0
        if r :
            for ingr in r.ingredients:
                for k in range(len(ingr.results)):
                    pos = getIngredientInstancePos(traj,ingr.name,k,frame)
                    autopack.helper.setTranslation(ingr.ipoly[k],pos)
                    indice+=1
                #k+=1
#        k=0
        #compartment ingr
        for orga in h.compartments:
            #compartment surface ingr
            rs =  orga.surfaceRecipe
            if rs :
                for ingr in rs.ingredients:
                    for k in range(len(ingr.results)):
                        pos = getIngredientInstancePos(traj,ingr.name,k,frame)
                        autopack.helper.setTranslation(ingr.ipoly[k],pos)
                        indice+=1
            #compartment matrix ingr
            ri =  orga.innerRecipe
            if ri :
                k=0
                for ingr in ri.ingredients:
                    for k in range(len(ingr.results)):
                        pos = getIngredientInstancePos(traj,ingr.name,k,frame)
                        autopack.helper.setTranslation(ingr.ipoly[k],pos)
                        indice+=1

def applyState_primitive_name(traj, h, frame):
        r =  h.exteriorRecipe
        indice = 0
        if r :
            for ingr in r.ingredients:
                print ingr.name
                if ingr.Type == "SingleSphere" or ingr.Type == "SingleCube":
                    for k in range(len(ingr.results)):
                        pos = traj.getIngredientInstancePos(ingr.name,k,frame)
                        autopack.helper.setTranslation(autopack.helper.getName(ingr.ipoly[k]),pos)
                        indice+=1
                elif ingr.Type == "MultiSphere":
                    level = ingr.maxLevel
                    nbPrim = len(ingr.radii[level])
                    for k in range(len(ingr.results)*nbPrim):
                        pos = traj.getIngredientInstancePos(ingr.name,k,frame)
                        print k
                        autopack.helper.setTranslation(autopack.helper.getName(ingr.isph[k]),pos)
                        indice+=1                
        for orga in h.compartments:
            #compartment surface ingr
            rs =  orga.surfaceRecipe
            if rs :
                for ingr in rs.ingredients:
                    if ingr.Type == "SingleSphere" or ingr.Type == "SingleCube":
                        for k in range(len(ingr.results)):
                            pos = traj.getIngredientInstancePos(ingr.name,k,frame)
                            autopack.helper.setTranslation(autopack.helper.getName(ingr.ipoly[k]),pos)
                            indice+=1
                    elif ingr.Type == "MultiSphere":
                        level = ingr.maxLevel
                        nbPrim = len(ingr.radii[level])
                        for k in range(len(ingr.results)*nbPrim):
                            pos = traj.getIngredientInstancePos(ingr.name,k,frame)
                            autopack.helper.setTranslation(autopack.helper.getName(ingr.isph[k]),pos)
                            indice+=1
            #compartment matrix ingr
            ri =  orga.innerRecipe
            if ri :
                k=0
                for ingr in ri.ingredients:
                    if ingr.Type == "SingleSphere" or ingr.Type == "SingleCube":
                        for k in range(len(ingr.results)):
                            pos = traj.getIngredientInstancePos(ingr.name,k,frame)
                            autopack.helper.setTranslation(autopack.helper.getName(ingr.ipoly[k]),pos)
                            indice+=1
                    elif ingr.Type == "MultiSphere":
                        level = ingr.maxLevel
                        nbPrim = len(ingr.radii[level])
                        for k in range(len(ingr.results)*nbPrim):
                            pos = traj.getIngredientInstancePos(ingr.name,k,frame)
                            autopack.helper.setTranslation(autopack.helper.getName(ingr.isph[k]),pos)
                            indice+=1