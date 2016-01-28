# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 17:15:12 2016

@author: ludo
"""
import sys
import numpy as np
from numpy import savez, dtype , fromfile 

#sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R17 Demo_BEED7D65/plugins/ePMV/mgl64/MGLToolsPckgs/")
#sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R17 Demo_BEED7D65/plugins/ePMV/mgl64/MGLToolsPckgs/PIL")

sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs")

from autopack import transformation as tr

#parse MolDB, then write again but convert all RH->LH using -X
#unity X pointing to the right, Y pointing up, and Z pointing away from the user, towards the screen 
#so invert X
filename="test_full_final_nctr.molb"
f=open(filename,"rb")


filename="test_full_final_nctr_fixed_tr.molb"
fout=open(filename,"wb")

#the header size to jump is 
#nmolecules instance ?
nframe=100
nr=0;
step_nr=0
mols=0
state_size=0
b = [0,0,0,0];
size = 8;
liste_m=[]
nr = 0
size = 8;
m = np.zeros((4,4))

ind = np.array([1,2,3,0])
mp = np.array([ 1, -1, -1,  1])

b=fromfile(f,'<i',count=4)
step_nr = fromfile(f,'<i',count=1)[0]
mols = fromfile(f,'<i',count=1)[0]#nb of instance total
state_size= fromfile(f,'<i',count=1)[0] #header ?
header_size=60
f.seek(state_size,1)
f.seek(size,1)            
f.seek(size*nr*4,1)

#p=[]
#q=[]
pos = fromfile(f,'<d',count=mols*4).reshape(mols,4)*np.array([-1,1,1,1])
rot = fromfile(f,'<d',count=mols*9).reshape(mols,3,3)

#maybe the reshape doesnt give the correct 3x3

qq=[tr.quaternion_from_matrix(r.transpose()).tolist() for r in rot]
quat=np.array(qq)[:,ind]*mp


xz
#euler=[tr.euler_from_matrix(r.transpose()) for r in rot]
#quat =[tr.quaternion_from_euler(e[0],-e[1],-e[2],axes='szxy') for e in euler]

#try different axis and use euler->quaternion
#quat =[MayaRotationToUnity(e) for e in euler]
#conversion
#p.extend(pos.tolist())
#q.extend(quat[:])

np.array(pos,'f').tofile(fout)
np.array(quat,'f').tofile(fout)
for i in range(1,nframe):
    f.seek(header_size,1)
    pos = fromfile(f,'<d',count=mols*4).reshape(mols,4)*np.array([-1,1,1,1])
    rot = fromfile(f,'<d',count=mols*9).reshape(mols,3,3)

    qq=[tr.quaternion_from_matrix(r.transpose()).tolist() for r in rot]
    quat=np.array(qq)[:,ind]*mp

#    euler=[tr.euler_from_matrix(r.transpose()) for r in rot]
#    quat =[tr.quaternion_from_euler(e[0],-e[1],-e[2],axes='szxy') for e in euler]
#    p.extend(pos.tolist())
#    q.extend(quat)
    np.array(pos,'f').flatten().tofile(fout)
    np.array(quat,'f').flatten().tofile(fout)
#np.array(p,'f').flatten().tofile(fout)
#np.array(q,'f').flatten().tofile(fout)   
fout.close()
f.close()

#alterntive to check
#ind = np.array([0,2,3,1])
#mp = np.array([ 1, -1, -1,  1])
#qq = [tr.quaternion_from_matrix(r).tolist() for r in rot]
#nqq=np.array(qq)[:,i]*m


#there is 4*mols double
#m[3][:3] = fromfile(f,'<d',count=4)[:3]#size*4 4 double ?
#f.seek( size*(mols-1-nr)*4, 1 )
#f.seek( size*nr*9, 1 )
#m[:3,:3] = fromfile(f,'<d',count=9).reshape(3,3)#.transpose()#size*4 4 double ?
#f.seek( size*(mols-1-nr)*9, 1)  
#liste_m.append(m)
#if (log) : print step_nr,m,mols
#f.close()



#public static Vector3 euler_from_matrix(Matrix4x4 M)
#    {
#        /*
#         * _AXES2TUPLE = {
#            'sxyz': (0, 0, 0, 0), 'sxyx': (0, 0, 1, 0), 'sxzy': (0, 1, 0, 0),
#            'sxzx': (0, 1, 1, 0), 'syzx': (1, 0, 0, 0), 'syzy': (1, 0, 1, 0),
#            'syxz': (1, 1, 0, 0), 'syxy': (1, 1, 1, 0), 'szxy': (2, 0, 0, 0),
#            'szxz': (2, 0, 1, 0), 'szyx': (2, 1, 0, 0), 'szyz': (2, 1, 1, 0),
#            'rzyx': (0, 0, 0, 1), 'rxyx': (0, 0, 1, 1), 'ryzx': (0, 1, 0, 1),
#            'rxzx': (0, 1, 1, 1), 'rxzy': (1, 0, 0, 1), 'ryzy': (1, 0, 1, 1),
#            'rzxy': (1, 1, 0, 1), 'ryxy': (1, 1, 1, 1), 'ryxz': (2, 0, 0, 1),
#            'rzxz': (2, 0, 1, 1), 'rxyz': (2, 1, 0, 1), 'rzyz': (2, 1, 1, 1)}
#            */
#        float _EPS = 8.8817841970012523e-16f;
#        Vector4 _NEXT_AXIS = new Vector4(1, 2, 0, 1);
#        Vector4 _AXES2TUPLE = new Vector4(0, 0, 0, 0);//sxyz
#        int firstaxis = (int)_AXES2TUPLE[0];
#        int parity = (int)_AXES2TUPLE[1];
#        int repetition = (int)_AXES2TUPLE[2];
#        int frame = (int)_AXES2TUPLE[3];
#
#        int i = firstaxis;
#        int j = (int)_NEXT_AXIS[i + parity];
#        int k = (int)_NEXT_AXIS[i - parity + 1];
#        float ax = 0.0f;
#        float ay = 0.0f;
#        float az = 0.0f;
#        //Matrix4x4 M = matrix;// = numpy.array(matrix, dtype=numpy.float64, copy=False)[:3, :3]
#        if (repetition == 1)
#        {
#            float sy = Mathf.Sqrt(M[i, j] * M[i, j] + M[i, k] * M[i, k]);
#            if (sy > _EPS)
#            {
#                ax = Mathf.Atan2(M[i, j], M[i, k]);
#                ay = Mathf.Atan2(sy, M[i, i]);
#                az = Mathf.Atan2(M[j, i], -M[k, i]);
#            }
#            else
#            {
#                ax = Mathf.Atan2(-M[j, k], M[j, j]);
#                ay = Mathf.Atan2(sy, M[i, i]);
#                az = 0.0f;
#            }
#        }
#        else
#        {
#            float cy = Mathf.Sqrt(M[i, i] * M[i, i] + M[j, i] * M[j, i]);
#            if (cy > _EPS)
#            {
#                ax = Mathf.Atan2(M[k, j], M[k, k]);
#                ay = Mathf.Atan2(-M[k, i], cy);
#                az = Mathf.Atan2(M[j, i], M[i, i]);
#            }
#            else
#            {
#                ax = Mathf.Atan2(-M[j, k], M[j, j]);
#                ay = Mathf.Atan2(-M[k, i], cy);
#                az = 0.0f;
#            }
#        }
#        if (parity == 1)
#        {
#            ax = -ax;
#            ay = -ay;
#            az = -az;
#        }
#        if (frame == 1)
#        {
#            ax = az;
#            az = ax;
#        }
#        //radians to degrees
#        //array([-124.69515353,  -54.90319877,  108.43494882])
#        //to -56.46 °110.427 °126.966 °
#        //Vector3 euler = new Vector3 (ay * Mathf.Rad2Deg, az * Mathf.Rad2Deg, -ax * Mathf.Rad2Deg);
#        Vector3 euler = new Vector3(ax * Mathf.Rad2Deg, ay * Mathf.Rad2Deg, az * Mathf.Rad2Deg);
#        //if (euler.x < 0) euler.x += 360.0f;
#        //if (euler.y < 0) euler.y += 360.0f;
#        //if (euler.z < 0) euler.z += 360.0f;
#        return euler;
#    }
#var mat = MyUtility.quaternion_matrix(rotation);
#            var euler = MyUtility.euler_from_matrix(mat);
#            rotation = MyUtility.MayaRotationToUnity(euler)

