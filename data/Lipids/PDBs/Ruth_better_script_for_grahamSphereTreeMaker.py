#!/usr/bin/env python

# Graham wrote this on Feb 26 2013
# It is used to create a sphereTree directly from a PDB file of a small molecule
# It is ideal for lipids where the clusterGUI algorithm is unhappy with sparse atom sets
# It simply takes every other atom and makes a spheretree with a radius of 1.5 for each atom position.

import glob, os

def saveSphereModel(fileRead,filename, minR=-1, maxR=-1 ):
    
    minR = 1.50 #self.radg #added by Graham 4/4/11
    maxR = 2.00 #self.radm #added by Graham 4/4/11
    r=3.0
    centers=(19,7,4)

    fw = open(filename, 'w')
    ptr = open(fileRead)
    ctr = 0
    lines = ptr.readlines()
    for l in lines:
        if l.find("ATOM")==0 or l.find("HETA")==0:
            ctr += 1
    fw.write("# rmin rmax\n")
    fw.write("%6.2f  %6.2f\n"%(minR, maxR))
    fw.write("\n")
    fw.write("# number of levels\n")
    fw.write("1\n")
    fw.write("\n")
    fw.write("# number of spheres in level 1\n")
    fw.write("%d\n"%int(ctr/2))  #len(ctr))
    fw.write("\n")
    fw.write("# x y z r of spheres in level 1\n")
    ctr2=0
    for l in lines:
        if l.find("ATOM")==0 or l.find("HETA")==0:
            xcoord = float(l[30:38])
            ycoord = float(l[38:46])
            zcoord = float(l[46:54])
            if ctr2%2==0:
                fw.write("%6.2f %6.2f %6.2f %6.2f\n"%(xcoord,ycoord,zcoord,r))
            ctr2 += 1
    fw.close()


filelist = glob.glob("*.pdb")
for f in filelist:
    ctr = 1
    name = os.path.splitext(os.path.basename(f))[0]
    outputfilename = name+ '.sph'
    saveSphereModel(f,outputfilename, -1, -1 )