# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 09:43:18 2017

@author: ludov
"""

#with given recipe open, show one ingredient sphereTree + PDB
import c4d
env = c4d.af.values()[0].histoVol.values()[0]
name = "mpn529" #HU
ingr =env.getIngrFromName(name)
pdb = ingr.source["pdb"]
# load in pmv ?
sphTree_coords = ingr.positions
sphTree_radii = ingr.radii
#build points clouds
import upy
helper = upy.getHelperClass()()
helper.PointCloudObject("sph"+name,vertices=sphTree_coords[0])