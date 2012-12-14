# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 15:05:07 2012

@author: -
"""
import sys
sys.path.append("/Users/ludo/DEV/Autofill_svn_test/autofill/trunk")
sys.path.append("/Users/ludo/DEV/MGLTOOLS/mgl32/MGLToolsPckgs/")
sys.path.append("/Users/ludo/DEV/MGLTOOLS/mgl32/MGLToolsPckgs/PIL")
import numpy
from AutoFill.Organelle import Organelle
o1 = Organelle("vesicle_Mesh",None, None, None,
               filename="/Users/ludo/DEV/Autofill_svn_test/autofill/trunk/AutoFill/cache_organelles/vesicle")
print "ok o1"
inner, surface = o1.getSurfaceInnerPointsSDF(o1.bb,20.0)
#print inner
# 
#from UTpackages.UTsdf import utsdf
#dim = 64
#dim1=dim+1
#size = dim1*dim1*dim1
# #can be 16,32,64,128,256,512,1024
#verts = numpy.array(o1.vertices)
#tris = numpy.array(o1.faces,"int")
#utsdf.setParameters(dim,0,1,[0,0,0,0,0,0])#size, bool isNormalFlip, bool insideZero,bufferArr
#datap = utsdf.computeSDF(verts,tris)

#try :
#    datap = utsdf.computeSDF(verts,tris)
#except :
#    print "what"
