import upy
import numpy
import scipy
# http://docs.scipy.org/doc/scipy/reference/spatial.distance.html
from scipy.spatial import distance
helper = upy.getHelperClass()()

parent = helper.getObject("parent")
childs_pos = numpy.array([helper.ToVec(helper.getTranslation(ch)) for ch in helper.getChilds(parent)])
pairwise_dist = distance.pdist(childs_pos)
#this is a condensed distance matrix, to get the full matrice
all_dist = distance.squareform(pairwise_dist)
#look at the min,max, avg distance for each row/proteins then analyze the radial distribution function
#http://www.shocksolution.com/microfluidics-and-biotechnology/calculating-the-pair-correlation-function-in-python/
# apply the code from https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation
