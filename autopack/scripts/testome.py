from aicsimage import io#, processing
#img = processing.aicsImage.AICSImage("C:\\Users\\ludov\\Downloads\\AICS-12_881_7.ome.tif")
img = io.OmeTifReader("C:\\Users\\ludov\\Downloads\\AICS-12_881_7.ome.tif")
data = img.load()

import numpy
import scipy
import scipy.ndimage

def expand(datazyx, scalezyx):
    mindim = min(scalezyx[0], min(scalezyx[1], scalezyx[2]))
    factors = scalezyx/mindim
    # factors = maxdim / scalexyz
    print(factors)
    return scipy.ndimage.zoom(datazyx, factors, None, 1)


# Z - Y - X
# mydata = numpy.ndarray(shape=(3,2,2), dtype=float, buffer=numpy.array([[[1.0,1.0],[1.0,1.0]], [[0.33,0.33],[0.33,0.33]], [[0.0,0.0],[0.0,0.0]]]))
#img.data is T C Z Y X
print (expand(img.data[0][1], numpy.array([4, 1, 1])))


#note onoptimising the mesh in meshlab:
#1-clean->remove duplicate Faces
#2-clean->remove duplicate Vertices
#3-clean->remove unreferenced Vertices
#4-remeshing->Screened Poisson reconstruction
#5-smoothing->Taubin smoothing
#6-simplification->clustering decimation
