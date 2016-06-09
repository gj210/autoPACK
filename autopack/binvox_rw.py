#  Copyright (C) 2012 Daniel Maturana
#  This file is part of binvox-rw-py.
#
#  binvox-rw-py is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  binvox-rw-py is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with binvox-rw-py. If not, see <http://www.gnu.org/licenses/>.
#

"""
Binvox to Numpy and back.

binvox argument
Usage: binvox [-d <voxel dimension>] [-t <voxel file type>] [-c] [-v] <model filespec>
  -license: show software license
  -d: specify voxel grid size (default 256, max 1024)(no max when using -e)
  -t: specify voxel file type (default binvox, also supported: hips, mira, vtk, raw, schematic, msh)
  -c: z-buffer based carving method only
  -dc: dilated carving, stop carving 1 voxel before intersection
  -v: z-buffer based parity voting method only (default is both -c and -v)
  -e: exact voxelization (any voxel with part of a triangle gets set)(does not use graphics card)
Additional parameters:
  -bb <minx> <miny> <minz> <maxx> <maxy> <maxz>: force a different input model bounding box
  -ri: remove internal voxels
  -cb: center model inside unit cube
  -rotx: rotate object 90 degrees ccw around x-axis before voxelizing
  -rotz: rotate object 90 degrees cw around z-axis before voxelizing
    both -rotx and -rotz can be used multiple times
  -aw: _also_ render the model in wireframe (helps with thin parts)
  -fit: only write the voxels in the voxel bounding box
  -bi <id>: when converting to schematic, use block ID <id>
  -mb: when converting using -e from .obj to schematic, parse block ID from material spec 'usemtl blockid_<id>' (ids 1-255 only)
  -down: downsample voxels by a factor of 2 in each dimension (can be used multiple times)
  -dmin <nr>: when downsampling, destination voxel is on if >= <nr> source voxels are (default 4)
Supported 3D model file formats:
  VRML V2.0: almost fully supported
  UG, OBJ, OFF, DXF, XGL, POV, BREP, PLY, JOT: only polygons supported
Example:
binvox -c -d 200 -t mira plane.wrl


>>> import numpy as np
>>> import binvox_rw
>>> with open('chair.binvox', 'rb') as f:
...     m1 = binvox_rw.read_as_3d_array(f)
...
>>> m1.dims
[32, 32, 32]
>>> m1.scale
41.133000000000003
>>> m1.translate
[0.0, 0.0, 0.0]
>>> with open('chair_out.binvox', 'wb') as f:
...     m1.write(f)
...
>>> with open('chair_out.binvox', 'rb') as f:
...     m2 = binvox_rw.read_as_3d_array(f)
...
>>> m1.dims==m2.dims
True
>>> m1.scale==m2.scale
True
>>> m1.translate==m2.translate
True
>>> np.all(m1.data==m2.data)
True

>>> with open('chair.binvox', 'rb') as f:
...     md = binvox_rw.read_as_3d_array(f)
...
>>> with open('chair.binvox', 'rb') as f:
...     ms = binvox_rw.read_as_coord_array(f)
...
>>> data_ds = binvox_rw.dense_to_sparse(md.data)
>>> data_sd = binvox_rw.sparse_to_dense(ms.data, 32)
>>> np.all(data_sd==md.data)
True
>>> # the ordering of elements returned by numpy.nonzero changes with axis
>>> # ordering, so to compare for equality we first lexically sort the voxels.
>>> np.all(ms.data[:, np.lexsort(ms.data)] == data_ds[:, np.lexsort(data_ds)])
True
"""

import numpy as np

class Voxels(object):
    """ Holds a binvox model.
    data is either a three-dimensional numpy boolean array (dense representation)
    or a two-dimensional numpy float array (coordinate representation).

    dims, translate and scale are the model metadata.

    dims are the voxel dimensions, e.g. [32, 32, 32] for a 32x32x32 model.

    scale and translate relate the voxels to the original model coordinates.

    To translate voxel coordinates i, j, k to original coordinates x, y, z:

    x_n = (i+.5)/dims[0]
    y_n = (j+.5)/dims[1]
    z_n = (k+.5)/dims[2]
    x = scale*x_n + translate[0]
    y = scale*y_n + translate[1]
    z = scale*z_n + translate[2]

    """

    def __init__(self, data, dims, translate, scale, axis_order):
        self.data = data
        self.dims = dims
        self.translate = translate
        self.scale = scale
        assert (axis_order in ('xzy', 'xyz'))
        self.axis_order = axis_order
        nx, ny, nz = dims
#        self.ijk = self.cartesian([range(nx), range(ny), range(nz)])
#        self.index = np.apply_along_axis(self.getIndex,1,self.ijk)
#        self.ordered_data_sparse = np.apply_along_axis(self.getIndexData,1,self.ijk)
#        self.data = self.data[self.index].reshape(self.dims)
        
    def clone(self):
        data = self.data.copy()
        dims = self.dims[:]
        translate = self.translate[:]
        return Voxels(data, dims, translate, self.scale, self.axis_order)

    def write(self, fp):
        write(self, fp)

    def ijkToxyz(self):
        #depends on axis_order
        sz = self.dims[0]*self.dims[0]*self.dims[0]
        if self.axis_order == "xyz":
            x_n = (self.data[0]+0.5)/self.dims[0]
            y_n = (self.data[1]+0.5)/self.dims[1]
            z_n = (self.data[2]+0.5)/self.dims[2]
            x = self.scale*x_n + self.translate[0]
            y = self.scale*y_n + self.translate[1]
            z = self.scale*z_n + self.translate[2]
        else :
            x_n = (self.data[0]+.5)/self.dims[0]
            y_n = (self.data[2]+.5)/self.dims[2]
            z_n = (self.data[1]+.5)/self.dims[1]
            x = self.scale*x_n + self.translate[0]
            y = self.scale*y_n + self.translate[2]
            z = self.scale*z_n + self.translate[1]           
        return np.vstack((x, y, z)).transpose()

    
    def xyzToijk(self, xyz):
        xyz = xyz.transpose()        
        x=xyz[0]
        y=xyz[1]
        z=xyz[2]
        x_n = (x-self.translate[0])/self.scale
        y_n = (y-self.translate[1])/self.scale
        z_n = (z-self.translate[2])/self.scale
        i = (x_n*self.dims[0])-0.5
        j = (y_n*self.dims[1])-0.5
        k = (z_n*self.dims[2])-0.5
        return np.vstack((i, j, k)).transpose()
        
    def ijkToIndex(self, ijk):
        ijk = ijk.transpose()
        wxh = self.dims[0]*self.dims[2]
        width = self.dims[0]
        index = ijk[0] * wxh + ijk[2] * width + ijk[1];  #// wxh = width * height = d * d
        return index
        
    def getIndex(self, ijk):
        wxh = self.dims[0]*self.dims[2]
        width = self.dims[0]
        index = ijk[0] * wxh + ijk[2] * width + ijk[1];  #// wxh = width * height = d * d
        return index

#    def getIJK(self, index):
#        wxh = self.dims[0]*self.dims[2]
#        width = self.dims[0]
#        index = xyz[0] * wxh + xyz[2] * width + xyz[1];  #// wxh = width * height = d * d
#        x=index / wxh
#        zp=index % wxh
#        z = zp / self.dims[0]
#        y = zp % self.dims[0]
#        x = nz_voxels / (dims[0]*dims[2])
#        zwpy = nz_voxels % (dims[0]*dims[2]) # z*w + y
#        z = zwpy / dims[0]
#        y = zwpy % dims[0]
#        return index

    def getIndexData(self, ijk):
        wxh = self.dims[0]*self.dims[2]
        width = self.dims[0]
        index = ijk[0] * wxh + ijk[2] * width + ijk[1];  #// wxh = width * height = d * d
        value = self.data[index]
        return index

    def cartesian(self, arrays, out=None):
            """
            #http://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
            Generate a cartesian product of input arrays.
        
            Parameters
            ----------
            arrays : list of array-like
                1-D arrays to form the cartesian product of.
            out : ndarray
                Array to place the cartesian product in.
        
            Returns
            -------
            out : ndarray
                2-D array of shape (M, len(arrays)) containing cartesian products
                formed of input arrays.
        
            Examples
            --------
            >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
            array([[1, 4, 6],
                   [1, 4, 7],
                   [1, 5, 6],
                   [1, 5, 7],
                   [2, 4, 6],
                   [2, 4, 7],
                   [2, 5, 6],
                   [2, 5, 7],
                   [3, 4, 6],
                   [3, 4, 7],
                   [3, 5, 6],
                   [3, 5, 7]])
        
            """
    
            arrays = [np.asarray(x) for x in arrays]
            dtype = arrays[0].dtype
    
            n = np.prod([x.size for x in arrays])
            if out is None:
                out = np.zeros([n, len(arrays)], dtype=dtype)
    
            m = n / arrays[0].size
            out[:, 0] = np.repeat(arrays[0], m)
            if arrays[1:]:
                self.cartesian(arrays[1:], out=out[0:m, 1:])
                for j in xrange(1, arrays[0].size):
                    out[j * m:(j + 1) * m, 1:] = out[0:m, 1:]
            return out         
            
def read_header(fp):
    """ Read binvox header. Mostly meant for internal use.
    """
    line = fp.readline().strip()
    if not line.startswith(b'#binvox'):
        raise IOError('Not a binvox file')
    dims = list(map(int, fp.readline().strip().split(b' ')[1:]))
    translate = list(map(float, fp.readline().strip().split(b' ')[1:]))
    scale = list(map(float, fp.readline().strip().split(b' ')[1:]))[0]
    line = fp.readline()
    return dims, translate, scale

def read(fp):
    """ Read binary binvox format as array.

    Returns the model with accompanying metadata.

    Voxels are stored in a three-dimensional numpy array, which is simple and
    direct, but may use a lot of memory for large models. (Storage requirements
    are 8*(d^3) bytes, where d is the dimensions of the binvox model. Numpy
    boolean arrays use a byte per element).

    Doesn't do any checks on input except for the '#binvox' line.
    """
    dims, translate, scale = read_header(fp)
    raw_data = np.frombuffer(fp.read(), dtype=np.uint8)
    #  The first byte of each pair is the value byte and is either 0 or 1 
    # (1 signifies the presence of a voxel). The second byte is the count byte 
    sz = np.prod(dims)
    data = np.empty(sz, dtype=np.bool)
    index, end_index = 0, 0
    i = 0
    while i < len(raw_data):
        value, count = map(int, (raw_data[i], raw_data[i+1]))
        end_index = index+count
        data[index:end_index] = value
        index = end_index
        i += 2
    # data = data.reshape(dims)
    return Voxels(data, dims, translate, scale, 'xzy'),raw_data
    
def read_as_3d_array(fp, fix_coords=True):
    """ Read binary binvox format as array.

    Returns the model with accompanying metadata.

    Voxels are stored in a three-dimensional numpy array, which is simple and
    direct, but may use a lot of memory for large models. (Storage requirements
    are 8*(d^3) bytes, where d is the dimensions of the binvox model. Numpy
    boolean arrays use a byte per element).

    Doesn't do any checks on input except for the '#binvox' line.
    """
    dims, translate, scale = read_header(fp)
    raw_data = np.frombuffer(fp.read(), dtype=np.uint8)
    # if just using reshape() on the raw data:
    # indexing the array as array[i,j,k], the indices map into the
    # coords as:
    # i -> x
    # j -> z
    # k -> y
    # if fix_coords is true, then data is rearranged so that
    # mapping is
    # i -> x
    # j -> y
    # k -> z
    values, counts = raw_data[::2], raw_data[1::2]
    data = np.repeat(values, counts).astype(np.bool)
    data = data.reshape(dims)
    if fix_coords:            
        data = np.transpose(data, (0, 2, 1))
        axis_order = 'xyz'
    else:
        axis_order = 'xzy'
    return Voxels(data, dims, translate, scale, axis_order)

def read_as_coord_array(fp, fix_coords=True):
    """ Read binary binvox format as coordinates.

    Returns binvox model with voxels in a "coordinate" representation, i.e.  an
    3 x N array where N is the number of nonzero voxels. Each column
    corresponds to a nonzero voxel and the 3 rows are the (x, z, y) coordinates
    of the voxel.  (The odd ordering is due to the way binvox format lays out
    data).  Note that coordinates refer to the binvox voxels, without any
    scaling or translation.

    Use this to save memory if your model is very sparse (mostly empty).

    Doesn't do any checks on input except for the '#binvox' line.
    """
    dims, translate, scale = read_header(fp)
    raw_data = np.frombuffer(fp.read(), dtype=np.uint8)

    values, counts = raw_data[::2], raw_data[1::2]

    sz = np.prod(dims)
    index, end_index = 0, 0
    end_indices = np.cumsum(counts)
    indices = np.concatenate(([0], end_indices[:-1])).astype(end_indices.dtype)

    values = values.astype(np.bool)
    indices = indices[values]
    end_indices = end_indices[values]

    nz_voxels = []
    for index, end_index in zip(indices, end_indices):
        nz_voxels.extend(range(index, end_index))
    nz_voxels = np.array(nz_voxels)
    # TODO are these dims correct?
    # according to docs,
    # index = x * wxh + z * width + y; // wxh = width * height = d * d

    x = nz_voxels / (dims[0]*dims[2])
    zwpy = nz_voxels % (dims[0]*dims[2]) # z*w + y
    z = zwpy / dims[0]
    y = zwpy % dims[0]
    if fix_coords:
        data = np.vstack((x, y, z))
        axis_order = 'xyz'
    else:
        data = np.vstack((x, z, y))
        axis_order = 'xzy'

    #return Voxels(data, dims, translate, scale, axis_order)
    return Voxels(np.ascontiguousarray(data), dims, translate, scale, axis_order)

def dense_to_sparse(voxel_data, dtype=np.int):
    """ From dense representation to sparse (coordinate) representation.
    No coordinate reordering.
    """
    if voxel_data.ndim!=3:
        raise ValueError('voxel_data is wrong shape; should be 3D array.')
    return np.asarray(np.nonzero(voxel_data), dtype)

def sparse_to_dense(voxel_data, dims, dtype=np.bool):
    if voxel_data.ndim!=2 or voxel_data.shape[0]!=3:
        raise ValueError('voxel_data is wrong shape; should be 3xN array.')
    if np.isscalar(dims):
        dims = [dims]*3
    dims = np.atleast_2d(dims).T
    # truncate to integers
    xyz = voxel_data.astype(np.int)
    # discard voxels that fall outside dims
    valid_ix = ~np.any((xyz < 0) | (xyz >= dims), 0)
    xyz = xyz[:,valid_ix]
    out = np.zeros(dims.flatten(), dtype=dtype)
    out[tuple(xyz)] = True
    return out

#def get_linear_index(x, y, z, dims):
    #""" Assuming xzy order. (y increasing fastest.
    #TODO ensure this is right when dims are not all same
    #"""
    #return x*(dims[1]*dims[2]) + z*dims[1] + y

def write(voxel_model, fp):
    """ Write binary binvox format.

    Note that when saving a model in sparse (coordinate) format, it is first
    converted to dense format.

    Doesn't check if the model is 'sane'.

    """
    if voxel_model.data.ndim==2:
        # TODO avoid conversion to dense
        dense_voxel_data = sparse_to_dense(voxel_model.data, voxel_model.dims)
    else:
        dense_voxel_data = voxel_model.data

    fp.write('#binvox 1\n')
    fp.write('dim '+' '.join(map(str, voxel_model.dims))+'\n')
    fp.write('translate '+' '.join(map(str, voxel_model.translate))+'\n')
    fp.write('scale '+str(voxel_model.scale)+'\n')
    fp.write('data\n')
    if not voxel_model.axis_order in ('xzy', 'xyz'):
        raise ValueError('Unsupported voxel model axis order')

    if voxel_model.axis_order=='xzy':
        voxels_flat = dense_voxel_data.flatten()
    elif voxel_model.axis_order=='xyz':
        voxels_flat = np.transpose(dense_voxel_data, (0, 2, 1)).flatten()

    # keep a sort of state machine for writing run length encoding
    state = voxels_flat[0]
    ctr = 0
    for c in voxels_flat:
        if c==state:
            ctr += 1
            # if ctr hits max, dump
            if ctr==255:
                fp.write(chr(state))
                fp.write(chr(ctr))
                ctr = 0
        else:
            # if switch state, dump
            fp.write(chr(state))
            fp.write(chr(ctr))
            state = c
            ctr = 1
    # flush out remainders
    if ctr > 0:
        fp.write(chr(state))
        fp.write(chr(ctr))

if __name__ == '__main__':
    import doctest
    doctest.testmod()

