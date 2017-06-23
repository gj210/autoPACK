# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 13:18:40 2016

@author: ludov
"""


from collada import Collada
from collada import material
from collada import source
from collada import geometry
from collada import scene
import numpy

def oneMaterial(name,collada_xml,color=None):
    if color == None :
        color = numpy.random.random(4)
        color[3] = 1.0
    if len(color) == 3 :
        color = [color[0],color[1],color[2],1.0]
    effect = material.Effect("effect"+name, [], "phong", 
                                 diffuse=color)
    mat = material.Material("material"+name, name+"_material", effect)
    matnode = scene.MaterialNode("material"+name, mat, inputs=[])
    collada_xml.effects.append(effect)
    collada_xml.materials.append(mat)
    return matnode

def buildMeshGeom(name,vertices, faces,vnormals, collada_xml,matnode):
    iname = name
    vertxyz = numpy.array(vertices)# * numpy.array([1,1,-1])
    vert_src = source.FloatSource(iname+"_verts-array", vertxyz.flatten(), ('X', 'Y', 'Z'))
    geom = geometry.Geometry(collada_xml, "geometry"+iname, iname, [vert_src])
    input_list = source.InputList()
    input_list.addInput(0, 'VERTEX', "#"+iname+"_verts-array")
    fi=numpy.array(faces,int)#[:,::-1]
    triset = geom.createTriangleSet(fi.flatten(), input_list, iname+"materialref")
    geom.primitives.append(triset)
    collada_xml.geometries.append(geom)
    master_geomnode = scene.GeometryNode(geom, [matnode])
    master_node = scene.Node("node_"+iname, children=[master_geomnode,])#,transforms=[tr,rz,ry,rx,s])
    return master_node    
        
def coarseMolSurface(coords, radii, XYZd =[32,32,32],isovalue=6.0,resolution=-0.1,padding=0.0,
                         name='CoarseMolSurface',geom=None):
    from UTpackages.UTblur import blur
    if radii is None :
        radii = np.ones(len(coords))*1.8
    #volarr, origin, span = blur.generateBlurmap(coords, radii, XYZd,resolution, padding = 0.0)
    volarr, origin, span = blur.generateBlurmap(np.ascontiguousarray(coords).tolist(), radii.tolist(), XYZd,resolution, padding = 0.0)
    volarr.shape = (XYZd[0],XYZd[1],XYZd[2])
    volarr = np.ascontiguousarray(np.transpose(volarr), 'f')
    weights =  np.ones(len(radii),"f")
    h = {}
    from Volume.Grid3D import Grid3DF
    maskGrid = Grid3DF( volarr, origin, span , h)
    h['amin'], h['amax'],h['amean'],h['arms']= maskGrid.stats()
    from UTpackages.UTisocontour import isocontour
    isocontour.setVerboseLevel(0)
    data = maskGrid.data
    origin = np.array(maskGrid.origin).astype('f')
    stepsize = np.array(maskGrid.stepSize).astype('f')
    # add 1 dimension for time steps amd 1 for multiple variables
    if data.dtype.char!= np.float32:
#            print 'converting from ', data.dtype.char
        data = data.astype('f')#Numeric.Float32)
    newgrid3D = np.ascontiguousarray(np.reshape( np.transpose(data),
                                          (1, 1)+tuple(data.shape) ), data.dtype.char)
#        print "ok"       
    ndata = isocontour.newDatasetRegFloat3D(newgrid3D, origin, stepsize)
#        print "pfff"
    isoc = isocontour.getContour3d(ndata, 0, 0, isovalue,
                                       isocontour.NO_COLOR_VARIABLE)
    vert = np.zeros((isoc.nvert,3)).astype('f')
    norm = np.zeros((isoc.nvert,3)).astype('f')
    col = np.zeros((isoc.nvert)).astype('f')
    tri = np.zeros((isoc.ntri,3)).astype('i')
    isocontour.getContour3dData(isoc, vert, norm, col, tri, 0)
    #print vert
    if maskGrid.crystal:
        vert = maskGrid.crystal.toCartesian(vert)
    return vert, norm, tri        
    
collada_xml = Collada()
collada_xml.assetInfo.unitname="centimeter"
collada_xml.assetInfo.unitmeter=0.01
collada_xml.assetInfo.upaxis="Y_UP"

root_env=scene.Node(env.name)
myscene = scene.Scene(env.name+"_Scene", [root_env])
collada_xml.scenes.append(myscene)
collada_xml.scene = myscene                

name = "myMolecule"
matnode=oneMaterial(name,collada_xml)
master_node=buildMeshGeom(name,v,f,collada_xml,matnode)
collada_xml.nodes.append(master_node)

collada_xml.write("test.dae")
