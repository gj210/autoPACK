# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 10:32:41 2012

@author: LUDOVIC AUTIN
#execfile("/Users/ludo/DEV/Autofill_svn_test/autofill/trunk/AutoFillClean/testGamer.py")
"""
import gamer
import numpy as np
from time import time
import upy
helperClass = upy.getHelperClass()
helper = helperClass()

def timeFunction(function,args,kw):
    """
    Mesure the time for performing the provided function. 

    @type  function: function
    @param function: the function to execute
    @type  args: liste
    @param args: the liste of arguments for the function
    

    @rtype:   list/array
    @return:  the center of mass of the coordinates
    """
     
    t1=time()
    if len(kw):
        res = function(*args,**kw)
    else :
        res = function(*args)
    print(("time "+function.__name__, time()-t1))
    return res


def gamer_to_host(gmesh, boundaries, mesh_name="gamer_improved",
                  switch_layer=True,create_new_mesh=True):
    #myprint("gamer_to_host ",gmesh)
    # Check arguments
#    if not isinstance(gmesh, gamer.SurfaceMesh):
#        self.drawError(errormsg="expected a SurfaceMesh") 
    
    # Get scene
    scn = helper.getCurrentScene()
#    self.waitingCursor(1)
    
    verts = [(gvert.x, gvert.y, gvert.z) for gvert in gmesh.vertices()]
    un_selected_vertices = [i for i, gvert in enumerate(gmesh.vertices())
                            if not gvert.sel]
    faces = [(gface.a, gface.b, gface.c) for gface in gmesh.faces()]
    markers = [(i, gface.m) for i, gface in enumerate(gmesh.faces()) \
               if gface.m != -1]

    # If we create a new mesh we copy the boundaries to a dict
#    if self.getVal(self.gparams["create_new_mesh"]):
#        new_boundaries = {}
#        for boundary_name in boundaries.keys():
#            boundary = boundaries[boundary_name]
#            new_boundaries[boundary_name] = dict(
#                marker=boundary["marker"], r=boundary["r"], g=boundary["g"], \
#                b=boundary["b"], faces={})
#        
#        # Do not copy the faces information grab that from the gamer mesh
#        boundaries = new_boundaries
    
    # Create marker to boundary map
#    face_markers = {}
#    for boundary in boundaries.values():
#        face_markers[boundary["marker"]] = []
#
#    # Gather all faces of a marker
#    for face, marker in markers:
#        if marker in face_markers:
#            face_markers[marker].append(face)
#
#    # Set the faces of the corresponding boundary
#    for boundary in boundaries.values():
#        self._set_boundary_faces(boundary, face_markers[boundary["marker"]])
#    
#    # Ensure editmode is off
#    editmode = self.helper.toggleEditMode()

    if create_new_mesh:

        # Create new mesh
        #createsNmesh(self,name,vertices,vnormals,faces,color=[1,0,0],
        #                material=None,smooth=True,proxyCol=False, **kw)
        obj, bmesh = helper.createsNmesh(mesh_name, verts, None, \
                                              faces, smooth=0)
#        # If not generating a totally new mesh
#        switch_to_layer = scn.getLayers()[-1]
#        if switch_layer:
#            # Switch to another layer
#            switch_to_layer += 1
#            switch_to_layer = 1 if switch_to_layer > 20 else switch_to_layer
#
#        scn.setLayers([switch_to_layer])
#        obj.layers = [switch_to_layer]
#        
#        self.helper.addObjectToScene(scn, obj)
#        self.helper.ObjectsSelection([obj])
#
#        # Set the property dictionary
#        # FIXME: Is this safe? Is boundaries always something I can use?
#        self.helper.setProperty(obj, "boundaries", boundaries)
    
    else:
        # Get selected mesh
        # Update present mesh
        obj = helper.getObject(mesh_name)
        helper.updateMesh(obj, vertices=verts, faces=faces)
    
#    #myprint("un_selected_vertices ", len(un_selected_vertices))
#    self.helper.selectVertices(obj, un_selected_vertices, False)
#
#    # Restore editmode
#    self.helper.restoreEditMode(editmode)
#
#    # Repaint boundaries if there were markers in the GAMer data
#    if markers:
#        self._repaint_boundaries(obj)
#    
#    self.waitingCursor(0)
#    self.updateViewer()


def host_to_gamer(obj=None, check_for_vertex_selection=True,compId = 1):
    "Transfer the active mesh to a GAMer surface mesh"
    # Take the first one
    #myprint("host_to_gamer ", obj)
    # Get selected mesh
#    if obj is None:
#        obj = self._get_selected_mesh()
    #myprint("_get_selected_mesh ", obj)
#    if obj is None:
#        return None, None

#    self.waitingCursor(1)
    
    # Grab vertices and Faces
    vertices, selected_vertices = helper.getMeshVertices(obj, selected=True)
    vertices = helper.getMeshVertices(obj)
    faces = helper.getMeshFaces(obj)
    translation = helper.ToVec(helper.getTranslation(obj))
    # Init gamer mesh
    gmesh = gamer.SurfaceMesh(len(vertices), len(faces))
    def setVert(co, gvert, sel):
        gvert.x = co[0] + translation[0]
        gvert.y = co[1] + translation[1]
        gvert.z = co[2] + translation[2]
        gvert.sel = sel
    
    selected_vertices = np.ones(len(vertices), dtype=bool)
    
    [setVert(*args) for args in zip(vertices, gmesh.vertices(), selected_vertices)]

    # Transfere data from blender mesh to gamer mesh
    for face, gface in zip(faces, gmesh.faces()):
        gface.a, gface.b, gface.c = face
        gface.m = compId#comp Id 
    return gmesh

def makeFaces(cells):
    # faces are counter clockwise in the file
    faces = {} # dict used to avoid duplicates
    allFaceList = [] # list of triangular faces
    allTetFaces = [] # list of face indices for faces of each tet
                  # negative index means reverse indices to make clockwise
    nbf = 0
    nbt = 0
    n = 0
    allTetFaces = []
    for i,j,k,l in cells:
        allTetFaces.extend([ [i,k,j], [k,l,j], [l,k,i] , [l,i,j] ])
    return allTetFaces

#get an object, pass it to gamer
#build tetrahedre
obj = helper.getObject("Cube")
gm1 = host_to_gamer(obj=obj,compId = 1)
gm1.use_volume_contraint = True
gm1.volume_contraint = 0.01
gm1.marker = 1
gm1.as_hole = False

obj2 = helper.getObject("Sphere")
gm2 = host_to_gamer(obj=obj2,compId = 2)
gm2.use_volume_contraint = True
gm2.volume_contraint =0.01
gm2.marker = 2
gm2.as_hole = False

gem_mesh = gamer.GemMesh([gm1, gm2])
#res = timeFunction(function,args,kw)
#gem_mesh.vertex
#gem_mesh.vertices
#gem_mesh.num_vertices
#gem_mesh.num_cells
#gem_mesh.cell fa,fb,fc,fd,grp,id,mat,na,nb,nc,nd
#gem_mesh.cells
def makeObject(gem_mesh):
    verts = [(gvert.x, gvert.y, gvert.z) for gvert in gem_mesh.vertices()]
    cids =  [gvert.chrt for gvert in gem_mesh.vertices()]
    cells = [(gface.na, gface.nb, gface.nc,gface.nd) for gface in gem_mesh.cells()]
    markers = [(i, gface.mat) for i, gface in enumerate(gem_mesh.cells())]
    faces = makeFaces(cells)
    
    obj, bmesh = helper.PointCloudObject("tetraPoint", vertices=verts)
    obj, bmesh = helper.createsNmesh("tetra", verts, None,faces, smooth=0)

def spheres(gem_mesh):
    p = helper.getObject("debug")
    if p is None :
        p = helper.newEmpty("debug")
    meshsphere = helper.getObject("base_sphere")
    verts1 = [(gvert.x, gvert.y, gvert.z) for gvert in gem_mesh.vertices() if gvert.chrt == 1]
    verts2 = [(gvert.x, gvert.y, gvert.z) for gvert in gem_mesh.vertices() if gvert.chrt == 0]    
    result = helper.instancesSphere("marker1",verts1,[2.0],meshsphere,[[1,0,0]],None,parent=p)
    result = helper.instancesSphere("marker0",verts2,[2.0],meshsphere,[[0,0,1]],None,parent=p)
   