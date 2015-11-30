from pandac.PandaModules import NodePath, GeomVertexFormat, GeomVertexWriter, \
  GeomVertexData, Geom, GeomTriangles, GeomNode

# ------------------------------------------------------------------------------
# GEOMETRY DATA GENERATORS
# ------------------------------------------------------------------------------
import math

class GeometryData:
  def __init__(self):
    pass
  
  def getVertex(self, i):
    return self.vertices[i]
  
  def getNormal(self, i):
    return self.normal[i]
  
  def getTriangle(self, i):
    return self.triangles[i]
  
  def getUv(self, i):
    return self.uv[i]

# a icosaeder but unfortunately not lying on the points
class IcosaederData(GeometryData):
  phi = (1.0 + math.sqrt(5.0)) / 2.0
  a = 1.0/2.0 / 0.587785243988
  b = 1.0/(2.0 * phi) / 0.587785243988
  vertices = [
      (-b,  a, 0),
      ( b,  a, 0),
      ( 0,  b, a),
      ( a,  0, b),
      ( 0, -b, a),
      (-a,  0, b),
      ( 0,  b,-a),
      (-a,  0,-b),
      ( 0, -b,-a),
      (-b, -a, 0),
      ( b, -a, 0),
      ( a,  0,-b),
    ]
  
  triangles = (
      ( 2,  1,  0),
      ( 2,  3,  1),
      ( 2,  5,  4),
      ( 2,  4,  3),
      ( 2,  0,  5),
      ( 8,  7,  6),
      ( 8, 10,  9),
      ( 8, 11, 10),
      ( 8,  9,  7),
      ( 8,  6, 11),
      ( 1, 11,  6),
      ( 6,  0,  1),
      ( 0,  6,  7),
      ( 7,  5,  0),
      ( 5,  7,  9),
      ( 9,  4,  5),
      ( 4,  9, 10),
      (10,  3,  4),
      ( 3, 10, 11),
      (11,  1,  3),
    )
  
  lines = {
      0: [ 1, 2, 5, 6, 7],
      1: [ 0, 2, 3, 6,11],
      2: [ 0, 1, 3, 4, 5],
      3: [ 1, 2, 4,10,11],
      4: [ 2, 3, 5, 9,10],
      5: [ 0, 2, 4, 7, 9],
      6: [ 0, 1, 7, 8,11],
      7: [ 0, 5, 6, 8, 9],
      8: [ 6, 7, 9,10,11],
      9: [ 4, 5, 7, 8,10],
     10: [ 3, 4, 8, 9,11],
     11: [ 1, 3, 6, 8,10]
    }
  
  def __init__(self, radius=1.0):
    GeometryData.__init__(self)
    self.radius = radius
  
  def getVertex(self, i):
    #v = Vec3(self.vertices[i][0] * self.radius, self.vertices[i][1] * self.radius, self.vertices[i][2] * self.radius)
    #print "ico", v.length()
    return [self.vertices[i][0] * self.radius, self.vertices[i][1] * self.radius, self.vertices[i][2] * self.radius]
  
  def getNormal(self, i):
    return self.vertices[i]
  
  def getUv(self, i):
    x,y,z = self.vertices[i]
    u = -((math.atan2(x,y)) / math.pi) / 2.0 + 0.5
    v = (z / self.b) / 2.0 + 0.5
    return [u,v]

class Cube2Data(GeometryData):
  triangles = (
      # front
      ( 0*3+1, 4*3+1, 5*3+1),
      ( 0*3+1, 5*3+1, 1*3+1),
      # back
      ( 6*3+1, 2*3+1, 3*3+1),
      ( 6*3+1, 3*3+1, 7*3+1),
      # left
      ( 2*3+0, 0*3+0, 1*3+0),
      ( 2*3+0, 1*3+0, 3*3+0),
      # right
      ( 4*3+0, 6*3+0, 7*3+0),
      ( 4*3+0, 7*3+0, 5*3+0),
      # top
      ( 1*3+2, 5*3+2, 7*3+2),
      ( 1*3+2, 7*3+2, 3*3+2),
      # bottom
      ( 2*3+2, 6*3+2, 4*3+2),
      ( 2*3+2, 4*3+2, 0*3+2)
    )
  uv = (
      [1,0], [0,0], [1,0], # 0
      [0,0], [0,1], [1,1], # 1
      [1,1], [1,0], [0,0], # 2
      [0,1], [0,0], [0,1], # 3
      [0,0], [1,0], [1,1], # 4
      [0,1], [1,1], [1,0], # 5
      [1,0], [1,1], [0,1], # 6
      [1,1], [0,1], [0,0], # 7
    )
  
  
  def __init__(self, width, height, depth):
    GeometryData.__init__(self)
    self.vertices = list()
    self.normal = list()
    for x in [0,width]:
      for y in [0,height]:
        for z in [0,depth]:
          self.vertices.append([x,y,z])
          self.vertices.append([x,y,z])
          self.vertices.append([x,y,z])
          
          nx = (x-width/2.)/(width/2.0)
          self.normal.append([nx,0,0])
          ny = (y-height/2.)/(height/2.0)
          self.normal.append([0,ny,0])
          nz = (z-depth/2.)/(depth/2.0)
          self.normal.append([0,0,nz])
  
  def getUv(self, i):
    return self.uv[i]


class SphereData(GeometryData):
  def __init__(self, radius, widthSegments=32, heightSegments=17, angle=360, inverted=False):
    GeometryData.__init__(self)
    self.radius = radius
    assert(heightSegments>=3)
    assert(widthSegments>=4)
    self.vertices = list()
    self.uv = list()
    for j in xrange(heightSegments-1):
      for i in xrange(widthSegments):
        thetaS = float(j)/(heightSegments-2)
        theta = thetaS * math.pi * (angle/360.)
        phiS = float(i)/(widthSegments-1 )
        phi  = phiS * math.pi * 2
        x = math.sin(theta) * math.cos(phi)
        y = math.cos(theta)
        z = -math.sin(theta) * math.sin(phi)
        self.vertices.append( (x,z,y) )
        self.uv.append( (thetaS, phiS))
    
    self.triangles = list()
    for j in xrange(heightSegments-2):
      for i in xrange(widthSegments-1):
        v0 = (j  )*widthSegments + i
        v1 = (j+1)*widthSegments + i+1
        v2 = (j  )*widthSegments + i+1
        if not inverted: self.triangles.append( (v2,v1,v0) )
        else:            self.triangles.append( (v0,v1,v2) )
        v0 = (j  )*widthSegments + i  
        v1 = (j+1)*widthSegments + i  
        v2 = (j+1)*widthSegments + i+1
        if not inverted: self.triangles.append( (v2,v1,v0) )
        else:            self.triangles.append( (v0,v1,v2) )
  
  def getVertex(self, i):
    return [self.vertices[i][0] * self.radius, self.vertices[i][1] * self.radius, self.vertices[i][2] * self.radius]
  def getUv(self, i):
    return self.uv[i]
  def getNormal(self, i):
    return self.vertices[i]

class TubeWallData(GeometryData):
  def __init__(self, radius=1.0, length=1.0, radiusSegments=32, inverted=False):
    GeometryData.__init__(self)
    assert(radiusSegments>=4)
    self.vertices = list()
    self.normal = list()
    self.uv = list()
    for i in xrange(radiusSegments):
      partial = float(i)/(radiusSegments-1)
      alpha = partial * math.pi *2
      xn = math.cos(alpha)
      x = xn * radius
      yn = math.sin(alpha)
      y = yn * radius
      #print "v", Vec3(x,y,0).length()
      #print alpha, x, y
      self.vertices.append( (x,y,0) )
      if not inverted:
        self.normal.append( (xn,yn,0) )
      else:
        self.normal.append( (-xn,-yn,0) )
      self.uv.append( (partial, 0) )
      self.vertices.append( (x,y,length) )
      if not inverted:
        self.normal.append( (xn,yn,0) )
      else:
        self.normal.append( (-xn,-yn,0) )
      self.uv.append( (partial, 1.) )
    
    self.triangles = list()
    if not inverted:
      for i in xrange(0,radiusSegments-1):
        self.triangles.append( [i*2+2,i*2+1,i*2  ] )
        self.triangles.append( [i*2+3,i*2+1,i*2+2] )
    else:
      for i in xrange(0,radiusSegments-1):
        self.triangles.append( [i*2  ,i*2+1,i*2+2] )
        self.triangles.append( [i*2+2,i*2+1,i*2+3] )
  
  def getUv(self, i):
    return self.uv[i]
  
  def getNormal(self, i):
    return self.normal[i]


class CircleData(GeometryData):
  def __init__(self, outerRadius=1.0, innerRadius=0.0, radiusSegments=32, inverted=False):
    GeometryData.__init__(self)
    self.inverted = inverted
    self.vertices = list()
    self.uv = list()
    for i in xrange(radiusSegments):
      alpha = float(i)/(radiusSegments-1) * math.pi *2
      x = math.cos(alpha)
      y = math.sin(alpha)
      self.vertices.append( (x*outerRadius,y*outerRadius,0) )
      self.uv.append([x/2.+.5,y/2.+.5])
      self.vertices.append( (x*innerRadius,y*innerRadius,0) )
      self.uv.append([x*innerRadius/outerRadius/2.+.5,y*innerRadius/outerRadius/2.+.5])
    
    self.triangles = list()
    if not inverted:
      for i in xrange(0,radiusSegments-1):
        self.triangles.append( [i*2+2,i*2+1,i*2  ] )
        self.triangles.append( [i*2+3,i*2+1,i*2+2] )
    else:
      for i in xrange(0,radiusSegments-1):
        self.triangles.append( [i*2  ,i*2+1,i*2+2] )
        self.triangles.append( [i*2+2,i*2+1,i*2+3] )
  
  def getUv(self, i):
    return self.uv[i]
  def getNormal(self, i):
    if not self.inverted:
      return [0,0,1]
    else:
      return [0,0,-1]


# ------------------------------------------------------------------------------
# DEFAULT GEOMETRY GENERATOR
# ------------------------------------------------------------------------------
class Geometry(NodePath):
  def __init__(self):
    NodePath.__init__(self, 'Geometry')
  
  def addGeometry(self, geomData):
    debugGui = dict()
    
    format = GeomVertexFormat.getV3n3t2()
    vdata = GeomVertexData('name', format, Geom.UHStatic)
    vertex = GeomVertexWriter(vdata, 'vertex')
    normal = GeomVertexWriter(vdata, 'normal')
    texcoord = GeomVertexWriter(vdata, 'texcoord')
    prim = GeomTriangles(Geom.UHStatic)
    
    postphonedTriangles = list()
    vtxTargetId0 = vtxTargetId1 = vtxTargetId2 = None
    vtxDataCounter = 0
    for vtxSourceId0, vtxSourceId1, vtxSourceId2 in geomData.triangles:
      vx0,vy0,vz0 = v0 = geomData.getVertex(vtxSourceId0)
      vx1,vy1,vz1 = v1 = geomData.getVertex(vtxSourceId1)
      vx2,vy2,vz2 = v2 = geomData.getVertex(vtxSourceId2)
      # prepare the vertices
      uvx0, uvy0 = uv0 = geomData.getUv(vtxSourceId0)
      uvx1, uvy1 = uv1 = geomData.getUv(vtxSourceId1)
      uvx2, uvy2 = uv2 = geomData.getUv(vtxSourceId2)
      #
      n0 = geomData.getNormal(vtxSourceId0)
      n1 = geomData.getNormal(vtxSourceId1)
      n2 = geomData.getNormal(vtxSourceId2)
      
      # make it wrap nicely
      if min(uvx0,uvx1,uvx2) < .25 and max(uvx0,uvx1,uvx2) > 0.75:
        if uvx0 < 0.25: uvx0 += 1.0
        if uvx1 < 0.25: uvx1 += 1.0
        if uvx2 < 0.25: uvx2 += 1.0
      
      vertex.addData3f(*v0)
      normal.addData3f(*n0)
      texcoord.addData2f(*uv0)
      vtxTargetId0 = vtxDataCounter
      vtxDataCounter += 1
    
      vertex.addData3f(*v1)
      normal.addData3f(*n1)
      texcoord.addData2f(*uv1)
      vtxTargetId1 = vtxDataCounter
      vtxDataCounter += 1
    
      vertex.addData3f(*v2)
      normal.addData3f(*n2)
      texcoord.addData2f(*uv2)
      vtxTargetId2 = vtxDataCounter
      vtxDataCounter += 1
      
      prim.addVertex(vtxTargetId0)
      prim.addVertex(vtxTargetId1)
      prim.addVertex(vtxTargetId2)
      prim.closePrimitive()
      
      if False:
        if vtxSourceId0 not in debugGui:
          i = InfoTextBillaboarded(render)
          i.setScale(0.05)
          i.billboardNodePath.setPos(Vec3(x0,y0,z0)*1.1)
          i.setText('%i: %.1f %.1f %.1f\n%.1f %.1f' % (vtxSourceId0, x0,y0,z0, nx0, ny0))
          debugGui[vtxSourceId0] = i
        if vtxSourceId1 not in debugGui:
          i = InfoTextBillaboarded(render)
          i.setScale(0.05)
          i.billboardNodePath.setPos(Vec3(x1,y1,z1)*1.1)
          i.setText('%i: %.1f %.1f %.1f\n%.1f %.1f' % (vtxSourceId1, x1,y1,z1, nx1, ny1))
          debugGui[vtxSourceId1] = i
        if vtxSourceId2 not in debugGui:
          i = InfoTextBillaboarded(render)
          i.setScale(0.05)
          i.billboardNodePath.setPos(Vec3(x2,y2,z2)*1.1)
          i.setText('%i: %.1f %.1f %.1f\n%.1f %.1f' % (vtxSourceId2, x2,y2,z2, nx2, ny2))
          debugGui[vtxSourceId2] = i
    
    geom = Geom(vdata)
    geom.addPrimitive(prim)
    
    node = GeomNode('gnode')
    node.addGeom(geom)
    
    nodePath = self.attachNewNode(node)
    return nodePath


# ------------------------------------------------------------------------------
# CUSTOM GEOMETRY GENERATORS
# ------------------------------------------------------------------------------
class Cube(Geometry):
  def __init__(self, height, width, depth):
    Geometry.__init__(self)
    self.addGeometry(Cube2Data(height, width, depth))

class Icosaeder(Geometry):
  def __init__(self, radius=1.0):
    Geometry.__init__(self)
    self.addGeometry(IcosaederData(radius))

class Cylinder(Geometry):
  def __init__(self, radius, length, radiusSegments):
    Geometry.__init__(self)
    side   = self.addGeometry(TubeWallData(radius=radius, length=length, radiusSegments=radiusSegments))
    bottom = self.addGeometry(CircleData(outerRadius=radius, innerRadius=0.0, radiusSegments=radiusSegments,inverted=True))
    top    = self.addGeometry(CircleData(outerRadius=radius, innerRadius=0.0, radiusSegments=radiusSegments))
    top.setZ(length)

class Tube(Geometry):
  def __init__(self, outerRadius=1.0, innerRadius=0.5, length=1.0, radiusSegments=32):
    Geometry.__init__(self)
    cylinderTop       = self.addGeometry(CircleData(outerRadius=outerRadius, innerRadius=innerRadius, radiusSegments=radiusSegments))
    cylinderTop.setZ(length)
    cylinderBottom    = self.addGeometry(CircleData(outerRadius=outerRadius, innerRadius=innerRadius, radiusSegments=radiusSegments,inverted=True))
    cylinderInnerWall = self.addGeometry(TubeWallData(radius=innerRadius, length=length, radiusSegments=radiusSegments,inverted=True))
    cylinderOuterWall = self.addGeometry(TubeWallData(radius=outerRadius, length=length, radiusSegments=radiusSegments))

class Capsule(Geometry):
  def __init__(self, radius=1.0, length=1.0, radiusSegments=32, heightSegments=16):
    Geometry.__init__(self)
    side = self.addGeometry(TubeWallData(radius=radius, length=length, radiusSegments=radiusSegments, inverted=False))
    topSphere = self.addGeometry(SphereData(radius=radius, widthSegments=radiusSegments, heightSegments=heightSegments, angle=180))
    topSphere.setZ(length)
    bottomSphere = self.addGeometry(SphereData(radius=radius, widthSegments=radiusSegments, heightSegments=heightSegments, angle=180))
    bottomSphere.setR(180)
    bottomSphere.setH(180)

class Pyramid(Geometry):
  def __init__(self, radius=1.0, sides=4):
    Geometry.__init__(self)
    topSphere = self.addGeometry(SphereData(radius=radius, widthSegments=sides+1, heightSegments=3, angle=180))
    bottom    = self.addGeometry(CircleData(outerRadius=radius, innerRadius=0.0, radiusSegments=sides+1, inverted=True))

class Sphere(Geometry):
  def __init__(self, radius=.5, segements=16):
    Geometry.__init__(self)
    sphere = self.addGeometry(SphereData(radius=radius, widthSegments=segements, heightSegments=segements, angle=360))

