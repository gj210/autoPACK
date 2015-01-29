# import sys
# sys.path.append('C:\Users\Kevin\AppData\Roaming\MAXON\CINEMA 4D R15 Student_0471F205\plugins\ePMV\mgl64\MGLToolsPckgs\AutoFill')

import c4d
from math import sqrt, ceil
# import AutoFill
import autopack
from autopack import ingr_ui as UI
# helper = AutoFill.helper
helper = None
if helper is None:
    import upy
    helperClass = upy.getHelperClass()
    helper = helperClass()
    # AutoFill.helper = helper
import numpy
import heapq
"""
1. Build grid using radius of gyration of smallest guy we want to pack.
   - SphereTreeMaker (Rmin is radius of gyration, Rmax is the encapsulating radius).
2. Walk through triangles (that make up A,B,C) that make up the containers, assign triangles
   to closest grid points at the resolution of the grid.
    for triangle in polyhedron:
        if lengthOfAnySideLongerThanGridSpace:
            subdivide
                # Break any edge longer into segments max size of unit
                # Marching cubes?
                # Coarse projetion of our polyhedron onto our grid (at resolution of grid)


Plot computation time vs number of vertices/polygons in the original polyhedron

"""

class gridPoint:
    def __init__(self,i,globalC,isPolyhedron):
        self.index = int(i)
        self.isOutside = None
        self.minDistance = 99999 # Only store a number here if within certain distance from polyhedron
        self.representsPolyhedron = isPolyhedron
        self.closeFaces = []
        self.closestFaceIndex = 0
        self.testedEndpoint = None
        self.allDistances = [] # Stores a tuple list of distances to all points. (point,distance) = (5,2.5)
        self.globalCoord = numpy.array(globalC) # Stores the global coordinate associated with this point

def visualizeAtomArray(name = 'Test', vertexCoordinates = [(0,0,0)], sphereRadius = 5):
    """
    Given coordinates, visualize them in Cinema 4D using an Atom Array Object.
    """
    atom = c4d.BaseObject(c4d.Oatomarray)
    atom.SetName(name)
    atom[1000] = 3 #radius cylinder
    atom[1001] = sphereRadius #radius sphere
    atom[1002] = 5 #subdivision
    helper.addObjectToScene(helper.getCurrentScene(),atom,parent=None)
    
    parent = atom
    coords = vertexCoordinates
    nface = 0
    obj= c4d.PolygonObject(len(coords), nface)
    obj.SetName(name + '_child')
    cd4vertices = map(helper.FromVec,coords)
    map(obj.SetPoint,range(len(coords)),cd4vertices)
    helper.addObjectToScene(helper.getCurrentScene(),obj,parent=parent)

def setupBBOX(object_target):
    """
    Sets up a bounding box that encompasses a selected object in Cinema 4D. Returns the bounding box
    object.
    """
    mesh = object_target

    faces,vertices,vnormals,fn = helper.DecomposeMesh(mesh,edit=False,copy=False,tri=True,transform=True,fn=True)
    
    from autopack.Compartment import Compartment
    o1 = Compartment(helper.getName(mesh),vertices, faces, vnormals,fnormals=fn)
    o1.ref_obj = mesh
    o1.number = 1

    b = helper.Box("BBOX", cornerPoints=o1.bb)
    bb=o1.bb
    return o1

def makeGrid(lowerLeft,upperRight,spacing, inclusive = True):
    """
    Given lower left back corner, an upper right front corner, and a grid spacing, construct a list of
    coordinates that representsly an evenly spaced grid spanning the entire volume. Returned coordinates
    preserve the right/lefthandedness of the input coordinates, making this function work equally well
    in both 3D coordinate systems.
    """
    pointCoords = []
    xSize,ySize,zSize = upperRight[0] - lowerLeft[0], upperRight[1] - lowerLeft[1], upperRight[2] - lowerLeft[2]
    if inclusive:
        nx,ny,nz = xSize / spacing + 1, ySize / spacing + 1, zSize / spacing + 1
    else:
        nx,ny,nz = xSize / spacing, ySize / spacing, zSize / spacing
    nx,ny,nz = round(nx),round(ny),round(nz)
    x = numpy.linspace(lowerLeft[0],upperRight[0],nx)
    y = numpy.linspace(lowerLeft[1],upperRight[1],ny)
    z = numpy.linspace(lowerLeft[2],upperRight[2],nz)
    for a in x:
        for b in y:
            for c in z:
                pointCoords.append((a,b,c))
    # Order of above had to change to correctly produce a viable BBOX
    return pointCoords, (nx,ny,nz)

def getPointFrom3D(pt3d, gridSpacing, gridPtsPerEdge, BBOX):
    """
    Starts from orthogonal grid. Has a lower left, upper right coordinate from bounding box.
    Given lower left and known spacing of grid, we can provide any vector in 3D space and it will
    algebraically return index of closest grid point.

    Given a RIGHT-HANDED fine coordinate, it will return the index of the closest point in the LEFT-HANDED coarse grid.
    """
    x, y, z = pt3d  # Continuous 3D point to be discretized
    spacing1 = 1./gridSpacing  # Grid spacing = diagonal of the voxel determined by smalled packing radius
    NZ, NY, NX = gridPtsPerEdge  # vector = [length, height, depth] of grid, units = gridPoints
    
    # Added to fix some weird edge behavior
    #NZ, NY, NX = round(NZ),round(NY),round(NX)
    #print(NX,NY,NZ)
    OZ, OY, OX = BBOX.getBoundingBox()[0] # origin of Pack grid. OX and OZ were switched to compensate for coordinate systems
    #print(OX,OY,OZ)
    
    # Algebra gives nearest gridPoint ID to pt3D
    i = min( NX-1, max( 0, round((x-OX)*spacing1)))
    j = min( NY-1, max( 0, round((y-OY)*spacing1)))
    k = min( NZ-1, max( 0, round((z-OZ)*spacing1)))
    result = k*NX*NY + j*NX + i
    return int(result)
    #return int(k*NX*NY + j*NX + i)

def makeMarchingCube(gridSpacing,r):
    """
    Create a numpy array that represents the precomputed distances to each point
    for the cube of points surrounding our center point.
    """
    def _pythagorean(*edgeLengths):
        from math import sqrt
        total = 0
        for length in edgeLengths:
            total += length * length
        distance = sqrt(total)
        return distance
    from math import ceil
    pointsForRadius = ceil(r/gridSpacing) # Number of grid points required to represent our radius, rounded up
    pointsInEdge = 2 * pointsForRadius + 1 # Number of points in one edge of our cube

    center = pointsForRadius # The index if our center point
    cube = numpy.zeros(shape=(pointsInEdge,pointsInEdge,pointsInEdge))
    distX = numpy.zeros(shape=(pointsInEdge,pointsInEdge,pointsInEdge))
    distY = numpy.zeros(shape=(pointsInEdge,pointsInEdge,pointsInEdge))
    distZ = numpy.zeros(shape=(pointsInEdge,pointsInEdge,pointsInEdge))
    for a in range(pointsInEdge):
        lenX = a - center
        for b in range(pointsInEdge):
            lenY = b - center
            for c in range(pointsInEdge):
                lenZ = c - center
                cube[a][b][c] = _pythagorean(lenX,lenY,lenZ) * gridSpacing
                distX[a][b][c] = lenX
                distY[a][b][c] = lenY
                distZ[a][b][c] = lenZ
    return cube, distX, distY, distZ

def f_dot_product(vector1,vector2):
    """
    Return the dot product of two 3D vectors.
    """
    dottedVectors = [vector1[i] * vector2[i] for i in range(len(vector1))]
    return sum(dottedVectors)

def vcross(v1,v2):
    """
    Returns the cross-product of two 3D vectors.
    """
    x1,y1,z1 = v1
    x2,y2,z2 = v2
    return (y1*z2-y2*z1, z1*x2-z2*x1, x1*y2-x2*y1)

def vlen(v1):
    """
    Returns the length of a 3D vector.
    """
    x1,y1,z1 = v1
    return sqrt(x1*x1 + y1*y1 + z1*z1)

def f_ray_intersect_polyhedron(pRayStartPos, pRayEndPos, faces,vertices, pTruncateToSegment):
    """This function returns TRUE if a ray intersects a triangle.
    It also calculates and returns the UV coordinates of said colision as part of the intersection test,
    Makes sure that we are working with arrays

    * TAKES IN pRayStartPos as LEFT-HANDED COORDINATE (Z,Y,X)
    * TAKES IN pRayEndPos AS LEFT-HANDED COORDINATE (Z,Y,X)
    * faces is a list of vertex indices that defines the polyhedron
    * vertices is the list of global coordinates that faces refers to
    * pTruncateToSegment decides whether or not the segment terminates at the end position, or keeps on
      going forever
    """
    pRayStartPos = numpy.array(pRayStartPos)
    pRayEndPos = numpy.array(pRayEndPos)

    vLineSlope = pRayEndPos - pRayStartPos #This line segment defines an infinite line to test for intersection
    vPolyhedronPos = numpy.array((0,0,0))
    vTriPoints = vertices
    vTriPolys = faces
    vBackface = None
    vBackfaceFinal = None
    maxDistance = 99999
    vHitCount = 0
    
    # Globalize the polyhedron
    # Present in original coffee script, but unnecessary because DecomposeMesh gives global coords
    #for i in range(len(vTriPoints)):
    #   vTriPoints[i] = vTriPoints[i] + vPolyhedronPos

    vEpsilon = 0.00001
    vBreakj = False
    vCollidePos = None
    
    # Walk through each polygon in a polyhedron
    for testingFace in vTriPolys:
        #  Loop through all the polygons in an input polyhedron
        #vQuadrangle = 1  #  Default says polygon is a quadrangle.
        #vLoopLimit = 2  #  Default k will loop through polygon assuming its a quad.
        #if (vTriPolys[j+3] == vTriPolys[j+2])  #  Test to see if quad is actually just a triangle.
        #   {
        vQuadrangle = 0  #  Current polygon is not a quad, it's a triangle.
        vLoopLimit = 1  #  Set k loop to only cycle one time.
        
        for k in range(vLoopLimit):
            vTriPt0 = numpy.array(vTriPoints[testingFace[0]])  # Always get the first point of a quad/tri
            vTriPt1 = numpy.array(vTriPoints[testingFace[1]])  # Get point 1 for a tri and a quad's first pass, but skip for a quad's second pass
            vTriPt2 = numpy.array(vTriPoints[testingFace[2]])  # Get point 2 for a tri and a quad's first pass, but get point 3 only for a quad on its second pass.
    
            vE1 = vTriPt1 - vTriPt0  #  Get the first edge as a vector.
            vE2 = vTriPt2 - vTriPt0  #  Get the second edge.
            h = vcross(vLineSlope, vE2)
            
            a = f_dot_product(vE1,h)  #  Get the projection of h onto vE1.
            if a > -vEpsilon and a < vEpsilon:
                continue# If the ray is parallel to the plane then it does not intersect it, i.e, a = 0 +/- given rounding slope.
                #  If the polygon is a quadrangle, test the other triangle that comprises it.

            F = 1.0/a       
            s = pRayStartPos - vTriPt0  #  Get the vector from the origin of the triangle to the ray's origin.
            u = F * f_dot_product(s,h)
            if u < 0.0 or u > 1.0:
                continue
                #/* Break if its outside of the triangle, but try the other triangle if in a quad.
                #U is described as u = : start of vE1 = 0.0,  to the end of vE1 = 1.0 as a percentage.  
                #If the value of the U coordinate is outside the range of values inside the triangle,
                #then the ray has intersected the plane outside the triangle.*/

            q = vcross(s, vE1)
            v = F * f_dot_product(vLineSlope,q)
            if v <0.0 or u+v > 1.0:
                continue  
                #/*  Break if outside of the triangles v range.
                #If the value of the V coordinate is outside the range of values inside the triangle,
                #then the ray has intersected the plane outside the triangle.
                #U + V cannot exceed 1.0 or the point is not in the triangle.           
                #If you imagine the triangle as half a square this makes sense.  U=1 V=1 would be  in the 
                #lower left hand corner which would be in the second triangle making up the square.*/

            vCollidePos = vTriPt0 + u*vE1 + v*vE2  #  This is the global collision position.
            assert len(vCollidePos) == 3

            #  The ray is hitting a triangle, now test to see if its a triangle hit by the ray.
            vBackface = False
            if f_dot_product(vLineSlope, vCollidePos - pRayStartPos) > 0:  #  This truncates our infinite line to a ray pointing from start THROUGH end positions.
                vHitCount += 1
                #print(testingFace)
                if pTruncateToSegment and vlen(vLineSlope) < vlen(vCollidePos - pRayStartPos):
                    print('broken')
                    break # This truncates our ray to a line segment from start to end positions.

                d = squaredTwoPointDistance(pRayStartPos,vCollidePos)
                if d >= maxDistance:
                    continue
                if a < vEpsilon:  #  Test to see if the triangle hit is a backface.
                    #set master grid to organelle->getname inside
                    vBackface = True
                    #  This stuff is specific to our Point inside goals.
                    vBreakj = True  #  To see if a point is inside, I can stop at the first backface hit.
                    #break
                else:
                    vBreakj = True
                vBackfaceFinal = vBackface
                maxDistance = d
    #vBackfaceFinal = vBackfaces[distances.index(min(distances))]
    return vHitCount,vBackfaceFinal

def findPointsCenter(*args):
    # Average down the column, such that we're averaging across all measurements in one dimension
    center = numpy.mean(args[0], axis = 0)
    return center

def squaredTwoPointDistance(point1,point2):
    """Computes and returns the distance between two points in 3D space"""
    dists = [numpy.float64(point1[i] - point2[i]) for i in range(len(point1))]
    distsSquare = [x**2 for x in dists]
    dist = sum(distsSquare)
    return dist

def projectPolyhedronToGrid2(gridSpacing, radius = None):
    """
    Takes a polyhedron, and builds a grid. In this grid:
        - Projects the polyhedron to the grid.
        - Determines which points are inside/outside the polyhedron
        - Determines point's distance to the polyhedron.
    Usage:
    Select desired object, run this program with desired grid spacing resolution. By default,
    radius is 2 * gridSpacing. superFine provides the option doing a super leakproof test when
    determining which points are inside or outside. This usually not necessary, because the
    built-in algorithm has no known leakage cases. It is simply there as a safeguard.
    """
    from time import time
    if radius == None:
        radius = gridSpacing

    startTime = time()
    from math import ceil
    # Get the object, and decompose it
    object_target = helper.getCurrentSelection()[0]
    faces,verticesLH,vnormals,faceNormals = helper.DecomposeMesh(object_target,edit=False,copy=False,tri=True,transform=True,fn=True)
    # Setup a bounding box for our polyhedron
    b = setupBBOX(object_target)
    corners = b.getBoundingBox()
    
    # Create grid points to fill in our bbox with grid points at regular distances
    points, gridPtsPerEdge = makeGrid(corners[0],corners[1],gridSpacing) # These are left-handed

    # Create a gridPoint object for every single point we have in our grid. This pre-allocates all memory, so that future
    # computations are faster.
    gridPoints = []
    i = 0
    for point in points:
        gridPoints.append(gridPoint(i,point,isPolyhedron = False))
        i += 1
    assert len(gridPoints) == len(points)

    # Make a precomputed distance cube, where 
    distanceCube,distX,distY,distZ = makeMarchingCube(gridSpacing,radius)
    # Flatten and combine these arrays. This makes it easier to iterate over.
    distanceCubeF,distXF,distYF,distZF = distanceCube.flatten(),distX.flatten(),distY.flatten(),distZ.flatten()
    zippedNumbers = zip(distanceCubeF,distXF,distYF,distZF)

    NZ,NY,NX = gridPtsPerEdge
    OZ, OY, OX = b.getBoundingBox()[0]
    spacing1 = 1./gridSpacing # Inverse of the spacing. We compute this here, so we don't have to recompute it repeatedly
    allCoordinates = [] # Tracker for all the fine coordiantes that we have interpolated for the faces of the polyhedron
    pointsToTestInsideOutside = set()
    # Walk through the faces, projecting each to the grid and marking immediate neighbors so we can test said 
    # neighbors for inside/outside later.
    for face in faces:
        # Dist is a small function that calculates the distance between two points.
        # dist = lambda x,y: vlen([y[i] - x[i] for i in range(len(x))])
        # Get the vertex coordinates and conver to numpy arrays...
        triCoords = [numpy.array(verticesLH[i]) for i in face]
        thisFaceFineCoords = list(triCoords)
        allCoordinates.extend(triCoords)
        # We will use these u/v vectors to interpolate points that reside on the face
        pos = triCoords[0]
        u = triCoords[1] - pos
        v = triCoords[2] - pos
        # Smetimes the hypotenuse isn't fully represented. To remedy this, we will use an additional w vector
        # to interpolate points on the hypotenuse
        w = triCoords[2] - triCoords[1]

        # If either u or v is greater than the grid spacing, then we need to subdivide it
        # We will use ceil: if we have a u of length 16, and grid spacing of 5, then we want
        # a u at 0, 5, 10, 15 which is [0, 1, 2, 3] * gridSpacing. We need ceil to make this happen.
        
        # This is a much more efficient solution than the previous method of artificially increasing
        # the densith of the polygon mesh by 4, then mapping each of those small points to the grid.
        # However, it is also leaks more, for (as of yet) unknown reasons.
        
        # It was discovered that by using the default gridspacing, some faces will produce leakage. One solution is to
        # temporarily use a somewhat denser gridspacing to interpolate, and then project back onto our original spacing.
        # We'll decrease the gridspacing by 25% (so that it's 75% of the original). This seems to fix the issue without
        # any appreciable increase in computation time.
        gridSpacingTempFine = gridSpacing / 3
        # Determine the number of grid spacing-sized points we can fit on each vector.
        # Minimum is one because range(1) gives us [0]
        uSubunits, vSubunits, wSubunits = 1, 1, 1
        if vlen(u) > gridSpacingTempFine:
            uSubunits = ceil(vlen(u)/gridSpacingTempFine) + 1
        if vlen(v) > gridSpacingTempFine:
            vSubunits = ceil(vlen(v)/gridSpacingTempFine) + 1
        if vlen(w) > gridSpacingTempFine:
            wSubunits = ceil(vlen(w)/gridSpacingTempFine) + 1
        # Because we have observed leakage, maybe we want to try trying a denser interpolation, using numpy's linspace?
        for uSub in range(uSubunits):
            percentU = uSub * gridSpacingTempFine / vlen(u)
            percentU = min(percentU, 1.0) # Make sure that we have not stepped outside of our original u vector
            # h represents the height of the hypotenuse at this u. Naturally, we cannot go past the hypotenuse, so this will be
            # our upper bound.
            h = percentU * u + (1 - percentU) * v
            for vSub in range(vSubunits):
                percentV = vSub * gridSpacingTempFine / vlen(v)
                percentV = min(percentV, 1.0) # Make sure that we have not stepped oustide of our original v vector.
                interpolatedPoint = percentU * u + percentV * v
                # The original if: statement asks if the distance from the origin to the interpolated point is less than
                # the distance from the origin to the hypotenuse point, as such:
                # if vlen(interpolatedPoint) < vlen(h):
                # Wouldn't it be a better idea to measure distance to the u position instead? This is implemented below.
                if (vlen(interpolatedPoint - percentU * u) < vlen(h - percentU * u)):
                    allCoordinates.append(interpolatedPoint + pos)
                    thisFaceFineCoords.append(interpolatedPoint + pos)
                else:
                    break
        # Because the above only interpolates the face, and may not completely capture the hypotenuse, let's separately
        # interpolate points on the hypotenuse (the w vector). This way we can be 100% sure that our final solution is 
        # completely leak-proof
        for wSub in range(wSubunits):
            # Apply the same proceudre we did above for u/v, just for w (for hypotenuse interpolation)
            percentW = wSub * gridSpacingTempFine / vlen(w)
            percentW = min(percentW, 1.0)
            interpolatedPoint = percentW * w
            allCoordinates.append(interpolatedPoint + triCoords[1])
            thisFaceFineCoords.append(interpolatedPoint + triCoords[1])
        projectedIndices = set()
        # Once we have interpolated the face, let's project each interpolated point to the grid.
        for coord in thisFaceFineCoords:
            projectedPointIndex = getPointFrom3D(coord[::-1],gridSpacing,gridPtsPerEdge,b)
            projectedIndices.add(projectedPointIndex)

        # Walk through each grid point that our face spans, gather its closest neighbors, and annotate the neighbors with
        # minimum distance and closest faces, and flag them for testing inside/outside later.
        for P in list(projectedIndices):
            # if P == 53144:
            #     print('We have projected to the point that causes issues!!!!!!!!!!!')
            # Get the point object corresponding to the index, and set its polyhedron attribute to true
            g = gridPoints[P]
            g.representsPolyhedron = True
            # Get the coordinates of the point, and convert them to grid units
            zTemp,yTemp,xTemp = g.globalCoord
            i,j,k = round((xTemp-OX)*spacing1), round((yTemp-OY)*spacing1), round((zTemp-OZ)*spacing1)
            # Let's step through our distance cube, and assign faces/closest distances to each
            for d,x,y,z in zippedNumbers:
                # Get the grid indices for the point we're considering, and pass if we're stepping oustide the boundaries
                newI, newJ, newK = i + x, j + y, k + z
                if newI < 0 or newI > (NX-1) or newJ < 0 or newJ > (NY-1) or newK < 0 or newK > (NZ-1):
                    continue
                # Get the point index that this coordinate corresponds to.
                desiredPointIndex = int(round(newK*NX*NY + newJ*NX + newI))
                # if desiredPointIndex == 54199:
                #     print('HIT with face ' + str(face) + ' and polygon point ' + str(P))
                desiredPoint = gridPoints[desiredPointIndex]
                if desiredPoint.representsPolyhedron == True:
                    continue
                # Add the current face to the its list of closest faces
                if face not in desiredPoint.closeFaces:
                    desiredPoint.closeFaces.append(face)
                # Add the distance to the point's list of distances, and overwrite minimum distance if appropriate
                desiredPoint.allDistances.append((v,d))
                if d < desiredPoint.minDistance:
                    desiredPoint.minDistance = d
                # Later down the road, we want to test as few points as possible for inside/outside. Therefore,
                # we will only test points that are 
                # if abs(x) <= 1 and abs(y) <= 1 and abs(z) <= 1:
                #     pointsToTestInsideOutside.add(desiredPointIndex)
    # visualizeAtomArray('fineCoords',allCoordinates,5)
    projectedIndices = [x.index for x in gridPoints if x.representsPolyhedron == True]
    # projectedIndices = set()
    # for coord in allCoordinates:
    #     projectedPointIndex = getPointFrom3D(coord[::-1],gridSpacing,gridPtsPerEdge,b)
    #     projectedIndices.add(projectedPointIndex)
    # visualizeAtomArray('projectedCoords',[points[i] for i in projectedIndices],5)
    timeFinishProjection = time()
    print('Projecting polyhedron to grid took ' + str(timeFinishProjection - startTime) + ' seconds.')
    
    # Let's start filling in inside outside. Here's the general algorithm:
    # Walk through all the points in our grid. Once we encounter a point that has closest faces, 
    # then we know we need to test it for inside/outside. Once we test that for inside/outside, we
    # fill in all previous points with that same inside outisde property. To account for the possible
    # situation that there is a surface that is only partially bound by the bbox, then we need to
    # reset the insideOutsideTracker every time we have a change in more than 1 of the 3 coordinates
    # because that indicates we're starting a new row/column of points.

    # This tracks our inside/outside.
    isOutsideTracker = None
    # This tracks the points that we've iterated over which we do not know if inside/outside. Remember
    # to reset every time we find an inside/outside.
    emptyPointIndicies = []
    mismatchCounter = 0
    for g in gridPoints:        
        # Check if we've started a new line. If so, then we reset everything.
        # This test should precede all other test, because we don't want old knowldge
        # to carry over to the new line, since we don't know if the polygon is only partially encapsulated by the bounding box.
        if g.index > 0: # We can't check the first element, so we can skip it. 
            coordDiff = g.globalCoord - gridPoints[g.index - 1].globalCoord
            coordDiffNonzero = [x != 0 for x in coordDiff]
            if sum(coordDiffNonzero) > 1:
                # assert len(emptyPointIndicies) == 0 # When starting a new line, we shouldn't have any unknowns from the previous line
                isOutsideTracker = None
                emptyPointIndicies = []

        # There's no point testing inside/outside for points that are on the surface. So if
        # we get to one, we just skip over it.
        if g.representsPolyhedron == True:
            g.isOutside = None
            continue

        if len(g.closeFaces) == 0:
            # If it's not close to any faces, and we don't know if this row is inside/outside, then
            # we have to wait till later to figure it out
            if isOutsideTracker == None:
                emptyPointIndicies.append(g.index)
            # However, if we do know , we can just use the previous one to fill
            else:
                g.isOutside = isOutsideTracker
        
        # If there are close faces attached to it, then we need to test it for inside/outside.
        else:
            # Find centroid of all the vertices of all the close faces. This will be our endpoint
            # when casting a ray for collision testing. 
            uniquePoints = []
            # This takes just the first face and projects to the center of it.
            # [uniquePoints.append(x) for x in g.closeFaces[0] if x not in uniquePoints]
            [uniquePoints.append(x) for x in g.closeFaces[g.closestFaceIndex] if x not in uniquePoints]
            uniquePointsCoords = verticesLH[uniquePoints]
            endPoint = findPointsCenter(uniquePointsCoords)
            g.testedEndpoint = endPoint

            # Draw a ray to that point, and see if we hit a backface or not
            numHits, thisBackFace = f_ray_intersect_polyhedron(g.globalCoord,endPoint,g.closeFaces,verticesLH,False)
            
            # We can check the face as well if we want to be super precise. If they dont' agree, we then check against the entire polyhedron.
            # if superFine == True:
            #     if len(g.closeFaces) > 1:
            #         uniquePoints2 = []
            #         [uniquePoints2.append(x) for x in g.closeFaces[1] if x not in uniquePoints2]
            #         uniquePointsCoords2 = verticesLH[uniquePoints2]
            #         endPoint2 = findPointsCenter(uniquePointsCoords2)
            #         numHits2, thisBackFace2 = f_ray_intersect_polyhedron(g.globalCoord,endPoint2,g.closeFaces,verticesLH,False)
            #     if len(g.closeFaces) == 1 or thisBackFace != thisBackFace2:
            #         mismatchCounter += 1
            #         numHits, thisBackFace = f_ray_intersect_polyhedron(g.globalCoord,numpy.array([0.0,0.0,0.0]),faces,verticesLH,False)
            
            # Fill in inside outside attribute for this point, as pRayStartPos, pRayEndPos, faces, vertices, pTruncateToSegmentll as for any points before it
            g.isOutside = not thisBackFace
            isOutsideTracker = not thisBackFace
            for i in emptyPointIndicies:
                gridPoints[i].isOutside = isOutsideTracker
            # Because we have filled in all the unknowns, we can reset that counter.
            emptyPointIndicies = []

    # Final pass through for quality/sanity checks.
    for g in gridPoints:
        if g.representsPolyhedron == True:
            assert g.isOutside == None
            continue
        else:
            if g.isOutside == None:
                g.isOutside = True
    print('Flood filling grid inside/outside took ' + str(time() - timeFinishProjection) + ' seconds.')
    insidePoints = [g.index for g in gridPoints if g.isOutside == False]
    outsidePoints = [g.index for g in gridPoints if g.isOutside == True]
    surfacePoints = [g.index for g in gridPoints if g.representsPolyhedron == True]

    visualizeAtomArray('insidePoints',[points[i] for i in insidePoints],0.5)
    visualizeAtomArray('outsidePoints',[points[i] for i in outsidePoints],0.5)
    visualizeAtomArray('surfacePoints',[points[i] for i in surfacePoints],0.5)
    # if superFine == True:
    #     print('Superfine was ' + str(superFine) + ' and there were ' + str(mismatchCounter) + ' mismatches.')
    print('Grid construction took ' + str(time() - startTime) + ' seconds for ' + str(len(faces)) + ' faces and ' + str(len(gridPoints)) + ' points.')
    
    # print('The bad point is at global index ' + str(insidePoints[15090]))
    # visualizeAtomArray('badPoint', [gridPoints[insidePoints[15090]].globalCoord])
    # visualizeAtomArray('badPointEndPoint',[gridPoints[insidePoints[15090]].testedEndpoint])
    # visualizeAtomArray('badPointCloseFaces',verticesLH[gridPoints[insidePoints[15090]].closeFaces])
    # visualizeAtomArray('badPointMasterPoint', [gridPoints[53144].globalCoord])

    return gridPoints

projectPolyhedronToGrid2(5)
