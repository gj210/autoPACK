from AutoFill.ingr_ui import SphereTreeUI

anyNestedOrganelles = 0 # variable user can set to properly handle nested organelles, but expensive so we'll turn off by default and offer in GUI.

for organelle in organelles:
    # Calculate encapsulatingSpheres of each organelle- these will be used for optimization. Can use code from SphereTreeUI
    organelle.encapsulatingSphere = calculateEncapsulatingSphere #(use Cluster code from SphereTreeUI to calculate
                                                                 # center of gravity and 1 encapsulating sphere for)
    organelle.data = makeMeshData(organelle.mesh)
    
    # In many (most) cases, the meshes need to be subdivided or there are situations that can return the "wrong" closest vertex 
    # from BHtree. This happens when an organelle folds back with a pseudopod or fimbrae that almost touches its own surface
    # OPTIMIZATION... even though we calculate a subdivided mesh, we won't always need it (encapsulating sphere tests will
    #   check for need of subdivide vs original mesh, so we are wise to keep two meshes for each subdivided organelle if not too much RAM
    vOrganellesCount = len(organelles)
    if vOrganellesCount > 1:  # first optimization, we only need refined grids if more than one organelle surface or if concave
        organelle.subdividedMesh = subdivideMesh(organelle) # no triangle u or v can be longer than gridSpacing
                                                            # Michel made a subdivide function for this somewhere- needs to be uPy safe
                                                            # His code or new code needs to walk through every triangle, if u or v
                                                            #  if u>gridSpacing: split in half (add new vert and replace affected faces
                                                            #  repeat for v
                                                            # This shoudl return, verts, faces, fnormals, pnormals and BHtree for 
        organelle.subdividedMesh.data = makeMeshData(organelle.subdividedMesh)
    if len(organelles) == 1:
        concave = testMeshForConcavity() # this may be more expensive than subdividing the mesh, but we can only know by testing... should 
                               # depend on the details of the mesh itself, so worth building either way
        if concave:
            organelle.subdividedMesh = subdivideMesh(organelle) # if a mesh is concave, it has a chance of having a vertex      
                                                                # closer to a gridVertext during BHtree test that can cause false closest point test   
#INSIDE POINTS LOOP
for p in masterGridPoints:
    distToClosestOrganelleProxy = 0
    organelleDistsToProxies = []
    for organelle in organelles:
        distToClosestOrganelleProxyCurrent = distance(p.position(), organelle.encapsulatingSphere.radius) #distance function returns distance to 
                                    # surface of current organelle) 
        organelleDistsToProxies[organelle] = distToClosestOrganelleProxyCurrent
        if distToClosestOrganelleProxyCurrent < distToClosestOrganelleProxy:
            distToClosestOrganelleProxy = distToClosestOrganelleProxyCurrent
            closestOrganelleRadius = organelle.encapsulatingSphere.radius
    safeDistanceLimit = distToClosestOrganelleProxy + closestOrganelleDiameter # (precalculate diameter and store with organelle to avoid
                                                                                # doing muliplication every time)
    organellesToTest = []
    for oraganelle in organelleDistances:  # we don't need to check for organelles that are safely too far away!
        # e.g. if one organelle is a cube with only 8 very far apart points and an organelle that is a fine blobby mesh
        # like an amoeba comes near the surface of that cube, when getting the points and polygons to use for testing,
        # we need to get the far away cube points, so points near the inside surface of the cube will test against the cube
        # and will not test only against the adjacent blobby organelle that has closer mesh points... however, if that organelle
        # is farther away than the encapsulating radius of the cube, we don't need to test it!  Because the BHTree will return
        # The cube vertices as the closest points and we can get the correct polygons to test from that for a very efficient test
        # The same optimization will apply to convex meshes... we don't need to use the subdivided mesh if the surface is 100% convex
        if organelleDistsToProxies(organelle) <= safeDistanceLimit:
            organellesToTest.append(organelle)
    for organelle in organellesToTest:
#        if len(organellesToTest) == 1 && organelle.concave == 0: # organelle is convex, so we can use simpler mesh, not subdivided... 
#                                            #  eliminate this test if we find testing for concavity is to expensive (see above)
#            #doOriginal dot product test against closest point normal... you can even use the closest face normal for this
#            #and coudl add an optimization up front to say that if surface is convex we don't need to calculate point normals
#            pointNormals = organelle.data.pnormals
#
#            modifiedVersionOfOriginalDotProductAlgorithm(p.pos, organelle.BHTree, pointNormals) #return inside/outside and 
#                                                                                                    #distance (negative if inside)
# Acutally, this faster option won't work for a very thin organelle... need to just use the test below... 
# could let user say use this if they deem no thin sections, but turning it off for now
        if len(organellesToTest) > 1 
            pointNormals = organelle.subdividedMesh.data.pnormals # this is only reliable if you use subdivided meshes (see large faced
                                                                  # cube near amoeboid object example described above)
            nearestPointsEachOrganelle = []
            nearestPointsOrg = []
            for organelle in organellesToTest:
                nearestPointsOrg.append(organelle.subdividedMesh.data.BHtree.vertices.pos[0] )
                if organelle.concave == 0
                    return
                else: # Need to collect patches that may be within range, but not closest to test for point inside of longer polygon, but nearer
                    # to polygon from polyhedron that is folded back onto itself
                    i=0
                    while distance(p.pos, organelle.subdividedMesh.data.BHtree.vertices[i].pos) <= smallestProteinRadius
                        nearestPointsOrg.append(organelle.subdividedMesh.data.BHtree.vertices.pos[i] )
                        i++
            if nearestPointsOrg :
                organelle.nearestPointsOrg.append(nearestPointsOrg)
            elif anyNestedOrganelles == 0: 
                pop.organellesToTest(organelle)
    polygonsToRayTest = []
    for organelle in organellesToTest:  # (These must be in nested order from outside to inside... they come in properly from histovol.organelles, 
                                        #  be sure to retain that order)
        for closePoint in organelle.nearestPointsOrg: # Optimization loop: Don't test adjacent patches!  
                                                      # Only multiple patches on polyhedra that fold back onto themselves
            polygonsToRayTestTemp = getClosePolygons(closePoint)
            for polygon1 in polygonsToRayTestTemp:
                for polygon2 in polygonsToRayTest:
                    if polygon1 != polygon2
                        polygonsToRayTest.append(polygon2, organelleID)  # somehow need to associate polygon with the organelle it comes from e.g., next line:
        organelle.polygonsToRayTest = polygonsToRayTest
        for closePoint in organelle.nearestPointsOrg:
            vRayEndPos = closePoint.pnormals * 0.0000000001
            for polygon in organelle.polygonsToRayTest:
                organelle.hitList = ray_collide_triangle(vRayStartPos(p.pos, vRayEndPos, polygon) # return list of hit polygons, hitDistace = closest point of hit, 
                                                                        # frontface or backface
            closestHit = sort(hitlist, byDistanceOfRayCollision)
                        #need to keep a closestHit for every organelle
            if closestHit:
                if abs(hitDistance) < abs(masterGridPoints[p].distance)
                    masterGridPoints[p].distance = hitDistance
                    masterGridPoints[p].organelleList[len(masterGridPoints[p].organelleList)-1].append[hitDistance] 
                if abs(hitDistace) == abs(masterGridPoints[p].distance)
                    masterGridPoints[p].organelleList.append[hitDistance] #trying to indicate multiDimensinal array where outside point is close to 
                                                                     # multiple organelles at once.
            if hitDistance < (-bilayerHalfThickness) # point is inside
                organelle.insidePoints.append(p) # need to get the correct organelle ID here... maybe use 2D array for polygonsToRayTest with organelleID?
                gridPtID[p] = -organelle.ID
            if hitDistance > 0 && hitDistace <= (bilayerThickness) # point is in bilayer
                organelle.surfacePoints.append[p]
                gridPtID[p] = organelle.ID
        # This will recursively change ptIDs for nested organelles and works from outside to Inside!


##########################################
# Pseudocode for missing functions:
#########################################
                                                         
def makeMeshData(organelle)
    verts = organelle.vertices # make sure all proper uPy right-handed
    faces = organelle.faces # make sure all proper uPy right-handed
    fnormals = organelle.normals # if you think these are broken, make a uPy safe (left vs right handed input converted
    # to right if geometry is coming from host using my COFFEE code below for
    # f_calculate_array_of_surface_normals
    organelle.pnormals = f_calculate_array_of_point_normals(verts, faces, fnormals)
    organelle.BHtree = BHtree(organelle, points, faces)


def testMeshForConcavity(organelle) #test every face against every face for angle between normals >90deg
    concave = 0
    points = []
    neighborFaces = []
    for face in organelle.faces:
        points = face.points # list of vertices that make up the triangle (should keep universal for quads too)
        for p in points: # need to efficiently get list of faces touching the current triangle
            for face in organelle.faces
                if face.vertices = p:
                    neighborFaces.append(face)
        for faceN in neighborFaces:
            if dot(face.normal, faceN.normal) >= 1 # whatever dot product test for 90 degrees is)
                concave = 1
                return concave  # end test- this mesh is concave (don't need to test more)- must subdivide mesh
    return concave
        
def getClosePolygons(ptID,organelle)
    closePolygons = []
    for polygons in organelle.faces:
        for pt in polygons:
            if pt = ptID:
                closePolygons.append(polygons)
    return closePolygons
                                                         

                                                         
                                                         
                                                         
##########################################
# COFFEE functions below need to be converted to Python and made uPy safe if planning 
#  to use (I can't recall if my cross products are left or right handded, etc.
##########################################

f_get_polygon_normal(pPointArray)
	{
	var vE1 = pPointArray[1] - pPointArray[0];  //  Get the first edge as a vector.
	var vE2 = pPointArray[2] - pPointArray[0];  //  Get the second edge.
	var vPolygonNormal = vcross(vE1, vE2);
	return vPolygonNormal;///vlen(vPolygonNormal);
}


f_calculate_array_of_point_normals(pPointPositions, pPolygonPositions)
	{
	/*For a basic vertex normal calculation create an array for the normals of the same size as the point array. Then loop thru all polygons and add 
	the surface normal for each point in the array. After that you normalize each normal in the array. Note that this will not give phongbreaks 
	because for this you need a vertex normal for each vertex for each polygon.*/
	
	var vArrayOfSurfaceNormals = f_calculate_array_of_surface_normals(pPointPositions, pPolygonPositions);
	var vArrayOfPointNormals = new (array, sizeof(pPointPositions));
	var i, j, k;
	for (i = 0; i < sizeof(pPointPositions); i++)
		{
		var vPointIndex = i;
		var vPointNormal = vector(0,0,0);
		for (j = 0; j < sizeof(pPolygonPositions)/4; j++)
			{
			var vLoopLimit = 4;  //  Default k will loop through polygon assuming its a quad.
			if (pPolygonPositions[j*4+3] == pPolygonPositions[j*4+2])  //  Test to see if quad is actually just a triangle.
				vLoopLimit = 3;  //  Set k loop to only cycle one time.
            
			for (k = 0; k < vLoopLimit; k++)
				{
				if (pPolygonPositions[k+j*4] == vPointIndex)
					{
					vPointNormal = vPointNormal + vArrayOfSurfaceNormals[j][0];//+= vArrayOfSurfaceNormals[j];
					break;
                }
            }
        }
		vArrayOfPointNormals[i] = vnorm(vPointNormal);
    }
	return vArrayOfPointNormals;
}

f_calculate_array_of_surface_normals(pPointPositions, pPolygonPositions)
	{
	var vArrayOfSurfaceNormals = new (array, sizeof(pPolygonPositions)/4, 2);
	var i;
	for (i = 0; i < sizeof(pPolygonPositions)/4; i++)
		{
		var j = i * 4;
		var vPointArray = new(array, 3);
		vPointArray[0] = pPointPositions[ pPolygonPositions[j] ];
		vPointArray[1] = pPointPositions[ pPolygonPositions[j+1] ];
		vPointArray[2] = pPointPositions[ pPolygonPositions[j+2] ];
		var vNormal = f_get_polygon_normal(vPointArray);
		var vLen = vlen(vNormal);
		vArrayOfSurfaceNormals[i][0] = vNormal/vLen;
		vArrayOfSurfaceNormals[i][1] = vLen/2;
		}
		return vArrayOfSurfaceNormals;
}