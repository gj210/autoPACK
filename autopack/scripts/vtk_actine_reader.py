# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 17:27:53 2016

@author: ludov

STED http://bmcneurosci.biomedcentral.com/articles/10.1186/1471-2202-12-16

https://sites.cns.utexas.edu/sites/default/files/synapseweb/files/2007_els_bourne_harris_dendritic_spines.pdf
http://www.redheracles.net/media/upload/research/pdf/248541201412672579.pdf
http://www.cell.com/neuron/pdf/S0896-6273(14)00251-7.pdf
GFP fused protein
The amount of postsynaptic protein in the spine head is proportional to the spine volume. 
neurone dendritic spine
surface :
Na+,K+-ATPase  (NKA or sodium pump), 
The resulting STED images showed that neuronal NKA is located in discernible pools in the spine-heads and spine-necks as well as within the connecting dendritic structures.

NMDA type glutamate receptor
AMPA type glutamate receptor
inside : 
G-actin
F-actin
calmodulin depedant protein kinase II
cofilin

neurotramsitter receptor (GluA1 subunit of AMPAR)
signal transduction molecules (alpha and beta subunit CaMKII)
PSD scaffolding protein (PSD-95, Homer1b, Shank1b, SAP97)
Within the PSD, there are over 300 individual shank molecules, roughly 5% of the total protein molecules within the PSD
actin and regulatory protein (cofilin-1, actin interacting protein 1, p21 subunit of Arp2/3, profilin IIA)
F-actin cross-link (drebin A, alpha-actinin2,CaMKIIbeta )
dendritic structural protein (septin7)


parse branche.vtk
 the points are defined first after
DATASET POLYDATA 
POINTS n dataType 
p 0x p 0y p 0z 
p 1x p 1y p 1z 
then the line indices
LINES n size 
numPoints 0 , i 0 ,j 0 ,k 0 , ... 
numPoints 1 , i 1 ,j 1 ,k 1 , ... 
... 
numPoints n-1 , i n-1 ,j n-1 ,k n-1 , ... 
then some attributes 
SCALARS dataName dataType numComp 
LOOKUP_ TABLE tableName 
s 0 
s 1 
... 
s n-1 
Pierrick string set super fast 
"""
from scipy.spatial import cKDTree
from heapq import heappush, heappop
import numpy as np
import pdb

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def dist(x,y):   
    return np.sqrt(np.sum((x-y)**2))
    
def traverse(xyz, tree, indice,distance_pt, segments, segments_key):
    query_res = tree.query(xyz[indice],k=10)
    res = query_res[1][query_res[0] < distance_pt]
    #print query_res
    for i in res :
        if i == indice or i in segments[segments_key] : continue
        segments[segments_key].append(i)
        segments = traverse(xyz, tree, i ,distance_pt, segments, segments_key)
    return segments


#
#def checkExist(indice, segments):
#    found = False
#    found = indice in segments
#    if not found :
#        for key in segments:
#            found = indice in segments[key]
#            if found :
#                break
#    return found
#  
def traverse_edge(xyz, tree, indice, segments, segments_key, branch = 0, start = False):
    degree = tree.degree(indice) # branching ?
    res = tree.edge[indice].keys() #previous - next and branch
    res.sort()
    if degree == 2 and not start:
        if res[1] in segments[segments_key][branch]:
            res = [res[0]]
        else :
            res = [res[1]]
    if (degree == 2 and start) or ( degree > 2 ):
        if not start :
            print "branching ",indice, degree,len( segments[segments_key]),segments_key, start
            #change the segments_key for the branching i+1  
            #BRANCH segments_key IS RES [0]
            if int(res[0]) not in segments[segments_key][0]:
                segments = traverse_edge(xyz, tree, int(res[0]) , segments, segments_key, branch = 0)
            if (indice not in segments):
                segments[indice]=[]
                for i in range(1,degree):
                    segments[indice].append( [indice] )
            for i in range(1,degree):
                if i >= len( segments[indice]):
                    segments[indice].append([int(res[i])])
                else :
                    if int(res[i]) not in segments[indice][i-1]:
                        segments[indice][i-1].append(int(res[i]))
                    else :
                        break    
                pdb.set_trace()
                segments = traverse_edge(xyz, tree, int(res[i]) , segments, indice, branch = i-1)
            return segments
        # go through each branch
        for i in range(degree):
            if i >= len( segments[segments_key]):
                segments[segments_key].append([int(res[i])])
            else :
                if int(res[i]) not in segments[segments_key][i]:
                    segments[segments_key][i].append(int(res[i]))
                else :
                    continue
            segments = traverse_edge(xyz, tree, int(res[i]) , segments, segments_key, branch = i)
        return segments
    i=res[0]
    if int(i) == indice or int(i) in segments[segments_key][branch] : 
        return segments
    segments[segments_key][branch].append(int(i))
    segments = traverse_edge(xyz, tree, int(i) , segments, segments_key, branch = branch)

    return segments

def simplifyGraph(G):
    H = G.copy()
    Nnode= G.number_of_nodes()
    toremove=[]
    #allpts = range(Nnode)  
    #Nedge = G.number_of_edges()
    #start = 0
    for i in H.node.keys():#range(Nnode):
        degree = H.degree(i)
        edges = H.edge[i].keys()
        if degree == 2 :
            #d=dist(xyz[i],xyz[start])
            #not a branch
            #distance from start
            #if d < distance :           
            H.remove_node(i)
            H.add_edge(edges[0],edges[1])
       #     toremove.append(i)
            #H.remove_node(i)
            #H.add_edge(edges[0],edges[1])
        elif degree <1:     
            print i
#            H.remove_node(i)
        else :
            continue
    #for n in toremove:
    #    edges = H.edge[i].keys()
    #    H.add_edge(edges[0],edges[1])
    #    H.remove_node(i)
    return H

def interpolate(liste_pts,coords,distance):
    s=[]
    s.append(coords[liste_pts[0]].tolist())
    prev_points = coords[liste_pts[0]]
    #traverse the segment and interpolate
    for i in range(1,len(liste_pts)-1):
        d=dist(coords[liste_pts[i]],prev_points)#coords[liste_pts[i-1]])
        if d < distance :
            continue
        #else interpolate at the exact position
        direction = coords[liste_pts[i]]-coords[liste_pts[i-1]]
        ndirection = unit_vector(direction)
        new_points = coords[liste_pts[i]]#prev_points+ndirection*distance
        s.append(new_points.tolist())
        prev_points = new_points
    s.append(coords[liste_pts[-1]].tolist())
    return s
    
def rec_traverse_merge():
    pass

def makeFaceFromEdge(graph):
    nodes = np.array(graph.node.keys(),int)
    faces=[]
    search=set()
    count = 0
    total = graph.number_of_edges()
    for i in range(len(nodes)):
        edgs = graph.edge[nodes[i]].keys()
        i_edgs = np.searchsorted(nodes,np.array(edgs,int))
        if (count % 2000)==0 :
             print count,total,len(edgs)
        for j in range(len(edgs)):
             f1=[i,i_edgs[j]]
             f2=[i_edgs[j],i]
             if str(f2[0])+"_"+str(f2[1]) in search :
                continue
             faces.append(f1)  
             search.add(str(f1[0])+"_"+str(f1[1]))
        count+=1
    return faces


def getAllStartingPoint(graph):
    res=[]
    for n in graph.degree_iter():
        if n[1]==1:
            res.append(n[0])
    return res

def findNext(xyz,edges,current,prev):
    v1 = xyz[current]-xyz[prev]
    angles = []
    for i in range(len(edges)):
        if edges[i] == current :
            angles.append(-1)
            continue
        if edges[i] == prev :
            angles.append(-1)
            continue
        v2=xyz[current]-xyz[edges[i]]
        angles.append(angle_between(v1,v2))
    max_angle = max(angles)
    indice = angles.index(max_angle)
    if (max_angle < np.pi/2):
        indice = -1
    return indice,angles

def checkExist(listeSegments, query):
    #check if already visit the query
    for s in listeSegments:
        if query in s:
            return True
    return False
    
def mergeSegment(prev_node, current_node, segments,depth,allsegments):
    if depth >100:
        return depth,segments
    nedgs=H.edge[current_node].keys()
    if H.degree(current_node) == 1:
        if nedgs[0] not in segments :
            segments.append(nedgs[0])
        #print "return degree 1 ", current_node, depth
        return depth,segments
    elif H.degree(current_node) == 2:
        toadd = nedgs[0]
        if (nedgs[0]==prev_node):
            toadd = nedgs[1]
        if toadd not in segments :
            segments.append(toadd)
        #print "recursif degree 2 ", current_node, depth
        depth,segments = mergeSegment(current_node,toadd, segments,depth,allsegments)
    else :
        ind,angles = findNext(xyz,nedgs,current_node,prev_node)
        #print ind, max(angles)
        if ind != -1:
            if checkExist(allsegments,nedgs[ind]):
                #print "return used already ", current_node, depth
                return depth,segments
            else :   
                segments.append(nedgs[ind])
                depth+=1
                depth,segments = mergeSegment(current_node, nedgs[ind], segments,depth,allsegments)
        else :
            #print "return bad angle ", current_node, depth
            return depth,segments
    return [depth,segments]
    
        
def getSegments(graphS,graph):
    elemid = []
    search =set()
    for ed in graphS.edge :
        edgs = graphS.edge[ed].keys()
        for i in range(len(edgs)):
            if graphS.degree(ed) == 1 and  graphS.degree(edgs[i]) == 1:
                search.add(str(int(ed))+"_"+str(int(edgs[i])))
                search.add(str(int(edgs[i]))+"_"+str(int(ed)))
                continue
            if str(int(edgs[i]))+"_"+str(int(ed)) in search :
                continue
            if str(int(ed))+"_"+str(int(edgs[i])) in search :
                continue
            li = nx.shortest_path(graph,int(ed),int(edgs[i]))#give back all the point for this segment
            search.add(str(int(ed))+"_"+str(int(edgs[i])))
            search.add(str(int(edgs[i]))+"_"+str(int(ed)))
            if len(li) <= 2:
                continue
            elemid.append(li)
    return elemid
    
def showGraph(graph,xyz):
    import upy
    helper = upy.getHelperClass()()
    nodes = np.array(graph.node.keys(),int)
    faces = makeFaceFromEdge(graph)
    test = helper.createsNmesh("tesT", xyz[nodes], None, faces)

def showSpline(graphS,graph,xyz):
    import upy
    helper = upy.getHelperClass()()
    segments = getSegments(graphS,graph)
    a=0
    for s in segments :
        if len(s) > 1 :
            helper.spline("test_"+str(a), xyz[s])
            a+=1

def toBinary(filename,points,branch):
    from struct import pack
    #format is N points followed by float*3
    fptr = open(filename, "wb")
    data = pack('i', len(points))
    fptr.write(data)
    np.array(points, 'f').flatten().tofile(fptr)  # 4float position
    data = pack('i', len(branch))
    fptr.write(data)
    for b in branch:
        data = pack('i', len(b))
        fptr.write(data)
        np.array(b, int).flatten().tofile(fptr)  # 4flaot quaternion
    fptr.close()
            
#def simplifyLine(G,distance):
#    #resample and interpolate according distance_interval
#    # one segment is between two node with branch > 2
#    
import sys
sys.path.append("C:\\Users\\ludov\\Anaconda3\\envs\\py27\\Lib\\site-packages\\")
import pyvtk
filename ='C:/Dev/flexpack_dev_1.0/data/branches.vtk'
vtk = pyvtk.VtkData('C:/Dev/flexpack_dev_1.0/data/branches.vtk')
nSegments = len(vtk.structure.lines)
xyz = np.array(vtk.structure.points,float)
import networkx as nx
G=nx.Graph()
for seg in vtk.structure.lines:
    G.add_path(seg)

#reduce the graph ?
H = simplifyGraph(G)
nodes = np.array(H.node.keys(),int)

coords = xyz[nodes]
H_kdtree = cKDTree(coords)

show = False
export = False
reconnect = True
merge_longest = False
if reconnect:
    for i in range(len(coords)):
        indices = H_kdtree.query_ball_point(coords[i],28)[1:]
        #add edge
        for l in indices :
            if nodes[l] not in H.edge[nodes[i]].keys():
                if (int(l)==i): continue
                #H.add_edge(nodes[i],nodes[int(l)])
                G.add_edge(nodes[i],nodes[int(l)])#here as well?
    H = simplifyGraph(G)
if merge_longest :
    s=getAllStartingPoint(H)
    #check closest ?
    #traverse from here
    
    
    branch = getSegments(H,G)
    #branch start and end are node in H
    #can mesure angle and continuity ?
    
#build the simplified network

if show :
    showGraph(H,xyz)
    showSpline(H,G,xyz)
#write the points
branch = getSegments(H,G)
#write the somplified graph
if export :
    toBinary("C:/Dev/flexpack_dev_1.0/data/cache_raw_new.bin",xyz,branch)
    nx.write_gml(H,'C:/Dev/flexpack_dev_1.0/data/graph_node.gml')
    nx.write_gml(G,'C:/Dev/flexpack_dev_1.0/data/graph_node_all.gml')


actine_radius = 27.6/2.0
distance_pt=(27.6/100.0)#0.31#

if False:    
    filename ='C:/Dev/flexpack_dev_1.0/data/cache.bin_points_raw.txt'
    data = np.loadtxt(filename,delimiter=" " )
    x,y,z,r = data.transpose()
    xyz = np.column_stack([x,y,z])
    
    filename ='C:/Dev/flexpack_dev_1.0/data/cache.bin_graph_raw.txt'
    data_graph = np.loadtxt(filename,delimiter=" " )
    
    import networkx as nx
    G=nx.Graph()
    G.add_nodes_from(range(len(xyz)))
    G.add_edges_from(data_graph.astype(int))
    dG=nx.DiGraph()
    dG.add_edges_from(data_graph)
    
    #need to resample to reduce the size, per
    H = simplifyGraph(G,0,xyz,0)
    #H is only the branching node
    nodes = np.array(H.node.keys(),int)
    
    
    #    
    #for ed in H.edge :
    #     ind = np.searchsorted(nodes, int(ed))
    #     edgs = H.edge[ed].keys()
    #     if (count % 1000)==0 :
    #         print count,total,len(edgs)
    #     for i in range(len(edgs)):
    #         f1=[int(ind),np.searchsorted(nodes,int(edgs[i]))]
    #         f2=[f1[1],f1[0]]
    #         if f1 in faces or f2 in faces :
    #             continue
    #         faces.append(f1)
    #     count+=1
    
    xyz = xyz * 100.0
    coords = xyz[nodes]
    tree_node = cKDTree(coords)
    
    
    #reconnet when close to each other < 28
    for i in range(len(coords)):
        indices = tree_node.query_ball_point(coords[i],28)[1:]
        #add edge
        for l in indices :
            if nodes[l] not in H.edge[nodes[i]].keys():
                if (int(l)==i): continue
                H.add_edge(nodes[i],nodes[int(l)])
                G.add_edge(nodes[i],nodes[int(l)])
                #print i,nodes[i],l,nodes[l]
    
    nx.write_gml(H,'C:/Dev/flexpack_dev_1.0/data/graph_node.gml')
    nx.write_gml(G,'C:/Dev/flexpack_dev_1.0/data/graph_node_all.gml')
    
        

    H2 = simplifyGraph(H,0,xyz,0)
    nodes = np.array(H2.node.keys(),int)
    coords = xyz[nodes]
    
    #gather starting node with degree 1
    starting_node=[]
    for key in H.nodes_iter():
        if (H.degree(key) == 1) :
            starting_node.append(key)
    
    def traverse(node,segments):
        edgs = H.edge[node].keys()
        for e in edgs:
            traverse(e)
        return edgs
    
    def findNext(xyz,edges,current,start):
        v1 = xyz[current]-xyz[start]
        angles = []
        for i in range(len(edges)):
            if edges[i] == current :
                angles.append(-1)
                continue
            v2=xyz[current]-xyz[edges[i]]
            angles.append(angle_between(v1,v2))
        max_angle = max(angles)
        indice = angles.index(max_angle)
        if (max_angle < np.pi/2):
            indice = -1
        return indice,angles
    
    def checkExist(listeSegments, query):
        #check if already visit the query
        for s in listeSegments:
            if query in s:
                return True
        return False
        
    def mergeSegment(prev_node, current_node, segments,depth,allsegments):
        if depth >100:
            return depth,segments
        nedgs=H.edge[current_node].keys()
        if H.degree(current_node) == 1:
            if nedgs[0] not in segments :
                segments.append(nedgs[0])
            #print "return degree 1 ", current_node, depth
            return depth,segments
        elif H.degree(current_node) == 2:
            toadd = nedgs[0]
            if (nedgs[0]==prev_node):
                toadd = nedgs[1]
            if toadd not in segments :
                segments.append(toadd)
            #print "recursif degree 2 ", current_node, depth
            depth,segments = mergeSegment(current_node,toadd, segments,depth,allsegments)
        else :
            ind,angles = findNext(xyz,nedgs,current_node,prev_node)
            #print ind, max(angles)
            if ind != -1:
                if checkExist(allsegments,nedgs[ind]):
                    #print "return used already ", current_node, depth
                    return depth,segments
                else :   
                    segments.append(nedgs[ind])
                    depth+=1
                    depth,segments = mergeSegment(current_node, nedgs[ind], segments,depth,allsegments)
            else :
                #print "return bad angle ", current_node, depth
                return depth,segments
        return [depth,segments]
    
       
    nodes_done=[]
    segments=[]
    allsegments=[]
    done = False
    c=0
    total_node=0
    while not done :
        #traverse
        if checkExist(allsegments, starting_node[c]):
            c+=1
            continue
        edgs = H.edge[starting_node[c]].keys()
        # check angle
        segments.append(starting_node[c])
        segments.append(edgs[0])
        depth,segments=mergeSegment(starting_node[c], edgs[0], segments,0,allsegments)
        total_node=total_node+len(segments)
        allsegments.append(segments)
        segments=[]
        c+=1
        if c >= len(starting_node):
            done = True
       
    #nx.write_gml(G,"test.gml")
    
    #faces=[]
    #search=set()
    #count = 0
    #total = H2.number_of_edges()
    #for i in range(len(nodes)):
    #    edgs = H2.edge[nodes[i]].keys()
    #    i_edgs = np.searchsorted(nodes,np.array(edgs,int))
    #    if (count % 2000)==0 :
    #         print count,total,len(edgs)
    #    for j in range(len(edgs)):
    #         f1=[i,i_edgs[j]]
    #         f2=[i_edgs[j],i]
    #         if str(f2[0])+"_"+str(f2[1]) in search :
    #            continue
    #         faces.append(f1)  
    #         search.add(str(f1[0])+"_"+str(f1[1]))
    #    count+=1
    #resample
    #create subdivide segment interpolated ?
    
    segments = []
    search=set()
    new_xyz=[]
    new_indices=[]
    old_indices=[]
    indices = set()
    correct_segments=[]
    for ed in H2.edge :
        edgs = H2.edge[ed].keys()
        for i in range(len(edgs)):
            if str(int(edgs[i]))+"_"+str(int(ed)) in search :
                continue
            li = nx.shortest_path(G,int(ed),int(edgs[i]))#give back all the point for this segment
            asegment=[]
            #first and last shouldnt change
            #add only if not there already
            if int(ed) in indices:
                asegment.append(old_indices.index(int(ed)))
                old_indices.append(-1)
            else :
                new_xyz.append(xyz[ed])#start
                asegment.append(len(new_xyz))
                old_indices.append(ed)
                indices.add(ed)
            start = len(new_xyz)
            new_points =  interpolate (li, xyz, 27.6)
            new_xyz.extend(new_points)
            asegment.extend(np.array(range(len(new_points)))+start)     
            if int(edgs[i]) in indices:
                asegment.append(old_indices.index(int(edgs[i])))
                old_indices.append(-1)
            else :
                new_xyz.append(xyz[int(edgs[i])])#start
                new_indices.append(len(new_xyz))
                old_indices.append(edgs[i])  
                indices.add(edgs[i])
            segments.append(asegment)
            search.add(str(int(ed))+"_"+str(int(edgs[i])))
            search.add(str(int(edgs[i]))+"_"+str(int(ed)))
    #how deal with connecting segments, whoch should share the position
    #new_xyz
    #indices
            
    #create a new structure from here
    #with xyz
    #segment IDs
    
    #segments=[]
    #for s in allsegments:
    #    points = list(s)
    #    segs=[]
    #    for i in range(len(points)-1):
    #        li = nx.shortest_path(G,int(points[i]),int(points[i+1]))
    #        new_points =  interpolate (li, xyz, 27.6)
    #        segs.extend(new_points)
    #    segments.append(segs)

    segments=[]
    elemid = []
    search =set()
    for ed in H.edge :
        edgs = H.edge[ed].keys()
        for i in range(len(edgs)):
            if str(int(edgs[i]))+"_"+str(int(ed)) in search :
                continue
            li = nx.shortest_path(G,int(ed),int(edgs[i]))#give back all the point for this segment
            new_points =  interpolate (li, xyz, 27.6)
            #conectivity?
            
    import upy
    helper = upy.getHelperClass()()
    
    #tree_node = cKDTree(xyz[list_branching_point])
    #t = helper.PointCloudObject("pts_cloud",vertices=xyz[nodes])
    #t = helper.createsNmesh("tesT", coords, None, faces)
    a=0
    for s in segments :
        if len(s) > 1 :
            helper.spline("test_"+str(a), s)
            a+=1
    #resample ?
    
    #test = list(nx.bfs_edges(G,0))
    #traverse the graph
    #allpts = range(len(data))    
    #segments={}
    #done = False
    #cutoff=1
    #count=0
    #scount=0
     #node with nore than 2 edges, or starting point with two edges
    #b
    #while not done:
    #    degree = G.degree(allpts[0])
    #    segments[allpts[0]] = []
    #    for i in range(degree):
    #         segments[allpts[0]].append([allpts[0]])
    #    segments=traverse_edge(xyz,G,allpts[0],segments,allpts[0],start=True)
    #    for i in segments[allpts[0]] :
    #        for j in i :
    #            if j in allpts:
    #                allpts.remove(j)
    #    count+=1
    #    if count >= cutoff:
    #        break
    
    #list_branching_point =[]
    #for s in segments :
    #    if s not in list_branching_point :
    #        list_branching_point.append(s)
    #    for branch in segments[s]:
    #        end = branch[-1]
    #        if end not in list_branching_point :
    #            list_branching_point.append(end)
    #need first and last point of  each segments. make a tree and check if should merge
    #import upy
    #helper = upy.getHelperClass()()
    #
    #
    ##tree_node = cKDTree(xyz[list_branching_point])
    ##t = helper.PointCloudObject("pts_cloud",vertices=xyz[nodes])
    #t = helper.createsNmesh("tesT", coords, None, faces)
    #a=0
    #for s in segments :
    #    for branch in segments[s]:
    #        if len(branch) > 1 :
    #            helper.spline("test_"+str(s)+"_"+str(a), xyz[branch])
    #            a+=1
    #a=0
    #for s in segments :
    #    for branch in segments[s]:
    #        if len(branch) > 1 :
    #            helper.spline("test_"+str(s)+"_"+str(a), xyz[branch])
    #            a+=1
    
    
    #test = nx.bfs_successors(G,0)
    #res=list(nx.edge_dfs(G, 0, orientation='ignore'))
    
    #
    #allpoints=[res[0][0]]
    #for r in res:
    #    allpoints.append(r[1])
    #done = False
    #cutoff=50
    #elems=[[0]]
    #i=0
    #count=0
    #scount=0
    #allpts.remove(0)
    #icount=0
    #found = True
    #while not done:
    #    if i not in test:
    ##        break
    #        icount=0
    #        while i not in test :
    #            #pick another lead
    #            if icount >= len(allpts):
    #                found = False
    #                break
    #            i = int( allpts[icount])
    #            icount+=1
    #            found = True
    #        elems.append([i])
    #        scount+=1
    #        allpts.remove(i)
    #        if not found :
    #            break   
    #        test = nx.bfs_successors(G,i)
    #    #do it for each branch ?
    #    next_elem = int(test[i][0])
    #    if next_elem not in elems[scount]:
    #        elems[scount].append(next_elem)
    #    if next_elem in allpts :
    #        allpts.remove(next_elem)
    #    i = next_elem
    #    count+=1
    #    if count > cutoff:
    #        break
    #    if not len(allpts):
    #        break
    #
    #
    #tree = cKDTree(xyz)
    #Npoints = len(xyz)
    #ordered_pair = []
    #segments={}
    #finished = False
    #i=0
    #j=0
    #
    #segments[0] = [0]    
    #
    #test = traverse(xyz, tree, 0,distance_pt, segments, 0)
    #
    #while not finished:
    #    segments = traverse(xyz, tree, i,distance_pt, segments, i)
    #    found = False
    #    while not found:
    #        i+=1
    #        found = not checkExist(i, segments)
    #    segments[i]=[i]
    #    if i > cutoff:
    #        finished = True
#    


#helper.spline("test", xyz[allpoints])

               
#    res = tree.query_ball_point(xyz[i],distance_pt)
#    if len(res) > 1:
#        for j in res:
#            if (i==j) continue
#            new_res = tree.query_ball_point(xyz[j],distance_pt)
#            if len(new_res) > 2:
#                
#        res = tree.query(xyz[i],3)[1][1]
#    if (i in ordered_pair):
#        print i
#        i = tree.query(xyz[i],3)[1][2]
#        print i
#    if len(ordered_pair) > cutoff:
#        break
#execfile('C:/Users/ludov/OneDrive/Documents/autoPACK/autopack/scripts/vtk_actine_reader.py')
