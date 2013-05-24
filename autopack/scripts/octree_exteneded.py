### Pymel Expanding Object-Oriented Octree v.1
### http://zoomy.net/2010/02/07/pymel-oooctree/
### based on Ben Harling's Python Octree v.1
### http://code.activestate.com/recipes/498121/

#from pymel import *
#import pymel.core.datatypes as dt
import random
import upy

global MAX_OBJECTS_PER_NODE, MINIMUM_SIZE, DIRLOOKUP
MAX_OBJECTS_PER_NODE = 10 # max number of objs allowed per node  
MINIMUM_SIZE = 1 # prevents infinitely small nodes
DIRLOOKUP = {"3":0, "2":1, "-2":2, "-1":3, "1":4, "0":5, "-4":6, "-3":7} 
# used by findOctant() to get the correct octant index

### START CLASS OCTNODE

class OctNode:
  def __init__(self, position, size, objects):
    self.position = position
    self.size = size

    # store our objects
    self.objects = []
    if len(objects) > 1: self.objects.extend[objects]

    # give it some empty subnodes
    self.subnodes = [None, None, None, None, None, None, None, None]

    # empty subnodes don't count as subnodes
    self.hasSubnodes = False

    # once it's been subdivided, don't accept any more objects
    self.hasSubdivided = False

    # node's bounding coordinates
    self.ldb = (position[0] - (size / 2), position[1] - (size / 2), position[2] - (size / 2))
    self.ruf = (position[0] + (size / 2), position[1] + (size / 2), position[2] + (size / 2))

    self.bounds = [self.ldb, self.ruf]
### END CLASS OCTNODE



### START CLASS OCTREE
      
class Octree:
    def __init__(self, initSize,helper = None):
        # init the octree's root cube at the world origin
        self.root = self.addNode([0,0,0], initSize, [])
        self.helper = helper
        if self.helper is None :
            self.helper = upy.getHelperClass()()
        
    def addNode(self, position, size, objects):
        # creates an OctNode
        return OctNode(position, size, objects)
        
    # is obj bb entirely inside node bb? returns True or False
    def enclosesObj(self, node, obj,bb=None):
        if bb is None :
            if hasattr(obj,"bb"):
                bb=obj.bb
            else :
                bb = self.helper.getBoundingBox(obj)
        nbb = node.bounds
        # min and max of the bb
        bbmin = bb[0]
        bbmax = bb[1]
        nbbmin = nbb[0]
        nbbmax = nbb[1]
        return bbmin[0] > nbbmin[0] and \
               bbmin[1] > nbbmin[1] and \
               bbmin[2] > nbbmin[2] and \
               bbmax[0] < nbbmax[0] and \
               bbmax[1] < nbbmax[1] and \
               bbmax[2] < nbbmax[2]

    # is obj bb inside node bb at all? returns True or False
    def containsObj(self, node, obj,bb=None):
        if bb is None :
            if hasattr(obj,"bb"):
                bb=obj.bb
            else :
                bb = self.helper.getBoundingBox(obj)
        nbb = node.bounds
        bbmin = bb[0]
        bbmax = bb[1]
        nbbmin = nbb[0]
        nbbmax = nbb[1]
        return bbmin[0] < nbbmax[0] and \
               bbmin[1] < nbbmax[1] and \
               bbmin[2] < nbbmax[2] and \
               bbmax[0] > nbbmin[0] and \
               bbmax[1] > nbbmin[1] and \
               bbmax[2] > nbbmin[2]

    def findOctantCenter(self, size, pos, octant):
        newCenter = (0,0,0) # initialize vector
        offset = size / 2.0
        if octant == 0:
          # left down back
          newCenter = (pos[0]-offset, pos[1]-offset, pos[2]-offset )
        elif octant == 1:
          # left down forwards
          newCenter = (pos[0]-offset, pos[1]-offset, pos[2]+offset )
        elif octant == 2:
          # right down forwards
          newCenter = (pos[0]+offset, pos[1]-offset, pos[2]+offset )
        elif octant == 3:
          # right down back
          newCenter = (pos[0]+offset, pos[1]-offset, pos[2]-offset )
        elif octant == 4:
          # left up back
          newCenter = (pos[0]-offset, pos[1]+offset, pos[2]-offset )
        elif octant == 5:
          # left up forward
          newCenter = (pos[0]-offset, pos[1]+offset, pos[2]+offset )
        elif octant == 6:
          # right up forward
          newCenter = (pos[0]+offset, pos[1]+offset, pos[2]+offset )
        elif octant == 7:
          # right up back
          newCenter = (pos[0]+offset, pos[1]+offset, pos[2]-offset )
        return newCenter    

    # find center of a node's octant
    def findCenter(self, size, node, obj):
        octant = self.findOctant(node, obj)
        # new center = parent position + (octant direction * offset)
        pos = node.position
        return self.findOctantCenter(size, pos, octant)

    # superdivide the root toward obj
    def superdivide(self, obj):
        node = self.root
        # make a new root node twice as large as the old
        size = node.size
        newCenter = self.findCenter(size, node, obj)
        newRoot = self.addNode(newCenter, size*2, [])
        # make node a subnode of newRoot - octant is the corner of newRoot 
        #toward node, relative to obj
        octant = self.findOctant(newRoot, node)
        oldRoot = self.root
        newRoot.subnodes[octant] = oldRoot
        self.root = newRoot
        newRoot.hasSubnodes = True # default is False, override

    # split node into subnodes
    def subdivide(self, node):
        size = node.size
        node.hasSubdivided = True
        for i, item in enumerate(node.subnodes):
            # replace any "None" placeholders with nodes
            if item == None:
                # Find the subnode's center point based on octant
                pos = node.position
                newCenter = self.findOctantCenter(size/2, pos, i)
                newNode = self.addNode(newCenter, size/2, [])
                node.subnodes[i] = newNode
    
        # reassign contents of the node into the new nodes
        for obj in node.objects:
            inNodes = self.findContainingNodes(obj, node)
            self.insertNode(obj, inNodes)
        node.objects = []

    # attempt to assign obj to nodes in list of containing nodes
    def insertNode(self, obj, inNodes):
        global MAX_OBJECTS_PER_NODE, MINIMUM_SIZE
        # if tree does not contain obj, expand.
#        while not self.enclosesObj(self.root, obj):
#            self.superdivide(obj)
#            inNodes = [self.root]
        returnNodes = []
        for node in inNodes:
            if node.hasSubdivided == False:#leaf?
                # update node obj list
                if not obj in node.objects:
                    node.objects.append(obj)
                    returnNodes.append(node)
            # if too full, subdivide
            if len(node.objects) >= MAX_OBJECTS_PER_NODE:
                if node.size > MINIMUM_SIZE:
                    self.subdivide(node)
        return returnNodes
  
    # finds node's octant that contains or is closest to obj
    def findOctant(self, node, obj):
        if node == None: return None
    
        vec1 = node.position
        try: vec2 = obj.position # if obj is an octree node
        except: vec2 = self.ToVec(self.helper.getTranslation(obj)) # if obj is a transform
    
        result = 0
        # Equation created by adding nodes with known octant directions into the tree, 
        # and comparing results. See DIRLOOKUP above for the corresponding return 
        # values and octant indices
        # I don't currently understand this.
        for i in range(3):
            if vec1[i] <= vec2[i]:
                result += (-4 / (i + 1) / 2)
            else:
                result += (4 / (i + 1) / 2)
        result = DIRLOOKUP[str(result)]
        return result

    def findContainingNode(self, obj, node):
        # Basic collision lookup that finds the leaf node containing the specified position
        # Returns the child objects of the leaf, or None if the leaf is empty or none
        if node == None:
            return None
        elif node.hasSubdivided == False:
            return node
        else:
            octant = self.findOctant(node, obj)
            return self.findContainingNode(node.subnodes[octant], obj)
            
    # returns list of nodes containing obj's bb
    def findContainingNodes(self, obj, node,bb=None):
        if not self.containsObj(self.root, obj,bb=bb): return []
        if node.hasSubdivided == False:
            verified = [node]
        else:
            verified = []
        toCheck = [node]
        while len(toCheck) > 0:
            # don't update the same list we're iterating through
            checking = toCheck[:]
            for node in checking:
                toCheck.remove(node)
                # if it has subnodes, recurse
                for subnode in node.subnodes:
                    if subnode != None and self.containsObj(subnode, obj,bb=bb):
                        if subnode.hasSubdivided == False:
                            verified.append(subnode)
                        toCheck.append(subnode)
                        if self.enclosesObj(subnode, obj,bb=bb) and node in verified:
                            verified.remove(node)
        return list(set(verified))
      
    # checks if obj collides with any objs in nodes
    def findCollisions(self, obj, nodes):
        suspects = set()
        collisions = []
        for node in nodes:
            for x in node.objects:
                suspects.add(x)
        obb = obj.getBoundingBox()
        for x in list(suspects):
            if obb.intersects(x.getBoundingBox()):
                collisions.append(x)
        return list(set(collisions))
  
### END CLASS OCTREE
#if __name__ == "__main__":
#
    ### Object Insertion Test ###
    #from AutoFill.octree_exteneded import Octree
    # So lets test the adding:
#    import random
#    import time
#    helper = upy.getHelperClass()())
#    # Create a new octree, size of world
#    # size of the grid world X Y Z 
#    myTree = Octree(1500.0000,helper=helper)
#
#    # Number of collisions we're going to test
#    NUM_COLLISION_LOOKUPS = 20
#
#    #should use currentSelection ?
#    listeObj = helper.getCurrentSelection()
#    # Insert some random objects and time it
#    Start = time.time()
#    for o in listeObj:
#        myTree.insertNode(myTree.root, 15000.000, myTree.root, o)
#    End = time.time() - Start
#
#    # print some results.
#    print str(len(listeObj)) + "-Node Tree Generated in " + str(End) + " Seconds"
#    print "Tree Leaves contain a maximum of " + str(MAX_OBJECTS_PER_CUBE) + " objects each."
#
#    ### Lookup Tests ###
#
#    # Look up some random positions and time it
#    Start = time.time()
#    for o in listeObj:
#        pos = self.ToVec(self.helper.getTranslation(obj))
#        #result = myTree.findPosition(myTree.root, pos)
#        nodes = myTree.findContainingNodes(o, myTree.root)
#        collisions = myTree.findCollisions(obj, nodes)
#        if len(collisions) > 0: print "ok ",o
#
#        ##################################################################################
#        # This proves that results are being returned - but may result in a large printout
#        # I'd just comment it out and trust me :)
#        #print "Results for test at: " + str(pos)
#        #if result != None:
#        #    for i in result:
#        #        print i.name, i.position
#        # print
#        ##################################################################################
#        
#    End = time.time() - Start
#
#    # print some results.
#    print str(NUM_COLLISION_LOOKUPS) + " Collision Lookups performed in " + str(End) + " Seconds"
#    print "Tree Leaves contain a maximum of " + str(MAX_OBJECTS_PER_CUBE) + " objects each."
#
#    x = raw_input("Press any key (Wheres the any key?):")
#
#### END OOOCTREE v0.1