# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 22:21:27 2013

@author: ludo
"""
from AutoFill import octree_exteneded as octree
from AutoFill.octree_exteneded import Octree
import random
import time
import upy


        
helper = upy.getHelperClass()()
def getBoundingBox(o,**kw):
        if o is None :
            return
        if type(o) is str:
            o = helper.getObject(o)
        m = helper.getMesh(o)
        r= m.GetRad()#size on X,Y,and Z
        c= m.GetMp()*o.GetMg()
        bb1=c-r
        bb2=c+r
        return [helper.ToVec(bb1),helper.ToVec(bb2)]
def displaysubnode(parentnode,i):
    if not parentnode.hasSubnodes and not parentnode.hasSubdivided: return
    for subnode in parentnode.subnodes:
        b=helper.box("node"+str(i),center=subnode.position,size=[subnode.size,]*3)
        print "node"+str(i),len(subnode.objects)
        for io in subnode.objects :
            print helper.getName(io)
        displaysubnode(subnode,i)
        i+=1

helper.getBoundingBox=getBoundingBox
octree.MINIMUM_SIZE = 26
# Create a new octree, size of world
# size of the grid world X Y Z 
myTree = Octree(1500.0000,helper=helper)
listeObj = helper.getCurrentSelection()
Start = time.time()
for o in listeObj: 
    myTree.insertNode(o,[myTree.root])#,15000.000,myTree.root,o)
    b=helper.box(helper.getName(o)+"box",cornerPoints=getBoundingBox(o))
End = time.time() - Start
#display the octree by creatin a box at each leaf ? or for each node ?

# print some results.
print str(len(listeObj)) + "-Node Tree Generated in " + str(End) + " Seconds"
print "Tree Leaves contain a maximum of " + str(octree.MAX_OBJECTS_PER_NODE) + " objects each."

#listeObj2 = helper.getCurrentSelection()
#o=listeObj2[0]
#nodes = myTree.findContainingNodes(o, myTree.root)
#helper.setCurrentSelection(nodes[0].objects)

root = helper.box("root",center=myTree.root.position,size=[myTree.root.size,]*3)
def displaysubnode(parentnode,i):
    if not parentnode.hasSubnodes and not parentnode.hasSubdivided: return
    for subnode in parentnode.subnodes:
        if subnode is None : continue
        b=helper.box("node"+str(i),center=subnode.position,size=[subnode.size,]*3)
        print "node"+str(i),len(subnode.objects)
        for io in subnode.objects :
            print helper.getName(io)
        displaysubnode(subnode,i)
        i+=1
displaysubnode(myTree.root,0)
#import c4d;h=c4d.af.values()[0].histoVol.values()[0];h.afviewer.displaysubnode(h.octree.root,0)
#from AutoFill.Ingredient import IngredientInstanceDrop
#jtrans,rotMat,ingr,ptid = h.molecules[0]
#dropedObject = IngredientInstanceDrop(ptid, jtrans, rotMat, ingr)
#h.octree.insertNode(dropedObject,[h.octree.root]) 
#h.afviewer.displaysubnode(h.octree.root,0)
#for subnode in root.subnodes:
#    b=helper.box("root",center=subnode.position,size=[subnode.size,]*3)
#execfile("/Users/ludo/DEV/autofill_svn/trunk/AutoFillClean/octree_test_1.py")