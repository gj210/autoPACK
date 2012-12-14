# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 11:05:22 2012

@author: -
"""
import sys
try :
    import panda3d
except :    
    p="/Developer/Panda3D/lib"
    sys.path.append(p)
    import panda3d
from panda3d.core import loadPrcFileData

loadPrcFileData("", "window-type none" ) 
# Make sure we don't need a graphics engine 
#(Will also prevent X errors / Display errors when starting on linux without X server)
loadPrcFileData("", "audio-library-name null" ) # Prevent ALSA errors 

import direct.directbase.DirectStart
    
from panda3d.bullet import BulletWorld
from panda3d.bullet import BulletDebugNode
from panda3d.bullet import BulletPlaneShape
from panda3d.bullet import BulletBoxShape
from panda3d.bullet import BulletSphereShape
from panda3d.bullet import BulletRigidBodyNode
from panda3d.core import Vec3
from panda3d.core import Vec4
from panda3d.core import Point3
from panda3d.core import TransformState
from panda3d.core import BitMask32
from panda3d.core import NodePath

def delRB(world, node):
    world.removeRigidBody(node)
    np = NodePath(node)
    np.removeNode()
        
def addRBSphere(world,pos, worldNP):
    shape = BulletSphereShape(0.6)
    inodenp = worldNP.attachNewNode(BulletRigidBodyNode("Sphere"))
    inodenp.node().setMass(0.0)
#    inodenp.node().addShape(shape)
    inodenp.node().addShape(shape, TransformState.makePos(pos))
#        spherenp.setPos(-2, 0, 4)
    inodenp.setCollideMask(BitMask32.allOn())
    world.attachRigidBody(inodenp.node())
    return inodenp.node()

def addRBCube(world,pos, worldNP):
    shape = BulletBoxShape(Vec3(12.5, 12.5, 12.5))
    inodenp = worldNP.attachNewNode(BulletRigidBodyNode("Box"))
    inodenp.node().setMass(1.0)
#    inodenp.node().addShape(shape)
    inodenp.node().addShape(shape, TransformState.makePos(pos))
#        spherenp.setPos(-2, 0, 4)
    inodenp.setCollideMask(BitMask32.allOn())
    world.attachRigidBody(inodenp.node())
    return inodenp.node()    

def collide(world,node1, node2=None):
    result1 = world.contactTest(node1)
    print "node 1 contacts ",node1.getName(), result1.getNumContacts ()    
    if node2 is not None :
        result2 = world.contactTestPair(node1, node2)
        print "node 1 node 2 contacts ",node2.getName(), result2.getNumContacts () 
        
if __name__ == "__main__":
    worldNP = render.attachNewNode('World')
    
    world = BulletWorld()
    world.setGravity(Vec3(0, 0, 0))

    # Box
    box1 =addRBCube(world,Point3(2, 0, 4), worldNP)
    box2 =addRBCube(world,Point3(0, 0, 0), worldNP)
    collide(world,box1, node2=box2)
    NodePath(box2).setPos(12, 0, 4)#this is additive to the rbnoe transform!
    collide(world,box2, node2=box1)
    NodePath(box2).setPos(100, 0, 4)#this is additive to the rbnoe transform!
    collide(world,box2, node2=box1)

    node1 = addRBSphere(world,Point3(-2, 0, 4), worldNP)
    collide(world,node1)
    box3 =addRBCube(world,Point3(-10.24, 0, 4), worldNP)
    collide(world,box3,node2=node1)
#    print "sphere contact ",result.getContacts()
#    print "sphere 1 contact ",node1, result2.getNumContacts ()    
#    if result2.getNumContacts ()  :
#        delRB(world, node1)
#    else :
#        for contact in result2.getContacts():
#              cp = contact.getManifoldPoint()
#              node0 = contact.getNode0()
#              node1 = contact.getNode1()
#              print node0.getName(), node1.getName(), cp.getDistance() 
##    dt= globalClock.getDt()    
##    world.doPhysics(0.0000001, 10, 0.008)#dt, subset n, subset time
#              
#    #################
#    node2 = addRB(world,Point3(-2, 0, 4), worldNP)
#    result2 = world.contactTest(node2)
##    print "sphere contact ",result.getContacts()
#    print "sphere 2 contact ",node2, result2.getNumContacts ()  
#    result = world.contactTestPair(node1, node2)
#    print "sphere 1&2 contact ", result.getNumContacts ()  
#    print result.getContacts() 
#    if result2.getNumContacts ()  :
#        delRB(world, node2)
#    else :
#        for contact in result2.getContacts():
#              cp = contact.getManifoldPoint()
#              node0 = contact.getNode0()
#              node1 = contact.getNode1()
#              print node0.getName(), node1.getName(), cp.getDistance()        
#    dt= globalClock.getDt()    
#    world.doPhysics(0.0000001, 10, 0.008)#dt, subset n, subset time
#
#    node3 = addRB(world,Point3(-2, 0, 4), worldNP)
#    result2 = world.contactTest(node3)
#    result = world.contactTestPair(node3, node2)
#    print "sphere 3&2 contact ", result.getNumContacts ()  
##    print "sphere contact ",result.getContacts()
#    print "sphere 3 contact ",node3, result2.getNumContacts ()    
#    if result2.getNumContacts ()  :
#        delRB(world, node3)
#    else :
#        for contact in result2.getContacts():
#              cp = contact.getManifoldPoint()
#              node0 = contact.getNode0()
#              node1 = contact.getNode1()
#              print node0.getName(), node1.getName(), cp.getDistance()            
        
#    box.setActive(True)
#    sphere.setActive(True)
#    #check contact
#    result = world.contactTestPair(sphere, box)
#    print "result1 ",result.getContacts()
#    result = world.contactTest(box)
#    print "box contact ",result.getContacts()
#    result = world.contactTest(sphere)
#    print "sphere contact ",result.getContacts()
#    spherenp.setPos(2.1, 0, 4)
##    boxnp.setPos(-2, 0, 4)
##    box.setTransform(TransformState.makePos(Point3(1.246, 0,4)))
##    sphere.setTransform(TransformState.makePos(Point3(1.246, 0,4)))
#    dt= globalClock.getDt()    
#    world.doPhysics(0.0000001, 10, 0.008)#dt, subset n, subset time
#    result = world.contactTestPair(sphere, box)
#    print "result2 ",result.getContacts()
#    result2 = world.contactTest(box)
#    print "box contact ",result2.getContacts() #gave nothing ???
#    result3 = world.contactTest(sphere)
#    print "sphere contact ",result3.getContacts() #give the contact on the Box ....
#    #move sphere
#    print box.checkCollisionWith(sphere)
#    
#    spherenp.setPos(10.1, 0, 4)
##    boxnp.setPos(-2, 0, 4)
##    box.setTransform(TransformState.makePos(Point3(1.246, 0,4)))
##    sphere.setTransform(TransformState.makePos(Point3(1.246, 0,4)))
#    dt= globalClock.getDt()    
#    world.doPhysics(0.0000001, 10, 0.008)#dt, subset n, subset time
#    result = world.contactTestPair(sphere, box)
#    print "result2 ",result.getContacts()
#    result2 = world.contactTest(box)
#    print "box contact ",result2.getContacts() #gave nothing ???
#    result3 = world.contactTest(sphere)
#    print "sphere contact ",result3.getContacts() #give the contact on the Box ....
#    #move sphere
#    print box.checkCollisionWith(sphere)

#    for contact in result2.getContacts():
#          cp = contact.getManifoldPoint()
#          node0 = contact.getNode0()
#          node1 = contact.getNode1()
#          print node0.getName(), node1.getName(), cp.getDistance()
#          print cp.getPositionWorldOnA ()
#          print cp.getPositionWorldOnB ()
##          print cp.getAppliedImpulse()
##          print cp.getLocalPointA()
##          print cp.getLocalPointB()          