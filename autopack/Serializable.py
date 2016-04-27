# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 15:55:27 2016

@author: ludo
"""
#serializabl;e class for auopack recipe
import json


#toDO class Assambly for nesting of protein quaternary structure.
#class sAssambly
#    static_id=0
#    def __init__(self,name,**kwds):
#        self.local_id=0; #// The id of the group relative only to the parent
#        self.unique_id=sAssambly.static_id;  #// It  does not matter how this one is obtained as long as it is unique for each group
#        self.name=name;
#        self.nb=0
#        self.Ingredients=[]
#        self.Assambly=[]
#        for k in kwds:
#            setattr(self,k,kwds[k])
#        sAssambly.static_id+=1

class sCompartment(object):
    static_id = 0
    def __init__(self, name, **kwds):
        self.local_id = 0  # // The id of the compartment relative only to the parent
        self.unique_id = sCompartment.static_id  # // It  does not matter how this one is obtained as long as it is unique for each group
        self.name = name
        self.Compartments = []
        self.IngredientGroups = []
        sCompartment.static_id += 1
        for k in kwds:
            setattr(self, k, kwds[k])
            
    def addCompartment(self, compartment):
        self.Compartments.append(compartment)
        compartment.local_id = len(self.Compartments)-1

    def addIngredientGroup(self, ingrgroup):
        self.IngredientGroups.append(ingrgroup)
        ingrgroup.local_id = len(self.IngredientGroups)-1

    def to_JSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, 
            sort_keys=True, indent=4)


class sIngredientGroup(object):
    static_id=0
    def __init__(self, name, groupType, **kwds):
        self.local_id = 0  # // The id of the group relative only to the parent
        self.unique_id = sIngredientGroup.static_id  # // It  does not matter how this one is obtained as long as it is unique for each group
        self.name=name
        self.Ingredients = []
        self.groupType = groupType #currently 0 1 2 protein, fiber, lipids
        for k in kwds:
            setattr(self, k, kwds[k])
        sIngredientGroup.static_id += 1
            
    def addIngredient(self,ingredient):
        self.Ingredients.append(ingredient)

    def to_JSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, 
            sort_keys=True, indent=4)


class sIngredient(object):
    static_id = [0, 0, 0]
    def __init__(self, name, groupType, **kwds):
        self.ingredient_id = sIngredient.static_id[groupType]
        self.name = name
        self.path = ""  # // The path is made out of the local id of each parent nodes starting from the root
        for k in kwds:
            setattr(self, k, kwds[k])
        sIngredient.static_id[groupType] += 1
        
    def to_JSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, 
            sort_keys=True, indent=4)


class sIngredientFiber(object):
    static_id=0
    def __init__(self, name, **kwds):
        self.ingredient_id = sIngredientFiber.static_id
        self.name = name
        self.path = ""  # // The path is made out of the local id of each parent nodes starting from the root
        for k in kwds:
            setattr(self, k, kwds[k])
        sIngredientFiber.static_id += 1
        
    def to_JSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, 
            sort_keys=True, indent=4)


#import sys
#import os
##import c4d
#
#sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/")
#sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/PIL/")
#
#from autopack.Serializable import *
#
#
#c1=sCompartment("test1")
#c2=sCompartment("test2")
#
#g1=sIngredientGroup("g1")
#g2=sIngredientGroup("g2")
#
#i1=sIngredient("i1")
#i2=sIngredient("i2")
#
#g1.addIngredient(i1)
#g2.addIngredient(i2)
#
#c1.addIngredientGroup(g1)
#c2.addIngredientGroup(g2)
#
#c1.addCompartment(c2)
