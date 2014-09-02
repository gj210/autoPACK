# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 11:44:29 2014

@author: ludo
"""

import c4d
import autopack
def getIngredientInstancePos(traj, name,instance_id,frame):
        indice = traj.mapping[name][1][instance_id]
        pos = traj.getPosAt(indice,frame)
        return pos
        
def applyState(traj, h, frame):
        r =  h.exteriorRecipe
        indice = 0
        if r :
            for ingr in r.ingredients:
                for k in range(len(ingr.results)):
                    pos = getIngredientInstancePos(traj,ingr.name,k,frame)
                    autopack.helper.setTranslation(ingr.ipoly[k],pos)
                    indice+=1
                #k+=1
#        k=0
        #compartment ingr
        for orga in h.compartments:
            #compartment surface ingr
            rs =  orga.surfaceRecipe
            if rs :
                for ingr in rs.ingredients:
                    for k in range(len(ingr.results)):
                        pos = getIngredientInstancePos(traj,ingr.name,k,frame)
                        autopack.helper.setTranslation(ingr.ipoly[k],pos)
                        indice+=1
            #compartment matrix ingr
            ri =  orga.innerRecipe
            if ri :
                k=0
                for ingr in ri.ingredients:
                    for k in range(len(ingr.results)):
                        pos = getIngredientInstancePos(traj,ingr.name,k,frame)
                        autopack.helper.setTranslation(ingr.ipoly[k],pos)
                        indice+=1
def main():
    h=c4d.af.values()[0].histoVol.values()[0]
    cframe = doc.GetTime().GetFrame(doc.GetFps())
    print "cframe",cframe
    applyState(h.traj,h,cframe)