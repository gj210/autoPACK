# -*- coding: utf-8 -*-
"""
Created on Fri Apr 4 2014

@author: Ludovic Autin
"""
import c4d
import autopack
import numpy
import json
#Welcome to the world of Python
import upy
helper = upy.getHelperClass()()
def main():
    helper.doc = doc    
    PS = doc.GetParticleSystem()
    ob = op.GetObject()
    mx = ob.GetMg()
    cframe = doc.GetTime().GetFrame(doc.GetFps())
    doit = True
    PS = doc.GetParticleSystem()
    output=[]
#    if hasattr(h,"has_particle") :
#        doit = not h.has_particle
    if cframe == 0 :#and doit:
        liste_particles_coords=[]#coordinates
#        liste_particles_radius=[]#coordinates
        group_name=ob.GetName()
        if ob.GetType() == c4d.Opolygon :
            liste_particles_coords=helper.getMeshVertices(ob,transform=True)
        else :
            childrens=helper.getChilds(ob)
            if len(childrens):
#                liste_particles_coords=[helper.ToVec(helper.getTranslation(ch)) for ch in childrens]
                for ch in childrens:
                    coord = helper.ToVec(helper.getTranslation(ch))
                    rot = helper.getMatRotation(ch)
                    liste_particles_coords.append(coord)
                    output.append([coord,rot.tolist(),group_name,1,0])
        indice_group=0
        N=len(liste_particles_coords)
        g = helper.particle("PS_"+group_name,liste_particles_coords,group_name="PS_"+group_name)#,radius=liste_particles_radius)
        #fluogrp = helper.createGroup("PSfluo",[1.,1.,1.])
        #g = helper.createGroup("PS_"+group_name,color=[float(i)/float(indice_group),0,0],PS=PS)
        tpg = helper.newTPgeometry("gPS_"+group_name,group = g)
        #array is ( jtrans, rotMatj, self, ptInd )
        numpy.savetxt("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/"+group_name+"_pos.txt",numpy.array(liste_particles_coords),delimiter=",")
        #numpy.savetxt("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/"+group_name+"_autopack.txt",numpy.array(output),delimiter=",")
        with open("/Users/ludo/Downloads/autopack_ouput/hiv_experiment/"+group_name+"_autopack.json", 'w') as fp :#doesnt work with symbol link ?
            json.dump(output,fp)
        #map(h.ps.SetGroup,range(N),liste_particles_group_id[:N])
        #helper.setParticlProperty("group",range(N),[g,]*N,PS=PS)
        print "OK"
    else :
        #update ?
        pass