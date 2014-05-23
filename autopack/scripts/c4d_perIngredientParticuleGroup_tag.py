# -*- coding: utf-8 -*-
"""
Created on Fri Apr 4 2014

@author: Ludovic Autin
"""
import c4d
import autopack
#Welcome to the world of Python
h=c4d.af.values()[0].histoVol.values()[0]
helper = autopack.helper 
def main():
    helper.doc = doc
    PS = doc.GetParticleSystem()
    ob = op.GetObject()
    mx = ob.GetMg()
    cframe = doc.GetTime().GetFrame(doc.GetFps())
    doit = True
    if hasattr(h,"has_particle") :
        doit = not h.has_particle
    if cframe == 0 and doit:
        liste_particles_coords=[]#coordinates
        liste_particles_radius=[]#coordinates
        liste_particles_group=[]#indice of the group
        group_dic={}
        group_name=[]
        indice_group=0
        h.collectResultPerIngredient()
        r =  h.exteriorRecipe
        if r :
            for ingr in r.ingredients:
                if not len(ingr.results) : continue
                group_name.append(ingr.name)
                for r in ingr.results:  
                    if hasattr(r[0],"tolist"):
                        r[0]=r[0].tolist()#position
                        liste_particles_coords.append(r[0])
                        liste_particles_group.append(indice_group)
                        liste_particles_radius.append(ingr.encapsulatingRadius)
                indice_group+=1
        #compartment ingr
        for orga in h.compartments:
            #compartment surface ingr
            rs =  orga.surfaceRecipe
            if rs :
                for ingr in rs.ingredients:
                    if not len(ingr.results) : continue
                    group_name.append(ingr.name)
                    for r in ingr.results:  
                        if hasattr(r[0],"tolist"):
                            r[0]=r[0].tolist()#position
                            liste_particles_coords.append(r[0])
                            liste_particles_group.append(indice_group)
                            liste_particles_radius.append(ingr.encapsulatingRadius)
                    indice_group+=1
            #compartment matrix ingr
            ri =  orga.innerRecipe
            if ri :
                for ingr in ri.ingredients:
                    if not len(ingr.results) : continue
                    group_name.append(ingr.name)
                    for r in ingr.results:  
                        if hasattr(r[0],"tolist"):
                            r[0]=r[0].tolist()#position
                            liste_particles_coords.append(r[0])
                            liste_particles_group.append(indice_group)
                            liste_particles_radius.append(ingr.encapsulatingRadius)
                    indice_group+=1
        N=len(liste_particles_coords)
        h.ps = helper.particle("all",liste_particles_coords,radius=liste_particles_radius)
        fluogrp = helper.createGroup("PSfluo",[1.,1.,1.])
        for i in range(indice_group):#nb of ingredient
            g = helper.createGroup("PS_"+group_name[i],color=[float(i)/float(indice_group),0,0],parent = fluogrp,PS=h.ps)
            tpg = helper.newTPgeometry("gPS_"+str(i),group = g)
            group_dic[i]=g
            #print g,i
        liste_particles_group_id =[ group_dic[i] for i in liste_particles_group ]
        #map(h.ps.SetGroup,range(N),liste_particles_group_id[:N])
        helper.setParticlProperty("group",range(N),
                                      liste_particles_group_id[:N],PS=h.ps)
        h.has_particles = True
        print "OK"
        print doit