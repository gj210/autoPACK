# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 22:07:48 2014

@author: ludo
"""
import os
import numpy as np
import math
from autopack.transformation import euler_from_matrix,matrixToEuler
#some progress report would be nice

class bd_box:
    def __init__(self,name="box",bounding_box=[[0.0, 0.0, 0.], [190.0, 190.0, 190]]):
        #define whats is need for bd_box simulation
        self.base_name=name   
        self.bounding_box=bounding_box
        self.params_file = self.base_name+".prm"
        self.params_string=""
        self.order=[]
        self.params={}
        self.str_string=""
        self.str_bead_id=1
        self.bd_binary=""
        self.bond_dictionary={}
        self.bond_length_dictionary={}#should be distance between sphere
        self.setup()     
        self.makePrmFile()

    def makePrmFile(self,):
        self.params_string=""
        for k in self.order :
            if type(self.params[k]) == dict :
                if self.params[k]["type"]=="bool":
                    if self.params[k]["value"] :
                        self.params_string+=k+" yes\n"
                    else :self.params_string+=k+" no\n"
                else :
                    self.params_string+=k+" "+str(self.params[k]["value"])+"\n"
            else :
                self.params_string+=k+" "+str(self.params[k])+"\n"

    def write(self,):
        f=open(self.params_file,"w")
        f.write(self.params_string)
        f.close()
        
    def setBinary(self,path):
        self.bd_binary=path
        
    def run(self,):
        #call os ? or pipe ? or thread ?
        return
#        import os
#        os.exec(self.bd_binary+" "+self.params_file)
    
        
class flex_box(bd_box):
    def setup(self,):
        #define whats is need for bd_box simulation
        self.order=["dt","T","visc","vfactor","bdsteps","save_dcd_freq","save_rst_freq",
                    "save_enr_freq","save_xyz_freq","dcd_filename","enr_filename","str_filename",
                    "out_filename","pqr_filename","xyz_filename","alpha_lj","lj_6_term","elec",
                    "gamma_c","kappa_c","move_attempts","cutoff_lj","cutoff_c",
                    "bond_lj_scale","bond_c_scale","check_overlap","bc","xbox",
                    "ybox","zbox","E_ext","E_magn","E_factor","E_dir1","E_dir2",
                    "E_freq","E_type","ewald_recip","ewald_real","ewald_method",
                    "cuda_block","cuda_devices","algorithm","rand_seed","hydro",
                    "epsilon_c","e_collision","nb_list"]
        self.params["dt"] = 10.0
        self.params["T"] =  298.15
        self.params["visc"] =  0.0102
        self.params["vfactor"] =  14.4
        self.params["bdsteps"] =  50000
        self.params["save_dcd_freq"] =  10
        self.params["save_rst_freq"] =  10
        self.params["save_enr_freq"] =  10
        self.params["save_xyz_freq"] =  10
        self.params["dcd_filename"] =  self.base_name+".dcd"
        self.params["enr_filename"] =  self.base_name+".enr"
        self.params["str_filename"] =  self.base_name+".str"
        self.params["out_filename"] =  self.base_name+".out"
        self.params["pqr_filename"] =  self.base_name+".pqr"
        self.params["xyz_filename"] =  self.base_name+".xyz"
        self.params["alpha_lj"] =  4.0
        self.params["lj_6_term"] =  "yes"
        self.params["elec"]= "yes" 
        self.params["gamma_c"] =  332.4
        self.params["kappa_c"] =  0.1
        self.params["epsilon_c"] =  78.54
        self.params["move_attempts"] =  1000000
        self.params["cutoff_lj"] =  -1
        self.params["cutoff_c"] =  200.0
        self.params["bond_lj_scale"] =  1.0
        self.params["bond_c_scale"] =  0.0
        self.params["check_overlap"] =  "yes"
        self.params["bc"] =  "pbc"
        self.params["xbox"] =  self.bounding_box[1][0]-self.bounding_box[0][0]
        self.params["ybox"] =  self.bounding_box[1][1]-self.bounding_box[0][1]
        self.params["zbox"] =  self.bounding_box[1][2]-self.bounding_box[0][2]
        self.params["E_ext"] =  "yes"
        self.params["E_magn"] =  1e+08
        self.params["E_factor"] =  0.2e-08
        self.params["E_dir1"] =  "x"
        self.params["E_dir2"] =  "y"
        self.params["E_freq"] =  0.01
        self.params["E_type"] =  "RF"
        self.params["ewald_recip"] =  0#3
        self.params["ewald_real"] =  0#3
        self.params["ewald_method"] =  "smith"
        self.params["cuda_block"] =  32
        self.params["cuda_devices"] =  "1 0"
        self.params["algorithm"] =  "ermak_const"
        self.params["rand_seed"] =  7239847
        self.params["hydro"] =  "cholesky"
        self.params["e_collision"] =  1#1 (perfectly elastic collisions) 0 (perfectly inelastic collisions)
        self.params["nb_list"] =  "spatial"
#        self.bd_binary="/Users/ludo/DEV/BDBOX/bin/bd_flex"  #v1.0
        self.bd_binary="/Users/ludo/DEV/BD_BOX2/bin/bd_flex" #v2.0
        
    def addSubBead(self,name,x,y,z,t,Q,R,LJ,m,bid=None,bond=False,bond_id=0,):
        #name, x, y, z, Q, LJ, m, bid=None
        ##sub name id x y z t Q 2R LJ m
        #sub ext__ix400_n25 11 1024.1435 804.5214 700.7810 25.0000 -5.0000 50.0000 0.5922 1.0000
        #no negative position!
        if bid == None :
            bid = self.str_bead_id
            self.str_bead_id+=1
        self.str_string+="sub %s %d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n" %(name,
                                bid, x,   y,   z,    t,   Q,  R,  LJ, m)
        if bond :
            if name not in self.bond_dictionary:
                self.bond_dictionary[name]={}
            if bond_id not in self.bond_dictionary[name]:
                self.bond_dictionary[name][bond_id]={}
                self.bond_dictionary[name][bond_id]["id"]=[]
                self.bond_dictionary[name][bond_id]["pos"]=[]
            self.bond_dictionary[name][bond_id]["id"].append(bid)
            self.bond_dictionary[name][bond_id]["pos"].append([x,   y,   z])
            
    def addAutoPackIngredient(self, ingredient):
        #name should change to three letter code for the xyz formats
        bi=0
        for r in ingredient.results:  
            if hasattr(r[0],"tolist"):
                r[0]=r[0].tolist()#position
            if hasattr(r[1],"tolist"):
                r[1]=r[1].tolist()#rotation that should be apply to sphereTree 
            if ingredient.Type == "SingleSphere" or ingredient.Type == "SingleCube":
                self.addSubBead(ingredient.name,r[0][0],r[0][1],r[0][2],ingredient.encapsulatingRadius,-5.0,
                            ingredient.encapsulatingRadius,0.5922,1)
            elif ingredient.Type == "MultiSphere":
                #one line per sphere, position came from rotation+need bond
                level = ingredient.maxLevel
                px = ingredient.transformPoints(r[0],r[1], ingredient.positions[level])
                for ii in range(len(ingredient.radii[level])):
                    self.addSubBead(ingredient.name,px[ii][0],px[ii][1],px[ii][2],ingredient.radii[level][ii],-5.0,
                                ingredient.radii[level][ii],0.5922,1,bond=True,bond_id=bi) 
                bi+=1
                    
    def makeBond(self):
        #bond id(1) id(2) ro rmax H
        #eq_bond_length=0
        #max_bond_length=0
        H=1
        for iname in self.bond_dictionary:
            for bondid in self.bond_dictionary[iname]:
                liste_bond = self.bond_dictionary[iname][bondid]["id"]
                xv, yv = np.meshgrid(liste_bond,liste_bond)
                for i in range(len(liste_bond)-1):
                    for j in range(1,len(liste_bond)):
                        # treat xv[j,i], yv[j,i]
                        if xv[j,i] == yv[j,i] :
                            continue
                        p1=self.bond_dictionary[iname][bondid]["pos"][i]
                        p2=self.bond_dictionary[iname][bondid]["pos"][j]
                        d = math.sqrt((p2[0] - p1[0]) ** 2 +
                                    (p2[1] - p1[1]) ** 2 +
                                    (p2[2] - p1[2]) ** 2)
                        self.str_string+="bond %d %d %.4f %.4f %.4f\n" %(
                                        xv[j,i],yv[j,i],d,10000000.0000,H)

    def write(self,):
        bd_box.makePrmFile()
        self.makeBond()        
        f=open(self.params["str_filename"],"w")
        f.write(self.str_string)
        f.close()
    
class rigid_box(bd_box):
    def setup(self):
        self.order=["dt","min_dt","move_attempts","T","bdsteps","save_enr_freq",
        "save_dcd_freq","save_molb_freq","save_rst_freq","rst_filename",
        "molb_filename","dcd_filename","enr_filename","out_filename","pqr_filename",
        "alpha","repul_B","lj_attract","lj_repul","elec","lj_cutoff","beta",
        "eps","epsin","gamma","kappa","check_overlaps","overlaps","overlaps_removal",
        "points_per_sphere","probe_radius","srf_ratio","hdb_surf","hdb_cutoff",
        "hdb_a","hdb_b","hdb_phi","elec_cutoff","elec_desolv","bc","xbox",
        "ybox","zbox","bd_algorithm","benchmark","rand_seed","nb_list","pos_min","ext_coor_file"]

        self.initial_coordinates_file=self.base_name+".txt"
        self.params["out_filename"]={"name":"out_filename","value":self.base_name+".log","default":"err.log",
                                           "type":"str","description":"plain text log file (output)",
                                           "width":30}
        self.params["dt"]={"name":"dt","value":5,"default":5,
                                           "type":"int","description":"timestep",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["T"]={"name":"dt","value":298.15,"default":298.15,
                                           "type":"float","description":"temperature",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["eps"]={"name":"eps","value":78.54,"default":78.54,
                                           "type":"float","description":"solvent dielectric constant",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["epsin"]={"name":"epsin","value":1,"default":1,
                                           "type":"float","description":"solute dielectric constant",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["kappa"]={"name":"kappa","value":0.1,"default":0.1,
                                           "type":"float","description":"inverse of the Debye screening length",
                                           "mini":1,"maxi":1000,
                                           "width":30} 
        self.params["elec_cutoff"]={"name":"elec_cutoff","value":45.0,"default":45.0,
                                           "type":"float","description":"cutoff for electrostatic interactions",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["bdsteps"]={"name":"bdsteps","value":45000000,"default":45000000,
                                           "type":"int","description":"number of simulation's steps",
                                           "mini":1,"maxi":1000,
                                           "width":30} 
        self.params["gamma"]={"name":"gamma","value":331.842,"default":331.842,
                                           "type":"float","description":"scaling factor for electrostatic interactions",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["alpha"]={"name":"alpha","value":4,"default":4,
                                           "type":"int","description":"scaling factor for L-J interactions",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["beta"]={"name":"beta","value":1,"default":1,
                                           "type":"int","description":"scaling factor for cavity (desolvation) terms",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["lj_cutoff"]={"name":"lj_cutoff","value":10,"default":10,
                                           "type":"int","description":"cutoff for L-J interactions",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["lj_attract"]={"name":"lj_attract","value":6,"default":6,
                                           "type":"int","description":"power, attractive term in the L-J potential",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["lj_repul"]={"name":"lj_repul","value":12,"default":12,
                                           "type":"int","description":"power, repulsive term in the L-J potential",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["save_dcd_freq"]={"name":"save_dcd_freq","value":100000,"default":1,
                                           "type":"int","description":"frequency (number of steps) for writing to the dcd trajectory",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["save_rst_freq"]={"name":"save_rst_freq","value":500,"default":1,
                                           "type":"int","description":"frequency (number of steps) for writing to the restart file",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["save_enr_freq"]={"name":"save_enr_freq","value":50,"default":1,
                                           "type":"int","description":"frequency (number of steps) for writing to the energy file",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["save_molb_freq"]={"name":"save_molb_freq","value":50,"default":1,
                                           "type":"int","description":"frequency (number of steps) for writing to the molb trajectory",
                                           "mini":1,"maxi":1000,
                                           "width":30} 
        self.params["xbox"]={"name":"xbox","value":self.bounding_box[1][0]-self.bounding_box[0][0],
                                            "default":0,
                                           "type":"float","description":"computational box, x-size",
                                           "mini":1,"maxi":10000,
                                           "width":30} 
        self.params["ybox"]={"name":"ybox","value":self.bounding_box[1][1]-self.bounding_box[0][1],
                                            "default":0,
                                           "type":"float","description":"computational box, y-size",
                                           "mini":1,"maxi":10000,
                                           "width":30}
        self.params["zbox"]={"name":"zbox","value":self.bounding_box[1][2]-self.bounding_box[0][2],
                                            "default":0,
                                           "type":"float","description":"computational box, z-size",
                                           "mini":1,"maxi":10000,
                                           "width":30}
        self.params["elec"]={"name":"elec","value":True,"default":True,
                                           "type":"bool","description":"switch electrostatic interactions on/off",
                                           "width":30}
        self.params["bd_algorithm"]={"name":"bd_algorithm","value":"ermak_const","default":"ermak_const",
                                           "type":"str","description":"BD algorithm to be used",
                                           "mini":1,"maxi":10000,
                                           "width":30}#list ?
        self.params["dcd_filename"]={"name":"dcd_filename","value":self.base_name+".dcd","default":"",
                                           "type":"str","description":"dcd–formatted trajectory file (output)",
                                           "mini":1,"maxi":10000,
                                           "width":30}
        self.params["enr_filename"]={"name":"enr_filename","value":self.base_name+".enr","default":"",
                                           "type":"str","description":"plain-text file with energies (output)",
                                           "mini":1,"maxi":10000,
                                           "width":30} 
        self.params["rst_filename"]={"name":"rst_filename","value":self.base_name+".rst","default":"",
                                           "type":"str","description":"name of the restart file",
                                           "mini":1,"maxi":10000,
                                           "width":30}
        self.params["pqr_filename"]={"name":"pqr_filename","value":self.base_name+".pqr","default":"",
                                           "type":"str","description":"pqr-formatted template file (output)",
                                           "mini":1,"maxi":10000,
                                           "width":30}
        self.params["molb_filename"]={"name":"molb_filename","value":self.base_name+".molb","default":"",
                                           "type":"str","description":"molb–formatted trajectory file (output)",
                                           "mini":1,"maxi":10000,
                                           "width":30}
        self.params["rand_seed"]={"name":"rand_seed","value":12345,"default":12345,
                                           "type":"int","description":"random generator seed",
                                           "mini":1,"maxi":100000,
                                           "width":30}
        self.params["E_ext"]={"name":"E_ext","value":False,"default":False,
                                           "type":"bool","description":"yes/no - switch the external electric field on/off",
                                           "width":30}
        self.params["E_magn"]={"name":"E_magn","value":1,"default":1,
                                           "type":"int","description":"magnitude of the external electric field",
                                           "mini":1,"maxi":100000,
                                           "width":30}
        self.params["E_type"]={"name":"E_type","values":["DC","AC","RF"],"value":"DC","default":"DC",
                                           "type":"liste","description":"DC/AC/RF - a choice between AC, DC or RF electric fields",
                                           "width":30} 
        self.params["E_freq"]={"name":"E_freq","value":0,"default":0,
                                           "type":"int","description":"frequency of the external electric field, applicable in case of AC or RF",
                                           "mini":1,"maxi":100000,
                                           "width":30}
        self.params["E_dir1"]={"name":"E_dir1","values":["x","y","z"],"value":"x","default":"x",
                                           "type":"liste","description":"x/y/z - the direction in which the external electric field is applied",
                                           "width":30} 
        self.params["E_dir2"]={"name":"E_dir2","values":["x","y","z"],"value":"y","default":"y",
                                           "type":"liste","description":"x/y/z - the second direction in which the external electric field (only for RF) is applied",
                                           "width":30}
        self.params["E_factor"]={"name":"E_factor","value":1,"default":1,
                                           "type":"int","description":"electric field, units conversion factor",
                                           "mini":1,"maxi":100000,
                                           "width":30}
        self.params["bc"]={"name":"bc","value":"pbc","default":"pbc",
                                           "type":"str","description":"boundary conditions, currently none, a rectangular box or a reflective sphere",
                                           "width":30}
        self.params["sphere_radius"]={"name":"sphere_radius","value":0,"default":0,
                                           "type":"float","description":"radius of the reflective sphere enclosing the studied system",
                                           "mini":1,"maxi":100000,
                                           "width":30}
        self.params["restart"]={"name":"restart","value":"","default":"",
                                           "type":"str","description":"apply restart from filename",
                                           "width":30} 
        self.params["move_attempts"]={"name":"move_attempts","value":10,"default":10,
                                           "type":"int","description":"number of attempts to draw a new random vector in case of an overlap between particles",
                                           "mini":1,"maxi":100000,
                                           "width":30}
        self.params["cuda_devices"]={"name":"cuda_devices","value":"0 1","default":"0 1",
                                           "type":"str","description":"Cuda devices to be used",
                                           "width":30}
        self.params["sboundary"]={"name":"sboundary","value":False,"default":False,
                                           "type":"bool","description":"yes/no - switch the bounding sphere field on/off",
                                           "width":30}
        self.params["sboundary_A"]={"name":"sboundary_A","value":1,"default":1,
                                           "type":"int","description":"the magnitude of the bounding sphere force",
                                           "mini":1,"maxi":100000,
                                           "width":30} 
        self.params["sboundary_n"]={"name":"sboundary_n","value":2,"default":2,
                                           "type":"int","description":"the power of the radial distance dependence of the bounding sphere force",
                                           "mini":1,"maxi":100000,
                                           "width":30}
        self.params["sboundary_cutoff"]={"name":"sboundary_cutoff","value":0,"default":0,
                                           "type":"int","description":"the bounding sphere force is applied outside this cutoff radius",
                                           "mini":1,"maxi":100000,
                                           "width":30}
        self.params["check_overlaps"]={"name":"check_overlaps","value":True,"default":True,
                                           "type":"bool","description":"whether to detect overlaps after each step",
                                           "width":30} 
        self.params["overlaps"]={"name":"overlaps","value":"trees","default":"trees",
                                           "type":"str","description":"an algorithm used to check for overlaps",
                                           "width":30}
        self.params["overlaps_removal"]={"name":"overlaps_removal","value":"none","default":"none",
                                           "type":"str","description":"whether to remove overlaps between particles using an iterative procedure, yes/no",
                                           "width":30} 
        self.params["cuda_block"]={"name":"cuda_block","value":256,"default":256,
                                           "type":"int","description":"The dimension of the thread block",
                                           "mini":1,"maxi":2048,
                                           "width":30}
        self.params["elec_desolv"]={"name":"elec_desolv","value":True,"default":True,
                                           "type":"bool","description":"switch electrostatic interactions - desolvation forces on/off",
                                           "width":30}
        self.params["min_dt"]={"name":"min_dt","value":-1,"default":-1,
                                           "type":"int","description":"minimal time step in the variable-step BD algorithm",
                                           "mini":-1,"maxi":100000,
                                           "width":30} 
        self.params["probe_radius"]={"name":"probe_radius","value":1.4,"default":1.4,
                                           "type":"float","description":"spherical probe with a given radius",
                                           "mini":1,"maxi":10000,
                                           "width":30} 
        self.params["srf_ratio"]={"name":"srf_ratio","value":0.9,"default":0.9,
                                           "type":"float","description":"predefined threshold",
                                           "mini":0.001,"maxi":10000,
                                           "width":30}
        self.params["points_per_sphere"]={"name":"points_per_sphere","value":1000,"default":1000,
                                           "type":"int","description":"number of a set of points",
                                           "mini":1,"maxi":10000,
                                           "width":30} 
        self.params["nb_list"]={"name":"nb_list","value":"brute","default":"trees",
                                           "type":"str","description":"generation of nonbonded interactions lists algorithm",
                                           "width":30}
        self.params["hdb_surf"]={"name":"hdb_surf","value":"none","default":"none",
                                           "type":"str","description":"whether to evaluate hydrophobic interactions",
                                           "width":30} 
        self.params["hdb_cutoff"]={"name":"hdb_cutoff","value":10,"default":10,
                                           "type":"float","description":"cutoff for hydrophobic (SASA based) interactions",
                                           "mini":0.001,"maxi":10000,
                                           "width":30}
        self.params["pos_iter"]={"name":"pos_iter","value":100,"default":100,
                                           "type":"float","description":"number of iterations during which the program tries to place randomly oriented molecules inside a primary simulation cell",
                                           "mini":0.001,"maxi":10000,
                                           "width":30} 
        self.params["pos_min"]={"name":"pos_min","value":0,"default":0,
                                           "type":"float","description":"min distance beetwen surfaces of molecules",
                                           "mini":0,"maxi":10000,
                                           "width":30}
        self.params["ext_coor_file"]={"name":"ext_coor_file","value":self.initial_coordinates_file,"default":self.base_name+".txt",
                                           "type":"str","description":"text file with coordinates of molecules that can be used to initatiate a simulation",
                                           "width":30}
#        self.params["ext_coor_molb"]={"name":"ext_coor_molb","value":"","default":"",
#                                           "type":"str","description":"molb file that can be used to initiate a simulation",
#                                           "width":30}
#        self.params["ext_coor_molb_step"]={"name":"ext_coor_molb_step","value":0,"default":0,
#                                           "type":"int","description":"number of frame to read",
#                                           "mini":0,"maxi":10000,
#                                           "width":30} 
        self.params["repul_B"]={"name":"repul_B","value":1,"default":1,
                                           "type":"int","description":"magnitude of the repulsive part in the L-J potential",
                                           "mini":0,"maxi":10000,
                                           "width":30}
        self.params["dm_alg"]={"name":"dm_alg","value":"atomdec","default":"atomdec",
                                           "type":"str","description":"type of distributed memory algorithm for computing nonbonded interactions",
                                           "width":30}
        self.params["rep"]={"name":"rep","values":["lj","dcd","pqr","LJ","Q"],"value":"lj","default":"lj",
                                           "type":"liste","description":"lj,dcd and pqr molecular representation, LJ or Q from mstr-file",
                                           "width":30}
        self.params["verlet_count"]={"name":"verlet_count","value":1,"default":1,
                                           "type":"int","description":"frequency of nonbonded neighbor list regeneration",
                                           "mini":0,"maxi":10000,
                                           "width":30}
        self.params["verlet_roff"]={"name":"verlet_roff","value":10,"default":10,
                                           "type":"int","description":"outer cutoff in the double-cutoff scheme",
                                           "mini":0,"maxi":10000,
                                           "width":30}
        self.params["benchmark"]={"name":"benchmark","value":False,"default":False,
                                           "type":"bool","description":"flag used to activate forces benchmarking",
                                           "width":30} 
        self.params["hdb_phi"]={"name":"hdb_phi","value":0.5,"default":0.5,
                                           "type":"float","description":"",
                                           "mini":0,"maxi":10000,
                                           "width":30} 
        self.params["hdb_a"]={"name":"hdb_a","value":3.1,"default":3.1,
                                           "type":"float","description":"",
                                           "mini":0,"maxi":10000,
                                           "width":30} 
        self.params["hdb_b"]={"name":"hdb_b","value":4.35,"default":4.35,
                                           "type":"float","description":"",
                                           "mini":0,"maxi":10000,
                                           "width":30}
        
        #object string integer - the name of the object (molecule) to be simulated and the number of objects of this
#kind in the system - the program will look for appropriate string.mstr and string.dt 
#les; multiple object lines
##can be present in the parameter 
#le - their order determines the ordering of molecules in all output 
#les as well
##as the order of objects in the external coordinate 
#le de
#ned with the ext coor 
#le key
        self.objects=[]#string integer
        self.objects_nb=[]#string integer
        self.objects_str_string=[]
        self.objects_dt_string=[]
        # A separate *.mstr  file is required for each kind of a molecule that 
        # is present in the simulated system
        self.mstr_files={}
        # A separate *.dt  file is required for each kind of a molecule that 
        # is present in the simulated system.
        self.dt_files={}
        self.initial_coordinates_str=""
        self.bd_binary="/Users/ludo/DEV/BD_BOX2/bin/bd_rigid"

#    def addObjToPrmFile(self,):
#       #add the object parameters
#        for i in range(len(self.objects)):
#            self.params_string+="object "+self.objects[i]+" "+self.objects_nb+"\n"
            
    def addObject(self,name,number,pos,radius):
        self.objects.append(name)
        self.objects_nb.append(number)
        aStr=""
        input_str="\n"#comment line
        input_str+="0\n" #computations are to be performed for a coarse-grained model
        input_str+=str(len(pos))+"\n" #number of beads in a model
        for i in range(len(pos)):
            aStr+=self.addBead(pos[i],radius[i])
            input_str+="%f %f %f %f\n" % (pos[i][0],pos[i][1],pos[i][2],radius[i])
        input_str+="0.01 298.15\n"#viscosity temperature - viscosity of the solvent (Poisse) and temperature (K)
        input_str+="0.0\n"       #solvent radius - radius of solvent molecules (typically set to 0.0A in case of coarse-grained models)
        self.objects_dt_string.append(input_str)
        self.objects_str_string.append(aStr)
        self.params_string+="object "+name+" "+str(number)+"\n"
        
    def addCharge(self,pos,charge,radius):
        """Q x1 y1 z1 charge1 radius1
where x, y, z denote Cartesian coordinates of a dielectric sphere in the molecular frame, 
charge denotes its central charge and radius denotes its radius. """
        aStr="Q %.4f %.4f %.4f %.4f %.4f\n" %(pos[0],pos[1],pos[2], radius, charge)
        return aStr
        
    def addBead(self,pos,radius,e=0.5922, d11=0, d21=0, c01=0, c11=0, c21=0,
                c31=0,mu1=0,SASAO1=0):
        """LJ x1 y1 z1 radius1 e1 d11 d21 c01 c11 c21 c31 mu1 SASAo1
where x, y, z denote Cartesian coordinates of spheres in the molecular frame, 
radius denotes radius used in evaluations of excluded volume and hydrophobic interactions, 
e1 is the well depth of excluded volume interactions,
parameters d and c are for the FACTS model (equations 6, 7), 
SASAo denotes the solvent accessible surface
area of a given sphere in the model (Equation 9) and 
mu1 (Equations 5, 13) denotes the surface tension parameter.
ex: LJ 5.26316 -12.7961 -4.13284 4.24053 0.5922 0 0 0 0 0 0 0 0
""" 
        aStr="LJ %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n" %(
            pos[0],pos[1],pos[2], radius,e, d11, d21, c01, c11, c21, c31, mu1, SASAO1)
        return aStr

    def addAutoPackIngredient(self, ingredient):
        #name should change to three letter code for the xyz formats
        #add the object, create the mstr        
        level = ingredient.maxLevel
        if len(ingredient.results):
            self.addObject(ingredient.name,len(ingredient.results),ingredient.positions[level],ingredient.radii[level])
        #add the coordinate
        for r in ingredient.results:  
            if hasattr(r[0],"tolist"):
                r[0]=r[0].tolist()#position
            if hasattr(r[1],"tolist"):
                r[1]=r[1].tolist()#rotation that should be apply to sphereTree
            #self.addCoordinate(r[0],euler_from_matrix(r[1]))
            self.addCoordinate(r[0],matrixToEuler(r[1]))
            
    def addCoordinate(self,pos,rot):
        self.initial_coordinates_str+="%.4f %.4f %.4f %.4f %.4f %.4f\n" %(
                                pos[0],pos[1],pos[2],rot[0],rot[1],rot[2])
        
    def writeInitialCoordinates(self,):
        """Initial coordinates of objects in the studied system (in the laboratory frame) can be speci
ed by the user in an
external text 
le, which should adhere to the following format
vx,vy,vz,omx,omy,omz (translation / rotation) for ith molecule
"""
        f=open(self.initial_coordinates_file,"w")
        f.write(self.initial_coordinates_str)
        f.close()

    def computeDmatrix(self,filename):
        executable="/Users/ludo/DEV/bd_box-2.1/tools/dmatrix/dmatrix"#comment line
        from subprocess import call
        fileName, fileExtension = os.path.splitext(filename)
        call(['mv',filename,'input.txt'])
        return_code = call([executable, filename])
        call(['mv','input.txt',filename])
        call(['mv','tensor.txt',fileName+".dt"])
        print return_code
        aStr="""  .218365E-01  .000000E+00  .000000E+00  .000000E+00  .000000E+00  .000000E+00
  .000000E+00  .218365E-01  .000000E+00  .000000E+00  .000000E+00  .000000E+00
  .000000E+00  .000000E+00  .218365E-01  .000000E+00  .000000E+00  .000000E+00
  .000000E+00  .000000E+00  .000000E+00  .163774E-03  .000000E+00  .000000E+00
  .000000E+00  .000000E+00  .000000E+00  .000000E+00  .163774E-03  .000000E+00
  .000000E+00  .000000E+00  .000000E+00  .000000E+00  .000000E+00  .163774E-03
"""
        f=open(fileName+".dt","w")
        f.write(aStr)
        f.close()
       
    def write(self,):
        bd_box.write(self)
        #write one mstr file per objects too         
        for i in range(len(self.objects)):
            f=open(self.objects[i]+".mstr","w")
            f.write(self.objects_str_string[i])
            #divides by zero, then releases the kraken!
            f.close()
            f=open(self.objects[i]+".beads","w")
            f.write(self.objects_dt_string[i])
            f.close()
            self.computeDmatrix(self.objects[i]+".beads")
        #dt file ? make on the fly?
        self.writeInitialCoordinates()
        