# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 12:06:42 2014

@author: ludo
"""

import os
import numpy as np
import math
from autopack.transformation import euler_from_matrix,matrixToEuler,euler_matrix
#some progress report would be nice
#EMAN2?
class tem_sim:
    def __init__(self,name="box",bounding_box=[[0.0, 0.0, 0.], [190.0, 190.0, 190]]):
        #define whats is need for bd_box simulation
        #The following three columns are Euler angles for rotation around
        #the z axis, then around the x axis, and again around the z axis. Coordinates are in
        #nmunits, and angles in degrees.
        self.base_name=name 
        self.base_directory=os.path.dirname(name)
        self.bounding_box=np.array(bounding_box)*1.5
        self.params_file = self.base_name+"_input.txt"
        self.params_string=""
        self.order=[]
        self.params={}
        self.str_string=""
        self.str_bead_id=1
        self.bd_binary=""
        self.bond_dictionary={}
        self.bond_length_dictionary={}#should be distance between sphere
        self.setup()     
        #self.makePrmFile()

    def setup(self):
        self.objects={}#string integer
        self.objects_nb=[]#string integer
        self.objects_str_order={}
        self.initial_coordinates_str={}
        self.objects_str_order["pdb"]=["source","use_imag_pot","famp","voxel_size",
                                "pdb_file_in","map_file_re_out","map_file_im_out",
                                "map_file_format",
                                "map_file_byte_order","map_axis_order"]
        self.objects_str_order["map"]=["source","use_imag_pot","famp","nx","ny","nz",
                                "map_file_re_out","map_file_im_out",
                                "map_file_im_in","map_file_im_in","map_file_format",
                                "map_file_byte_order","map_axis_order"]
        self.objects_str_order["random"]=["source","use_imag_pot","famp","voxel_size","nx","ny","nz",
                                "contrast_re","constrast_im","smoothness","make_positive",
                                "map_file_re_out","map_file_im_out",
                                ]  
        self.objects_set_str_order=["particle_type","particle_coords",
                                "num_particles","where","coord_file_in","coord_file_out"]
        self.order_p=["simulation","sample"]
        self.order={"simulation":["generate_micrographs","generate_volumes","generate_particle_maps",
                                  "log_file","rand_seed"],
                    "sample":["diameter","thickness_center","thickness_edge",
                              "offset_x","offset_y"],
                    "geometry" : ["gen_tilt_data","tilt_axis","ntilts","theta_start",
                                  "theta_incr","geom_errors"],
                    "electronbeam":["acc_voltage","energy_spread","gen_dose",
                                    "dose_per_im"],
                    "optics":["magnification","cs","cc","aperture","focal_length",
                    "cond_ap_angle","gen_defocus","defocus_nominal"],
                    "detector":["det_pix_x","det_pix_y","pixel_size","gain",
                                "use_quantization","dqe","mtf_a","mtf_b","mtf_c",
                                "mtf_alpha","mtf_beta","image_file_out"]
                    }

        self.initial_coordinates_file=self.base_name+".txt"
        #SIMULATION
        self.params["generate_micrographs"]={"name":"generate_micrographs","value":True,"default":True,
                                           "type":"bool","description":"micrographs should be generated",
                                           "width":30}
        self.params["generate_volumes"]={"name":"generate_volumes","value":False,"default":False,
                                           "type":"bool","description":"micrographs should be generated",
                                           "width":30}
        self.params["generate_particle_maps"]={"name":"generate_volumes","value":False,"default":False,
                                           "type":"bool","description":"micrographs should be generated",
                                           "width":30}
        self.params["log_file"]={"name":"log_file","value":self.base_name+".log","default":"input.log",
                                           "type":"str","description":"a log file is saved showing all input parameters",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["rand_seed"]={"name":"rand_seed","value":12345,"default":12345,
                                           "type":"int","description":"rand_seed",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        #SAMPLE
        self.params["diameter"]={"name":"diameter","value":1200,"default":1200,
                                           "type":"float","description":"diameter",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["thickness_center"]={"name":"thickness_center","value":50,"default":50,
                                           "type":"float","description":"thickness_center",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["thickness_edge"]={"name":"thickness_edge","value":150,"default":150,
                                           "type":"float","description":"thickness_edge",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["offset_x"]={"name":"offset_x","value":0,"default":0,
                                           "type":"float","description":"offset_x",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.params["offset_y"]={"name":"offset_y","value":0,"default":0,
                                           "type":"float","description":"offset_y",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        

    def makePrmFile(self,):
        self.params_string="# To run the simulation, execute the command\n"
        self.params_string+="# <path to executable>/TEM-simulator "+self.params_file+"\n"        
#        self.params_string+="\n=== simulation ===\n"
#        self.params_string+="# The simulation component specifies what kind of computations should be done.\n"
        
        for i in self.order_p:
            self.params_string+="\n=== "+i+" ===\n"
            for k in self.order[i] :
                if type(self.params[k]) == dict :
                    if self.params[k]["type"]=="bool":
                        if self.params[k]["value"] :
                            self.params_string+=k+" = yes\n"
                        else :self.params_string+=k+" = no\n"
                    else :
                        self.params_string+=k+" = "+str(self.params[k]["value"])+"\n"
                else :
                    self.params_string+=k+" = "+str(self.params[k])+"\n"

        for obj in self.objects:
            self.params_string+="\n=== particle "+obj+" ===\n"
            for k in self.objects_str_order[self.objects[obj]["source"]["value"]]:
                    if type(self.objects[obj][k]) == dict :
                        if self.objects[obj][k]["type"]=="bool":
                            if self.objects[obj][k]["value"] :
                                self.params_string+=k+" = yes\n"
                            else :self.params_string+=k+" = no\n"
                        else :
                            self.params_string+=k+" = "+str(self.objects[obj][k]["value"])+"\n"
                    else :
                        self.params_string+=k+" = "+str(self.objects[obj][k])+"\n"
        for obj in self.objects:   
            self.params_string+="\n=== particleset ===\n"                        
            for k in self.objects_set_str_order:
                    if type(self.objects[obj][k]) == dict :
                        if self.objects[obj][k]["type"]=="bool":
                            if self.objects[obj][k]["value"] :
                                self.params_string+=k+" = yes\n"
                            else :self.params_string+=k+" = no\n"
                        else :
                            self.params_string+=k+" = "+str(self.objects[obj][k]["value"])+"\n"
                    else :
                        self.params_string+=k+" = "+str(self.objects[obj][k])+"\n"
        self.params_string+="""
=== geometry ===
gen_tilt_data = yes
tilt_axis = 0
ntilts = 1
theta_start = 0
theta_incr = 0
geom_errors = none

=== electronbeam ===
acc_voltage = 200
energy_spread = 1.3
gen_dose = yes
dose_per_im = 4000

=== optics ===
magnification = 25000
cs = 2.1
cc = 2.2
aperture = 40
focal_length = 2.7
cond_ap_angle = 0.1
gen_defocus = yes
defocus_nominal = 3.0

=== detector ===
det_pix_x = 2000
det_pix_y = 2000
pixel_size = 16
gain = 80
use_quantization = yes
dqe = 0.4
mtf_a = 0.7
mtf_b = 0.2
mtf_c = 0.1
mtf_alpha = 10
mtf_beta = 40
image_file_out = output_tem.mrc

=== detector ===
det_pix_x = 2000
det_pix_y = 2000
pixel_size = 16
gain = 80
use_quantization = no
dqe = 0.4
mtf_a = 0.7
mtf_b = 0.2
mtf_c = 0.1
mtf_alpha = 10
mtf_beta = 40
image_file_out = output_tem_nonoise.mrc
"""      

    def addGoldMarker(self):
        self.addObject("marker","random",None,30)
        
    def addObject(self,name,source,file_in,number):
        self.objects[name]={}
        
        self.objects[name]["source"]={"name":"source","value":source,"default":"pdb",
                                           "type":"str","description":"source for particle, e.g pdb,map,random",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.objects[name]["use_imag_pot"]={"name":"use_defocus_corr","value":True,"default":True,
                                           "type":"bool","description":"use_defocus_corr",
                                           "mini":0.001,"maxi":1000,
                                           "width":30}
        #if use_imag_pot is no
        self.objects[name]["famp"]={"name":"famp","value":0.3,"default":0.3,
                                           "type":"float","description":"thickness_edge",
                                           "mini":0.00,"maxi":1.0,
                                           "width":30}

        self.objects[name]["pdb_file_in"]={"name":"pdb_file_in","value":file_in,"default":self.base_name+".log",
                                           "type":"str","description":"pdb file if source is pdb",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.objects[name]["use_defocus_corr"]={"name":"use_defocus_corr","value":False,"default":False,
                                           "type":"bool","description":"use_defocus_corr",
                                           "mini":0.001,"maxi":1000,
                                           "width":30}
       # The PDB file is first converted into a volume map of the electrostatic potential
       # in the molecule. The size of the voxels in this map is set to 0.1 nm.
        self.objects[name]["voxel_size"]={"name":"voxel_size","value":0.1,"default":0.1,
                                           "type":"float","description":"thickness_edge",
                                           "mini":0.001,"maxi":100,
                                           "width":30}
        self.objects[name]["nx"]={"name":"voxel_size","value":1,"default":1,
                                           "type":"int","description":"thickness_edge",
                                           "mini":0,"maxi":100,
                                           "width":30}
        self.objects[name]["ny"]={"name":"voxel_size","value":1,"default":1,
                                           "type":"int","description":"thickness_edge",
                                           "mini":0,"maxi":100,
                                           "width":30}
        self.objects[name]["nz"]={"name":"voxel_size","value":1,"default":1,
                                           "type":"int","description":"thickness_edge",
                                           "mini":0,"maxi":100,
                                           "width":30}
 
       # The volume maps are saved to MRC files. In subsequent simulations, the maps
        # can be read directly from these files to save time.
        self.objects[name]["map_file_re_out"]={"name":"map_file_re_out","value":self.base_directory+os.sep+name+"_map.mrc","default":name+"_map.mrc",
                                           "type":"str","description":"pdb file if source is pdb",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.objects[name]["map_file_im_out"]={"name":"map_file_im_out","value":self.base_directory+os.sep+name+"_apbs_map.mrc","default":name+"_apbs_map.mrc",
                                           "type":"str","description":"pdb file if source is pdb",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.objects[name]["map_file_re_in"]={"name":"map_file_re_in","value":self.base_directory+os.sep+name+"_map.mrc","default":name+"_map.mrc",
                                           "type":"str","description":"pdb file if source is pdb",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.objects[name]["map_file_im_in"]={"name":"map_file_im_in","value":self.base_directory+os.sep+name+"_apbs_map.mrc","default":name+"_apbs_map.mrc",
                                           "type":"str","description":"pdb file if source is pdb",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.objects[name]["map_file_format"]={"name":"map_file_format","value":"mrc","default":"mrc",
                                           "type":"str","description":"pdb file if source is pdb",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.objects[name]["map_file_byte_order"]={"name":"map_file_byte_order","value":"native","default":"native",
                                           "type":"str","description":"pdb file if source is pdb",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.objects[name]["map_axis_order"]={"name":"map_axis_order","value":"xyz","default":"xyz",
                                           "type":"str","description":"map_axis_order",
                                           "mini":1,"maxi":1000,
                                           "width":30} 
        self.objects[name]["contrast_re"]={"name":"contrast_re","value":20,"default":20,
                                           "type":"float","description":"thickness_edge",
                                           "mini":0.001,"maxi":100,
                                           "width":30}
        self.objects[name]["constrast_im"]={"name":"constrast_im","value":20,"default":20,
                                           "type":"float","description":"thickness_edge",
                                           "mini":0.001,"maxi":100,
                                           "width":30}
        self.objects[name]["smoothness"]={"name":"smoothness","value":8,"default":8,
                                           "type":"float","description":"thickness_edge",
                                           "mini":0.001,"maxi":100,
                                           "width":30}
        self.objects[name]["make_positive"]={"name":"make_positive","value":True,"default":True,
                                           "type":"bool","description":"thickness_edge",
                                           "mini":0.001,"maxi":100,
                                           "width":30}
        
        self.objects[name]["particle_type"]={"name":"particle_type","value":name,"default":name,
                                           "type":"str","description":"the particle name id",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        #file,random,grid
        self.objects[name]["particle_coords"]={"name":"particle_coords","value":"random","default":"random",
                                           "type":"str","description":"particle_coords",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.objects[name]["num_particles"]={"name":"num_particles","value":number,"default":number,
                                           "type":"int","description":"pdb file if source is pdb",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        #for random generation, where volume,surface top or bottom
        self.objects[name]["where"]={"name":"where","value":"volume","default":"volume",
                                           "type":"str","description":"pdb file if source is pdb",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.objects[name]["coord_file_out"]={"name":"coord_file_out","value":self.base_directory+os.sep+name+"_coordinates.txt","default":name+"_coordinates.txt",
                                           "type":"str","description":"pdb file if source is pdb",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        self.objects[name]["coord_file_in"]={"name":"coord_file_in","value":self.base_directory+os.sep+name+"_coordinates.txt","default":name+"_coordinates.txt",
                                           "type":"str","description":"pdb file if source is pdb",
                                           "mini":1,"maxi":1000,
                                           "width":30}
        

        self.initial_coordinates_str[name]=""## File created for TEM-simulator, version 1.3.\n"

    def addAutoPackIngredient(self, ingredient):
        #name should change to three letter code for the xyz formats
        #add the object, create the mstr        
        level = ingredient.maxLevel
        if len(ingredient.results):
            #name,source,file_in,number
            self.addObject(ingredient.name,"pdb",ingredient.pdb,len(ingredient.results))
        #add the coordinate
        for r in ingredient.results:  
            if hasattr(r[0],"tolist"):
                r[0]=r[0].tolist()#position
            if hasattr(r[1],"tolist"):
                r[1]=r[1].tolist()#rotation that should be apply to sphereTree
            euler = euler_from_matrix(r[1])#Or angle from software ?
            euler = [math.degrees(euler[0]),math.degrees(euler[1]),math.degrees(euler[2])]
            self.addCoordinate(ingredient.name,r[0],euler)
            
    def addCoordinate(self,name,pos,rot):
        #position in nm, angle in degree
        self.initial_coordinates_str[name]+="%.4f %.4f %.4f %.4f %.4f %.4f\n" %(
                                pos[0]/10.0,pos[1]/10.0,pos[2]*0.0,rot[0],rot[1],rot[2])

        
    def write(self,):
        self.makePrmFile()
        f=open(self.params_file,"w")
        f.write(self.params_string)
        f.close()
        for ob in self.objects:
            print ("write TEM-simulator  coord input files for object "+ob)
            f=open(self.base_directory+os.sep+ob+"_coordinates.txt","w")
            f.write("# File created for TEM-simulator, version 1.3.\n")
            f.write(" "+str(self.objects[ob]["num_particles"]["value"])+" 6\n")
            f.write(self.initial_coordinates_str[ob])
            f.close()
        
    def setBinary(self,path):
        self.bd_binary=path
        
    def run(self,):
        #call os ? or pipe ? or thread ?
        return
#        import os
#        os.exec(self.bd_binary+" "+self.params_file)

    def setupFromEnv(self,env):
        r =  env.exteriorRecipe
        if r :
            for ingr in r.ingredients:
                self.addAutoPackIngredient(ingr)

        #compartment ingr
        for orga in env.compartments:
            #compartment surface ingr
            rs =  orga.surfaceRecipe
            if rs :
                for ingr in rs.ingredients:
                    self.addAutoPackIngredient(ingr)
            #compartment matrix ingr
            ri =  orga.innerRecipe
            if ri :
                for ingr in ri.ingredients:
                    self.addAutoPackIngredient(ingr)
                   
       


