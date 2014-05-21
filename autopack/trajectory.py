# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:50:35 2014

@author: ludo
"""

#check if MDAnalysis is available we can parse DCD file, otherwise parse the txt
import re
import autopack
import numpy as np
from numpy import savez, dtype , fromfile 
from os.path import getsize
from time import clock

MDAnalysis = None
try :
    import MDAnalysis
except :
    MDAnalysis = None
    
class Trajectory:
    traj_type=""
    reg=r"[-+]?[0-9]*\.?[0-9]+"
    def __init__(self, filename, nbIngredients):
        #need to parse
        self.data={}
        self.nconf=0
        self.filename = filename
        #self.f = open(self.filename,"r")
        self.nbIngredients=nbIngredients
        self.parse(filename=filename,nbMol=nbIngredients)
        print ("ok parsed") 
        self.current_indice = 0
        self.sub_id = 0
        self.mapping={}

    def makeIngrMapping(self,name,nbInstance):
        self.mapping[name]=[self.sub_id]
        sub_id_instance = []
        for i in range(nbInstance): #
            sub_id_instance.append(self.current_indice)
            self.current_indice += 1
        self.mapping[name].append(sub_id_instance)
        self.sub_id+=1
        
    
    def getIngredientInstancePos(self, name,instance_id,frame):
        indice = self.mapping[name][1][instance_id]
        pos = self.getPosAt(indice,frame)
        return pos
        
    #callback that should update the scene with the trajectory
    #just position....no rotation ?
    def generalApply(self,obj,obj_id,frame):
        matrice = self.getPosAt(obj_id,frame)
        autopack.helper.setObjectMatrix(autopack.helper.getName(obj),matrice)
        
    def applyState(self, h, frame):
        r =  h.exteriorRecipe
        indice = 0
        if r :
            for ingr in r.ingredients:
                for k in range(len(ingr.results)):
                    pos = self.getIngredientInstancePos(ingr.name,k,frame)
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
                        pos = self.getIngredientInstancePos(ingr.name,k,frame)
                        autopack.helper.setTranslation(ingr.ipoly[k],pos)
                        indice+=1
            #compartment matrix ingr
            ri =  orga.innerRecipe
            if ri :
                k=0
                for ingr in ri.ingredients:
                    for k in range(len(ingr.results)):
                        pos = self.getIngredientInstancePos(ingr.name,k,frame)
                        autopack.helper.setTranslation(ingr.ipoly[k],pos)
                        indice+=1

    def applyState_name(self, h, frame):
        r =  h.exteriorRecipe
        indice = 0
        if r :
            for ingr in r.ingredients:
                for k in range(len(ingr.results)):
                    pos = self.getIngredientInstancePos(ingr.name,k,frame)
                    autopack.helper.setTranslation(autopack.helper.getName(ingr.ipoly[k]),pos)
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
                        pos = self.getIngredientInstancePos(ingr.name,k,frame)
                        autopack.helper.setTranslation(autopack.helper.getName(ingr.ipoly[k]),pos)
                        indice+=1
            #compartment matrix ingr
            ri =  orga.innerRecipe
            if ri :
                k=0
                for ingr in ri.ingredients:
                    for k in range(len(ingr.results)):
                        pos = self.getIngredientInstancePos(ingr.name,k,frame)
                        autopack.helper.setTranslation(autopack.helper.getName(ingr.ipoly[k]),pos)
                        indice+=1

    def applyState_primitive_name(self, h, frame):
        r =  h.exteriorRecipe
        indice = 0
        if r :
            for ingr in r.ingredients:
                if ingr.Type == "SingleSphere" or ingr.Type == "SingleCube":
                    for k in range(len(ingr.results)):
                        pos = self.getIngredientInstancePos(ingr.name,k,frame)
                        autopack.helper.setTranslation(autopack.helper.getName(ingr.ipoly[k]),pos)
                        indice+=1
                elif ingr.Type == "MultiSphere":
                    level = ingr.maxLevel
                    nbPrim = len(ingr.radii[level])
                    for k in range(len(ingr.results)*nbPrim):
                        pos = self.getIngredientInstancePos(ingr.name,k,frame)
                        autopack.helper.setTranslation(autopack.helper.getName(ingr.isph[k]),pos)
                        indice+=1                
        for orga in h.compartments:
            #compartment surface ingr
            rs =  orga.surfaceRecipe
            if rs :
                for ingr in rs.ingredients:
                    if ingr.Type == "SingleSphere" or ingr.Type == "SingleCube":
                        for k in range(len(ingr.results)):
                            pos = self.getIngredientInstancePos(ingr.name,k,frame)
                            autopack.helper.setTranslation(autopack.helper.getName(ingr.ipoly[k]),pos)
                            indice+=1
                    elif ingr.Type == "MultiSphere":
                        level = ingr.maxLevel
                        nbPrim = len(ingr.radii[level])
                        for k in range(len(ingr.results)*nbPrim):
                            pos = self.getIngredientInstancePos(ingr.name,k,frame)
                            autopack.helper.setTranslation(autopack.helper.getName(ingr.isph[k]),pos)
                            indice+=1
            #compartment matrix ingr
            ri =  orga.innerRecipe
            if ri :
                k=0
                for ingr in ri.ingredients:
                    if ingr.Type == "SingleSphere" or ingr.Type == "SingleCube":
                        for k in range(len(ingr.results)):
                            pos = self.getIngredientInstancePos(ingr.name,k,frame)
                            autopack.helper.setTranslation(autopack.helper.getName(ingr.ipoly[k]),pos)
                            indice+=1
                    elif ingr.Type == "MultiSphere":
                        level = ingr.maxLevel
                        nbPrim = len(ingr.radii[level])
                        for k in range(len(ingr.results)*nbPrim):
                            pos = self.getIngredientInstancePos(ingr.name,k,frame)
                            autopack.helper.setTranslation(autopack.helper.getName(ingr.isph[k]),pos)
                            indice+=1

    def applyState_cb(self, h, frame,cb):
        r =  h.exteriorRecipe
        indice = 0
        if r :
            for ingr in r.ingredients:
                for k in range(len(ingr.results)):
                    pos = self.getIngredientInstancePos(ingr.name,k,frame)
                    cb(ingr,pos)
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
                        pos = self.getIngredientInstancePos(ingr.name,k,frame)
                        cb(ingr,pos)
                        indice+=1
            #compartment matrix ingr
            ri =  orga.innerRecipe
            if ri :
                k=0
                for ingr in ri.ingredients:
                    for k in range(len(ingr.results)):
                        pos = self.getIngredientInstancePos(ingr.name,k,frame)
                        cb(ingr,pos)
                        indice+=1
                        
    def completeMapping(self,h):
        self.current_indice = 0
        self.sub_id = 0
        self.mapping={}
        r =  h.exteriorRecipe
        if r :
            for ingr in r.ingredients:
                if not len(ingr.results) : continue
                if ingr.Type == "SingleSphere" or ingr.Type == "SingleCube" or self.traj_type=="molb":
                    self.makeIngrMapping(ingr.name,len(ingr.results))
                elif ingr.Type == "MultiSphere":
                    level = ingr.maxLevel
                    nbPrim = len(ingr.radii[level])
                    self.makeIngrMapping(ingr.name,len(ingr.results)*nbPrim)
                    
        #compartment ingr
        for orga in h.compartments:
            #compartment surface ingr
            rs =  orga.surfaceRecipe
            if rs :
                for ingr in rs.ingredients:
                    if not len(ingr.results) : continue
                    if ingr.Type == "SingleSphere" or ingr.Type == "SingleCube"  or self.traj_type=="molb":
                        self.makeIngrMapping(ingr.name,len(ingr.results))
                    elif ingr.Type == "MultiSphere":
                        level = ingr.maxLevel
                        nbPrim = len(ingr.radii[level])
                        self.makeIngrMapping(ingr.name,len(ingr.results)*nbPrim)
            #compartment matrix ingr
            ri =  orga.innerRecipe
            if ri :
                for ingr in ri.ingredients:
                    if not len(ingr.results) : continue
                    if ingr.Type == "SingleSphere" or ingr.Type == "SingleCube" or self.traj_type=="molb":
                        self.makeIngrMapping(ingr.name,len(ingr.results))
                    elif ingr.Type == "MultiSphere":
                        level = ingr.maxLevel
                        nbPrim = len(ingr.radii[level])
                        self.makeIngrMapping(ingr.name,len(ingr.results)*nbPrim)
                                        
class dcdTrajectory(Trajectory): 
    traj_type="dcd"
    dcd=None
    def parse(self,filename=None,**kwds):
        if MDAnalysis == None :
            return
        if filename==None:
            filename = self.filename
        self.dcd = MDAnalysis.coordinates.DCD.DCDReader(filename)

    def getPosAt(self,indice,frame):
        return self.dcd[frame][indice]
                
class xyzTrajectory(Trajectory): 
    #NMR Type file type
    traj_type="xyz"
    def parse(self,filename=None,**kwds):
        #brutforce parsing
        if filename==None:
            filename = self.filename
        self.f = open(self.filename,"r")
        self.line_offset = []
        self.nLines=0
        offset = 0
        for line in self.f:
            self.line_offset.append(offset)
            offset += len(line)
            self.nLines+=1
        # Now, to skip to line n (with the first line being line 0), just do
        #file.seek(line_offset[n])
        #'ATOM     39 ext__i200_n50 sub    39     -21.404-192.392  -1.462 -5.00 25.00\n'
        self.nconf = self.nLines/float(self.nbIngredients+2)
             
    def getPosAt(self,indice,frame):
        real_indice = (self.nbIngredients+3)*frame + indice
        return self.getPosLine(real_indice)
        
    def getPosLine(self, line):
        self.f.seek(self.line_offset[line])
        line = self.f.readline()
        matches = re.findall(self.reg, line)
        return [matches[-5],matches[-4],matches[-3]]

class molbTrajectory(Trajectory): 
    #NMR Type file type
    traj_type="molb"
    def parse(self,filename=None,nbMol=0,log=False):
        #brutforce parsing
        if filename==None:
            filename = self.filename
        #for all molecules grab the result ? nbMol is the total number of molecules
        #but how make difference between 2 molecules    
        #dicrionary indice in moldb-> indice in input ?
        for i in range(nbMol):
           mat = self.parse_one_mol(filename,i,log=log)
           self.data[i]=mat
               
    def parse_one_mol(self,filename,instance_id, ftype="double",log=False ):
        #start_time = clock()
        #progresss bard ? but it slow down
        if filename==None:
            filename = self.filename
        f=open(filename,"rb")
        nr=0;
        step_nr=0
        mols=0
        state_size=0
        b = [0,0,0,0];
        size = 8;
        liste_m=[]
        nr = instance_id
        size = 8;
        while 1:
            m = np.zeros((4,4))
            b=fromfile(f,'<i',count=4)
            if (log) : print b
            if not len(b) : break
            step_nr = fromfile(f,'<i',count=1)[0]
            mols = fromfile(f,'<i',count=1)[0]
            if autopack.helper is not None :
                autopack.helper.progressBar(label="parsin moldb "+str(step_nr)+" "+str(mols)+" "+str(nr))
            if ( mols <= nr ):
                print("Number of molecules (%d) is less of equal than the number given (%d)\n"% (mols,nr));
            if mols == 0 :
                break
            state_size= fromfile(f,'<i',count=1)[0] 
            f.seek(state_size,1)
            f.seek(size,1)            
            f.seek(size*nr*4,1)
            m[3][:3] = fromfile(f,'<d',count=4)[:3]#size*4 4 double ?
            f.seek( size*(mols-1-nr)*4, 1 )
            f.seek( size*nr*9, 1 )
            m[:3,:3] = fromfile(f,'<d',count=9).reshape(3,3)#.transpose()#size*4 4 double ?
            f.seek( size*(mols-1-nr)*9, 1)  
            liste_m.append(m)
            if (log) : print step_nr,m,mols
        f.close()
        #end_time = clock()
        #print 'It took',end_time - start_time,'seconds'
        return np.array(liste_m)
             
    def getPosAt(self,indice,frame):        
        return self.data[indice][frame]

    def applyState_name(self, h, frame):
        r =  h.exteriorRecipe
        indice = 0
        if r :
            for ingr in r.ingredients:
                for k in range(len(ingr.results)):
                    pos = self.getIngredientInstancePos(ingr.name,k,frame)
                    autopack.helper.setTransformation(autopack.helper.getName(ingr.ipoly[k]),mat=pos)
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
                        pos = self.getIngredientInstancePos(ingr.name,k,frame)
                        #pos
                        newp = np.identity(4)
                        newp[3,:3] = pos[3,:3]
                        newp[:3,:3] = pos.transpose()[:3,:3]
                        autopack.helper.setTransformation(autopack.helper.getName(ingr.ipoly[k]),mat=pos)
                        indice+=1
            #compartment matrix ingr
            ri =  orga.innerRecipe
            if ri :
                k=0
                for ingr in ri.ingredients:
                    for k in range(len(ingr.results)):
                        pos = self.getIngredientInstancePos(ingr.name,k,frame)
                        autopack.helper.setTransformation(autopack.helper.getName(ingr.ipoly[k]),mat=pos)
                        indice+=1

    def applyState_primitive_name(self, h, frame):
        return self.applyState_name(h,frame)          
#from autopack.trajectory import dcdTrajectory
#dcd = dcdTrajectory("/Users/ludo/Downloads/autopack_ouput/spheres.dcd",len(h.result))