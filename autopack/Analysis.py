# -*- coding: utf-8 -*-
"""
Created on Mon May  6 22:58:44 2013

@author: ludo
"""
import numpy
import csv

class AnalyseAP:
    def __init__(self, env=None, viewer=None, result_file=None):
        self.env=env
        self.afviewer = viewer
        self.helper = self.afviewer.vi
        self.resutl_file = result_file

    def getPositionsFromResFile(self,):
        #could actually restore file using histoVol.
        #or not
        #need to parse apr file here anyway
        return []
        
    def getPositionsFromObject(self,parents): 
        positions=[]
        for parent in parents:
            obparent = self.helper.getObject(parent)
            childs = self.helper.getChilds(obparent)
            for ch in childs:
                ingr_name = self.helper.getName(ch)
                meshp = self.helper.getObject("Meshs_"+ingr_name.split("_")[0])
                if meshp is None :
                    c = self.helper.getChilds(ch)
                    if not len(c) :
                        continue                        
                    meshpchilds = self.helper.getChilds(c[0])#continue #should get sphere/cylnder parent ?
                else :
                    meshpchilds = self.helper.getChilds(meshp)
                for cc in meshpchilds:
                    pos = self.helper.ToVec(self.helper.getTranslation(cc))
                    positions.append(pos)
        return positions
        
    def getDistanceFrom(self,target,parents=None,**options):
        """
        target : name or host object target or target position
        parent : name of host parent object for the list of object to measre distance from
        objects : list of object or list of points
        """
        #get distance from object to the target.
        #all object are in h.molecules and orga.molecules
        #get options
        targetPos = [0,0,0]
        usePoint = False
        threshold = 99999.
        if "usePoint" in options:
            usePoint = options["usePoint"]        
        if "threshold" in options:
            threshold = options["threshold"]
        if type(target) == list or type(target) == tuple:
            targetPos = target
        elif type(target) == unicode or type(target) == str : 
            o = self.helper.getObject(target)
            if o is not None :
                targetPos = self.helper.ToVec(self.helper.getTranslation(o)) #hostForm
        else :
            o = self.helper.getObject(target)
            if o is not None :
                targetPos = self.helper.ToVec(self.helper.getTranslation(o)) #hostForm
        listeObjs=[]
        listeDistances = []
        listeCenters =[]
        if self.resutl_file is None :
            if parents is None and self.resutl_file is None:
                listeParent = [self.env.name+"_cytoplasm"]
                for o in self.env.organelles :
                    listeParent.append(o.name+"_Matrix")
                    listeParent.append(o.name+"_surface")
            elif parents is not None and self.resutl_file is None:       
                listeParent = parents
            listeCenters =self.getPositionsFromObject(listeParent)
        else :
            #use data from file
            listeCenters =self.getPositionsFromResFile(listeParent) 

        delta = numpy.array(listeCenters)-numpy.array(targetPos)
        delta *= delta
        distA = numpy.sqrt( delta.sum(1) )           
        return distA

    def getClosestDistance(self,parents=None,**options):
        if self.resutl_file is None :
            if parents is None and self.resutl_file is None:
                listeParent = [self.env.name+"_cytoplasm"]
                for o in self.env.organelles :
                    listeParent.append(o.name+"_Matrix")
                    listeParent.append(o.name+"_surface")
            elif parents is not None and self.resutl_file is None:       
                listeParent = parents
            listeCenters =self.getPositionsFromObject(listeParent)
        else :
            #use data from file
            listeCenters =self.getPositionsFromResFile(listeParent) 
        #is the distance in the result array ?
        listeDistance=numpy.zeros(len(listeCenters))+99999
        for i in range(len(listeCenters) ):
            for j in range(i+1,len(listeCenters)):
                #should use point
                d = self.helper.measure_distance(listeCenters[i],listeCenters[j])
                if d < listeDistance[i] :
                    listeDistance[i] = d
        return listeDistance
                            
    #should take any type of list...
    def save_csv(self,data,filename=None):
        if filename is None :
            filename = "output.csv"
        resultFile = open(filename,'wb')
        wr = csv.writer(resultFile, dialect='excel')
        #wr.writerows(data) list of list ?
        #resultFile.close()
        for item in data: wr.writerow([item,]) 
        resultFile.close()

#from bhtree import bhtreelib
#bht = bhtreelib.BHtree( verts, None, 10)
#closest = bht.closestPointsArray(tuple(grdPos), diag, returnNullIfFail)
        