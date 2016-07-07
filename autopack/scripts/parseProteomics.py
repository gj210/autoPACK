# -*- coding: utf-8 -*-
"""
Created on Wed Apr 06 10:44:55 2016

@author: ludov
"""

#class protein:
#    def __init__(self,header,mweight,nbmol):
#        self.header = header
#        self.mweight = mweight
#        self.nbmol = nbmol
#
#class recipe:
#    list_prot=[]
#    def addProt(self,header,mweight,nbmol):
#        self.list_prot.append(protein(header,mweight,nbmol))
        
import csv
import math
import numpy

#f="D:\\Data\\cellPAC_data\\cellPACK_database_1.1.0\\proteomics\\M.mycoidesSV.csv"
f="D:\\Data\\cellPAC_data\\cellPACK_database_1.1.0\\proteomics\\M.mycoidesSV_identified.csv"
all_data=[]
cell_radius = 0.15 #um
cell_density = 1.07#g/cc
protein_content_fraction=0.163#155.0/950.0#0.163#by weight
cell_volume = 4.0*math.pi*(math.pow(cell_radius,3.0)/3.0)#cu um
cell_mass = cell_volume*cell_density*math.pow(10,-12)#g
protein_mass = cell_mass*protein_content_fraction#g
total_I=[0,0,0,0]
total_IBAQ=[0,0,0,0]
total_LFQ=[0,0,0,0]
col_id=[42,47,51]
total = [total_I,total_IBAQ,total_LFQ]
oneMol=[total_I,total_IBAQ,total_LFQ]
with open(f, 'r') as csvfile:
    spamreader = csv.reader(csvfile)
    for row in spamreader:
        all_data.append(row)
        if len(all_data) == 1 : continue
        for i in range(3):
            for j in range(4):
                intensity = row[col_id[i]+j]
                if intensity == "NaN" :
                    intensity = 0
                else :
                    intensity = intensity.replace(",",".")
                total[i][j]+=float(intensity)
nbMol=[]
avgMols=[]
mWeight=[]
names=[]
for d in range(1,len(all_data)):
    #compute nbMol
    w = all_data[d][25].replace(",",".")
    mWeight.append(float(w))
    oneMol=[[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    avgMol=[0,0,0]
    for i in range(3):
        for j in range(4):
            intensity = all_data[d][col_id[i]+j]
            if intensity == "NaN" :
                intensity = 0
            else :
                intensity = intensity.replace(",",".")
            oneMol[i][j]=int(round((float(intensity)/total[i][j])*(protein_mass/(mWeight[-1]*1000.0))*6.022e23))
        nb=numpy.array(oneMol[i])
        avgMol[i]=numpy.average(nb[numpy.nonzero(nb)]).astype('int')
    nbMol.append(oneMol)  
    avgMols.append(avgMol)
    names.append(all_data[d][5])
nbmol = numpy.array(avgMols)
#mWeight 25 kDa
#intensity 42,43,44,45
#IBAQ intensity  47,48,49,50
#LFQ intensity 51,52,53,54
#fasta header 5
#ribosome 30-40/100nm3 #100000A3
#2265 diameter is 30*10 (300nm)
#The average of experimentally determined partial specific volumes for soluble, 
#globular proteins is approximately 0.73 cm3/g
#(1.21 x MW) A3/molecule
#nbBeads can be prot_volume/beads_volume
#ABC transporter
#preprotein translocase
#"ATP synthase"
#"phenylalanyl-tRNA synthetase"
def oneProtein(name):
    nb=[0,0,0]
    na=0
    mw=0
    for n in range(len(names)) :
        if names[n].find(name) != -1 :
            print (n,names[n],nbmol[n],mWeight[n])
            for i in range(3):
                nb[i]+=nbmol[n][i]
            na+=1
            mw+=mWeight[n]
    return nb,na,mw/na
