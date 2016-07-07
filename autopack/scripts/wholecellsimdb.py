# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 12:09:35 2016

@author: ludov
"""
# lot of information are available on their website
# http://wholecelldb.stanford.edu
# http://www.wholecellkb.org
import h5py
import json
import numpy as np

def getRawData(frame=0, ignore_null=False):
    #name,count,PDBid,id in alldata
    data={}
    if os.path.isfile("C:\\Users\\ludov\\Downloads\\rawrecipe_frame_%i.json" % frame):    
        f = h5py.File("C:\\Users\\ludov\\Downloads\\simulation-1.h5")
        dataMono = f.get("states/ProteinMonomer/counts/data")
        labelMono = f.get("states/ProteinMonomer/counts/labels/0")
    #    labelCompartments = f.get("states/ProteinMonomer/counts/labels/1")
        dataComplex = f.get("states/ProteinComplex/counts/data")
        labelComplex = f.get("states/ProteinComplex/counts/labels/0")
        for i in range(6):
            if ignore_null:
                names = labelMono[dataMono[frame,i,:] > 0]
                counts = dataMono[frame,i,:][dataMono[frame,i,:] > 0]
            else :
                names = labelMono
                counts = dataMono[frame,i,:]#[dataMono[frame,i,:] > 0]
            mono = np.column_stack([names,counts])           
            if ignore_null:
                names = labelComplex[dataComplex[:,i,frame] > 0]
                counts = dataComplex[:,i,frame][dataComplex[:,i,frame] > 0]
            else:
                names = labelComplex
                counts = dataComplex[:,i,frame]
            compl = np.column_stack([names,counts])

            data[i] = np.vstack([mono,compl]).tolist()
        with open("C:\\Users\\ludov\\Downloads\\rawrecipe_frame_%i.json" % frame, 'w') as fp:  # doesnt work with symbol link ?
            json.dump(data,fp)
    else :
        with open("C:\\Users\\ludov\\Downloads\\rawrecipe_frame_%i.json" % frame, 'r') as fp:  # doesnt work with symbol link ?
            data = json.load(fp)        
    return data

data=getRawData(frame=0,ignore_null=True)#ignore_null means if the count is 0 we ignore the entry
#print (data[0])
#there is 6 compartments in this data
#c Cytosol 0
#d DNA     1
#e Extracellular Space  2
#m Membrane  3
#tc Terminal Organelle Cytosol  4
#tm Terminal Organelle Membrane  5
localisation={"c":0,"d":1,"e":2,"m":3,"tc":4,"tm":5}    
#this should be cross with the wholecellKB data that can be download as json file
#you could filter per localisation, pathway etc.
with open("C:\\Users\\ludov\\Downloads\\wholecellmyco.json", 'r') as fp:  # doesnt work with symbol link ?
    celldic = json.load(fp)

mapping={"Chromosomes":[0,1],
        "ChromosomesFeatures":[1,2305], #2305
        "Compartment":[2306,6],      #6
        "Gene":[2312,525],             #525
        "Metabolite":[2837,722],
        "Note":[3559,35],
        "Parameter":[3594,154],
        "Pathway":[3748,16],
        "Process":[3764,29],
        "ProteinComplex":[3793,201],
        "ProteinMonomer":[3994,482],
        "Reaction":[4476,1857],
        "Reference":[6333,930],
        "Species":[7261,1],
        "State":[7262,16],
        "Stimulus":[7278,10],
        "TranscriptionUnit":[7288,335],
        "TranscriptionalRegulation":[7623,30],
        "Type":[7653,140]
        }
def getGenePosIn(GenId,Data):
    for i in range(len(Data)):
        if Data[i]["wid"] == GenId :
            return i
    return -1
    
verbose = 0
gen_data = celldic["data"][mapping["Gene"][0]:mapping["Gene"][0]+mapping["Gene"][1]]
N=mapping["ProteinMonomer"][1]
for i in range(N):
    elem = celldic["data"][mapping["ProteinMonomer"][0]+i]
    gen = elem["gene"]
#    print elem['wid'] 
    uni_id = ''
    gid = getGenePosIn(gen, gen_data)
    if len(gen_data[gid]["cross_references"]) >=4 :
        uni_id = gen_data[gid]["cross_references"][3]['xid']
    if verbose :
        print elem['wid'], uni_id, elem["name"],elem["localization"]
N=mapping["ProteinComplex"][1]
for i in range(N):
    elem = celldic["data"][mapping["ProteinComplex"][0]+i]
    
