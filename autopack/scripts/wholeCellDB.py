# -*- coding: utf-8 -*-
"""
Created on Tue May 10 15:50:52 2016

@author: ludov
"""
import os
import json
import numpy as np
import h5py
import urllib2
import time
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Blast.NCBIWWW import qblast
from Bio.Blast import NCBIXML 
from Bio import Entrez
Entrez.email = "autin@scripps.edu"
entrez_db = Entrez.read(Entrez.einfo())

taxonid = 243273
handle = Entrez.efetch("taxonomy", id="243273", retmode="xml")
taxonomy = Entrez.read(handle)
handle = Entrez.esearch(db="protein", term="taxonomy:'Mycoplasma genitalium (strain ATCC 33530 / G-37 / NCTC 10195) [243273]'")

#http://www.uniprot.org/uniprot/?query=taxonomy%3A%22Mycoplasma%20genitalium%20(strain%20ATCC%2033530%20%2F%20G-37%20%2F%20NCTC%2010195)%20%5B243273%5D%22&fil=reviewed%3Ayes&columns=id%2Centry%20name%2Creviewed%2Cprotein%20names%2Cgenes%2Corganism%2Ccomment(DEVELOPMENTAL%20STAGE)%2Ccomment(SUBCELLULAR%20LOCATION)%2C3d%2Cdatabase(ProteinModelPortal)%2Cdatabase(PDB)%2Cdatabase(DisProt)%2Cdatabase(SMR)%2Cdatabase(ExpressionAtlas)%2Cdatabase(Bgee)%2Cdatabase(CleanEx)%2Cdatabase(CollecTF)%2Cdatabase(Genevisible)%2Cdatabase(PaxDb)%2Cdatabase(PeptideAtlas)%2Cdatabase(ProMEX)%2Cdatabase(TopDownProteomics)%2Cdatabase(PRIDE)%2Cdatabase(MaxQB)%2Cdatabase(EPD)%2Ccomment(MASS%20SPECTROMETRY)%2Cmass

#organism:"Mycoplasma genitalium (strain ATCC 33530 / G-37 / NCTC 10195) [243273]"
#http://www.uniprot.org/proteomes/UP000000807
# proteomecomponent:chromosome AND organism:"Mycoplasma genitalium (strain ATCC 33530 / G-37 / NCTC 10195) [243273]" AND proteome:up000000807
# proteomecomponent:chromosome AND organism:"Mycoplasma genitalium (strain ATCC 33530 / G-37 / NCTC 10195) [243273]" AND proteome:up000000807
from bioservices import WSDbfetch
from bioservices import UniProt
from bioservices import PDB
import json
#http://pax-db.org/api/search?q=ygiF&species=511145

columns_name=["id",
              "entry name",
              "protein names",
              "genes",
              "organism",
              "length",
              "sequence",
              "mass",
            "comment(SUBCELLULAR LOCATION)",
            "go(cellular component)",
            "go(molecular function)",
            "go(biological process)",
            "3d",
            "database(PDB)",
            "database(ModBase)",
            "database(ProteinModelPortal)",
            "database(Bgee)",
            "keywords",
            "database(DIP)",
            "existence",
            "go(molecular function)",
            "feature(INTRAMEMBRANE)",
            "feature(TRANSMEMBRANE)",
            "database(PaxDb)"]
            
def getPDB(pdbid):
    s = PDB()
    res = s.get_file("1FBV", "pdb")
    return res
#    import tempfile
#    fh = tempfile.NamedTemporaryFile()
#    fh.write(res)

def getAbundancy(uniprot_id):
    headers = {"Accept" : "application/vnd.paxdb.search+json;version=4"}
    req = urllib2.Request("http://pax-db.org/api/search?q="+uniprot_id, None, headers)
    #req = urllib2.Request("http://pax-db.org/api/search?q=cdc", None, headers)
    response = urllib2.urlopen(req)
    data = response.read()
    response.close()
    result = json.loads(data)
    if int(result["searchResponse"]["hits"]["@totalHits"]) == 0 :
        print "didnt found "+uniprot_id
        return 0.0
    proteins = result["searchResponse"]["hits"]["protein"]
    if int(result["searchResponse"]["hits"]["@totalHits"]) == 1 :
        pid = proteins["@id"]#, proteins[1]["@name"]    
    else :
        pid = proteins[1]["@id"]#, proteins[1]["@name"]
    #then gather the data for the id http://pax-db.org/api/proteins?ids=6883672
    headers = {"Accept" : "application/vnd.paxdb.proteins+json;version=4"}
    req = urllib2.Request("http://pax-db.org/api/proteins?ids="+pid, None, headers)
    #req = urllib2.Request("http://pax-db.org/api/search?q=cdc", None, headers)
    response = urllib2.urlopen(req)
    data = response.read()
    response.close()
    result = json.loads(data)
    #print result["searchResponse"]["hits"]["@totalHits"]
    #proteins = result["searchResponse"]["hits"]["protein"]
    #print proteins[1]["@id"], proteins[1]["@name"]
    #the avg ppm is 
    return float(result["proteinsResponse"]["proteins"]['abundances'][-1]['formattedAbundance'].split()[0])
  
#mycoplasma
#http://www.uniprot.org/uniprot/?query=taxonomy:mycoplasma mycoides&sort=score&columns=id,entry name,reviewed,protein names,genes,organism,length,genes(ALTERNATIVE),comment(MASS SPECTROMETRY),comment(SUBCELLULAR LOCATION),3d,database(Bgee),go(cellular component),keywords,database(PDB),database(ModBase),database(ProteinModelPortal),database(DIP),existence,go(molecular function),feature(INTRAMEMBRANE),feature(TRANSMEMBRANE),database(PaxDb),database(MaxQB),database(PRIDE),database(PeptideAtlas),database(ProMEX),go-id,go(biological process),mass,comment(CAUTION),comment(TISSUE SPECIFICITY),comment(INDUCTION),comment(DEVELOPMENTAL STAGE),database(ExpressionAtlas)
#ecoli
#http://www.uniprot.org/uniprot/?query=taxonomy:83333&columns=id,entry name,reviewed,protein names,genes,organism,length,genes(ALTERNATIVE),comment(MASS SPECTROMETRY),comment(SUBCELLULAR LOCATION),3d,database(Bgee),go(cellular component),keywords,database(PDB),database(ModBase),database(ProteinModelPortal),database(DIP),existence,go(molecular function),feature(INTRAMEMBRANE),feature(TRANSMEMBRANE),database(PaxDb),database(MaxQB),database(PRIDE),database(PeptideAtlas),database(ProMEX),go-id

#u.mapping("ID", "PDB_ID", "P43403")
#u.get_fasta_sequence("P43403")
#Accession via entry name (e.g., ZAP70_HUMAN) is faster than by Entry (e.g., P43403)
def getUniprotInfo(uni_id):
    u = UniProt()#verbose=False)
    frmt="tab"
    columns = ','.join(columns_name)
    alldata = u.search(uni_id,frmt=frmt,columns=columns)
    dataline = alldata.split("\n")
    data = [l.split("\t") for l in dataline[1:]]
    header = dataline[0].split("\t")
    dic_data = []
    for j in  range(len(data)-1):
        dic={}
        for i,key in enumerate(columns_name):
            dic[key] = data[j][i]
        dic_data.append(dic)
    return dic_data,data,header

def getTaxonomyProtein(taxonomy,format="tab"):
    u = UniProt()#verbose=False)
    query="taxonomy:"+taxonomy
    frmt="tab"
    columns = ','.join(columns_name)
    #get all entry_name as a data_frame    
#    entry_name=u.search(query,frmt=frmt,columns="entry name")
#    entry_name_1 = str(entry_name).split("\n")
#    enrty_name = entry_name_1[1:-1]
    #this is no enought informtion 
    #get using the seach
    alldata = u.search(query,frmt=format,columns=columns)
    dataline = alldata.split("\n")
    data = [l.split("\t") for l in dataline[1:]]
    header = dataline[0].split("\t")
    return alldata,data,header

def gatherInfo(data,header):
    #if protein is not "Inferred from homology" or "Predicted"
    result=[]
    r={}
    for i in range(len(data)) :
        if len(data[i])==1 : continue
        if data[i][header.index('Protein existence')] != "Predicted" and data[i][header.index('Protein existence')] != "Uncertain":#"Evidence at protein level":
            r={"id":i}
            ipdb=header.index('Cross-reference (PDB)')
            if data[i][ipdb] !='':
                r["pdb"] = data[i][ipdb]
#            elif data[i][12] !='':
#                a=getAbundancy(data[i][12][:-1])
#                r["abundancy"] = data[i][ipdb]
            else :
                continue   
            result.append(r)
        else :
            continue
    return result

#data,header = getTaxonomyProteins("mycoplasma mycoides")
#res = gatherInfo(data,header)
#http://www.uniprot.org/uniprot/?query=proteomecomponent%3achromosome&fil=organism%3a%22Mycoplasma+genitalium+(strain+ATCC+33530+%2f+G-37+%2f+NCTC+10195)+%5b243273%5d%22+AND+proteome%3aup000000807&offset=25
taxonomy = "C:\\Users\\ludov\\Downloads\\uniprot-proteomecomponent_chromosome.xml"
#p.print_p_contents(text, lookforData=None, lookforTag="sequence",lookforAttr=["mass",] )

#Mycoplasma genitalium G37
#data = urllib.urlopen("http://www.uniprot.org/uniprot/" + code + ".txt").read()
#MG_131
#filename = C:\\Users\\ludov\\Downloads\\wholecellmyco.json
columns_name=["id","entry name","protein names","genes","organism","length","mass"
            "comment(SUBCELLULAR LOCATION)","go(cellular component)","go(molecular function)",
            "go(biological process)",
            "3d","database(PDB)","database(ModBase)","database(ProteinModelPortal)",
            "database(Bgee)","keywords",
            "database(DIP)","existence","go(molecular function)","feature(INTRAMEMBRANE)",
            "feature(TRANSMEMBRANE)","database(PaxDb)"]
            
#from bioservices import UniProt
#u = UniProt(verbose=False)
#data = u.search("zap70+and+taxonomy:9606", frmt="tab", limit=3,
#                 columns="entry name,length,id, genes")

def blastIt(sequence):
    res = qblast("blastp", "pdb", sequence, auto_format=None, 
           composition_based_statistics=None, db_genetic_code=None, 
           endpoints=None, entrez_query='(none)', expect=10.0, filter=None, 
           gapcosts=None, genetic_code=None, hitlist_size=50, i_thresh=None, 
           layout=None, lcase_mask=None, matrix_name=None, nucl_penalty=None, 
           nucl_reward=None, other_advanced=None, perc_ident=None, 
           phi_pattern=None, query_file=None, query_believe_defline=None, 
           query_from=None, query_to=None, searchsp_eff=None, service=None, 
           threshold=None, ungapped_alignment=None, word_size=None, 
           alignments=500, alignment_view=None, descriptions=500, 
           entrez_links_new_window=None, expect_low=None, expect_high=None, 
           format_entrez_query=None, format_object=None, format_type='XML', 
           ncbi_gi=True, results_file=None, show_overview=None, megablast=None)
    data = NCBIXML.parse(res)  
    for record in data:
        if record.alignments :  #skip queries with no matches      
            print "QUERY: %s" % record.query[:60]      
            for align in record.alignments:         
                for hsp in align.hsps:     
#                    if hsp.expect < E_VALUE_THRESH:        
                        print "MATCH: %s " % align.title[:60]        
                        print hsp.expect  
    return res          
           
def getMWFromSequence(sequence):
    X = ProteinAnalysis(sequence)
    mw = X.molecular_weight()  
    return mw
    
def getWholeCellKB(filename):
    if not os.path.isfile(filename):#"C:\\Users\\ludov\\Downloads\\wholecellkb_genitalium.json"):
        request = urllib2.Request("http://www.wholecellkb.org/export/Mgenitalium/?all_model_types=True&format=json")
        response = urllib2.urlopen(request)
        data = response.read()
        time.sleep(1)
        f=open(filename,"w")
        f.write(data)
        f.close()
        return json.loads(data)
    else :
        f=open(filename,"w")
        data = json.load(f)
        f.close()
        return data

def getGenePosIn(GenId,Data):
    for i in range(len(Data)):
        if Data[i]["wid"] == GenId :
            return i
    return -1
    
def getUniprotIdFromGenId(GenId,Data):
    for entry in Data:
        if entry["wid"] == GenId :
            if len(entry["cross_references"]) >= 4:
                return entry["cross_references"][3]['xid']
            else :
                return ""
    return ""
    
def getUniprotData(EntryID):
    url = "http://www.uniprot.org/uniprot/" + EntryID + ".txt"
    #try except :
    response = urllib2.urlopen(url)
    data = response.read()
    return data
    
def parseUniprotData(prot_data,sequence = False):
    lines = prot_data.split("\n")
    PDBid=""
    Sequence=""
    Model=""
    mw = 0
    for i in range(len(lines)) :
        elem = lines[i].split()
        if not len(elem):
            continue
        if elem[0] == "DR":
            if elem[1] == "PDB;":
                PDBid = elem[2][:-1]
            if elem[1] == "ProteinModelPortal;":
                Model = elem[2][:-1]
        if elem[0] == "SQ":
            #get until //
            i=i+1
            while not lines[i].startswith("//"):
                Sequence += lines[i].replace(" ","")
                i=i+1
            mw = getMWFromSequence(Sequence)
    if PDBid == "" and Model != "":
        #grab the model ?
        # print "found ProteinModelPortal id ", Model
        PDBid = Model
    if sequence :
        return PDBid,mw,Sequence
    return PDBid, mw
                
def getDataWID(EntryWID):
    url = 'http://www.wholecellkb.org/detail/Mgenitalium/'+EntryWID+'?format=json'
    # print url
    request = urllib2.Request(url)
    response = urllib2.urlopen(request)
    data = response.read()
#    time.sleep(1)
    return data

def getStructure(uni_id,gen):
#    uni_id = getUniprotIdFromGenId(gene,gen_data)
    m = 0
    if uni_id == '':
        return "null", m
    prot_data = getUniprotData(uni_id)
    p, m = parseUniprotData(prot_data)
    if p == '':
        return "null", m
    return p, m

def sumUpMWComplex(elem,cross):
    #'biosynthesis'
    bio = elem['biosynthesis' ]
    n=len(bio)
    mw=0.0
    for i in range(n-1):
        nb = abs(eval(bio[i]['coefficient']))
        wid = bio[i]['molecule']
        #what if wid not yet in cross ?
        if wid not in cross:
            continue
        mw += nb * cross[wid][-1]
    return mw
    
def composeDataCross(all_data,mapping):
    # get  dictionary {wid:[pdbid,name,molecularweight,gene,access in all_data]}
    # start with monomer
    cross={}
    if not os.path.isfile("C:\\Users\\ludov\\Downloads\\rawcross.json"):
        gen_data = celldic["data"][mapping["Gene"][0]:mapping["Gene"][0]+mapping["Gene"][1]]
        N=mapping["ProteinMonomer"][1]
        for i in range(N):
            elem = celldic["data"][mapping["ProteinMonomer"][0]+i]
            gen = elem["gene"]
            if elem['wid'] not in cross :
                gid = getGenePosIn(gen,gen_data)
                uni_id = ''
                if len(gen_data[gid]["cross_references"]) >=4 :
                    uni_id = gen_data[gid]["cross_references"][3]['xid']
#                pdbId, mw = getStructure(uni_id,gen)
#                if pdbId == "":
#                    pdbId = ""
#                    if len(elem['cross_references']):
#                        for j in range(len(elem['cross_references'])):
#                            pdbId+=elem['cross_references'][j]["xid"]+" " 
                cross[elem['wid']]=[uni_id,mapping["Gene"][0]+gid,mapping["ProteinMonomer"][0]+i,elem["name"],0.0]
        N=access["ProteinComplex"][1]
        for i in range(N):
            elem = celldic["data"][mapping["ProteinComplex"][0]+i]
            pdbId=""
            if len(elem['cross_references']):
                for j in range(len(elem['cross_references'])):
                    pdbId+=elem['cross_references'][j]["xid"]+" "
            if elem['wid'] not in cross :
                mw = sumUpMWComplex(elem, cross)
                cross[elem['wid']]=[pdbId,-1,mapping["ProteinComplex"][0]+i,elem["name"],mw]   
        with open("C:\\Users\\ludov\\Downloads\\rawcross.json", 'w') as fp:  # doesnt work with symbol link ?
            json.dump(cross,fp)
    else :
        with open("C:\\Users\\ludov\\Downloads\\rawcross.json", 'r') as fp:  # doesnt work with symbol link ?
            cross = json.load(fp)       
    return cross

def getRawRecipe(cross,frame=0):
    #name,count,PDBid,id in alldata
    Recipe={}
    if os.path.isfile("C:\\Users\\ludov\\Downloads\\rawrecipe_frame_%i.json" % frame):    
        f = h5py.File("C:\\Users\\ludov\\Downloads\\simulation-1.h5")
        dataMono = f.get("states/ProteinMonomer/counts/data")
        labelMono = f.get("states/ProteinMonomer/counts/labels/0")
    #    labelCompartments = f.get("states/ProteinMonomer/counts/labels/1")
        dataComplex = f.get("states/ProteinComplex/counts/data")
        labelComplex = f.get("states/ProteinComplex/counts/labels/0")
        for i in range(6):
            index = np.nonzero(dataMono[frame,i,:])[0]
            names = labelMono[dataMono[frame,i,:] > 0]
            counts = dataMono[frame,i,:][dataMono[frame,i,:] > 0]
            monoPDBids = [cross[name.split("-")[0]][0] for name in names]
            monoNames = [cross[name.split("-")[0]][-1] for name in names]
    #        monoNames = [widnames[name.split("-")[0]] for name in names]
            mono = np.column_stack([names,index,monoPDBids,counts,monoNames])
            index = np.nonzero(dataComplex[:,i,frame])[0]
            names = labelComplex[dataComplex[:,i,frame] > 0]
            counts = dataComplex[:,i,frame][dataComplex[:,i,frame] > 0]
            complPDBids = [cross[name.split("-")[0]][0] for name in names]
            complNames = [cross[name.split("-")[0]][-1] for name in names]#mw
    #        complNames = [widnames[name.split("-")[0]] for name in names]
            compl = np.column_stack([names,index,complPDBids,counts,complNames])
            Recipe[i] = np.vstack([mono,compl]).tolist()
        with open("C:\\Users\\ludov\\Downloads\\rawrecipe_frame_%i.json" % frame, 'w') as fp:  # doesnt work with symbol link ?
            json.dump(Recipe,fp)
    else :
        with open("C:\\Users\\ludov\\Downloads\\rawrecipe_frame_%i.json" % frame, 'r') as fp:  # doesnt work with symbol link ?
            Recipe = json.load(fp)        
    return Recipe
    
#gather PDB from names and other information
#open json wholecellKB
with open("C:\\Users\\ludov\\Downloads\\wholecell.json", 'r') as fp:  # doesnt work with symbol link ?
    jsondic = json.load(fp)
    
with open("C:\\Users\\ludov\\Downloads\\wholecellmyco.json", 'r') as fp:  # doesnt work with symbol link ?
    celldic = json.load(fp)
    
access={"Chromosomes":[0,1],
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
        
print "uniprot and gene"      
cross = composeDataCross(celldic, access)
print "done with uniprot and gene"  
recipe = getRawRecipe(cross,frame=0)
print "done with raw recipe"  

#mda P47254 ~/Desktop/MDA percent 60 group false color blast 
# P47619 AMD P47254
#1XZP(2) + 3CP2(2) = 3ces
# P47619 MG379 
# PDB gave some Protein Stoichiometry information
#1 chromosomes
#2305 ch features
#http://www.wholecellkb.org/detail/SpeciesWID/EntryWID?format=Format 
#for every entry access the gen ID and get the uniprot code-> PDB or Sequence
#gen_data = celldic["data"][access["Gene"][0]:access["Gene"][0]+access["Gene"][1]]
#gen_data also contain info about metabolite
#
#cross={}
#widnames={}
#no_structure_list=[]
#for entry in jsondic["data"]:
#    # check 'cross_references'
#    if len(entry['cross_references']):
#        # print entry['wid'],entry['cross_references'][0]['xid']
#        cross[entry['wid']] = entry['cross_references'][0]['xid']
#    else :
#        cross[entry['wid']] = "None"
#    widnames[entry['wid']] = entry['name']
#    if entry["model"] == "ProteinMonomer":
#        gene = entry["gene"]
#        uni_id = getUniprotIdFromGenId(gene,gen_data)
#        if uni_id == "" :
#            print "no structure for ",entry['wid'],uni_id,p
#            #predict structure ?
#            no_structure_list.append(entry['wid'])
#            continue
#        # gene_data = json.loads(getDataWID(gene))
#        # uni_id = gene_data['data'][0]["cross_references"][3]['xid']
#        prot_data = getUniprotData(uni_id)
#        p, s = parseUniprotData(prot_data)
#        cross[entry['wid']] = p
#        if p == '':
#            #no structure
#            print "no structure for ",entry['wid'],uni_id,p
#            #predict structure ?
#            no_structure_list.append(entry['wid'])
#    else :
#        if cross[entry['wid']] == "None":
#            print "no structure for ", entry['wid'],entry["model"],cross[entry['wid']]
#            # predic structure ?
#            no_structure_list.append(entry['wid'])
#6 compartments
#percompartments,which is localisation
#pname,i,count
#add Molecular Weight
#or just access it ?
#Recipe={}
#for i in range(6):
#    index = np.nonzero(dataMono[0,i,:])[0]
#    names = labelMono[dataMono[0,i,:] > 0]
#    counts = dataMono[0,i,:][dataMono[0,i,:] > 0]
#    monoPDBids = [cross[name.split("-")[0]] for name in names]
#    monoNames = [widnames[name.split("-")[0]] for name in names]
#    mono = np.column_stack([names,index,monoPDBids,counts,monoNames])
#    index = np.nonzero(dataComplex[:,i,0])[0]
#    names = labelComplex[dataComplex[:,i,0] > 0]
#    counts = dataComplex[:,i,0][dataComplex[:,i,0] > 0]
#    complPDBids = [cross[name.split("-")[0]] for name in names]
#    complNames = [widnames[name.split("-")[0]] for name in names]
#    compl = np.column_stack([names,index,complPDBids,counts,complNames])
#    Recipe[i] = np.vstack([mono,compl]).tolist()
#with open("C:\\Users\\ludov\\Downloads\\rawrecipe.json", 'w') as fp:  # doesnt work with symbol link ?
#    jsondic = json.dump(Recipe,fp)


#make a recipe out of it


#build a cellPACK recipe 
#c Cytosol 0
#d DNA     1
#e Extracellular Space  2
#m Membrane  3
#tc Terminal Organelle Cytosol  4
#tm Terminal Organelle Membrane  5
## only do c,d,e,m
#MG_237 -> 1td6_A template
#MG_200 4DCZ
#MG_289 3eki
#MG_032 2i1kA align with dna binding protein , but say outside in exterior ?
#MG_074 4ZDH 
#