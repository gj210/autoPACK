# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 11:19:57 2015

@author: ludo
"""
import urllib2
from bioservices import WSDbfetch
from bioservices import UniProt
from bioservices import PDB
import json
#http://pax-db.org/api/search?q=ygiF&species=511145

columns_name=["id","entry name","protein names","genes","organism","length","mass"
            "comment(SUBCELLULAR LOCATION)","go(cellular component)","go(molecular function)",
            "go(biological process)",
            "3d","database(PDB)","database(ModBase)","database(ProteinModelPortal)",
            "database(Bgee)","keywords",
            "database(DIP)","existence","go(molecular function)","feature(INTRAMEMBRANE)",
            "feature(TRANSMEMBRANE)","database(PaxDb)"]
            
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

def getTaxonomyProtein(taxonomy):
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
    alldata = u.search(query,frmt='xls',columns=columns)
    dataline = alldata.split("\n")
    data = [l.split("\t") for l in dataline[1:]]
    header = dataline[0].split("\t")
    return data,header

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

data,header = getTaxonomyProtein("mycoplasma mycoides")
res = gatherInfo(data,header)
#u = UniProt()#verbose=False)
#proteomeID="UP000000318"
#taxonomy="83333 - Escherichia coli (strain K12)"
#query="taxonomy:83333"
#frmt="tab"
#columns="id,entry name,reviewed,protein names,genes,organism,length,genes(ALTERNATIVE),comment(MASS SPECTROMETRY),comment(SUBCELLULAR LOCATION),3d,database(Bgee),go(cellular component),keywords,database(PDB),database(ModBase),database(ProteinModelPortal),database(DIP),existence,go(molecular function),feature(INTRAMEMBRANE),feature(TRANSMEMBRANE),database(PaxDb)"
##res=u.search(query,frmt=frmt,columns=columns)
#
#entry_name=u.search(query,frmt=frmt,columns="entry name")
#entry_name_1 = str(entry_name).split("\n")
#enrty_name = entry_name_1[1:-1]
#
##alldata = u.search(query,frmt=frmt,columns="id,entry name,protein names, genes,subcellular locations,go,3d,database(PDB),database(ModBase),database(ProteinModelPortal)existence,go(molecular function),feature(INTRAMEMBRANE),feature(TRANSMEMBRANE),database(PaxDb)")
#alldata = u.search(query,frmt='xls',columns="id,entry name,protein names, genes,subcellular locations,3d,go(cellular component),go(molecular function),go(biological process),database(PDB),database(ModBase),database(ProteinModelPortal),feature,database(PaxDb)")
#dataline = alldata.split("\n")
#
#header = dataline[0].split("\t")
#data = [l.split("\t") for l in dataline[1:]]
#pdbc=[]
##abd=[]
#for i in range(len(data)) :
#    ipdb=header.index('Cross-reference (PDB)')
#    if data[i][ipdb] !='':
#        pdbc.append([i,data[i][ipdb]])
##    elif data[i][12] !='':
##        a=getAbundancy(data[i][12][:-1])
##        abd.append([i,a])
#    else :
#        continue