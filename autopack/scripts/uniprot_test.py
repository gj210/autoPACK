# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 11:19:57 2015

@author: ludo
"""
import csv
import os
import math
import urllib2
import bioservices
from bioservices import WSDbfetch
from bioservices import UniProt
from bioservices import PDB
import json
from Bio.Blast.NCBIWWW import qblast
from Bio.Blast import NCBIXML 
import numpy as np
from difflib import SequenceMatcher
import requests
from collections import OrderedDict
#http://pax-db.org/api/search?q=ygiF&species=511145
#http://www.uniprot.org/uniprot/?query=taxonomy%3A%22Mycoplasma%20pneumoniae%20(strain%20ATCC%2029342%20%2F%20M129)%20%5B272634%5D%22&fil=reviewed%3Ayes&columns=id%2Centry%20name%2Creviewed%2Cprotein%20names%2Cgenes%2Corganism%2C3d%2Cdatabase(ProteinModelPortal)%2Cdatabase(PDB)%2Cdatabase(DisProt)%2Cdatabase(SMR)%2Cdatabase(ExpressionAtlas)%2Cdatabase(Bgee)%2Cdatabase(CleanEx)%2Cdatabase(CollecTF)%2Cdatabase(Genevisible)%2Cdatabase(PaxDb)%2Cdatabase(PeptideAtlas)%2Cdatabase(ProMEX)%2Cdatabase(TopDownProteomics)%2Cdatabase(PRIDE)%2Cdatabase(MaxQB)%2Cdatabase(EPD)%2Ccomment(MASS%20SPECTROMETRY)%2Cmass%2Cannotation%20score%2Cfragment%2Ccomment(SUBUNIT%20STRUCTURE)%2Ccomment(SUBCELLULAR%20LOCATION)%2Cfeature(INTRAMEMBRANE)%2Cfeature(TOPOLOGICAL%20DOMAIN)%2Cfeature(TRANSMEMBRANE)%2Cdatabase(PDBsum)%2Cdatabase(IntAct)%2Cdatabase(Reactome)%2Cexistence%2Cgenes(PREFERRED)%2Cgenes(ALTERNATIVE)%2Cgenes(OLN)%2Cgenes(ORF)

#collada
from collada import Collada
from collada import material
from collada import source
from collada import geometry
from collada import scene

import sys
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs\lib-tk")


global pdb_dic
pdb_dic={}
global protein_names_dic
protein_names_dic={}

from HTMLParser import HTMLParser

class MyHTMLParser(HTMLParser):
    def print_p_contents(self, html, lookforTag = "a", lookforData = "[Show]", lookforAttr=None, gatherData=False):
        self.lookforTag  = lookforTag
        self.lookforData = lookforData
        self.lookforAttr = lookforAttr
        self.gatherData = gatherData
        self.tag_stack = False
        self.tag_attr = []
        self.stored = []
        self.stored_attr = []
        self.feed(html)
        
    def handle_starttag(self, tag, attrs):
        if tag == self.lookforTag:
            self.tag_stack = True
            self.tag_attr = []
            # print "Encountered the beginning of a %s tag" % tag
            # print attrs           
            if len(attrs):
                if self.lookforAttr is not None :
                    for atr in self.lookforAttr :
                        for attr in attrs :
                            if attr[0] == atr :
                                self.tag_attr.append(attr[1])
                else :
                    self.tag_attr = attrs[0][1]
            
    def handle_endtag(self, tag):
        self.tag_stack = False
        #self.tag_attr = None
        pass
#        print "Encountered the end of a %s tag" % tag

    def handle_data(self, data):
        if self.tag_stack:
            if self.lookforData is None :
                if self.gatherData:
                    self.stored.append(data)
                else :
                    if len(self.tag_attr):
                        self.stored.append(self.tag_attr)                
            else :
                if data == self.lookforData:
                    #print "Encountered data %s" % data
                    #print self.tag_attr 
                    self.stored.append(self.tag_attr)
                
                
columns_name=["id","entry name","reviewed","protein names","genes","genes(OLN)","organism","3d",
         "database(ProteinModelPortal)","database(PDB)","database(DisProt)","database(SMR)",
         "database(ExpressionAtlas)","database(Bgee)","database(CleanEx)","database(CollecTF)",
         "database(Genevisible)","database(PaxDb)","database(PeptideAtlas)","database(ProMEX)",
         "database(TopDownProteomics)","database(PRIDE)","database(MaxQB)","database(EPD)",
         "comment(MASS SPECTROMETRY)","mass","annotation score","fragment","comment(SUBUNIT STRUCTURE)",
         "comment(SUBCELLULAR LOCATION)","feature(INTRAMEMBRANE)","feature(TOPOLOGICAL DOMAIN)",
         "feature(TRANSMEMBRANE)","database(PDBsum)","database(IntAct)","database(Reactome)","existence"]#,"sequence"]   

for c in columns_name :
    if c not in bioservices.uniprot.UniProt._valid_columns:
        bioservices.uniprot.UniProt._valid_columns.append(c)
  


from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBList import PDBList
data_folder = "D:\\Data\\cellPACK_data\\cellPACK_database_1.1.0\\other"
geom_folder = "D:\\Data\\cellPACK_data\\cellPACK_database_1.1.0\\geometries\\"

fetch = PDBList(pdb=data_folder)
p = PDBParser(PERMISSIVE=1,QUIET=True)

k=["ID","PDB","BU","ORG","SCORE","LOC","STOICH","Name"]
summary = OrderedDict()
if os.path.isfile("D:\\Data\\cellPACK_data\\Mycoplasma\\SummaryProteome.csv"):
    csvfile = open('D:\\Data\\cellPACK_data\\Mycoplasma\\SummaryProteome.csv')
    reader = csv.DictReader(csvfile)  
    for row in reader:
        if row['Name']!='' :
            summary[row['ID'].lower()] = OrderedDict(row)
    csvfile.close()

def colladaMesh(name,v,n,f,collada_xml,matnode=None):
    vertxyz = np.array(v)# * numpy.array([1,1,-1])
    # vertxyz = numpy.vstack([x,y,z]).transpose()* numpy.array([1,1,-1])        
    vert_src = source.FloatSource(name+"_verts-array", vertxyz.flatten(), ('X', 'Y', 'Z'))
    input_list = source.InputList()
    input_list.addInput(0, 'VERTEX', "#"+name+"_verts-array")
    if type(n) != type(None) and len(n)  :
        norzyx=np.array(n)
        nz,ny,nx=norzyx.transpose()
        norxyz = np.vstack([nx,ny,nz]).transpose()#* numpy.array([1,1,-1])
        normal_src = source.FloatSource(name+"_normals-array", norxyz.flatten(), ('X', 'Y', 'Z'))
        geom = geometry.Geometry(collada_xml, "geometry"+name, name, [vert_src,normal_src])
        input_list.addInput(0, 'NORMAL', "#"+name+"_normals-array")
    else :
        geom = geometry.Geometry(collada_xml, "geometry"+name, name, [vert_src])        
    fi=np.array(f,int)        
    triset = geom.createTriangleSet(fi.flatten(), input_list, name+"materialref")
    geom.primitives.append(triset)
    collada_xml.geometries.append(geom)
    master_geomnode = scene.GeometryNode(geom, [matnode])
    master_node = scene.Node("node_"+name, children=[master_geomnode,])#,transforms=[tr,rz,ry,rx,s])
    return master_node        
    
def exportCollada(name, v,f,n,filename):

    collada_xml = Collada()
    collada_xml.assetInfo.unitname="centimeter"
    collada_xml.assetInfo.unitmeter=0.01
    collada_xml.assetInfo.upaxis="Y_UP"

    root_env=scene.Node(name)
    myscene = scene.Scene(name+"_Scene", [root_env])
    collada_xml.scenes.append(myscene)
    collada_xml.scene = myscene                

    effect = material.Effect("effect0", [], "phong", diffuse=(1,0,0), specular=(0,1,0))
    mat = material.Material("material0", "mymaterial", effect)
    collada_xml.effects.append(effect)
    collada_xml.materials.append(mat)

    node = colladaMesh(name,v,n,f,collada_xml,matnode=mat)
    root_env.children.append(node)  
    
    collada_xml.write(filename)
    return collada_xml


def exportOBJ(name, v,face,n,filename):
    f = open(filename,"w")
    for i in range(len(v)):
        f.write("v %f %f %f\n"%(v[i][0],v[i][1],v[i][2]))
    for i in range(len(n)):
        f.write("vn  %f %f %f\n"%(n[i][0],n[i][1],n[i][2]))
    for i in range(len(face)):
        #f v1//vn1 v2//vn2 v3//vn3 ...
        f.write("f %i// %i// %i//\n"%(face[i][0]+1,face[i][1]+1,face[i][2]+1))
    f.close()
    
def coarseMolSurface(coords, radii, XYZd =[16,16,16],isovalue=6.0,resolution=-0.1,padding=0.0,
                         name='CoarseMolSurface',geom=None):
    from UTpackages.UTblur import blur
    print "res",resolution
    #if radii is None :
    radii = np.ascontiguousarray(np.ones(len(coords))*2.6).tolist()
    volarr, origin, span = blur.generateBlurmap(np.ascontiguousarray(coords).tolist(), radii, XYZd,resolution, padding = 0.0)
    volarr.shape = (XYZd[0],XYZd[1],XYZd[2])
#        print volarr
    volarr = np.ascontiguousarray(np.transpose(volarr), 'f')
    #weights =  np.ones(len(radii), typecode = "f")
    h = {}
    from Volume.Grid3D import Grid3DF
    maskGrid = Grid3DF( volarr, origin, span , h)
    h['amin'], h['amax'],h['amean'],h['arms']= maskGrid.stats()
    #(self, grid3D, isovalue=None, calculatesignatures=None, verbosity=None)
    from UTpackages.UTisocontour import isocontour
    isocontour.setVerboseLevel(0)

    data = maskGrid.data

    origin = np.array(maskGrid.origin).astype('f')
    stepsize = np.array(maskGrid.stepSize).astype('f')
    # add 1 dimension for time steps amd 1 for multiple variables
    if data.dtype.char!=np.float32:
#            print 'converting from ', data.dtype.char
        data = data.astype('f')#Numeric.Float32)

    newgrid3D = np.ascontiguousarray(np.reshape( np.transpose(data),
                                          (1, 1)+tuple(data.shape) ), data.dtype.char)
#        print "ok"       
    ndata = isocontour.newDatasetRegFloat3D(newgrid3D, origin, stepsize)

#        print "pfff"
    isoc = isocontour.getContour3d(ndata, 0, 0, isovalue,
                                       isocontour.NO_COLOR_VARIABLE)
    vert = np.zeros((isoc.nvert,3)).astype('f')
    norm = np.zeros((isoc.nvert,3)).astype('f')
    col = np.zeros((isoc.nvert)).astype('f')
    tri = np.zeros((isoc.ntri,3)).astype('i')
    isocontour.getContour3dData(isoc, vert, norm, col, tri, 0)
    #print vert
    if maskGrid.crystal:
        vert = maskGrid.crystal.toCartesian(vert)
    return vert, norm, tri
    
def computeOPM(filename,pdbId):
    if not os.path.isfile(data_folder+os.sep+pdbId+"_mb.pdb"):  
        url = "http://sunshine.phar.umich.edu/upload_file.php"
        payload = {"submit":"Submit","inout":"in","yesno":"no"}#,"userfile":filename}
        r = requests.post(url, data=payload,files={"userfile":open(filename,"r")},timeout=60.0)
        p = MyHTMLParser()
        p.print_p_contents(r.content,lookforData=None, lookforTag="a",lookforAttr=['href'])
        if not len(p.stored) :
            return ""
        url_download = "http://sunshine.phar.umich.edu/"+p.stored[0][0]
        response = urllib2.urlopen(url_download)
        res = response.read()
        f = open(data_folder+os.sep+pdbId+"_mb.pdb","w")
        f.write(res)
        f.close()
    return data_folder+os.sep+pdbId+"_mb.pdb"
    
def getOPM(pdbId):
    if not os.path.isfile(data_folder+os.sep+pdbId+"_mb.pdb"):    
        search_url = "http://opm.phar.umich.edu/protein.php?search="+pdbId#1l7v
        response = urllib2.urlopen(search_url)
        res = response.read()
        #<a href="pdb/1l7v.pdb">Download Coordinates</a>
        p = MyHTMLParser()
        p.print_p_contents(res,lookforData="Download Coordinates", lookforTag="a",lookforAttr=['href'])
        if not len(p.stored) :
            #check for reference
            pattern = "protein.php?pdbid="
            ind = res.find(pattern)
            if ind != -1 :
                #do it again on new ref pdb
                newpdb = res[ind+len(pattern):ind+len(pattern)+4]
                return getOPM(newpdb)
            else :
                print "problem ",pdbId
                return "",pdbId
        url = "http://opm.phar.umich.edu/"+p.stored[0][0]
        print url
        #direct access
        try :
            response = urllib2.urlopen(url)
            res = response.read()
            f = open(data_folder+os.sep+pdbId+"_mb.pdb","w")
            f.write(res)
            f.close()
            return data_folder+os.sep+pdbId+"_mb.pdb",pdbId
        except:
            #not found
            #build
            print ("problem accessing " + url)
            return "",pdbId
    else :
        return data_folder+os.sep+pdbId+"_mb.pdb",pdbId


def gunzip(gzpath, path):
    import gzip
    gzf = gzip.open(gzpath)
    f = open(path, 'wb')
    f.write(gzf.read())
    f.close()
    gzf.close()

#get bioassamble n
def fetch_biounit(pdbId,bio):
    suffix='.pdb%d'%bio
    site='ftp.wwpdb.org'
    url_pattern='ftp://%s/pub/pdb/data/biounit/coordinates/all/%s.gz'
    name=pdbId+suffix
    file_url=url_pattern%(site,name.lower())
    response=urllib2.urlopen(file_url)
    res=response.read()
    f=open(data_folder+os.sep+name+".gz","wb")
    f.write(res)
    f.close()
    gunzip(data_folder+os.sep+name+".gz",data_folder+os.sep+pdbId+"_"+str(bio)+".pdb")


def getBioAssamblyInfo(pdbid,bio=0):  
    #pdb Id + bu
    url = "http://www.rcsb.org/pdb/json/symmetryOrientation?pdbID="+pdbid+"&bioassembly="+str(bio)
    response = urllib2.urlopen(url)
    res = response.read()
    allinfo = json.loads(res)    
    return allinfo
    #resAU['symmetries'][0]['stoichiometry']
    
#http://pir.georgetown.edu/cgi-bin/pairwise.pl?seq_data1=&option=Submit &seq_id1=P53039&seq_data2=&seq_id2=Q6FQ69
#align2uniEntry("P53039","Q6FQ69")
def align2uniEntry(entry1,entry2):
    url="http://pir.georgetown.edu/cgi-bin/pairwise.pl"
    payload = {"submit":"Submit","inout":"in","yesno":"no","seq_id1":entry1,"seq_id2":entry2}#,"userfile":filename}
    r = requests.post(url, data=payload, timeout=60.0)
    indice = r.content.find("Smith-Waterman score")
    Score = r.content[indice:indice+100].split("\n")[0]
    return Score
#    <PRE>
#    >>Q6FQ69Q6FQ69_CANGA  ProteinYIP (247 aa)
#     s-w opt: 1194  Z-score: 1466.6  bits: 278.9 E(): 7e-80
#    Smith-Waterman score: 1194; 70.8% identity (91.2% similar) in 250 aa overlap (1-248:1-247)    
    
def getSource(pdbId):
    url = "http://www.rcsb.org/pdb/rest/customReport.xml?pdbids=%s&customReportColumns=source,uniprotAcc,experimentalTechnique"%pdbId
    response = urllib2.urlopen(url)
    res = response.read()
    #import xml.etree.ElementTree
    #e = xml.etree.ElementTree.parse('thefile.xml').getroot()
    from xml.dom import minidom
    xmldoc = minidom.parseString(res)
    objects = xmldoc.getElementsByTagName('dimStructure.experimentalTechnique')        
    types = ""
    if len(objects) != 0 :
        types = objects[0].firstChild.nodeValue
    if types.lower().find("model") != -1 :
        source = "model"
    else :
        objects = xmldoc.getElementsByTagName('dimEntity.source')        
        source = ""
        if len(objects) != 0 :
            source = objects[0].firstChild.nodeValue
#    objects = xmldoc.getElementsByTagName('dimEntity.uniprotAcc')  
#    #need the uniprot of the proper chain    
#    uniProt = ""
#    if len(objects) != 0 :
#        uniProt =objects[0].firstChild.nodeValue
    return source#,uniProt
    
def getPDBDescription(pdbid):
    #http://www.rcsb.org/pdb/rest/describeMol?structureId=1aon
    #http://www.rcsb.org/pdb/rest/describeMol?structureId=1we3
    #need a xml parser
    #print "grab " +pdbid
    if pdbid not in pdb_dic:
        url = "http://www.rcsb.org/pdb/rest/describeMol?structureId="+pdbid
        response = urllib2.urlopen(url)
        res = response.read()
        #import xml.etree.ElementTree
        #e = xml.etree.ElementTree.parse('thefile.xml').getroot()
        from xml.dom import minidom
        xmldoc = minidom.parseString(res)
        polymers_objects = xmldoc.getElementsByTagName('polymer')
        result = []
        for p in polymers_objects:
            polymer=""
            name = ""
            accession=""
            chains = p.getElementsByTagName("chain")
            descriptions = p.getElementsByTagName('polymerDescription')
            if len(descriptions) : 
                polymer = descriptions[0].getAttribute("description")
            synonyms = p.getElementsByTagName('synonym')
            if len(synonyms) : 
                name = synonyms[0].getAttribute("name")
            accessions = p.getElementsByTagName('accession')
            if len(accessions) : 
                accession = accessions[0].getAttribute("id")
            result.append([polymer,name,len(chains),accession])
        #count number of chains per polymer ?
        pdb_dic[pdbid] = result
        #get the bioassambly nb
        url = "http://www.rcsb.org/pdb/rest/getEntityInfo?structureId="+pdbid
        response = urllib2.urlopen(url)
        res = response.read()
        xmldoc = minidom.parseString(res)
        pdb_elem = xmldoc.getElementsByTagName('PDB')
        ass = pdb_elem[0].getAttribute("bioAssemblies")
        nAssambly = 0
        if ass != "" :
            nAssambly = int(pdb_elem[0].getAttribute("bioAssemblies"))
        info = None
        if nAssambly != 0 :
            info = getBioAssamblyInfo(pdbid)
        pdb_dic[pdbid]=[result,nAssambly,info]#if not 0 it means biomt, also different biomolecule in the file
        result = [result,nAssambly,info]
    else :
        result = pdb_dic[pdbid]
    return result

def compareDescriptionLabel(descr,label):
    descr=descr.lower()
    label=label.lower()
    s=SequenceMatcher(None, descr,label)
    match = s.find_longest_match(0,len(descr),0,len(label))
    return [s.ratio(),match.size]

def pickAmodel(length):
    maxLength = 0
    longest = -1
    for i in range(len(length)):
        ln = eval(length[i][1])-eval(length[i][0])
        if ln > maxLength :
            maxLength = ln
            longest = i
    return longest

def pickAModelOnDescription(templates, label):
    choice = [0,0]
    maxscore = 0
    for i in range(len(templates)):
        #gettemplateDescription
        if len(templates[i]) == 0 :
            continue
        descr = getPDBDescription(templates[i][:4])#[polymer,name,len(chains)]
        #check all description for this given PDB, store in a dictionary of PDB: polymer-nbChain/stochio
        nbPoly = len(descr[0])
        #do we have the proper name?
        for j,d in enumerate(descr[0]):
            score = compareDescriptionLabel(d[0],label)
            if score[0] > 0.5 :
                if score[0]> maxscore : 
                    maxscore = score[0]
                    choice = [i,j]
                #print j,templates[i],d[0],score, nbPoly
#    print "should use ",label,choice
#    print templates[choice[0]]
#    print pdb_dic[templates[choice[0]][:4]]
    return templates[choice[0]]    

def getPolymerAccession(pdbId,label):
    choice = 0
    maxscore = -1
    accession = ""
    descr = getPDBDescription(pdbId)
    nbPoly = len(descr[0])
    #do we have the proper name?
    for j,d in enumerate(descr[0]):  
        score = compareDescriptionLabel(d[0],label)
        if score[0] > 0.5 :
            if score[0]> maxscore : 
                maxscore = score[0]
                choice = j
                accession = d[3]
    if maxscore == -1:
        accession = descr[0][0][3]
    return accession,j,maxscore

def GetScoreAndSource(label,uniId,pdbId):
    accession,j,maxscore = getPolymerAccession(pdbId,label)
    score = align2uniEntry(uniId,accession)
    source = getSource(pdbId)
    return score, source
    
def filterTemplatesComplex(complexe,uniData):
    #loop through template until finding one that has the subunit
    choice = [0,0]
    maxscore = 0
    gscores={}
    name = complexe[0]
    subunitgene = complexe[2]
    #each sunbuti gene can have a template (that may come from something else than protein portal model)
    all_templates=[]
    for akey in subunitgene :    
        if akey not in uniData :
            continue
        templates = uniData[akey]['Template'].split(",")
        all_templates.extend(templates)
        gscores[akey]=[0,[0,0,""]]
    #find one template that has all proteins
    nGene = len(subunitgene)
    if len(all_templates) == 0 :
        return gscores
    for i in range(len(all_templates)):
        if len(all_templates[i]) == 0 :
            continue
        if all_templates[i] == "None" :
            continue
        descr = getPDBDescription(all_templates[i][:4])
        nbPoly = len(descr[0])
        
        #p = float(nGene)/float(nbPoly)
        if nbPoly == 0 :
            print "no desc", all_templates[i][:4],descr
            
#        if nbPoly == 2:
#            descr = [descr,]
        de = nGene - nbPoly
        if de < 2 :
            for j,d in enumerate(descr[0]):
                #score agains all gene
                if type(descr[0][j]) == int :
                    print "pb", descr,all_templates[i][:4]
                    continue
                for akey in subunitgene : 
                    if akey not in uniData :
                         continue
                    label = uniData[akey]["Name"].lower()
                    score = compareDescriptionLabel(descr[0][j][0],label)
                    if score[0] > 0.5 :
                        if score[0]> gscores[akey][0] : 
                            gscores[akey][0] = score[0]
                            gscores[akey][1] = [i,j,all_templates[i][:4]]
    return gscores
                        
def getListPArtner(interactId):
    #should return list of gene_name/protein name of partner
    return []
    
def getPMPTemplates(uid):
    p = MyHTMLParser()
    url = "http://www.proteinmodelportal.org/query/up/"+uid
    response = urllib2.urlopen(url)
    res = response.read()
    #find the Show button
    p.print_p_contents(res,lookforData="[Show]", lookforTag="a",lookforAttr=['href'])
    #first one should be the best
    if len(p.stored) == 0 :
        return "None"
    else :
        pdbs=""
        for i in range(len(p.stored)):
            pdbs+=p.stored[i][0].split("template=")[1].split("&pmpuid")[0]+","
        return pdbs[:-1]
  
def queryPMP(uid,name):
    print "get ",uid,name
    if not os.path.isfile(data_folder+os.sep+name+".pdb"):
        p = MyHTMLParser()
        url = "http://www.proteinmodelportal.org/query/up/"+uid
        response = urllib2.urlopen(url)
        res = response.read()

        #find the Show button
        p.print_p_contents(res,lookforData="[Show]", lookforTag="a",lookforAttr=['href'])
        #p.print_p_contents(res)
        urls = p.stored[:]
        if not len(p.stored) :
            return ""
        p.print_p_contents(res,lookforData=None, lookforTag="span",lookforAttr=['data-pmpseqstart','data-pmpseqend'])
        length = p.stored[:]
        i = pickAmodel(length)
        url = "http://www.proteinmodelportal.org/"+urls[i][0]
        response = urllib2.urlopen(url)
        res = response.read()
        p.print_p_contents(res,lookforTag= "a",lookforData="[ download ]",lookforAttr=['href'])
        url = p.stored[0][0]
        response = urllib2.urlopen(url)
        res = response.read()
        #check if there is tag        
        p.print_p_contents(res,lookforData=None, lookforTag="content",lookforAttr=None,gatherData=True)
        if len(p.stored):
            res = p.stored[0][1:-2]
        f = open(data_folder+os.sep+name+".pdb","w")
        f.write(res)
        f.close()
    return data_folder+os.sep+name+".pdb"
       
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
    #alldata = u.search(query,frmt='xls',columns=columns)
    alldata = u.search(query,frmt='tab',columns=columns)
    dataline = alldata.split("\n")
    data = [l.split("\t") for l in dataline[1:]]
    header = dataline[0].split("\t")
    return data,header

def listToDic(row, header):
    rdic={}
    for i in range(len(header)):
        rdic[header[i]] = row[i]
    return rdic
    
def gatherInfo(data,header):
    #if protein is not "Inferred from homology" or "Predicted"
    result={}
    r={}
    for i in range(len(data)) :
        if len(data[i])==1 : continue
#        if data[i][header.index('Protein existence')] != "Predicted" and data[i][header.index('Protein existence')] != "Uncertain":#"Evidence at protein level":
        r={"id":i}
        r["name"] = str(data[i][header.index('Gene names  (ordered locus )')]).lower().replace("_","")
        r["label"] = data[i][header.index('Protein names')]#Gene names (ordered locus )
        r["udic"] = listToDic(data[i],header)
        r["pdb"] = "None" 
        ipdb=header.index('Cross-reference (PDB)')
        if data[i][ipdb] !='':
            r["pdb"] = data[i][ipdb]
        elif data[i][header.index('Cross-reference (ProteinModelPortal)')]!='':
            print "found model ",data[i][header.index('Cross-reference (ProteinModelPortal)')]#templates for PMP
            pmpid = data[i][header.index('Cross-reference (ProteinModelPortal)')]
            r["pdb"] = getPMPTemplates(pmpid)
        else :
            #get the sequence ?
            r["pdb"] = "None"   
        if data[i][header.index('Cross-reference (PaxDb)')] !='':
            #a=getAbundancy(data[i][12][:-1])
            r["abundancy"] = data[i][ipdb]               
        result[r["name"]]=r
#        else :
#            print "xistence :",data[i][header.index('Protein existence')]
#            continue
    return result

def getRadiusFromMW(mw):
    V = mw * 1.21    
    return math.pow((3*V)/(4*math.pi),1.0/3.0)


def parse_PDB_BIOMT(filename):
    #row major matrix
    #REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B, C, D, E
    #REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
    #REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
    #REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
    dic_mat={}#chain_group, liste_matrix
    chains = []
    for line in open(filename, 'r'):
        if line.startswith("REMARK 350 APPLY THE FOLLOWING TO") :
            chains=line.split("TO CHAINS: ")[1][:-1]#.split(",")
            if chains not in dic_mat:
                dic_mat[chains]={}
        if line.startswith("REMARK 350   BIOMT") :
            spl = line.split()
            symOpNum = int(spl[3])
            symx = float(spl[4])  # float(l[23:33])
            symy = float(spl[5])  # float(l[33:43])
            symz = float(spl[6])  # float(l[43:53])
            tr = float(spl[7])  # float(l[53:])
            if symOpNum not in dic_mat[chains]:
                dic_mat[chains][symOpNum] = []
            dic_mat[chains][symOpNum].append([symx, symy, symz, tr])
    return dic_mat
 
def ApplyMatrix(coords,mat):
    mat = np.array(mat)
    coords = np.array(coords)
    one = np.ones( (coords.shape[0], 1), coords.dtype.char )
    c = np.concatenate( (coords, one), 1 )
    return np.dot(c, np.transpose(mat))[:, :3]


def getBU(pdb,dstch):
    bu = 0
    if len(pdb) == 4 and pdb.lower() !="none":
        descr = getPDBDescription(pdb)
        nBU=descr[1]
        for i in range(nBU+1):
            info = getBioAssamblyInfo(pdb,bio=i)
            if len(info['symmetries']):
                pdb_stch=1
                pdb_stch_str = info['symmetries'][0]['stoichiometry'].split("A")[1].split("B")[0]#problem with A2B2
                #print pdb,i,pdb_stch_str,DData[ingr.o_name]["stoich"],info['symmetries'][0]['stoichiometry']
                if pdb_stch_str!='' and len(pdb_stch_str) != 0:
                    pdb_stch = int(pdb_stch_str)
#                print pdb_stch,pdb_stch_str
                if pdb_stch == int(dstch) :
                    bu = i
                    return bu
    return bu
    
def cluster(coords):
    #center should be the surface offset ?
    #pcpalVector ?
    proxy_radius=11.85
    Vproxy = 4*math.pi*(proxy_radius*proxy_radius*proxy_radius)/3.0
    V = len(coords) * 10.0 * 1.21
    nProxy = int(round(V/Vproxy))
    if nProxy < 10:
        nProxy = 10
    print "ncluster ", nProxy, len(coords)
    from sklearn.cluster import KMeans#MiniBatchKMeans
    k_means = KMeans(init='k-means++', n_clusters=nProxy, n_init=10)
    k_means.fit(coords)
    return k_means.cluster_centers_    
    
def build_proxy_biomt(pdbId, filename):
    lmat = parse_PDB_BIOMT(filename)
    s = p.get_structure(pdbId, filename)
    chains_coords={}
    chains_proxy={}
    for m in s.get_models():
        for ch in  m.get_chains():
            chains_coords[ch.id] = [atom.coord.tolist() for atom in ch.get_atoms() if atom.parent.resname != "DUM"]
            chains_proxy[ch.id] = cluster(chains_coords[ch.id])
    #build proxy for all chain, apply the transformation, save
    all_cluster_pos=[]
    all_atoms=[]
    for ch_key in lmat:
        for ch in ch_key.split(","):
            print ch,len(lmat[ch_key])
            #apply the trasnformation
            for mat in lmat[ch_key]:
                #transform the coordinate
                amat = lmat[ch_key][mat]
                if len(amat) == 3:
                    amat.append([0., 0., 0., 1.])  # ?
                amat = np.array(amat)
                new_coords = ApplyMatrix(chains_proxy[ch.strip()],amat)
                new_atoms =  ApplyMatrix(chains_coords[ch.strip()],amat)
                all_cluster_pos.extend(new_coords) 
                all_atoms.extend(new_atoms.tolist()) 
    np.savetxt(data_folder+os.sep+pdbId+"_cl.txt",all_cluster_pos)
    mesh_file = geom_folder+os.sep+pdbId+"_cms.dae"
    if not os.path.isfile(mesh_file) :
        vert, norm, tri = coarseMolSurface(all_atoms,None)
        dae = exportCollada(pdbId, vert,tri,norm,mesh_file)
        exportOBJ(pdbId, vert,tri,norm,geom_folder+os.sep+pdbId+"_cms.obj")
    return all_cluster_pos
      
def buildProxy(PDBid,proxy_radius,name,mw,stch=0,surface=False, overwrite=False, biomt=False):
    # read the pdb
    # build cluster
    # pass pos/radii
    # check for CA only ?
    # should return two level of sphere, encapsuilating then cluster
    fname = PDBid
    overwrite = False
    name = name.split("-")[0]
    mesh = []
    filename=""
    structure_id = PDBid
    center=[0,0,0]
    usebu=False
    
    #print "build ",PDBid,surface
    R = getRadiusFromMW(mw)
    if PDBid == "none":
        print "PDBid is none"
        return [[[0,0,0]]],[[R]],R,[],"None",center
    if len(PDBid) > 4 :
        #file should be in data folder
        filename = data_folder+os.sep+PDBid+".pdb"   
        if not os.path.isfile(filename):  
            print "?? ",PDBid,name,filename
            return [[[0,0,0]]],[[R]],R,[],"None",center
        #if surface :
        #    filename = computeOPM(filename,PDBid)
        fname = PDBid
    elif len(PDBid) == 4 :
        #check stoichiometry?
        if surface :
            filename,PDBid = getOPM(PDBid)       
            if filename == "":
                temp_filename = data_folder+os.sep+PDBid+".pdb"
                if not os.path.isfile(temp_filename):  
                    temp_filename = fetch.retrieve_pdb_file(PDBid)
                #opm doesnt accept file name as pdbPDBI.ent
                filename = computeOPM(temp_filename,PDBid)
            
        else :
            #test if it exist in the othe folder
            filename = data_folder+os.sep+PDBid+".pdb"
            bu = getBU(PDBid,stch) 
            print "********** use bu ",bu,stch,PDBid
            if bu != 0 :
               filename = data_folder+os.sep+PDBid+"_"+str(bu)+".pdb"
               if not os.path.isfile(filename): 
                   fetch_biounit(PDBid,int(bu)) 
               fname = PDBid+"_"+str(bu)
               usebu = True
            else :
                if not os.path.isfile(filename):  
                    #data_folder
                    filename = fetch.retrieve_pdb_file(structure_id)
        fname = PDBid
    elif len(PDBid) > 0 :
        #parse from protein model portal
        filename = data_folder+os.sep+PDBid+".pdb"
        if not os.path.isfile(filename):  
            print "?? ",PDBid,name,filename
            return [[[0,0,0]]],[[R]],R,[],"None",center
#        if surface :
#            filename = computeOPM(filename,PDBid)
        fname = PDBid
    else:
        print "?? ",PDBid,name
        return [[[0,0,0]]],[[R]],R,[],"None",center
        #filename = queryPMP(PDBid,name)
        #fname = name
        #if surface and filename != "":
        #   filename = computeOPM(filename,name)
        #http://salilab.org/modbase/retrieve/modbase?databaseID=P21812
    if filename == "" :
        #gather the radius at least
        return [[[0,0,0]]],[[R]],R,[],"None",center
    if len(fname) == 4 and surface :
        fname = fname+"_mb"
    print "clustering ",fname
    mesh_file = geom_folder+os.sep+fname+"_cms.dae" # should be pointer to mesh file, dae ?
    if not os.path.isfile(data_folder+os.sep+fname+"_cl.txt") or overwrite :   
        if fname.find("biomt") != -1 :
             centroids = build_proxy_biomt(fname, filename)
             nProxy = len(centroids)
             R = np.linalg.norm(centroids,axis=1).max()+proxy_radius
             center = [0,0,0]
        else :
            s = p.get_structure(structure_id, filename)
#            if usebu :
            atoms_coord=[]
            for m in s.get_models():
                coords = [atom.coord.tolist() for atom in m.get_atoms() if atom.parent.resname != "DUM"]
                atoms_coord.extend(coords)
#            else :
#                for m in s.get_models():
#                    atoms_coord = [atom.coord.tolist() for atom in m.get_atoms() if atom.parent.resname != "DUM"]
#                    break
            center = np.average(atoms_coord,axis=0)
            atoms_coord_centerd = np.array(atoms_coord) - center
            #center should be the surface offset ?
            #pcpalVector ?
            R = np.linalg.norm(atoms_coord_centerd,axis=1).max()
            Vproxy = 4*math.pi*(proxy_radius*proxy_radius*proxy_radius)/3.0
            V = len(atoms_coord) * 10.0 * 1.21
            nProxy = int(round(V/Vproxy))
            if nProxy < 10:
                nProxy = 10
            print "ncluster ", nProxy, len(atoms_coord),structure_id
            # ~(1.21 x MW) A**3/molecule
            # V=4*math.pi*(r*r*r)/3.0
            r =  math.pow((3*V)/(4*math.pi),1.0/3.0)
            # r = math.pow(r3,1.0/3.0)
            #cluster
        #	0.73 cm**3/g  x  10**24 A**3/cm**3  x  molecular weight g/mole
        #	--------------------------------------------------------------
        #			 6.02 x 10**23 molecules/mole
            from sklearn.cluster import KMeans#MiniBatchKMeans
    #        from sklearn.cluster import AffinityPropagation
    #        from sklearn import metrics
    #        from sklearn.datasets.samples_generator import make_blobs
    
            k_means = KMeans(init='k-means++', n_clusters=nProxy, n_init=10)
            k_means.fit(atoms_coord_centerd)
            #k_means_labels = k_means.labels_
            #k_means_cluster_centers = k_means.cluster_centers_
            #k_means_labels_unique = np.unique(k_means_labels)
            #np.savetxt(pdb_directory+os.sep+pdbname+"_kmeans15.txt",
            #       k_means.cluster_centers_, fmt='%f')
            centroids = k_means.cluster_centers_
            if not os.path.isfile(mesh_file) or overwrite :
                vert, norm, tri = coarseMolSurface(atoms_coord_centerd.tolist(),None)
                dae = exportCollada(fname, vert,tri,norm,mesh_file)
                exportOBJ(fname, vert,tri,norm,geom_folder+os.sep+fname+"_cms.obj")
            #msms ?
            np.savetxt(data_folder+os.sep+fname+"_cl.txt",k_means.cluster_centers_)
            center = center.tolist()
            R = np.linalg.norm(centroids,axis=1).max()+proxy_radius
    else :
        if fname.find("biomt") != -1 :
            center = [0,0,0]
        else :
            s = p.get_structure(structure_id, filename)
            atoms_coord=[]
            for m in s.get_models():
                coords = [atom.coord.tolist() for atom in m.get_atoms() if atom.parent.resname != "DUM"]
                atoms_coord.extend(coords)
            center = np.average(atoms_coord,axis=0)
            if not os.path.isfile(mesh_file) or overwrite :
                atoms_coord_centerd = np.array(atoms_coord) - center
                vert, norm, tri = coarseMolSurface(atoms_coord_centerd.tolist(),None)
                dae = exportCollada(fname, vert,tri,norm,mesh_file)
                exportOBJ(fname, vert,tri,norm,geom_folder+os.sep+fname+"_cms.obj")
                #write as dae
            center = center.tolist()
        centroids = np.loadtxt(data_folder+os.sep+fname+"_cl.txt")
        nProxy = len(centroids)
        R = np.linalg.norm(centroids,axis=1).max()+proxy_radius
    #return [centroids.tolist(),], [(np.ones(nProxy)*proxy_radius).tolist(),], R, mesh_file, fname, center
    return [[[0,0,0],],centroids.tolist()], [[R],(np.ones(nProxy)*proxy_radius).tolist()], R, mesh_file, fname, center
    
def addIngredient(name, label, pdbId, surface,mw,count,stch,uniId):
    from autopack.Ingredient import MultiSphereIngr
    bu = 0 
    testpdb = pdbId.split("_")[0]
#    if name in summary:
#       bu = int(summary[name]["BU"])
#    else :
#       bu = getBU(testpdb,stch)  
    if not surface :
        bu = getBU(testpdb,stch)  
        if len(pdbId) == 4 and bu != 0:
            pdbId = pdbId+"_"+str(bu)
    print "build ingr ",name, pdbId, stch#, bu
    score = 1
    source = ""
    testpdb = pdbId.split("_")[0]
    if len(testpdb) == 4 and testpdb!= "none":
        score, source = GetScoreAndSource(label,uniId,pdbId.split("_")[0])
    pos,rad,R, mesh_file, pdbId, center  = buildProxy(pdbId,11.85,name,mw,stch,surface=surface,biomt=False) 
    if surface :
        if len(pdbId) == 4:
            pdbId=pdbId+"_mb"
#        if len(testpdb) == 4 and pdbId.lower() !="none" and pdbId.find("_mb")==-1:
#            pdbId=pdbId+"_mb"
    #check biomt statuts
#    if pdbId != "none":
#        if len(pdbId)>=4:
#            getPDBDescription(pdbId[:4])
#            print pdbId,stch,pdb_dic[pdbId[:4]][1]
    print "build ingr ",name, pdbId, stch,mesh_file
    ingr = MultiSphereIngr(name=name, pdb=pdbId, nbMol=count, positions=pos,
                               radii=rad,encapsulatingRadius=R)
    ingr.encapsulatingRadius=R
    ingr.principalVector = [0,0,1]
    if surface :
        ingr.offset = [center[0],center[1],center[2]]
        ingr.jitterMax = [0.5,0.5,0.00]
        #print ingr.offset
    ingr.description = label
    ingr.organism = source
    ingr.score = score
    ingr.bu = bu
    ingr.sphereFile = data_folder+os.sep+pdbId+"_cl.txt"
    ingr.meshFile = mesh_file
    ingr.meshName = pdbId
    #ingr.source["biomt"] = False
    return ingr

def changeSphFile(ingr):
    ingr.sphereFile=None               
#

def specialCase(env):
    special_props={}
    Hu_props = {}
    Hu_props["name"] = "mpn529"
    Hu_props["properties"]={"pairs":[[2],[1,2,3],[0,1],[3,4],[1,2,3],[2]],"beadsin":0,"beadsout":7,"range":[2,0],"marge_out":[82,86],"diehdral":[40,73],"marge_in":[74,94],"length":45.859,"bend":True,"pt4":[5.931,11.516,99.092],"pt3":[1.775,25.153,1.866],"pt2":[-2.677,-22.192,-4.602],"pt1":[-16.547,-73.938,77.388]}
    Hu_props["sphereFile"]=None
    Hu_props["positions"] = [(5.637, 20.943, 13.172),
                            (5.957, 7.007, -4.044),
                            (-4.645, -3.305, -17.916),
                            (-7.46, -4.494, 1.825),
                            (0.723, -11.892, 20.555)]
    Hu_props["radii"] =[11.85,11.85,11.85,11.85,11.85]
    Hu_props["partners_name"] = ["DNA"]
    
    special_props["mpn529"] = Hu_props
    
    RNAP_props = {}
    RNAP_props["name"] = "mpn191"
    RNAP_props["properties"]={"pairs":[[1],[1],[1],[2],[2],[2]],"st_pt2":[38.786,22.058,27.949],"st_pt1":[3.087,22.485,13.395],"b_pt1":[0],"marge_out":[60,70],"diehdral":[70,90],"marge_in":[0,5],"length":63.423051,"bend":True,"pt4":[-52.616,74.174,60.161],"pt3":[-1.663,11.179,-1.802],"pt2":[-45.827,-18.367,32.825],"pt1":[-119.96,-59.814,89.309],"st_ingr":"mRNA"}
#    RNAP_props["sphereFile"]="2o5i_1.pdb_kmeans3.txt"
    RNAP_props["partners_name"] = ["DNA"]
    special_props["mpn191"] = RNAP_props
    
    RIB_props = {}
    RIB_props["name"] = "mpn208"
    RIB_props["properties"]={"pairs":[[1],[1],[1],[2],[2],[2]],"b_pt1":[0],"st_pt2":[-17.79,95.035,-90.413],"st_pt1":[-3.29,20.055,-26.37],"marge_out":[24,28],"diehdral":[67,71],"marge_in":[48,50],"length":28.348752,"bend":True,"pt4":[72.775,-54.869,-78.361],"pt3":[-10.252,-39.096,-24.904],"pt2":[-38.275,-39.914,-20.695],"pt1":[-100.506,-118.001,-15.25],"st_ingr":"peptides"}    
#    RIB_props["sphereFile"]="1TWT_1TWV.pdb_kmeans3.txt"
    RIB_props["partners_name"] = ["mRNA"]
    special_props["mpn208"] = Hu_props
    
    for k in special_props:
        ingr = env.getIngrFromName(k)
        ingr.properties =   special_props[k]["properties"]      
#        ingr.sphereFile =   special_props[k]["sphereFile"]  
        ingr.partners_name =   special_props[k]["partners_name"]  
        if "positions" in special_props[k]:
            if len(ingr.positions)==2:
                ingr.positions[1] =   special_props[k]["positions"]   
                ingr.radii[1] =   special_props[k]["radii"]
            else :
                ingr.positions[0] =   special_props[k]["positions"]   
                ingr.radii[0] =   special_props[k]["radii"]   
        if "sphereFile"in special_props[k]:
            ingr.sphereFile =   special_props[k]["sphereFile"]   
        #ingr.nbMol = -ingr.nbMol # this should let cellpack skip it ? or manually do it
    
    #FIBER
    from autopack.Ingredient import GrowIngrediant
    #DNA
    kw={"packingMode":"random",
       "partners_position":[[[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]],[[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]]],
       "orientation":[1,0,0],"weight":0.2,"color":[1,0.498,0.314],"coordsystem":"left","nbJitter":10,"proba_not_binding":0.5,
       "radii":[[11.85]],"cutoff_boundary":55,"encapsulatingRadius":51,"positions2":[[[51,0,0]]],"rejectionThreshold":30,
       "isAttractor":False,"Type":"Grow","useLength":False,"uLength":102,
        "source":{"mWeight":799723980,"pdb":"dna_single_base","transform":{"center":False}},
        "gradient":"","use_rbsphere":51,"jitterMax":[1,1,1],"packingPriority":-200,
        "orientBiasRotRangeMax":-3.1415927,"closed":True,"useHalton":True,"overwrite_nbMol_value":30,
        "biased":0.5,"molarity":0,"useOrientBias":False,"rotRange":6.2831,"orientBiasRotRangeMin":-3.1415927,
        "marge":25,"perturbAxisAmplitude":0.1,"nbMol":1,"use_mesh_rb":False,"principalVector":[1,0,0],"properties":{"lnSeg":3338},
        "partners_name":[Hu_props["name"],RNAP_props["name"]],"name":"DNA",
        "positions":[[[-51,0,0]]],"excluded_partners_name":[],"placeType":"jitter","cutoff_surface":55,"length":5126954,
        "walkingMode":"sphere","constraintMarge":False,"compMask":[-2],"proba_binding":0.5,"pdb":"dna_single_base",
        "startingMode":"file","useRotAxis":False}
    ingr = GrowIngrediant(**kw)
    ingr.score =""
    ingr.organism = ""
    ingr.bu=0
    env.compartments[0].innerRecipe.addIngredient(ingr)
    #mRNA
    kw={"packingMode":"random","partners_position":[[[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]]],"orientation":[1,0,0],"weight":1,
        "color":[1,0.498,0.314],"coordsystem":"left","nbJitter":1,"proba_not_binding":0.5,"radii":[[6]],
        "cutoff_boundary":55,"encapsulatingRadius":50,"positions2":[[[50,0,0]]],"rejectionThreshold":300,
        "isAttractor":False,"Type":"Grow","useLength":False,"uLength":100,
        "source":{"mWeight":0.66,"pdb":"RNA_G_Base","transform":{"center":False}},"gradient":"","use_rbsphere":50,"jitterMax":[1,1,1],
        "packingPriority":-20,"orientBiasRotRangeMax":-3.1415927,"closed":False,"useHalton":True,"overwrite_nbMol_value":10,
        "biased":0.5,"molarity":0,"useOrientBias":False,"rotRange":6.2831,"orientBiasRotRangeMin":-3.1415927,"marge":80,
        "perturbAxisAmplitude":0.1,"nbMol":0,"use_mesh_rb":False,"principalVector":[1,0,0],"properties":{"lnSeg":-1,"ratioSeg":0.333333},
        "partners_name":[RIB_props["name"]],"name":"mRNA","positions":[[[-50,0,0]]],"excluded_partners_name":[],"placeType":"jitter",
        "cutoff_surface":25,"length":3000,"walkingMode":"sphere","constraintMarge":False,"compMask":[-2],
        "proba_binding":0.5,"pdb":"RNA_G_Base","startingMode":"line","useRotAxis":False}
    ingr = GrowIngrediant(**kw)
    ingr.score =""
    ingr.organism = ""
    ingr.bu=0
    env.compartments[0].innerRecipe.addIngredient(ingr)
    #peptide
    kw={"packingMode":"random","partners_position":[],"orientation":[1,0,0],"weight":0.2,"color":[1,0.498,0.314],
            "coordsystem":"left","nbJitter":1,"proba_not_binding":0.5,"radii":[[3]],"cutoff_boundary":55,"encapsulatingRadius":10,"positions2":[[[10,0,0]]],
            "rejectionThreshold":300,"isAttractor":False,"Type":"Grow","useLength":False,"uLength":20,
            "source":{"mWeight":0,"pdb":"alanine","transform":{"center":True}},"gradient":"","use_rbsphere":10,"jitterMax":[1,1,1],
            "packingPriority":-1,"orientBiasRotRangeMax":-3.1415927,"closed":False,"useHalton":True,"overwrite_nbMol_value":5,
            "biased":0.5,"molarity":0,"useOrientBias":False,"rotRange":6.2831,"orientBiasRotRangeMin":-3.1415927,"marge":120,
            "perturbAxisAmplitude":0.1,"nbMol":0,"use_mesh_rb":False,"principalVector":[1,0,0],"properties":{},"partners_name":[],
            "name":"peptides","positions":[[[-10,0,0]]],"excluded_partners_name":[],"placeType":"jitter","cutoff_surface":25,
            "length":500,"walkingMode":"sphere","constraintMarge":False,"compMask":[-2],"proba_binding":0.5,"startingMode":"line","useRotAxis":False}
    ingr = GrowIngrediant(**kw)
    ingr.score =""
    ingr.organism = ""
    ingr.bu=0
    env.compartments[0].innerRecipe.addIngredient(ingr)
    return env
          
def buildRecipe(dData):
    from upy import hostHelper
    import autopack
    autopack.helper = hostHelper.Helper()
    autopack.helper.host = "none"
    autopack.forceFetch = False

    from autopack.Environment import Environment
    from autopack.Compartment import Compartment
    from autopack.Recipe import Recipe
    from autopack.Ingredient import MultiSphereIngr
    env = Environment(name="Mpn")
    #two compartments - 1 surface, 1 interior 
    osurf = Compartment("MycoplasmaSurface", None, None, None)
    print ("added compartment ", "MycoPnsurface")
    osurf.ghost = True
    osurf.filename = "D:\\Data\\cellPACK_data\\Mycoplasma\\mpn_surface_center.dae"
    osurf.name = "mpn_center"
    osurf.ref_obj = "mpn_surface_center"
    env.addCompartment(osurf)
    osurf.ghost = False
#    oin = Compartment("MycoplasmaInner", None, None, None)
#    print ("added compartment ", "MycoPnInner")
#    oin.ghost = True
#    oin.filename = "D:\\Data\\cellPACK_data\\Mycoplasma\\mpn_surface_inner.dae"
#    oin.name = "mpn_surface_inner"
#    oin.ref_obj = "mpn_surface_inner"
#    env.addCompartment(oin)
#    oin.ghost = False
    ofoot = Compartment("MycoplasmaFoot", None, None, None)
    print ("added compartment ", "MycoPnFoot")
    ofoot.ghost = True
    ofoot.filename = "D:\\Data\\cellPACK_data\\Mycoplasma\\mpn_foot.dae"
    ofoot.name = "mpn_foot"
    ofoot.ref_obj = "mpn_foot"
    env.addCompartment(ofoot)
    ofoot.ghost = False
    #surface    
    rSurf = Recipe(name="surface")
#    o.setSurfaceRecipe(rSurf)    
    
    rMatrix = Recipe(name="interior")
#    o.setInnerRecipe(rMatrix)
    
    #loop through the entry , entry is the geneId
    for entry in dData :
        name = entry
        uniId = dData[entry]["uniprot"]
        label = dData[entry]["Name"]
#        cname = dData[entry]["CompName"]
#        if cname != "":
#            name = cname
        stch = dData[entry]["stoich"]
        if stch == "" :
            continue
        if int(eval(stch)) < 0 :
            continue
        location = dData[entry]["location"]
        surface = False
        if location == "TM" or location == "LIPO":
            surface = True
        #if location == "LIPO":
        #    continue
        partner = dData[entry]["partners"]
        pdb = dData[entry]["pdb"].strip()
        count =  int(eval(dData[entry]["count"]))
        mw = float(dData[entry]["protein mass/cell (Mda)"])
        if pdb == "":
            continue
        #if len(pdb) > 4 :
        #    print "problem pdb ",pdb
        #    continue
        ingr = addIngredient(entry, label, pdb, surface, mw, count,stch, uniId)
        if surface :
            rSurf.addIngredient(ingr)
        else :
            rMatrix.addIngredient(ingr)
    #check lefty
    osurf.setSurfaceRecipe(rSurf)    
    osurf.setInnerRecipe(rMatrix)
    
    env.boundingBox=[[0,0,0],[2000,2000,2000]]
    env.name = "MycoPn"
   
    env = specialCase(env)
   
    #save serialized
    from autopack.IOutils import serializedRecipe,saveResultBinary
    djson, all_pos, all_rot = serializedRecipe(env,False,True)#transpose, use_quaternion, result=False, lefthand=False
    with open("D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_1.0_serialized.json","w") as f:
        f.write(djson)
    with open("C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_1.0_serialized.json","w") as f:
        f.write(djson)       

    env.saveRecipe("D:\\Data\\cellPACK_data\\Mycoplasma\\MpnFlex_1.0.json",useXref=False,format_output="json")   
        
    env.loopThroughIngr(changeSphFile)
    env.saveRecipe("D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_1.0.json",useXref=False,format_output="json")   
    
    return env
#
#data,header = getTaxonomyProtein("mycoplasma mycoides")
#res = gatherInfo(data,header)

#data,header = getTaxonomyProtein('"Influenza A virus (strain A/Udorn/307/1972 H3N2) [381517]" AND reviewed:yes')
#res = gatherInfo(data,header)

DHeader = ["Stable ID","Protein name", "Control abundance (copies/cell)","protein mass/cell (Mda)","coord source","stoich gene","partners","location","complexes	models"]
import csv
csvfile = open('D:\\Data\\cellPACK_data\\Mycoplasma\\MyMpn_proteome_afterHomomerCheck.csv')
reader = csv.DictReader(csvfile)
DData = OrderedDict()
for row in reader:
    if row['Name']!='' :
        DData[row['ID'].lower()] = OrderedDict(row)
        DData[row['ID'].lower()]["partOf"] = ""
csvfile.close()

#read the complex
csvfile = open('D:\\Data\\cellPACK_data\\Mycoplasma\\BookComplexs.csv')
reader = csv.DictReader(csvfile)
complexDData = {}
mapping={}
for row in reader:
    #split the raw using []
    #if int(row["Attachment"]) == 1 : continue
    elem = row["Genes"].split()
    gnames = []
    for el in elem :
        if el.startswith("MPN") :
            gname =  el.strip().lower()
            gnames.append(gname)
            if gname in DData:
                if DData[gname]['Name'].find(row['Complex name']) != -1 :
                    continue
            if gname not in mapping:
                mapping[gname] = [[row['Complex name'],row['Type']],]
            else :
                mapping[gname].append([row['Complex name'],row['Type']])
    if row['Complex name'] not in complexDData:
        complexDData[row['Complex name']]=[row['Complex name'],row['Type'],gnames]
    else :
        complexDData[row['Complex name']][2].extend(gnames)
csvfile.close()

#if os.path.isfile("D:\\Data\\cellPACK_data\\Mycoplasma\\pdb_dic.json"):
#    with open("D:\\Data\\cellPACK_data\\Mycoplasma\\pdb_dic.json", 'r') as fp:  # doesnt work with symbol link ?
#        pdb_dic = json.load(fp)
#gather the PDB template from PMP
#id, name, pdb, organism, score

#compare
import json
if not os.path.isfile("D:\\Data\\cellPACK_data\\Mycoplasma\\uniprot_data.json"):  
    taxonomyid = 272634
    taxonomy='"Mycoplasma pneumoniae (strain ATCC 29342 / M129) [272634]" AND reviewed:yes'
    data,header = getTaxonomyProtein(taxonomy)
    res = gatherInfo(data,header)
    #write csv compare to david csv
    with open("D:\\Data\\cellPACK_data\\Mycoplasma\\uniprot_data.json","w") as fp:  # doesnt work with symbol link ?
         json.dump(res, fp, indent=1, separators=(',', ':'))  # ,indent=4, separators=(',', ': ')
    uniData = res
#    f = open("D:\\Data\\cellPACK_data\\Mycoplasma\\resume.txt","w")
#    #header
#    f.write("Id;Name;PDB;Count;Template;stochiometryD;components;localisation;loca_u\n")
#    for d in DData :
#        if (DData[d]["Control abundance (copies/cell)"]==''): continue
#        #if int(DData[d]["Control abundance (copies/cell)"]) < 20 : continue
#        pdb = DData[d]['PDB']
#        if pdb =="":
#            pdb = "none"
#        stochiometryU = "None"
#        stochiometryD = DData[d]['stoich']
#        stochiometryD = stochiometryD if stochiometryD !="" else "None"
#        localisationD = DData[d]["location"] if DData[d]["location"] !="" else "None"
#        localisationU = "None"
#        if d in res:
#            stochiometryU = res[d]["udic"][header.index('Subunit structure [CC]')]
#            stochiometryU = stochiometryU if stochiometryU !="" else "None"
#            localisationU = res[d]["udic"][header.index("Subcellular location [CC]")] if res[d]["udic"][header.index("Subcellular location [CC]")]!="" else "None"
#            res[d]["udic"][header.index('Subunit structure [CC]')]
#            #print d," PDB ",DData[d]['PDB']," Count ",DData[d]["Control abundance (copies/cell)"], res[d]['pdb']
#            f.write(d+";"+res[d]["label"].replace(";",",")+";"+pdb+";"+str(DData[d]["Control abundance (copies/cell)"])+";"+res[d]['pdb'].replace(";",",")+";"+stochiometryD+";"+stochiometryU+";"+localisationD+";"+localisationU+"\n")
#        else :
#            #print d," PDB ",DData[d]['PDB']," Count ",DData[d]["Control abundance (copies/cell)"], 'None'
#            f.write(d+";"+DData[d]["Name"].replace(";",",")+";"+pdb+";"+str(DData[d]["Control abundance (copies/cell)"])+";None;"+stochiometryD+";"+stochiometryU+";"+localisationD+";"+localisationU+"\n")
#    f.close()
    
else :
    with open("D:\\Data\\cellPACK_data\\Mycoplasma\\uniprot_data.json", 'r') as fp:  # doesnt work with symbol link ?
        uniData = json.load(fp)
#    csvfile2 = open("D:\\Data\\cellPACK_data\\Mycoplasma\\resume.txt","r")
#    reader2 = csv.DictReader(csvfile2,delimiter=';',)
#    uniData = {}
#    dowrite = True
#    field = None
#    for row in reader2:
#        #[None, 'Id;Name;PDB;Count;Template;stochiometryD;components;localisation;loca_u']
#        uniData[row['Id']]=row
#        field = row.keys()
#        #add the complex information
#        if 'complex' not in uniData[row['Id']]:
#            uniData[row['Id']]['complex']=[]
#            if row['Id'] in mapping:
#                dowrite = False
#                uniData[row['Id']]['complex']=mapping[row['Id']]
#    csvfile2.close()
#    del reader2
    
    
    #loop over the complex
##allres={}
##for c in complexDData:
##    allres[c] = filterTemplatesComplex(complexDData[c],uniData)
##
#def testOneTemplates(akey):
#    if akey not in uniData:
#        return
#    templates = uniData[akey]['Template'].split(",")
#    label = uniData[akey]["Name"].lower()
#    return pickAModelOnDescription(templates, label)
        

##corss and verify data with uniprot e.g. localisation and template
for entry in DData:
    if entry not in uniData:
        print entry," not in unidata"
    name = entry
    label = DData[entry]["Name"]
    stch = DData[entry]["stoich"]
    DData[entry]["uniprot"] = ""
    if entry in mapping :
        DData[entry]["partOf"] = mapping[entry]
    if stch == "" :
        continue
    if int(eval(stch)) < 0 :
        continue
    location = DData[entry]["location"]
    if entry in uniData:
        localisation = uniData[entry]['udic']['Subcellular location [CC]']
        if location == "":
            if localisation != "" :
                t = uniData[entry]['udic']['Subcellular location [CC]'].lower().split("note")[0].find("membrane")
                if t != -1 :
                    print "****************",entry,uniData[entry]['udic']['Subcellular location [CC]']
                    location = "TM"   
                    DData[entry]["location"] = "TM"
        surface = False
        DData[entry]["uniprot"] = uniData[entry]['udic']['Entry']
    if location == "TM" or location == "LIPO":
        surface = True
    #if location == "LIPO":
    #    continue
    partner = DData[entry]["partners"]
    pdb = DData[entry]["pdb"].strip()
    count =  int(eval(DData[entry]["count"]))
    mw = float(DData[entry]["protein mass/cell (Mda)"])

    #check in uniprot
    if pdb == "" :
        if entry in uniData:
            templates = uniData[entry]['pdb'].lower().split(",")
            label = uniData[entry]["label"].lower()
            if len(templates)!=0:
                pdb = templates[0]#pickAModelOnDescription(templates, label)
                DData[entry]["pdb"] = pdb
#    description = ""
#    if pdb != "" and pdb.lower() != "none":
#        if len(pdb.split("_")[0])==4:        
#             description = getPDBDescription(pdb[:4])#polymer,name,len(chains),accession->result,nAssambly,info
#    DData[entry]["pdb_description"] = description
    

#write the DDAta
env = buildRecipe(DData)

if not os.path.isfile("D:\\Data\\cellPACK_data\\Mycoplasma\\pdb_dic.json"):
    with open("D:\\Data\\cellPACK_data\\Mycoplasma\\pdb_dic.json", 'w') as fp:  # doesnt work with symbol link ?
        json.dump(pdb_dic,fp,indent=1, separators=(',', ':'))
    
    
for i in env.compartments[0].surfaceRecipe.ingredients:
    if i.source["pdb"] == "None":
        print i.name, i.source
for i in env.compartments[0].innerRecipe.ingredients:
    if i.source["pdb"] == "None":
        print i.name, i.source

#save a simple file with

def printSTOICH(ingr):
    pdb = ingr.source["pdb"]
    ingr.bu = 0
    if len(pdb) == 4 and pdb.lower() !="none":
        nBU=pdb_dic[pdb][1]
        print pdb
        for i in range(nBU+1):
            info = getBioAssamblyInfo(pdb,bio=i)
            if len(info['symmetries']):
                pdb_stch=1
                pdb_stch_str = info['symmetries'][0]['stoichiometry'].split("A")[1].split("B")[0]#problem with A2B2
                #print pdb,i,pdb_stch_str,DData[ingr.o_name]["stoich"],info['symmetries'][0]['stoichiometry']
                if pdb_stch_str!='' and len(pdb_stch_str) != 0:
                    pdb_stch = int(pdb_stch_str)
#                print pdb_stch,pdb_stch_str
                if pdb_stch == int(DData[ingr.o_name]["stoich"]) :
                    print "use bio ",i," for ", pdb, pdb_stch_str,pdb_stch,DData[ingr.o_name]["stoich"] 
                    ingr.bu = i
                    return
 

#id, name, pdb, organism, score
k=["ID","PDB","BU","ORG","SCORE","LOC","STOICH","Name"]
with open('D:\\Data\\cellPACK_data\\Mycoplasma\\SummaryProteome.csv','wb') as fou:
    writer = csv.DictWriter(fou,k,delimiter=",")
    writer.writeheader()
    for i in env.compartments[0].surfaceRecipe.ingredients:
        writer.writerow({"ID":i.o_name,"PDB":i.source["pdb"],"BU":i.bu, "ORG" : i.organism,"SCORE": i.score,"LOC":"membrane","STOICH":DData[i.o_name]["stoich"], "Name":DData[i.o_name]["Name"]})
    for i in env.compartments[0].innerRecipe.ingredients:
        writer.writerow({"ID":i.o_name,"PDB":i.source["pdb"],"BU":i.bu, "ORG" : i.organism,"SCORE": i.score, "LOC":"cytoplasm","STOICH":DData[i.o_name]["stoich"],"Name":DData[i.o_name]["Name"]})

k=['ID',
 'Name',
 'count',
 'label',
 'location',
 'name gene',
 'partners',
 'pdb',
 'protein mass/cell (Mda)',
 'source_coord',
 'stoich',
 'to do',
 'unit Control abundance (copies/cell)',
 'w',
 'partOf',
 "uniprot"
]
with open('D:\\Data\\cellPACK_data\\Mycoplasma\\MyMpn_proteome_L.csv','wb') as fou:
    writer = csv.DictWriter(fou,k,delimiter=",")
    writer.writeheader()
    for d in DData:
        writer.writerow(DData[d])
#for gname in complexDData[c][2]:
        #get the PDB ? the count ?
    #    print gname
    #    testOneTemplates(gname)       
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
# after runing autopack in c4d we can run 
#                from autopack.IOutils import serializedFromResult, toBinary, gatherResult
#d,p,r = serializedFromResult(result,False,True,result=True,lefthand=False)
#
#ap, ar = gatherResult(posrot, True, True,  lefthand=True)
#
#
#from autopack.IOutils import serializedRecipe, saveResultBinary, toBinary
#env = c4d.af.values()[0].histoVol.values()[0]
#djson, all_pos, all_rot = serializedRecipe(env,True,True,True,True)#env, transpose, use_quaternion, result=False, lefthand=False
#with open("D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_1.0_serialized.json","w") as f:
#    f.write(djson)
#if resultfile!= "" :
#    saveResultBinary(env,"D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_1.0_serialized.bin",False,True, lefthand=True)   
#from autopack.Serializable import sCompartment,sIngredientGroup,sIngredient,sIngredientFiber
#sCompartment.static_id = 0;sIngredientFiber.static_id = 0;sIngredient.static_id = [0, 0, 0];sIngredientGroup.static_id= 0;djson, all_pos, all_rot = serializedRecipe(env,True,True,True,True);toBinary(all_pos, all_rot,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_1.0_serializedTR_L.bin")     
#sCompartment.static_id = 0;sIngredientFiber.static_id = 0;sIngredient.static_id = [0, 0, 0];sIngredientGroup.static_id= 0;djson, all_pos, all_rot = serializedRecipe(env,False,True,True,True);toBinary(all_pos, all_rot,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_1.0_serialized_L.bin")     
#sCompartment.static_id = 0;sIngredientFiber.static_id = 0;sIngredient.static_id = [0, 0, 0];sIngredientGroup.static_id= 0;djson, all_pos, all_rot = serializedRecipe(env,False,True,True,False);toBinary(all_pos, all_rot,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_1.0_serialized.bin")     
#sCompartment.static_id = 0;sIngredientFiber.static_id = 0;sIngredient.static_id = [0, 0, 0];sIngredientGroup.static_id= 0;djson, all_pos, all_rot = serializedRecipe(env,True,True,True,False);toBinary(all_pos, all_rot,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_1.0_serializedTR.bin")     