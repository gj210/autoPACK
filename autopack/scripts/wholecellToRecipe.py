# -*- coding: utf-8 -*-
"""
Created on Wed May 11 13:23:47 2016

@author: ludov
"""
import h5py
import json
import sys
import os
import math
import urllib2
import requests

import numpy as np
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs\lib-tk")

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
                
f = h5py.File("C:\\Users\\ludov\\Downloads\\simulation-1.h5")
dataMono = f.get("states/ProteinMonomer/counts/data")
labelMono = f.get("states/ProteinMonomer/counts/labels/0")
labelCompartments = f.get("states/ProteinMonomer/counts/labels/1")
dataComplex = f.get("states/ProteinComplex/counts/data")
labelComplex = f.get("states/ProteinComplex/counts/labels/0")

frame = 0
with open("C:\\Users\\ludov\\Downloads\\rawrecipe_frame_%i.json" % frame, 'r') as fp:  # doesnt work with symbol link ?
    recipe = json.load(fp)
with open("C:\\Users\\ludov\\Downloads\\rawcross.json", 'r') as fp:  # doesnt work with symbol link ?
    cross = json.load(fp)

data_folder = "D:\\Data\\cellPACK_data\\cellPACK_database_1.1.0\\other"
import autopack
from autopack.Environment import Environment
from autopack.Compartment import Compartment
from autopack.Recipe import Recipe
from autopack.Ingredient import MultiSphereIngr

from upy import hostHelper

autopack.helper = hostHelper.Helper()
autopack.helper.host = "none"
autopack.forceFetch = False

from scipy.cluster.vq import kmeans,vq
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBList import PDBList
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa

fetch = PDBList(pdb=data_folder)
p = PDBParser(PERMISSIVE=1)

def getMWFromSequence(sequence):
    X = ProteinAnalysis(sequence)
    mw = X.molecular_weight()  
    return mw

def getSequenceStructure(s):
    seq = ""
    for r in s.get_residues():
        if is_aa(r.get_resname(), standard=True):
            seq+=three_to_one(r.get_resname())
        else:
            seq+="G"  
    return seq
    
def getRadiusFromMW(mw):
    V = mw * 1.21    
    return math.pow((3*V)/(4*math.pi),1.0/3.0)

def coarseMolSurface(coords, radii, XYZd =[32,32,32],isovalue=6.0,resolution=-0.1,padding=0.0,
                         name='CoarseMolSurface',geom=None):
    from UTpackages.UTblur import blur
#        print "res",resolution
    if radii is None :
        radii = np.ones(len(coords))*1.8
    volarr, origin, span = blur.generateBlurmap(coords, radii, XYZd,resolution, padding = 0.0)
    volarr.shape = (XYZd[0],XYZd[1],XYZd[2])
#        print volarr
    volarr = np.ascontiguousarray(np.transpose(volarr), 'f')
    weights =  np.ones(len(radii), typecode = "f")
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
    if data.dtype.char!=np.Float32:
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
    
def pickAmodel(length):
    maxLength = 0
    longest = -1
    for i in range(len(length)):
        ln = eval(length[i][1])-eval(length[i][0])
        if ln > maxLength :
            maxLength = ln
            longest = i
    return longest
    
#1KZN
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
    #http://www.proteinmodelportal.org/?pid=modelDetail&provider=MODBASE&template=1i4jA&pmpuid=1000812646454&range_from=1&range_to=144&ref_ac=P47402&mapped_ac=P47402&zid=async
    #<a href="?pid=modelDetail&provider=MODBASE&template=1i4jA&pmpuid=1000812646454&range_from=1&range_to=144&ref_ac=P47402&mapped_ac=P47402" target="_blank" class="externlink">[Show]</a>
#    http://www.proteinmodelportal.org/?pid=modelDetail&provider=MODBASE&template=1i4jA&pmpuid=1000812646454&range_from=1&range_to=144&ref_ac=P47402&mapped_ac=P47402&zid=async


#filename=data_folder+"\\MG_101_MONOMER-bound.pdb"
#from ghost import Ghost
#g = Ghost(wait_timeout=20)
#session = g.start()
#def browsePPM():
#    url = "http://sunshine.phar.umich.edu/server.php"
#    page, resources = session.open(url)
#    payload = {"submit":"Submit","inout":"in","yesno":"no","userfile":filename}
#    session.wait_timeout = 999
#    #value, resources = session.evaluate( 
#    #       'document.getElementById("inout").checked') 
#    value, resources = session.evaluate('validate()')
#    values = {"userfile":filename}
#    session.fill('form',values)
#    session.call('form','submit',expect_loading=True)    
#    p,r = session.wait_for_page_loaded()
#    #document.getElementById('myform').submit();
#    value, resources = session.evaluate("document.main.submit();")
#server.php?
    
def computeOPM(filename,pdbId):
    if not os.path.isfile(data_folder+os.sep+pdbId+"_mb.pdb"):  
        url = "http://sunshine.phar.umich.edu/upload_file.php"
        payload = {"submit":"Submit","inout":"in","yesno":"no"}#,"userfile":filename}
        r = requests.post(url, data=payload,files={"userfile":open(filename,"r")},timeout=120.0)
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
        print p.stored    
        if not len(p.stored) :
            return ""
        url = "http://opm.phar.umich.edu/"+p.stored[0][0]
        print url
        #direct access
        try :
            response = urllib2.urlopen(url)
            res = response.read()
            f = open(data_folder+os.sep+pdbId+"_mb.pdb","w")
            f.write(res)
            f.close()
            return data_folder+os.sep+pdbId+"_mb.pdb"
        except:
            #not found
            #build
            print ("problem accessing " + url)
            return ""
    else :
        return data_folder+os.sep+pdbId+"_mb.pdb"

def buildProxy(PDBid,proxy_radius,name,surface=False, overwrite=False):
    # read the pdb
    # build cluster
    # pass pos/radii
    # check for CA only ?
    overwrite = True
    name = name.split("-")[0]
    mesh = []
    structure_id = PDBid
    center=[0,0,0]
    if len(PDBid) == 4 :
        if surface :
            filename = getOPM(pdbId)       
            if filename == "":
                temp_filename = fetch.retrieve_pdb_file(structure_id)
                filename = computeOPM(temp_filename,PDBid)
        else :
            filename = fetch.retrieve_pdb_file(structure_id)
        fname = PDBid
    else :
        #parse from protein model portal
        filename = queryPMP(PDBid,name)
        fname = name
        if surface and filename != "":
           filename = computeOPM(filename,name)
        #http://salilab.org/modbase/retrieve/modbase?databaseID=P21812
    if filename == "" :
        #gather the radius at least
        mw = cross[name.split("-")[0]][-1]
        R = getRadiusFromMW(mw)
        return [[[0,0,0]]],[[R]],R,[],"",center
    if not os.path.isfile(data_folder+os.sep+fname+"_cl.txt") or overwrite :            
        s = p.get_structure(structure_id, filename)
        for m in s.get_models():
            atoms_coord = [atom.coord.tolist() for atom in m.get_atoms() if atom.parent.resname != "DUM"]
            break
        center = np.average(atoms_coord,axis=0)
        atoms_coord_centerd = np.array(atoms_coord) - center
        R = np.linalg.norm(atoms_coord_centerd,axis=1).max()
        Vproxy = 4*math.pi*(proxy_radius*proxy_radius*proxy_radius)/3.0
        V = len(atoms_coord) * 10.0 * 1.21
        nProxy = int(round(V/Vproxy))
        print "cluster ", nProxy, len(atoms_coord)
        if nProxy == 0:
            nProxy = 1
        # ~(1.21 x MW) A**3/molecule
        # V=4*math.pi*(r*r*r)/3.0
        r =  math.pow((3*V)/(4*math.pi),1.0/3.0)
        # r = math.pow(r3,1.0/3.0)
        #cluster
    #	0.73 cm**3/g  x  10**24 A**3/cm**3  x  molecular weight g/mole
    #	--------------------------------------------------------------
    #			 6.02 x 10**23 molecules/mole
        centroids,_ = kmeans(atoms_coord_centerd,nProxy)
        mesh = []#coarseMolSurface(atoms_coord_centerd.tolist(),None)
        #msms ?
        np.savetxt(data_folder+os.sep+fname+"_cl.txt",centroids)
        center = center.tolist()
    else :
        centroids = np.loadtxt(data_folder+os.sep+fname+"_cl.txt")
        nProxy = len(centroids)
        R = np.linalg.norm(centroids,axis=1).max()+proxy_radius
    return [centroids.tolist(),], [(np.ones(nProxy)*proxy_radius).tolist(),], R, mesh, fname, center
    
#X = ProteinAnalysis("MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV")
#print(X.count_amino_acids())
#print(X.get_amino_acids_percent())
#print(X.molecular_weight())
left=[]
env = Environment(name="MG")
#c Cytosol 
#d DNA 
#e Extracellular Space 
#m Membrane 
#tc Terminal Organelle Cytosol 
#tm Terminal Organelle Membrane 

#rCyto = Recipe()
## sorted(numbers, key=str.lower)
#for ing_name in sorted(ingrs_dic, key=unicode.lower):  # ingrs_dic:
#    # either xref or defined
#    ing_dic = ingrs_dic[ing_name]
#    ingr = io_ingr.makeIngredientFromJson(inode=ing_dic, recipe=env.name)
#    rCyto.addIngredient(ingr)
#    # setup recipe
#env.setExteriorRecipe(rCyto)

#put a score of confidence on structure
#Xray,NMR,cryo, modeling (score)

#c Cytosol 0
#d DNA     1
#e Extracellular Space  2
#m Membrane  3
#tc Terminal Organelle Cytosol  4
#tm Terminal Organelle Membrane  5

nIngrInterior = len(recipe["0"])
nIngrSurface = len(recipe["2"])+len(recipe["3"])
#need a filename or v,f,n
#o = Compartment("MycoGenit", None, None, None, filename=None,
#                object_name=None, object_filename=None)
o = Compartment("MycoGenit", None, None, None)
print ("added compartment ", "MycoGenit")
o.ghost = True
o.filename = "MycoplasmaGenitalium.dae"
o.name = "MycoplasmaGenitalium"
o.gname = "MycoplasmaGenitalium"
env.addCompartment(o)

rCyto = Recipe()
# oustide or surface ?
for i in range(len(recipe["2"])):
    pos=None
    rad=None
    R = 20.0
    wid, wnum, pdbId, count, full_name =  recipe["2"][i]
    if pdbId != "null" and pdbId != "" :
        pos,rad,R, mesh, pdbId, center = buildProxy(pdbId,10.0,wid,surface=True)
    else :
        left.append([wid, full_name,"2",i])
        print ("no information for ",wid,pdbId,count,wnum)
        continue
    ingr = MultiSphereIngr(name=wid, pdb=pdbId, nbMol=int(eval(count)), positions=pos,
                           radii = rad,encapsulatingRadius=R,principalVector=(0,0,1),
                            offset = center, jitterMax=[1,1,0])
    rCyto.addIngredient(ingr)
    print pdbId
env.setExteriorRecipe(rCyto)

#surface    
rSurf = Recipe(name="surf_" + str(len(env.compartments) - 1))
for i in range(len(recipe["3"])):
    pos=None
    rad=None
    R = 20.0
    wid, wnum, pdbId, count, full_name =  recipe["3"][i]
    if pdbId != "null" and pdbId != "" :
        pos,rad,R, mesh, pdbId, center  = buildProxy(pdbId,10.0,wid,surface=True) 
    else :
        left.append([wid, full_name,"3",i])
        print ("no information for ",wid)
        continue
    ingr = MultiSphereIngr(name=wid, pdb=pdbId, nbMol=int(eval(count)), positions=pos,
                           radii=rad,encapsulatingRadius=R)
    rSurf.addIngredient(ingr)
    print pdbId
o.setSurfaceRecipe(rSurf)    

rMatrix = Recipe(name="int_" + str(len(env.compartments) - 1))
for i in range(len(recipe["0"])):
    pos=None
    rad=None
    R = 20.0
    wid, wnum, pdbId, count, full_name =  recipe["0"][i]
    if pdbId != "null" and pdbId != "" :
        pos,rad,R, mesh, pdbId, center  = buildProxy(pdbId,10.0,wid)
    else :
        left.append([wid, full_name, "0",i])
        print ("no information for ",wid)
        continue      
    #read the pdb
    #build proxy and mesh
    ingr = MultiSphereIngr(name=wid, pdb=pdbId, nbMol=int(eval(count)), positions=pos,
                           radii=rad,encapsulatingRadius=R)
    rMatrix.addIngredient(ingr)
    print pdbId
o.setInnerRecipe(rMatrix)

#check lefty

env.boundingBox=[[0,0,0],[2000,2000,2000]]
env.name = "MycoplasmaGenitalium"
env.saveRecipe("C:\\Users\\ludov\\Downloads\\MG_1.0.json",useXref=False,format_output="json")
#save serialized
from autopack.IOutils import serializedRecipe,saveResultBinary
djson = serializedRecipe(env)
f=open("C:\\Users\\ludov\\Downloads\\MG_1.0_serialized.json","w")
f.write(djson)
f.close()
saveResultBinary(env,"C:\\Users\\ludov\\Downloads\\MG_1.0_serialized.bin",False,True)

#MG4 - nascent 753
#MG4 - processed I 753
#MG4 - processed II 753
#MG4 - signal sequence 753
#MG4 - folded 776
#MG4 - mature 776
#MG4 - inactivated 776
#MG4 - bounded 776
#MG4 - misfolded 776
#MG4 - damaged 776

#MG3 - folded 776
#MG3 - nascent 753
#DNA_GYRASE-nascent 0-1
#DNA_GYRASE-mature 1-4
#DNA_GYRASE-inactivated 0
#DNA_GYRASE-bounded 87 #8H 50 after 1H
#DNA_GYRASE-misfolded 0
#DNA_GYRASE-damaged 0
#    
#MG14-15. 500 monomer -> 40 dimer ?
#    
#for i in range(len(recipe["3"])):
#    
#if "surface" in comp_dic:
#    snode = comp_dic["surface"]
#    ingrs_dic = snode["ingredients"]
#    if len(ingrs_dic):
#        rSurf = Recipe(name="surf_" + str(len(env.compartments) - 1))
#        #                        rSurf = Recipe(name=o.name+"_surf")
#        for ing_name in sorted(ingrs_dic, key=unicode.lower):  # ingrs_dic:
#            # either xref or defined
#            ing_dic = ingrs_dic[ing_name]
#            ingr = io_ingr.makeIngredientFromJson(inode=ing_dic, recipe=env.name)
#            rSurf.addIngredient(ingr)
#            # setup recipe
#        o.setSurfaceRecipe(rSurf)
#if "interior" in comp_dic:
#    snode = comp_dic["interior"]
#    ingrs_dic = snode["ingredients"]
#    if len(ingrs_dic):
#        #                        rMatrix = Recipe(name=o.name+"_int")
#        rMatrix = Recipe(name="int_" + str(len(env.compartments) - 1))
#        for ing_name in sorted(ingrs_dic, key=unicode.lower):  # ingrs_dic:
#            # either xref or defined
#            ing_dic = ingrs_dic[ing_name]
#            ingr = io_ingr.makeIngredientFromJson(inode=ing_dic, recipe=env.name)
#            rMatrix.addIngredient(ingr)
#            # setup recipe
#        o.setInnerRecipe(rMatrix)