import sys
import os
import math
#import c4d
import numpy as np
#append MGLTools
#WINDOWS
sys.path.append("C:\\Users\\ludov\\Downloads\\blender-2.77-windows64\\MGLToolsPckgs")
#sys.path.append("C:\\Users\\ludo\\Downloads\\mgltools_win_amd64_latest")
#sys.path.append("C:\\Users\\ludo\\AppData\\Roaming\\MAXON\\CINEMA\ 4D\ R17\ Demo_E0A949BC\\plugins\\ePMV\\mgl64\\MGLToolsPckgs")
#LINUX
#sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/")
#sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/PIL/")
from upy.dejavuTk.dejavuHelper import dejavuHelper as helperClass
import autopack

helper = autopack.helper
if helper is None :
    import upy
    #helperClass = upy.getHelperClass()
    helper = helperClass(vi="nogui")
autopack.helper = helper

from sklearn.cluster import MiniBatchKMeans, KMeans
from sklearn.cluster import AffinityPropagation
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs

try :
    import prody
except:
    prody = None
    
from autopack.Environment import Environment
#from autopack.Graphics import AutopackViewer as AFViewer
#
pathToCP="/opt/data/dev/cellPACK/cellPACK_data_git/cellPACK_database_1.1.0/"
pathToCP="D:\\DATA\\cellPACK_data\\cellPACK_database_1.1.0\\"
#D:\Data\cellPACK_data\cellPACK_database_1.1.0\other

# filename = "/home/ludo/hivexp/BloodHIV1.0.json"
filename = "/opt/data/dev/cellPACK/cellPACK_data_git/cellPACK_database_1.1.0/recipes/Mycoplasma1.6.json"
filename = "/opt/data/dev/cellPACK/cellPACK_data_git/cellPACK_database_1.1.0/recipes/DNAplectoneme.1.0.json"
filename = pathToCP+os.sep+"recipes"+os.sep+"Mycoplasma1.6_full.json"
#path = "/home/ludo/hivexp/"
fileName, fileExtension = os.path.splitext(filename)
n=os.path.basename(fileName)
h = Environment(name=n)
##h.helper = helper
#recipe=n
h.loadRecipe(filename)

#for each ingredient, PDB, make a sphere Tree using SBL
sbl_binary_U="sbl-ballcovor-pdb-U.exe"
sbl_binary_A="sbl-ballcovor-pdb-A.exe"
sbl_binary_C="sbl-ballcovor-pdb-C.exe"
sbl_binary_T="sbl-ballcovor-txt.exe"

sbl_arg_file=" -f "
sbl_arg_nballs=" --n-inner-balls "
sbl_arg_outer=" --outer "
sbl_arg_interpolated=" --interpolated "
sbl_arg_verbose=" -v "

pdb_directory=pathToCP+os.sep+"other"+os.sep

LOD_LEVELES=np.array([0.15,0.1,0.05,0.01,0.001,0.0])#0 is encapsulating radius

try :
    import prody
except:
    prody= None
    
import glob

def sql_cg(pdbname):
    if len(pdbname) == 4: #PDBID
        prody.fetchPDB(pdbname, compressed=False)
        pdbname+=".pdb"
    mol = prody.parsePDB(pdbname)
    na=mol.numAtoms()
    print ("num atoms is ",na)
    LOD_NATOMS = (na*LOD_LEVELES).astype(int)
    filename = pdb_directory+pdbname
    binary = sbl_binary_U
    sblname="sbl-ballcovor-pdb-U__inner_approximation.txt"
    vmdscript="*.vmd"
    for i in range(4):
        #os.system(binary+sbl_arg_file+filename+sbl_arg_nballs+str(LOD_NATOMS[i])+sbl_arg_outer+sbl_arg_interpolated+sbl_arg_verbose)
        os.system(binary+sbl_arg_file+filename+sbl_arg_nballs+str(LOD_NATOMS[i])+sbl_arg_verbose)
        data=np.loadtxt(sblname)
        #sqrt the radius and save
        if len(data)==0:
            print ("********************************")
            print (pdbname)
            print ("********************************")
            continue
        data[:,3] = np.sqrt(data[:,3])
        ofilename = pdbname+"_"+str(LOD_LEVELES[i])+'.txt'
        np.savetxt(ofilename, data, fmt='%g')
        print ("success ", ofilename)
        os.system("rm "+sblname)
        os.system("mkdir "+pdbname+"_sph")
        os.system("mv *.vmd "+pdbname+"_sph/.")



def prodyLoad(pdbname,biomt=False):
    #biomt?
    if len(pdbname) == 4: #PDBID
        prody.fetchPDB(pdbname.lower(), compressed=False)
        pdbname = pdbname.lower() + ".pdb"
    else :
        if pdbname[-4:] != ".pdb":
            pdbname += ".pdb"

    if biomt:
        mol,header = prody.parsePDB(pdbname, header=True)        
        if len(header['biomoltrans']):
            mol = prody.buildBiomolecules( header, mol)
    else:
        mol = prody.parsePDB(pdbname, header=False)
        
    na=mol.numAtoms()
    c=mol.getCoords()
    center_c = c - np.average(c,0)
    np.savetxt(pdb_directory+os.sep+pdbname+"_cl.txt",
               center_c, fmt='%f')
    return center_c
    
def doIngredientSPH(ingredient):
    if "pdb" in ingredient.source:
        print (ingredient.source["pdb"])
        sql_cg(ingredient.source["pdb"])

def unique_rows(data):
    uniq = np.unique(data.view(data.dtype.descr * data.shape[1]))
    return uniq.view(data.dtype).reshape(-1, data.shape[1])

def voxelize(pdbname,spacing=10.0,padding=0.0):
    #spacing=1/5.0
    #padding=5.0
    if len(pdbname) == 4: #PDBID
        prody.fetchPDB(pdbname, compressed=False)
        pdbname+=".pdb"
    mol = prody.parsePDB(pdbname)
    na=mol.numAtoms()
    c=mol.getCoords()
    center_c = c - np.average(c,0)
    bot=np.min(center_c,0)+padding
    top=np.max(center_c,0)+padding
    ijk = (1/spacing * (center_c)).astype(int)
    n=np.max(ijk)
    ijku=unique_rows(ijk)
    out_coords = (ijku*spacing)
    np.savetxt(pdbname+".vox", out_coords, delimiter=' ', fmt='%f')
    print ("success ", pdbname+".vox")

def voxelize_avg(pdbname,spacing=20.0,padding=20.0):
    #spacing=1/5.0
    #padding=5.0
    if len(pdbname) == 4: #PDBID
        prody.fetchPDB(pdbname, compressed=False)
        pdbname+=".pdb"
    mol = prody.parsePDB(pdbname)
    na=mol.numAtoms()
    c=mol.getCoords()
    center_c = c - np.average(c,0)
    bot=np.min(center_c,0)+padding
    top=np.max(center_c,0)+padding
    #new_center = top-bot

    ind = np.array((1/spacing*(center_c - bot)), 'int')
    maxi = np.max(ind, 0)
    mask = np.zeros( maxi+1 )
    #ind1 = [tuple(x.tolist()) for x in ind]
    mask[ [ind[:,0],ind[:,1],ind[:,2]] ] = 1

    ijk = (1/spacing * (center_c-bot)).astype(int)
    n=np.max(ijk)
    ijku=unique_rows(ijk)
    out_coords = (ijku*spacing)+bot
    np.savetxt(pdbname+".xyz", out_coords, delimiter=' ', fmt='%f')
    avg=[ijk[0].tolist(),]
    coords=[[center_c[0]],]
    for i in range(1,len(ijk)):
        found =False
        for j in range(len(avg)) :
            if ijk[i].tolist()==avg[j] :
                coords[j].append(center_c[i].tolist())
                found = True
                break
        if not found :
            avg.append(ijk[i].tolist())
            coords.append([center_c[i].tolist()])

    cavg=[]
    for c in coords:
       cavg.append(np.average(c,0).tolist())

    np.savetxt(pdbname+"_avg.xyz", cavg, delimiter=' ', fmt='%f')
    print ("success ", pdbname+".xyz")


def doCluster(pdbname,spacing=20.0,padding=20.0, percentile=0.001,biomt=False):
    if len(pdbname) == 4: #PDBID
        pdbname = pdbname.lower() + ".pdb"
    else :
        if pdbname[-4:] != ".pdb":
            pdbname += ".pdb"
    if prody is not None :
        center_c = prodyLoad(pdbname,biomt=biomt)
    else :
        data=np.loadtxt(pdb_directory+os.sep+pdbname+"_cl.txt")
        center_c = data - np.average(data,0)
    V1=6970.17465149966
    V2=20.5795262761155
    natom=V1/V2
    ncluster = int(math.ceil(len(center_c)/natom))
    print ("name ",pdbname,ncluster,len(center_c))
    if ncluster < 10:
        ncluster = 10
        doKMeans(pdbname,center_c,ncluster)
    #Kmean of affinity etc
    #doAffinity(pdbname,center_c)
#    doKMeans(pdbname,center_c,ncluster)
    
def doAffinity(pdbname,points):
    # Compute Affinity Propagation
    af = AffinityPropagation(preference=-50).fit(points)
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_
    n_clusters_ = len(cluster_centers_indices)
    print (n_clusters_)
    print (labels)
    np.savetxt(pdbname+"_kmeans2.txt",
               af.cluster_centers_, fmt='%f')
    np.savetxt(pdbname+"_cl.txt",
               points, fmt='%f')

def doKMeans(pdbname,points,ncluster):
    
    k_means = KMeans(init='k-means++', n_clusters=ncluster, n_init=10)
    k_means.fit(points)
    #k_means_labels = k_means.labels_
    #k_means_cluster_centers = k_means.cluster_centers_
    #k_means_labels_unique = np.unique(k_means_labels)
    np.savetxt(pdb_directory+os.sep+pdbname+"_kmeans3.txt",
               k_means.cluster_centers_, fmt='%f')
            
def doKmeans(pdbname,spacing=20.0,padding=20.0, percentile=0.001):
    if len(pdbname) == 4: #PDBID
        prody.fetchPDB(pdbname.lower(), compressed=False)
        pdbname = pdbname.lower() + ".pdb"
    else :
        if pdbname[-4:] != ".pdb":
            pdbname += ".pdb"
    #mol,header = prody.parsePDB(pdbname, header=True)
    mol = prody.parsePDB(pdbname, header=False)
    #if len(header['biomoltrans']):
    #    mol = prody.buildBiomolecules( header, mol)
    na=mol.numAtoms()

    c=mol.getCoords()
    center_c = c - np.average(c,0)
    #print int(round(len(center_c)*0.008))
    ncluster = int(round(len(center_c)*percentile))
    if mol.numAtoms('ca') == mol.numAtoms():
        ncluster *= 5
    if ncluster == 0:
        ncluster = int(na/10.0)
    if ncluster <= 6:
        ncluster *= 2
    if ncluster == 0:
        print (pdbname,"no cluster")
        return
    print (ncluster)

    k_means = KMeans(init='k-means++', n_clusters=ncluster, n_init=10)
    k_means.fit(center_c)
    k_means_labels = k_means.labels_
    k_means_cluster_centers = k_means.cluster_centers_
    k_means_labels_unique = np.unique(k_means_labels)
    np.savetxt(pdbname+"_kmeans2.txt",
               k_means.cluster_centers_, fmt='%f')
    np.savetxt(pdbname+"_cl.txt",
               center_c, fmt='%f')

def doIngredientVOX(ingredient):
    if "pdb" in ingredient.source:
        if ingredient.source["pdb"] is not None and ingredient.source["pdb"] != "None":
            print (ingredient.source["pdb"])
            biomt = False
            if "biomt" in ingredient.source :
                biomt= ingredient.source["biomt"]
            doCluster(ingredient.source["pdb"],biomt=biomt)

h.loopThroughIngr(doIngredientVOX)

#from autopack.IOutils import serializedRecipe, saveResultBinary
#djson = serializedRecipe(h, False, True, lefthand = True)
#f=open("/opt/data/dev/cellPACK/cellPACK_data_git/cellPACK_database_1.1.0/recipes/Mycoplasma1.6_serialized.json", "w")
#f.write(djson[0])
#f.close()

#h.saveRecipe("/home/ludo/hivexp/BloodHIV1.0.full.json", useXref=False, mixed=True, kwds=["source","name"], result=False, grid=True, packing_options=True, indent=False, quaternion=True)
#h.saveRecipe("/opt/data/dev/cellPACK/cellPACK_data_git/cellPACK_database_1.1.0/recipes/Mycoplasma1.6_full.json", useXref=False, mixed=True, result=False, grid=True, packing_options=True, indent=False, quaternion=True)

#dataavg=numpy.loadtxt("/opt/data/dev/cellPACK/cellPACK_data_git/cellPACK_database_1.1.0/other/3ghg_ABCDEFMNOP.pdb_avg.vox")
#helper.instancesSphere("test",data,[10.0,]*len(data),helper.getObject("Sphere"),[[1,0,0]],None)