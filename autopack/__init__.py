# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 23:53:00 2012
###############################################################################
#
# autoPACK Authors: Graham T. Johnson, Mostafa Al-Alusi, Ludovic Autin, Michel Sanner
#   Based on COFFEE Script developed by Graham Johnson between 2005 and 2010
#   with assistance from Mostafa Al-Alusi in 2009 and periodic input
#   from Arthur Olson's Molecular Graphics Lab
#
# __init__.py Authors: Ludovic Autin with minor editing/enhancement from Graham Johnson
#
# Copyright: Graham Johnson ©2010
#
# This file "__init__.py" is part of autoPACK.
#
#    autoPACK is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    autoPACK is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with autoPACK (See "CopyingGNUGPL" in the installation.
#    If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
Name: 'autoPACK'
Define here some usefull variable and setup filename path that facilitate
AF 
@author: Ludovic Autin with editing by Graham Johnson
"""
packageContainsVFCommands = 1
import sys
import os
import re
import shutil
from os import path, environ
use_json_hook=True
#in casesimplejson not there
try :
    import simplejson as json
    #we can use the hookup at loading
except :
    import json
    if sys.version[0:3] <= "2.6" :
        use_json_hook=False
try :
    from collections import OrderedDict
except :
    from ordereddict import OrderedDict
try :
    import urllib.request as urllib# , urllib.parse, urllib.error
except :
    import urllib

afdir = os.path.abspath(__path__[0])

#==============================================================================
# #Setup autopack data directory.
#==============================================================================
#the dir will have all the recipe + cache.
APPNAME = "autoPACK"
if sys.platform == 'darwin':
    #from AppKit import NSSearchPathForDirectoriesInDomains
    # http://developer.apple.com/DOCUMENTATION/Cocoa/Reference/Foundation/Miscellaneous/Foundation_Functions/Reference/reference.html#//apple_ref/c/func/NSSearchPathForDirectoriesInDomains
    # NSApplicationSupportDirectory = 14
    # NSUserDomainMask = 1
    # True for expanding the tilde into a fully qualified path
    #appdata = path.join(NSSearchPathForDirectoriesInDomains(14, 1, True)[0], APPNAME)
    appdata = os.path.expanduser("~")+"/Library/Application Support/autoPACK"
elif sys.platform == 'win32':
    appdata = path.join(environ['APPDATA'], APPNAME)
else:
    appdata = path.expanduser(path.join("~", "." + APPNAME))
if not os.path.exists(appdata):
    os.makedirs(appdata)
    print("autoPACK data dir created")
    print(appdata)

#==============================================================================
# setup Panda directory
#==============================================================================
PANDA_PATH=""
if sys.platform == 'darwin':
    PANDA_PATH=afdir+os.sep+".."+os.sep+"Panda3D"     
    sys.path.append("/Developer/Panda3D/lib/")#in case already installed
    #TODO need to fix the dependancy that are locally set to /Developer/Panda3D/lib/
elif sys.platform == 'win32':
    PANDA_PATH=afdir+os.sep+".."+os.sep+"Panda3d-1.9.0-x64"
    PANDA_PATH_BIN=PANDA_PATH+os.sep+"bin"
    try:
        if PANDA_PATH_BIN not in os.environ.get('PATH', ''):
            os.environ['PATH'] = os.pathsep.join((PANDA_PATH_BIN, os.environ.get('PATH', '')))
    except Exception:
        pass
    sys.path.append(PANDA_PATH_BIN)
    sys.path.append(PANDA_PATH)
elif sys.platform == "linux2" :#linux ? blender and maya ?
    PANDA_PATH="/usr/lib/python2.7/dist-packages/"
    PANDA_PATH_BIN="/usr/lib/panda3d/"
    #sys.path.append(PANDA_PATH)
else :
    pass
sys.path.append(PANDA_PATH+os.sep+"lib")

def checkURL(URL):
    try :
        response = urllib.urlopen(URL)
    except :
        return False
    return response.code != 404



#==============================================================================
# setup the cache directory inside the app data folder
#==============================================================================
cache_results = appdata+os.sep+"cache_results"
if not os.path.exists(cache_results):
    os.makedirs(cache_results)
cache_geoms = appdata+os.sep+"cache_geometries"
if not os.path.exists(cache_geoms):
    os.makedirs(cache_geoms)
cache_sphere = appdata+os.sep+"cache_collisionTrees"
if not os.path.exists(cache_sphere):
    os.makedirs(cache_sphere)
#cacheo = appdata+os.sep+"cache_organelles"
#if not os.path.exists(cacheo):
#    os.makedirs(cacheo)
cache_recipes=appdata+os.sep+"cache_recipes"
if not os.path.exists(cache_recipes):
    os.makedirs(cache_recipes)
    
preferences = appdata+os.sep+"preferences"
if not os.path.exists(preferences):
    os.makedirs(preferences)
#we can now use some json/xml file for storing preferences and options.
#need others ?
cache_dir={
"geometries":cache_geoms,
"results":cache_results,
"collisionTrees":cache_sphere,
"recipes":cache_recipes,
"prefs":preferences,
}

#    
#autopack_cache_data (e.g. recipe_available.json)
#or call this autopack_cache_recipelists
#autopack_cache_recipes
#autopack_cache_results
#autopack_cache_analysis
#autopack_cache_geometries
#autopack_cache_paths

usePP = False
helper = None
#try :
#    from panda3d.core import Mat4
#    LISTPLACEMETHOD =["jitter","spring","rigid-body","pandaBullet","pandaBulletRelax"]
#except:
#    LISTPLACEMETHOD =  ["jitter","spring","rigid-body"]

LISTPLACEMETHOD =["jitter","spring","rigid-body","pandaBullet","pandaBulletRelax","pandaDev","RAPID"]

ncpus = 2
#forceFetch is for any file not only recipe/ingredient etc...
forceFetch = False
checkAtstartup = True
testPeriodicity = False
biasedPeriodicity = None#[1,1,1]
fixpath = False
verbose = 0 
messag = '''Welcome to autoPACK.
Please update to the latest version under the Help menu.
'''

#we have to change the name of theses files. and decide how to handle the 
#currated recipeList, and the dev recipeList
#same for output and write theses file see below for the cache directories
#all theses file will go in the pref folder ie cache_path
recipe_web_pref_file = preferences+os.sep+"recipe_available.json"
recipe_user_pref_file = preferences+os.sep+"user_recipe_available.json"
recipe_dev_pref_file = preferences+os.sep+"autopack_serverDeveloper_recipeList.json"
autopack_path_pref_file = preferences+os.sep+"path_preferences.json"
autopack_user_path_pref_file = preferences+os.sep+"path_user_preferences.json"

#Default values    
legacy_autoPACKserver="https://autofill.googlecode.com/git/autoPACK_database_1.0.0"#XML
autoPACKserver="https://raw.githubusercontent.com/mesoscope/cellPACK_data/master/cellPACK_database_1.1.0"
autoPACKserver_default="https://raw.githubusercontent.com/mesoscope/cellPACK_data/master/cellPACK_database_1.1.0"#XML
autoPACKserver_alt="http://mgldev.scripps.edu/projects/autoPACK/data/cellPACK_data/cellPACK_database_1.1.0"
filespath = "https://raw.githubusercontent.com/mesoscope/cellPACK_data/master/autoPACK_filePaths.json"
recipeslistes = autoPACKserver+"/autopack_recipe.json"

autopackdir=str(afdir)#copy
def checkPath(autopack_path_pref_file):
    fname = filespath#autoPACKserver+"/autoPACK_filePaths.json"
    if fname.find("http") != -1 or fname.find("ftp")!= -1 :
        try :
            import urllib.request as urllib# , urllib.parse, urllib.error
        except :
            import urllib
        if checkURL(fname):
            urllib.urlretrieve(fname, autopack_path_pref_file)
        else :
            print ("problem accessing path "+fname)
    else :
        autopack_path_pref_file = fname
        
#get user / default value 
if not os.path.isfile(autopack_path_pref_file):
    print (autopack_path_pref_file+" file is not found")
    checkPath(autopack_path_pref_file)

doit=False
if os.path.isfile(autopack_user_path_pref_file):
    f=open(autopack_user_path_pref_file,"r")  
    doit=True
elif os.path.isfile(autopack_path_pref_file):
    f=open(autopack_path_pref_file,"r")
    doit=True
if doit :
    pref_path = json.load(f)
    f.close()
    if "autoPACKserver" not in pref_path :
        print ("problem with autopack_path_pref_file ",autopack_path_pref_file)
        print ("reset to default")
    else :
        autoPACKserver=pref_path["autoPACKserver"]
        if "filespath" in pref_path:
            if pref_path["filespath"] != "default" :
                filespath =pref_path["filespath"]
        if "recipeslistes" in pref_path:
            if pref_path["recipeslistes"] != "default" :
                recipeslistes =pref_path["recipeslistes"]
        if "autopackdir" in pref_path:
            if pref_path["autopackdir"] != "default" :
                autopackdir=pref_path["autopackdir"]

replace_autoPACKserver=["autoPACKserver",autoPACKserver]
replace_autopackdir=["autopackdir",autopackdir]
replace_autopackdata=["autopackdata",appdata]

replace_path=[replace_autoPACKserver,replace_autopackdir,replace_autopackdata]
global current_recipe_path
current_recipe_path=appdata
#we keep the file here, it come with the distribution 
#wonder if the cache shouldn use the version like other appDAta
#ie appData/AppName/Version/etc...
if not os.path.isfile(afdir+os.sep+"version.txt"):
    f=open(afdir+os.sep+"version.txt","w")
    f.write("0.0.0")
    f.close()
f = open(afdir+os.sep+"version.txt","r")
__version__ = f.readline()
f.close()

#should we check filespath

info_dic = ["setupfile","resultfile","wrkdir"]
#change the setupfile access to online in recipe_available.xml
#change the result access to online in recipe_available.xml

#hard code recipe here is possible
global RECIPES
RECIPES = OrderedDict()
# = {
#"Test_CylindersSpheres2D":{
#    "1.0":
#    {
#    "setupfile":afdir+os.sep+"autoFillRecipeScripts"+os.sep+"2DcylinderSphereFill"+os.sep+"2DCylindersSpheres_setup_recipe.py",
#    "resultfile":afdir+os.sep+"autoFillRecipeScripts"+os.sep+"2DcylinderSphereFill"+os.sep+"results"+os.sep+"CylSpherefillResult.afr.txt",
#    "wrkdir":afdir+os.sep+"autoFillRecipeScripts"+os.sep+"2DcylinderSphereFill"
#    }
#}


USER_RECIPES={}

def resetDefault():
    if os.path.isfile(autopack_user_path_pref_file):
        os.remove(autopack_user_path_pref_file)
    autoPACKserver=autoPACKserver_default#"http://autofill.googlecode.com/git/autoPACK_database_1.0.0"
    filespath = "https://raw.githubusercontent.com/mesoscope/cellPACK_data/master/autoPACK_filePaths.json"
    recipeslistes = autoPACKserver+"/autopack_recipe.json"
    

def revertOnePath(p):
    for v in replace_path:
        p=p.replace(v[1],v[0])
    return p

def checkErrorInPath(p,toreplace):
    #if in p we already have part of the replace path
    part = p.split(os.sep)
    newpath=""
    for i,e in enumerate(part) :
        f=re.findall('{0}'.format(re.escape(e)), toreplace)
        if not len(f):
            newpath+=e+"/"
    if part[0] == "http:":
        newpath="http://"+newpath[6:]
    return newpath[:-1]

def fixOnePath(p):
    for v in replace_path:
        #fix before
        if fixpath and re.findall('{0}'.format(re.escape(v[0])), p):
            p=checkErrorInPath(p,v[1])
            #check for legacyServerautoPACK_database_1.0.0
            p=checkErrorInPath(p,"autoPACK_database_1.0.0")           
        p=p.replace(v[0],v[1])
    return p

def updateReplacePath(newPaths):
    for w in newPaths :    
        found=False
        for i,v in enumerate(replace_path):
            if v[0]==w[0]:
                replace_path[i][1]=w[1]
                found=True
        if not found:
            replace_path.append(w)
                
def retrieveFile(filename,destination="",cache="geometries",force=None):
#    helper = autopack.helper
    if force is None :
        force = forceFetch
    if filename.find("http") == -1 and filename.find("ftp")== -1 :
        filename=fixOnePath(filename)
    print ("autopack retrieve file ",filename)
    if filename.find("http") != -1 or filename.find("ftp")!= -1 :
        #check if usin autoPACKserver
        useAPServer = False
        if filename.find(autoPACKserver) != -1 :
            useAPServer = True
        reporthook = None
        if helper is not None:        
            reporthook=helper.reporthook
        name = filename.split("/")[-1]#the recipe name
        tmpFileName = cache_dir[cache]+os.sep+destination+name
        if not os.path.exists(cache_dir[cache]+os.sep+destination):
            os.makedirs(cache_dir[cache]+os.sep+destination)
        #check if exist first
        #print ("isfile ",tmpFileName)
        if not os.path.isfile(tmpFileName) or force :
            if checkURL(filename):
                try:
                    urllib.urlretrieve(filename, tmpFileName,reporthook=reporthook)
                except:
                    if useAPServer :
                        print ("try alternate server")
                        urllib.urlretrieve(autoPACKserver_alt+"/"+cache+"/"+name, tmpFileName,reporthook=reporthook)                        
            else :
                if not os.path.isfile(tmpFileName)  :
                    print ("not isfile ",tmpFileName)
                    return  None
        filename = tmpFileName
        print ("autopack return grabbed ",filename)
        #check the file is not an error            
        return filename
    print ("autopack search ",filename)
    #if no folder provided, use the current_recipe_folder
    if not os.path.isfile(filename):#use the cache system ? geom / recipes / ingredients / others ?
        print ("search in cache",cache_dir[cache]+os.sep+filename)
        print ("search on server",autoPACKserver+"/"+cache+"/"+filename)
        print ("search on alternate_backup_server",autoPACKserver_alt+"/"+cache+"/"+filename)
        print ("search in curr_dir",current_recipe_path+os.sep+filename)
        if  os.path.isfile(cache_dir[cache]+os.sep+filename):
            return cache_dir[cache]+os.sep+filename
        elif checkURL(autoPACKserver+"/"+cache+"/"+filename):
            reporthook = None
            if helper is not None:        
                reporthook=helper.reporthook
            name = filename.split("/")[-1]#the recipe name
            tmpFileName = cache_dir[cache]+os.sep+destination+name
            try:
                urllib.urlretrieve(autoPACKserver+"/"+cache+"/"+filename, tmpFileName,reporthook=reporthook)
                return tmpFileName
            except:
                print ("try alternate server")
                urllib.urlretrieve(autoPACKserver_alt+"/"+cache+"/"+filename, tmpFileName,reporthook=reporthook)
                #check the file is not an error
                return tmpFileName
        elif checkURL(autoPACKserver_alt+"/"+cache+"/"+filename):
            reporthook = None
            if helper is not None:        
                reporthook=helper.reporthook
            name = filename.split("/")[-1]#the recipe name
            tmpFileName = cache_dir[cache]+os.sep+destination+name
            try:
                urllib.urlretrieve(autoPACKserver_alt+"/"+cache+"/"+filename, tmpFileName,reporthook=reporthook)
            except:
                print ("not on alternate server")
                return None
            #check the file is not an error
            return tmpFileName
        elif os.path.isfile(current_recipe_path+os.sep+filename):
            return current_recipe_path+os.sep+filename
        else :
            print (filename," not found ")
    return filename

def fixPath(adict):#, k, v):
    for key in list(adict.keys()):
        if type(adict[key]) is dict or type(adict[key]) is OrderedDict:
            fixPath(adict[key])
        else :
#        if key == k:
            adict[key]=fixOnePath(adict[key])
#            if type(v) is list or type(v) is tuple:
#                adict[key]=adict[key].replace(v[0],v[1])
#            else :
#                adict[key] = v
#        elif type(adict[key]) is dict:
#            fixPath(adict[key], k, v)

def updatePathJSON():
    if not os.path.isfile(autopack_path_pref_file):
        print (autopack_path_pref_file+" file is not found")
        return
    if os.path.isfile(autopack_user_path_pref_file):
        f=open(autopack_user_path_pref_file,"r")    
    else :
        f=open(autopack_path_pref_file,"r")
    pref_path = json.load(f)
    f.close()
    autoPACKserver=pref_path["autoPACKserver"]
    replace_autoPACKserver[1]=autoPACKserver
    filespath = autoPACKserver+"/autoPACK_filePaths.json"
    if "filespath" in pref_path:
        if pref_path["filespath"] != "default" :
            filespath =pref_path["filespath"]
    recipeslistes = autoPACKserver+"/autopack_recipe.json"
    if "recipeslistes" in pref_path:
        if pref_path["recipeslistes"] != "default" :
            recipeslistes =pref_path["recipeslistes"]
    if "autopackdir" in pref_path:
        if pref_path["autopackdir"] != "default" :
            autopackdir=pref_path["autopackdir"]
            replace_autopackdir[1]=pref_path["autopackdir"]
            
def updatePath():
    #now get it
    fileName, fileExtension = os.path.splitext(autopack_path_pref_file)
    if fileExtension.lower() == ".xml":
        pass#updateRecipAvailableXML(recipesfile)
    elif fileExtension.lower() == ".json":
        updatePathJSON()
        
def checkRecipeAvailable():
#    fname = "http://mgldev.scripps.edu/projects/AF/datas/recipe_available.xml"
#    fname = "https://sites.google.com/site/autofill21/recipe_available/recipe_available.xml?attredirects=0&d=1"#revision2
    fname = fixOnePath(recipeslistes)#autoPACKserver+"/autopack_recipe.json"
    try :
        import urllib.request as urllib# , urllib.parse, urllib.error
    except :
        import urllib
    if checkURL(fname):
        urllib.urlretrieve(fname, recipe_web_pref_file)
    else :
        print ("problem accessing recipe "+fname)
        
def updateRecipAvailableJSON(recipesfile):
    if not os.path.isfile(recipesfile):
        print (recipesfile+" was not found")
        return
    #replace shortcut pathby hard path
    f=open(recipesfile,"r")
    if use_json_hook:
        recipes=json.load(f,object_pairs_hook=OrderedDict) 
    else :
        recipes=json.load(f)       
    f.close()
    RECIPES.update(recipes)
    print ("recipes updated "+str(len(RECIPES)))
    
def updateRecipAvailableXML(recipesfile):
    if not os.path.isfile(recipesfile):
        return
    from xml.dom.minidom import parse
    XML = parse(recipesfile) # parse an XML file by name
    res = XML.getElementsByTagName("recipe")
    for r in res :
        name = r.getAttribute("name")
        version = r.getAttribute("version")
        if name in RECIPES : #update te value
            if version in RECIPES[name]:
                for info in info_dic :
                    text = r.getElementsByTagName(info)[0].childNodes[0].data.strip().replace("\t","")
                    if text[0] != "/" and text.find("http") == -1:
                        text = afdir+os.sep+text
                    RECIPES[name][version][info] = str(text)
            else :
                RECIPES[name][version] = {}
                for info in info_dic :
                    text = r.getElementsByTagName(info)[0].childNodes[0].data.strip().replace("\t","")
                    if text[0] != "/" and text.find("http") == -1:
                        text = afdir+os.sep+text
                    RECIPES[name][version][info] = str(text)
        else : #append to the dictionary
            RECIPES[name]={}
            RECIPES[name][version] = {}
            for info in info_dic :
                text = r.getElementsByTagName(info)[0].childNodes[0].data.strip().replace("\t","")
                if text[0] != "/" and text.find("http") == -1:
                    text = afdir+os.sep+text
                RECIPES[name][version][info] = str(text)
    print ("recipes updated "+str(len(RECIPES)))
    
def updateRecipAvailable(recipesfile):
    if not os.path.isfile(recipesfile):
        return
    #check format xml or json
    fileName, fileExtension = os.path.splitext(recipesfile)
    if fileExtension.lower() == ".xml":
        updateRecipAvailableXML(recipesfile)
    elif fileExtension.lower() == ".json":
        updateRecipAvailableJSON(recipesfile)
    fixPath(RECIPES)
#    fixPath(RECIPES,"wrkdir")#or autopackdata
#    fixPath(RECIPES,"resultfile")
    print ("recipes updated and path fixed "+str(len(RECIPES)))
    
def saveRecipeAvailable(recipe_dictionary,recipefile):
    from xml.dom.minidom import getDOMImplementation
    impl = getDOMImplementation()
    XML = impl.createDocument(None, "autoPACK_recipe", None)
    root = XML.documentElement
    for k in recipe_dictionary:     
        for v in recipe_dictionary[k]:
            relem=XML.createElement("recipe")
            relem.setAttribute("name",k)
            relem.setAttribute("version",v)
            root.appendChild(relem)
            for l in recipe_dictionary[k][v]:
                node = XML.createElement(l)
                data = XML.createTextNode(recipe_dictionary[k][v][l])
                node.appendChild(data)
                relem.appendChild(node)
    f = open(recipefile,"w")        
    XML.writexml(f, indent="\t", addindent="", newl="\n")
    f.close()

def saveRecipeAvailableJSON(recipe_dictionary,filename):
    with open(filename, 'w') as fp :#doesnt work with symbol link ?
        json.dump(recipe_dictionary,fp,indent=1, separators=(',', ': '))#,indent=4, separators=(',', ': ')

def clearCaches(*args):
    #can't work if file are open!
    for k in cache_dir:
        try :
            shutil.rmtree(cache_dir[k])
            os.makedirs(cache_dir[k])
        except:
            print ("problem cleaning ",cache_dir[k])
   
#we should read a file to fill the RECIPE Dictionary so we can add some and write/save setup 
#afdir  or user_pref

if checkAtstartup :
    checkPath(autopack_path_pref_file)
updatePathJSON()
print ("path are updated ")

if checkAtstartup :
    #get from server the list of recipe
    #recipe_web_pref_file
    checkRecipeAvailable()
updateRecipAvailable(recipe_web_pref_file)
updateRecipAvailable(recipe_user_pref_file)
updateRecipAvailable(recipe_dev_pref_file)

print ("currently nb recipes is "+str(len(RECIPES)))
#check cach directory create if doesnt exit.abs//should be in user pref?
#?
#need a distinction between autopackdir and cachdir
wkr = afdir
#in the preefined working directory