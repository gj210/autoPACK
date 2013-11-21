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
import os
import json
try :
    import urllib.request as urllib# , urllib.parse, urllib.error
except :
    import urllib
usePP = False
helper = None
#try :
#    from panda3d.core import Mat4
#    LISTPLACEMETHOD =["jitter","spring","rigid-body","pandaBullet","pandaBulletRelax"]
#except:
#    LISTPLACEMETHOD =  ["jitter","spring","rigid-body"]

LISTPLACEMETHOD =["jitter","spring","rigid-body","pandaBullet","pandaBulletRelax","pandaDev","RAPID"]
afdir = os.path.abspath(__path__[0])
ncpus = 2
forceFetch = False
checkAtstartup = True
verbose = 0 
messag = '''Welcome to autoPACK.
Please update to the latest version under the Help menu.
'''
autoPACKserver="http://autofill.googlecode.com/git"
replace_autoPACKserver=["autoPACKserver","http://autofill.googlecode.com/git"]
autopackdir=str(afdir)#copy
replace_autopackdir=["autopackdir",afdir]
replace_path=[replace_autoPACKserver,replace_autopackdir]

#we have to change the name of theses files. and decide how to handle the 
#currated recipeList, and the dev recipeList
#same for output and write theses file see below for the cache directories
recipe_web_pref_file = afdir+os.sep+"recipe_available.json"
recipe_user_pref_file = afdir+os.sep+"user_recipe_available.json"
recipe_dev_pref_file = afdir+os.sep+"autopack_serverDeveloper_recipeList.json"
autopack_path_pref_file = afdir+os.sep+"path_preferences.json"

if not os.path.isfile(afdir+os.sep+"version.txt"):
    f=open(afdir+os.sep+"version.txt","w")
    f.write("0.0.0")
    f.close()
f = open(afdir+os.sep+"version.txt","r")
__version__ = f.readline()
f.close()


info_dic = ["setupfile","resultfile","wrkdir"]
#change the setupfile access to online in recipe_available.xml
#change the result access to online in recipe_available.xml

#hard code recipe here is possible
global RECIPES
RECIPES = {}
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

def checkURL(URL):
    try :
        response = urllib.urlopen(URL)
    except :
        return False
    return response.code != 404
 
def retrieveFile(filename,destination=os.sep):
    filename=autopack.fixOnePath(filename)
    if filename.find("http") != -1 or filename.find("ftp")!= -1 :
        name = filename.split("/")[-1]#the recipe name
        tmpFileName = afdir+os.sep+"autoFillRecipeScripts"+os.sep+destination+name
        if not os.path.exists(afdir+os.sep+"autoFillRecipeScripts"+os.sep+destination):
		    os.makedirs(afdir+os.sep+"autoFillRecipeScripts"+os.sep+destination)
        #check if exist first
        if not os.path.isfile(tmpFileName) or forceFetch :
            if checkURL(filename):
                urllib.urlretrieve(filename, tmpFileName,reporthook=None)
            else :
                if not os.path.isfile(tmpFileName)  :
                    return  None
        filename = tmpFileName
        return filename
    return filename

def fixOnePath(p):
    for v in replace_path:
        p=p.replace(v[0],v[1])
    return p

def fixPath(adict, k, v):
    for key in list(adict.keys()):
        if key == k:
            if type(v) is list or type(v) is tuple:
                adict[key]=adict[key].replace(v[0],v[1])
            else :
                adict[key] = v
        elif type(adict[key]) is dict:
            fixPath(adict[key], k, v)

def updatePathJSON():
    if not os.path.isfile(autopack_path_pref_file):
        print (autopack_path_pref_file+" file is not found")
        return
    f=open(autopack_path_pref_file,"r")
    pref_path = json.load(f)
    f.close()
    autoPACKserver=pref_path["autoPACKserver"]
    replace_autoPACKserver[1]=autoPACKserver
    if "autopackdir" in  pref_path:
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
            
def checkPath():
    fname = autoPACKserver+"/autoPACK_filePaths.json"
    try :
        import urllib.request as urllib# , urllib.parse, urllib.error
    except :
        import urllib
    if checkURL(fname):
        urllib.urlretrieve(fname, autopack_path_pref_file)
    else :
        print ("problem accessing "+fname)
       
def checkRecipeAvailable():
#    fname = "http://mgldev.scripps.edu/projects/AF/datas/recipe_available.xml"
#    fname = "https://sites.google.com/site/autofill21/recipe_available/recipe_available.xml?attredirects=0&d=1"#revision2
    fname = autoPACKserver+"/autopack_recipe.json"
    try :
        import urllib.request as urllib# , urllib.parse, urllib.error
    except :
        import urllib
    if checkURL(fname):
        urllib.urlretrieve(fname, recipe_web_pref_file)
    else :
        print ("problem accessing "+fname)
        
def updateRecipAvailableJSON(recipesfile):
    if not os.path.isfile(recipesfile):
        print (recipesfile+" was not found")
        return
    #replace shortcut pathby hard path
    f=open(recipesfile,"r")
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
    fixPath(RECIPES,"setupfile",replace_autoPACKserver)
    fixPath(RECIPES,"wrkdir",replace_autopackdir)
    fixPath(RECIPES,"resultfile",replace_autoPACKserver)
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
        json.dump(recipe_dictionary,fp,indent=4, separators=(',', ': '))#,indent=4, separators=(',', ': ')
    
#we should read a file to fill the RECIPE Dictionary so we can add some and write/save setup 
#afdir  or user_pref
if checkAtstartup :
    checkPath()
updatePathJSON()
print ("path are updated ")

if checkAtstartup :
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
cache = wkr+os.sep+"cache_results"
if not os.path.exists(cache):
    os.makedirs(cache)
cachei = wkr+os.sep+"cache_ingredients"
if not os.path.exists(cachei):
    os.makedirs(cachei)
cache_sphere = wkr+os.sep+"cache_ingredients"+os.sep+"sphereTree"
if not os.path.exists(cache_sphere):
    os.makedirs(cache_sphere)
cacheo = wkr+os.sep+"cache_organelles"
if not os.path.exists(cacheo):
    os.makedirs(cacheo)
