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
try :
    import urllib.request as urllib# , urllib.parse, urllib.error
except :
    import urllib

helper = None
autoPACKserver="http://autofill.googlecode.com/svn/data/"
#try :
#    from panda3d.core import Mat4
#    LISTPLACEMETHOD =["jitter","spring","rigid-body","pandaBullet","pandaBulletRelax"]
#except:
#    LISTPLACEMETHOD =  ["jitter","spring","rigid-body"]

LISTPLACEMETHOD =["jitter","spring","rigid-body","pandaBullet","pandaBulletRelax"]
afdir = os.path.abspath(__path__[0])

forceFetch = False
checkAtstartup = True
messag = '''Welcome to autoPACK.
Please update to the latest version under the Help menu.
'''

recipe_web_pref_file = afdir+os.sep+"recipe_available.xml"
recipe_user_pref_file = afdir+os.sep+"user_recipe_available.xml"

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
RECIPES = {
#"Test_CylindersSpheres2D":{
#    "1.0":
#    {
#    "setupfile":afdir+os.sep+"autoFillRecipeScripts"+os.sep+"2DcylinderSphereFill"+os.sep+"2DCylindersSpheres_setup_recipe.py",
#    "resultfile":afdir+os.sep+"autoFillRecipeScripts"+os.sep+"2DcylinderSphereFill"+os.sep+"results"+os.sep+"CylSpherefillResult.afr.txt",
#    "wrkdir":afdir+os.sep+"autoFillRecipeScripts"+os.sep+"2DcylinderSphereFill"
#    }
#}
}

USER_RECIPES={}

def checkURL(URL):
    try :
        response = urllib.urlopen(URL)
    except :
        return False
    return response.code != 404
    
def checkRecipeAvailable():
    fname = "http://autofill.googlecode.com/svn/data/recipe_available.xml"
    try :
        import urllib.request as urllib# , urllib.parse, urllib.error
    except :
        import urllib
    if checkURL(fname):
        urllib.urlretrieve(fname, recipe_web_pref_file)
    
def updateRecipAvailable(recipesfile):
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
    
#we should read a file to fill the RECIPE Dictionary so we can add some and write/save setup 
#afdir  or user_pref
if checkAtstartup :
    checkRecipeAvailable()
updateRecipAvailable(recipe_web_pref_file)
updateRecipAvailable(recipe_user_pref_file)

#check cach directory create if doesnt exit.abs//should be in user pref?
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