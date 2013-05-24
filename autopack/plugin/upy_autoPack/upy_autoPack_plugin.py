# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 23:53:00 2012
###############################################################################
#
# autoPACK Authors: Graham T. Johnson, Ludovic Autin, Mostafa Al-Alusi, Michel Sanner
#   Based on COFFEE Script developed by Graham Johnson between 2005 and 2010
#   with assistance from Mostafa Al-Alusi in 2009 and periodic input
#   from Arthur Olson's Molecular Graphics Lab
#
# AFGui.py Authors: Ludovic Autin with minor editing/enhancement from Graham Johnson
#
# Copyright: Graham Johnson ©2010
#
# This file "upy_autoPack_plugin.py" is part of autoPACK, cellPACK.
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
    
Name: 'autopack/AutoCell GUI'
@author: Ludovic Autin with minor editing/enhancement by Graham Johnson
"""


#this should be part of the adaptor?
__author__ = "Ludovic Autin, Graham Johnson"
__url__ = [""]
__version__="0.0.0.1"
__doc__ = "AF v"+__version__
__doc__+"""\
AUTOFILL by Graham Jonhson,Ludovic Autin,Michel Sanner.
Develloped in the Molecular Graphics Laboratory directed by Arthur Olson.
Develloped @UCSF.
"""
# -------------------------------------------------------------------------- 
# ***** BEGIN GPL LICENSE BLOCK ***** 
# 
# This program is free software; you can redistribute it and/or 
# modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation; either version 2 
# of the License, or (at your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License 
# along with this program; if not, write to the Free Software Foundation, 
# Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 
# 
# ***** END GPL LICENCE BLOCK ***** 
# -------------------------------------------------------------------------- 
#=======
#should be universal
import os,sys
import platform

#SHOULD HANDLE SOME PYTHON SYSTEM PATH HERE

import sys
import os
from time import time
	
#be sure to use a unique ID obtained from www.plugincafe.com
#from Blender import Draw
PLUGIN_ID = 5555551#need a plug id

#from maya import cmds
#prefpath=cmds.internalVar(userPrefDir=True)
#os.chdir(prefpath)
#os.chdir(".."+os.sep)
#softdir = os.path.abspath(os.curdir)
#mgldirfile=prefpath+os.sep+"mgltoolsdir"
#local = False
#localpath = softdir+os.sep+"plug-ins"+os.sep+"MGLToolsPckgs"
#
###what about the path ?
##import c4d
##prefpath=c4d.storage.GeGetC4DPath(1)
##os.chdir(prefpath)
#os.chdir(".."+os.sep+"plugins"+os.sep+"autopack")
#plugdir = os.path.abspath(os.curdir)
#
#sys.path.append(plugdir)

#upy UIadaptor stuff
import upy
upy.setUIClass() #set the class
#get the pluginClass
plugTypeClass,opType = upy.getPluginClass(plug="command")#= operator in blender

from autopack import AFGui

class af_Dialog(plugTypeClass):
    plugin_name =  "autopack"
    plugin_id = PLUGIN_ID
    plugin_tooltip = "This is autopack"
    hasGui = True

    def setgui(self,dname):
        self.gui = AFGui.AFGui()
        self.gui.setup(rep=dname,host=upy.host)
        self.hasGui = True
        self.gui.display()
        print ("gui setup")
        
    def resetgui(self,dname):
        self.gui = AFGui.AFGui()
        self.gui.setup(rep=dname,host=upy.host)
        self.hasGui = True
        self.gui.display()
        
af_plugin = af_Dialog(name="autopack",pluginId=PLUGIN_ID,
                              tooltip="This is autopack",
                              hasGui=True)

af_plugin.setIcon(image_name="autoCell.tif")

if "__res__" in locals() :
    print __res__
    af_plugin.register(af_Dialog,Object=af_plugin,menuadd={"head":None,"mt":None},res=__res__)
#else :
#    af_plugin.register(af_Dialog,Object=af_plugin,menuadd={"head":None,"mt":None})

def register():
    print (__name__)
def unregister():
    pass

#Maya function
def initializePlugin(mobject):
    af_plugin.register(af_Dialog,mobject=mobject,menuadd={"head":None,"mt":None})
def uninitializePlugin(mobject):
    af_plugin.deregister(object=mobject)

if __name__ == '__main__':
    af_plugin.register(af_Dialog,Object=af_plugin,menuadd={"head":None,"mt":None})
