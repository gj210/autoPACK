# upyAutoPACKPlugin
import os
import sys
import upy
from PySide import QtGui, QtCore, shiboken

found = True
MGL_ROOT=upy.__path__[0]+"\\..\\"

sys.path.append(MGL_ROOT+os.sep+"PIL")
sys.path.append(MGL_ROOT+os.sep+"lib-tk")

upy.setUIClass()

from autopack import Gui

#if 3dsMax not in the PATH add it
if "C:\\Program Files\\Autodesk\\3ds Max 2015" not in os.environ['PATH']:
    path = os.environ['PATH']
    os.environ['PATH'] = "C:\\Program Files\\Autodesk\\3ds Max 2015" + os.pathsep +  path

def upyAP_Execute():
    apgui = Gui.AutoPackGui() 
    apgui.setup(rep="autopack",host='3dsmax')
    apgui.display()
    return True
    
#app = QtGui.QApplication.instance()
#if not app:#
#	app = QtGui.QApplication([])
	
def main():		
    #MaxPlus.FileManager.Reset(True)#?
    upyAP_Execute(  )
    
if __name__ == '__main__':
	main()