# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 15:46:19 2013

@author: ludo
"""
import sys
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14 Student_7B992864/plugins/ePMV/mgl64/MGLToolsPckgs")
sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R14 Student_7B992864/plugins/ePMV/mgl64/MGLToolsPckgs/PIL/")

from AutoFill.GeometryTools import GeometriTools,Rectangle
g=GeometriTools()
rect=Rectangle(0,1000,0,1000)#top,bottom, right, left
m=[500.0,500.0]
r=710.0
ch=g.check_rectangle_oustide(rect,m,r)
if ch :
    leftBound,rightBound = g.getBoundary(rect,m,r)
    print leftBound,rightBound
    area = g.get_rectangle_cercle_area(rect,m,r,rightBound,leftBound)
    print (area)
