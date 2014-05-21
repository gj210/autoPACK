# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 15:27:52 2013

@author: ludovic Autin

2D integration of circle rectangle intersection area from 
http://www.eex-dev.net/index.php?id=100
Copyright © Emanuel Jöbstl, 18.06.2011
reative Commons Attribution-ShareAlike 

3D integration sphee cube intersection volume from 
http://crowsandcats.blogspot.com/2013/04/cube-sphere-intersection-volume.html
Copyright ©  Andrew White 
"""
import math
import numpy as np


try :
    from scipy.integrate import quad
except :
    quad = None
    print ("no scipy installed")
    
from math import pi
 
class Rectangle:
    def __init__(self,top,bottom, right, left):
        self.Top=top
        self.Bottom = bottom
        self.Right = right
        self.Left = left
        
class GeometriTools:  
    Resolution = 0.01

#==============================================================================
# 2D circle - rectangle intersection area
#==============================================================================

#    // rect is a System.Drawing.RectangleF structure 
#    which defines our rectangle  
#    // x is the x coordinate to calculate f(x) for  
    def  UpperRectangleFunction(self,rect, x)  :
        return rect.Top 
      
    def LowerRectangleFunction(self,rect, x) : 
        return rect.Bottom
      
#    // m is a System.Drawing.PointF which represents 
#      the center of our circle  
#    // r is the radius of our circle  
#    // x is the x coordinate to calculate f(x) for  
#   // circle is x**2+y**2=r
#      //(x − a)2 + (y − b)2 = r2
#      //top coordinate is py+r?
#      //y = cy + sqrt[r*r - (cx-4)**2]
    def UpperCircleFunction(self,m, r,  x):  
        return m[1] + math.sqrt((r * r) - pow((x - m[0]), 2))  
      
    def LowerCircleFunction(self,m, r, x) : 
        return m[1] - math.sqrt((r * r) - pow((x - m[0]), 2))  

    def GetDistance(self,p1,p2):
        return math.sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2));  

    def check_sphere_inside(self,rect,m,r):
        """check if sphere completly inside rect
        check if the four corner are more than radius away from circle center        
        """
        #1000,0,1000,0 bottom top right left
        if (abs(rect.Bottom - m[1]) > r and abs(rect.Left - m[0]) > r and
             abs(m[1] - rect.Top) > r and abs(m[0] - rect.Right) > r) :
                 #sphere inside 
                 return 0
        return 1 
        
#//Check whether the rectangle lies completely outside of the circle.   
#//Note: It is easier to check if a rectangle is outside another rectangle or  
#//circle than to check whether it is inside.  
    def check_rectangle_oustide(self,rect,m,r):
        if ((rect.Bottom < m[1] and rect.Right < m[0])  
             and (self.GetDistance([rect.Bottom, rect.Right], m) > r) or  
           (rect.Top > m[1] and rect.Right < m[0])  
             and (self.GetDistance([rect.Top, rect.Right], m) > r) or  
           (rect.Bottom < m[1] and rect.Left > m[0])  
             and (self.GetDistance([rect.Bottom, rect.Left], m) > r) or  
          (rect.Top > m[1] and rect.Left > m[0])  
             and (self.GetDistance([rect.Top, rect.Left], m) > r))  :
             return 0  # //Terminate fast 
        return 1 

    def getBoundary(self,rect,m,r):
        #//A variable storing the nearest horizontal 
        #edge of the rectangle.   
        nearestRectangleEdge = 0  
          
        #//Determine what is nearer to the circle center - 
        #  the rectangle top edge or the rectangle bottom edge  
        if (abs(m[1] - rect.Top) > abs(m[1] - rect.Bottom))  :
            nearestRectangleEdge = rect.Bottom;  
        else : 
            nearestRectangleEdge = rect.Top;   
          
        #//The bounds of our integration  
        leftBound = 0;  
        rightBound = 0;  
#        print m[1] >= rect.Bottom and m[1] <= rect.Top
        if (m[1] >= rect.Bottom and m[1] <= rect.Top):  
            #//Take care if the circle's center lies 
            #within the rectangle.   
            leftBound =  max(-r + m[0], rect.Left);  
            rightBound = min(r + m[0], rect.Right);  
#            print leftBound,rightBound
            return leftBound,rightBound
        elif (r >= abs(nearestRectangleEdge - m[1])) : 
            #//If the circle's center lies outside of the rectangle, we can choose optimal bounds.  
            leftBound = max(-math.sqrt(r * r - abs(pow(nearestRectangleEdge - m[1], 2))) + m[0], rect.Left);  
            rightBound = min(math.sqrt(r * r - abs(pow(nearestRectangleEdge - m[1], 2))) + m[0], rect.Right);  
            return leftBound,rightBound
        return leftBound,rightBound
    
    def get_rectangle_cercle_area(self,rect,m,r,leftBound,rightBound):
#        self.Resolution = 0.01#M;  
          
        #//Loop trough the intersection area and sum up 
        # the area 
        a=0.0
        i=leftBound + self.Resolution
#        print a,i,self.Resolution,leftBound,rightBound
        while (i <= rightBound):
            upperBound = min(self.UpperRectangleFunction(rect, i - self.Resolution / 2.0), 
                             self.UpperCircleFunction(m, r, i - self.Resolution / 2.0));  
            lowerBound = max(self.LowerRectangleFunction(rect, i - self.Resolution / 2.0), 
                             self.LowerCircleFunction(m, r, i - self.Resolution / 2.0));  
#            print upperBound,lowerBound
            a += (upperBound - lowerBound ) * self.Resolution;
            i+=self.Resolution
        return a;  

#==============================================================================
#    3d Sphere Cube intersection     
#==============================================================================
    #Require Scipy
    #square regions
    def region_1(self,rho, d):
        return 0.5 * (d**2 * np.sqrt(rho**2 - 2 * d ** 2))
     
    def region_1_2_theta(self,rho, d):
        return np.arccos(d / np.sqrt(rho ** 2 - d ** 2))
     
    #curved square regions
    def region_2_integrand(self,theta, rho, d):
        return np.sqrt(np.cos(theta) ** 2 - (d / rho)**2) / (np.cos(theta) ** 3)
     
    def region_2(self,rho, d):
        #require scipy
        i4 = d**3 / 6. * (rho**2 / d **2 - 1) * (pi / 4 - self.region_1_2_theta(rho, d))
        i3 = d ** 2 * rho / 3. * quad(self.region_2_integrand, self.region_1_2_theta(rho,d), pi / 4, args=(rho, d))[0]
        return i4 + i3    
     
    #spherical region
    def region_3_integrand(self,theta, rho, d):
        return np.sqrt(np.cos(theta) ** 2 - (d / rho)**2) / np.cos(theta)
     
    def region_3(self,rho, d):
        #require scipy
        return rho ** 3 / 3. * (d / rho * (pi / 4 - self.region_1_2_theta(rho, d)) - quad(self.region_3_integrand, self.region_1_2_theta(rho, d), pi / 4, args=(rho, d))[0])
     
    def calc_volume(self,rho, d):         
        alpha = rho / d         
        if(alpha <= 1):
            return 4./3 * pi * (rho) ** 3         
        if(alpha <= np.sqrt(2)):
            return 4. / 3 * pi * (rho) ** 3 - 6. * pi * (2 * (rho)**3 - 3 * d * (rho)**2 + d**3) / 3.         
        if(alpha < np.sqrt(3)):
            return 16. * (self.region_1(rho,d) + self.region_2(rho, d) + self.region_3(rho,d))        
        return 8. * d ** 3
