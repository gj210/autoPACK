# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 19:00:44 2014

@author: ludo
"""
import c4d

h=c4d.af.values()[0].histoVol.values()[0]
frame=0

def main():
    #depending time -> frane
    #first frame everytinh
    frame+=1
    h.traj.applyState(h,frame)
    