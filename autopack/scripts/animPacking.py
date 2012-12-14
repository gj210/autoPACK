# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 21:39:41 2012

@author: -
"""

import c4d
#Welcome to the world of Python


def main():
    f="/Users/ludo/DEV/Autofill_svn_test/autofill/trunk/AutoFillClean/cache/test.afr"#path t result file
    afui = c4d.af.values()[0]
    h = afui.histoVol.values()[0]
    afui.helper.doc = doc
    if not hasattr(h,"ingrnumber"):
        h.ingrnumber={}
    if not len(h.organelles[0].molecules):
        result,orgaresult,freePoint=h.load(resultfilename=f,restore_grid=False) 
        h.organelles[0].molecules = orgaresult[0]
    orga = h.organelles[0]
    r = h.organelles[0].molecules        
    ob = op.GetObject()
    cframe = doc.GetTime().GetFrame(doc.GetFps())
    #cframe gave which ingredient to toggle
    if cframe == 0 :
        for e in r :
            pos,rot,name,compNum,ptInd = e
            if name not in h.ingrnumber :
                h.ingrnumber[name]=0
            instancename = "%sF%i%s%i" % (orga.name,h.cFill,name,h.ingrnumber[name])
            print instancename
            o=afui.helper.getObject(instancename)
            afui.helper.toggleDisplay(o,False)
            h.ingrnumber[name]+=1
        h.ingrnumber={}
    else :
        pos,rot,name,compNum,ptInd = r[cframe-1]
        if name not in h.ingrnumber :
            h.ingrnumber[name]=0
        instancename = "%sF%i%s%i" % (orga.name,h.cFill,name,h.ingrnumber[name])
        o=afui.helper.getObject(instancename)
        afui.helper.toggleDisplay(o,True)        
        h.ingrnumber[name]+=1
        