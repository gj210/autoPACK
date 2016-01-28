import sys
import os
#import c4d

sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/")
sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/PIL/")

import autopack

helper = autopack.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper = helperClass(vi="nogui")
autopack.helper = helper


from autopack.Environment import Environment
#from autopack.Graphics import AutopackViewer as AFViewer

filename = "/home/ludo/hivexp/BloodHIV1.0.json"
path = "/home/ludo/hivexp/"
fileName, fileExtension = os.path.splitext(filename)
n=os.path.basename(fileName)
h = Environment(name=n)
#h.helper = helper
recipe=n
h.loadRecipe(filename)
h.placeMethod = "pandaBullet"
#h.resultfile = "/home/ludo/hivexp/pack_hiv_from_ncfix"
#h.helper = helper
#afviewer=None
#print h,helper
#setattr(h,"helper",helper)
#afviewer = AFViewer(ViewerType=h.helper.host,helper=h.helper)
#afviewer.SetHistoVol(h,20.0,display=False)
#h.host=h.helper.host
#previousresult = "/home/ludo/hivexp/BloodHIV1.0_mixed.json"
#previousresult = "/home/ludo/hivexp/BloodHIV1.0_centered.json"
#previousresult = "/home/ludo/hivexp/HIVCAfix.json"#good in c4d. transpose
previousresult = "/home/ludo/hivexp/pack_hiv_from_ncfix.json"
#previousresult = "/home/ludo/hivexp/pack_hiv_from_ncfix_tr.json"
##prepare a file with only the CA
## the ingredient order is not keep ?
r=h.loadResult(previousresult,transpose=False)#for c4d transpose
ingredients = h.restore(*r)

from autopack.IOutils import serializedRecipe,saveResultBinary
djson = serializedRecipe(h)
f=open("/home/ludo/hivexp/pack_hiv_from_ncfix_serialized.json","w")
f.write(djson)
f.close()
saveResultBinary(h,"/home/ludo/hivexp/pack_hiv_from_ncfix_serialized.bin",False,True)