import sys
import os
#import c4d

sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/")
sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/PIL/")
#windows path
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17_8DE13DAD\plugins\ePMV\mgl64\MGLToolsPckgs")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17_8DE13DAD\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")

import autopack

helper = autopack.helper
if helper is None :
    import upy
    helperClass = upy.getHelperClass()
    helper = helperClass(vi="nogui")
autopack.helper = helper


from autopack.Environment import Environment
#from autopack.Graphics import AutopackViewer as AFViewer

workingdir = "C:\Dev\IntegraseProject\\"
recipefile = "HIV_IN_XP.json"
resultfile = "INT_50_random_tr.json"

workingdir = "C:\\Users\\ludov\\OneDrive\\Documents\\myRecipes\\"
recipefile = "rbc1.json"
resultfile = "rbc_result_tr.json"

workingdir = "D:\\Data\\cellPAC_data\\cellPACK_database_1.1.0\\"
recipefile ="recipes\\Mycoplasma1.6_full.json"
resultfile = "results\\Mycoplasma1.6_full_result_1.json"
		
#workingdir ="C:\Users\ludov\OneDrive\Documents\cellVIEW-i\Data\\"
#recipefile="christmas.json"

fileName, fileExtension = os.path.splitext(recipefile)
n=os.path.basename(fileName)
h = Environment(name=n)
#h.helper = helper
recipe=n
h.loadRecipe(workingdir+recipefile)

r=h.loadResult(workingdir+resultfile,transpose=False)#for c4d transpose
ingredients = h.restore(*r)

from autopack.IOutils import serializedRecipe, saveResultBinary
djson, all_pos, all_rot = serializedRecipe(h,False,True)#transpose, use_quaternion, result=False, lefthand=False
with open(workingdir+fileName+"_serialized.json","w") as f:
    f.write(djson)
saveResultBinary(h,workingdir+fileName+"_serialized.bin",False,True, lefthand=True)        
