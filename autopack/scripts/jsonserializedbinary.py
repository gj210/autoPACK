import sys
import os
#import c4d

sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/")
sys.path.append("/home/ludo/Tools/mgltools_x86_64Linux2_latest/MGLToolsPckgs/PIL/")
#windows path
#sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17_8DE13DAD\plugins\ePMV\mgl64\MGLToolsPckgs")
#sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17_8DE13DAD\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")

sys.path.insert(0,"C:\Users\ludov\AppData\Roaming\MAXON\Cinema 4D R19 Demo_B31EE4CE\plugins\ePMV\mgl64\MGLToolsPckgs")
sys.path.insert(1,"C:\Users\ludov\AppData\Roaming\MAXON\Cinema 4D R19 Demo_B31EE4CE\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")

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

workingdir = "D:\\Data\\cellPAC_data\\cellPACK_database_1.1.0\\"
recipefile ="recipes\\Mycoplasma1.6_full.json"
resultfile = "results\\Mycoplasma1.6_full_result_1.json"
		
  
recipefile ="recipes\\BloodHIVMycoRB.1.0.json"
resultfile = ""

recipefile ="recipes\\Mycoplasma1.7.json"
resultfile = ""

recipefile ="recipes\\DNAplectoneme.1.0.json"

recipefile ="recipes\\DNAplectoneme.1.0.json"

recipefile ="D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_1.0.json"
tr=""#"_tr"
transpose = True
resultfile ="D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_Results_H"+tr+".json"# or Tr
workingdir ="D:\\Data\\cellPACK_data\\Mycoplasma\\"

workingdir = "D:\\Data\\cellPACK_data_git\\influenza\Influenza\\"
recipefile = "Influenza_3.json"
resultfile = ""

#workingdir ="C:\Users\ludov\OneDrive\Documents\cellVIEW-i\Data\\"
#recipefile="christmas.json"

fileName, fileExtension = os.path.splitext(recipefile)
n=os.path.basename(fileName)
h = Environment(name=n)
#h.helper = helper
recipe=n
h.loadRecipe(workingdir+recipefile)

# save a full recipe no Xref
#"nbMol": nbmol,  
#"molarity" : ingr.molarity,
#"source": ingr.source, 
#"positions":ingr.positions, 
#"radii":ingr.radii}

rname=workingdir+"recipes\\BloodHIVMycoRB_full.1.0.json"
rname=workingdir+"recipes\\Mycoplasma1.7.json"
rname=workingdir+recipefile
#h.saveRecipe(rname,useXref=False,mixed=True,
#                     kwds=["source","name","nbMol","molarity","positions","radii"],result=True,
#                   grid=False,packing_options=False,indent=False,quaternion=True)  
if resultfile!= "" :
    r=h.loadResult(resultfile,transpose=transpose)#for c4d transpose
    ingredients = h.restore(*r)
    h.saveRecipe("D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_Results_H_1.json",useXref=False,mixed=True, kwds=["source","name"],result=True, transpose = True,  grid=False,packing_options=False,indent=False,quaternion=True) 
    h.saveRecipe("D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_Results_H_2.json",useXref=False,mixed=True, kwds=["source","name"],result=True, transpose = False, grid=False,packing_options=False,indent=False,quaternion=True)                            

#                           
##gagpol [0,0,-143.687] offset, pcpal 0,0,1
from autopack.IOutils import serializedRecipe, saveResultBinary, toBinary
djson, all_pos, all_rot = serializedRecipe(h,False,True,True,True)#transpose, use_quaternion, result=False, lefthand=False
with open(workingdir+fileName+"_serialized.json","w") as f:
    f.write(djson)
if resultfile!= "" :
    toBinary(all_pos, all_rot,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_1.0_serialized_L01"+tr+".bin")     
#    #saveResultBinary(h,fileName+"_serialized.bin",False,True, lefthand=True)        
#from autopack.Serializable import sCompartment,sIngredientGroup,sIngredient,sIngredientFiber
#sCompartment.static_id = 0;sIngredientFiber.static_id = 0;sIngredient.static_id = [0, 0, 0];sIngredientGroup.static_id= 0;djson, all_pos, all_rot = serializedRecipe(h,True,True,True,True);toBinary(all_pos, all_rot,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_1.0_serializedTR_L"+tr+".bin")     
#sCompartment.static_id = 0;sIngredientFiber.static_id = 0;sIngredient.static_id = [0, 0, 0];sIngredientGroup.static_id= 0;djson, all_pos, all_rot = serializedRecipe(h,False,True,True,True);toBinary(all_pos, all_rot,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_1.0_serialized_L"+tr+".bin")     
#sCompartment.static_id = 0;sIngredientFiber.static_id = 0;sIngredient.static_id = [0, 0, 0];sIngredientGroup.static_id= 0;djson, all_pos, all_rot = serializedRecipe(h,False,True,True,False);toBinary(all_pos, all_rot,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_1.0_serialized"+tr+".bin")     
#sCompartment.static_id = 0;sIngredientFiber.static_id = 0;sIngredient.static_id = [0, 0, 0];sIngredientGroup.static_id= 0;djson, all_pos, all_rot = serializedRecipe(h,True,True,True,False);toBinary(all_pos, all_rot,"C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_1.0_serializedTR"+tr+".bin")     
#
#
#with open("C:\\Users\\ludov\\OneDrive\\Documents\\OnlinePacking_Tobias\\cellVIEW-OP\\Data\\Mpn_1.0_serialized.json","w") as f:
#    f.write(djson)