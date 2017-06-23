# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 13:44:05 2016

@author: ludov
"""
import numpy as np
import json
import sys
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs\PIL")
sys.path.append("C:\Users\ludov\AppData\Roaming\MAXON\CINEMA 4D R17 Demo_E0A949BC\plugins\ePMV\mgl64\MGLToolsPckgs\lib-tk")

mainPATH = "C:\\Users\\ludov\\Downloads\\"
mainPATH = "C:\\Dev\\Brett\\"
#fin = "D:\\Data\\HIV\\blood_hiv_immature_split.json"
#fin = "C:\\Users\\ludov\\Downloads\\MG_results_tr.json"
fin = mainPATH+"RECIPE-syn1.0_complexes_wipdb_wholeRIB_tip_shifted_res_tr.json"
fout = mainPATH+"RECIPE-syn1.0_complexes_wipdb_wholeRIB_tip_shifted_res_tr_S"
fin = "C:\\Users\\ludov\\OneDrive\\Documents\\cellVIEW\\Data\\packing_results\\HIV-1_0.1.6-8_mixed_pdb.cpr"
fout = "C:\\Users\\ludov\\OneDrive\\Documents\\cellVIEW\\Data\\packing_results\\HIVcapsid.bin"
with open( fin ,"r") as f :
    result = json.load(f)

#fout = "C:\\Users\\ludov\\Downloads\\MG_1.0_serialized"

from autopack.IOutils import serializedFromResult, toBinary, gatherResult
d,p,r = serializedFromResult(result,False,True,result=True,lefthand=False)

ap, ar = gatherResult(posrot, True, True,  lefthand=True)

#toBinary(p, r,fout+".bin")
#djson,p,r = serializedFromResult(result,True,True,result=True, lefthand=False)
#toBinary(p, r,fout+"_tr.bin")
#d,p,r = serializedFromResult(result,True,True,result=True, lefthand=True)
#toBinary(p, r,fout+"_tr_lh.bin")
#d,p,r = serializedFromResult(result,False,True,result=True, lefthand=True)
#toBinary(p, r,fout+"_lh.bin")
#
#with open( fout+".json","w") as f :
#    f.write(djson)
#build color palette ?
#color_palette={}
##go through everything and generate color
#def getColor(dic,parentname=""):
#    parentname+=dic["name"]
##    print parentname
#    if len(dic["IngredientGroups"]):
#        for group in dic["IngredientGroups"]:
##            print group["name"]
#            for ingr in group["Ingredients"]:
#                print parentname+"."+ingr["name"]
#                c = np.random.random(3)*255.0
#                color_palette[parentname+"."+ingr["name"]]={"x":c[0],"y":c[1],"z":c[2]}
#    for comp in dic["Compartments"]:
#        getColor(comp,parentname=parentname+".")
#    
#new_recipe = json.loads(djson)
#
#with open( "D:\\Data\\HIV\\blood_hiv_immature_palette.json","w") as f :
#    json.dump(color_palette,f)
#
#p,r=saveResultBinaryDic(result,"D:\\Data\\HIV\\blood_hiv_immature_split_serialized_tr",True,True,False)
#p,r=saveResultBinaryDic(result,"D:\\Data\\HIV\\blood_hiv_immature_split_serialized_tr_lh",True,True,True)
#p,r=saveResultBinaryDic(result,"D:\\Data\\HIV\\blood_hiv_immature_split_serialized",False,True,False)
#p,r=saveResultBinaryDic(result,"D:\\Data\\HIV\\blood_hiv_immature_split_serialized_lh",False,True,True)
