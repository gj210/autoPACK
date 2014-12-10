# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 09:04:10 2013

@author: Ludovic Autin
"""
import os
import numpy
import sys
import autopack
from autopack.Ingredient import GrowIngrediant,ActinIngrediant,KWDS


from xml.dom.minidom import getDOMImplementation 
try :
    import simplejson as json
    from simplejson import encoder    
except :
    import json
    from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.8g')
try :
    from collections import OrderedDict
except :
    from ordereddict import OrderedDict

#if python 2.6 need to convert keyword from unicode to string
def flatten_unicode_keys(d):
    if sys.version[0:3] > "2.6" :
        return d
    for k in d:
        if isinstance(k, unicode):
            v = d[k]
            del d[k]
            d[str(k)] = v
    return d
            
def getValueToXMLNode(vtype,node,attrname):
#        print "getValueToXMLNode ",attrname
    value = node.getAttribute(attrname)
#        print "value " , value
    value = str(value)
    if not len(value):
        return None
    if vtype not in ["liste","filename","string"] :
        value=eval(value)
    else :
        value=str(value)
    return value
           
def setValueToXMLNode(value,node,attrname):
    if value is None:
        print (attrname, " is None !")
        return
    if attrname == "color" :
        if type(value) != list and type(value) != tuple :
            if autopack.helper is not None : 
                value=autopack.helper.getMaterialProperty(value,["color"])[0]
            else :
                value = [1.,0.,0.]
    if type (value) == numpy.ndarray :
        value = value.tolist()
    elif type(value) == list or type(value) == tuple:
        for i,v in enumerate(value) :
            if type(v) == numpy.ndarray :
                value[i] = v.tolist()
            elif type(v) == list  or type(v) == tuple:
                for j,va in enumerate(v) :
                    if type(va) == numpy.ndarray :
                        v[j] = va.tolist()                        
    node.setAttribute(attrname,str(value))

def setValueToJsonNode(value,attrname):
    vdic=OrderedDict()
    if value is None:
        print (attrname, " is None !")
        return {attrname:None}
    if attrname == "color" :
        if type(value) != list and type(value) != tuple :
            if autopack.helper is not None : 
                value=autopack.helper.getMaterialProperty(value,["color"])[0]
            else :
                value = [1.,0.,0.]
    if type (value) == numpy.ndarray :
        value = value.tolist()
    elif type(value) == list or type(value) == tuple:
        for i,v in enumerate(value) :
            if type(v) == numpy.ndarray :
                value[i] = v.tolist()
            elif type(v) == list  or type(v) == tuple:
                for j,va in enumerate(v) :
                    if type(va) == numpy.ndarray :
                        v[j] = va.tolist()                        
    #node.setAttribute(attrname,str(value))
    vdic[attrname]=value
    return vdic
    
def setValueToPythonStr(value,attrname):  
    if value is None:
        print (attrname, " is None !")
        return 
    if attrname == "color" :
        if type(value) != list and type(value) != tuple :
            if autopack.helper is not None : 
                value=autopack.helper.getMaterialProperty(value,["color"])[0]
            else :
                value = [1.,0.,0.]
    if type (value) == numpy.ndarray :
        value = value.tolist()
    elif type(value) == list :
        for i,v in enumerate(value) :
            if type(v) == numpy.ndarray :
                value[i] = v.tolist()
            elif type(v) == list :
                for j,va in enumerate(v) :
                    if type(va) == numpy.ndarray :
                        v[j] = va.tolist()                        
#        print ("setValueToXMLNode ",attrname,value,str(value))  
    if type (value) == str:
        return "%s = '%s'" % (attrname,str(value))
    else :
        return "%s = %s" % (attrname,str(value))  

def getStringValueOptions(value,attrname):
    """
    Helper function to return the given environment option as a string to
    be write in the xml file.
    """
    if value is None:
        return "None"
    if attrname == "color" :
        if type(value) != list and type(value) != tuple :
            if autopack.helper is not None : 
                value=autopack.helper.getMaterialProperty(value,["color"])[0]
            else :
                value = [1.,0.,0.]
    if type (value) == numpy.ndarray :
        value = value.tolist()
    elif type(value) == list :
        for i,v in enumerate(value) :
            if type(v) == numpy.ndarray :
                value[i] = v.tolist()
            elif type(v) == list :
                for j,va in enumerate(v) :
                    if type(va) == numpy.ndarray :
                        v[j] = va.tolist() 
    if type(value) == str :
        value = '"'+value+'"'
    return str(value)     

class GrabResult(object):
    """Class for callbacks
    """

    def __init__(self,env):
        self.collision = []
        #self.lock = thread.allocate_lock()
    
    def reset(self):
        self.collision = []
        
    def grab(self, value):
        """
        the callback function
        """
        # we must use lock here because += is not atomic
        #self.lock.acquire()
        self.collision.append(value)
        #self.lock.release()

class IOingredientTool(object):
    #parser that can return an ingredient
    def __init__(self,env=None):
        super(IOingredientTool,self)
        self.env=env

    def read(self,filename):
        fileName, fileExtension = os.path.splitext(filename)
        if fileExtension == '.xml':     
            pass#self.load_XML(setupfile)
        elif fileExtension == '.py':   #execute ?  
            pass#return IOutils.load_Python(env,setupfile)
        elif fileExtension == '.json':
            pass#return IOutils.load_Json(env,setupfile) 
            
    def write(self,ingr,filename,ingr_format="xml"):
        if ingr_format == "json" :
            ingdic = self.ingrJsonNode(ingr)
            with open(filename+".json", 'w') as fp :#doesnt work with symbol link ?
                json.dump(ingdic,fp,indent=1,separators=(',', ':'))#,indent=4, separators=(',', ': ')            
        elif ingr_format == "xml" :
            ingrnode,xmldoc = self.ingrXmlNode(ingr)
            f = open(filename+".xml","w")        
            xmldoc.writexml(f, indent="\t", addindent="", newl="\n")
            f.close()
        elif ingr_format == "python" :
            ingrnode = self.ingrPythonNode(ingr)
            f = open(filename+".py","w")        
            f.write(ingrnode)
            f.close()            
        elif ingr_format == "all" :
            ingdic = self.ingrJsonNode(ingr)
            with open(filename+".json", 'w') as fp :#doesnt work with symbol link ?
                json.dump(ingdic,fp,indent=4, separators=(',', ': '))#,indent=4, separators=(',', ': ')            
            ingrnode,xmldoc = self.ingrXmlNode(ingr)
            f = open(filename+".xml","w")        
            xmldoc.writexml(f, indent="\t", addindent="", newl="\n")
            f.close()
            ingrnode = self.ingrPythonNode(ingr)
            f = open(filename+".py","w")        
            f.write(ingrnode)
            f.close()

    def makeIngredientFromXml(self,inode=None,filename=None, recipe="Generic"):
        print ("makeIngredientFromXml",inode,filename,recipe)
        if filename is None and inode is not None :
            f=str(inode.getAttribute("include"))
            if f != '':
                filename = str(f)
        if filename is not None :
            #filter the filename or pass the custom-path?
            filename=autopack.retrieveFile(filename,
                            #destination = recipe+os.sep+"recipe"+os.sep+"ingredients"+os.sep,
                            cache="recipes")
            from xml.dom.minidom import parse
            print ("parsing an ingredient xml",filename)
            xmlingr = parse(filename) # parse an XML file by name
            ingrnode = xmlingr.documentElement
        elif inode is not None:
            ingrnode = inode
        else :
            print ("filename is None")
            return None
        kw=self.parseIngrXmlNode(ingrnode)
#        print ("all KW",kw)
        #check for overwritten parameter
        overwrite_node = inode.getElementsByTagName("overwrite")
        if len(overwrite_node) :
            kwo=self.parseIngrXmlNode(overwrite_node[0])
            kw.update(kwo)
        name = str(ingrnode.getAttribute("name"))
        kw.update({"name":name})  
        ingre = self.makeIngredient(**kw)                    
        return ingre

    def parseIngrXmlNode(self,ingrnode):
#        print ("parseIngrXmlNode,",ingrnode)
#        name = str(ingrnode.getAttribute("name"))
#        print ("name",name)
        kw = {}
        for k in KWDS:
            v=getValueToXMLNode(KWDS[k]["type"],ingrnode,k)
#example of debugging...
#            if k == "rejectionThreshold" :
#                print "rejectionThreshold",KWDS[k]["type"],v,v is not None
#                print "rejectionThreshold",ingrnode.getAttribute(k)
            if v is not None :
                kw[k]=v                   
        #create the ingredient according the type
#        ingre = self1.makeIngredient(**kw)
#        kw.update({"name":name})                    
        return kw 
                
    def ingrXmlNode(self,ingr,xmldoc=None):
        rxmldoc=False
        if xmldoc is None : 
            rxmldoc = True
            impl = getDOMImplementation()
            #what about afviewer
            xmldoc = impl.createDocument(None, "ingredient", None)
            ingrnode = xmldoc.documentElement
            ingrnode.setAttribute("name",str(ingr.name))
        else :
            ingrnode = xmldoc.createElement("ingredient")
            ingrnode.setAttribute("name",str(ingr.name))
        for k in ingr.KWDS:
            v = getattr(ingr,k)
            setValueToXMLNode(v,ingrnode,k)
        if rxmldoc :
            return ingrnode,xmldoc
        else :
            return ingrnode

    def makeIngredientFromJson(self,inode=None,filename=None, recipe="Generic"):
        overwrite_dic={}
        ingr_dic={}
        if filename is None and inode is not None :
            if "include" in inode:
                filename = inode["include"]
            if "overwrite" in inode :
                overwrite_dic = inode["overwrite"]
        if filename is not None :
            filename=autopack.retrieveFile(filename,
                            #destination = recipe+os.sep+"recipe"+os.sep+"ingredients"+os.sep,
                            cache="recipes")
            with open(filename, 'r') as fp :#doesnt work with symbol link ?
                ingr_dic =json.load(fp)               
        elif inode is not None:
            ingr_dic = inode
        else :
            print ("filename is None and not ingredient dictionary provided")
            return None
        kw=ingr_dic
        #check for overwritten parameter
        if len(overwrite_dic) :
            kw.update(overwrite_dic)
#        print ("OK ",kw)
        kw=flatten_unicode_keys(kw)
        ingre = self.makeIngredient(**kw)                    
        return ingre
        
    def ingrJsonNode(self,ingr):
        ingdic={}
        for k in ingr.KWDS:
            v = getattr(ingr,k)
            if hasattr(v,"tolist"):
                v=v.tolist()
            ingdic[k] = v
        ingdic["results"]=[] 
        for r in ingr.results:  
            if hasattr(r[0],"tolist"):
                r[0]=r[0].tolist()
            if hasattr(r[1],"tolist"):
                r[1]=r[1].tolist()
            ingdic["results"].append([r[0],r[1]])
        if isinstance(ingr, GrowIngrediant) or isinstance(ingr, ActinIngrediant):
            ingdic["nbCurve"]=ingr.nbCurve
            for i in range(ingr.nbCurve):
                lp = numpy.array(ingr.listePtLinear[i])
                ingr.listePtLinear[i]=lp.tolist()                 
                ingdic["curve"+str(i)] = ingr.listePtLinear[i]
        ingdic["name"]=ingr.o_name
        return ingdic


    def ingrPythonNode(self,ingr,recipe="recipe"):
        inrStr="#include as follow : execfile('pathto/"+ingr.name+".py',globals(),{'recipe':recipe_variable_name})\n"
        if ingr.Type == "MultiSphere":
            inrStr+="from autopack.Ingredient import SingleSphereIngr, MultiSphereIngr\n"
            inrStr+=ingr.name+"= MultiSphereIngr( \n"  
        if ingr.Type == "MultiCylinder":
            inrStr+="from autopack.Ingredient import MultiCylindersIngr\n"
            inrStr+=ingr.name+"= MultiCylindersIngr( \n"                
        for k in ingr.KWDS:
            v = getattr(ingr,k)
            aStr = setValueToPythonStr(v,k)
            if aStr is not None :
                inrStr+=aStr+",\n" 
        inrStr+=")\n"
        inrStr+=recipe+".addIngredient("+ingr.name+")\n"        
        return inrStr

    def makeIngredient(self,**kw):
        from autopack.Ingredient import SingleSphereIngr, MultiSphereIngr,SingleCubeIngr
        from autopack.Ingredient import MultiCylindersIngr, GrowIngrediant
        ingr = None
        if kw["Type"]=="SingleSphere":
            kw["position"] = kw["positions"][0][0]
            kw["radius"]=kw["radii"][0][0]
            del kw["positions"]
            del kw["radii"]
            ingr = SingleSphereIngr(**kw)
        elif kw["Type"]=="MultiSphere":
            ingr = MultiSphereIngr(**kw)                    
        elif kw["Type"]=="MultiCylinder":
            ingr = MultiCylindersIngr(**kw)                    
        elif kw["Type"]=="SingleCube":
            kw["positions"]=[[[0,0,0],[0,0,0],[0,0,0],]]
            kw["positions2"]=None
            ingr = SingleCubeIngr(**kw)                    
        elif kw["Type"]=="Grow":
            ingr = GrowIngrediant(**kw)                    
        elif kw["Type"]=="Actine":
            ingr = ActinIngrediant(**kw)       
        if "gradient" in kw and kw["gradient"] != "" and kw["gradient"]!= "None":
            ingr.gradient = kw["gradient"]           
        return ingr    

    def set_recipe_ingredient(self,xmlnode,recipe):
        #get the defined ingredient
        ingrnodes = xmlnode.getElementsByTagName("ingredient")
        for ingrnode in ingrnodes:
            ingre = self.makeIngredientFromXml(inode = ingrnode)# , recipe=self.name)
            if ingre : recipe.addIngredient(ingre) 
            else : print ("PROBLEM creating ingredient from ",ingrnode,) 
        #check for includes 
        ingrnodes_include = xmlnode.getElementsByTagName("include")
        for inclnode in ingrnodes_include:
            xmlfile = str(inclnode.getAttribute("filename"))
            ingre = self.makeIngredientFromXml(filename = xmlfile)#, recipe=self.name)
            if ingre : recipe.addIngredient(ingre)
            else : print ("PROBLEM creating ingredient from ",ingrnode)
            #look for overwritten attribute
        
    
#put here the export/import ?
def save_asJson(env,setupfile,useXref=True):
        """
        Save the current environment setup as an json file.
        env is the environment / recipe to be exported.
        """
        io_ingr = IOingredientTool(env=env)
        env.setupfile = setupfile#+".json"provide the server?
        #the output path for this recipes files
        if env.setupfile.find("http") != -1 or env.setupfile.find("ftp")!= -1 :
            pathout=os.path.dirname(os.path.abspath(autopack.retrieveFile( env.setupfile)))
        else :
            pathout=os.path.dirname(os.path.abspath(env.setupfile))
        if env.version is None :
            env.version = "1.0"
        env.jsondic = OrderedDict({"recipe":{"name":env.name,"version":env.version}})
        if env.custom_paths :
            #this was the used path at loading time
            env.jsondic["recipe"]["paths"]=env.custom_paths
        env.jsondic["options"]={}
        for k in env.OPTIONS:
            v = getattr(env,k)
            if k == "gradients" :
                v = env.gradients.keys()
#            elif k == "runTimeDisplay"
            env.jsondic["options"].update(setValueToJsonNode(v,k))
        #add the boundin box
        env.jsondic["options"].update(setValueToJsonNode(env.boundingBox,"boundingBox"))

        #grid path information
        if env.grid is not None :
            if env.grid.filename is not None or env.grid.result_filename is not None:
                env.jsondic["grid"]={"grid_storage":str(env.grid.filename),
                                    "grid_result":str(env.grid.result_filename)}

        #gradient information
        if len(env.gradients):
            env.jsondic["gradients"]={}
            for gname in env.gradients:
                g = env.gradients[gname]
                env.jsondic["gradients"][str(g.name)]={}
                for k in g.OPTIONS:
                    v = getattr(g,k)
                    env.jsondic["gradients"][str(g.name)].update(setValueToJsonNode(v,k))      
        
        r =  env.exteriorRecipe
        if r :
            env.jsondic["cytoplasme"]={}
            env.jsondic["cytoplasme"]["ingredients"]={}
            for ingr in r.ingredients:                
                if useXref :
                    #write the json file for this ingredient
                    io_ingr.write(ingr,pathout+os.sep+ingr.o_name,ingr_format="json")
                    #use reference file : str(pathout+os.sep+ingr.o_name+".json")
                    ing_filename = ingr.o_name+".json"#autopack.revertOnePath(pathout+os.sep+ingr.o_name+".json")
                    env.jsondic["cytoplasme"]["ingredients"][ingr.o_name]={"name":ingr.o_name,"include":ing_filename}
                else :
                    env.jsondic["cytoplasme"]["ingredients"][ingr.o_name]={"name":ingr.o_name}
                    for k in ingr.KWDS:
                        v = getattr(ingr,k)
    #                    print ingr.name+" keyword ",k,v
                        env.jsondic["cytoplasme"]["ingredients"][ingr.o_name].update(setValueToJsonNode(v,k)) 
                    env.jsondic["cytoplasme"]["ingredients"][ingr.o_name]["name"]=ingr.o_name
        if len(env.compartments):
            env.jsondic["compartments"]=OrderedDict()
        for o in env.compartments:
            env.jsondic["compartments"][str(o.name)]=OrderedDict()
            env.jsondic["compartments"][str(o.name)]["geom"]=str(o.filename)#should point to the used filename
            env.jsondic["compartments"][str(o.name)]["name"]=str(o.ref_obj)
            if o.representation is not None :
                fileName, fileExtension = os.path.splitext(o.representation_file)
                env.jsondic["compartments"][str(o.name)]["rep"]=str(o.representation)#None
                env.jsondic["compartments"][str(o.name)]["rep_file"]=str(fileName)
#            else :
#                fileName = None
            rs = o.surfaceRecipe
            if rs :
                env.jsondic["compartments"][str(o.name)]["surface"]={}
                env.jsondic["compartments"][str(o.name)]["surface"]["ingredients"]={}
                for ingr in rs.ingredients: 
                    if useXref :
                        #write the json file for this ingredient
                        io_ingr.write(ingr,pathout+os.sep+ingr.o_name,ingr_format="json")
                        #use reference file
                        env.jsondic["compartments"][str(o.name)]["surface"]["ingredients"][ingr.o_name]={"name":ingr.o_name,"include":str(ingr.o_name+".json")}
                    else :
                        env.jsondic["compartments"][str(o.name)]["surface"]["ingredients"][ingr.o_name]={"name":ingr.o_name}
                        for k in ingr.KWDS:
                            v = getattr(ingr,k)
        #                    print ingr.name+" keyword ",k,v
                            env.jsondic["compartments"][str(o.name)]["surface"]["ingredients"][ingr.o_name].update(setValueToJsonNode(v,k)) 
                        env.jsondic["compartments"][str(o.name)]["surface"]["ingredients"][ingr.o_name]["name"]=ingr.o_name
            ri = o.innerRecipe
            if ri :
                env.jsondic["compartments"][str(o.name)]["interior"]={}
                env.jsondic["compartments"][str(o.name)]["interior"]["ingredients"]={}
                for ingr in ri.ingredients: 
                    if useXref :
                        #write the json file for this ingredient
                        io_ingr.write(ingr,pathout+os.sep+ingr.o_name,ingr_format="json")
                        #use reference file
                        env.jsondic["compartments"][str(o.name)]["interior"]["ingredients"][ingr.o_name]={"name":ingr.o_name,"include":str(ingr.o_name+".json")}
                    else :
                        env.jsondic["compartments"][str(o.name)]["interior"]["ingredients"][ingr.o_name]={"name":ingr.o_name}
                        for k in ingr.KWDS:
                            v = getattr(ingr,k)
        #                    print ingr.name+" keyword ",k,v
                            env.jsondic["compartments"][str(o.name)]["interior"]["ingredients"][ingr.o_name].update(setValueToJsonNode(v,k)) 
                        env.jsondic["compartments"][str(o.name)]["interior"]["ingredients"][ingr.o_name]["name"]=ingr.o_name
        with open(setupfile, 'w') as fp :#doesnt work with symbol link ?
            json.dump(env.jsondic,fp,indent=1,separators=(',', ':'))#,indent=4, separators=(',', ': ')
        print ("recipe saved to ",setupfile)
        
def save_asXML(env,setupfile,useXref=True):
        """
        Save the current environment setup as an xml file.
        env is the environment / recipe to be exported.
        """
        io_ingr = IOingredientTool(env=env)
#        env.setupfile = setupfile+".xml"
        pathout=os.path.dirname(os.path.abspath(env.setupfile))
        #export all information as xml
        #histovol is a tag, option are attribute of the tag
        from xml.dom.minidom import getDOMImplementation
        impl = getDOMImplementation()
        #what about afviewer
        env.xmldoc = impl.createDocument(None, "autopackSetup", None)
        root = env.xmldoc.documentElement
        root.setAttribute("name",str(env.name))
        if env.custom_paths :
            setValueToXMLNode(env.custom_paths,root,"paths")
        options=env.xmldoc.createElement("options")
        for k in env.OPTIONS:
            v = getattr(env,k)
            if k == "gradients" :
                v = env.gradients.keys()
#            elif k == "runTimeDisplay"
            setValueToXMLNode(v,options,k)
        #add the boundin box
        setValueToXMLNode(env.boundingBox,options,"boundingBox")
        setValueToXMLNode(env.version,options,"version")#version?
        root.appendChild(options)
        
        if len(env.gradients):
            gradientsnode=env.xmldoc.createElement("gradients")
            root.appendChild(gradientsnode)
            for gname in env.gradients:
                g = env.gradients[gname]
                grnode = env.xmldoc.createElement("gradient")
                gradientsnode.appendChild(grnode)
                grnode.setAttribute("name",str(g.name))
                for k in g.OPTIONS:
                    v = getattr(g,k)
                    setValueToXMLNode(v,grnode,k)      

        #grid path information
        if env.grid.filename is not None or env.grid.result_filename is not None:
            gridnode=env.xmldoc.createElement("grid")
            root.appendChild(gridnode)
            gridnode.setAttribute("grid_storage",str(env.grid.filename))
            gridnode.setAttribute("grid_result",str(env.grid.result_filename))
        
        r =  env.exteriorRecipe
        if r :
            rnode=env.xmldoc.createElement("cytoplasme")
            root.appendChild(rnode)
            for ingr in r.ingredients:                
                if useXref :
                    io_ingr.write(ingr,pathout+os.sep+ingr.name,ingr_format="xml")
                    ingrnode = env.xmldoc.createElement("ingredient")
                    rnode.appendChild(ingrnode)
                    ingrnode.setAttribute("include",str(pathout+os.sep+ingr.name+".xml"))                    
                else :
                    ingrnode = env.xmldoc.createElement("ingredient")
                    rnode.appendChild(ingrnode)
                    ingrnode.setAttribute("name",str(ingr.name))
                    for k in ingr.KWDS:
                        v = getattr(ingr,k)
    #                    print ingr.name+" keyword ",k,v
                        setValueToXMLNode(v,ingrnode,k)
        for o in env.compartments:
            onode=env.xmldoc.createElement("compartment")
            root.appendChild(onode)
            onode.setAttribute("name",str(o.name))
            onode.setAttribute("geom",str(o.filename))#should point to the used filename
            onode.setAttribute("rep",str(o.representation))#None
            if o.representation is not None :
                fileName, fileExtension = os.path.splitext(o.representation_file)
            else :
                fileName = None
            onode.setAttribute("rep_file",str(fileName))#None
            rs = o.surfaceRecipe
            if rs :
                onodesurface=env.xmldoc.createElement("surface")
                onode.appendChild(onodesurface)
                for ingr in rs.ingredients: 
                    if useXref :
                        io_ingr.write(ingr,pathout+os.sep+ingr.name,ingr_format="xml")
                        ingrnode = env.xmldoc.createElement("ingredient")
                        onodesurface.appendChild(ingrnode)
                        ingrnode.setAttribute("include",str(pathout+os.sep+ingr.name+".xml")) 
                    else :
                        ingrnode = env.xmldoc.createElement("ingredient")
                        onodesurface.appendChild(ingrnode)
                        ingrnode.setAttribute("name",str(ingr.name))                       
                        for k in ingr.KWDS:
                            v = getattr(ingr,k)
                            setValueToXMLNode(v,ingrnode,k)
            ri = o.innerRecipe
            if ri :
                onodeinterior=env.xmldoc.createElement("interior")
                onode.appendChild(onodeinterior)             
                for ingr in ri.ingredients: 
                    if useXref :
                        io_ingr.write(ingr,pathout+os.sep+ingr.name,ingr_format="xml")
                        ingrnode = env.xmldoc.createElement("ingredient")
                        onodeinterior.appendChild(ingrnode)
                        ingrnode.setAttribute("include",str(pathout+os.sep+ingr.name+".xml")) 
                    else :
                        ingrnode = env.xmldoc.createElement("ingredient")
                        onodeinterior.appendChild(ingrnode)
                        ingrnode.setAttribute("name",str(ingr.name))                       
                        for k in ingr.KWDS:
                            v = getattr(ingr,k)
                            setValueToXMLNode(v,ingrnode,k)
        f = open(setupfile,"w")        
        env.xmldoc.writexml(f, indent="\t", addindent="", newl="\n")
        f.close()

def save_asPython(env,setupfile,useXref=True):
    """
    Save the current environment setup as a python script file.
    """
    io_ingr = IOingredientTool(env=env)
    env.setupfile = setupfile
    pathout=os.path.dirname(os.path.abspath(env.setupfile))
    #add the import statement
    setupStr="""
import sys
import os
#autopack
import autopack
localdir = wrkDir = autopack.__path__[0]
from autopack.Ingredient import SingleSphereIngr, MultiSphereIngr
from autopack.Ingredient import MultiCylindersIngr,GrowIngrediant,ActinIngrediant
from autopack.Compartment import Compartment
from autopack.Recipe import Recipe
from autopack.Environment import Environment
from autopack.Graphics import AutopackViewer as AFViewer
#access the helper
helper = autopack.helper
if helper is None :
import upy
helperClass = upy.getHelperClass()
helper =helperClass()
#create the viewer
ViewerType=autopack.helper.host    
afviewer = AFViewer(ViewerType=helper.host,helper=helper)#long ?
#make some option here     
afviewer.doPoints = False
afviewer.doSpheres = False
afviewer.quality = 1 #lowest quality for sphere and cylinder
afviewer.visibleMesh = True #mesh default visibility 
#create the env
h1 = Environment()
"""
    setupStr+="h1.name='"+env.name+"'\n"    
    for k in env.OPTIONS:
        v = getattr(env,k)
        if k == "gradients" :
            v = env.gradients.keys()
        vstr=getStringValueOptions(v,k)#env.setValueToXMLNode(v,options,k)
        setupStr+="h1.%s=%s\n" % (k,vstr)
    #add the boundin box
    vstr=getStringValueOptions(env.boundingBox,"boundingBox")#env.setValueToXMLNode(v,options,k)
    setupStr+="h1.%s=%s\n" % ("boundingBox",vstr)
    vstr=getStringValueOptions(env.version,k)#env.setValueToXMLNode(v,options,k)
    setupStr+="h1.%s=%s\n" % ("version",vstr)
    
#TODO : GRADIENT
#        if len(env.gradients):
#            gradientsnode=env.xmldoc.createElement("gradients")
#            root.appendChild(gradientsnode)
#            for gname in env.gradients:
#                g = env.gradients[gname]
#                grnode = env.xmldoc.createElement("gradient")
#                gradientsnode.appendChild(grnode)
#                grnode.setAttribute("name",str(g.name))
#                for k in g.OPTIONS:
#                    v = getattr(g,k)
#                    env.setValueToXMLNode(v,grnode,k)      
#
#        grid path information
#        if env.grid.filename is not None or env.grid.result_filename is not None:
#            gridnode=env.xmldoc.createElement("grid")
#            root.appendChild(gridnode)
#            gridnode.setAttribute("grid_storage",str(env.grid.filename))
#            gridnode.setAttribute("grid_result",str(env.grid.result_filename))
#        
    r =  env.exteriorRecipe
    if r :
        setupStr+="cytoplasme = Recipe()\n"
        for ingr in r.ingredients:                
            if useXref :
                io_ingr.write(ingr,pathout+os.sep+ingr.name,ingr_format="python")
                setupStr+="execfile('"+pathout+os.sep+ingr.name+".py',globals(),{'recipe':cytoplasme})\n"
            else :
                ingrnode = io_ingr.ingrPythonNode(ingr,recipe="cytoplasme")
                setupStr+=ingrnode
        setupStr+="h1.setExteriorRecipe(cytoplasme)\n"                    
    for o in env.compartments:
        setupStr+=o.name+" = Compartment('"+o.name+"',None, None, None,\n"
        setupStr+="         filename='"+o.filename+"',\n"
        if o.representation is not None:
            setupStr+="         object_name ='"+o.representation+"',\n"
            setupStr+="         object_filename ='"+o.representation_file+"'\n"
        setupStr+="         )\n"
        setupStr+="h1.addCompartment("+o.name+")\n"
        rs = o.surfaceRecipe
        if rs :
            setupStr+=o.name+"_surface = Recipe(name='"+o.name+"_surf')\n"
            for ingr in rs.ingredients:                
                if useXref :
                    io_ingr.write(ingr,pathout+os.sep+ingr.name,ingr_format="python")
                    setupStr+="execfile('"+pathout+os.sep+ingr.name+".py',globals(),{'recipe':"+o.name+"_surface})\n"
                else :
                    ingrnode = io_ingr.ingrPythonNode(ingr,recipe=o.name+"_surface")
                    setupStr+=ingrnode 
            setupStr+=o.name+".setSurfaceRecipe("+o.name+"_surface)\n"
        ri = o.innerRecipe
        if ri :
            setupStr+=o.name+"_inner = Recipe(name='"+o.name+"_int')\n"
            for ingr in rs.ingredients:                
                if useXref :
                    io_ingr.write(ingr,pathout+os.sep+ingr.name,ingr_format="python")
                    setupStr+="execfile('"+pathout+os.sep+ingr.name+".py',globals(),{'recipe':"+o.name+"_inner})\n"
                else :
                    ingrnode = io_ingr.ingrPythonNode(ingr,recipe=o.name+"_inner")
                    setupStr+=ingrnode   
            setupStr+=o.name+".setInnerRecipe("+o.name+"_inner)\n"
    setupStr+="afviewer.SetHistoVol(h1,0,display=False)\n"
    setupStr+="afviewer.displayPreFill()\n"
    setupStr+="bbox = afviewer.helper.getObject('histvolBB')\n"
    setupStr+="if bbox is None : bbox = afviewer.helper.box('histvolBB',cornerPoints=h1.boundingBox)\n"
    setupStr+="helper = afviewer.helper\n"
    setupStr+="noGUI = False\n"
    setupStr+="try :\n"
    setupStr+="    print ('try')\n"
    setupStr+="    AFGui.Set('"+env.name+"',helper=afviewer.helper,afviewer=afviewer,histoVol=h1,bbox=bbox)\n"
    setupStr+="except:\n"
    setupStr+="    print ('no GUI')\n"
    setupStr+="    noGUI = True\n"
    f = open(setupfile,"w")        
    f.write(setupStr)
    f.close()

def load_XML(env,setupfile):
    """
    Setup the environment according the given xml file. 
    """
    env.setupfile = setupfile
#    from autopack import Ingredient as ingr
    io_ingr = IOingredientTool(env=env)
    from xml.dom.minidom import parse
    env.xmldoc = parse(setupfile) # parse an XML file by name
    root = env.xmldoc.documentElement
    env.name = str(root.getAttribute("name"))
    env.custom_paths = getValueToXMLNode("g",root,"paths")
    env.current_path=os.path.dirname(os.path.abspath(env.setupfile))   
    if env.custom_paths:
#        autopack.replace_path.extend(env.custom_paths)#keyWordPAth,valuePath
        autopack.updateReplacePath(env.custom_paths)
    autopack.current_recipe_path=env.current_path
    options=root.getElementsByTagName("options")
    if len(options) :
        options=options[0]
        for k in env.OPTIONS:
            if k == "gradients": 
                continue
            v=getValueToXMLNode(env.OPTIONS[k]["type"],options,k)
            if v is not None :
                setattr(env,k,v)
        v=getValueToXMLNode("vector",options,"boundingBox")
        env.boundingBox = v
        v=getValueToXMLNode("string",options,"version")
        env.version = v
        
    gradientsnode=root.getElementsByTagName("gradients")
    if len(gradientsnode) :
        gradientnode=gradientsnode[0]
        grnodes = gradientnode.getElementsByTagName("gradient")
        for grnode in grnodes:
            name = str(grnode.getAttribute("name"))
            mode = str(grnode.getAttribute("mode"))
            weight_mode = str(grnode.getAttribute("weight_mode"))
            pick_mode = str(grnode.getAttribute("pick_mode"))
            direction = str(grnode.getAttribute("direction"))#vector
            description=str(grnode.getAttribute("description"))
            radius=float(str(grnode.getAttribute("radius")))
#                print "weight_mode",weight_mode
            env.setGradient(name=name,mode=mode, direction=eval(direction),
                        weight_mode=weight_mode,description=description,
                        pick_mode=pick_mode,radius=radius)

    gridnode=root.getElementsByTagName("grid")
    if len(gridnode) :
        gridn=gridnode[0]
        env.grid_filename = str(gridn.getAttribute("grid_storage"))
        env.grid_result_filename = str(gridn.getAttribute("grid_result"))
    
    from autopack.Recipe import Recipe
    rnode=root.getElementsByTagName("cytoplasme")
    if len(rnode) :
        rCyto = Recipe()
        rnode=rnode[0]
        #check for include list of ingredients
        ingredients_xmlfile = str(rnode.getAttribute("include"))
        if ingredients_xmlfile :#open the file and parse the ingredient:
            #check if multiple include filename, aumngo',' in the path
            liste_xmlfile = ingredients_xmlfile.split(",")
            for xmlf in liste_xmlfile :                    
                xmlfile = autopack.retrieveFile(xmlf,
                        #destination = self.name+os.sep+"recipe"+os.sep,
                        cache="recipes")
                if xmlfile :
                    xmlinclude = parse(xmlfile).documentElement
                    io_ingr.set_recipe_ingredient(xmlinclude,rCyto)
            
        io_ingr.set_recipe_ingredient(rnode,rCyto)                
        #setup recipe
        env.setExteriorRecipe(rCyto)
        
    onodes = root.getElementsByTagName("compartment")#Change to Compartment
    if not len(onodes) :
        #backward compatibility
        onodes = root.getElementsByTagName("organelle")#Change to Compartment            
    from autopack.Compartment import Compartment
    for onode in onodes:
        name = str(onode.getAttribute("name"))
        geom = str(onode.getAttribute("geom"))
        rep =  str(onode.getAttribute("rep"))
        rep_file=str(onode.getAttribute("rep_file"))
        #print ("is it working")
        print ("xml parsing ",name,geom,rep,rep_file)
        #print (rep,rep_file,len(rep),rep == '',rep=="",rep != "None",rep != "None" or len(rep) != 0)
        if (rep != "None" and len(rep) != 0) :
            rname =  rep_file.split("/")[-1]
            fileName, fileExtension = os.path.splitext(rname)
            if fileExtension == "" :
                fileExtension = autopack.helper.hext
                if fileExtension == "" :
                    rep_file = rep_file+fileExtension
                else :
                    rep_file = rep_file+"."+fileExtension   
        else :
            rep=None
            rep_file=None
            print ("no representation found")
        print ("add compartment ",name,geom,rep,rep_file)
        o = Compartment(name,None, None, None,filename=geom,object_name=rep,object_filename=rep_file)
        print ("added compartment ",name)
        env.addCompartment(o)
        rsnodes = onode.getElementsByTagName("surface")
        from autopack.Recipe import Recipe
        if len(rsnodes) :
            rSurf = Recipe(name=o.name+"_surf")
            rsnodes=rsnodes[0]
            ingredients_xmlfile = str(rsnodes.getAttribute("include"))
            if ingredients_xmlfile :#open the file and parse the ingredient:
                #check if multiple include filename, aumngo',' in the path
                liste_xmlfile = ingredients_xmlfile.split(",")
                for xmlf in liste_xmlfile :                    
                    xmlfile = autopack.retrieveFile(xmlf,
                            #destination = self.name+os.sep+"recipe"+os.sep,
                            cache="recipes")
                    if xmlfile :
                        xmlinclude = parse(xmlfile).documentElement
                        io_ingr.set_recipe_ingredient(xmlinclude,rSurf)
            io_ingr.set_recipe_ingredient(rsnodes,rSurf)                
            o.setSurfaceRecipe(rSurf)                
        rinodes = onode.getElementsByTagName("interior")
        from autopack.Recipe import Recipe
        if len(rinodes) :
            rMatrix = Recipe(name=o.name+"_int")
            rinodes=rinodes[0]
            ingredients_xmlfile = str(rinodes.getAttribute("include"))
            if ingredients_xmlfile :#open the file and parse the ingredient:
                #check if multiple include filename, aumngo',' in the path
                liste_xmlfile = ingredients_xmlfile.split(",")
                for xmlf in liste_xmlfile :
                    xmlfile = autopack.retrieveFile(xmlf,
                            #destination = self.name+os.sep+"recipe"+os.sep,
                        cache="recipes")
                    if xmlfile :
                        xmlinclude = parse(xmlfile).documentElement
                        io_ingr.set_recipe_ingredient(xmlinclude,rMatrix)
            io_ingr.set_recipe_ingredient(rinodes,rMatrix)                
            o.setInnerRecipe(rMatrix)
    #Go through all ingredient and setup the partner  
    env.loopThroughIngr(env.set_partners_ingredient)
#        if self.placeMethod.find("panda") != -1 :
#            self.setupPanda()
        


def load_Json(env,setupfile):
    """
    Setup the environment according the given xml file. 
    """
    if setupfile == None:
        setupfile = env.setupfile
    with open(setupfile, 'r') as fp :#doesnt work with symbol link ?
        env.jsondic=json.load(fp,object_pairs_hook=OrderedDict)#,indent=4, separators=(',', ': ')
    env.current_path=os.path.dirname(os.path.abspath(env.setupfile))   
    from autopack import Ingredient as ingr
    io_ingr = IOingredientTool(env=env)
    env.name = env.jsondic["recipe"]["name"]
    env.version = env.jsondic["recipe"]["version"]
    #is there any cutoms paths
    if "paths" in env.jsondic["recipe"]:
        env.custom_paths = env.jsondic["recipe"]["paths"]#list(env.jsondic["recipe"]["paths"].items())
#        autopack.replace_path.extend(env.custom_paths)#keyWordPAth,valuePath
        autopack.updateReplacePath(env.custom_paths)
    autopack.current_recipe_path=env.current_path    
    options=env.jsondic["options"]
    if len(options) :
        for k in env.OPTIONS:
            if k == "gradients": 
                continue
            setattr(env,k,options[k])
        env.boundingBox = options["boundingBox"]
    if "gradients" in env.jsondic:
        gradientsnode=env.jsondic["gradients"]
        if len(gradientsnode) :#number of gradients defined
            for g_name in gradientsnode:
                g_dic = gradientsnode[g_name]
                env.setGradient(name=g_name,mode=g_dic["mode"], 
                            direction=g_dic["direction"],weight_mode=g_dic["weight_mode"],
                            description=g_dic["description"],
                            pick_mode=g_dic["pick_mode"],radius=g_dic["radius"])

    if "grid" in env.jsondic:
        gridnode=env.jsondic["grid"]
        if len(gridnode) :
            env.grid_filename = str(gridnode["grid_storage"])
            env.grid_result_filename = str(gridnode["grid_result"])
    from autopack.Recipe import Recipe
    if "cytoplasme" in env.jsondic:
        rnode=env.jsondic["cytoplasme"]
        ingrs_dic = env.jsondic["cytoplasme"]["ingredients"]
        if len(ingrs_dic) :
            rCyto = Recipe()
            for ing_name in ingrs_dic:
                #either xref or defined
                ing_dic = ingrs_dic[ing_name]
                ingr=io_ingr.makeIngredientFromJson(inode=ing_dic,recipe=env.name)
                rCyto.addIngredient(ingr) 
            #setup recipe
            env.setExteriorRecipe(rCyto)
    
    if  "compartments" in env.jsondic:    
        if len(env.jsondic["compartments"]):
            from autopack.Compartment import Compartment
            for cname in env.jsondic["compartments"]:
                comp_dic = env.jsondic["compartments"][cname]
                name = str(comp_dic["name"])
                geom = str(comp_dic["geom"])
                rep=""
                if "rep" in comp_dic :
                    rep =  str(comp_dic["rep"])
                rep_file=""
                if "rep_file" in comp_dic:
                    rep_file=str(comp_dic["rep_file"])
                print ("rep ?",name,geom,rep,rep_file,(rep != "None" and len(rep) != 0 and rep != '' and rep == ""))
#                print (len(rep),rep == '',rep=="",rep != "None",rep != "None" or len(rep) != 0)
                if rep != "None" and len(rep) != 0 and rep != '' and rep != "":
                    rname =  rep_file.split("/")[-1]
                    fileName, fileExtension = os.path.splitext(rname)
                    if fileExtension == "" :
                        fileExtension = autopack.helper.hext
                        if fileExtension == "" :
                            rep_file = rep_file+fileExtension
                        else :
                            rep_file = rep_file+"."+fileExtension   
                else :
                    rep=None
                    rep_file=None
                    print ("NONENE")
                print ("add compartment ",name,geom,rep,rep_file)
                o = Compartment(name,None, None, None,filename=geom,
                                object_name=rep,object_filename=rep_file)
                print ("added compartment ",name)
                env.addCompartment(o)
                if "surface" in comp_dic:
                    snode=comp_dic["surface"]
                    ingrs_dic = snode["ingredients"]
                    if len(ingrs_dic) :
                        rSurf = Recipe(name="surf_"+str(len(env.compartments)-1))
#                        rSurf = Recipe(name=o.name+"_surf")
                        for ing_name in ingrs_dic:
                            #either xref or defined
                            ing_dic = ingrs_dic[ing_name]
                            ingr=io_ingr.makeIngredientFromJson(inode=ing_dic,recipe=env.name)
                            rSurf.addIngredient(ingr) 
                        #setup recipe
                        o.setSurfaceRecipe(rSurf)                
                if "interior" in comp_dic:
                    snode=comp_dic["interior"]
                    ingrs_dic = snode["ingredients"]
                    if len(ingrs_dic) :
#                        rMatrix = Recipe(name=o.name+"_int")
                        rMatrix = Recipe(name="int_"+str(len(env.compartments)-1))
                        for ing_name in ingrs_dic:
                            #either xref or defined
                            ing_dic = ingrs_dic[ing_name]
                            ingr=io_ingr.makeIngredientFromJson(inode=ing_dic,recipe=env.name)
                            rMatrix.addIngredient(ingr) 
                        #setup recipe
                        o.setInnerRecipe(rMatrix)                
    #Go through all ingredient and setup the partner  
    env.loopThroughIngr(env.set_partners_ingredient)
#        if env.placeMethod.find("panda") != -1 :
#            env.setupPanda()
        
