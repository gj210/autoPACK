/*
###############################################################################
#
# autoPACK Authors: Graham T. Johnson, Mostafa Al-Alusi, Ludovic Autin, Michel Sanner
#   Based on COFFEE Script developed by Graham Johnson between 2005 and 2010 
#   with assistance from Mostafa Al-Alusi in 2009 and periodic input 
#   from Arthur Olson's Molecular Graphics Lab
#
# autopack.cpp Authors: Ludovic Autin
#
# Translation from Python initiated March 15, 2010 by Ludovic Autin
#
#
# Copyright: Graham Johnson Ludovic Autin ©2010
#
# This file "autopack.cpp" is part of autoPACK, cellPACK.
#    
#    autoPACK is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    autoPACK is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with autoPACK (See "CopyingGNUGPL" in the installation.
#    If not, see <http://www.gnu.org/licenses/>.
#
#
###############################################################################
@author: Graham Johnson, Ludovic Autin, & Michel Sanner
*/

#include <vector>
#include <map>
#include <pugixml.hpp>

#include "BigGrid.h"
#include "Sphere.h"
#include "Types.h"
#include "IngradientsFactory.h"

/* XML CODE */
/* 
parsing information from the autopack setup file as well as the collada mesh file 
*/

double getRadii(std::string str){
    //std::string str(input);  
    std::replace(str.begin(), str.end(), '[', ' ');  
    std::replace(str.begin(), str.end(), ']', ' ');  
    std::replace(str.begin(), str.end(), '(', ' ');  
    std::replace(str.begin(), str.end(), ')', ' ');  
    return atof(str.c_str());
}

std::vector<float> getBox(std::string str){
    std::replace(str.begin(), str.end(), '[', ' ');  
    std::replace(str.begin(), str.end(), ']', ' ');  
    std::replace(str.begin(), str.end(), '(', ' ');  
    std::replace(str.begin(), str.end(), ')', ' ');  
    std::replace(str.begin(), str.end(), ',', ' ');
    // Build an istream that holds the input string
    std::istringstream iss(str);
    //std::cout << str << std::endl;
    // If possible, always prefer std::vector to naked array
    std::vector<float> v;
    // Iterate over the istream, using >> to grab floats
    // and push_back to store them in the vector
    std::copy(std::istream_iterator<float>(iss),
        std::istream_iterator<float>(),
        std::back_inserter(v));
    return v;
}

std::vector<int> getInts(std::string str){
    // Build an istream that holds the input string
    std::istringstream iss(str);
    //std::cout << str << std::endl;
    // If possible, always prefer std::vector to naked array
    std::vector<int> v;
    // Iterate over the istream, using >> to grab floats
    // and push_back to store them in the vector
    std::copy(std::istream_iterator<int>(iss),
        std::istream_iterator<int>(),
        std::back_inserter(v));
    return v;
}

std::vector<openvdb::Vec3f> getPositions(std::vector<float> pos){
    std::vector<float>::size_type i = 0; 
    std::vector<openvdb::Vec3f> positions;   
    while (i<pos.size()-2){
         openvdb::Vec3f p(pos[i],pos[i+1],pos[i+2]);
         positions.push_back(p);
         i = i+3;
    }
    return positions;
}
std::vector<openvdb::Vec3s> getPositionsS(std::vector<float> pos){
    std::vector<float>::size_type i = 0; 
    std::vector<openvdb::Vec3s> positions;   
    while (i<pos.size()-2){
         openvdb::Vec3f p(pos[i],pos[i+1],pos[i+2]);
         positions.push_back(p);
         i = i+3;
    }
    return positions;
}

std::vector<openvdb::Vec3I> getPositionsInt(std::vector<int> pos){
    std::vector<int>::size_type i = 0; 
    std::vector<openvdb:: Vec3I > positions;   
    while (i<pos.size()-2){
         //openvdb::Vec3i p(pos[i],pos[i+1],pos[i+2]);
         openvdb:: Vec3I  p(pos[i],pos[i+1],pos[i+2]);//, openvdb::util::INVALID_IDX);
         positions.push_back(p);
         i = i+3;
    }
    return positions;
}

openvdb::Vec3f getArray(std::string str){
    std::replace(str.begin(), str.end(), '[', ' ');  
    std::replace(str.begin(), str.end(), ']', ' ');  
    std::replace(str.begin(), str.end(), '(', ' ');  
    std::replace(str.begin(), str.end(), ')', ' ');  
    std::replace(str.begin(), str.end(), ',', ' ');
    // Build an istream that holds the input string
    std::istringstream iss(str);
    //std::cout << str << std::endl;
    // If possible, always prefer std::vector to naked array
    std::vector<float> v;
    // Iterate over the istream, using >> to grab floats
    // and push_back to store them in the vector
    std::copy(std::istream_iterator<float>(iss),
        std::istream_iterator<float>(),
        std::back_inserter(v));
    openvdb::Vec3f p(v[0],v[1],v[2]);
    return p;
}
/* 
assimp test, but doesnt work for all colldada file

//use Open Asset Import BSD Licence
//see here for repo woth xcode project https://bitbucket.org/sherief/open-asset-import/src
void getMeshs_assimp_node(const struct aiScene* scene,const struct aiNode* nd,
                            struct aiMatrix4x4* trafo,
                            std::vector<mesh>* meshs){
    std::cout << "# getNode " << std::endl;    
    struct aiMatrix4x4 prev;
    unsigned int n = 0, t;
    prev = *trafo;
    *trafo = prev * nd->mTransformation;
    //aiMultiplyMatrix4(trafo,&nd->mTransformation);
    std::cout << "# mNumMeshes is " <<  nd->mNumMeshes << std::endl;
    for (; n < nd->mNumMeshes; ++n) {
        mesh mesh3d;  
        std::vector<openvdb::Vec3s> vertices;
        std::vector<openvdb::Vec3I> faces;
        std::vector<openvdb::Vec4I> quads;
        const struct aiMesh* m = scene->mMeshes[nd->mMeshes[n]];
        for (t = 0; t < m->mNumVertices; ++t) {
            struct aiVector3D tmp = m->mVertices[t];//(*trafo);
            tmp*=*trafo;
            //aiTransformVecByMatrix4(&tmp,trafo);
            vertices.push_back(openvdb::Vec3s(tmp.x,tmp.y,tmp.z));
        // get the vertices and fill mesh         
        }
        for (t = 0; t < m->mNumFaces; ++t) {
            struct aiFace f = m->mFaces[t];
            if (f.mNumIndices == 3) faces.push_back(openvdb::Vec3I(f.mIndices[0],f.mIndices[1],f.mIndices[2]));
            else if (f.mNumIndices == 3) quads.push_back(openvdb::Vec4I(f.mIndices[0],f.mIndices[1],f.mIndices[2],f.mIndices[3]));
        }
        mesh3d.vertices = vertices; 
        mesh3d.faces = faces;
        mesh3d.quads = quads;
        meshs->push_back(mesh3d);
        std::cout << "# nmesh is " <<  meshs->size() << std::endl;   
   }
    std::cout << "# nchildren is " <<  nd->mNumChildren << std::endl;
    for (n = 0; n < nd->mNumChildren; ++n) {
        getMeshs_assimp_node(scene,nd->mChildren[n],trafo,meshs);
    }
    *trafo = prev;
}

std::vector<mesh> getMeshs_assimp(std::string path){
    unsigned int n = 0, t;
    std::vector<mesh> meshs;
    const aiScene* scene = NULL;
    std::cout << "# scene to load  " <<  path << std::endl; 
     // Create an instance of the Importer class
    //Assimp::Importer imp;
    //importer.SetExtraVerbose(true);
    // And have it read the given file with some example postprocessing
    // Usually - if speed is not the most important aspect for you - you'll
    // propably to request more postprocessing than we do in this example.
    //importer.GetOrphanedScene() ;    
    //importer.FreeScene( );
    const unsigned int flags = 
        aiProcess_Triangulate |
        //aiProcess_JoinIdenticalVertices |
        //aiProcess_GenSmoothNormals |
        //aiProcess_ValidateDataStructure |
        //aiProcess_RemoveRedundantMaterials |
        //aiProcess_SortByPType |
        //aiProcess_FindDegenerates |
        //aiProcess_FindInvalidData |;
        aiProcess_GenUVCoords;// |
        //aiProcess_OptimizeMeshes |
        //aiProcess_OptimizeGraph;    
    scene = importer.ReadFile(path.c_str(),flags);
    scene = importer.GetScene();
    if( !scene)
    {
        std::cout << "#X" << importer.GetErrorString() << "X " << std::endl;
        return meshs;
    }
    //scene = Assimp::aiImportFile(path.c_str(),flags);//aiProcessPreset_TargetRealtime_Quality);//require const char*
    std::cout << "# scene loaded " << std::endl;
    //std::cout << " scene loaded  and has mesh ? " <<  scene->HasMeshes() << std::endl;    //segfault ?
    //for every node ?recursively from scene->mRootNode
    //struct aiNode* nd = scene->mRootNode;
    struct aiMatrix4x4 trafo;
    //aiIdentityMatrix4(&trafo);
    std::cout << "# scene loaded2 " << std::endl;
    getMeshs_assimp_node(scene,scene->mRootNode,&trafo,&meshs);
    //importer.GetOrphanedScene() ;    
    //importer.FreeScene( );  
    //delete scene;
    //Assimp::aiReleaseImport(scene); 
    return meshs;
}
*/
//
openvdb::math::Mat4d getNodeTransformation(pugi::xml_node nd){
    openvdb::math::Transform::Ptr transfo =
            openvdb::math::Transform::createLinearTransform();
    if (std::string(nd.name()) == "visual_scene")
        return transfo->baseMap()->getAffineMap()->getMat4();//should be identity    
    std::string stringT( nd.child("translate").text().get() );
    std::vector<float> pos = getBox(stringT);    
    if (DEBUG) std::cout << "# trans  " <<stringT<< std::endl; 
    //std::vector<openvdb::Vec3s> translate = getPositionsS(pos);
    openvdb::math::Mat4d T;
    T.setIdentity();
    openvdb::Vec3d trans(pos[0],pos[1],pos[2]);
    if (trans.length() == 0.0) T.setIdentity();
    else T.setTranslation(trans);
    openvdb::math::Mat4d rX;
    openvdb::math::Mat4d rY;
    openvdb::math::Mat4d rZ;
    rX.setIdentity();
    rY.setIdentity();
    rZ.setIdentity();
    //rotation 
    for (pugi::xml_node rotation = nd.child("rotate"); rotation; rotation = rotation.next_sibling("rotate")){
        std::string stringR( rotation.text().get() );
        std::vector<float> r = getBox(stringR);  
        if (DEBUG) std::cout << "#rotation  " <<stringR<< std::endl;           
        if (std::string(rotation.attribute("sid").value()) =="rotateY"){
            rY.setToRotation(openvdb::Vec3d(r[0],r[1],r[2]),(r[3]/180.0)*M_PI);
        }
        else if (std::string(rotation.attribute("sid").value()) =="rotateX"){
            rX.setToRotation(openvdb::Vec3d(r[0],r[1],r[2]),(r[3]/180.0)*M_PI);        
        }
        else if (std::string(rotation.attribute("sid").value()) =="rotateZ"){
            rZ.setToRotation(openvdb::Vec3d(r[0],r[1],r[2]),(r[3]/180.0)*M_PI);
        }
    }    
    //transfo->postTranslate(trans);
    //transfo->postMult(T*rY*rX*rZ); 
    transfo->preMult(rZ*rX*rY*T); 
    return transfo->baseMap()->getAffineMap()->getMat4();
}
void getMeshs_nodew(pugi::xml_node nd,
                            openvdb::math::Mat4d* trafo,
                            std::map<std::string, mesh> scene_meshs,
                            std::vector<mesh>* meshs){
    if (DEBUG) std::cout << "#getNode " <<nd.attribute("id").value()<< std::endl;    
    openvdb::math::Mat4d prev;
    unsigned int n = 0, t;
    prev = *trafo;
    if (DEBUG) std::cout << "#getNodeTransformation * prev " <<prev<< std::endl;
    //get node transformation  
    *trafo = prev * getNodeTransformation(nd) ;
    //if (DEBUG) std::cout << "trafo= " << *trafo << std::endl;
    pugi::xml_node igeom_node = nd.child("instance_geometry");
    if (igeom_node) {
        if (DEBUG) std::cout << "#geomnode " << std::endl;
        std::string idsource(igeom_node.attribute("url").value());
        //remove the #
        idsource.erase(std::remove(idsource.begin(), idsource.end(), '#'), idsource.end()); 
        //gett he corresponding mesh and apply the matrice
        mesh mesh3d;  
        if (DEBUG) std::cout << "trafo= " << *trafo << std::endl;
        std::vector<openvdb::Vec3s> vertices;
        std::vector<openvdb::Vec3I> faces;
        std::vector<openvdb::Vec4I> quads;
        mesh m = scene_meshs[idsource];
        for (t = 0; t < m.vertices.size(); ++t) {
            //openvdb::Vec3s tmp =  m.vertices[t];//*(*trafo);
            // tmp=*trafo;
            //aiTransformVecByMatrix4(&tmp,trafo);
            vertices.push_back(m.vertices[t]*(*trafo));
        // get the vertices and fill mesh         
        }
        mesh3d.vertices = vertices; 
        mesh3d.faces = m.faces;
        mesh3d.quads = m.quads;
        meshs->push_back(mesh3d);
        if (DEBUG) std::cout << "#nmesh is " <<  meshs->size() << std::endl;   
    }
    for (pugi::xml_node node = nd.child("node"); node; node = node.next_sibling("node")){
        std::cout << "#" <<nd.attribute("id").value() << "child " << node.attribute("id").value() << std::endl;
        getMeshs_nodew(node,trafo,scene_meshs,meshs);
    }
    trafo = &prev;
}

void getMeshs_node(pugi::xml_node nd,
                            openvdb::math::Mat4d trafo,
                            std::map<std::string, mesh> scene_meshs,
                            std::vector<mesh>* meshs){
    if (DEBUG) std::cout << "#getNode " <<nd.attribute("id").value()<< std::endl;    
    openvdb::math::Mat4d prev;
    unsigned int n = 0, t;
    prev = trafo;
    if (DEBUG) std::cout << "prev =" <<prev<< std::endl;
    //get node transformation  
    trafo = prev * getNodeTransformation(nd) ;
    //if (DEBUG) std::cout << "trafo= " << *trafo << std::endl;
    pugi::xml_node igeom_node = nd.child("instance_geometry");
    if (igeom_node) {
        if (DEBUG) std::cout << "#geomnode " << std::endl;
        std::string idsource(igeom_node.attribute("url").value());
        //remove the #
        idsource.erase(std::remove(idsource.begin(), idsource.end(), '#'), idsource.end()); 
        //gett he corresponding mesh and apply the matrice
        mesh mesh3d;  
        if (DEBUG) std::cout << "trafo= " << trafo << std::endl;
        std::vector<openvdb::Vec3s> vertices;
        std::vector<openvdb::Vec3I> faces;
        std::vector<openvdb::Vec4I> quads;
        mesh m = scene_meshs[idsource];
        for (t = 0; t < m.vertices.size(); ++t) {
            //openvdb::Vec3s tmp =  m.vertices[t];//*(*trafo);
            // tmp=*trafo;
            //aiTransformVecByMatrix4(&tmp,trafo);
            vertices.push_back(m.vertices[t]*trafo);
        // get the vertices and fill mesh         
        }
        mesh3d.vertices = vertices; 
        mesh3d.faces = m.faces;
        mesh3d.quads = m.quads;
        meshs->push_back(mesh3d);
        if (DEBUG) std::cout << "#nmesh is " <<  meshs->size() << std::endl;   
    }
    for (pugi::xml_node node = nd.child("node"); node; node = node.next_sibling("node")){
        std::cout << "#" <<nd.attribute("id").value() << "child " << node.attribute("id").value() << std::endl;
        getMeshs_node(node,trafo,scene_meshs,meshs);
    }
    trafo = prev;
}
std::vector<mesh> getMeshs(std::string path){
    //need a more efficient parser that will get node-transformation-mesh
    std::map<std::string, mesh> scene_meshes;
    openvdb::math::Mat4d trafo;
    
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(path.c_str());
    std::cout << "#Load result: " << result.description() << " library_geometries " << doc.child("COLLADA").child("library_geometries") << std::endl;
    int nbGeom = 0;  
    std::vector<mesh> meshs;  
    for (pugi::xml_node geometry = doc.child("COLLADA").child("library_geometries").first_child(); geometry; geometry = geometry.next_sibling())
    {
        if (std::string(geometry.name()) != "geometry") continue;
        std::string idgeom(geometry.attribute("id").value());
        std::cout << "#geometry: " << idgeom << " " << nbGeom << std::endl;
        mesh mesh3d;   
        //find source geom
        pugi::xml_node meshnode = geometry.child("mesh");
        if (DEBUG)std::cout << "#meshnode " << meshnode << std::endl;
        pugi::xml_node vnode = meshnode.child("vertices");
        if (DEBUG)std::cout << "#vnode " << vnode << " x " << vnode.child("input") << std::endl;
        //get the ID of the array of float ie input source
        std::string idsource(vnode.child("input").attribute("source").value());
        //remove the #
        //std::replace(idsource.begin(), idsource.end(), '#', '');
        idsource.erase(std::remove(idsource.begin(), idsource.end(), '#'), idsource.end()); 
        if (DEBUG)std::cout << "#idsource " << idsource  << std::endl;
        pugi::xml_node varray;
        int vcount=0;
        //find_child_by_attribute(const char_t* name, const char_t* attr_name, const char_t* attr_value) const
        for (pugi::xml_node source = meshnode.child("source"); source; source = source.next_sibling("source")){
            std::string strid(source.attribute("id").value());
            if (DEBUG)std::cout << "#id " << strid  << std::endl;
            if ( strid == idsource){
                 varray =  source.child("float_array");  
                 vcount =  varray.attribute("count").as_int();
                 break;
            }
        }
        if (DEBUG)std::cout << "#vcount " << vcount  << std::endl;
        //std::string varraysource(varray.value());
        //std::cout << "#varray " << varray.text().get()  << std::endl;
        std::string strarray( varray.text().get() );
        std::vector<float> pos = getBox(strarray);     
        std::vector<openvdb::Vec3s> vertices = getPositionsS(pos);
        //now the face
        std::vector<openvdb:: Vec3I > faces;
        std::vector<openvdb:: Vec4I > quads;
        pugi::xml_node fnode = meshnode.child("polylist");//or triangles
        if (fnode) {
            //std::cout << "#fnode " << fnode.child("p").text().get()  << std::endl;
            //get the vertices fload array ID and parse it
            std::string strintfcount( fnode.child("vcount").text().get() );
            std::vector<int> fcount = getInts(strintfcount);
            std::string strintarray( fnode.child("p").text().get() );
            std::vector<int> f = getInts(strintarray);
            if (DEBUG)std::cout << "#fcount polylist " << fcount.size() << std::endl;
            for (std::vector<int>::size_type i = 0; i< fcount.size();i++){
                  if (fcount[i] == 3){
                    faces.push_back(openvdb:: Vec3I(f[i],f[i+1],f[i+2]));                
                    }
                  else if (fcount[i] == 4){
                    quads.push_back(openvdb:: Vec4I(f[i],f[i+1],f[i+2],f[i+3]));  
                    }       
            }
        }
        else {
            fnode = meshnode.child("triangles");//or triangles
            std::string strintarray( fnode.child("p").text().get() );
            std::vector<int> f = getInts(strintarray);
            //actually need every third position not all
            int counter = 0 ;
            int face[3];
            for (std::vector<int>::size_type i = 0; i< f.size();i++){
                 if ((i % 3) == 0)  {
                    face[counter] = f[i];
                    counter++;
                    if (counter == 3) {
                        counter = 0;
                         openvdb:: Vec3I  p(face[0],face[1],face[2]);//, openvdb::util::INVALID_IDX);
                         faces.push_back(p);
                    }
                }
            }    
            //faces = getPositionsInt(f);
            quads.resize(0);
            if (DEBUG)std::cout << "#fcount triangle " <<faces.size() << " " << f.size() << std::endl;
        }
        
        //same for the triangle
        mesh3d.vertices = vertices; 
        mesh3d.faces = faces;
        mesh3d.quads = quads;
        nbGeom+=1; 
        //meshs.push_back(mesh3d); 
        scene_meshes[idgeom] = mesh3d;
    }

    //we have the mesh now create the actual geom with a recursive function that will 
    //apply the transformation matrices
    pugi::xml_node scene_node = doc.child("COLLADA").child("library_visual_scenes").first_child();
    //go through node until finding instance geom
    trafo.setIdentity();
    getMeshs_node(scene_node,trafo,scene_meshes,&meshs);
    return meshs;
}

mesh getMesh(std::string path){
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(path.c_str());
    std::cout << "#Load result: " << result.description() << " library_geometries " << doc.child("COLLADA").child("library_geometries") << std::endl;
    mesh mesh3d;    
    //whatabout different geoms in one file
    //how many geometry node there is ?
    pugi::xml_node gnodes = doc.child("COLLADA").child("library_geometries").child("geometry"); //is this gave the number of geom ?    
    std::cout << "#gnodes " << gnodes << std::endl;//
    pugi::xml_node meshnode = doc.child("COLLADA").child("library_geometries").child("geometry").child("mesh");
    std::cout << "#meshnode " << meshnode << std::endl;
    pugi::xml_node vnode = meshnode.child("vertices");
    std::cout << "#vnode " << vnode << " x " << vnode.child("input") << std::endl;
    //get the ID of the array of float ie input source
    std::string idsource(vnode.child("input").attribute("source").value());
    //remove the #
    //std::replace(idsource.begin(), idsource.end(), '#', '');
    idsource.erase(std::remove(idsource.begin(), idsource.end(), '#'), idsource.end()); 
    std::cout << "#idsource " << idsource  << std::endl;
    pugi::xml_node varray;
    int vcount=0;
    //find_child_by_attribute(const char_t* name, const char_t* attr_name, const char_t* attr_value) const
    for (pugi::xml_node source = meshnode.child("source"); source; source = source.next_sibling("source")){
        std::string strid(source.attribute("id").value());
        std::cout << "#id " << strid  << std::endl;
        if ( strid == idsource){
             varray =  source.child("float_array");  
             vcount =  varray.attribute("count").as_int();
             break;
        }
    }
    std::cout << "#vcount " << vcount  << std::endl;
    //std::string varraysource(varray.value());
    //std::cout << "#varray " << varray.text().get()  << std::endl;
    std::string strarray( varray.text().get() );
    std::vector<float> pos = getBox(strarray);     
    std::vector<openvdb::Vec3s> vertices = getPositionsS(pos);
    //now the face
    std::vector<openvdb:: Vec3I > faces;
    std::vector<openvdb:: Vec4I > quads;
    pugi::xml_node fnode = meshnode.child("polylist");//or triangles
    if (fnode) {
        //std::cout << "#fnode " << fnode.child("p").text().get()  << std::endl;
        //get the vertices fload array ID and parse it
        std::string strintfcount( fnode.child("vcount").text().get() );
        std::vector<int> fcount = getInts(strintfcount);
        std::string strintarray( fnode.child("p").text().get() );
        std::vector<int> f = getInts(strintarray);
        for (std::vector<int>::size_type i = 0; i< fcount.size();i++){
              if (fcount[i] == 3){
                faces.push_back(openvdb:: Vec3I(f[i],f[i+1],f[i+2]));                
                }
              else if (fcount[i] == 4){
                quads.push_back(openvdb:: Vec4I(f[i],f[i+1],f[i+2],f[i+3]));  
                }       
        }
        //faces = getPositionsInt(f);
    }
    else {
        fnode = meshnode.child("triangles");//or triangles
        std::string strintarray( fnode.child("p").text().get() );
        std::vector<int> f = getInts(strintarray);
        faces = getPositionsInt(f);
        quads.resize(0);
    }
    //same for the triangle
    mesh3d.vertices = vertices; 
    mesh3d.faces = faces;
    mesh3d.quads = quads;
    //grid work
    //openvdb::FloatGrid::Ptr grid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
    //    *openvdb::math::Transform::createLinearTransform(), vertices, faces);
    //OR which work too
    //openvdb::tools::internal::MeshVoxelizer<openvdb::FloatTree>
    //    voxelizer(vertices, faces);
    //voxelizer.runParallel();
    return mesh3d;
}


//we load the autopack xml setup and create the grid class object
big_grid load_xml(std::string path,int _mode,unsigned _seed){
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(path.c_str());
    //std::cout << path << std::endl;
    std::cout << "#Load result: " << result.description() << ", AutoFillSetup name: " << doc.child("AutoFillSetup").attribute("name").value() << std::endl;
    //get the option
    const double smallestObjectSize = atof(doc.child("AutoFillSetup").child("options").attribute("smallestProteinSize").value());
    const unsigned seed = _seed;
    const int mode = _mode;
    stepsize = smallestObjectSize * 1.1547;
    std::cout <<"#step " << stepsize <<std::endl;
    std::cout <<"#seed " <<seed<<std::endl;
    std::string strbb(doc.child("AutoFillSetup").child("options").attribute("boundingBox").value());
    std::vector<float> bb = getBox(strbb);

    bool pickWeightedIngr =doc.child("AutoFillSetup").child("options").attribute("pickWeightedIngr").as_bool();
    bool pickRandPt =doc.child("AutoFillSetup").child("options").attribute("pickRandPt").as_bool();

    //only cytoplsme for now, should parse, organelle and gradient as well
    //could we use openvdb do compute compute/prepare the gradient?
    std::vector<sphere> _ingredients;
    //need to add organelle as well...organelle ingredient are inside organelle levelSet.
    pugi::xml_node cytoplasme = doc.child("AutoFillSetup").child("cytoplasme");
    //float radius, int mode, float concentration, 
    //     float packingPriority,int nbMol,std::string name, point color    
    for (pugi::xml_node ingredient = cytoplasme.child("ingredient"); ingredient; ingredient = ingredient.next_sibling("ingredient")){
        if (DEBUG) std::cout << "#ingredient " << ingredient.attribute("name").value() << std::endl;
        std::string str(ingredient.attribute("radii").value()); 
        std::vector<float> radii = getBox(str);     
        std::string posstr(ingredient.attribute("positions").value()); 
        std::vector<float> pos = getBox(posstr);     
        std::vector<openvdb::Vec3f> positions = getPositions(pos);
        //float r = getRadii(str); //different radii... as well as the position...  
        //std::cout << r << std::endl;
        float mol = ingredient.attribute("molarity").as_float();
        if (DEBUG)std::cout << "#molarity "<<mol << std::endl;
        float priority = ingredient.attribute("packingPriority").as_float();
        if (DEBUG)std::cout << "#priority "<< priority << std::endl;        
        int nMol = ingredient.attribute("nbMol").as_int();
        if (DEBUG)std::cout << "#nmol "<< nMol << std::endl; 
        std::string iname(ingredient.attribute("name").value()); 
        openvdb::Vec3f color(1,0,0);
        if (ingredient.attribute("color")){        
            std::string strcol(ingredient.attribute("color").value());
            openvdb::Vec3f color = getArray(strcol);
        }
        unsigned nbJitter = ingredient.attribute("nbJitter").as_int();
        if (DEBUG)std::cout << "#color "<< color << std::endl;
        std::string strjitter(ingredient.attribute("jitterMax").value());
        openvdb::Vec3f jitter =  getArray(strjitter);
        if (DEBUG)std::cout << "#jitter "<< jitter << std::endl;     
        std::string straxe(ingredient.attribute("principalVector").value());
        openvdb::Vec3f principalVector =  getArray(straxe);
        if (DEBUG)std::cout << "#principalVector "<< principalVector << std::endl;
        bool fSphere=forceSphere;
        if (ingredient.attribute("forceSphere"))
            fSphere = ingredient.attribute("forceSphere").as_bool(); 
        //also need packingMode,perturbAxisAmplitude
        //mesh file ... should use multiSphere function or makeMesh  
        //sphere ingr = makeSphere(r,_mode,mol,priority,nMol,iname,
        //                        color,nbJitter,jitter);
        //get the meshFile and get vertices+faces
        std::string strmeshFile(ingredient.attribute("meshFile").value());
        //if dae can use the xmlparser
        std::cout << "# mesh file ? " << strmeshFile << " x " << strmeshFile.empty() << std::endl;
        //can be none
        mesh mesh3d;
        std::vector<mesh> meshs;
        sphere ingr ;
        if ((!strmeshFile.empty())&&(!forceSphere)){
            //mesh3d = getMesh(strmeshFile);
            //meshs = getMeshs_assimp(strmeshFile);
            meshs = getMeshs(strmeshFile);
            std::cout << "# mesh nb " << meshs.size()  << std::endl;
            if (meshs.size()==1) {
                ingr = makeMeshIngredient(radii,_mode,mol,priority,nMol,iname,
                    color,nbJitter,jitter,meshs[0]);
            }
            else {
                ingr = makeMeshesIngredient(radii,_mode,mol,priority,nMol,iname,
                    color,nbJitter,jitter,meshs);
            }
            
        }
        else {
            ingr = makeMultiSpheres(radii,_mode,mol,priority,nMol,iname,
                    color,nbJitter,jitter,positions);
        }
        ingr.filename = strmeshFile;
        ingr.principalVector=principalVector;
        ingr.useRotAxis = false;
        if (ingredient.attribute("useRotAxis")){
            ingr.useRotAxis = ingredient.attribute("useRotAxis").as_bool();
            ingr.rotRange = ingredient.attribute("rotRange").as_float();
            if (ingr.useRotAxis){
                std::string strRaxe(ingredient.attribute("rotAxis").value());
                ingr.rotAxis =  getArray(strRaxe);
            }
        }
        ingr.perturbAxisAmplitude = 0.1f;
        if (ingredient.attribute("perturbAxisAmplitude")){
            ingr.perturbAxisAmplitude = ingredient.attribute("perturbAxisAmplitude").as_float();
        }
        //packing mode overwrite from xml file
        ingr.packingMode = std::string(ingredient.attribute("packingMode").value());
        _ingredients.push_back(ingr);
        float volume= 1;
        std::cout << "#ingredient done!" << std::endl;
    }
    //make the grid and return it
    const openvdb::Vec3d ibotleft(bb[0],bb[1],bb[2]);
    const openvdb::Vec3d iupright(bb[3],bb[4],bb[5]);    
    if (DEBUG)std::cout << "#box "<<ibotleft<<" "<<iupright << std::endl;

    //should create the grid from a xml file...easier setup
    big_grid g(stepsize,ibotleft,iupright,seed);
    g.vRangeStart=0.0;
    g.pickWeightedIngr=pickWeightedIngr;
    g.pickRandPt =pickRandPt; 
    //g.setMode(0);//1-close packing 0-random
    g.setIngredients(_ingredients);
    return g; 
}