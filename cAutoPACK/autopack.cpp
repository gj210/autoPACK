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
"""
*/

/*
openvdb includes
*/
#include <openvdb/openvdb.h>
#include <openvdb/Types.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tree/Tree.h>
#include <openvdb/tools/GridTransformer.h>

/*
general include
*/
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <random>

/*
xml parser for collada import, dont forgot the lib in the folder
*/
#include "../../cautoPACK/pugixml-1.2/src/pugixml.hpp" 


/*some general option for autpoack to run*/

float stepsize = 15*1.1547;         //grid step size ie smallest ingredients radius
const float dmax = 999999.0;        //maximum value assign to the grid for intialisation
const int rejectionThresholdIngredient = 300;//30 by default, number of rejection before stoppin a ingredient
const bool DEBUG = false;                   //DEBUG mode for prining some information
const float MaxRadius = (2.0*1.1547)*20.0;  //Maximum radius allowed, this is used for the voxelisation
const float spherewidth = 10.0;             //default or more ? close packing need to increase this number
const std::string packing = "random";       //packing mode, can be random or close
bool forceSphere = true;
//deprecated use now the packingMode from the setup file
//easy to overwrite but need to recompile

//a simple point, deprecated with openvdb::Vec3
struct point {
    float x,y,z;
};
/* 
a box, let define the grid size in x y an z
again deprectaed with openvdb
*/
struct box {
    unsigned x,y,z;
};

/* 
a simple mesh struct
need to add normal (surface packing and organelle definition)
- faceNormal
- verticeNormal ?
*/
struct mesh {
    std::vector<openvdb::Vec3s> vertices;
    std::vector<openvdb::Vec3I > faces;
    std::vector<openvdb::Vec4I > quads;
};

// Define a local function that retrieves a and b values from a CombineArgs
// struct and then sets the result member to the maximum of a and b.
// what if we want to keep the level set ? this apply to all value
struct Local {
    static inline void min(openvdb::CombineArgs<float>& args) {
        if (args.b() < args.a()) {
            if ((args.b() < 0.0 )||args.bIsActive()){
                // Transfer the B value and its active state if active
                args.setResult(args.b());
                args.setResultIsActive(false);//args.bIsActive()); set false == occupied
            }
            else {
                //also set B value but active is true == empty
                args.setResult(args.a());
                args.setResultIsActive(true);//args.aIsActive());
            }
        } else {
            // Preserve the A value and its active state.
            args.setResult(args.a());
            args.setResultIsActive(args.aIsActive());
        }
    }
    static inline void rmin(openvdb::CombineArgs<float>& args) {
        if (args.b() < args.a()) {
                // Transfer the B value and its active state 
                args.setResult(args.b());
                args.setResultIsActive(args.bIsActive());
        } else {
            // Preserve the A value and its active state.
            args.setResult(args.a());
            args.setResultIsActive(args.aIsActive());
        }
    }
};

/*
a container // ie organelle
not use yet
*/

struct container{
    std::string name;
    int id;
    openvdb::FloatGrid::Ptr grid;
    openvdb::CoordBBox bbox;
    //why no use the mesh struct.
    std::vector<openvdb::Vec3s> vertices;
    std::vector<openvdb::Vec3s> normals;
    std::vector<openvdb::Vec3I > faces;
    std::vector<openvdb::Vec4I > quads;
    //list of indice too ?
};

/* 
a SingleSphere ingredient with the radius and the packing mode
Need to  be rename to Ingredient
Should it bea class instead of a struct ?
*/
struct sphere {
    float radius;
    std::vector<float> radii;
    std::vector<openvdb::Vec3f> positions;
    float minRadius;
    float maxRadius;
    float stepsize;
    //point trans;
    int mode;
    int nbMol;
    float completion;
    float packingPriority;
    std::string name;
    std::string packingMode;
    int counter;
    int rejectionCounter;
    int rejectionThreshold;
    unsigned nbJitter;
    bool active;
    openvdb::Vec3f color;
    openvdb::Vec3f trans;
    openvdb::Vec3f jitterMax;
    openvdb::Vec3f principalVector;
    openvdb::FloatGrid::Ptr gsphere;
    openvdb::CoordBBox bbox;
    float miniVal;
    float maxiVal;
    mesh mesh3d;
    std::string filename;
    //rotation paramters
    bool useRotAxis;
    float rotRange;
    openvdb::Vec3f rotAxis;
    float perturbAxisAmplitude;
};

/* the comparison function are strict translatio from the python code */
//sort function for ingredient//
//The value returned indicates whether the element passed as first argument 
//is considered to go before the second in the specific strict weak ordering it defines.
//can weuse template here ? so ca accept any ingredient type..
bool ingredient_compare1(sphere* x, sphere* y){
    /*
    """
    sort ingredients using decreasing priority and decreasing radii for
    priority ties and decreasing completion for radii ties
    
    """
    */
    
    float p1 = x->packingPriority;
    float p2 = y->packingPriority;
    if (p1 < p2) //# p1 > p2
        return true;
    else if (p1==p2){ //# p1 == p1
       float r1 = x->minRadius;
       float r2 = y->minRadius;
       if (r1 > r2) //# r1 < r2
           return true;
       else if (r1==r2){ //# r1 == r2
           float c1 = x->completion;
           float c2 = y->completion;
           if (c1 >= c2) //# c1 > c2
               return true;
           else
               return false;
       }else
           return false;
    }else
       return false;
}

bool ingredient_compare0(sphere* x, sphere* y){
    /*
    """
    sort ingredients using decreasing priority and decreasing radii for
    priority ties and decreasing completion for radii ties
    
    """
    */
    
    float p1 = x->packingPriority;
    float p2 = y->packingPriority;
    if (p1 > p2) //# p1 > p2
        return true;
    else if (p1==p2){ //# p1 == p1
       float r1 = x->minRadius;
       float r2 = y->minRadius;
       if (r1 > r2) //# r1 < r2
           return true;
       else if (r1==r2){ //# r1 == r2
           float c1 = x->completion;
           float c2 = y->completion;
           if (c1 >= c2) //# c1 > c2
               return true;
           else
               return false;
       }else
           return false;
    }else
       return false;
}

bool ingredient_compare2(sphere* x, sphere* y){
    /*
    """
    sort ingredients using decreasing radii and decresing completion
    for radii matches
    
    """
    */
    
    float r1 = x->minRadius;
    float r2 = y->minRadius;
   if (r1 < r2) //# r1 < r2
       return true;
   else if (r1==r2){ //# r1 == r2
       float c1 = x->completion;
       float c2 = y->completion;
       if (c1 >= c2) //# c1 > c2
           return true;
       else
           return false;
   }
   else
    return false;
}



/*
openvdb::CoordBBox getBB(float radius,  openvdb::Vec3f pos){
    const openvdb::Vec3d ibotleft(pos-radius,pos-radius,pos-radius);//(0,0,0)
    const openvdb::Vec3d iupright(pos+radius,pos+radius,pos+radius);//(1000,1000,10);
    openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(stepsize);
    openvdb::Vec3d botleft=transform->worldToIndex(ibotleft);
    openvdb::Vec3d upright=transform->worldToIndex(iupright);
    openvdb::Coord left(openvdb::tools::local_util::roundVec3(botleft));//(openvdb::Int32)botleft.x(),(openvdb::Int32)botleft.y(),(openvdb::Int32)botleft.z());
    openvdb::Coord right(openvdb::tools::local_util::roundVec3(upright));//(openvdb::Int32)upright.x(),(openvdb::Int32)upright.y(),(openvdb::Int32)upright.z());
    //define the active region that will be our boundary. set to max everywhere
    openvdb::CoordBBox bbox = openvdb::CoordBBox(left,right);//min and max
    return bbox;
}
*/


//exapl of function applied on grid data
struct SetMaxToDefault {
    float _max;
    SetMaxToDefault(float max) {_max=max;}
    inline void operator()(const openvdb::FloatGrid::ValueAllIter& iter) const {
       if (iter.getValue() == _max) iter.setValue(dmax);
    }
};


//helper to create a singleSphere ingredient given a radius, and some options
inline sphere makeSphere(float radius, int mode, float concentration, 
         float packingPriority,int nbMol,std::string name, openvdb::Vec3f color,
        unsigned nbJitter,openvdb::Vec3f jitterMax){
    sphere sp;
    sp.radius = radius;
    sp.mode = mode;
    sp.minRadius = radius;
    sp.maxRadius = radius;
    sp.nbMol=nbMol;
    sp.completion = 0.0;
    sp.packingMode = packing;
    sp.packingPriority = packingPriority;
    sp.name = name;
    sp.counter = 0;
    sp.rejectionCounter = 0;
    sp.rejectionThreshold  = rejectionThresholdIngredient;
    sp.active = true;
    sp.color = color;
    sp.nbJitter = nbJitter;
    sp.trans = openvdb::Vec3f(0,0,0);
    sp.radii.resize(1);
    sp.radii[0] = radius;
    sp.positions.resize(1);
    sp.positions[0] = openvdb::Vec3f(0,0,0);
    sp.jitterMax = jitterMax;
    //build the grid
    float size=stepsize;
    if (stepsize > sp.minRadius) size = sp.minRadius/2.0;
    sp.stepsize = size;

    sp.gsphere = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
              /*radius=*/radius, /*center=*/sp.trans,//where is the sphere is  index or worl
               /*voxel size=*/size, /*width=*/25.0);//width could be maxRadius / stepsize (int)MaxRadius/stepsize
    //sp.gsphere->signedFloodFill();
        
    //sphere bounding box
    const openvdb::Vec3d ibotleft(-radius,-radius,-radius);//(0,0,0)
    const openvdb::Vec3d iupright(radius,radius,radius);//(1000,1000,10);
    openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(/*voxel size=*/stepsize);
    openvdb::Vec3d botleft=transform->worldToIndex(ibotleft);
    openvdb::Vec3d upright=transform->worldToIndex(iupright);
    
    openvdb::Coord left(openvdb::tools::local_util::roundVec3(botleft));//(openvdb::Int32)botleft.x(),(openvdb::Int32)botleft.y(),(openvdb::Int32)botleft.z());
    openvdb::Coord right(openvdb::tools::local_util::roundVec3(upright));//(openvdb::Int32)upright.x(),(openvdb::Int32)upright.y(),(openvdb::Int32)upright.z());
    //define the active region that will be our boundary. set to max everywhere
    sp.bbox = openvdb::CoordBBox(left,right);//min and max                
    sp.gsphere->evalMinMax(sp.miniVal,sp.maxiVal);
    // Apply the functor to all active values.
    //std::cout << "#mini " << sp.miniVal << " maxi " << sp.maxiVal << std::endl;
    //openvdb::tools::foreach(sp.gsphere->beginValueAll(), SetMaxToDefault(sp.maxiVal));
    //translate (const Coord &t)
    sp.useRotAxis=false;
    return sp;
}

//helper to create a multiSpheres ingredient given a list of radii and positions
//if only one radius and one position given we build a uniq sphere.
inline sphere makeMultiSpheres(std::vector<float> radii, int mode, float concentration, 
         float packingPriority,int nbMol,std::string name, openvdb::Vec3f color,
        unsigned nbJitter,openvdb::Vec3f jitterMax,std::vector<openvdb::Vec3f> positions){
    sphere sp;
    sp.radii = radii;
    sp.positions = positions;
    sp.radius = radii[0];
    sp.mode = mode;
    sp.minRadius = 9999999.0;
    sp.maxRadius = 0.0;
    for (int i =0 ; i < radii.size();i++){
        if (radii[i] < sp.minRadius)sp.minRadius = radii[i];
        if (radii[i] > sp.maxRadius)sp.maxRadius = radii[i];
    }
    if (DEBUG)std::cout << "#miniR " << sp.minRadius << " maxiR " << sp.maxRadius <<std::endl;
    sp.nbMol=nbMol;
    sp.completion = 0.0;
    sp.packingMode = packing;
    sp.packingPriority = packingPriority;
    sp.name = name;
    sp.counter = 0;
    sp.rejectionCounter = 0;
    sp.rejectionThreshold  = rejectionThresholdIngredient;
    sp.active = true;
    sp.color = color;
    sp.nbJitter = nbJitter;
    sp.trans = openvdb::Vec3f(0,0,0);
    sp.jitterMax = jitterMax;
    //build the grid
    //need to create as many grid as sphere, then combine then in one uniq ie union?

    float size=stepsize;
    if (stepsize > sp.minRadius) size = stepsize;//sp.minRadius/2.0;
    sp.stepsize = stepsize;
    
    if (radii.size() == 1 ){
         sp.gsphere = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
              /*radius=*/sp.radius, /*center=*/positions[0],//where is the sphere is  index or worl
               /*voxel size=*/size, /*width=*/spherewidth); 
        //sphere bounding box
        const openvdb::Vec3d ibotleft(-sp.radius,-sp.radius,-sp.radius);//(0,0,0)
        const openvdb::Vec3d iupright(sp.radius,sp.radius,sp.radius);//(1000,1000,10);
        openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(/*voxel size=*/stepsize);
        openvdb::Vec3d botleft=transform->worldToIndex(ibotleft);
        openvdb::Vec3d upright=transform->worldToIndex(iupright);
        
        openvdb::Coord left(openvdb::tools::local_util::roundVec3(botleft));//(openvdb::Int32)botleft.x(),(openvdb::Int32)botleft.y(),(openvdb::Int32)botleft.z());
        openvdb::Coord right(openvdb::tools::local_util::roundVec3(upright));//(openvdb::Int32)upright.x(),(openvdb::Int32)upright.y(),(openvdb::Int32)upright.z());
        //define the active region that will be our boundary. set to max everywhere
        sp.bbox = openvdb::CoordBBox(left,right);//min and max                
    }
    else {
        std::vector<openvdb::FloatGrid::Ptr> gspheres;
        sp.gsphere = openvdb::FloatGrid::create(dmax);
        sp.gsphere->setTransform(
            openvdb::math::Transform::createLinearTransform(/*voxel size=*/size)); 
        gspheres.resize(radii.size());
        for (int i =0 ; i < radii.size();i++){
            //is the position in xyz or ijk ?
            if (DEBUG)std::cout << "#r " << radii[i] << " pos " <<  positions[i] << std::endl;
            gspheres[i] = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
                  /*radius=*/radii[i], /*center=*/positions[i],//where is the sphere is  index or worl
                   /*voxel size=*/size, /*width=*/(int)spherewidth);
            //union with gsphere
            sp.gsphere->tree().combineExtended(gspheres[i]->tree(), Local::rmin);
            //openvdb::tools::csgUnion(sp.gsphere->tree(),gspheres[i]->tree());
        }
        sp.gsphere->prune();
        //sp.gsphere->signedFloodFill();
        //sphere bounding box
        //openvdb::Coord dimgrid = sp.gsphere->evalActiveVoxelDim();//ninjnk
        //sphere bounding box
        sp.bbox = sp.gsphere->evalActiveVoxelBoundingBox();
        sp.useRotAxis=false;
    }
    /*
    const openvdb::Vec3d ibotleft(-radius,-radius,-radius);//(0,0,0)
    const openvdb::Vec3d iupright(radius,radius,radius);//(1000,1000,10);
    openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(stepsize);
    openvdb::Vec3d botleft=transform->worldToIndex(ibotleft);
    openvdb::Vec3d upright=transform->worldToIndex(iupright);
    
    openvdb::Coord left(openvdb::tools::local_util::roundVec3(botleft));//(openvdb::Int32)botleft.x(),(openvdb::Int32)botleft.y(),(openvdb::Int32)botleft.z());
    openvdb::Coord right(openvdb::tools::local_util::roundVec3(upright));//(openvdb::Int32)upright.x(),(openvdb::Int32)upright.y(),(openvdb::Int32)upright.z());
    //define the active region that will be our boundary. set to max everywhere
    sp.bbox = openvdb::CoordBBox(left,right);//min and max 
    */               
    sp.gsphere->evalMinMax(sp.miniVal,sp.maxiVal);
    // Apply the functor to all active values.
    std::cout << "#mini " << sp.miniVal << " maxi " << sp.maxiVal << " bb " << sp.bbox<<std::endl;
    //openvdb::tools::foreach(sp.gsphere->beginValueAll(), SetMaxToDefault(sp.maxiVal));
    //translate (const Coord &t)
    return sp;
}

//helper to create an ingredient given a 3d mesh triangles or quads
inline sphere makeMeshIngredient(std::vector<float> radii, int mode, float concentration, 
         float packingPriority,int nbMol,std::string name, openvdb::Vec3f color,
        unsigned nbJitter,openvdb::Vec3f jitterMax, mesh mesh3d){
    sphere sp;
    sp.radii = radii;
    sp.positions.resize(1);
    sp.positions[0] = openvdb::Vec3f(0,0,0);
    //sp.positions = positions;
    sp.radius = radii[0];
    sp.mode = mode;
    sp.minRadius = sp.radius ;
    sp.maxRadius = 0.0;
    sp.nbMol=nbMol;
    sp.completion = 0.0;
    sp.packingMode = packing;
    sp.packingPriority = packingPriority;
    sp.name = name;
    sp.counter = 0;
    sp.rejectionCounter = 0;
    sp.rejectionThreshold  = rejectionThresholdIngredient;
    sp.active = true;
    sp.color = color;
    sp.nbJitter = nbJitter;
    sp.trans = openvdb::Vec3f(0,0,0);
    sp.jitterMax = jitterMax;
    //build the grid
    //need to create as many grid as sphere, then combine then in one uniq ie union?
    sp.gsphere = openvdb::FloatGrid::create(dmax);
    float size=stepsize;
    if (stepsize > sp.minRadius) size = stepsize;//sp.minRadius/2.0;
    sp.stepsize = stepsize;
    //levelSet
    if (DEBUG) std::cout << "#mesh have v "<< mesh3d.vertices.size() << " f " << mesh3d.faces.size() << " q " << mesh3d.quads.size() << std::endl;
    
    if ((mesh3d.quads.size() != 0)&&(mesh3d.faces.size() != 0)){
        sp.gsphere = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
            *openvdb::math::Transform::createLinearTransform(size), 
                mesh3d.vertices, mesh3d.faces, mesh3d.quads);
    }
    else if ((mesh3d.quads.size() != 0)&&(mesh3d.faces.size() == 0)){
        sp.gsphere = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
            *openvdb::math::Transform::createLinearTransform(size), 
            mesh3d.vertices, mesh3d.quads);
    }
    else if ((mesh3d.quads.size() == 0)&&(mesh3d.faces.size() != 0)){
        sp.gsphere = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
            *openvdb::math::Transform::createLinearTransform(size), 
            mesh3d.vertices, mesh3d.faces);
    }
    
    //sp.gsphere->signedFloodFill();
    //sp.gsphere = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(
    //    *openvdb::math::Transform::createLinearTransform(stepsize), mesh3d.vertices, mesh3d.faces,mesh3d.quads,(int)spherewidth,(int)spherewidth);
    //what about sdf ?
    //meshToSignedDistanceField (const openvdb::math::Transform &xform, const std::vector< Vec3s > &points, const std::vector< Vec3I > &triangles, const std::vector< Vec4I > &quads, float exBandWidth, float inBandWidth)
    //openvdb::tools::internal::MeshVoxelizer<openvdb::FloatTree>
    //   voxelizer(mesh3d.vertices, mesh3d.faces);
    //voxelizer.runParallel();
    //sp.gsphere = openvdb::FloatGrid::create(voxelizer);//use 4I ?
    sp.bbox = sp.gsphere->evalActiveVoxelBoundingBox();
    sp.gsphere->evalMinMax(sp.miniVal,sp.maxiVal);//distance
    // Apply the functor to all active values.
    sp.useRotAxis=false;
    std::cout << "#radius is "<< sp.minRadius <<" mini " << sp.miniVal << " maxi " << sp.maxiVal << " bb " << sp.bbox<<std::endl;
    //openvdb::tools::foreach(sp.gsphere->beginValueAll(), SetMaxToDefault(sp.maxiVal));
    //translate (const Coord &t)
    return sp;
}

//helper to create an ingredient given a different 3d mesh triangles or quads
inline sphere makeMeshesIngredient(std::vector<float> radii, int mode, float concentration, 
         float packingPriority,int nbMol,std::string name, openvdb::Vec3f color,
        unsigned nbJitter,openvdb::Vec3f jitterMax, std::vector<mesh> meshs){
    sphere sp;
    sp.radii = radii;
    sp.positions.resize(1);
    sp.positions[0] = openvdb::Vec3f(0,0,0);
    //sp.positions = positions;
    sp.radius = radii[0];
    sp.mode = mode;
    sp.minRadius = sp.radius ;
    sp.maxRadius = 0.0;
    sp.nbMol=nbMol;
    sp.completion = 0.0;
    sp.packingMode = packing;
    sp.packingPriority = packingPriority;
    sp.name = name;
    sp.counter = 0;
    sp.rejectionCounter = 0;
    sp.rejectionThreshold  = rejectionThresholdIngredient;
    sp.active = true;
    sp.color = color;
    sp.nbJitter = nbJitter;
    sp.trans = openvdb::Vec3f(0,0,0);
    sp.jitterMax = jitterMax;
    //build the grid
    //need to create as many grid as sphere, then combine then in one uniq ie union?
    sp.gsphere = openvdb::FloatGrid::create(dmax);
    float size=stepsize;
    if (stepsize > sp.minRadius) size = stepsize;//sp.minRadius/2.0;
    sp.stepsize = stepsize;
    //levelSet
    std::vector<openvdb::FloatGrid::Ptr> gspheres;
    sp.gsphere->setTransform(
            openvdb::math::Transform::createLinearTransform(/*voxel size=*/size));    
    gspheres.resize(meshs.size());
    if (DEBUG) std::cout << "#merge "<< meshs.size() << " voxelmesh " << std::endl;
    for (int i =0 ; i < meshs.size();i++){
        //is the position in xyz or ijk ?
        if (DEBUG) std::cout << "#mesh have v "<< meshs[i].vertices.size() << " f " << meshs[i].faces.size() << " q " << meshs[i].quads.size() << std::endl;
        if ((meshs[i].quads.size() != 0)&&(meshs[i].faces.size() != 0)){
            gspheres[i] = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
                *openvdb::math::Transform::createLinearTransform(size), 
                    meshs[i].vertices, meshs[i].faces, meshs[i].quads);
        }
        else if ((meshs[i].quads.size() != 0)&&(meshs[i].faces.size() == 0)){
            gspheres[i] = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
                *openvdb::math::Transform::createLinearTransform(size), 
                meshs[i].vertices, meshs[i].quads);
        }
        else if ((meshs[i].quads.size() == 0)&&(meshs[i].faces.size() != 0)){
            gspheres[i] = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(
                *openvdb::math::Transform::createLinearTransform(size), 
                meshs[i].vertices, meshs[i].faces);
        }
        if (DEBUG) std::cout <<  "combine " <<size <<std::endl;
        sp.gsphere->tree().combineExtended(gspheres[i]->tree(), Local::rmin);
        sp.gsphere->prune();
        //openvdb::tools::csgUnion(sp.gsphere->tree(),gspheres[i]->tree());
    }
    if (DEBUG) std::cout <<  "OK merged all of them and step is "<< size <<std::endl;
    //sp.gsphere->prune();
    sp.bbox = sp.gsphere->evalActiveVoxelBoundingBox();
    if (DEBUG) std::cout <<  "OK bbox "<< sp.bbox << std::endl;    
    sp.useRotAxis=false;
    sp.gsphere->evalMinMax(sp.miniVal,sp.maxiVal);//distance
    // Apply the functor to all active values.
    std::cout << "#radius is "<< sp.minRadius <<" mini " << sp.miniVal << " maxi " << sp.maxiVal << " bb " << sp.bbox<<std::endl;
    //openvdb::tools::foreach(sp.gsphere->beginValueAll(), SetMaxToDefault(sp.maxiVal));
    //translate (const Coord &t)
    return sp;
}

//code from the openvdb documentation
template<typename OpType>
void processTypedGrid(openvdb::GridBase::Ptr grid, OpType& op)
{
#define CALL_OP(GridType) \
    op.template operator()<GridType>(openvdb::gridPtrCast<GridType>(grid))

    if (grid->isType<openvdb::BoolGrid>())        CALL_OP(openvdb::BoolGrid);
    else if (grid->isType<openvdb::FloatGrid>())  CALL_OP(openvdb::FloatGrid);
    else if (grid->isType<openvdb::DoubleGrid>()) CALL_OP(openvdb::DoubleGrid);
    else if (grid->isType<openvdb::Int32Grid>())  CALL_OP(openvdb::Int32Grid);
    else if (grid->isType<openvdb::Int64Grid>())  CALL_OP(openvdb::Int64Grid);
    else if (grid->isType<openvdb::Vec3IGrid>())  CALL_OP(openvdb::Vec3IGrid);
    else if (grid->isType<openvdb::Vec3SGrid>())  CALL_OP(openvdb::Vec3SGrid);
    else if (grid->isType<openvdb::Vec3DGrid>())  CALL_OP(openvdb::Vec3DGrid);
    else if (grid->isType<openvdb::StringGrid>()) CALL_OP(openvdb::StringGrid);

#undef CALL_OP
}

// Define a functor that prunes the trees of grids of arbitrary type
// with a fixed pruning tolerance.
struct PruneOp {
    double tolerance;
    PruneOp(double t): tolerance(t) {}

    template<typename GridType>
    void operator()(typename GridType::Ptr grid) const
    {
        grid->tree().prune(typename GridType::ValueType(tolerance));
    }
    // Overload to handle string-valued grids
    void operator()(openvdb::StringGrid::Ptr grid) const
    {
        grid->tree().prune();
    }
};

//some expermental functions for the IJK<->U grid index convertion
inline unsigned getU(openvdb::Coord dim,openvdb::Coord ijk){
    return (int)(ijk.x()*dim.x()*dim.y() + ijk.y()*dim.x() + ijk.z());
}
inline void getIJK(int u,openvdb::Coord dim,int* i_ijk){
    // = {0,0,0};    
    //openvdb::Coord ijk(0,0,0);
    int nxnynz = dim.x()*dim.y()*dim.z();
    int nynz = dim.z()*dim.y();
    int nx = dim.x();
    int ny = dim.y();
    int nz = dim.z();
    int integer;
    float decimal;
    float fraction;
    int nysc;
    if (u < dim.z()){
        i_ijk[2] = u;
    }
    else if ((u < nynz)&&(u >= nxnynz)){
        //whats z
        fraction = (float)u/(float)dim.z();
        integer = (int) fraction;
        decimal = fraction - integer;
        i_ijk[2] = (int) round(decimal*dim.z());
        //whast y 
        i_ijk[1] = integer;  
    }
    else if ((u < nxnynz)&&(u >= nynz)){
        fraction = (float)u/(float)nynz;
        integer = (int) fraction;
        decimal = fraction - integer;
        nysc = ny * integer;
        //whast x 
        i_ijk[0] = integer;  
        fraction = (float)u/(float)nz;
        integer = (int) fraction;
        decimal = fraction - integer;
        //whats z        
        i_ijk[2] = (int) round(decimal*(float)nz);
        //whast y 
        //46867 
        //233 15477 201 77 603 77.7231
        //std::cout << integer << " " << nysc << " " << ny << " " << (int)((float)u/(float)nynz) << " " << nynz << " " << (float)u/(float)nynz << std::endl;
        i_ijk[1] = integer - (ny*(int)((float)u/(float)nynz));  
        //int (integer - (ny*int(float(u)/float(nynz))));
    }    
}

inline openvdb::Coord getIJKc(int u,openvdb::Coord dim){
    // = {0,0,0};    
    int i_ijk[3];
    i_ijk[0]=0;
    i_ijk[1]=0;
    i_ijk[2]=0;    
    //openvdb::Coord ijk(0,0,0);
    int nxnynz = dim.x()*dim.y()*dim.z();
    int nynz = dim.z()*dim.y();
    int nx = dim.x();
    int ny = dim.y();
    int nz = dim.z();
    int integer;
    float decimal;
    float fraction;
    int nysc;
    if (u < dim.z()){
        i_ijk[2] = u;
    }
    else if ((u < nynz)&&(u >= nxnynz)){
        //whats z
        fraction = (float)u/(float)dim.z();
        integer = (int) fraction;
        decimal = fraction - integer;
        i_ijk[2] = (int) round(decimal*dim.z());
        //whast y 
        i_ijk[1] = integer;  
    }
    else if ((u < nxnynz)&&(u >= nynz)){
        fraction = (float)u/(float)nynz;
        integer = (int) fraction;
        decimal = fraction - integer;
        nysc = ny * integer;
        //whast x 
        i_ijk[0] = integer;  
        fraction = (float)u/(float)nz;
        integer = (int) fraction;
        decimal = fraction - integer;
        //whats z        
        i_ijk[2] = (int) round(decimal*(float)nz);
        //whast y 
        //46867 
        //233 15477 201 77 603 77.7231
        //std::cout << integer << " " << nysc << " " << ny << " " << (int)((float)u/(float)nynz) << " " << nynz << " " << (float)u/(float)nynz << std::endl;
        i_ijk[1] = integer - (ny*(int)((float)u/(float)nynz));  
        //int (integer - (ny*int(float(u)/float(nynz))));
    } 
    return openvdb::Coord(i_ijk[0],i_ijk[1],i_ijk[2]);   
}


/* the main class that handle the packing
came from oleg trott code for the bidirectional array swapping.
extended to do the main autopack loop with ingredient and point picking
*/
struct big_grid { // needs 8*n bytes 
    //Lot of theses parameteres are deprecated and not used.
    //original wrote in hw.cc file
    std::vector<unsigned> all; //all point indice // should point to ijk
    std::vector<unsigned> empty;//available point indice 
    std::vector<point> data;    //the grid 3d coordintates
    std::vector<sphere> ingredients;//the list of sphere ingredient to pack
    std::vector<sphere*> activeIngr;
    std::vector<sphere*> activeIngr0;
    std::vector<sphere*> activeIngr12;
    std::vector<sphere*> droped_ingredient;
    std::vector<float> normalizedPriorities0;
    std::vector<float> normalizedPriorities;
    std::vector<float> thresholdPriorities; 
    std::vector<openvdb::Vec3f> rtrans;    //the grid 3d coordintates
    std::vector<openvdb::math::Mat4d> rrot;
    //float* distance; 
    std::vector<float> distance;    //the array of closest distances for each point
    unsigned num_empty;         //the number of free point available
    unsigned num_points;        //total number of point in the grid
    float space;                //spacing of the grid (unit depends on user)
    float maxradius;            //the biggest ingredient
    float minradius;            //the biggest ingredient
    float lowestPriority;       //priority the lowest after sorting
    float totalRadii;           //sum of all ingredient radii
    float totalPriorities;
    float vRangeStart;
    box nxyz;                   //the grid size in x, y and z
    point boundingBox0;         //the grid lower left corner coordinate
    point boundingBox1;         //the grid top right coordinate
    int mode;                   //the packing mode random or distance
    bool pickWeightedIngr = true;
    bool pickRandPt = true;
    bool use_gradient;
    int numActiveIngr;

    std::default_random_engine generator;
    std::uniform_real_distribution<float> uniform;// (0.0,1.0);
    std::normal_distribution<float> gauss;//(0.0,0.3);

    openvdb::FloatGrid::Ptr distance_grid;
    openvdb::CoordBBox bbox;
    openvdb::Coord dim;
    //openvdb::FloatGrid::Accessor accessor_distance;
    std::vector<openvdb::Coord> visited_rejected_coord;
    std::map<int, sphere*> results; 
    //the constructor that take as input the sizenor of the grid, the step, and the bouding box
    big_grid(unsigned nx, unsigned ny,unsigned nz, 
                float step, openvdb::Vec3d bot, openvdb::Vec3d up,float seed) : 
        all(nx*ny*nz),
        empty(nx*ny*nz),
        num_empty(nx*ny*nz),
        data(nx*ny*nz),
        distance(nx*ny*nz),
        uniform(0.0,1.0),
        gauss(0.0,0.3){
        
        generator.seed (seed);
        num_points = nx*ny*nz;      //initialize total number of point

        mode = 1;
        nxyz.x = nx;//+encapsulatingGrid
        nxyz.y = ny;//+encapsulatingGrid
        nxyz.z = nz;//+encapsulatingGrid 
        space = step;
        float xl=bot.x();
        float yl=bot.y();
        float zl=bot.z();
        
        distance_grid = openvdb::FloatGrid::create();//the indice ?
        distance_grid->setTransform(
            openvdb::math::Transform::createLinearTransform(/*voxel size=*/stepsize));
        //set active within the given bounding box
        // Define a coordinate with large signed indices.
        const openvdb::Vec3d ibotleft = bot;//(0,0,0)
        const openvdb::Vec3d iupright = up;//(1000,1000,10);
        openvdb::Vec3d botleft=distance_grid->worldToIndex(ibotleft);
        openvdb::Vec3d upright=distance_grid->worldToIndex(iupright);
        
        openvdb::Coord left(openvdb::tools::local_util::roundVec3(botleft));//(openvdb::Int32)botleft.x(),(openvdb::Int32)botleft.y(),(openvdb::Int32)botleft.z());
        openvdb::Coord right(openvdb::tools::local_util::roundVec3(upright));//(openvdb::Int32)upright.x(),(openvdb::Int32)upright.y(),(openvdb::Int32)upright.z());
        //define the active region that will be our boundary. set to max everywhere
        bbox = openvdb::CoordBBox(left,right);//min and max
        distance_grid->fill(bbox,dmax,true);//bbox, value, active
        distance_grid->prune();
        dim = distance_grid->evalActiveVoxelDim();
        
        num_points = dim.x()*dim.y()*dim.z();
        openvdb::FloatGrid::Accessor accessor_distance = distance_grid->getAccessor();
        std::cout << "#Testing distance access:" << std::endl;
        std::cout << "#Grid " << left << " "<< botleft << " = " << accessor_distance.getValue(left) << std::endl;
        std::cout << "#Grid " << right << " "<< upright << " = " << accessor_distance.getValue(right) << std::endl;
        std::cout << "#Grid " << bbox << std::endl;
        std::cout << "#Grid Npoints " << dim << " " << num_points << " "  << distance_grid->activeVoxelCount() << std::endl;
        unsigned i=0;           //counter for u indices             
        for(unsigned i = 0; i < num_points; ++i) { 
                all[i] = i;
                empty[i] = i;
            }//distance = new float[nx*ny*nz];
    }

    void setMode(unsigned _mode){
        //set the packing mode, random or distance
        mode = _mode;
    }

    void setIngredients(std::vector<sphere> _ingredients){
        //set the ingredients list to pack in the grid
        //retrieve the biggest one
        ingredients = _ingredients;
        maxradius = 0.0;
        minradius = 9999999.9;
        activeIngr.resize(ingredients.size());
        for(unsigned i = 0; i < ingredients.size(); ++i) { 
            if (ingredients[i].maxRadius > maxradius) 
                maxradius = ingredients[i].maxRadius;
            if (ingredients[i].minRadius < minradius) 
                minradius = ingredients[i].minRadius;
            activeIngr[i] = &ingredients[i];
        }
        numActiveIngr = ingredients.size();
        std::cout << "#min radius " << minradius << " stepsize " << stepsize << std::endl;
        //activeIngr = &ingredients;
    }
    
    void getMaxRadius(){
        maxradius = 0.0;
        //minradius = 9999.9;
        for(unsigned i = 0; i < numActiveIngr; ++i) { 
            if (activeIngr[i]->maxRadius > maxradius) 
                maxradius = activeIngr[i]->maxRadius;
            //if (activeIngr[i].minRadius < minradius) 
            //    minradius = activeIngr[i].minRadius;
        }
    }

    void getSortedActiveIngredients(){
        //pirorities120 is ot used ??
        std::vector<sphere*> ingr1;  // given priorities pass ptr ?
        std::vector<float> priorities1;
        std::vector<sphere*> ingr2;  // priority = 0 or none and will be assigned based on complexity
        std::vector<float> priorities2;
        std::vector<sphere*> ingr0;  // negative values will pack first in order of abs[packingPriority]
        std::vector<float> priorities0;
        sphere* lowestIng; 
        sphere* ing; 
        float r=0.0;
        float np=0.0;      
        for(unsigned i = 0; i < ingredients.size(); ++i) { 
            ing = &ingredients[i];
            if (ing->completion >= 1.0) continue;// # ignore completed ingredients
            if (ing->packingPriority == 0.0){
                ingr2.push_back(ing);
                priorities2.push_back(ing->packingPriority);
            }
            else if (ing->packingPriority > 0.0 ){
                ingr1.push_back(ing);
                priorities1.push_back(ing->packingPriority);
            }else{
                ingr0.push_back(ing);
                priorities0.push_back(ing->packingPriority);
            }
        }
        if (pickWeightedIngr){                      
            std::sort(ingr1.begin(),ingr1.end(),ingredient_compare1);
            std::sort(ingr2.begin(),ingr2.end(),ingredient_compare2);
            std::sort(ingr0.begin(),ingr0.end(),ingredient_compare0);
        }
        //#for ing in ingr3 : ing.packingPriority    = -ing.packingPriority
        //#GrahamAdded this stuff in summer 2011, beware!
        if (ingr1.size() != 0){
            lowestIng = ingr1[ingr1.size()-1];
            lowestPriority = lowestIng->packingPriority;
        }else 
            lowestPriority = 1.0;

        totalRadii = 0.0;
        for(unsigned i = 0; i < ingr2.size(); ++i) {
            ing = ingr2[i];
            r = ing->minRadius;
            totalRadii = totalRadii + r;
            if (r==0.0) { 
                //safety 
                totalRadii = totalRadii + 1.0;
             }
        }
        for(unsigned i = 0; i < ingr2.size(); ++i) {
            ing = ingr2[i];
            r = ing->minRadius;
            np = float(r)/float(totalRadii) * lowestPriority;
            //std::cout << "#packingPriority " << np << ' ' << r <<' '<< totalRadii<< ' ' << lowestPriority << std::endl;
            normalizedPriorities0.push_back(np);
            ing->packingPriority = np;
        }           
        
        activeIngr0 = ingr0;//#+ingr1+ingr2  #cropped to 0 on 7/20/10
        activeIngr12 = ingr1;//+ingr2
        activeIngr12.insert(activeIngr12.end(), ingr2.begin(), ingr2.end());
        //packingPriorities = priorities0;//+priorities1+priorities2
    }

    int prepareIngredient(){
        float pp;
        float np;
        sphere* ingr;
        getSortedActiveIngredients();
        std::cout << "#len(allIngredients) " << ingredients.size() << std::endl;
        std::cout << "#len(activeIngr0) " << activeIngr0.size() << std::endl;
        std::cout << "#len(activeIngr12) " << activeIngr12.size() << std::endl;
        //self.activeIngre_saved = self.activeIngr[:]
        totalPriorities = 0.0;// # 0.00001
        for(unsigned i = 0; i < activeIngr12.size(); ++i) {
            ingr = activeIngr12[i];
        //for priors in self.activeIngr12:
            pp = ingr->packingPriority;
            totalPriorities = totalPriorities + pp;
            std::cout << "#totalPriorities " << totalPriorities << ' ' << ingr->packingPriority <<std::endl;
        }
        
        
        float previousThresh = 0.0;
        //# Graham- Once negatives are used, if picked random# 
        //# is below a number in this list, that item becomes 
        //# the active ingredient in the while loop below
        for(unsigned i = 0; i < activeIngr0.size(); ++i) {
        //for priors in self.activeIngr0:
            normalizedPriorities.push_back(0.0);
            if (pickWeightedIngr)// #why ?
                thresholdPriorities.push_back(2);
        }        
        for(unsigned i = 0; i < activeIngr12.size(); ++i) {
        //for priors in self.activeIngr12:
            //#pp1 = 0
            pp = activeIngr12[i]->packingPriority;
            if (totalPriorities != 0)
                np = float(pp)/float(totalPriorities);
            else
                np=0.0;
            normalizedPriorities.push_back(np);
            std::cout << "#np is "<< np << " pp is "<< pp << " tp is " << np + previousThresh << std::endl;
            thresholdPriorities.push_back(np + previousThresh);
            previousThresh = np + float(previousThresh);
        }
        activeIngr = activeIngr0;// + self.activeIngr12
        activeIngr.insert(activeIngr.end(), activeIngr12.begin(), activeIngr12.end());

        int nls=0;
        int totalNumMols = 0;
        if (thresholdPriorities.size() == 0){
            for(unsigned i = 0; i < ingredients.size(); ++i) {
                totalNumMols += ingredients[i].nbMol;
                std::cout << "#nmol Fill5if is  for ingredient : "<< ingredients[i].name<< ' '  << ingredients[i].nbMol<< std::endl ;
            }
            std::cout << "#totalNumMols Fill5if = " << totalNumMols << std::endl;
        }else {                
            for(unsigned i = 0; i < thresholdPriorities.size(); ++i) {
                float threshProb = thresholdPriorities[i];
                sphere* ing = activeIngr[nls];
                std::cout << "#nmol Fill5else is for ingredient : "<< ingredients[i].name<< ' '  << ingredients[i].nbMol<< std::endl ;
                totalNumMols += ing->nbMol;
                nls++;
            }
            std::cout << "#totalNumMols Fill5else = " << totalNumMols << std::endl;
        } 
        return totalNumMols;   
    }

    void updatePriorities(sphere *ingr){
        vRangeStart = vRangeStart + normalizedPriorities[0];
        //ingr->completion = 1.0;
        int ind = 0 ;
        for (int i=0; i < activeIngr.size() ; i++){
            if (ingr->name == activeIngr[0]->name){
                ind = i;
                break;            
            }        
        }

        if (ind > 0){
            for (int j=0; j < ind;j++){         
                thresholdPriorities[j] = thresholdPriorities[j] + normalizedPriorities[ind];
            }
        }
        //remove ingredetien fomr activeIngr
        activeIngr.erase( activeIngr.begin()+ind );
        //# Start of massive overruling section from corrected thesis file of Sept. 25, 2012
        //#this function also depend on the ingr.completiion that can be restored ?
        getSortedActiveIngredients();        
        totalPriorities = 0; // 0.00001
        float pp;
        float np;
        for(unsigned i = 0; i < activeIngr12.size(); ++i) {
            pp = activeIngr12[i]->packingPriority;
            totalPriorities = totalPriorities + pp;
        }
        float previousThresh = 0;
        normalizedPriorities.clear();
        thresholdPriorities.clear();
        //# Graham- Once negatives are used, if picked random# 
        //# is below a number in this list, that item becomes 
        //#the active ingredient in the while loop below
        for(unsigned i = 0; i < activeIngr0.size(); ++i) {
            normalizedPriorities.push_back(0.0);
            if (pickWeightedIngr)
                thresholdPriorities.push_back(2.0);   
        } 
        for(unsigned i = 0; i < activeIngr12.size(); ++i) {
        //for priors in self.activeIngr12:
            //#pp1 = 0
            pp = activeIngr12[i]->packingPriority;
            if (totalPriorities != 0)
                np = float(pp)/float(totalPriorities);
            else
                np=0.0;
            normalizedPriorities.push_back(np);
            if (DEBUG)  std::cout << "#np is "<< np << " pp is "<< pp << " tp is " << np + previousThresh << std::endl;
            thresholdPriorities.push_back(np + previousThresh);
            previousThresh = np + float(previousThresh);
        }
        activeIngr = activeIngr0;// + self.activeIngr12
        activeIngr.insert(activeIngr.end(), activeIngr12.begin(), activeIngr12.end());
    }

    void dropIngredient(sphere *ingr){
        std::cout << "#drop ingredient " << ingr->name << " " << ingr->nbMol << " " << ingr->counter << " "<< ingr->rejectionCounter <<std::endl;
        int ingr_ind;
        bool found = false;
        for(unsigned i = 0; i < droped_ingredient.size(); ++i) {
            if (ingr->name == droped_ingredient[i]->name){
            ingr->active = false;
            ingr->completion = 1.0;
            std::cout << "#already drop ? \n";
            return;
            }
        }  
        for(unsigned i = 0; i < numActiveIngr; ++i) { 
            if (ingr->name == activeIngr[i]->name){ 
                ingr_ind = i;
                found = true;
                break;             
            }
        }
        if (found) {
            //swap  
            numActiveIngr--;
            std::swap(activeIngr[ingr_ind], activeIngr[numActiveIngr]);  
            //swap also the ingredient ? active Ingredient should be unsign int ? 
            //ingr->active = false;   
        }
        ingr->active = false; 
        getMaxRadius();
        //update thresholdproperties
        updatePriorities(ingr);
        droped_ingredient.push_back(ingr);
    }

    sphere* pickIngredient(){
        sphere* ingr;
        //float prob;
        //std::default_random_engine generator (seed);
        //std::default_random_engine generator(seed);
        //std::uniform_real_distribution<float> distribution(0.0,1.0);
        //std::normal_distribution<float> distribution(0.0,1.0);        
        int ingrInd;
        float threshProb;
        float r;
        int n ;
        if (pickWeightedIngr){ 
            if (thresholdPriorities[0] == 2){
                //# Graham here: Walk through -priorities first
                ingr = activeIngr[0];
            }else{
                //#prob = uniform(vRangeStart,1.0)  #Graham 9/21/11 This is wrong...vRangeStart is the point index, need active list i.e. thresholdPriority to be limited
                float prob = uniform(generator);
                ingrInd = 0;
                for(unsigned i = 0; i < thresholdPriorities.size(); ++i) {
                    threshProb = thresholdPriorities[i];
                    if (prob <= threshProb)
                        break;
                    ingrInd = ingrInd + 1;
                }
                //std::cout << "pick Ingredient "  << ingrInd << ' ' << prob <<' ' <<threshProb << std::endl;
                    
                if (ingrInd <  numActiveIngr)
                    ingr = activeIngr[ingrInd];
                else {
                    std::cout << "#error in histoVol pick Ingredient "  << ingrInd << ' ' <<numActiveIngr << std::endl;
                    ingr = activeIngr[0];
                }
            }
        }else {
            //pick random
            //ingr = activeIngr[rand() % numActiveIngr];
            ingr = activeIngr[(int) (uniform(generator) * numActiveIngr)];
            //uniform(generator)
        }
        return ingr;
    }

    inline sphere sample_ingredient() const {
        //randomly pick a ingredient from the list of ingredients
        return ingredients[rand() % numActiveIngr]; 
    }

    inline unsigned sample_empty() const {
        //randomly pick a point
        return empty[rand() % num_points];//num_empty]; 
    }
    
    inline unsigned sample_empty_distance(float radius) const {
        return 0;
    }

    inline unsigned sample_closest_empty(std::vector<float> Dist,std::vector<unsigned> PointID) const {
        return 0;
    }

    int getPointToDrop_i(sphere ingr, float radius,float jitter){
        //#should we take in account a layer? cuttof_boundary and cutoff_surface?
        unsigned ptInd;
        float cut = radius-jitter;
        float d;
        std::vector<unsigned> allIngrPts;
        std::vector<float> allIngrDist;
        if (ingr.packingMode=="close"){
            //for all inactive in the grid ?
            //get the value
            for(unsigned i = 0; i < num_empty; ++i) {
            //for pt in freePoints:#[:nbFreePoints]:
                d = distance[empty[i]];//#look up the distance
                if (d>=cut)//problem distance == max are actually spread along the grid...should be 999999
                    allIngrPts.push_back(empty[i]);
                    allIngrDist.push_back(d);
            }
        }
        else{
            if (mode == 1) {//distance
                for(unsigned i = 0; i < num_empty; ++i) { 
                    d=distance[empty[i]];
                    if (d>=cut) allIngrPts.push_back(empty[i]);
                }       

            }
            else if (mode == 0) {//pure random
                allIngrPts = empty;
            }
        }
        //std::cout << "allIngrPts " <<allIngrPts.size()<<std::endl;
        if (allIngrPts.size()==0){
            sphere tmp[] = {ingr}; 
            ingr.completion = 1.0;
            vRangeStart = vRangeStart + normalizedPriorities[0];
            std::cout << "#drop ingredent no more point for it " << std::endl;
            dropIngredient(&ingr);
            totalPriorities = 0; //# 0.00001
            return -1;                   
        }
        if (pickRandPt){
            //std::cout << "allIngrPts " <<allIngrPts.size()<<std::endl;
            if (ingr.packingMode=="close")
                ptInd = 0;//sample_closest_empty(allIngrDist,allIngrPts);
            else if (ingr.packingMode=="gradient") //&& (use_gradient)  
                ptInd =0;// self.gradients[ingr.gradient].pickPoint(allIngrPts) 
            else{
                ptInd = allIngrPts[rand() % allIngrPts.size()];
                if (ptInd > allIngrPts.size()) ptInd = allIngrPts[0];            
            }     
        }else {
            std:sort(allIngrPts.begin(),allIngrPts.end());//-(allIngrPts.size()-numActiveIngr)
            ptInd = allIngrPts[0];
        }
        return ptInd;
    }

    //use minRadius ?
    openvdb::Coord getPointToDropCoord(sphere* ingr, float radius,float jitter,int *emptyList){
        float cut = radius-jitter;//why - jitter ?
        float d;
        float mini_d=dmax;
        *emptyList = 0;
        if (DEBUG) std::cout << "#getPointToDropCoord " << cut << " " << mini_d <<std::endl;        
        openvdb::Coord cijk;
        openvdb::Coord mini_cijk;
        openvdb::FloatGrid::Accessor accessor_distance = distance_grid->getAccessor();
        std::vector<openvdb::Coord> allIngrPts;
        std::vector<float> allIngrDist;
        if (DEBUG) std::cout << "#retrieving available point from global grid " <<std::endl;  
        bool notfound = true;      
        for (openvdb::FloatGrid::ValueOnIter  iter = distance_grid->beginValueOn(); iter; ++iter) {
            d=iter.getValue();
            if (d>=cut){
                openvdb::Coord cc=iter.getCoord();
                //return getU(dim,cc);//return first found
                allIngrPts.push_back(cc);
                if (d < mini_d){
                    //if the point alread visisted and rejected
                    if (visited_rejected_coord.size() != 0) notfound =  (std::find(visited_rejected_coord.begin(), visited_rejected_coord.end(), cc) == visited_rejected_coord.end());
                    if (notfound) {
                		 // not found
                        mini_d = d;
                        mini_cijk = openvdb::Coord( cc.asVec3i() );
                    }               
                }
            }
        }
        
        if (DEBUG) std::cout << "#allIngrPts size " <<allIngrPts.size() << " nempty " << num_empty << " cutoff " << cut << " minid " << mini_d<<std::endl;
        if (allIngrPts.size()==0){
            //sphere* tmp[] = {ingr}; 
            ingr->completion = 1.0;
            vRangeStart = vRangeStart + normalizedPriorities[0];
            std::cout << "# drop no more point \n" ;
            dropIngredient(ingr); 
            //getSortedActiveIngredients();
            totalPriorities = 0; //# 0.00001
            *emptyList = 1;
            return openvdb::Coord(0,0,0);                   
        }
        if (pickRandPt){
            //std::cout << "allIngrPts " <<allIngrPts.size()<<std::endl;
            if (ingr->packingMode=="close"){
                
                if (mini_d == dmax){
                    cijk = allIngrPts[0];                    
                    }
                //want the smallest distance
                else {
                    cijk = mini_cijk;
                }
                std::cout << "#dist " << mini_d << " =dmax " << (mini_d == dmax) << " ptsijk " << cijk<<std::endl;
            }
            //    ptInd = 0;//sample_closest_empty(allIngrDist,allIngrPts);
            //else if (ingr->packingMode=="gradient") //&& (use_gradient)  
            //    ptInd =0;// self.gradients[ingr.gradient].pickPoint(allIngrPts) 
            //else{
            //else cijk = allIngrPts[rand() % allIngrPts.size()]; // is this working correctly?
            else cijk = allIngrPts[(int)(uniform(generator) * allIngrPts.size())];
            //if (ptInd > allIngrPts.size()) ptInd = allIngrPts[0];            
            //}     
        }else {
            std:sort(allIngrPts.begin(),allIngrPts.end());//-(allIngrPts.size()-numActiveIngr)
            cijk = allIngrPts[0];
        }
        return cijk;        

    }

    int getPointToDrop(sphere* ingr, float radius,float jitter){
        //#should we take in account a layer? cuttof_boundary and cutoff_surface?
        unsigned ptInd;
        float cut = radius-jitter;
        float d;
        //openvdb::Coord cijk;
        openvdb::FloatGrid::Accessor accessor_distance = distance_grid->getAccessor();
        std::vector<unsigned> allIngrPts;
        std::vector<float> allIngrDist;
        if (ingr->packingMode=="close"){
            //for all active in the grid ?
            for (openvdb::FloatGrid::ValueOnIter  iter = distance_grid->beginValueOn(); iter; ++iter) {
                d = iter.getValue();
                if (d>=cut) {
                    openvdb::Coord cc=iter.getCoord();
                    allIngrPts.push_back(getU(dim,cc));
                    allIngrDist.push_back(d);
                }                
            }
            //get the value
            /*for(unsigned i = 0; i < num_empty; ++i) {
                openvdb::Coord cijk = getIJKc(empty[i],dim);
                d=accessor_distance.getValue(cijk);
                if (d>=cut) {
                    allIngrPts.push_back(empty[i]);
                    allIngrDist.push_back(d);}
            } */      
        }
        else{
            for (openvdb::FloatGrid::ValueOnIter  iter = distance_grid->beginValueOn(); iter; ++iter) {
                d=iter.getValue();
                if (d>=cut){
                    openvdb::Coord cc=iter.getCoord();
                    //return getU(dim,cc);//return first found
                    allIngrPts.push_back(getU(dim,cc));
                }
            }
            //for(unsigned i = 0; i < num_empty; ++i) {
            //    openvdb::Coord cijk = getIJKc(empty[i],dim);
            //    d=accessor_distance.getValue(cijk);
            //    if (d>=cut)
            //        allIngrPts.push_back(empty[i]);
            //}
            //problem where no distance > cut
            //for (openvdb::FloatGrid::ValueAllIter  iter = distance_grid->beginValueAll(); iter; ++iter) {
            //    openvdb::Coord cc=iter.getCoord();
            //    if (bbox.isInside(cc)){
            //        if ((iter.getValue() == dmax )||(iter.getValue() >= cut))  {
            //            allIngrPts.push_back(getU(dim,cc));
            //        }            
            //    }       
            //}
            /*
            for(unsigned i = 0; i < num_empty; ++i) {
                openvdb::Coord cijk = getIJKc(empty[i],dim);
                d=accessor_distance.getValue(cijk);
                //std::cout << "#i " << i << " " << empty[i] << " " <<  cijk << " "  << d  << " " << dim << std::endl ;
                if (d>=cut) allIngrPts.push_back(empty[i]);
            } */      
        }
        //num_empty = allIngrPts.size();
        if (DEBUG) std::cout << "#allIngrPts " <<allIngrPts.size() << " " << num_empty << " " << cut << " " << dim <<std::endl;
        if (allIngrPts.size()==0){
            //sphere* tmp[] = {ingr}; 
            ingr->completion = 1.0;
            vRangeStart = vRangeStart + normalizedPriorities[0];
            std::cout << "# drop no more point \n" ; 
            dropIngredient(ingr); 
            //getSortedActiveIngredients();
            totalPriorities = 0; //# 0.00001
            return -1;                   
        }
        if (pickRandPt){
            //std::cout << "allIngrPts " <<allIngrPts.size()<<std::endl;
            if (ingr->packingMode=="close")
                ptInd = 0;//sample_closest_empty(allIngrDist,allIngrPts);
            else if (ingr->packingMode=="gradient") //&& (use_gradient)  
                ptInd =0;// self.gradients[ingr.gradient].pickPoint(allIngrPts) 
            else{
                ptInd = allIngrPts[rand() % allIngrPts.size()];
                if (ptInd > allIngrPts.size()) ptInd = allIngrPts[0];            
            }     
        }else {
            std:sort(allIngrPts.begin(),allIngrPts.end());//-(allIngrPts.size()-numActiveIngr)
            ptInd = allIngrPts[0];
        }
        return ptInd;
    }




    


    bool try_drop(unsigned pid,sphere *ingr)  {
        //main function that decide to drop an ingredient or not
        //std::cout  <<"test_data "<< ingr.name << ' ' << pid << std::endl;        
        int i_ijk[3];
        i_ijk[0]=0;
        i_ijk[1]=0;
        i_ijk[2]=0;
        getIJK(pid,dim,i_ijk);
        openvdb::Coord cijk(i_ijk[0],i_ijk[1],i_ijk[2]);
        return try_dropCoord(cijk,ingr);        
    }

    bool try_dropCoord(openvdb::Coord cijk,sphere *ingr)  {
        float rad = ingr->radius;        
        openvdb::FloatGrid::Accessor accessor_distance = distance_grid->getAccessor(); 
        float d = accessor_distance.getValue(cijk);
        //if (d < ingr.radius) {
        //    std::cout  <<"exit because radius "<< d << " "<< ingr.radius << " " << cijk << " " << pid << std::endl; 
        //    return true;
        // }
        openvdb::Vec3f center=distance_grid->indexToWorld(cijk);
        point px;   //the selected point where we want to drop
        px.x = center.x();
        px.y = center.y();        
        px.z = center.z();
        point target;           //the point with some jitter
        d = 0.0;          //distance to be computed
        unsigned nbJitter = ingr->nbJitter;  //nb of jitter
        //should be defined in ingredient
        float jx=ingr->jitterMax.x();           //jitter amount on x 
        float jy=ingr->jitterMax.y();            //jitter amount on y 
        float jz=ingr->jitterMax.z();            //jitter amount on z 
        //actuall jitter that will be apply to the point
        float dx=0.0;
        float dy=0.0;
        float dz=0.0;
        //square jitter
        float d2;
        float jitter = space;
        float jitter2 = jitter * jitter;
        bool collision;
        // -- Constructing a uniform, cell-centered transform --
        
        // The offset to cell-center points
        openvdb::math::Vec3f offset;//(delta/2., delta/2., delta/2.);
        openvdb::math::Mat4d rotMatj;
        rotMatj.setToRotation(openvdb::math::Vec3f(rand(),rand(),rand()),rand()*M_PI); // random value for axe and angle in radians
        //setTranslation>?
        
        //rotMatj=histoVol.randomRot.get() 

        //prepare the normal distribution for generating the jitter offset
        //std::default_random_engine generator;
        //std::normal_distribution<float> distribution(0.0,0.3);
        //std::cout  << "testData" << std::endl;
        unsigned totnbJitter=0;//total number of jitter
        for(unsigned i = 0; i < nbJitter; ++i) { 
            if (jitter2 > 0.0){
                bool found = false;//already found a good poisition ?
                while (!found){
                    //genereate the offset on x,y,z and the square value
                    dx = jx*jitter*gauss(generator);
                    dy = jy*jitter*gauss(generator);
                    dz = jz*jitter*gauss(generator);
                    d2 = dx*dx + dy*dy + dz*dz;
                    if (d2 < jitter2){
                        //if compNum > 0: # jitter less among normal
                        //    #if self.name=='2uuh C4 SYNTHASE':
                        //    #    import pdb
                        //    #    pdb.set_trace()
                        //    dx, dy, dz, dum = numpy.dot(rotMat, (dx,dy,dz,0))
                        target.x = px.x + dx;// = (tx+dx, ty+dy, tz+dz)
                        target.y = px.y + dy;//
                        target.z = px.z + dz;//
                        found = true;
                    }else{
                        //print('JITTER REJECTED', d2, jitter2)
                    }
                }
            }else{
                target.x = px.x;
                target.y = px.y;
                target.z = px.z;
                dx = dy = dz = 0.0;
            }
            dx = dy = dz = 0.0;
            ++totnbJitter;
            target.x = px.x + dx;// = (tx+dx, ty+dy, tz+dz)
            target.y = px.y + dy;//
            target.z = px.z + dz;//
            offset = openvdb::math::Vec3f(target.x,target.y,target.z);
            //rotMatj.setToRotation(openvdb::math::Vec3f(rand(),rand(),rand()),rand()%M_PI);
            if (ingr->useRotAxis){
                if (ingr->rotAxis.length() == 0.0)  rotMatj.setIdentity();
                else rotMatj.setToRotation(ingr->rotAxis,uniform(generator)*ingr->rotRange);
            }else {
                rotMatj.setToRotation(openvdb::math::Vec3f(uniform(generator),uniform(generator),uniform(generator)),uniform(generator)*2.0*M_PI);
            }
            //can get translate the grid correctly for some reason
            //jitterLength += dx*dx + dy*dy + dz*dz  //#Why is this expensive line needed?
            //jitterList.append( (dx,dy,dz) )      
            //check for collision at the given target point coordinate for the given radius      
            collision = checkSphCollisions(target,rotMatj,rad,ingr);
            //if (DEBUG) std::cout << "#" << rotMatj << "collide ? " << collision << std::endl;
            if (!collision) {
                openvdb::math::Vec3f offset(target.x,target.y,target.z);
                
                ingr->trans = openvdb::math::Vec3f(target.x,target.y,target.z);
//                std::cout << pid << " test data " << target.x <<' ' << target.y << ' ' << target.z <<' ' << collision << std::endl; 
//                std::cout << pid << " test data " << ingr.trans.x <<' ' << ingr.trans.y << ' ' << ingr.trans.z <<' ' << collision << std::endl; 
                rtrans.push_back(openvdb::math::Vec3f(target.x,target.y,target.z)); 
                rrot.push_back(openvdb::math::Mat4d(rotMatj));
                results[rtrans.size()-1] = ingr;

                if (DEBUG) std::cout << "#combine the grid "<< std::endl;

                // A linear transform with the correct spacing
                openvdb::math::Transform::Ptr sourceXform =
                    openvdb::math::Transform::createLinearTransform(ingr->stepsize);//should be stepsize or minRadius
                //which stepsize for the target transform
                openvdb::math::Transform::Ptr targetXform =
                    openvdb::math::Transform::createLinearTransform();//stepsize ?
                // Add the offset.
                openvdb::Vec3d woffset = distance_grid->worldToIndex(offset);//offset assume the stepsize
                targetXform->preMult(rotMatj);
                targetXform->postTranslate(woffset);

                openvdb::Vec3d vcc;
                openvdb::Vec3d cc;
                openvdb::Coord spcc;
                openvdb::Coord ci;
                openvdb::Vec3f spxyz;
                bool accepted = true;
                // Save copies of the two grids; compositing operations always
                // modify the A grid and leave the B grid empty.
                if (DEBUG) std::cout << "#duplicate ingredient grid "<< std::endl;
                openvdb::FloatGrid::Ptr copyOfGridSphere = openvdb::FloatGrid::create(dmax);
                //openvdb::FloatGrid::Ptr copyOfGridSphere = ingr.gsphere->deepCopy();
                copyOfGridSphere->setTransform(targetXform);

                openvdb::Mat4R xform =
                    sourceXform->baseMap()->getAffineMap()->getMat4() *
                    targetXform->baseMap()->getAffineMap()->getMat4().inverse();
                
                // Create the transformer.
                openvdb::tools::GridTransformer transformer(targetXform->baseMap()->getAffineMap()->getMat4());
                
                // Resample using nearest-neighbor interpolation.
                transformer.transformGrid<openvdb::tools::PointSampler, openvdb::FloatGrid>(
                    *ingr->gsphere, *copyOfGridSphere);
                copyOfGridSphere->tree().prune();
                //openvdb::tools::resampleToMatch<openvdb::tools::PointSampler>(ingr.gsphere,copyOfGridSphere);
                //openvdb::tools::compMin(*distance_grid, *ingr.gsphere);
                if (DEBUG) std::cout << "#combine grid "<< std::endl;
                distance_grid->tree().combineExtended(copyOfGridSphere->tree(), Local::min);//b is empty after
                if (DEBUG) std::cout << "#combine grid OK "<< std::endl;
                //problem doesnt fill everywhere....
                //maybe should use a dense fill instead of sparse.?
                //openvdb::tools::compMin(*distance_grid, *copyOfGridSphere);
                //ingr.gsphere = copyOfGridSphere->deepCopy();
                //update empty list
                //return collision;
                num_empty=0;
                for (openvdb::FloatGrid::ValueAllIter  iter = distance_grid->beginValueAll(); iter; ++iter) {
                    //create a sphere with color or radius dependant of value ?
                    ci=iter.getCoord();
                    if (bbox.isInside(ci)){
                        //std::cout << "inside \n";
                        //if (iter.getValue() < 0.0) {iter.setActiveState(false);}
                        //if ((iter.getValue() > 0.0)&&(iter.getValue() > 0.0) {iter.setActiveState(true);}
                        if (iter.getValue() > 0.0){
                             //std::cout << "#off "<<ci<<" "<<iter.getValue()<<std::endl;                             
                             //empty[num_empty] = getU(dim,ci);
                             //all[num_empty] = getU(dim,ci);
                             //num_empty++;
                             iter.setActiveState(true);                         
                        }
                        //if (iter.getValue() >= ingr->maxiVal){
                        //    iter.setValue(dmax);                        
                        // }
                        //if (iter.getValue() > 0.0) {iter.setValueOff();}
                        
                    }
                }
                num_empty=distance_grid->activeVoxelCount();
                if (DEBUG) std::cout << "#update num_empty "<< num_empty << " " << distance_grid->activeVoxelCount() << std::endl;
                //if (DEBUG) std::cout << "#num_empty "<< num_empty <<std::endl;
                // Compute the union (A u B) of the two level sets.
                //openvdb::tools::csgUnion(*distance_grid, *ingr.gsphere);
                /*
                openvdb::FloatGrid::Accessor accessor_distance = distance_grid->getAccessor();
                for (openvdb::FloatGrid::ValueAllCIter iter = ingr.gsphere->cbeginValueAll(); iter; ++iter) {
                    float dist = iter.getValue();//inside or outside // after / before the translation
                    spcc = iter.getCoord();//ijk or xyz
                    //std::cout << "#spcc " << spcc << " " <<ingr.bbox<< std::endl;
                    if  (ingr.bbox.isInside(spcc)){
                        //std::cout << "#spcc " << spcc << std::endl;
                        spxyz = transform->indexToWorld(spcc);
                        //spxyz = ingr.gsphere->indexToWorld(spcc);
                        //std::cout << "#spxyz " << spxyz << " " << std::endl;
                        //spxyz = spxyz + offset;
                        //spxyz = spcc*matrix;
                        //std::cout << "#spxyz " << spxyz << " " << std::endl;
                        cc=distance_grid->worldToIndex(spxyz);//spcc+woffset;//
                        //std::cout << "#cc " << cc << std::endl;
                        ci = openvdb::Coord((int)round(cc.x()),(int)round(cc.y()),(int)round(cc.z()));
                        //test if ci in bb ?
                        //std::cout << "#ci " << ci << std::endl;
                        float v  = accessor_distance.getValue(ci);
                        //std::cout << "#dist " << dist << " " << v << ci << std::endl;
                        if ((dist < v ))   {     
                            std::cout << "#dist " << dist << " " << v << ci << std::endl;
                            accessor_distance.setValue(ci,dist);
                            if (dist < 0.0) {
                                accessor_distance.setValueOff(ci);
                                //swap
                                if (bbox.isInside(ci)) set_filled(getU(dim,ci));
                                std::cout << "#set_filled_dist " << ci <<  getU(dim,ci) << " " << dist << " " << v << std::endl;
                                //--N;
                                }//only if inside sphere
                        }
                    }
                }  */
                //distance_grid->signedFloodFill();  
                //openvdb::FloatGrid::ConstPtr copyOfGridSphere = ingr.gsphere->deepCopy();
                //openvdb::tools::compMin(*distance_grid, *ingr.gsphere);maetof
                //ingr.gsphere = copyOfGridSphere->deepCopy();
                return collision;
            }
        }  
        return collision;
    }

    void set_filled(unsigned i) {
        //marke the given indice as occupide and remove it by swapping
        //Oleg contribution
        unsigned empty_index = all[i];
        if(empty_index < num_empty) { // really is empty
            --num_empty;
            unsigned filled_index = empty[num_empty]; 
            std::swap(empty[empty_index], empty[num_empty]); 
            std::swap(all[i], all[filled_index]);
        }
    }

    inline bool is_empty(unsigned i) const { 
        //check the statut of this point
        return all[i] < num_empty;
    }

    void updateDistance(point pos, float radii){       
    }

    bool checkSphCollisions(point pos,openvdb::math::Mat4d rotMatj, float radii, sphere* sp)  {
        openvdb::math::Vec3f offset(pos.x,pos.y,pos.z);//where the sphere should be.
        openvdb::Vec3d woffset = distance_grid->worldToIndex(offset);
        openvdb::math::Transform::Ptr transform =
            openvdb::math::Transform::createLinearTransform(sp->stepsize);//
        transform->postTranslate(woffset);
        //transform->worldToIndex
        //sp.gsphere->setTransform(transform);
        openvdb::Vec3d cc;
        openvdb::Coord spcc;
        openvdb::Coord ci;
        openvdb::Vec3f spxyz;
        bool collide = false;
        openvdb::FloatGrid::Accessor accessor_distance = distance_grid->getAccessor();

        //if intersect
        //should actually go through the sphere bounding box//ie object bounding 
        //value on is the shell...which could be the radius?
        float d;
        for (openvdb::FloatGrid::ValueAllIter iter = sp->gsphere->beginValueAll(); iter; ++iter) {
            spcc = iter.getCoord();//ijk or xyz
            d = iter.getValue();
            if  ((sp->bbox.isInside(spcc))&&(d<0) ){
                //std::cout << "#spcc " << spcc << std::endl;
                spxyz = sp->gsphere->indexToWorld(spcc);
                //std::cout << "#spxyz " << spxyz << " " << std::endl;
                //apply rotation
                spxyz =rotMatj.transform(spxyz);
                spxyz = spxyz + offset;
                //spxyz = transform->indexToWorld(spcc);
                //spxyz = spcc*matrix;
                //std::cout << "#spxyz " << spxyz << " " << std::endl;
                cc=distance_grid->worldToIndex(spxyz);//spcc+woffset;//
                //std::cout << "#cc " << cc << std::endl;
                ci = openvdb::Coord((int)round(cc.x()),(int)round(cc.y()),(int)round(cc.z()));
                //test if ci in bb ?
                //std::cout << "#ci " << ci << std::endl;
                //if  (!bbox.isInside(ci)){collide = true;return true;continue;}//or reject ?
                if  (!bbox.isInside(ci)){
                        //some point are not inside the grid but could update as well?
                        //continue;                       
                }//or reject ? or we can just look for collision outside ? could have an outside layer...
                //like two bounding box
                float v  = accessor_distance.getValue(ci);
                //std::cout << "#v " << v << " " << ci << std::endl;
                //float dist = iter.getValue();
                //std::cout << "dist " << dist << std::endl;
                //check in distance if already a inside value
                if (v < 0.0) {  
                    if (DEBUG) std::cout << "#sphere position " << ci << " " << offset << " " << woffset << " reject point at d " << v <<  std::endl;
                    if (DEBUG) std::cout << "#in sphere xyz " << sp->gsphere->indexToWorld(spcc) << " ijk " << spcc << "  to grid xyz " << spxyz << " ijk " << cc <<std::endl;
                    //reject point
                    //std::cout << "reject point" << std::endl;
                    collide = true;
                    //counterRej++;
                    return true;
                }
            }
        }
        return collide;
    }
};



/* XML CODE */
/* 
parsing information from the autopack setup file as well as the collada mesh file 
*/

float getRadii(std::string str){
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
    int i = 0; 
    std::vector<openvdb::Vec3f> positions;   
    while (i<pos.size()-2){
         openvdb::Vec3f p(pos[i],pos[i+1],pos[i+2]);
         positions.push_back(p);
         i = i+3;
    }
    return positions;
}
std::vector<openvdb::Vec3s> getPositionsS(std::vector<float> pos){
    int i = 0; 
    std::vector<openvdb::Vec3s> positions;   
    while (i<pos.size()-2){
         openvdb::Vec3f p(pos[i],pos[i+1],pos[i+2]);
         positions.push_back(p);
         i = i+3;
    }
    return positions;
}

std::vector<openvdb::Vec3I> getPositionsInt(std::vector<int> pos){
    int i = 0; 
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


std::vector<mesh> getMeshs(std::string path){
    //need a more efficient parser that will get node-transformation-mesh
    
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(path.c_str());
    std::cout << "#Load result: " << result.description() << " library_geometries " << doc.child("COLLADA").child("library_geometries") << std::endl;
    int nbGeom = 0;  
    std::vector<mesh> meshs;  
    for (pugi::xml_node geometry = doc.child("COLLADA").child("library_geometries").first_child(); geometry; geometry = geometry.next_sibling())
    {
        if (std::string(geometry.name()) != "geometry") continue;
        std::cout << "#geometry: " << nbGeom << std::endl;
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
            for (int i = 0; i< fcount.size();i++){
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
            for (int i = 0; i< f.size();i++){
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
        meshs.push_back(mesh3d);  
    }
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
        for (int i = 0; i< fcount.size();i++){
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
big_grid load_xml(std::string path,int _mode,float _seed){
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(path.c_str());
    //std::cout << path << std::endl;
    std::cout << "#Load result: " << result.description() << ", AutoFillSetup name: " << doc.child("AutoFillSetup").attribute("name").value() << std::endl;
    //get the option
    const float smallestObjectSize = atof(doc.child("AutoFillSetup").child("options").attribute("smallestProteinSize").value());
    const float seed = _seed;
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
            meshs = getMeshs(strmeshFile);
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
        ingr.perturbAxisAmplitude = 0.1;
        if (ingredient.attribute("perturbAxisAmplitude")){
            ingr.perturbAxisAmplitude = ingredient.attribute("perturbAxisAmplitude").as_float();
        }
        //packing mode overwrite from xml file
        ingr.packingMode = std::string(ingredient.attribute("packingMode").value());
        _ingredients.push_back(ingr);
        if (DEBUG)std::cout << "#ok ingredient "<< std::endl;
    }
    //make the grid and return it
    const openvdb::Vec3d ibotleft(bb[0],bb[1],bb[2]);
    const openvdb::Vec3d iupright(bb[3],bb[4],bb[5]);
    box nxyz;   
    nxyz.x = unsigned(ceil((iupright.x()-ibotleft.x())/stepsize));//+encapsulatingGrid
    nxyz.y = unsigned(ceil((iupright.y()-ibotleft.y())/stepsize));//+encapsulatingGrid
    nxyz.z = unsigned(ceil((iupright.z()-ibotleft.z())/stepsize));//+encapsulatingGrid
    if (DEBUG)std::cout << "#box "<<ibotleft<<" "<<iupright << std::endl;
    if (DEBUG)std::cout << "#box nxyz"<<nxyz.x<<" "<<nxyz.y<<" "<<nxyz.z << std::endl;
    //should create the grid from a xml file...easier setup
    big_grid g(nxyz.x,nxyz.y,nxyz.z,stepsize,ibotleft,iupright,seed);
    g.vRangeStart=0.0;
    g.pickWeightedIngr=pickWeightedIngr;
    g.pickRandPt =pickRandPt; 
    //g.setMode(0);//1-close packing 0-random
    g.setIngredients(_ingredients);
    return g; 
    //test<big_grid>(bb0,bb1,smallestObjectSize,_ingredients,_mode,_seed,pickWeightedIngr,pickRandPt); // runs in 0.22 s test<big_grid>(n); // runs in 0.06 s
}


void printIngredientGrid(sphere ingr){
    std::cout << "iname = \"" << ingr.name << "\"\n"; 
    std::cout << "inside=[]\n";
    int counter=0;
    //evalMinMax
    openvdb::Coord cc;
    for (openvdb::FloatGrid::ValueAllIter  iter = ingr.gsphere->beginValueAll(); iter; ++iter) {
        //create a sphere with color or radius dependant of value ?
        cc=iter.getCoord();
        if (ingr.bbox.isInside(cc)){
            counter++;
            float d = iter.getValue();
            float dv = d;
            openvdb::Vec3f pos=ingr.gsphere->indexToWorld(cc); //getValue? 
            /*if (iter.getValue() > 0.0 ){
                float d = iter.getValue();
                if (d == dmax) d=90.0;
                std::cout << "distances.append( "<< d <<")\n";
                std::cout << "all.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
            }*/
            //if (iter.getValue() > 0.0 ){
            //    if (d == dmax) dv=g.maxradius;
            //    std::cout << "Distances.append( "<< dv <<")\n";
            //    std::cout << "all.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
           // }
            if (d < 0.0 ) {
            //std::cout << "#" <<  cc.x() <<"," <<cc.y()<<"," << cc.z() <<std::endl;
                std::cout << "inside.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
            //std::cout << "s=helper.Sphere(\"sphere" << counter << "\",parent=parent,pos=pos)[0]\n";
            }
            else if ((d == dmax )||(d >= spherewidth*stepsize)) {
                if (d == dmax) d=spherewidth*stepsize;
                //std::cout << "background.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
                //std::cout << "distances.append( "<< d <<")\n";
            }
            else {
                //std::cout << "outside.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
        }
        }
    }
    //or pointCloud?
    std::cout << "isph=helper.PointCloudObject(iname+\"_inside\",vertices = inside,materials=[[0,0,1]],parent = parentHider)\n";
}

void printTheGrid(big_grid g){
    std::cout << "inside=[]\n";
    std::cout << "outside=[]\n";
    std::cout << "background=[]\n";
    std::cout << "distances=[]\n";
    std::cout << "Distances=[]\n";
    std::cout << "all=[]\n";
    int counter=0;
    openvdb::Coord cc;
    //the main grid, extract some distance information from it
    for (openvdb::FloatGrid::ValueAllIter  iter = g.distance_grid->beginValueAll(); iter; ++iter) {//g.distance_grid
    //for (openvdb::FloatGrid::ValueAllIter  iter = g.ingredients[1].gsphere->beginValueAll(); iter; ++iter) {
        //create a sphere with color or radius dependant of value ?
        cc=iter.getCoord();
        if (g.bbox.isInside(cc)){
        //if (g.ingredients[1].bbox.isInside(cc)){
            counter++;
            float d = iter.getValue();
            float dv = d;
            openvdb::Vec3f pos=g.distance_grid->indexToWorld(cc); //getValue? 
            /*if (iter.getValue() > 0.0 ){
                float d = iter.getValue();
                if (d == dmax) d=90.0;
                std::cout << "distances.append( "<< d <<")\n";
                std::cout << "all.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
            }*/
            //if (iter.getValue() > 0.0 ){
            //    if (d == dmax) dv=g.maxradius;
            //    std::cout << "Distances.append( "<< dv <<")\n";
            //    std::cout << "all.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
           // }
            if (d < 0.0 ) {
            //std::cout << "#" <<  cc.x() <<"," <<cc.y()<<"," << cc.z() <<std::endl;
                std::cout << "inside.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
            //std::cout << "s=helper.Sphere(\"sphere" << counter << "\",parent=parent,pos=pos)[0]\n";
            }
            else if ((d == dmax )||(d >= spherewidth*stepsize)) {
                if (d == dmax) d=spherewidth*stepsize;
                //std::cout << "background.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
                //std::cout << "distances.append( "<< d <<")\n";
            }
            else {
                //std::cout << "outside.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
        }
        }
    }
    std::cout << "isph=helper.PointCloudObject(\"inside\",vertices = inside,materials=[[0,0,1]])\n";
    std::cout << "isph=helper.PointCloudObject(\"background\",vertices = background)\n";
    std::cout << "isph=helper.PointCloudObject(\"outside\",vertices = outside,materials=[[0,1,0]])\n";
}
//main loop is here
int main(int argc, char* argv[])
{   
    float seed = atof(argv[2]);             //random seed
    forceSphere = (bool) atoi(argv[3]);     //force using sphere level set instead of mesh level set
    std::string filename = argv[1];         //xml setuo file
 
    bool ds_ingrgrid = (bool) atoi(argv[4]);//display level set grid for ingredient
    bool ds_grid = (bool) atoi(argv[5]);    //display global set grid
    //float stepsize = 5.0*1.1547;
    srand (seed);
    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    openvdb::initialize();
    
    //load and setup the pack
    big_grid g = load_xml(filename,0,seed);
    //return 0;

    int pt;
    int i_ijk[3]={0,0,0};
    int counterRej=0;
    int rejection=0;
    bool accepted = true;
    //int counterRej=0;
    openvdb::Coord cc;
    std::vector<openvdb::Vec3f> usedPoints;
    std::vector<float> radiis;
    std::vector<openvdb::Vec3f> colors;

    int totalNumMols = g.prepareIngredient();
    std::cout << "#prepare ingredient complete\n";
    //g.setMode(0);//random packing
    int counter = 0;
    int PlacedMols=0;
    int emptyList;
    openvdb::FloatGrid::Accessor accessor_distance = g.distance_grid->getAccessor();
    while(g.num_empty > 0) {     
    //for(unsigned i = 0; i < 2 ; ++i) {//nx*ny*nz
        //pick the object to drop
        if (DEBUG) std::cout << "#begin loop\n";
        if (g.numActiveIngr == 0) {
                std::cout << "#broken by no more ingredient Done!!!****\n";
                break;
        }
        if (g.vRangeStart>1.0){
                std::cout << "#broken by vRange and hence Done!!!****\n";
                break;
        }
        sphere* ingr = g.pickIngredient();
        //sphere ingr= g.sample_ingredient();//sampl using distance information as well ?
        if (!ingr->active) 
        {
                std::cout << "#inactive ingredient continue\n";
                //continue;        
        }
        if (ingr->completion >= 1.0) {
            std::cout << "#ingredient complete 1.0 \n";
            //g.dropIngredient(ingr);            
            //continue;
        }
        if (DEBUG) std::cout  << "#" << ingr->name << " c " << ingr->completion << " counter " << ingr->counter<<" nbmol " << ingr->nbMol <<  " active " << ingr->active <<std::endl;
        //unsigned s = g.sample_empty(); //get the next available point randomly
        //int s = g.getPointToDrop(ingr, ingr->radius,1.0);
        openvdb::Coord s = g.getPointToDropCoord(ingr, ingr->minRadius,1.0,&emptyList);
        if (DEBUG) std::cout  << "#" << s << " " <<emptyList << " " << accessor_distance.getValue(s) << " " << ingr->radius <<std::endl;
        if ( emptyList == -1 ){
            continue;
        }
        //unsigned s = g.sample_empty(); 
        
        //if (mode == 1) s=g.sample_empty_distance(ingr.radius); //get the next available point randomly  
        //else if (mode == 0)  s=g.sample_empty();    
        //unsigned s = g.sample_closest_empty(); 
        //bool collision = g.try_drop(s,ingr);
        bool collision = g.try_dropCoord(s,ingr);
        if (!collision){
            PlacedMols++;
            radiis.push_back(ingr->radius);
            colors.push_back(ingr->color);
            rejection=0;  
            ingr->counter++;
            ingr->completion = (float)ingr->counter/(float)ingr->nbMol;
            if (DEBUG) std::cout << "# main loop accepted c " <<ingr->completion<< " " << ingr->name << " counter " << ingr->counter << " ijk " << s << " nempty " << g.num_empty << " r " << rejection <<  " collide " << collision << std::endl;          
            if (ingr->completion >= 1.0) {
                std::cout << "#ingredient completion no more nmol to place "<< ingr->counter << std::endl;
                g.dropIngredient(ingr);            
            //continue;
                }
            }
        else {
            ingr->rejectionCounter++;
            if (ingr->rejectionCounter > ingr->rejectionThreshold){
                std::cout << "#ingredient completion too many rejection "<< ingr->rejectionCounter << std::endl;
                ingr->completion =1.0;//erase from list ingredient /
                g.dropIngredient(ingr); 
            }
            g.visited_rejected_coord.push_back(s);
            rejection++;  
            counterRej++; 
            if (DEBUG) std::cout << "# main loop rejected " << ingr->name << ' ' << s << ' ' << g.num_empty << ' ' << rejection <<  " collide " << collision << std::endl;            
                 
        }  
        counter++;
        //g.set_filled(s); 
        if   (rejection > 2000) break; 
        //if   (counterRej > 100) break; 
        //if   (counter > 10) break;        
    }
    //The packing is done, we will generate a python script executed in a 3d Host for the vizualisation
    //this can be replace by any output.

    std::cout << "#distance_grid->activeVoxelCount() " << g.num_empty << " " <<g.distance_grid->activeVoxelCount()<<std::endl;
    std::cout << "#main loop " << g.rtrans.size() << " on " << totalNumMols << std::endl;

    std::cout << "pts=[" << std::endl;
    openvdb::Vec3f pos;
    sphere * ingr;
    for(unsigned i = 0; i < g.rtrans.size(); ++i) { 
        ingr = g.results[i];
        openvdb::math::Transform::Ptr targetXform =
            openvdb::math::Transform::createLinearTransform();
        // Add the offset.
        targetXform->preMult(g.rrot[i]);
        targetXform->postTranslate(g.rtrans[i]);//should be woffset ? nope we apply on xyz not on ijk
        openvdb::math::Mat4d mat = targetXform->baseMap()->getAffineMap()->getMat4();
        for (int j = 0 ; j < g.results[i]->positions.size() ; j ++){
            pos = mat.transform(g.results[i]->positions[j]);
            std::cout << '[' << pos.x() <<',' << pos.y() << ',' << pos.z() << ']' << ',' <<std::endl; //rot? 
            //std::cout << mat << ',' << std::endl;
        }        
    }
    std::cout << "]" << std::endl; 
    std::cout << "matrices={}" << std::endl;
    for(unsigned i = 0; i < g.ingredients.size(); ++i) { 
        std::cout << "matrices[\"" << g.ingredients[i].name << "\"]=[]\n";
    }
    //openvdb::Vec3f pos;
    //sphere * ingr;
    for(unsigned i = 0; i < g.rtrans.size(); ++i) { 
        ingr = g.results[i];
        openvdb::math::Transform::Ptr targetXform =
            openvdb::math::Transform::createLinearTransform();
        // Add the offset.
        targetXform->preMult(g.rrot[i]);
        targetXform->postTranslate(g.rtrans[i]);
        openvdb::math::Mat4d mat = targetXform->baseMap()->getAffineMap()->getMat4();
        std::cout << "matrices[\""<< ingr->name <<"\"].append( " << mat  <<")"<< std::endl;
    }
    //std::cout << "]" << std::endl; 
    
    std::cout << "r=[" << std::endl  ;  
    for(unsigned i = 0; i < radiis.size(); ++i) { 
        //std::cout << i << ' ' << radiis[i] << std::endl;
        for (int j = 0 ; j < g.results[i]->radii.size() ; j ++){
            std::cout << g.results[i]->radii[j] << ',' <<std::endl; //rot? 
        }        
        //std::cout << 3.0 << ',' <<std::endl; 
        //std::cout << 5.0 << ',' <<std::endl;  
    }    
    std::cout << "]" << std::endl;

    std::cout << "color=[" << std::endl  ;  
    for(unsigned i = 0; i < colors.size(); ++i) { 
        //std::cout << i << ' ' << usedPoints[i] << std::endl;
        for (int j = 0 ; j < g.results[i]->radii.size() ; j ++){
            std::cout << '[' << colors[i].x() <<',' << colors[i].y() << ',' << colors[i].z() << ']' << ',' <<std::endl; 
        }
        //std::cout << '[' << colors[i].x() <<',' << colors[i].y() << ',' << colors[i].z() << ']' << ',' <<std::endl;  
    }  
     
    std::cout << "]" << std::endl;
    std::cout << "import upy\n" << "helper = upy.getHelperClass()()\n";
    std::cout << "pesph=helper.newEmpty(\"base_sphere\")\n";
    std::cout << "bsph=helper.Sphere(\"sphere\",parent=pesph)[0]\n";
    std::cout << "parent=helper.newEmpty(\"grid parent\")\n";
    std::cout << "parentHider=helper.newEmpty(\"hider parent\")\n";
    std::cout << "parentInstance=helper.newEmpty(\"instances parent\")\n";
    std::cout << "parentdistance=helper.newEmpty(\"grid distance\")\n";
    std::cout << "isph=helper.instancesSphere(\"cpp\",pts,r,pesph,color,None,parent=parent)\n";

    if (ds_grid) printTheGrid(g);
    
    std::cout << "#" << g.num_empty << std::endl;
    std::cout << "#" << g.distance_grid->activeVoxelCount() << std::endl;
    float mini=0.0;
    float maxi=0.0;
    g.distance_grid->evalMinMax(mini,maxi);
    std::cout << "#" << mini << " " <<maxi << std::endl;
    //mesh
    //foreac ingredient load the mesh?
    //if mesh
    //for each inredient get the original mesh and create instance for their different poisition in the grid
    for(unsigned i = 0; i < g.ingredients.size(); ++i) { //segmntaton fault here ?
        if (ds_ingrgrid) printIngredientGrid(g.ingredients[i]);
        //this will create the individual ingredient grid for debugging the voxelization
        //if (g.ingredients[i].filename.empty())
        //    continue;
        std::cout << "# ingr " << g.ingredients[i].radius << std::endl;
        std::cout << "iname = \"" << g.ingredients[i].name << "\"\n";    
        std::cout << "parent=helper.newEmpty(iname+\"_parent\",parent = parentHider)\n";
        std::cout << "helper.read(\"" << g.ingredients[i].filename << "\")\n"; 
        std::cout << "geom = helper.getObject(iname)\n";
        std::cout << "if geom is not None :\n";
        std::cout << "\thelper.rotateObj(geom,[0,math.pi/2.0,0.0])\n";//rotate on X in C4D
        std::cout << "\thelper.reParent(geom,parent)\n";
        std::cout << "\tiparent=helper.newEmpty(iname+\"_iparent\")\n";
        std::cout << "\taxis = " << g.ingredients[i].principalVector << "\n";
        std::cout << "\tipoly = helper.instancePolygon(iname+\"inst\",\n";
        std::cout << "                      matrices=matrices[iname],\n";
        std::cout << "                      mesh=parent,parent = iparent,\n";
        std::cout << "                      transpose= False,colors=["<<g.ingredients[i].color << "],\n";
        std::cout << "                      axis=axis)\n";
        //openvdb use convex - hull, art was suggesting decomposing the mesh to compensate
        //or just use the multisphere as in blood recipe. Using the mesh merging seems difficult because
        // of the transformation inside the collada file.
    }    
    //stdout the grid
}
//to compile assuming all dependancy present (openvdb, pugixml)
//sh ./compile_Example
//to run :
// time ./autopack <setup.xml> <seed> <forceSphere> <dsIngrG> <dsG> > ~/Dropbox/testCpp.py
//example : sh ./compile_Example;time ./autopack  CellScape1.0.xml 12 0 1 1 > ~/Dropbox/testCpp.py
//the command to run in python and visualize the result
//execfile("/Users/ludo/Dropbox/testCpp.py")