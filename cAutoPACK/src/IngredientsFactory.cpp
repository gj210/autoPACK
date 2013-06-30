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

#include "IngredientsFactory.h"


//helper to create an ingredient given a 3d mesh triangles or quads
Ingredient makeMeshIngredient(std::vector<double> radii, int mode, double concentration, 
                                 double packingPriority,int nbMol,std::string name, openvdb::Vec3d color,
                                 unsigned nbJitter,openvdb::Vec3d jitterMax, mesh mesh3d){
    Ingredient sp;
    sp.molarity=concentration;
    sp.radii = radii;
    sp.positions.resize(1);
    sp.positions[0] = openvdb::Vec3d(0,0,0);
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
    sp.color = color;
    sp.nbJitter = nbJitter;
    sp.jitterMax = jitterMax;
    //build the grid
    //need to create as many grid as sphere, then combine then in one uniq ie union?
    sp.gsphere = openvdb::DoubleGrid::create(dmax);
    double size=stepsize;
    if (stepsize > sp.minRadius) size = stepsize;//sp.minRadius/2.0;
    sp.stepsize = stepsize;
    //levelSet
    if (DEBUG) std::cout << "#mesh have v "<< mesh3d.vertices.size() << " f " << mesh3d.faces.size() << " q " << mesh3d.quads.size() << std::endl;

    if ((mesh3d.quads.size() != 0)&&(mesh3d.faces.size() != 0)){
        sp.gsphere = openvdb::tools::meshToLevelSet<openvdb::DoubleGrid>(
            *openvdb::math::Transform::createLinearTransform(size), 
            mesh3d.vertices, mesh3d.faces, mesh3d.quads);
    }
    else if ((mesh3d.quads.size() != 0)&&(mesh3d.faces.size() == 0)){
        sp.gsphere = openvdb::tools::meshToLevelSet<openvdb::DoubleGrid>(
            *openvdb::math::Transform::createLinearTransform(size), 
            mesh3d.vertices, mesh3d.quads);
    }
    else if ((mesh3d.quads.size() == 0)&&(mesh3d.faces.size() != 0)){
        sp.gsphere = openvdb::tools::meshToLevelSet<openvdb::DoubleGrid>(
            *openvdb::math::Transform::createLinearTransform(size), 
            mesh3d.vertices, mesh3d.faces);
    }

    //sp.gsphere->signedFloodFill();
    //sp.gsphere = openvdb::tools::meshToSignedDistanceField<openvdb::DoubleGrid>(
    //    *openvdb::math::Transform::createLinearTransform(stepsize), mesh3d.vertices, mesh3d.faces,mesh3d.quads,(int)spherewidth,(int)spherewidth);
    //what about sdf ?
    //meshToSignedDistanceField (const openvdb::math::Transform &xform, const std::vector< Vec3s > &points, const std::vector< Vec3I > &triangles, const std::vector< Vec4I > &quads, double exBandWidth, double inBandWidth)
    //openvdb::tools::internal::MeshVoxelizer<openvdb::FloatTree>
    //   voxelizer(mesh3d.vertices, mesh3d.faces);
    //voxelizer.runParallel();
    //sp.gsphere = openvdb::DoubleGrid::create(voxelizer);//use 4I ?
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
Ingredient makeMeshesIngredient(std::vector<double> radii, int mode, double concentration, 
                                   double packingPriority,int nbMol,std::string name, openvdb::Vec3d color,
                                   unsigned nbJitter,openvdb::Vec3d jitterMax, std::vector<mesh> meshs){
    Ingredient sp;
    sp.molarity=concentration;
    sp.radii = radii;
    sp.positions.resize(1);
    sp.positions[0] = openvdb::Vec3d(0,0,0);
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
    sp.color = color;
    sp.nbJitter = nbJitter;
    sp.jitterMax = jitterMax;
    //build the grid
    //need to create as many grid as sphere, then combine then in one uniq ie union?
    sp.gsphere = openvdb::DoubleGrid::create(dmax);
    double size=stepsize;
    if (stepsize > sp.minRadius) size = stepsize;//sp.minRadius/2.0;
    sp.stepsize = stepsize;
    //levelSet
    std::vector<openvdb::DoubleGrid::Ptr> gspheres;
    sp.gsphere->setTransform(
        openvdb::math::Transform::createLinearTransform(/*voxel size=*/size));    
    gspheres.resize(meshs.size());
    if (DEBUG) std::cout << "#merge "<< meshs.size() << " voxelmesh " << std::endl;
    for (std::vector<mesh>::size_type i =0 ; i < meshs.size();i++){
        //is the position in xyz or ijk ?
        if (DEBUG) std::cout << "#mesh have v "<< meshs[i].vertices.size() << " f " << meshs[i].faces.size() << " q " << meshs[i].quads.size() << std::endl;
        if ((meshs[i].quads.size() != 0)&&(meshs[i].faces.size() != 0)){
            gspheres[i] = openvdb::tools::meshToLevelSet<openvdb::DoubleGrid>(
                *openvdb::math::Transform::createLinearTransform(size), 
                meshs[i].vertices, meshs[i].faces, meshs[i].quads);
        }
        else if ((meshs[i].quads.size() != 0)&&(meshs[i].faces.size() == 0)){
            gspheres[i] = openvdb::tools::meshToLevelSet<openvdb::DoubleGrid>(
                *openvdb::math::Transform::createLinearTransform(size), 
                meshs[i].vertices, meshs[i].quads);
        }
        else if ((meshs[i].quads.size() == 0)&&(meshs[i].faces.size() != 0)){
            gspheres[i] = openvdb::tools::meshToLevelSet<openvdb::DoubleGrid>(
                *openvdb::math::Transform::createLinearTransform(size), 
                meshs[i].vertices, meshs[i].faces);
        }
        if (DEBUG) std::cout <<  "#combine " <<size <<std::endl;
        sp.gsphere->tree().combineExtended(gspheres[i]->tree(), Local::rmin);
        sp.gsphere->prune();
        //openvdb::tools::csgUnion(sp.gsphere->tree(),gspheres[i]->tree());
    }
    if (DEBUG) std::cout <<  "#OK merged all of them and step is "<< size <<std::endl;
    //sp.gsphere->prune();
    sp.bbox = sp.gsphere->evalActiveVoxelBoundingBox();
    if (DEBUG) std::cout <<  "#OK bbox "<< sp.bbox << std::endl;    
    sp.useRotAxis=false;
    sp.gsphere->evalMinMax(sp.miniVal,sp.maxiVal);//distance
    // Apply the functor to all active values.
    std::cout << "#radius is "<< sp.minRadius <<" mini " << sp.miniVal << " maxi " << sp.maxiVal << " bb " << sp.bbox<<std::endl;
    //openvdb::tools::foreach(sp.gsphere->beginValueAll(), SetMaxToDefault(sp.maxiVal));
    //translate (const Coord &t)
    return sp;
}


//helper to create a singleSphere ingredient given a radius, and some options
Ingredient makeSphere(double radius, int mode, double concentration, 
         double packingPriority,int nbMol,std::string name, openvdb::Vec3d color,
        unsigned nbJitter,openvdb::Vec3d jitterMax){
    Ingredient sp;
    sp.radius = radius;
    sp.molarity=concentration;
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
    sp.color = color;
    sp.nbJitter = nbJitter;
    sp.radii.push_back(radius);    
    sp.positions.push_back(openvdb::Vec3d(0,0,0));
    sp.jitterMax = jitterMax;
    //build the grid
    double size=stepsize;
    if (stepsize > sp.minRadius) size = sp.minRadius/2.0;
    sp.stepsize = size;

    sp.gsphere = openvdb::tools::createLevelSetSphere<openvdb::DoubleGrid>(
              /*radius=*/float(radius), /*center=*/openvdb::Vec3d(0,0,0),//where is the sphere is  index or worl
               /*voxel size=*/float(size), /*width=*/float(25.0));//width could be maxRadius / stepsize (int)MaxRadius/stepsize
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
Ingredient makeMultiSpheres(std::vector<double> radii, int mode, double concentration, 
         double packingPriority,int nbMol,std::string name, openvdb::Vec3d color,
        unsigned nbJitter,openvdb::Vec3d jitterMax,std::vector<openvdb::Vec3d> positions){
    Ingredient sp;
    sp.radii = radii;
    sp.positions = positions;
    sp.radius = radii[0];
    sp.mode = mode;
    sp.molarity=concentration;    

    auto minMaxIt = std::minmax_element(std::begin(radii), std::end(radii));
    sp.minRadius = *minMaxIt.first;
    sp.maxRadius = *minMaxIt.second;
    if (DEBUG)std::cout << "#miniR " << sp.minRadius << " maxiR " << sp.maxRadius <<std::endl;

    sp.nbMol=nbMol;
    sp.completion = 0.0;
    sp.packingMode = packing;
    sp.packingPriority = packingPriority;
    sp.name = name;
    sp.counter = 0;
    sp.rejectionCounter = 0;
    sp.rejectionThreshold  = rejectionThresholdIngredient;
    sp.color = color;
    sp.nbJitter = nbJitter;
    sp.jitterMax = jitterMax;
    //build the grid
    //need to create as many grid as sphere, then combine then in one uniq ie union?

    double size=stepsize;
    if (stepsize > sp.minRadius) size = stepsize;//sp.minRadius/2.0;
    sp.stepsize = stepsize;
    
    if (radii.size() == 1 ){
         sp.gsphere = openvdb::tools::createLevelSetSphere<openvdb::DoubleGrid>(
              /*radius=*/float(sp.radius), /*center=*/positions[0],//where is the sphere is  index or worl
               /*voxel size=*/float(size), /*width=*/float(spherewidth)); 
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
        std::vector<openvdb::DoubleGrid::Ptr> gspheres;
        sp.gsphere = openvdb::DoubleGrid::create(dmax);
        sp.gsphere->setTransform(
            openvdb::math::Transform::createLinearTransform(/*voxel size=*/size)); 
        gspheres.resize(radii.size());
        for (std::vector<double>::size_type i =0 ; i < radii.size();i++) {
            //is the position in xyz or ijk ?
            if (DEBUG)std::cout << "#r " << radii[i] << " pos " <<  positions[i] << std::endl;
            gspheres[i] = openvdb::tools::createLevelSetSphere<openvdb::DoubleGrid>(
                  /*radius=*/float(radii[i]), /*center=*/positions[i],//where is the sphere is  index or worl
                   /*voxel size=*/float(size), /*width=*/float(spherewidth));
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