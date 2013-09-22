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

//helper to create a multiSpheres ingredient given a list of radii and positions
//if only one radius and one position given we build a uniq sphere.
Ingredient makeMultiSpheres(std::vector<double> radii, int mode, double concentration, 
         double packingPriority,int nbMol,std::string name, openvdb::Vec3d color,
        unsigned nbJitter,openvdb::Vec3d jitterMax,std::vector<openvdb::Vec3d> positions){
    Ingredient sp;
    sp.radii = radii;
    sp.positions = positions;
    sp.radius = radii[0];
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

    if (radii.size() == 1 ){
         sp.gsphere = openvdb::tools::createLevelSetSphere<openvdb::DoubleGrid>(
              /*radius=*/float(sp.radius), /*center=*/positions[0],//where is the sphere is  index or worl
               /*voxel size=*/float(stepsize), /*width=*/float(spherewidth)); 
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
        
        sp.gsphere = openvdb::createLevelSet<openvdb::DoubleGrid>(stepsize, spherewidth);        
        for (std::vector<double>::size_type i =0 ; i < radii.size();i++) {
            //is the position in xyz or ijk ?
            if (DEBUG)std::cout << "#r " << radii[i] << " pos " <<  positions[i] << std::endl;
            openvdb::DoubleGrid::Ptr localSphere = openvdb::tools::createLevelSetSphere<openvdb::DoubleGrid>(
                  /*radius=*/float(radii[i]), /*center=*/positions[i],//where is the sphere is  index or worl
                   /*voxel size=*/float(stepsize), /*width=*/float(spherewidth)); 

            openvdb::tools::csgUnion(*(sp.gsphere), *localSphere);
        }

        for (auto iter = sp.gsphere->beginValueOn(); iter; ++iter) {
                if (iter.getValue() < 0.0){                    
                    iter.setActiveState(false);                         
                }
        }

        sp.gsphere->prune();
        //sp.gsphere->signedFloodFill();
        //sphere bounding box
        //openvdb::Coord dimgrid = sp.gsphere->evalActiveVoxelDim();//ninjnk
        //sphere bounding box
        sp.bbox = sp.gsphere->evalActiveVoxelBoundingBox();
        sp.useRotAxis=false;
    }
    return sp;
}