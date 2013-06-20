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

#pragma once

#include "Types.h"
#include "Ingredient.h"

//helper to create an ingredient given a 3d mesh triangles or quads
Ingredient makeMeshIngredient(std::vector<double> radii, int mode, double concentration, 
                                 double packingPriority,int nbMol,std::string name, openvdb::Vec3d color,
                                 unsigned nbJitter,openvdb::Vec3d jitterMax, mesh mesh3d);

//helper to create an ingredient given a different 3d mesh triangles or quads
Ingredient makeMeshesIngredient(std::vector<double> radii, int mode, double concentration, 
                                   double packingPriority,int nbMol,std::string name, openvdb::Vec3d color,
                                   unsigned nbJitter,openvdb::Vec3d jitterMax, std::vector<mesh> meshs);


//helper to create a singleSphere ingredient given a radius, and some options
Ingredient makeSphere(double radius, int mode, double concentration, 
         double packingPriority,int nbMol,std::string name, openvdb::Vec3d color,
        unsigned nbJitter,openvdb::Vec3d jitterMax);

//helper to create a multiSpheres ingredient given a list of radii and positions
//if only one radius and one position given we build a uniq sphere.
Ingredient makeMultiSpheres(std::vector<double> radii, int mode, double concentration, 
         double packingPriority,int nbMol,std::string name, openvdb::Vec3d color,
        unsigned nbJitter,openvdb::Vec3d jitterMax,std::vector<openvdb::Vec3d> positions);