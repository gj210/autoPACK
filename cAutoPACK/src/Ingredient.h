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
#include <vector>


class Ingredient {
public:
    double radius;
    std::vector<double> radii;
    std::vector<openvdb::Vec3d> positions;
    double minRadius;
    double maxRadius;
    double stepsize;
    //point trans;
    int mode;
    int nbMol;
    double molarity; 
    double completion;
    double packingPriority;
    std::string name;
    std::string packingMode;
    int counter;
    int rejectionCounter;
    int rejectionThreshold;
    unsigned nbJitter;
    openvdb::Vec3d color;
    openvdb::Vec3d trans;
    openvdb::Vec3d jitterMax;
    openvdb::Vec3d principalVector;
    openvdb::DoubleGrid::Ptr gsphere;
    openvdb::CoordBBox bbox;
    double miniVal;
    double maxiVal;
    mesh mesh3d;
    std::string filename;
    //rotation paramters
    bool useRotAxis;
    double rotRange;
    openvdb::Vec3d rotAxis;
    double perturbAxisAmplitude;

public:
    bool isActive() { return completion < 1.0; }
    void setCount(double volume);
};
