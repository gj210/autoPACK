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


class Ingradient {
public:
    float radius;
    std::vector<float> radii;
    std::vector<openvdb::Vec3f> positions;
    float minRadius;
    float maxRadius;
    float stepsize;
    //point trans;
    int mode;
    int nbMol;
    float molarity; 
    float completion;
    float packingPriority;
    std::string name;
    std::string packingMode;
    int counter;
    int rejectionCounter;
    int rejectionThreshold;
    unsigned nbJitter;
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

public:
    bool isActive() { return completion < 1.0; }
    void setCount(float volume);
};
