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

#ifdef _MSC_VER
    #pragma warning(disable:4146)
    #pragma warning(disable:4503) // OpenVDB "warning decorated name length exceeded, name was truncated"
#endif

#include <vector>

#ifdef _MSC_VER
    //Disable warnings from openvdb in Visual Studio
    #pragma warning(push, 0)   
#endif

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
#include <openvdb/tools/ParticlesToLevelSet.h>

#ifdef _MSC_VER
    //Disable warnings from openvdb in Visual Studio
    #pragma warning(pop)
#endif

const double dmax = 999999.0;        //maximum value assign to the grid for intialisation
const int rejectionThresholdIngredient = 300;//30 by default, number of rejection before stoppin a ingredient
const bool DEBUG = false;                   //DEBUG mode for prining some information
const double spherewidth = 12;             //default or more ? close packing need to increase this number
const std::string packing = "random";       //packing mode, can be random or close

extern double stepsize;         //grid step size ie smallest ingredients radius


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
    static inline void min(openvdb::CombineArgs<double>& args) {
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
    static inline void rmin(openvdb::CombineArgs<double>& args) {
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