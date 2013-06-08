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
#include <vector>
#include <random>

#include "Sphere.h"
#include "Types.h"

/* the main class that handle the packing
came from oleg trott code for the bidirectional array swapping.
extended to do the main autopack loop with ingredient and point picking
*/
struct big_grid { 
    
    //Temporarly must be first
    openvdb::FloatGrid::Ptr distance_grid;
    unsigned num_points;        //total number of point in the grid

    // needs 8*n bytes 
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
    
    float grid_volume;    
    float unit_volume;
    float space;                //spacing of the grid (unit depends on user)
    float maxradius;            //the biggest ingredient
    float minradius;            //the biggest ingredient
    float lowestPriority;       //priority the lowest after sorting
    float totalPriorities;
    float vRangeStart;
    point boundingBox0;         //the grid lower left corner coordinate
    point boundingBox1;         //the grid top right coordinate
    int mode;                   //the packing mode random or distance
    bool pickWeightedIngr;
    bool pickRandPt;
    bool use_gradient;
    int numActiveIngr;

    std::default_random_engine generator;
    std::uniform_real_distribution<float> uniform;// (0.0,1.0);
    std::normal_distribution<float> gauss;//(0.0,0.3);
    std::uniform_real_distribution<double> distribution;

    
    openvdb::CoordBBox bbox;
    openvdb::Coord dim;
    //openvdb::FloatGrid::Accessor accessor_distance;
    std::vector<openvdb::Coord> visited_rejected_coord;
    std::map<int, sphere*> results; 

    //the constructor that take as input the sizenor of the grid, the step, and the bouding box
    big_grid(float step, openvdb::Vec3d bot, openvdb::Vec3d up, unsigned seed);

    unsigned int initializeNumPointsCount();

    openvdb::FloatGrid::Ptr initializeDistanceGrid( openvdb::Vec3d bot, openvdb::Vec3d up );

    void setMode(unsigned _mode){
        //set the packing mode, random or distance
        mode = _mode;
    }

    void setIngredients(std::vector<sphere> _ingredients);
    
    void getMaxRadius();

    void getSortedActiveIngredients();

    int prepareIngredient();

    void updatePriorities(sphere *ingr);

    void dropIngredient(sphere *ingr);

    sphere* pickIngredient();
  
    int getPointToDrop_i(sphere ingr, float radius,float jitter);
    
    openvdb::Coord getPointToDropCoord(sphere* ingr, float radius,float jitter,int *emptyList);

    int getPointToDrop(sphere* ingr, float radius,float jitter);

    bool try_drop(unsigned pid,sphere *ingr);

    bool try_dropCoord(openvdb::Coord cijk,sphere *ingr);

    void set_filled(unsigned i);

    bool is_empty(unsigned i) const;

    bool checkSphCollisions(point pos,openvdb::math::Mat4d rotMatj, float radii, sphere* sp);
};