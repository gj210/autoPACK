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
#include <random>

#include "Ingredient.h"
#include "IngredientsDispatcher.h"
#include "ParticeList.h"


/* the main class that handle the packing
came from oleg trott code for the bidirectional array swapping.
extended to do the main autopack loop with ingredient and point picking
*/
struct big_grid { 
    
    //Temporarly must be first
    openvdb::CoordBBox bbox;
    openvdb::DoubleGrid::Ptr distance_grid;
    openvdb::Index64 num_points;        //total number of point in the grid

    IngradientsDispatcher ingredientsDipatcher;
   
    std::vector<openvdb::Vec3d> rtrans;    //the grid 3d coordintates
    std::vector<openvdb::math::Mat4d> rrot;
    std::vector<Ingredient*> results; 
    std::vector<openvdb::Vec3d> rpossitions;    //coordinations of molecules


    openvdb::Index64 num_empty;         //the number of free point available
       
    bool pickRandPt;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniform;// (0.0,1.0);
    std::uniform_real_distribution<double> half_uniform;// (-0.5,0.5);
    std::uniform_real_distribution<double> distribution;
    

    //the constructor that take as input the sizenor of the grid, the step, and the bouding box
    big_grid(std::vector<Ingredient> const & _ingredients, double step, openvdb::Vec3d bot, openvdb::Vec3d up, unsigned seed);

    openvdb::DoubleGrid::Ptr initializeDistanceGrid( openvdb::Vec3d bot, openvdb::Vec3d up, double voxelSize  );

    openvdb::Coord getPointToDropCoord(Ingredient* ingr, double radius, double jitter,int &emptyList);

    bool try_drop(unsigned pid,Ingredient *ingr);

    bool try_dropCoord(openvdb::Coord cijk,Ingredient *ingr);

    openvdb::math::Mat4d generateIngredientRotation( Ingredient const& ingr);

    bool checkSphCollisions(openvdb::Vec3d const& offset,openvdb::math::Mat4d rotMatj, double radii, Ingredient* sp);
    

private:
    openvdb::Vec3d generateRandomJitterOffset(openvdb::Vec3d const& center, openvdb::Vec3d const & ingrJitter );
    openvdb::Vec3d generateCenterJitterOffset( openvdb::Coord const& pointFrom,  openvdb::Coord const& pointToGo, openvdb::Vec3d const& ingrMaxJitter, Ingredient *ingr );
    double calculateValue( double i);
    openvdb::Coord findDirection(openvdb::Coord const& fromPoint, openvdb::Coord const& toCenter);
    void storePlacedIngradientInGrid( Ingredient * ingr, openvdb::Vec3d offset, openvdb::math::Mat4d rotMatj );
	double countDistance( std::vector<openvdb::Vec3d> const & rpossitions, Ingredient *ingr, openvdb::Vec3d const& offset, openvdb::math::Mat4d const& rotMatj );
	openvdb::Vec3d calculatePossition( Ingredient *ingr, openvdb::Vec3d const& offset, openvdb::math::Mat4d const& rotMatj );
	void printVectorPointsToFile(const std::vector<openvdb::Coord> &allIngrPts, const std::string &file, const std::string &radius);
    double countCurrentDistance( openvdb::Coord cijk, Ingredient *ingr );
    openvdb::Coord chooseTheBestPoint( const std::vector<openvdb::Coord> &allIngrPts, Ingredient *ingr );
    openvdb::Coord getGridMiddlePoint( );
    bool checkCollisionBasedOnGridValue( openvdb::math::Vec3d const& offset,openvdb::math::Mat4d rotMatj, Ingredient* sp );
    void placeSphereInTheGrid( openvdb::math::Vec3d const& offset,openvdb::math::Mat4d rotMatj, Ingredient* sp );
};