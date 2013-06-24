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

#include <numeric>
#include "BigGrid.h"

namespace {

//some expermental functions for the IJK<->U grid index convertion
//currently unused
/* 
unused but keep in case
inline unsigned getU(openvdb::Coord dim,openvdb::Coord ijk){
    return (int)(ijk.x()*dim.x()*dim.y() + ijk.y()*dim.x() + ijk.z());
}
*/

inline void getIJK(int u,openvdb::Coord dim,int* i_ijk){
    // = {0,0,0};    
    //openvdb::Coord ijk(0,0,0);
    int nxnynz = dim.x()*dim.y()*dim.z();
    int nynz = dim.z()*dim.y();
    //int nx = dim.x();
    int ny = dim.y();
    int nz = dim.z();
    int integer;
    double decimal;
    double fraction;
    int nysc;
    if (u < dim.z()){
        i_ijk[2] = u;
    }
    else if ((u < nynz)&&(u >= nxnynz)){
        //whats z
        fraction = (double)u/(double)dim.z();
        integer = (int) fraction;
        decimal = fraction - integer;
        i_ijk[2] = (int) round(decimal*dim.z());
        //whast y 
        i_ijk[1] = integer;  
    }
    else if ((u < nxnynz)&&(u >= nynz)){
        fraction = (double)u/(double)nynz;
        integer = (int) fraction;
        decimal = fraction - integer;
        nysc = ny * integer;
        //whast x 
        i_ijk[0] = integer;  
        fraction = (double)u/(double)nz;
        integer = (int) fraction;
        decimal = fraction - integer;
        //whats z        
        i_ijk[2] = (int) round(decimal*(double)nz);
        //whast y 
        //46867 
        //233 15477 201 77 603 77.7231
        //std::cout << integer << " " << nysc << " " << ny << " " << (int)((double)u/(double)nynz) << " " << nynz << " " << (double)u/(double)nynz << std::endl;
        i_ijk[1] = integer - (ny*(int)((double)u/(double)nynz));  
        //int (integer - (ny*int(double(u)/double(nynz))));
    }    
}

} //namespace

big_grid::big_grid( std::vector<Ingredient> const & _ingredients, double step, openvdb::Vec3d bot, openvdb::Vec3d up, unsigned seed ) :     
    distance_grid(initializeDistanceGrid(bot, up)),
    num_points(initializeNumPointsCount()),    
    ingredientsDipatcher(_ingredients, num_points, seed),
    num_empty(num_points),
    uniform(0.0,1.0),
    distribution(0.0,1.0),
    gauss(0.0,0.3),        
    pickRandPt(true),
    jitter(step),
    jitterSquare(step*step)
{
    generator.seed(seed);
}

openvdb::Index64 big_grid::initializeNumPointsCount()
{
    dim = distance_grid->evalActiveVoxelDim();    
    std::cout << "#Grid Npoints " << dim << distance_grid->activeVoxelCount() << std::endl;
    return distance_grid->activeVoxelCount();
}

openvdb::DoubleGrid::Ptr big_grid::initializeDistanceGrid( openvdb::Vec3d bot, openvdb::Vec3d up )
{
    openvdb::DoubleGrid::Ptr distance_grid;
    distance_grid = openvdb::DoubleGrid::create();
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

    openvdb::DoubleGrid::Accessor accessor_distance = distance_grid->getAccessor();
    std::cout << "#Testing distance access:" << std::endl;
    std::cout << "#Grid " << left << " "<< botleft << " = " << accessor_distance.getValue(left) << std::endl;
    std::cout << "#Grid " << right << " "<< upright << " = " << accessor_distance.getValue(right) << std::endl;
    std::cout << "#Grid " << bbox << std::endl;
    
    return distance_grid;
}


openvdb::Coord big_grid::getPointToDropCoord( Ingredient* ingr, double radius, double jitter, int *emptyList )
{
    const double cut = radius-jitter;//why - jitter ?
    double d;
    double mini_d=dmax;
    *emptyList = 0;
    if (DEBUG) std::cout << "#getPointToDropCoord " << cut << " " << mini_d <<std::endl;        
    openvdb::Coord cijk;
    openvdb::Coord mini_cijk;
    openvdb::DoubleGrid::Accessor accessor_distance = distance_grid->getAccessor();
    std::vector<openvdb::Coord> allIngrPts;
    std::vector<double> allIngrDist;
    if (DEBUG) std::cout << "#retrieving available point from global grid " <<std::endl;  
    bool notfound = true;      
    //only onValue or all value ? 
    /*
    It might surprise you that the Grid class doesn't directly provide access to 
    voxels via their $(i,j,k)$ coordinates. Instead, the recommended procedure is 
    to ask the grid for a "value accessor", which is an accelerator object that 
    performs bottom-up tree traversal using cached information from previous traversals. 
    Caching greatly improves performance, but it is inherently not thread-safe. 
    However, a single grid may have multiple value accessors, so each thread can 
    safely be assigned its own value accessor. Uncachedand therefore slower, 
    but thread-saferandom access is possible through a grid's tree, for example 
    with a call like grid.tree()->getValue(ijk).
    */
    //bottom-up tree traversal-> how to get real random after
    //this loop populate mostly with bottum up coordinate
    //may need a regular access with the 3 for loop and then grid.tree()->getValue(ijk).     
    for (openvdb::DoubleGrid::ValueOnCIter  iter = distance_grid->cbeginValueOn(); iter; ++iter) {
        //for (openvdb::DoubleGrid::ValueAllIter  iter = distance_grid->beginValueAll(); iter; ++iter) {
        //before getting value check if leaf or tile    
        d=iter.getValue();
        if (d>=cut){//the grid voxel is available and can receivethe given ingredient
            if (iter.isTileValue()){
                openvdb::CoordBBox bbox = iter.getBoundingBox();
                //iterate through it and add all the coordinates
                openvdb::Coord bbmini = bbox.min();
                openvdb::Coord bbmaxi = bbox.max();
                for (int k=bbmini.z();k<bbmaxi.z();k++){
                    for (int j=bbmini.y();j<bbmaxi.y();j++){
                        for (int i=bbmini.x();i<bbmaxi.x();i++){
                            openvdb::Coord nijk(i,j,k);
                            allIngrPts.push_back(nijk);
                            if (d < mini_d){
                                if (visited_rejected_coord.size() != 0) notfound =  (std::find(visited_rejected_coord.begin(), visited_rejected_coord.end(), nijk) == visited_rejected_coord.end());
                                if (notfound) {
                                    // not found
                                    mini_d = d;
                                    mini_cijk = openvdb::Coord( nijk.asVec3i() );
                                }               
                            }   
                        }
                    }
                }
            } 
            else {
                openvdb::Coord cc=iter.getCoord();
                //return getU(dim,cc);//return first found
                allIngrPts.push_back(cc);
                if (d < mini_d){
                    //if the point alread visisted and rejected 
                    //should be for this ingredient only
                    if (visited_rejected_coord.size() != 0) notfound =  (std::find(visited_rejected_coord.begin(), visited_rejected_coord.end(), cc) == visited_rejected_coord.end());
                    if (notfound) {
                        // not found
                        mini_d = d;
                        mini_cijk = openvdb::Coord( cc.asVec3i() );
                    }               
                }
            }           
        }
    }

    if (DEBUG) std::cout << "#allIngrPts size " <<allIngrPts.size() << " nempty " << num_empty << " cutoff " << cut << " minid " << mini_d<<std::endl;
    if (allIngrPts.size()==0){
        std::cout << "# drop no more point \n" ;
        ingredientsDipatcher.dropIngredient(ingr); 
        totalPriorities = 0; //# 0.00001
        *emptyList = 1;
        return openvdb::Coord(0,0,0);                   
    }
    if (pickRandPt){
        std::cout << "#allIngrPts " <<allIngrPts.size() << "mode " << ingr->packingMode <<std::endl;
        if (ingr->packingMode=="close"){
            //try to add here case where there is still space but need anoher starting point.
            if (mini_d == dmax){
                cijk = allIngrPts[0];                    
            }
            //want the smallest distance, but it is alway the same, so we get stuck here...
            //maybe use a weighting system based on distance, closed distance high prob.
            else {
                cijk = mini_cijk;
            }
            
			/*Daniel - doesn't work, always came here when rejectionCounter is equal to 0; */
			if (ingr->rejectionCounter != 0 && ingr->rejectionCounter % 300 == 0){
                cijk = allIngrPts[(int)(distribution(generator) * allIngrPts.size())];
                //increase the threshold ?
                //ingr->rejectionCounter = 0;//probably not enought....risk to never end...
            }
            std::cout << "#dist " << mini_d << " =dmax " << (mini_d == dmax) << " ptsijk " << cijk<<std::endl;
        }
        //    ptInd = 0;//sample_closest_empty(allIngrDist,allIngrPts);
        //else if (ingr->packingMode=="gradient") //&& (use_gradient)  
        //    ptInd =0;// self.gradients[ingr.gradient].pickPoint(allIngrPts) 
        //else{
        //else cijk = allIngrPts[rand() % allIngrPts.size()]; // is this working correctly?
        else {
            //or should I use a std::uniform_int_distribution<int> distribution(0,allIngrPts.size());
            cijk = allIngrPts[(int)(distribution(generator) * allIngrPts.size())];
            //rand tand to get small number first
            //cijk = allIngrPts[rand() % allIngrPts.size()];
            //why this would be on the edge and top. maybe the uniform distribution is not appropriate.
            //

        }
        //if (ptInd > allIngrPts.size()) ptInd = allIngrPts[0];            
        //}     
    }else {
        //ordered ?
        std::sort(allIngrPts.begin(),allIngrPts.end());//-(allIngrPts.size()-numActiveIngr)
        cijk = allIngrPts[0];
    }
    return cijk;
}

bool big_grid::try_drop( unsigned pid,Ingredient *ingr )
{
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

bool big_grid::try_dropCoord( openvdb::Coord cijk,Ingredient *ingr )
{
    const double rad = ingr->radius;        
    bool collision;
    openvdb::Vec3d center=distance_grid->indexToWorld(cijk);

    //actuall jitter that will be apply to the point
   
    openvdb::Vec3d offset;           //the point with some jitter
    for(unsigned i = 0; i < ingr->nbJitter; ++i) { 
        if (jitterSquare > 0.0){
            offset = center + generateRandomJitterOffset(ingr->jitterMax);
        }else{
            offset = center;
        }

        openvdb::math::Mat4d rotMatj;
        if (ingr->useRotAxis){
            if (ingr->rotAxis.length() == 0.0)  
                rotMatj.setIdentity();
            else 
                rotMatj.setToRotation(ingr->rotAxis,uniform(generator)*ingr->rotRange);
        }else {
            rotMatj.setToRotation(openvdb::math::Vec3d(uniform(generator),uniform(generator),uniform(generator)),uniform(generator)*2.0*M_PI);
        }

        //check for collision at the given target point coordinate for the given radius     
        collision = checkSphCollisions(offset,rotMatj,rad,ingr);
        //if (DEBUG) std::cout << "#" << rotMatj << "collide ? " << collision << std::endl;
        if (!collision) {
            ingr->trans = offset;

            rtrans.push_back(offset);
            rrot.push_back(rotMatj);
            results[rtrans.size()-1] = ingr;

            storePlacedIngradientInGrid(ingr, offset, rotMatj);

            if (DEBUG) std::cout << "#update num_empty "<< num_empty << " " << distance_grid->activeVoxelCount() << std::endl;
           
            return collision;
        }
    }  
    return collision;
}


bool big_grid::checkSphCollisions( openvdb::math::Vec3d const& offset,openvdb::math::Mat4d rotMatj, double radii, Ingredient* sp )
{
    openvdb::DoubleGrid::Accessor accessor_distance = distance_grid->getAccessor();

    //if intersect
    //should actually go through the sphere bounding box//ie object bounding 
    //value on is the shell...which could be the radius?
    //can we first test if boudning box hasOverlap ?
    for (openvdb::DoubleGrid::ValueAllCIter iter = sp->gsphere->cbeginValueAll(); iter; ++iter) {
        const openvdb::Coord sphereIndexCoord = iter.getCoord();//ijk or xyz
        const double d = iter.getValue();
        if  ( d<0 && sp->bbox.isInside(sphereIndexCoord) ) {
            //std::cout << "#sphereIndexCoord " << sphereIndexCoord << std::endl;
            openvdb::Vec3d sphereWorldCoord = sp->gsphere->indexToWorld(sphereIndexCoord);
            //std::cout << "#sphereWorldCoord " << sphereWorldCoord << " " << std::endl;
            //apply rotation
            sphereWorldCoord =rotMatj.transform(sphereWorldCoord);
            sphereWorldCoord = sphereWorldCoord + offset;
            //sphereWorldCoord = transform->indexToWorld(sphereIndexCoord);
            //sphereWorldCoord = sphereIndexCoord*matrix;
            //std::cout << "#sphereWorldCoord " << sphereWorldCoord << " " << std::endl;
            openvdb::Vec3d cc=distance_grid->worldToIndex(sphereWorldCoord);//sphereIndexCoord+woffset;//
            //std::cout << "#cc " << cc << std::endl;
            openvdb::Coord ci = openvdb::Coord(openvdb::tools::local_util::roundVec3(cc));
            //test if ci in bb ?
            //std::cout << "#ci " << ci << std::endl;
            //if  (!bbox.isInside(ci)){collide = true;return true;continue;}//or reject ?
            if  (!bbox.isInside(ci)){
                //some point are not inside the grid but could update as well?
                //continue;                       
            }//or reject ? or we can just look for collision outside ? could have an outside layer...
            //like two bounding box
            const double v  = accessor_distance.getValue(ci);
            //std::cout << "#v " << v << " " << ci << std::endl;
            //double dist = iter.getValue();
            //std::cout << "dist " << dist << std::endl;
            //check in distance if already a inside value
            if (v < 0.0) { 
                if (DEBUG) {
                    openvdb::Vec3d woffset = distance_grid->worldToIndex(offset);
                    std::cout << "#sphere position " << ci << " " << offset << " " << woffset << " reject point at d " << v <<  std::endl;
                    std::cout << "#in sphere xyz " << sp->gsphere->indexToWorld(sphereIndexCoord) << " ijk " << sphereIndexCoord << "  to grid xyz " << sphereWorldCoord << " ijk " << cc <<std::endl;
                }
                
                //reject point
                //std::cout << "reject point" << std::endl;
                //counterRej++;
                return true;
            }
        }
    }
    return false;
}

openvdb::Vec3d big_grid::generateRandomJitterOffset( openvdb::Vec3d const& ingrJitter )
{
    while (true) {
        const openvdb::Vec3d randomJitter( 
            jitter*gauss(generator)
            , jitter*gauss(generator)
            , jitter*gauss(generator));
        const openvdb::Vec3d deltaOffset (ingrJitter * randomJitter);
        if ( deltaOffset.lengthSqr() < jitterSquare )
            return deltaOffset;
    } 
}

void big_grid::storePlacedIngradientInGrid( Ingredient * ingr, openvdb::Vec3d offset, openvdb::math::Mat4d rotMatj )
{

    if (DEBUG) std::cout << "#combine the grid "<< std::endl;

    openvdb::math::Transform::Ptr targetXform =
        openvdb::math::Transform::createLinearTransform();//stepsize ?

    // Add the offset.
    const openvdb::Vec3d woffset = distance_grid->worldToIndex(offset);//offset assume the stepsize
    targetXform->preMult(rotMatj);
    targetXform->postTranslate(woffset);

    // Save copies of the two grids; compositing operations always
    // modify the A grid and leave the B grid empty.
    if (DEBUG) std::cout << "#duplicate ingredient grid "<< std::endl;

    openvdb::DoubleGrid::Ptr copyOfGridSphere = openvdb::DoubleGrid::create(dmax);            

    // Create the transformer.
    openvdb::tools::GridTransformer transformer(targetXform->baseMap()->getAffineMap()->getMat4());

    // Resample using nearest-neighbor interpolation.
    transformer.transformGrid<openvdb::tools::PointSampler, openvdb::DoubleGrid>(
        *ingr->gsphere, *copyOfGridSphere);
    copyOfGridSphere->tree().prune();

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
    for (openvdb::DoubleGrid::ValueAllIter  iter = distance_grid->beginValueAll(); iter; ++iter) {
        //create a sphere with color or radius dependant of value ?
        if (bbox.isInside(iter.getCoord())){
            //std::cout << "inside \n";                    
            if (iter.getValue() > 0.0){
                //std::cout << "#off "<<ci<<" "<<iter.getValue()<<std::endl;                             
                iter.setActiveState(true);                         
            }
        }
    }

    openvdb::Index64 result = distance_grid->activeVoxelCount();
    num_empty=unsigned(result);
}
