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


/* the comparison function are strict translatio from the python code */
//sort function for ingredient//
//The value returned indicates whether the element passed as first argument 
//is considered to go before the second in the specific strict weak ordering it defines.
//can weuse template here ? so ca accept any ingredient type..
bool ingredient_compare1(Ingredient* x, Ingredient* y){
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
           if (c1 > c2) //# c1 > c2
               return true;
           else
               return false;
       }else
           return false;
    }else
       return false;
}

bool ingredient_compare0(Ingredient* x, Ingredient* y){
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
           if (c1 > c2) //# c1 > c2
               return true;
           else
               return false;
       }else
           return false;
    }else
       return false;
}

bool ingredient_compare2(Ingredient* x, Ingredient* y){
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
       if (c1 > c2) //# c1 > c2
           return true;
       else
           return false;
   }
   else
    return false;
}


//some expermental functions for the IJK<->U grid index convertion
//currently unused
inline unsigned getU(openvdb::Coord dim,openvdb::Coord ijk){
    return (int)(ijk.x()*dim.x()*dim.y() + ijk.y()*dim.x() + ijk.z());
}

inline void getIJK(int u,openvdb::Coord dim,int* i_ijk){
    // = {0,0,0};    
    //openvdb::Coord ijk(0,0,0);
    int nxnynz = dim.x()*dim.y()*dim.z();
    int nynz = dim.z()*dim.y();
    //int nx = dim.x();
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


} //namespace

big_grid::big_grid( float step, openvdb::Vec3d bot, openvdb::Vec3d up, unsigned seed ) :     
    distance_grid(initializeDistanceGrid(bot, up)),
    num_points(initializeNumPointsCount()),    
    num_empty(num_points),
    uniform(0.0f,1.0f),
    distribution(0.0,1.0),
    gauss(0.0f,0.3f),
    pickWeightedIngr(true),    
    pickRandPt(true)
{
    generator.seed(seed);
    space = step;
}

unsigned int big_grid::initializeNumPointsCount()
{
    dim = distance_grid->evalActiveVoxelDim();    
    std::cout << "#Grid Npoints " << dim << " " << num_points << " "  << distance_grid->activeVoxelCount() << std::endl;
    return dim.x()*dim.y()*dim.z();
}

openvdb::FloatGrid::Ptr big_grid::initializeDistanceGrid( openvdb::Vec3d bot, openvdb::Vec3d up )
{
    openvdb::FloatGrid::Ptr distance_grid;
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

    openvdb::FloatGrid::Accessor accessor_distance = distance_grid->getAccessor();
    std::cout << "#Testing distance access:" << std::endl;
    std::cout << "#Grid " << left << " "<< botleft << " = " << accessor_distance.getValue(left) << std::endl;
    std::cout << "#Grid " << right << " "<< upright << " = " << accessor_distance.getValue(right) << std::endl;
    std::cout << "#Grid " << bbox << std::endl;
    

    return distance_grid;
}

void big_grid::setIngredients( std::vector<Ingredient> const & _ingredients )
{
    //set the ingredients list to pack in the grid
    //retrieve the biggest one
    ingredients = _ingredients;
    activeIngr.resize(ingredients.size());
    float unit_volume = pow(stepsize,3.0);
    float grid_volume =  num_points*unit_volume;
    std::cout << "#Grid Volume " << grid_volume << " unitVol " << unit_volume << std::endl;
    for(unsigned i = 0; i < ingredients.size(); ++i) { 
        activeIngr[i] = &ingredients[i];
        ingredients[i].setCount(grid_volume);
    }    
}

void big_grid::getSortedActiveIngredients()
{
    //pirorities120 is ot used ??
    std::vector<Ingredient*> ingr1;  // given priorities pass ptr ?
    std::vector<float> priorities1;
    std::vector<Ingredient*> ingr2;  // priority = 0 or none and will be assigned based on complexity
    std::vector<float> priorities2;
    std::vector<Ingredient*> ingr0;  // negative values will pack first in order of abs[packingPriority]
    std::vector<float> priorities0;
    Ingredient* lowestIng; 
    Ingredient* ing; 
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

    double totalRadii = 0.0;
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

void big_grid::prepareIngredient()
{
    getSortedActiveIngredients();
    std::cout << "#len(allIngredients) " << ingredients.size() << std::endl;
    std::cout << "#len(activeIngr0) " << activeIngr0.size() << std::endl;
    std::cout << "#len(activeIngr12) " << activeIngr12.size() << std::endl;

    calculateThresholdAndNormalizedPriorities();
}

void big_grid::updatePriorities()
{
    vRangeStart = vRangeStart + normalizedPriorities[0];
                
    //# Start of massive overruling section from corrected thesis file of Sept. 25, 2012    
    prepareIngredient();

}

void big_grid::dropIngredient( Ingredient *ingr )
{
    std::cout << "#drop ingredient " << ingr->name << " " << ingr->nbMol << " " << ingr->counter << " "<< ingr->rejectionCounter <<std::endl;
    
    ingr->completion = 1.0;
    //update priorities will regenerate activeIngr based on completion field
    updatePriorities();
}

Ingredient* big_grid::pickIngredient()
{
    Ingredient* ingr;
    //float prob;
    //std::default_random_engine generator (seed);
    //std::default_random_engine generator(seed);
    //std::uniform_real_distribution<float> distribution(0.0,1.0);
    //std::normal_distribution<float> distribution(0.0,1.0);        
    unsigned int ingrInd;
    float threshProb;
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

            if (ingrInd <  activeIngr.size())
                ingr = activeIngr[ingrInd];
            else {
                std::cout << "#error in histoVol pick Ingredient "  << ingrInd << ' ' <<activeIngr.size() << std::endl;
                ingr = activeIngr[0];
            }
        }
    }else {
        //pick random
        //ingr = activeIngr[rand() % numActiveIngr];
        ingr = activeIngr[(int) (uniform(generator) * activeIngr.size())];
        //uniform(generator)
    }
    return ingr;
}

openvdb::Coord big_grid::getPointToDropCoord( Ingredient* ingr, float radius,float jitter,int *emptyList )
{
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
    for (openvdb::FloatGrid::ValueOnIter  iter = distance_grid->beginValueOn(); iter; ++iter) {
        //for (openvdb::FloatGrid::ValueAllIter  iter = distance_grid->beginValueAll(); iter; ++iter) {
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
        vRangeStart = vRangeStart + normalizedPriorities[0];
        std::cout << "# drop no more point \n" ;
        dropIngredient(ingr); 
        totalPriorities = 0; //# 0.00001
        *emptyList = 1;
        return openvdb::Coord(0,0,0);                   
    }
    if (pickRandPt){
        std::cout << "#allIngrPts " <<allIngrPts.size() << "mode " << ingr->packingMode <<std::endl;
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
    rotMatj.setToRotation(openvdb::math::Vec3d(rand(),rand(),rand()),rand()*M_PI); // random value for axe and angle in radians
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

            // Save copies of the two grids; compositing operations always
            // modify the A grid and leave the B grid empty.
            if (DEBUG) std::cout << "#duplicate ingredient grid "<< std::endl;
            openvdb::FloatGrid::Ptr copyOfGridSphere = openvdb::FloatGrid::create(dmax);
            //openvdb::FloatGrid::Ptr copyOfGridSphere = ingr.gsphere->deepCopy();
            copyOfGridSphere->setTransform(targetXform);

            /*
            openvdb::Mat4R xform =
                sourceXform->baseMap()->getAffineMap()->getMat4() *
                targetXform->baseMap()->getAffineMap()->getMat4().inverse();
            */

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

            openvdb::Index64 result = distance_grid->activeVoxelCount();
            assert( result <= std::numeric_limits<unsigned int>::max() );
            num_empty=unsigned(result);

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


bool big_grid::checkSphCollisions( point pos,openvdb::math::Mat4d rotMatj, float radii, Ingredient* sp )
{
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
    //can we first test if boudning box hasOverlap ?
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

void big_grid::calculateThresholdAndNormalizedPriorities()
{
    normalizedPriorities.clear();
    thresholdPriorities.clear();
    //# Graham- Once negatives are used, if picked random# 
    //# is below a number in this list, that item becomes 
    //#the active ingredient in the while loop below
    normalizedPriorities.resize(activeIngr0.size(), 0.0);
    if (pickWeightedIngr)
            thresholdPriorities.resize(activeIngr0.size(), 2.0);   

    float totalPriorities = std::accumulate(activeIngr12.begin(), activeIngr12.end(), 0.0, 
        [](float accumulator , Ingredient* ingr) { return accumulator + ingr->packingPriority; });
    std::cout << "#totalPriorities " << totalPriorities << std::endl;
    
    if (totalPriorities == 0)
    {
        normalizedPriorities.resize(normalizedPriorities.size() + activeIngr12.size(), 0.0);
        thresholdPriorities.resize(thresholdPriorities.size() + activeIngr12.size(), 0.0);  
    }
    else
    {
        float previousThresh = 0;                
        for(Ingredient * ingr : activeIngr12) {
            const float np = ingr->packingPriority/totalPriorities;            
            if (DEBUG)  std::cout << "#np is "<< np << " pp is "<< ingr->packingPriority << " tp is " << np + previousThresh << std::endl;

            normalizedPriorities.push_back(np);
            previousThresh += np;
            thresholdPriorities.push_back(previousThresh);
        }    
    }

    activeIngr = activeIngr0;
    activeIngr.insert(activeIngr.end(), activeIngr12.begin(), activeIngr12.end());
}

int big_grid::calculateTotalNumberMols()
{
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
            //float threshProb = thresholdPriorities[i];
            Ingredient* ing = activeIngr[nls];
            std::cout << "#nmol Fill5else is for ingredient : "<< ingredients[i].name<< ' '  << ingredients[i].nbMol<< std::endl ;
            totalNumMols += ing->nbMol;
            nls++;
        }
        std::cout << "#totalNumMols Fill5else = " << totalNumMols << std::endl;
    }
    return totalNumMols;
}
