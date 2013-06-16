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

to run :
    > ./autopack <setup.xml> <seed> <forceSphere> <dsIngrGrid> <diplaysGrid> > ~/Dropbox/testCpp.py

    example : 
    > ./autopack  CellScape1.0.xml 12 0 1 1 > ~/Dropbox/testCpp.py

the command to run in python and visualize the result:
    > execfile("/Users/ludo/Dropbox/testCpp.py")
"""
*/

#include "Types.h"

/*
general include
*/
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <random>

/*
xml parser for collada import, dont forgot the lib in the folder
*/
#include <pugixml.hpp> 


#include "Ingredient.h"
#include "BigGrid.h"
#include "XMLLoader.h"

/*assimp include*/

//#include <assimp/assimp.h>// C importer interface
/* Deprecated use only pugixml
#include <assimp/assimp.hpp>// C++ importer interface
#include <assimp/aiPostProcess.h>// Post processing flags
#include <assimp/aiScene.h>// Output data structure
#include <assimp/DefaultLogger.h>
#include <assimp/LogStream.h>


// Create an instance of the Importer class
Assimp::Importer importer;
//const aiScene* scene;
*/

/*
a container // ie organelle
not use yet
*/

struct container{
    std::string name;
    int id;
    openvdb::FloatGrid::Ptr grid;
    openvdb::CoordBBox bbox;
    //why no use the mesh struct.
    std::vector<openvdb::Vec3s> vertices;
    std::vector<openvdb::Vec3s> normals;
    std::vector<openvdb::Vec3I > faces;
    std::vector<openvdb::Vec4I > quads;
    //list of indice too ?
};



float stepsize = float(15*1.1547);         //grid step size ie smallest ingredients radius
bool forceSphere = true;



/*
openvdb::CoordBBox getBB(float radius,  openvdb::Vec3f pos){
    const openvdb::Vec3d ibotleft(pos-radius,pos-radius,pos-radius);//(0,0,0)
    const openvdb::Vec3d iupright(pos+radius,pos+radius,pos+radius);//(1000,1000,10);
    openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(stepsize);
    openvdb::Vec3d botleft=transform->worldToIndex(ibotleft);
    openvdb::Vec3d upright=transform->worldToIndex(iupright);
    openvdb::Coord left(openvdb::tools::local_util::roundVec3(botleft));//(openvdb::Int32)botleft.x(),(openvdb::Int32)botleft.y(),(openvdb::Int32)botleft.z());
    openvdb::Coord right(openvdb::tools::local_util::roundVec3(upright));//(openvdb::Int32)upright.x(),(openvdb::Int32)upright.y(),(openvdb::Int32)upright.z());
    //define the active region that will be our boundary. set to max everywhere
    openvdb::CoordBBox bbox = openvdb::CoordBBox(left,right);//min and max
    return bbox;
}
*/


//exapl of function applied on grid data
struct SetMaxToDefault {
    float _max;
    SetMaxToDefault(float max) {_max=max;}
    inline void operator()(const openvdb::FloatGrid::ValueAllIter& iter) const {
       if (iter.getValue() == _max) iter.setValue(dmax);
    }
};




//code from the openvdb documentation
template<typename OpType>
void processTypedGrid(openvdb::GridBase::Ptr grid, OpType& op)
{
#define CALL_OP(GridType) \
    op.template operator()<GridType>(openvdb::gridPtrCast<GridType>(grid))

    if (grid->isType<openvdb::BoolGrid>())        CALL_OP(openvdb::BoolGrid);
    else if (grid->isType<openvdb::FloatGrid>())  CALL_OP(openvdb::FloatGrid);
    else if (grid->isType<openvdb::DoubleGrid>()) CALL_OP(openvdb::DoubleGrid);
    else if (grid->isType<openvdb::Int32Grid>())  CALL_OP(openvdb::Int32Grid);
    else if (grid->isType<openvdb::Int64Grid>())  CALL_OP(openvdb::Int64Grid);
    else if (grid->isType<openvdb::Vec3IGrid>())  CALL_OP(openvdb::Vec3IGrid);
    else if (grid->isType<openvdb::Vec3SGrid>())  CALL_OP(openvdb::Vec3SGrid);
    else if (grid->isType<openvdb::Vec3DGrid>())  CALL_OP(openvdb::Vec3DGrid);
    else if (grid->isType<openvdb::StringGrid>()) CALL_OP(openvdb::StringGrid);

#undef CALL_OP
}

// Define a functor that prunes the trees of grids of arbitrary type
// with a fixed pruning tolerance.
struct PruneOp {
    double tolerance;
    PruneOp(double t): tolerance(t) {}

    template<typename GridType>
    void operator()(typename GridType::Ptr grid) const
    {
        grid->tree().prune(typename GridType::ValueType(tolerance));
    }
    // Overload to handle string-valued grids
    void operator()(openvdb::StringGrid::Ptr grid) const
    {
        grid->tree().prune();
    }
};


inline openvdb::Coord getIJKc(int u,openvdb::Coord dim){
    // = {0,0,0};    
    int i_ijk[3];
    i_ijk[0]=0;
    i_ijk[1]=0;
    i_ijk[2]=0;    
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
    return openvdb::Coord(i_ijk[0],i_ijk[1],i_ijk[2]);   
}

void printIngredientGrid(Ingredient ingr){
    std::cout << "iname = \"" << ingr.name << "\"\n"; 
    std::cout << "inside=[]\n";
    int counter=0;
    //evalMinMax
    openvdb::Coord cc;
    for (openvdb::FloatGrid::ValueAllIter  iter = ingr.gsphere->beginValueAll(); iter; ++iter) {
        //create a sphere with color or radius dependant of value ?
        cc=iter.getCoord();
        if (ingr.bbox.isInside(cc)){
            counter++;
            float d = iter.getValue();
            //float dv = d;
            openvdb::Vec3f pos=ingr.gsphere->indexToWorld(cc); //getValue? 
            /*if (iter.getValue() > 0.0 ){
                float d = iter.getValue();
                if (d == dmax) d=90.0;
                std::cout << "distances.append( "<< d <<")\n";
                std::cout << "all.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
            }*/
            //if (iter.getValue() > 0.0 ){
            //    if (d == dmax) dv=g.maxradius;
            //    std::cout << "Distances.append( "<< dv <<")\n";
            //    std::cout << "all.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
           // }
            if (d < 0.0 ) {
            //std::cout << "#" <<  cc.x() <<"," <<cc.y()<<"," << cc.z() <<std::endl;
                std::cout << "inside.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
            //std::cout << "s=helper.Sphere(\"sphere" << counter << "\",parent=parent,pos=pos)[0]\n";
            }
            else if ((d == dmax )||(d >= spherewidth*stepsize)) {
                if (d == dmax) d=spherewidth*stepsize;
                //std::cout << "background.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
                //std::cout << "distances.append( "<< d <<")\n";
            }
            else {
                //std::cout << "outside.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
        }
        }
    }
    //or pointCloud?
    std::cout << "isph=helper.PointCloudObject(iname+\"_inside\",vertices = inside,materials=[[0,0,1]],parent = parentHider)\n";
}

void printTheGrid(big_grid g){
    std::cout << "inside=[]\n";
    std::cout << "outside=[]\n";
    std::cout << "background=[]\n";
    std::cout << "distances=[]\n";
    std::cout << "Distances=[]\n";
    std::cout << "all=[]\n";
    int counter=0;
    openvdb::Coord cc;
    //the main grid, extract some distance information from it
    for (openvdb::FloatGrid::ValueAllIter  iter = g.distance_grid->beginValueAll(); iter; ++iter) {//g.distance_grid
    //for (openvdb::FloatGrid::ValueAllIter  iter = g.ingredients[1].gsphere->beginValueAll(); iter; ++iter) {
        //create a sphere with color or radius dependant of value ?
        cc=iter.getCoord();
        if (g.bbox.isInside(cc)){
        //if (g.ingredients[1].bbox.isInside(cc)){
            counter++;
            float d = iter.getValue();
            //float dv = d;
            openvdb::Vec3f pos=g.distance_grid->indexToWorld(cc); //getValue? 
            /*if (iter.getValue() > 0.0 ){
                float d = iter.getValue();
                if (d == dmax) d=90.0;
                std::cout << "distances.append( "<< d <<")\n";
                std::cout << "all.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
            }*/
            //if (iter.getValue() > 0.0 ){
            //    if (d == dmax) dv=g.maxradius;
            //    std::cout << "Distances.append( "<< dv <<")\n";
            //    std::cout << "all.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
           // }
            if (d < 0.0 ) {
            //std::cout << "#" <<  cc.x() <<"," <<cc.y()<<"," << cc.z() <<std::endl;
                std::cout << "inside.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
            //std::cout << "s=helper.Sphere(\"sphere" << counter << "\",parent=parent,pos=pos)[0]\n";
            }
            else if ((d == dmax )||(d >= spherewidth*stepsize)) {
                if (d == dmax) d=spherewidth*stepsize;
                //std::cout << "background.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
                //std::cout << "distances.append( "<< d <<")\n";
            }
            else {
                //std::cout << "outside.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
        }
        }
    }
    std::cout << "isph=helper.PointCloudObject(\"inside\",vertices = inside,materials=[[0,0,1]])\n";
    std::cout << "isph=helper.PointCloudObject(\"background\",vertices = background)\n";
    std::cout << "isph=helper.PointCloudObject(\"outside\",vertices = outside,materials=[[0,1,0]])\n";
}

void generatePythonScript( big_grid &g, std::vector<float> &radiis, std::vector<openvdb::Vec3f> &colors, bool ds_grid, bool ds_ingrgrid ) 
{
    std::cout << "pts=[" << std::endl;
    openvdb::Vec3f pos;
    Ingredient * ingr;
    for(unsigned i = 0; i < g.rtrans.size(); ++i) { 
        ingr = g.results[i];
        openvdb::math::Transform::Ptr targetXform =
            openvdb::math::Transform::createLinearTransform();
        // Add the offset.
        targetXform->preMult(g.rrot[i]);
        targetXform->postTranslate(g.rtrans[i]);//should be woffset ? nope we apply on xyz not on ijk
        openvdb::math::Mat4d mat = targetXform->baseMap()->getAffineMap()->getMat4();
        for(openvdb::Vec3f const & position : g.results[i]->positions ) {        
            pos = mat.transform(position);
            std::cout << '[' << pos.x() <<',' << pos.y() << ',' << pos.z() << ']' << ',' <<std::endl; 
        }        
    }
    std::cout << "]" << std::endl; 
    std::cout << "matrices={}" << std::endl;
    for(unsigned i = 0; i < g.ingredients.size(); ++i) { 
        std::cout << "matrices[\"" << g.ingredients[i].name << "\"]=[]\n";
    }
    //openvdb::Vec3f pos;
    //sphere * ingr;
    for(unsigned i = 0; i < g.rtrans.size(); ++i) { 
        ingr = g.results[i];
        openvdb::math::Transform::Ptr targetXform =
            openvdb::math::Transform::createLinearTransform();
        // Add the offset.
        targetXform->preMult(g.rrot[i]);
        targetXform->postTranslate(g.rtrans[i]);
        openvdb::math::Mat4d mat = targetXform->baseMap()->getAffineMap()->getMat4();
        std::cout << "matrices[\""<< ingr->name <<"\"].append( " << mat  <<")"<< std::endl;
    }
    //std::cout << "]" << std::endl; 

    std::cout << "r=[" << std::endl  ;  
    for(std::vector<float>::size_type i = 0; i < radiis.size(); ++i) { 
        //std::cout << i << ' ' << radiis[i] << std::endl;
        for (std::vector<float>::size_type j = 0 ; j < g.results[i]->radii.size() ; j ++){
            std::cout << g.results[i]->radii[j] << ',' <<std::endl; //rot? 
        }        
        //std::cout << 3.0 << ',' <<std::endl; 
        //std::cout << 5.0 << ',' <<std::endl;  
    }    
    std::cout << "]" << std::endl;

    std::cout << "color=[" << std::endl  ;  
    for(std::vector<float>::size_type i = 0; i < colors.size(); ++i) { 
        //std::cout << i << ' ' << usedPoints[i] << std::endl;
        for (std::vector<float>::size_type j = 0 ; j < g.results[i]->radii.size() ; j ++){
            std::cout << '[' << colors[i].x() <<',' << colors[i].y() << ',' << colors[i].z() << ']' << ',' <<std::endl; 
        }
        //std::cout << '[' << colors[i].x() <<',' << colors[i].y() << ',' << colors[i].z() << ']' << ',' <<std::endl;  
    }  

    std::cout << "]" << std::endl;
    std::cout << "import upy\n" << "helper = upy.getHelperClass()()\n";
    std::cout << "pesph=helper.newEmpty(\"base_sphere\")\n";
    std::cout << "bsph=helper.Sphere(\"sphere\",parent=pesph)[0]\n";
    std::cout << "parent=helper.newEmpty(\"grid parent\")\n";
    std::cout << "parentHider=helper.newEmpty(\"hider parent\")\n";
    std::cout << "parentInstance=helper.newEmpty(\"instances parent\")\n";
    std::cout << "parentdistance=helper.newEmpty(\"grid distance\")\n";
    std::cout << "isph=helper.instancesSphere(\"cpp\",pts,r,pesph,color,None,parent=parent)\n";

    if (ds_grid) printTheGrid(g);

    std::cout << "#" << g.num_empty << std::endl;
    std::cout << "#" << g.distance_grid->activeVoxelCount() << std::endl;
    float mini=0.0;
    float maxi=0.0;
    g.distance_grid->evalMinMax(mini,maxi);
    std::cout << "#" << mini << " " <<maxi << std::endl;
    //mesh
    //foreac ingredient load the mesh?
    //if mesh
    //for each inredient get the original mesh and create instance for their different poisition in the grid
    for(unsigned i = 0; i < g.ingredients.size(); ++i) { //segmntaton fault here ?
        if (ds_ingrgrid) printIngredientGrid(g.ingredients[i]);
        //this will create the individual ingredient grid for debugging the voxelization
        //if (g.ingredients[i].filename.empty())
        //    continue;
        std::cout << "# ingr " << g.ingredients[i].radius << std::endl;
        std::cout << "iname = \"" << g.ingredients[i].name << "\"\n";    
        std::cout << "parent=helper.newEmpty(iname+\"_parent\",parent = parentHider)\n";
        std::cout << "helper.read(\"" << g.ingredients[i].filename << "\")\n"; 
        std::cout << "geom = helper.getObject(iname)\n";
        std::cout << "if geom is not None :\n";
        std::cout << "\thelper.rotateObj(geom,[0,math.pi/2.0,0.0])\n";//rotate on X in C4D
        std::cout << "\thelper.reParent(geom,parent)\n";
        std::cout << "\tiparent=helper.newEmpty(iname+\"_iparent\")\n";
        std::cout << "\taxis = " << g.ingredients[i].principalVector << "\n";
        std::cout << "\tipoly = helper.instancePolygon(iname+\"inst\",\n";
        std::cout << "                      matrices=matrices[iname],\n";
        std::cout << "                      mesh=parent,parent = iparent,\n";
        std::cout << "                      transpose= False,colors=["<<g.ingredients[i].color << "],\n";
        std::cout << "                      axis=axis)\n";
        //openvdb use convex - hull, art was suggesting decomposing the mesh to compensate
        //or just use the multisphere as in blood recipe. Using the mesh merging seems difficult because
        // of the transformation inside the collada file.
    } 
}

//main loop is here
int main(int argc, char* argv[])
{   
    clock_t beginRun = clock();
    const unsigned int seed = atoi(argv[2]);             //random seed
    forceSphere = atoi(argv[3])>0;          //force using sphere level set instead of mesh level set
    std::string filename = argv[1];         //xml setuo file
 
    bool ds_ingrgrid = atoi(argv[4])>0;//display level set grid for ingredient
    bool ds_grid = atoi(argv[5])>0;    //display global set grid
 
    srand (seed);
 
    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    openvdb::initialize();
    
    //load and setup the pack
    big_grid grid = load_xml(filename,0,seed);
    
    int counterRej=0;
    int rejection=0;

    openvdb::Coord cc;
    std::vector<openvdb::Vec3f> usedPoints;
    std::vector<float> radiis;
    std::vector<openvdb::Vec3f> colors;
    
    grid.prepareIngredient();

    int totalNumMols = grid.calculateTotalNumberMols();
    std::cout << "#prepare ingredient complete\n";

    int counter = 0;
    int PlacedMols=0;
    int emptyList;

    openvdb::FloatGrid::Accessor accessor_distance = grid.distance_grid->getAccessor();

    while(grid.num_empty > 0) {     
    //for(unsigned i = 0; i < 2 ; ++i) {//nx*ny*nz
        //pick the object to drop
        if (DEBUG) std::cout << "#begin loop\n";
        if (grid.activeIngr.size() == 0) {
                std::cout << "#broken by no more ingredient Done!!!****\n";
                break;
        }
        if (grid.vRangeStart>1.0){
                std::cout << "#broken by vRange and hence Done!!!****\n";
                break;
        }
        Ingredient* ingr = grid.pickIngredient();
        //sphere ingr= g.sample_ingredient();//sampl using distance information as well ?        
        if (!ingr->isActive()) {
            std::cout << "#ingredient not active, complete >= 1.0 \n";
            //g.dropIngredient(ingr);            
            //continue;
        }
        if (DEBUG) std::cout  << "#" << ingr->name << " c " << ingr->completion << " counter " << ingr->counter<<" nbmol " << ingr->nbMol <<std::endl;
        //unsigned s = g.sample_empty(); //get the next available point randomly
        //int s = g.getPointToDrop(ingr, ingr->radius,1.0);
        openvdb::Coord s = grid.getPointToDropCoord(ingr, ingr->minRadius,1.0,&emptyList);
        if (DEBUG) std::cout  << "#" << s << " " <<emptyList << " " << accessor_distance.getValue(s) << " " << ingr->radius <<std::endl;
        if ( emptyList == -1 ){
            continue;
        }
        //unsigned s = g.sample_empty(); 
        
        //if (mode == 1) s=g.sample_empty_distance(ingr.radius); //get the next available point randomly  
        //else if (mode == 0)  s=g.sample_empty();    
        //unsigned s = g.sample_closest_empty(); 
        //bool collision = g.try_drop(s,ingr);
        bool collision = grid.try_dropCoord(s,ingr);
        if (!collision){
            PlacedMols++;
            radiis.push_back(ingr->radius);
            colors.push_back(ingr->color);
            rejection=0;  
            ingr->counter++;
            ingr->completion = (float)ingr->counter/(float)ingr->nbMol;
            if (DEBUG) std::cout << "# main loop accepted c " <<ingr->completion<< " " << ingr->name << " counter " << ingr->counter << " ijk " << s << " nempty " << grid.num_empty << " r " << rejection <<  " collide " << collision << std::endl;          
            if (ingr->completion >= 1.0) {
                std::cout << "#ingredient completion no more nmol to place "<< ingr->counter << std::endl;
                grid.dropIngredient(ingr);            
            //continue;
                }
            }
        else {
            ingr->rejectionCounter++;
            if (ingr->rejectionCounter > ingr->rejectionThreshold){
                std::cout << "#ingredient completion too many rejection "<< ingr->rejectionCounter << std::endl;
                ingr->completion =1.0;//erase from list ingredient /
                grid.dropIngredient(ingr); 
            }
            grid.visited_rejected_coord.push_back(s);
            rejection++;  
            counterRej++; 
            if (DEBUG) std::cout << "# main loop rejected " << ingr->name << ' ' << s << ' ' << grid.num_empty << ' ' << rejection <<  " collide " << collision << std::endl;            
                 
        }  
        counter++;
        //g.set_filled(s); 
        if   (rejection > 2000) break; 
        //if   (counterRej > 100) break; 
        //if   (counter > 10) break;        
    }
    //The packing is done, we will generate a python script executed in a 3d Host for the vizualisation
    //this can be replace by any output.
    
    std::cout << "#distance_grid->activeVoxelCount() " << grid.num_empty << " " <<grid.distance_grid->activeVoxelCount()<<std::endl;
    std::cout << "#main loop " << grid.rtrans.size() << " on " << totalNumMols << std::endl;

    clock_t endRun = clock();
    std::cout << "#Running time: " << std::fixed << double(endRun-beginRun)/(60*1000) << std::defaultfloat << std::endl;
    generatePythonScript(grid, radiis, colors, ds_grid, ds_ingrgrid);


    //stdout the grid
}