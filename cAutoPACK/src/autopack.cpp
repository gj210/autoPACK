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

#include <boost/program_options.hpp>

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
    openvdb::DoubleGrid::Ptr grid;
    openvdb::CoordBBox bbox;
    //why no use the mesh struct.
    std::vector<openvdb::Vec3s> vertices;
    std::vector<openvdb::Vec3s> normals;
    std::vector<openvdb::Vec3I > faces;
    std::vector<openvdb::Vec4I > quads;
    //list of indice too ?
};



double stepsize = double(15*1.1547);         //grid step size ie smallest ingredients radius



/*
openvdb::CoordBBox getBB(double radius,  openvdb::Vec3d pos){
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
    double _max;
    SetMaxToDefault(double max) {_max=max;}
    inline void operator()(const openvdb::DoubleGrid::ValueAllIter& iter) const {
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
    else if (grid->isType<openvdb::DoubleGrid>())  CALL_OP(openvdb::DoubleGrid);
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
    return openvdb::Coord(i_ijk[0],i_ijk[1],i_ijk[2]);   
}

void printIngredientGrid(std::ofstream & resultFile, Ingredient const& ingr){
    resultFile << "iname = \"" << ingr.name << "\"\n"; 
    resultFile << "inside=[]\n";
    int counter=0;
    //evalMinMax
    openvdb::Coord cc;
    for (openvdb::DoubleGrid::ValueAllIter  iter = ingr.gsphere->beginValueAll(); iter; ++iter) {
        //create a sphere with color or radius dependant of value ?
        cc=iter.getCoord();
        if (ingr.bbox.isInside(cc)){
            counter++;
            double d = iter.getValue();
            //double dv = d;
            openvdb::Vec3d pos=ingr.gsphere->indexToWorld(cc); //getValue? 
            /*if (iter.getValue() > 0.0 ){
                double d = iter.getValue();
                if (d == dmax) d=90.0;
                resultFile << "distances.append( "<< d <<")\n";
                resultFile << "all.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
            }*/
            //if (iter.getValue() > 0.0 ){
            //    if (d == dmax) dv=g->maxradius;
            //    resultFile << "Distances.append( "<< dv <<")\n";
            //    resultFile << "all.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
           // }
            if (d < 0.0 ) {
            //resultFile << "#" <<  cc.x() <<"," <<cc.y()<<"," << cc.z() <<std::endl;
                resultFile << "inside.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
            //resultFile << "s=helper.Sphere(\"sphere" << counter << "\",parent=parent,pos=pos)[0]\n";
            }
            else if ((d == dmax )||(d >= spherewidth*stepsize)) {
                if (d == dmax) d=spherewidth*stepsize;
                //resultFile << "background.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
                //resultFile << "distances.append( "<< d <<")\n";
            }
            else {
                //resultFile << "outside.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
        }
        }
    }
    //or pointCloud?
    resultFile << "isph=helper.PointCloudObject(iname+\"_inside\",vertices = inside,materials=[[0,0,1]],parent = parentHider)\n";
}

void printTheGrid(std::ofstream & resultFile, std::shared_ptr<big_grid> g){
    resultFile << "inside=[]\n";
    resultFile << "outside=[]\n";
    resultFile << "background=[]\n";
    resultFile << "distances=[]\n";
    resultFile << "Distances=[]\n";
    resultFile << "all=[]\n";
    int counter=0;
    openvdb::Coord cc;
    //the main grid, extract some distance information from it
    for (openvdb::DoubleGrid::ValueAllIter  iter = g->distance_grid->beginValueAll(); iter; ++iter) {//g->distance_grid
    //for (openvdb::DoubleGrid::ValueAllIter  iter = g->ingredients[1].gsphere->beginValueAll(); iter; ++iter) {
        //create a sphere with color or radius dependant of value ?
        cc=iter.getCoord();
        if (g->bbox.isInside(cc)){
        //if (g->ingredients[1].bbox.isInside(cc)){
            counter++;
            double d = iter.getValue();
            //double dv = d;
            openvdb::Vec3d pos=g->distance_grid->indexToWorld(cc); //getValue? 
            /*if (iter.getValue() > 0.0 ){
                double d = iter.getValue();
                if (d == dmax) d=90.0;
                resultFile << "distances.append( "<< d <<")\n";
                resultFile << "all.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
            }*/
            //if (iter.getValue() > 0.0 ){
            //    if (d == dmax) dv=g->maxradius;
            //    resultFile << "Distances.append( "<< dv <<")\n";
            //    resultFile << "all.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
           // }
            if (d < 0.0 ) {
            //resultFile << "#" <<  cc.x() <<"," <<cc.y()<<"," << cc.z() <<std::endl;
                resultFile << "inside.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
            //resultFile << "s=helper.Sphere(\"sphere" << counter << "\",parent=parent,pos=pos)[0]\n";
            }
            else if ((d == dmax )||(d >= spherewidth*stepsize)) {
                if (d == dmax) d=spherewidth*stepsize;
                resultFile << "background.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
                //resultFile << "distances.append( "<< d <<")\n";
            }
            else {
                resultFile << "outside.append( ["<<pos.x()<<","<<pos.y()<<"," << pos.z() <<"])\n";
        }
        }
    }
    resultFile << "isph=helper.PointCloudObject(\"inside\",vertices = inside,materials=[[0,0,1]])\n";
    resultFile << "isph=helper.PointCloudObject(\"background\",vertices = background)\n";
    resultFile << "isph=helper.PointCloudObject(\"outside\",vertices = outside,materials=[[0,1,0]])\n";
}

void generatePythonScript( std::string file, std::shared_ptr<big_grid> &g, std::vector<double> &radiis, std::vector<openvdb::Vec3d> &colors, bool ds_grid, bool ds_ingrgrid ) 
{
    std::ofstream resultFile;
    resultFile.open (file);

    resultFile << "pts=[" << std::endl;
    openvdb::Vec3d pos;
    Ingredient * ingr;
    for(unsigned i = 0; i < g->rtrans.size(); ++i) { 
        ingr = g->results[i];
        openvdb::math::Transform::Ptr targetXform =
            openvdb::math::Transform::createLinearTransform();
        // Add the offset.
        targetXform->preMult(g->rrot[i]);
        targetXform->postTranslate(g->rtrans[i]);//should be woffset ? nope we apply on xyz not on ijk
        openvdb::math::Mat4d mat = targetXform->baseMap()->getAffineMap()->getMat4();
        for(openvdb::Vec3d const & position : g->results[i]->positions ) {        
            pos = mat.transform(position);
            resultFile << '[' << pos.x() <<',' << pos.y() << ',' << pos.z() << ']' << ',' <<std::endl; 
        }        
    }
    resultFile << "]" << std::endl; 
    resultFile << "matrices={}" << std::endl;
    for(Ingredient const & ingr: g->ingredientsDipatcher) { 
        resultFile << "matrices[\"" << ingr.name << "\"]=[]\n";
    }
    //openvdb::Vec3d pos;
    //sphere * ingr;
    for(unsigned i = 0; i < g->rtrans.size(); ++i) { 
        ingr = g->results[i];
        openvdb::math::Transform::Ptr targetXform =
            openvdb::math::Transform::createLinearTransform();
        // Add the offset.
        targetXform->preMult(g->rrot[i]);
        targetXform->postTranslate(g->rtrans[i]);
        openvdb::math::Mat4d mat = targetXform->baseMap()->getAffineMap()->getMat4();
        resultFile << "matrices[\""<< ingr->name <<"\"].append( " << mat  <<")"<< std::endl;
    }
    //resultFile << "]" << std::endl; 

    resultFile << "r=[" << std::endl  ;  
    for(std::vector<double>::size_type i = 0; i < radiis.size(); ++i) { 
        //resultFile << i << ' ' << radiis[i] << std::endl;
        for (std::vector<double>::size_type j = 0 ; j < g->results[i]->radii.size() ; j ++){
            resultFile << g->results[i]->radii[j] << ',' <<std::endl; //rot? 
        }        
        //resultFile << 3.0 << ',' <<std::endl; 
        //resultFile << 5.0 << ',' <<std::endl;  
    }    
    resultFile << "]" << std::endl;

    resultFile << "color=[" << std::endl  ;  
    for(std::vector<double>::size_type i = 0; i < colors.size(); ++i) { 
        //resultFile << i << ' ' << usedPoints[i] << std::endl;
        for (std::vector<double>::size_type j = 0 ; j < g->results[i]->radii.size() ; j ++){
            resultFile << '[' << colors[i].x() <<',' << colors[i].y() << ',' << colors[i].z() << ']' << ',' <<std::endl; 
        }
        //resultFile << '[' << colors[i].x() <<',' << colors[i].y() << ',' << colors[i].z() << ']' << ',' <<std::endl;  
    }  

    resultFile << "]" << std::endl;
    resultFile << "import upy\n" << "helper = upy.getHelperClass()()\n";
    resultFile << "pesph=helper.newEmpty(\"base_sphere\")\n";
    resultFile << "bsph=helper.Sphere(\"sphere\",parent=pesph)[0]\n";
    resultFile << "parent=helper.newEmpty(\"grid parent\")\n";
    resultFile << "parentHider=helper.newEmpty(\"hider parent\")\n";
    resultFile << "parentInstance=helper.newEmpty(\"instances parent\")\n";
    resultFile << "parentdistance=helper.newEmpty(\"grid distance\")\n";
    resultFile << "isph=helper.instancesSphere(\"cpp\",pts,r,pesph,color,None,parent=parent)\n";

    if (ds_grid) printTheGrid(resultFile, g);

    resultFile << "#" << g->num_empty << std::endl;
    resultFile << "#" << g->distance_grid->activeVoxelCount() << std::endl;
    double mini=0.0;
    double maxi=0.0;
    g->distance_grid->evalMinMax(mini,maxi);
    resultFile << "#" << mini << " " <<maxi << std::endl;
    //mesh
    //foreac ingredient load the mesh?
    //if mesh
    //for each inredient get the original mesh and create instance for their different poisition in the grid
    for(Ingredient const & ingr: g->ingredientsDipatcher) {
        if (ds_ingrgrid) printIngredientGrid(resultFile, ingr);
        //this will create the individual ingredient grid for debugging the voxelization
        //if (g->ingredients[i].filename.empty())
        //    continue;
        resultFile << "# ingr " << ingr.radius << std::endl;
        resultFile << "iname = \"" << ingr.name << "\"\n";    
        resultFile << "parent=helper.newEmpty(iname+\"_parent\",parent = parentHider)\n";
        resultFile << "helper.read(\"" << ingr.filename << "\")\n"; 
        resultFile << "geom = helper.getObject(iname)\n";
        resultFile << "if geom is not None :\n";
        resultFile << "\thelper.rotateObj(geom,[0,math.pi/2.0,0.0])\n";//rotate on X in C4D
        resultFile << "\thelper.reParent(geom,parent)\n";
        resultFile << "\tiparent=helper.newEmpty(iname+\"_iparent\")\n";
        resultFile << "\taxis = " << ingr.principalVector << "\n";
        resultFile << "\tipoly = helper.instancePolygon(iname+\"inst\",\n";
        resultFile << "                      matrices=matrices[iname],\n";
        resultFile << "                      mesh=parent,parent = iparent,\n";
        resultFile << "                      transpose= False,colors=["<<ingr.color << "],\n";
        resultFile << "                      axis=axis)\n";
        //openvdb use convex - hull, art was suggesting decomposing the mesh to compensate
        //or just use the multisphere as in blood recipe. Using the mesh merging seems difficult because
        // of the transformation inside the collada file.
    } 
}


void generateCoordFile( std::shared_ptr<big_grid> &g, std::vector<double> &radiis, std::vector<openvdb::Vec3d> &colors, bool ds_grid, bool ds_ingrgrid, const std::string &file ) 
{
    std::ofstream coordFile;
    coordFile.open (file);
    
    openvdb::Vec3d pos;
    Ingredient * ingr;
    for(unsigned i = 0; i < g->rtrans.size(); ++i) { 
        ingr = g->results[i];
        openvdb::math::Transform::Ptr targetXform =
            openvdb::math::Transform::createLinearTransform();
        // Add the offset.
        targetXform->preMult(g->rrot[i]);
        targetXform->postTranslate(g->rtrans[i]);//should be woffset ? nope we apply on xyz not on ijk
        openvdb::math::Mat4d mat = targetXform->baseMap()->getAffineMap()->getMat4();
        int j = 0;
        for(openvdb::Vec3d const & position : g->results[i]->positions ) {        
            pos = mat.transform(position);
            coordFile << pos.x() << ' ' << pos.y() << ' ' << pos.z() << ' ' << g->results[i]->radii[j] << std::endl;
            j++;
        }        
    }
    coordFile.close();
}

//main loop is here
int main(int argc, char* argv[])
{   
    namespace po = boost::program_options;
    po::options_description desc("Usage:\n\tautopack -i <xml recipe> <options>\n\n Options:");
    desc.add_options()
        ("help", "produce help message")
        ("recipe,i", po::value<std::string>(), "(required) Recipe XML file")
        ("output,o", po::value<std::string>()->default_value("results.py"), "Recipe XML file")
        ("randomSeed", po::value<unsigned int>()->default_value(12), "set random seed value")
        ("forceSphere", po::value<bool>()->default_value(false), "force using sphere level set instead of mesh level set")        
        ("format,f", po::value<std::string>()->default_value("PY"), "Output file format (PY | TXT)")
        ("level,l", po::value<bool>()->default_value(true), "display level set grid for ingredient")
        ("grid,g", po::value<bool>()->default_value(true), "display global set grid")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
    
    if (vm.count("help") || vm.count("recipe")==0) {
        std::cout << desc << "\n";
        return 1;
    }    

    clock_t beginRun = clock();
    const unsigned int seed = vm["randomSeed"].as<unsigned int>();
    const bool forceSphere = vm["forceSphere"].as<bool>();
    const std::string filename = vm["recipe"].as<std::string>();
 
    bool ds_ingrgrid = vm["level"].as<bool>();
    bool ds_grid = vm["grid"].as<bool>();
    const bool createResFile =  vm["format"].as<std::string>() == "TXT";
    const std::string outputFile (vm["output"].as<std::string>());
    srand (seed);
 
    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    openvdb::initialize();
    
    //load and setup the pack
    std::shared_ptr<big_grid> grid = load_xml(filename,0,seed, forceSphere);
    
    int rejection=0;

    openvdb::Coord cc;
    std::vector<openvdb::Vec3d> usedPoints;
    std::vector<double> radiis;
    std::vector<openvdb::Vec3d> colors;

    std::cout << "#prepare ingredient complete\n";

    int emptyList;

    openvdb::DoubleGrid::Accessor accessor_distance = grid->distance_grid->getAccessor();

    while(grid->num_empty > 0 && grid->ingredientsDipatcher.hasNextIngradient()) {     
        if (DEBUG) std::cout << "#begin loop\n";

        Ingredient* ingr = grid->ingredientsDipatcher.pickIngredient();
        if (!ingr->isActive()) {
            std::cout << "#ingredient not active, complete >= 1.0 \n";
        }

        if (DEBUG) std::cout  << "#" << ingr->name << " c " << ingr->completion << " counter " << ingr->counter<<" nbmol " << ingr->nbMol <<std::endl;

        openvdb::Coord s = grid->getPointToDropCoord(ingr, ingr->minRadius,1.0,emptyList);
        if (DEBUG) std::cout  << "#" << s << " " <<emptyList << " " << accessor_distance.getValue(s) << " " << ingr->radius <<std::endl;
        if ( emptyList == -1 ){
            continue;
        }

        bool collision = grid->try_dropCoord(s,ingr);
        if (!collision){
            radiis.push_back(ingr->radius);
            colors.push_back(ingr->color);
            rejection=0;  
            ingr->counter++;
            ingr->completion = (double)ingr->counter/(double)ingr->nbMol;
            if (DEBUG) std::cout << "# main loop accepted c " <<ingr->completion<< " " << ingr->name << " counter " << ingr->counter << " ijk " << s << " nempty " << grid->num_empty << " r " << rejection <<  " collide " << collision << std::endl;          
            if (ingr->completion >= 1.0) {
                std::cout << "#ingredient completion no more nmol to place "<< ingr->counter << std::endl;
                grid->ingredientsDipatcher.dropIngredient(ingr);            
            }
            std::cout << "## Ingredient placed, count:" << grid->rtrans.size() << std::endl;
        }            
        else {
            ingr->rejectionCounter++;
            if (ingr->rejectionCounter > ingr->rejectionThreshold){
                std::cout << "#ingredient completion too many rejection "<< ingr->rejectionCounter << std::endl;
                grid->ingredientsDipatcher.dropIngredient(ingr); 
            }
            ingr->visited_rejected_coord.push_back(s);
            rejection++;  
            std::cout << "## Ingredient rejected, already placed: " << grid->rtrans.size() << std::endl;
            if (DEBUG) std::cout << "# main loop rejected " << ingr->name << ' ' << s << ' ' << grid->num_empty << ' ' << rejection <<  " collide " << collision << std::endl;            
                 
        }  
        if   (rejection > 2000) break; 
    }
    
    //The packing is done, we will generate a python script executed in a 3d Host for the vizualisation
    //this can be replace by any output.
    
    std::cout << "#distance_grid->activeVoxelCount() " << grid->num_empty << " " <<grid->distance_grid->activeVoxelCount()<<std::endl;
    std::cout << "#main loop " << grid->rtrans.size() << std::endl;

    clock_t endRun = clock();
    std::cout << "#Running time: " << std::fixed << double(endRun-beginRun)/(60*1000) << std::defaultfloat << std::endl;
    
    if(createResFile)
        generateCoordFile( grid, radiis, colors, ds_grid, ds_ingrgrid, outputFile);
    else
        generatePythonScript(outputFile, grid, radiis, colors, ds_grid, ds_ingrgrid);

    //stdout the grid
}