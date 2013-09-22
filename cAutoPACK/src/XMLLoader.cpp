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

#include "Types.h"

#include <vector>
#include <map>
#include <pugixml.hpp>

#include "BigGrid.h"
#include "Ingredient.h"

#include "IngredientsFactory.h"

/* XML CODE */
/* 
parsing information from the autopack setup file as well as the collada mesh file 
*/

double getRadii(std::string str){
    //std::string str(input);  
    std::replace(str.begin(), str.end(), '[', ' ');  
    std::replace(str.begin(), str.end(), ']', ' ');  
    std::replace(str.begin(), str.end(), '(', ' ');  
    std::replace(str.begin(), str.end(), ')', ' ');  
    return atof(str.c_str());
}

std::vector<double> getBox(std::string str){
    std::replace(str.begin(), str.end(), '[', ' ');  
    std::replace(str.begin(), str.end(), ']', ' ');  
    std::replace(str.begin(), str.end(), '(', ' ');  
    std::replace(str.begin(), str.end(), ')', ' ');  
    std::replace(str.begin(), str.end(), ',', ' ');
    // Build an istream that holds the input string
    std::istringstream iss(str);
    //std::cout << str << std::endl;
    // If possible, always prefer std::vector to naked array
    std::vector<double> v;
    // Iterate over the istream, using >> to grab floats
    // and push_back to store them in the vector
    std::copy(std::istream_iterator<double>(iss),
        std::istream_iterator<double>(),
        std::back_inserter(v));
    return v;
}

openvdb::Vec3d moveMultisphereToGeomCenter( std::vector<openvdb::Vec3d> &positions )
{
    const auto geometricCenterOffset = Geometric::calculateGeometricCenter(std::begin(positions), std::end(positions));

    std::transform(positions.begin(), positions.end(), positions.begin(),
        [&geometricCenterOffset] (openvdb::Vec3d const& item) { return item - geometricCenterOffset; } );

    return geometricCenterOffset;
}


std::vector<openvdb::Vec3d> getPositions(std::vector<double> pos){
    std::vector<double>::size_type i = 0; 
    std::vector<openvdb::Vec3d> positions;   
    while (i<pos.size()-2){
         openvdb::Vec3d p(pos[i],pos[i+1],pos[i+2]);
         positions.push_back(p);
         i = i+3;
    }

    return positions;
}


openvdb::Vec3d getArray(std::string str){
    std::replace(str.begin(), str.end(), '[', ' ');  
    std::replace(str.begin(), str.end(), ']', ' ');  
    std::replace(str.begin(), str.end(), '(', ' ');  
    std::replace(str.begin(), str.end(), ')', ' ');  
    std::replace(str.begin(), str.end(), ',', ' ');
    // Build an istream that holds the input string
    std::istringstream iss(str);
    //std::cout << str << std::endl;
    // If possible, always prefer std::vector to naked array
    std::vector<double> v;
    // Iterate over the istream, using >> to grab floats
    // and push_back to store them in the vector
    std::copy(std::istream_iterator<double>(iss),
        std::istream_iterator<double>(),
        std::back_inserter(v));
    openvdb::Vec3d p(v[0],v[1],v[2]);
    return p;
}

//we load the autopack xml setup and create the grid class object
std::shared_ptr<big_grid> load_xml(std::string path,int _mode,unsigned _seed, bool forceSphere){
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(path.c_str());
    //std::cout << path << std::endl;
    std::cout << "#Load result: " << result.description() << ", AutoFillSetup name: " << doc.child("AutoFillSetup").attribute("name").value() << std::endl;
    //get the option
    const double smallestObjectSize = atof(doc.child("AutoFillSetup").child("options").attribute("smallestProteinSize").value());
    const unsigned seed = _seed;
    //const int mode = _mode; 
    stepsize = smallestObjectSize * 1.1547;
    std::cout <<"#step " << stepsize <<std::endl;
    std::cout <<"#seed " <<seed<<std::endl;
    std::string strbb(doc.child("AutoFillSetup").child("options").attribute("boundingBox").value());
    std::vector<double> bb = getBox(strbb);

    bool pickWeightedIngr =doc.child("AutoFillSetup").child("options").attribute("pickWeightedIngr").as_bool();
    bool pickRandPt =doc.child("AutoFillSetup").child("options").attribute("pickRandPt").as_bool();

    //only cytoplsme for now, should parse, organelle and gradient as well
    //could we use openvdb do compute compute/prepare the gradient?
    std::vector<Ingredient> _ingredients;
    //need to add organelle as well...organelle ingredient are inside organelle levelSet.
    pugi::xml_node cytoplasme = doc.child("AutoFillSetup").child("cytoplasme");
    //double radius, int mode, double concentration, 
    //     double packingPriority,int nbMol,std::string name, point color    
    for (pugi::xml_node ingredient = cytoplasme.child("ingredient"); ingredient; ingredient = ingredient.next_sibling("ingredient")){
        if (DEBUG) std::cout << "#ingredient " << ingredient.attribute("name").value() << std::endl;
        std::string str(ingredient.attribute("radii").value()); 
        std::vector<double> radii = getBox(str);     
        std::string posstr(ingredient.attribute("positions").value()); 
        std::vector<double> pos = getBox(posstr);   
        
        std::vector<openvdb::Vec3d> positions = getPositions(pos);
        auto geometricCenter = moveMultisphereToGeomCenter(positions);
        
        //double r = getRadii(str); //different radii... as well as the position...  
        //std::cout << r << std::endl;
        double mol = ingredient.attribute("molarity").as_float();
        if (DEBUG)std::cout << "#molarity "<<mol << std::endl;
        double priority = ingredient.attribute("packingPriority").as_float();
        if (DEBUG)std::cout << "#priority "<< priority << std::endl;        
        int nMol = ingredient.attribute("nbMol").as_int();
        if (DEBUG)std::cout << "#nmol "<< nMol << std::endl; 
        std::string iname(ingredient.attribute("name").value()); 
        openvdb::Vec3d color(1,0,0);
        if (ingredient.attribute("color")){        
            std::string strcol(ingredient.attribute("color").value());
            openvdb::Vec3d color = getArray(strcol);//TODO : pass to inredient
        }
        unsigned nbJitter = ingredient.attribute("nbJitter").as_int();
        if (DEBUG)std::cout << "#color "<< color << std::endl;
        std::string strjitter(ingredient.attribute("jitterMax").value());
        openvdb::Vec3d jitter =  getArray(strjitter);
        if (DEBUG)std::cout << "#jitter "<< jitter << std::endl;     
        std::string straxe(ingredient.attribute("principalVector").value());
        openvdb::Vec3d principalVector =  getArray(straxe);
        if (DEBUG)std::cout << "#principalVector "<< principalVector << std::endl;
        bool fSphere=forceSphere;
        if (ingredient.attribute("forceSphere"))
            fSphere = ingredient.attribute("forceSphere").as_bool(); 
        //also need packingMode,perturbAxisAmplitude
        //mesh file ... should use multiSphere function or makeMesh  
        //sphere ingr = makeSphere(r,_mode,mol,priority,nMol,iname,
        //                        color,nbJitter,jitter);
        //get the meshFile and get vertices+faces
        std::string strmeshFile(ingredient.attribute("meshFile").value());
        //if dae can use the xmlparser
        std::cout << "# mesh file ? " << strmeshFile << " x " << strmeshFile.empty() << std::endl;
        //can be none
        mesh mesh3d;
        std::vector<mesh> meshs;
        Ingredient ingr ;
        if ((!strmeshFile.empty())&&(!fSphere)){
            assert(!"Meshes are not  supported.");
        }
        else {
            ingr = makeMultiSpheres(radii,_mode,mol,priority,nMol,iname,
                    color,nbJitter,jitter,positions);
        }
        ingr.filename = strmeshFile;
        ingr.geometricCenter = geometricCenter;
        ingr.principalVector=principalVector;
        ingr.useRotAxis = false;
        if (ingredient.attribute("useRotAxis")){
            ingr.useRotAxis = ingredient.attribute("useRotAxis").as_bool();
            ingr.rotRange = ingredient.attribute("rotRange").as_double();
            if (ingr.useRotAxis){
                std::string strRaxe(ingredient.attribute("rotAxis").value());
                ingr.rotAxis =  getArray(strRaxe);
            }
        }
        ingr.perturbAxisAmplitude = 0.1;
        if (ingredient.attribute("perturbAxisAmplitude")){
            ingr.perturbAxisAmplitude = ingredient.attribute("perturbAxisAmplitude").as_double();
        }
        //packing mode overwrite from xml file
        ingr.packingMode = std::string(ingredient.attribute("packingMode").value());
        //rejectionThreshold        
        if (ingredient.attribute("rejectionThreshold")){      
            ingr.rejectionThreshold = ingredient.attribute("rejectionThreshold").as_int();
        }
        _ingredients.push_back(ingr);
        std::cout << "#ingredient done!" << std::endl;
    }
    
    //make the grid and return it
    const openvdb::Vec3d ibotleft(bb[0],bb[1],bb[2]);
    const openvdb::Vec3d iupright(bb[3],bb[4],bb[5]);    
    if (DEBUG)std::cout << "#box "<<ibotleft<<" "<<iupright << std::endl;

    //should create the grid from a xml file...easier setup    
    std::shared_ptr<big_grid> grid(new big_grid(_ingredients,stepsize,ibotleft,iupright,seed));
    grid->ingredientsDipatcher.pickWeightedIngr=pickWeightedIngr;
    grid->pickRandPt =pickRandPt; 
    //g.setMode(0);//1-close packing 0-random    
    return grid;
}