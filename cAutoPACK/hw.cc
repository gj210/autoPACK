/* 
###############################################################################
#
# autoPACK Authors: Graham T. Johnson, Mostafa Al-Alusi, Ludovic Autin, Michel Sanner
#   Based on COFFEE Script developed by Graham Johnson between 2005 and 2010 
#   with assistance from Mostafa Al-Alusi in 2009 and periodic input 
#   from Arthur Olson's Molecular Graphics Lab
#
# hw.cc Authors: Ludovic Autin
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
"""
*/

#include <vector>
#include <cstdlib>
#include <iostream>
#include <math.h> 
#include <algorithm>
#include <random>
#include <ctime>
#include <iterator>
#include <cassert>
#include <sstream>

#include "./pugixml-1.2/src/pugixml.hpp" 

#include <openvdb/openvdb.h>


/* a point definition
used for building the grid coordinates 
*/
struct point {
    float x,y,z;
};

/* 
a box, let define the grid size in x y an z
*/
struct box {
    unsigned x,y,z;
};

/* 
a SingleSphere ingredient with the radius and the packing mode
*/
struct sphere {
    float radius;
    float minRadius;
    //point trans;
    int mode;
    int nbMol;
    float completion;
    float packingPriority;
    std::string name;
    std::string packingMode="random";
    int counter;
    int rejectionCounter;
    int rejectionThreshold;
    unsigned nbJitter;
    bool active;
    point color;
    point trans;
    point jitterMax;
};

//sort function for ingredient//
//The value returned indicates whether the element passed as first argument 
//is considered to go before the second in the specific strict weak ordering it defines.
//can weuse template here ? so ca accept any ingredient type..
bool ingredient_compare1(sphere x, sphere y){
    /*
    """
    sort ingredients using decreasing priority and decreasing radii for
    priority ties and decreasing completion for radii ties
    
    """
    */
    
    float p1 = x.packingPriority;
    float p2 = y.packingPriority;
    if (p1 < p2) //# p1 > p2
        return true;
    else if (p1==p2){ //# p1 == p1
       float r1 = x.minRadius;
       float r2 = y.minRadius;
       if (r1 > r2) //# r1 < r2
           return true;
       else if (r1==r2){ //# r1 == r2
           float c1 = x.completion;
           float c2 = y.completion;
           if (c1 >= c2) //# c1 > c2
               return true;
           else
               return false;
       }else
           return false;
    }else
       return false;
}

bool ingredient_compare0(sphere x, sphere y){
    /*
    """
    sort ingredients using decreasing priority and decreasing radii for
    priority ties and decreasing completion for radii ties
    
    """
    */
    
    float p1 = x.packingPriority;
    float p2 = y.packingPriority;
    if (p1 > p2) //# p1 > p2
        return true;
    else if (p1==p2){ //# p1 == p1
       float r1 = x.minRadius;
       float r2 = y.minRadius;
       if (r1 > r2) //# r1 < r2
           return true;
       else if (r1==r2){ //# r1 == r2
           float c1 = x.completion;
           float c2 = y.completion;
           if (c1 >= c2) //# c1 > c2
               return true;
           else
               return false;
       }else
           return false;
    }else
       return false;
}

bool ingredient_compare2(sphere x, sphere y){
    /*
    """
    sort ingredients using decreasing radii and decresing completion
    for radii matches
    
    """
    */
    
    float r1 = x.minRadius;
    float r2 = y.minRadius;
   if (r1 < r2) //# r1 < r2
       return true;
   else if (r1==r2){ //# r1 == r2
       float c1 = x.completion;
       float c2 = y.completion;
       if (c1 >= c2) //# c1 > c2
           return true;
       else
           return false;
   }
   else
    return false;
}



/*
the main struct that define the environement
should it be a class ?
*/
struct big_grid { // needs 8*n bytes 
    std::vector<unsigned> all; //all point indice
    std::vector<unsigned> empty;//available point indice
    std::vector<point> data;    //the grid 3d coordintates
    std::vector<sphere> ingredients;//the list of sphere ingredient to pack
    std::vector<sphere> activeIngr;
    std::vector<sphere> activeIngr0;
    std::vector<sphere> activeIngr12;
    std::vector<float> normalizedPriorities0;
    std::vector<float> normalizedPriorities;
    std::vector<float> thresholdPriorities; 
    std::vector<point> rtrans;    //the grid 3d coordintates
    
    //float* distance; 
    std::vector<float> distance;    //the array of closest distances for each point
    unsigned num_empty;         //the number of free point available
    unsigned num_points;        //total number of point in the grid
    float space;                //spacing of the grid (unit depends on user)
    float maxradius;            //the biggest ingredient
    float lowestPriority;       //priority the lowest after sorting
    float totalRadii;           //sum of all ingredient radii
    float totalPriorities;
    float vRangeStart;
    box nxyz;                   //the grid size in x, y and z
    point boundingBox0;         //the grid lower left corner coordinate
    point boundingBox1;         //the grid top right coordinate
    int mode;                   //the packing mode random or distance
    bool pickWeightedIngr = true;
    bool pickRandPt = true;
    bool use_gradient;
    int numActiveIngr;

    std::default_random_engine generator;
    std::uniform_real_distribution<float> uniform;// (0.0,1.0);
    std::normal_distribution<float> gauss;//(0.0,0.3);

    //the constructor that take as input the sizenor of the grid, the step, and the bouding box
    big_grid(unsigned nx, unsigned ny,unsigned nz, 
                float step, point bb0, point bb1,float seed) : 
        all(nx*ny*nz),
        empty(nx*ny*nz),
        num_empty(nx*ny*nz),
        data(nx*ny*nz),
        distance(nx*ny*nz),
        uniform(0.0,1.0),
        gauss(0.0,0.3){
        generator.seed (seed);
        num_points = nx*ny*nz;      //initialize total number of point
        boundingBox0 = bb0;
        boundingBox1 = bb1;
        mode = 1;
        nxyz.x = nx;//+encapsulatingGrid
        nxyz.y = ny;//+encapsulatingGrid
        nxyz.z = nz;//+encapsulatingGrid 
        space = step;
        float xl=bb0.x;
        float yl=bb0.y;
        float zl=bb0.z;
        unsigned i=0;           //counter for u indices             
        //distance = new float[nx*ny*nz];
        for(unsigned zi = 0; zi < nz; ++zi) { 
            for(unsigned yi = 0; yi < ny; ++yi) { 
                for(unsigned xi = 0; xi < nx; ++xi) { 
                    //build the coordinate and initialize the array of indice
                    point px;
                    px.x = xl+(float)xi*space;
                    px.y = yl+(float)yi*space;
                    px.z = zl+(float)zi*space;
                    data[i]=px;
                    all[i] = i;
                    empty[i] = i;
                    distance[i] = 999999.0;
                    //if(xi + 10 >= nx) // print the last 10 indexes
                    //    std::cout << i << ' ' << px.x <<' ' << px.y << ' ' << px.z <<' ' << distance[i] << '\n';
                    ++i;
                }
            }
        }
    }

    void setMode(unsigned _mode){
        //set the packing mode, random or distance
        mode = _mode;
    }

    void setIngredients(std::vector<sphere> _ingredients){
        //set the ingredients list to pack in the grid
        //retrieve the biggest one
        ingredients = _ingredients;
        maxradius = 0.0;
        for(unsigned i = 0; i < ingredients.size(); ++i) { 
            if (ingredients[i].radius > maxradius) 
                maxradius = ingredients[i].radius;
        }
        numActiveIngr = ingredients.size();
        activeIngr = ingredients;
    }
    
    void getMaxRadius(){
        maxradius = 0.0;
        for(unsigned i = 0; i < numActiveIngr; ++i) { 
            if (activeIngr[i].radius > maxradius) 
                maxradius = activeIngr[i].radius;
        }
    }

    void getSortedActiveIngredients(){
        //pirorities120 is ot used ??
        std::vector<sphere> ingr1;  // given priorities
        std::vector<float> priorities1;
        std::vector<sphere> ingr2;  // priority = 0 or none and will be assigned based on complexity
        std::vector<float> priorities2;
        std::vector<sphere> ingr0;  // negative values will pack first in order of abs[packingPriority]
        std::vector<float> priorities0;
        sphere lowestIng; 
        sphere ing; 
        float r=0.0;
        float np=0.0;      
        for(unsigned i = 0; i < ingredients.size(); ++i) { 
            ing = ingredients[i];
            if (ing.completion >= 1.0) continue;// # ignore completed ingredients
            if (ing.packingPriority == 0.0){
                ingr2.push_back(ingredients[i]);
                priorities2.push_back(ing.packingPriority);
            }
            else if (ing.packingPriority > 0.0 ){
                ingr1.push_back(ingredients[i]);
                priorities1.push_back(ing.packingPriority);
            }else{
                ingr0.push_back(ingredients[i]);
                priorities0.push_back(ing.packingPriority);
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
            lowestPriority = lowestIng.packingPriority;
        }else 
            lowestPriority = 1.0;

        totalRadii = 0.0;
        for(unsigned i = 0; i < ingr2.size(); ++i) {
            ing = ingr2[i];
            r = ing.minRadius;
            totalRadii = totalRadii + r;
            if (r==0.0) { 
                //safety 
                totalRadii = totalRadii + 1.0;
             }
        }
        for(unsigned i = 0; i < ingr2.size(); ++i) {
            ing = ingr2[i];
            r = ing.minRadius;
            np = float(r)/float(totalRadii) * lowestPriority;
            //std::cout << "#packingPriority " << np << ' ' << r <<' '<< totalRadii<< ' ' << lowestPriority << '\n';
            normalizedPriorities0.push_back(np);
            ingr2[i].packingPriority = np;
        }           
        
        activeIngr0 = ingr0;//#+ingr1+ingr2  #cropped to 0 on 7/20/10
        activeIngr12 = ingr1;//+ingr2
        activeIngr12.insert(activeIngr12.end(), ingr2.begin(), ingr2.end());
        //packingPriorities = priorities0;//+priorities1+priorities2
    }

    void prepareIngredient(){
        float pp;
        float np;
        sphere ingr;
        getSortedActiveIngredients();
        std::cout << "#len(allIngredients) " << ingredients.size() << '\n';
        std::cout << "#len(activeIngr0) " << activeIngr0.size() << '\n';
        std::cout << "#len(activeIngr12) " << activeIngr12.size() << '\n';
        //self.activeIngre_saved = self.activeIngr[:]
        totalPriorities = 0.0;// # 0.00001
        for(unsigned i = 0; i < activeIngr12.size(); ++i) {
            ingr = activeIngr12[i];
        //for priors in self.activeIngr12:
            pp = ingr.packingPriority;
            totalPriorities = totalPriorities + pp;
            std::cout << "#totalPriorities " << totalPriorities << ' ' << ingr.packingPriority <<'\n';
        }
        
        
        float previousThresh = 0.0;
        //# Graham- Once negatives are used, if picked random# 
        //# is below a number in this list, that item becomes 
        //# the active ingredient in the while loop below
        for(unsigned i = 0; i < activeIngr0.size(); ++i) {
        //for priors in self.activeIngr0:
            normalizedPriorities.push_back(0.0);
            if (pickWeightedIngr)// #why ?
                thresholdPriorities.push_back(2);
        }        
        for(unsigned i = 0; i < activeIngr12.size(); ++i) {
        //for priors in self.activeIngr12:
            //#pp1 = 0
            pp = activeIngr12[i].packingPriority;
            if (totalPriorities != 0)
                np = float(pp)/float(totalPriorities);
            else
                np=0.0;
            normalizedPriorities.push_back(np);
            std::cout << "#np is "<< np << " pp is "<< pp << " tp is " << np + previousThresh << '\n';
            thresholdPriorities.push_back(np + previousThresh);
            previousThresh = np + float(previousThresh);
        }
        activeIngr = activeIngr0;// + self.activeIngr12
        activeIngr.insert(activeIngr.end(), activeIngr12.begin(), activeIngr12.end());

        int nls=0;
        int totalNumMols = 0;
        if (thresholdPriorities.size() == 0){
            for(unsigned i = 0; i < ingredients.size(); ++i) {
                totalNumMols += ingredients[i].nbMol;
            }
            std::cout << "#totalNumMols Fill5if = " << totalNumMols << '\n';
        }else {                
            for(unsigned i = 0; i < thresholdPriorities.size(); ++i) {
                float threshProb = thresholdPriorities[i];
                sphere ing = activeIngr[nls];
                std::cout << "#threshprop Fill5else is "<< threshProb << " for ingredient : "<< ing.name<< ' '  << ing.nbMol<< '\n' ;
                totalNumMols += ing.nbMol;
                std::cout << "#totalNumMols Fill5else = " << totalNumMols << '\n';
                nls++;
            }
        }    
    }

    void dropIngredient(sphere ingr){
        int ingr_ind;
        bool found = false;
        for(unsigned i = 0; i < numActiveIngr; ++i) { 
            if (ingr.name == activeIngr[i].name){ 
                ingr_ind = i;
                found = true;
                break;             
            }
        }
        if (found) {
            //swap  
            numActiveIngr--;
            std::swap(activeIngr[ingr_ind], activeIngr[numActiveIngr]);  
            //swap also the ingredient ? active Ingredient should be unsign int ? 
            ingr.active = false;   
        }
        getMaxRadius();
    }

    sphere pickIngredient(){
        sphere ingr;
        //float prob;
        //std::default_random_engine generator (seed);
        //std::default_random_engine generator(seed);
        //std::uniform_real_distribution<float> distribution(0.0,1.0);
        //std::normal_distribution<float> distribution(0.0,1.0);        
        int ingrInd;
        float threshProb;
        float r;
        int n ;
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
                //std::cout << "pick Ingredient "  << ingrInd << ' ' << prob <<' ' <<threshProb << '\n';
                    
                if (ingrInd <  numActiveIngr)
                    ingr = activeIngr[ingrInd];
                else {
                    std::cout << "#error in histoVol pick Ingredient "  << ingrInd << ' ' <<numActiveIngr << '\n';
                    ingr = activeIngr[0];
                }
            }
        }else {
            //pick random
            ingr = activeIngr[rand() % numActiveIngr];
        }
        return ingr;
    }

    inline sphere sample_ingredient() const {
        //randomly pick a ingredient from the list of ingredients
        return ingredients[rand() % numActiveIngr]; 
    }

    inline unsigned sample_empty() const {
        //randomly pick a point
        return empty[rand() % num_empty]; 
    }
    
    inline unsigned sample_empty_distance(float radius) const {
        //randomly pick a point that can accept an ingredient of a given radius
        float cut = radius - 1.0;
        std::vector<unsigned> freept;
        for(unsigned i = 0; i < num_empty; ++i) { 
            unsigned pt = empty[i];
            float d=distance[pt];
            if (d>=cut) freept.push_back(pt);
        }       
        return freept[rand() % freept.size()]; 
    }

    inline unsigned sample_closest_empty(std::vector<float> Dist,std::vector<unsigned> PointID) const {
        //we want the indice of the point that the smallest distance
        float cut  = 10.0;      
        float d;  
        float mini = 9999999999.9;
        unsigned ptId=0;
        for(unsigned i = 0; i < Dist.size(); ++i) { 
            d=Dist[empty[i]];
            if (d >= cut){
                if (d < mini) {
                    mini = d;
                    ptId = PointID[i];
                }
            }                  
        }
        return ptId;
    }

    int getPointToDrop(sphere ingr, float radius,float jitter){
        //#should we take in account a layer? cuttof_boundary and cutoff_surface?
        unsigned ptInd;
        float cut= radius-jitter;
        float d;
        std::vector<unsigned> allIngrPts;
        std::vector<float> allIngrDist;
        if (ingr.packingMode=="close"){
            for(unsigned i = 0; i < num_empty; ++i) {
            //for pt in freePoints:#[:nbFreePoints]:
                d = distance[empty[i]];//#look up the distance
                if (d>=cut)
                    allIngrPts.push_back(empty[i]);
                    allIngrDist.push_back(d);
            }
        }
        else{
            if (mode == 1) {//distance
                for(unsigned i = 0; i < num_empty; ++i) { 
                    d=distance[empty[i]];
                    if (d>=cut) allIngrPts.push_back(empty[i]);
                }       

            }
            else if (mode == 0) {//pure random
                allIngrPts = empty;
            }
        }
        //std::cout << "allIngrPts " <<allIngrPts.size()<<'\n';
        if (allIngrPts.size()==0){
            sphere tmp[] = {ingr}; 
            ingr.completion = 1.0;
            vRangeStart = vRangeStart + normalizedPriorities[0];
            dropIngredient(ingr); 
            getSortedActiveIngredients();
            //drop the ingredient and update
            //update the priorities ?
            //ind = activeIngr.index(ingr)
            /*std::vector<int>::iterator it;
            it = std::search(activeIngr.begin(), activeIngr.end(), tmp, tmp+1);
            int ind = (it-activeIngr.begin());
            vRangeStart = vRangeStart + normalizedPriorities[0];
            if (ind > 0){
                for(unsigned j = 0; j < ind; ++j) { 
                //for j in range(ind):                
                    thresholdPriorities[j] = thresholdPriorities[j] + normalizedPriorities[ind];
                }
            }*/
            //activeIngr.erase(it);
            //# Start of massive overruling section from corrected thesis file of Sept. 25, 2012
            //#this function also depend on the ingr.completiion that can be restored ?
            //getSortedActiveIngredients();//, (self.activeIngr,False))
            totalPriorities = 0; //# 0.00001
            //thresholdPriorities.pop(it)
            //normalizedPriorities.pop(ind)
            //if verbose:
            //        print(("time to reject the picking", time()-t))
            //print(("vRangeStart",vRangeStart))
            return -1;                   
            //seems duplicate from getSortedActiveIngredient
/*            for(unsigned i = 0; i < activeIngr12.size(); ++i) {
            //for priors in self.activeIngr12:
                float pp = activeIngr12[i].packingPriority;
                totalPriorities = totalPriorities + pp;
                //std::cout << "totalPriorities = "<< totalPriorities<< '\n';
            }
            previousThresh = 0
            self.normalizedPriorities = []
            self.thresholdPriorities = [] 
            # Graham- Once negatives are used, if picked random# 
            # is below a number in this list, that item becomes 
            #the active ingredient in the while loop below
            for priors in self.activeIngr0:
                self.normalizedPriorities.append(0)
                if self.pickWeightedIngr :
                    self.thresholdPriorities.append(2)
            for priors in self.activeIngr12:
                #pp1 = 0
                pp = priors.packingPriority
                if self.totalPriorities != 0:
                    np = float(pp)/float(self.totalPriorities)
                else:
                    np=0.
                self.normalizedPriorities.append(np)
                if verbose :
                    print ('np is ', np, ' pp is ', pp, ' tp is ', np + previousThresh)
                self.thresholdPriorities.append(np + previousThresh)
                previousThresh = np + float(previousThresh)
            self.activeIngr = self.activeIngr0 + self.activeIngr12
            
*/
        }
        if (pickRandPt){
            //std::cout << "allIngrPts " <<allIngrPts.size()<<'\n';
            if (ingr.packingMode=="close")
                ptInd = sample_closest_empty(allIngrDist,allIngrPts);
            else if (ingr.packingMode=="gradient") //&& (use_gradient)  
                ptInd =0;// self.gradients[ingr.gradient].pickPoint(allIngrPts) 
            else{
                ptInd = allIngrPts[rand() % allIngrPts.size()];
                if (ptInd > allIngrPts.size()) ptInd = allIngrPts[0];            
            }     
            /*if ptInd is None :
                t=time()
                if verbose:
                    print('No point left for ingredient %s %f minRad %.2f jitter %.3f in component %d'%(
                    ingr.name, ingr.molarity, radius, jitter, compNum))
                ingr.completion = 1.0
                ind = self.activeIngr.index(ingr)
                #if ind == 0:
                vRangeStart = vRangeStart + self.normalizedPriorities[0]
                if ind > 0:
                    #j = 0
                    for j in range(ind):                
                        self.thresholdPriorities[j] = self.thresholdPriorities[j] + self.normalizedPriorities[ind]
                self.activeIngr.pop(ind)
                if verbose:
                    print('popping this gradient ingredient array must be redone using Sept 25, 2011 thesis version as above for nongraient ingredients, TODO: July 5, 2012')
                self.thresholdPriorities.pop(ind)
                self.normalizedPriorities.pop(ind)
                if verbose:
                        print(("time to reject the picking", time()-t))
                print(("vRangeStart",vRangeStart))
                return False,vRangeStart                    

#            print(("time to random pick a point", time()-t2))
*/
        }else {
            std:sort(allIngrPts.begin(),allIngrPts.end());//-(allIngrPts.size()-numActiveIngr)
            ptInd = allIngrPts[0];
        }
        return ptInd;
    }

    inline float ManhattanDistance( point c1, point c2 ) const
    {
        //distance approximation
        float dx = abs(c2.x - c1.x);
        float dy = abs(c2.y - c1.y);
        float dz = abs(c2.z - c1.z);   
        return dx+dy+dz;
    }

    inline float Distance( point c1, point c2 ) const
    {
        //real distance between two point in 3d
        float s=0.0;
        point d;
        d.x = (c2.x - c1.x);
        d.y = (c2.y - c1.y);
        d.z = (c2.z - c1.z);   
        s = d.x*d.x+d.y*d.y+d.z*d.z;
        return sqrt(s);
    }

    bool test_data(unsigned pid,sphere ingr)  {
        //main function that decide to drop an ingredient or not
        //std::cout  <<"test_data "<< ingr.name << ' ' << pid << '\n';        
        float rad = ingr.radius;        
        point px = data[pid];   //the selected point where we want to drop
        point target;           //the point with some jitter
        float d = 0.0;          //distance to be computed
        unsigned nbJitter = ingr.nbJitter;  //nb of jitter
        //should be defined in ingredient
        float jx=ingr.jitterMax.x;           //jitter amount on x 
        float jy=ingr.jitterMax.y;            //jitter amount on y 
        float jz=ingr.jitterMax.z;            //jitter amount on z 
        //actuall jitter that will be apply to the point
        float dx=0.0;
        float dy=0.0;
        float dz=0.0;
        //square jitter
        float d2;
        float jitter = space;
        float jitter2 = jitter * jitter;
        bool collision;
        //prepare the normal distribution for generating the jitter offset
        //std::default_random_engine generator;
        //std::normal_distribution<float> distribution(0.0,0.3);
        //std::cout  << "testData" << '\n';
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
            ++totnbJitter;
            target.x = px.x + dx;// = (tx+dx, ty+dy, tz+dz)
            target.y = px.y + dy;//
            target.z = px.z + dz;//
            //jitterLength += dx*dx + dy*dy + dz*dz  //#Why is this expensive line needed?
            //jitterList.append( (dx,dy,dz) )      
            //check for collision at the given target point coordinate for the given radius      
            collision = checkSphCollisions(target,rad);
            if (!collision) {
                ingr.trans.x = target.x;
                ingr.trans.y = target.y;
                ingr.trans.z = target.z;
//                std::cout << pid << " test data " << target.x <<' ' << target.y << ' ' << target.z <<' ' << collision << '\n'; 
//                std::cout << pid << " test data " << ingr.trans.x <<' ' << ingr.trans.y << ' ' << ingr.trans.z <<' ' << collision << '\n'; 
                rtrans.push_back(target);                
                break;
            }
        }   
        return collision;
    }

    void set_filled(unsigned i) {
        //marke the given indice as occupide and remove it by swapping
        //Oleg contribution
        unsigned empty_index = all[i];
        if(empty_index < num_empty) { // really is empty
            --num_empty;
            unsigned filled_index = empty[num_empty]; 
            std::swap(empty[empty_index], empty[num_empty]); 
            std::swap(all[i], all[filled_index]);
        }
    }

    inline bool is_empty(unsigned i) const { 
        //check the statut of this point
        return all[i] < num_empty;
    }


    std::vector<unsigned> getPointsInCube(point bb0,point bb1,point pos,float radii) const {
        //for a given boundin box retrieve the indice of the grid point that are inside
        std::vector<unsigned> ptIndices;
        
        const float spacing1 = 1.0 / space;
        const int NX = nxyz.x; 
        const int NY = nxyz.y; 
        const int NZ = nxyz.z; 
        const float OX = boundingBox0.x;
        const float OY = boundingBox0.y;
        const float OZ = boundingBox0.z; //# origin of fill grid-> bottom lef corner not origin
        //std::cout << "bb0 " << bb0.x << ' ' << bb0.y << ' ' << bb0.z<< '\n';
        //std::cout << "bb1 " << bb1.x << ' ' << bb1.y << ' ' << bb1.z<< '\n';
        
        int i0 = std::max(0, (int)floor((bb0.x-OX)*spacing1));
        int i1 = std::min(NX, (int)((bb1.x-OX)*spacing1)+1);

        int j0 = std::max(0, (int)floor((bb0.y-OY)*spacing1));
        int j1 = std::min(NY, (int)((bb1.y-OY)*spacing1)+1);

        int k0 = std::max(0, (int)floor((bb0.z-OZ)*spacing1));
        int k1 = std::min(NZ, (int)((bb1.z-OZ)*spacing1)+1);

        //std::cout << "i0i1... " << i0 << ' ' << i1 << ' ' << j0 << ' ' << j1 << ' ' << k0 << ' ' << k1 << '\n';
        int zPlaneLength = NX*NY;
        int offz =0;
        int off = 0;
        for(int z = k0; z < k1; ++z) { 
            offz = z*zPlaneLength ;          
            for(int y = j0; y < j1; ++y) { 
                off = y*NX + offz;
                for(int x = i0; x < i1; ++x) { 
                    ptIndices.push_back( x + off);
                }
            }
        }
        return ptIndices;
    }

    void updateDistance(point pos, float radii){
        //update the distance array
        //std::cout << "updateDistance" << '\n';
        float d=0.0;
        //do it only in the bounding box covering maxradius
        float x = pos.x;
        float y = pos.y;
        float z = pos.z;
        point bb0 ;
        bb0.x = x-radii-(maxradius*2.0);
        bb0.y = y-radii-(maxradius*2.0);
        bb0.z = z-radii-(maxradius*2.0);
        point bb1 ;
        bb1.x = x+radii+(maxradius*2.0);
        bb1.y = y+radii+(maxradius*2.0);
        bb1.z = z+radii+(maxradius*2.0);
        std::vector<unsigned> pointsInCube = getPointsInCube(bb0,bb1, pos, radii);
        for(unsigned i = 0; i < pointsInCube.size(); ++i) { 
            d = Distance(data[pointsInCube[i]],pos); 
            if (d < distance[i]) distance[i] = d; 
        }       
    }

    bool checkSphCollisions(point pos, float radii)  {
            //check sphere collision
            //look first in the sphere bounding box
            //then the point indice inside the sphere
            //if point are free continue
            //if found a point that is occupied -> colision
            bool collision = false;            
            float x = pos.x;
            float y = pos.y;
            float z = pos.z;
            point bb0 ;
            bb0.x = x-radii;
            bb0.y = y-radii;
            bb0.z = z-radii;
            point bb1 ;
            bb1.x = x+radii;
            bb1.y = y+radii;
            bb1.z = z+radii;
            std::vector<unsigned> pointsInCube = getPointsInCube(bb0,bb1, pos, radii);//#indices
            std::vector<unsigned> ptsInSphereId;
            std::vector<float> dist;
            //std::cout << " ptinCube " << pointsInCube.size() << '\n';
            //mesure distance grid points coordinate of indice pointsInCube from pos
            //point in sphere are point where distance less than radii
            float d=0.0;
            for(unsigned i = 0; i < pointsInCube.size(); ++i) { 
                d = Distance(data[pointsInCube[i]],pos);
                if (d < radii){
                    ptsInSphereId.push_back(pointsInCube[i]);
                    dist.push_back(d-radii);
                }
                //std::cout << i << ' ' <<  " distance " << d <<" radiu " << radii << '\n';
                //std::cout << pointsInCube[i] << " coord " << data[pointsInCube[i]].x <<' ' << data[pointsInCube[i]].y << ' ' << data[pointsInCube[i]].z << '\n';  
            }
            //std::cout << " ptinSphere " << ptsInSphereId.size() << '\n';
            for(unsigned i = 0; i < ptsInSphereId.size(); ++i) { 
                //std::cout << i << " distance " << distance[ptsInSphereId[i]] << " empty " << is_empty(ptsInSphereId[i]) << '\n';
                //d = dist-raii;
                if (mode==1){
                    if (distance[ptsInSphereId[i]]<-0.0001) {   
                        collision = true;  
                        break; 
                    }
                }
                if (!is_empty(ptsInSphereId[i])){
                    collision = true;  
                    break;                 
                }
                //if (d+distance[pt] < histoVol.maxColl:)
                //    histoVol.maxColl = d+distance[pt]
                //        #print("in collision histovol.maxColl if")
                //return True
            }
            //std::cout << "colliding ? " << collision  << '\n';
            if (!collision) {
                for(unsigned i = 0; i < ptsInSphereId.size(); ++i) { 
                    //d = dist[i]-radii;    
                    //if (dist[i] < radii)  //# point is inside dropped sphere
                    set_filled(ptsInSphereId[i]);
                    //if (dist[i] < distance[ptsInSphereId[i]]) //# point in region of influence
                    //        distance[ptsInSphereId[i]] = dist[i];//??
                    //std::cout << "update " << ptsInSphereId[i] << '\n';
                    //std::cout << i << " distance " << distance[ptsInSphereId[i]] << " di " << dist[i] << '\n';
                
                            //std::swap(distance[ptsInSphereId[i]], dist[i]);
                            //std::swap_ranges(distance.begin()+ptsInSphereId[i], distance.begin()+ptsInSphereId[i], dist.begin()+i);
                } 
                //distance for all point?
                if (mode==1) updateDistance(pos,radii);           
            }
            return collision; 
        }
};

inline box computeGridNumberOfPoint(point bb0,point bb1, float space){
    //compute number of points  given the bounding box and the stepsize
    box nxyz;   
    nxyz.x = unsigned(ceil((bb1.x-bb0.x)/space));//+encapsulatingGrid
    nxyz.y = unsigned(ceil((bb1.y-bb0.y)/space));//+encapsulatingGrid
    nxyz.z = unsigned(ceil((bb1.z-bb0.z)/space));//+encapsulatingGrid
    return nxyz;
}

template<typename T> void test(point bb0,point bb1, float smallestObjectSize,
            std::vector<sphere> _ingredients, int mode, float seed,
            bool pickWeightedIngr, bool pickRandPt) { 
    //Main Pack function that will go until no more free point
    //or broke due to too much rejection.
    //should ad stop by completion
    float space = smallestObjectSize*1.1547;
    box nxyz = computeGridNumberOfPoint(bb0,bb1,space); 
    std::cout << "#grid step " << space << " grid Npoint " << nxyz.x*nxyz.y*nxyz.z << '\n';
    std::cout << "#"<< nxyz.x<< ' ' <<nxyz.y<< ' ' <<nxyz.z << '\n';
    
    clock_t start;
    double diff;
    start = clock();
    T g(nxyz.x,nxyz.y,nxyz.z,space,bb0,bb1,seed);
    diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
    std::cout<<"#time: "<< diff <<'\n';
    //mainLoop for packing
    std::vector<point> usedPoints;
    std::vector<float> radiis;
    std::vector<point> colors;
    int rejection=0;
    g.vRangeStart=0.0;
    int PlacedMols = 0;
    g.pickWeightedIngr=pickWeightedIngr;
    g.pickRandPt =pickRandPt; 
    g.setMode(mode);//1-close packing 0-random
    g.setIngredients(_ingredients);
    g.prepareIngredient();
    //g.setMode(0);//random packing
    while(g.num_empty > 1) {     
    //for(unsigned i = 0; i < 2 ; ++i) {//nx*ny*nz
        //pick the object to drop
        if (g.vRangeStart>1.0){
                std::cout << "#broken by vRange and hence Done!!!****\n";
                break;
        }
        sphere ingr = g.pickIngredient();
        //sphere ingr= g.sample_ingredient();//sampl using distance information as well ?
        if (!ingr.active) continue;        
        if (ingr.completion >= 1.0) {
            g.dropIngredient(ingr);            
            continue;
        }
        //std::cout  << ingr.name << '\n';
        unsigned s = g.sample_empty(); //get the next available point randomly
        //int s = g.getPointToDrop(ingr, ingr.radius,1.0);
        if ( s == -1 ){
            continue;
        }
        //unsigned s = g.sample_empty(); 
        //std::cout  << s << '\n';
        //if (mode == 1) s=g.sample_empty_distance(ingr.radius); //get the next available point randomly  
        //else if (mode == 0)  s=g.sample_empty();    
        //unsigned s = g.sample_closest_empty();        
        if (!g.test_data(s,ingr)){
            //drop ingredient
            //usedPoints.push_back(ingr.trans);//oups the jitter!!
            //std::cout  << ingr.name << ' ' << g.rtrans.back().x <<' ' << g.rtrans.back().y << ' ' << g.rtrans.back().z << '\n'; 
            radiis.push_back(ingr.radius);
            colors.push_back(ingr.color);
            rejection=0;  
            ingr.counter++;
            ingr.completion = float(ingr.counter)/float(ingr.nbMol);
        }
        else {
            ingr.rejectionCounter++;
            if (ingr.rejectionCounter > ingr.rejectionThreshold)
                ingr.completion =1.0;//erase from list ingredient /
            rejection++;            
        }      
        //g.set_filled(s); 
        //std::cout << "# main loop " << ingr.name << ' ' << s << ' ' << g.num_empty << ' ' << rejection << '\n';
        if   (rejection > 3000) break;
        //when stop
        //fill the point if distance to something correct? 
        //if( i % 100 == 0) // print the last 10 indexes
        //    std::cout << i << ' ' << s << '\n';
    }
    //std::cout << "#main loop " << '\n';
    
    //fancy python ouput for easy visualization in 3d Host
    std::cout << "#main loop " << g.rtrans.size() << '\n';
    std::cout << "pts=[" << '\n';
    for(unsigned i = 0; i < g.rtrans.size(); ++i) { 
        //std::cout << i << ' ' << usedPoints[i] << '\n';
        std::cout << '[' << g.rtrans[i].x <<',' << g.rtrans[i].y << ',' << g.rtrans[i].z << ']' << ',' <<'\n';  
    } 
    std::cout << "]" << '\n'; 
  
    std::cout << "r=[" << '\n'  ;  
    for(unsigned i = 0; i < radiis.size(); ++i) { 
        //std::cout << i << ' ' << usedPoints[i] << '\n';
        std::cout << radiis[i] << ',' <<'\n';  
    }    
    std::cout << "]" << '\n';

    std::cout << "color=[" << '\n'  ;  
    for(unsigned i = 0; i < colors.size(); ++i) { 
        //std::cout << i << ' ' << usedPoints[i] << '\n';
        std::cout << '[' << colors[i].x <<',' << colors[i].y << ',' << colors[i].z << ']' << ',' <<'\n';  
    }    
    std::cout << "]" << '\n';
    std::cout << "import upy\n" << "helper = upy.getHelperClass()()\n";
    std::cout << "pesph=helper.newEmpty(\"base_sphere\")\n";
    std::cout << "bsph=helper.Sphere(\"sphere\",parent=pesph)[0]\n";
    std::cout << "isph=helper.instancesSphere(\"cpp\",pts,r,pesph,color,None)\n";
//#execfile("/Users/ludo/DEV/testCpp.py")
//#5sec step 5 mode 1
   
}

//helper to create a singleSphere ingredient given a radius
inline sphere makeSphere(float radius, int mode, float concentration, 
         float packingPriority,int nbMol,std::string name, point color,
        unsigned nbJitter,point jitterMax){
    sphere sp;
    sp.radius = radius;
    sp.mode = mode;
    sp.minRadius = radius;
    sp.nbMol=nbMol;
    sp.completion = 0.0;
    sp.packingPriority = packingPriority;
    sp.name = name;
    sp.counter = 0;
    sp.rejectionCounter = 0;
    sp.rejectionThreshold  = 30;
    sp.active = true;
    sp.color = color;
    sp.nbJitter = nbJitter;
    sp.trans.x = 0.0;
    sp.trans.y = 0.0;
    sp.trans.z = 0.0;
    sp.jitterMax = jitterMax;
    return sp;
}

inline point makePoint(float a,float b, float c ){
    point p;
    p.x = a;
    p.y = b;
    p.z = c;
    return p;
}

float getRadii(std::string str){
    //std::string str(input);  
    std::replace(str.begin(), str.end(), '[', ' ');  
    std::replace(str.begin(), str.end(), ']', ' ');  
    std::replace(str.begin(), str.end(), '(', ' ');  
    std::replace(str.begin(), str.end(), ')', ' ');  
    return atof(str.c_str());
}

std::vector<float> getBox(std::string str){
    std::replace(str.begin(), str.end(), '[', ' ');  
    std::replace(str.begin(), str.end(), ']', ' ');  
    std::replace(str.begin(), str.end(), '(', ' ');  
    std::replace(str.begin(), str.end(), ')', ' ');  
    std::replace(str.begin(), str.end(), ',', ' ');
    // Build an istream that holds the input string
    std::istringstream iss(str);
    //std::cout << str << '\n';
    // If possible, always prefer std::vector to naked array
    std::vector<float> v;
    // Iterate over the istream, using >> to grab floats
    // and push_back to store them in the vector
    std::copy(std::istream_iterator<float>(iss),
        std::istream_iterator<float>(),
        std::back_inserter(v));
    return v;
}

point getArray(std::string str){
    std::replace(str.begin(), str.end(), '[', ' ');  
    std::replace(str.begin(), str.end(), ']', ' ');  
    std::replace(str.begin(), str.end(), '(', ' ');  
    std::replace(str.begin(), str.end(), ')', ' ');  
    std::replace(str.begin(), str.end(), ',', ' ');
    // Build an istream that holds the input string
    std::istringstream iss(str);
    //std::cout << str << '\n';
    // If possible, always prefer std::vector to naked array
    std::vector<float> v;
    // Iterate over the istream, using >> to grab floats
    // and push_back to store them in the vector
    std::copy(std::istream_iterator<float>(iss),
        std::istream_iterator<float>(),
        std::back_inserter(v));
    point p;
    p.x = v[0];
    p.y = v[1];
    p.z = v[2];
    return p;
}
void load_xml(std::string path,int _mode,float _seed){
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(path.c_str());
    //std::cout << path << '\n';
    std::cout << "#Load result: " << result.description() << ", AutoFillSetup name: " << doc.child("AutoFillSetup").attribute("name").value() << std::endl;
    //get the option
    //std::cout <<doc.child("AutoFillSetup").child("options").attribute("smallestProteinSize").value()<<'\n';
    const float smallestObjectSize = atof(doc.child("AutoFillSetup").child("options").attribute("smallestProteinSize").value());
    const float seed = _seed;
    const int mode = _mode;
    //std::cout <<smallestObjectSize<<'\n';
    //std::cout <<seed<<'\n';
    std::string strbb(doc.child("AutoFillSetup").child("options").attribute("boundingBox").value());
    std::vector<float> bb = getBox(strbb);
    point bb0 ;
    bb0.x = bb[0];
    bb0.y =  bb[1];
    bb0.z =  bb[2];
    point bb1 ;
    bb1.x =  bb[3];
    bb1.y =  bb[4];
    bb1.z =  bb[5];

    bool pickWeightedIngr =doc.child("AutoFillSetup").child("options").attribute("pickWeightedIngr").as_bool();
    bool pickRandPt =doc.child("AutoFillSetup").child("options").attribute("pickRandPt").as_bool();

    //for (unsigned i=0;i < bb.size(); i++){
    //    std::cout << '#' << bb[i] << '\n';
    //}
    //get the ingredient.
    //only cytoplsme for now
    std::vector<sphere> _ingredients;
    pugi::xml_node cytoplasme = doc.child("AutoFillSetup").child("cytoplasme");
    //float radius, int mode, float concentration, 
    //     float packingPriority,int nbMol,std::string name, point color    
    for (pugi::xml_node ingredient = cytoplasme.child("ingredient"); ingredient; ingredient = ingredient.next_sibling("ingredient")){
        std::cout << "#ingredient " << ingredient.attribute("name").value() << '\n';
        std::string str(ingredient.attribute("radii").value());         
        float r = getRadii(str);   
        //std::cout << r << '\n';
        float mol = ingredient.attribute("molarity").as_float();
        //std::cout << mol << '\n';
        float priority = ingredient.attribute("packingPriority").as_float();
        //std::cout << priority << '\n';        
        int nMol = ingredient.attribute("nbMol").as_int();
        //std::cout << nMol << '\n'; 
        std::string iname(ingredient.attribute("name").value()); 
        //std::cout << iname << '\n'; 
        std::string strcol(ingredient.attribute("color").value());
        point color = getArray(strcol);
        unsigned nbJitter = ingredient.attribute("nbJitter").as_int();
        //std::cout << color.x << ' '<< color.y << ' ' << color.z << '\n';
        std::string strjitter(ingredient.attribute("jitterMax").value());
        point jitter =  getArray(strjitter);       
        sphere ingr = makeSphere(r,_mode,mol,priority,nMol,iname,
                                color,nbJitter,jitter);
        _ingredients.push_back(ingr);
    }
    test<big_grid>(bb0,bb1,smallestObjectSize,_ingredients,_mode,_seed,pickWeightedIngr,pickRandPt); // runs in 0.22 s test<big_grid>(n); // runs in 0.06 s
}

//1000,1000,10 / 15
int main(int argc, char* argv[]) {
    //for a grid we need a boundin box
    //[[-500.,-500,-0.5],[500,500,0.5]]
  

    point bb0 ;
    bb0.x = -500.0;
    bb0.y = -500.0;
    bb0.z = -0.5;
    point bb1 ;
    bb1.x = 500.0;
    bb1.y = 500.0;
    bb1.z = 0.5;
    //const float smallestObjectSize = atoi(argv[1]);
    int mode = atoi(argv[2]);
    float seed = atof(argv[3]);

    std::string filename = argv[1];    
    load_xml(filename,mode,seed);
    /*std::vector<sphere> _ingredients;
    _ingredients.resize(5);
    _ingredients[0] = makeSphere(200.0,0,0.001,0,4,"5_n200",makePoint(0.827, 0.827, 0.827));
    _ingredients[1] = makeSphere(100.0,0,0.01,0,8,"10_n100",makePoint(0.498, 0.498, 0.498));
    _ingredients[2] = makeSphere(50.0,0,0.1,0,16,"200_n50",makePoint(0.306, 0.451, 0.816));
    _ingredients[3] = makeSphere(25.0,0,1.0,0,32,"x400_n25",makePoint(0.784, 0.204, 0.204));
    _ingredients[4] = makeSphere(25.0,0,10.0,0,32,"y400_n25",makePoint(0.467, 0.239, 0.972));
    */
    //and a stepsize 
    //const unsigned nx = 58;//58; 
    //const unsigned ny = 58;//58; 
    //const unsigned nz = 1;//1 
    //should use xml setup from autoPack
    //and save result as apr or python
    //test<big_grid>(bb0,bb1,smallestObjectSize,_ingredients,mode,seed); // runs in 0.22 s test<big_grid>(n); // runs in 0.06 s
}

//still problem ingredient picking
//pointdropping picking
//os.system("/Users/ludo/DEV/autopack /Users/ludo/DEV/autofill_svn/trunk/AutoFillClean/autoFillRecipeScripts/2DsphereFill/Spheres2D.xml 0 10 > /Users/ludo/DEV/testCpp.py");execfile("/Users/ludo/DEV/testCpp.py")
