
#include <numeric>
#include "IngredientsDispatcher.h"

namespace {
/* the comparison function are strict translation from the python code */
//sort function for ingredient//
//The value returned indicates whether the element passed as first argument 
//is considered to go before the second in the specific strict weak ordering it defines.
//can weuse template here ? so ca accept any ingredient type..
bool ingredient_compare1(Ingredient* x, Ingredient* y){
    /*
    """
    sort ingredients using decreasing priority and decreasing radii for
    priority ties and decreasing completion for radii ties
    used for initial priority > 0
    """
    */
    
    double p1 = x->packingPriority;
    double p2 = y->packingPriority;
    if (p1 < p2) //# p1 > p2
        return true;
    else if (p1==p2){ //# p1 == p1
       double r1 = x->minRadius;
       double r2 = y->minRadius;
       if (r1 > r2) //# r1 < r2
           return true;
       else if (r1==r2){ //# r1 == r2
           double c1 = x->completion;
           double c2 = y->completion;
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
    used for initial priority < 0
    """
    */
    
    double p1 = x->packingPriority;
    double p2 = y->packingPriority;
    if (p1 > p2) //# p1 > p2
        return true;
    else if (p1==p2){ //# p1 == p1
       double r1 = x->minRadius;
       double r2 = y->minRadius;
       if (r1 > r2) //# r1 < r2
           return true;
       else if (r1==r2){ //# r1 == r2
           double c1 = x->completion;
           double c2 = y->completion;
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
    used for initial priority == 0
    """
    */
    
    double r1 = x->minRadius;
    double r2 = y->minRadius;
    //ask graham about this issue
   if (r1 < r2) //# is r1 < r2 - this should be r1 > r2
       return true;
   else if (r1==r2){ //# r1 == r2
       double c1 = x->completion;
       double c2 = y->completion;
       if (c1 > c2) //# c1 > c2
           return true;
       else
           return false;
   }
   else
    return false;
}

}

IngradientsDispatcher::IngradientsDispatcher( std::vector<Ingredient> const & _ingredients, openvdb::Index64 num_points, unsigned int seed) 
    :uniform(0.0, 1.0)
    , vRangeStart(0.0)
    , generator(std::default_random_engine(seed))
    , pickWeightedIngr(true)
{
    //set the ingredients list to pack in the grid
    //retrieve the biggest one
    ingredients = _ingredients;
    activeIngr.resize(ingredients.size());
    double unit_volume = pow(stepsize,3.0);
    double grid_volume =  num_points*unit_volume;
    std::cout << "#Grid Volume " << grid_volume << " unitVol " << unit_volume << std::endl;
    for(unsigned i = 0; i < ingredients.size(); ++i) { 
        activeIngr[i] = &ingredients[i];
        ingredients[i].setCount(grid_volume);
    }    

    prepareIngredient();

    int totalNumberOfMols = std::accumulate(ingredients.cbegin(), ingredients.cend(), 0, [](int acc, Ingredient const& ingr) { return acc+=ingr.nbMol; } ); 
    std::cout << "# Total number of mols to pack: " << totalNumberOfMols << std::endl;
}


Ingredient* IngradientsDispatcher::pickIngredient()
{
    Ingredient* ingr;      
    unsigned int ingrInd;
    double threshProb;
    if (pickWeightedIngr){ 
        if (thresholdPriorities[0] == 2){
            //# Graham here: Walk through -priorities first
            ingr = activeIngr[0];
        }else{
            //#prob = uniform(vRangeStart,1.0)  #Graham 9/21/11 This is wrong...vRangeStart is the point index, need active list i.e. thresholdPriority to be limited
            double prob = uniform(generator);
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

void IngradientsDispatcher::dropIngredient( Ingredient *ingr )
{
    ingr->completion = 1.0;
    vRangeStart = vRangeStart + normalizedPriorities[0];
    std::cout << "#drop ingredient " << ingr->name << " " << ingr->nbMol << " " << ingr->completion << " " << ingr->counter << " "<< ingr->rejectionCounter <<std::endl;
    //update priorities will regenerate activeIngr based on completion field
    updatePriorities();
}

bool IngradientsDispatcher::hasNextIngradient() const {
    if (activeIngr.empty()) { 
        std::cout << "#broken by no more ingredient Done!!!****\n";
        return false;
    }
    if (vRangeStart > 1.0) {
        std::cout << "#broken by vRange and hence Done!!!****\n";
        return false;
    }
    return true;
}

void IngradientsDispatcher::prepareIngredient()
{
    getSortedActiveIngredients();
    std::cout << "#len(allIngredients) " << ingredients.size() << std::endl;
    std::cout << "#len(activeIngr0) " << activeIngr0.size() << std::endl;
    std::cout << "#len(activeIngr12) " << activeIngr12.size() << std::endl;

    calculateThresholdAndNormalizedPriorities();
}

void IngradientsDispatcher::updatePriorities()
{
    vRangeStart = vRangeStart + normalizedPriorities[0];
    //# Start of massive overruling section from corrected thesis file of Sept. 25, 2012    
    prepareIngredient();
}

void IngradientsDispatcher::calculateThresholdAndNormalizedPriorities()
{
    normalizedPriorities.clear();
    thresholdPriorities.clear();
    //# Graham- Once negatives are used, if picked random# 
    //# is below a number in this list, that item becomes 
    //#the active ingredient in the while loop below
    normalizedPriorities.resize(activeIngr0.size(), 0.0);
    if (pickWeightedIngr)
        thresholdPriorities.resize(activeIngr0.size(), 2.0);   

    double totalPriorities = std::accumulate(activeIngr12.begin(), activeIngr12.end(), 0.0, 
        [](double accumulator , Ingredient* ingr) { return accumulator + ingr->packingPriority; });
    std::cout << "#totalPriorities " << totalPriorities << std::endl;

    if (totalPriorities == 0)
    {
        normalizedPriorities.resize(normalizedPriorities.size() + activeIngr12.size(), 0.0);
        thresholdPriorities.resize(thresholdPriorities.size() + activeIngr12.size(), 0.0);  
    }
    else
    {
        double previousThresh = 0;                
        for(Ingredient * ingr : activeIngr12) {
            const double np = ingr->packingPriority/totalPriorities;            
            if (DEBUG)  std::cout << "#np is "<< np << " pp is "<< ingr->packingPriority << " tp is " << np + previousThresh << std::endl;

            normalizedPriorities.push_back(np);
            previousThresh += np;
            thresholdPriorities.push_back(previousThresh);
        }    
    }

    activeIngr = activeIngr0;
    activeIngr.insert(activeIngr.end(), activeIngr12.begin(), activeIngr12.end());
}

void IngradientsDispatcher::getSortedActiveIngredients()
{
    //pirorities120 is ot used ??
    std::vector<Ingredient*> ingr1;  // given priorities pass ptr ?
    std::vector<double> priorities1;
    std::vector<Ingredient*> ingr2;  // priority = 0 or none and will be assigned based on complexity
    std::vector<double> priorities2;
    std::vector<Ingredient*> ingr0;  // negative values will pack first in order of abs[packingPriority]
    std::vector<double> priorities0;
    Ingredient* lowestIng; 
    Ingredient* ing; 
    double r=0.0;
    double np=0.0;      
    for(unsigned i = 0; i < ingredients.size(); ++i) { 
        ing = &ingredients[i];
        if (ing->completion >= 1.0) continue;// # ignore completed ingredients
        //std::cout << "ingr " << ing->name << " " << ing->completion << std::endl;
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
    double lowestPriority;
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
        np = double(r)/double(totalRadii) * lowestPriority;
        //std::cout << "#packingPriority " << np << ' ' << r <<' '<< totalRadii<< ' ' << lowestPriority << std::endl;
        ing->packingPriority = np;
    }           

    activeIngr0 = ingr0;
    activeIngr12 = ingr1;
    activeIngr12.insert(activeIngr12.end(), ingr2.begin(), ingr2.end());
}
