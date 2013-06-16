#pragma once

#include "Ingredient.h"
#include <random>
#include <numeric>


class IngradientsDispatcher {
public: //Temporary public, should be private
    typedef std::vector<Ingredient>::const_iterator const_iterator;
    bool pickWeightedIngr;

    std::vector<Ingredient*> activeIngr;
    std::vector<Ingredient*> activeIngr0;
    std::vector<Ingredient*> activeIngr12;

    std::vector<float> normalizedPriorities;
    float vRangeStart;

    std::vector<float> thresholdPriorities; 

    std::vector<Ingredient> ingredients;//the list of sphere ingredient to pack


    std::uniform_real_distribution<float> uniform;
    std::default_random_engine & generator;

public:
    IngradientsDispatcher(std::vector<Ingredient> const & _ingredients,  unsigned num_points, unsigned int seed);    

    Ingredient* pickIngredient();

    void dropIngredient( Ingredient *ingr );

    bool hasNextIngradient() const;

    const_iterator begin() const { return ingredients.cbegin();}

    const_iterator end() const { return ingredients.cend();}

private:
    IngradientsDispatcher(const IngradientsDispatcher & that);
    void prepareIngredient();

    void updatePriorities();

    void calculateThresholdAndNormalizedPriorities();

    void getSortedActiveIngredients();

};