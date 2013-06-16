/*
        # Graham here- In the new setup, priority is infinite with abs[priority] increasing (+)
        # An ingredients with (-) priority will pack from greatest abs[-priority] one at a time
        #     to lease abs[-priority]... each ingredient will attempt to reach its molarity
        #     before moving on to the next ingredient, and all (-) ingredients will try to
        #     deposit before other ingredients are tested.
        # An ingredient with (+) priority will recieve a weighted value based on its abs[priority]
        #     e.g. an ingredient with a priority=10 will be 10x more likely to be picked than
        #     an ingredient with a priority=1.
        # An ingredient with the default priority=0 will recieve a weighted value based on its
        #     complexity. (currently complexity = minRadius), thus a more 'complex' ingredient
        #     will more likely try to pack before a less 'complex' ingredient.
        #     IMPORTANT: the +priority list does not fully mix with the priority=0 list, but this
        #     should be an option... currently, the priority=0 list is normalized against a range
        #     up to the smallest +priority ingredient and appended to the (+) list
        # TODO: Add an option to allow + ingredients to be weighted by assigned priority AND complexity
        #     Add an option to allow priority=0 ingredients to fit into the (+) ingredient list
        #       rather than appending to the end.
        #     Even better, add an option to set the max priority for the 0list and then plug the results
        #       into the (+) ingredient list.
        #     Get rid of the (-), 0, (+) system and recreate this as a new flag and a class function
        #        so we can add multiple styles of sorting and weighting systems.
        #     Make normalizedPriorities and thresholdPriorities members of Ingredient class to avoid
        #        building these arrays.

*/
#define BOOST_TEST_MODULE PrioritiesTesting

#if !defined( WIN32 )
    #define BOOST_TEST_DYN_LINK
#endif
#include <boost/test/unit_test.hpp>

#include <string>
#include "../src/IngredientsDispatcher.h"
#include "../src/Types.h"

float stepsize = float(15*1.1547); 

static Ingredient makeIngredient(std::string name, float minRadius, float packingPriority, float completition) {
    Ingredient ing1;    
    ing1.minRadius=minRadius;
    ing1.packingPriority=packingPriority;
    ing1.completion=completition;
    ing1.name = name;
    return ing1;
};

BOOST_AUTO_TEST_CASE( PrioritiestSortingTest ) {
        
    std::vector<Ingredient> ingredients;    
    ingredients.push_back(makeIngredient("r1p0.1c0" , 1,  0.1, 0));
    ingredients.push_back(makeIngredient("r1p0.2c0" , 1,  0.2, 0));
    ingredients.push_back(makeIngredient("r1p0c0"   , 1,    0, 0));
    ingredients.push_back(makeIngredient("r1p-0.1c0", 1, -0.1, 0));
    ingredients.push_back(makeIngredient("r1p-0.2c0", 1, -0.2, 0));
    
    ingredients.push_back(makeIngredient("r1p0.1c0.1" , 1,  0.1, 0.1));
    ingredients.push_back(makeIngredient("r1p0.2c0.1" , 1,  0.2, 0.1));
    ingredients.push_back(makeIngredient("r1p0c0.1"   , 1,    0, 0.1));
    ingredients.push_back(makeIngredient("r1p-0.1c0.1", 1, -0.1, 0.1));
    ingredients.push_back(makeIngredient("r1p-0.2c0.1", 1, -0.2, 0.1));
    
    ingredients.push_back(makeIngredient("r2p0.1c0" , 2,  0.1, 0));
    ingredients.push_back(makeIngredient("r2p0.2c0" , 2,  0.2, 0));
    ingredients.push_back(makeIngredient("r2p0c0"   , 2,    0, 0));
    ingredients.push_back(makeIngredient("r2p-0.1c0", 2, -0.1, 0));
    ingredients.push_back(makeIngredient("r2p-0.2c0", 2, -0.2, 0));

    IngradientsDispatcher dispatcher(ingredients, 300 /*not important in test*/, /* seed */ 123);    

    //assertion for active inrgedent size
    BOOST_CHECK_EQUAL(dispatcher.ingredients.size(),ingredients.size());
    BOOST_CHECK_EQUAL(dispatcher.activeIngr0.size(),6) ;    //negative packingPriority
    BOOST_CHECK_EQUAL(dispatcher.activeIngr12.size(),9) ;   //zero or positive packingPriority
    BOOST_CHECK_EQUAL(dispatcher.activeIngr.size(),ingredients.size());//initially all ingredient should be active

    std::vector<Ingredient> expected; //excpectation after initialization

    expected.push_back(makeIngredient("r2p-0.1c0"  , 2,-0.1 ,   0));
    expected.push_back(makeIngredient("r1p-0.1c0.1", 1,-0.1 , 0.1));
    expected.push_back(makeIngredient("r1p-0.1c0"  , 1,-0.1 ,   0));
    expected.push_back(makeIngredient("r2p-0.2c0"  , 2,-0.2 ,   0));
    expected.push_back(makeIngredient("r1p-0.2c0.1", 1,-0.2 , 0.1));
    expected.push_back(makeIngredient("r1p-0.2c0"  , 1,-0.2 ,   0));
    expected.push_back(makeIngredient("r2p0.1c0"   , 2, 0.1 ,   0));
    expected.push_back(makeIngredient("r1p0.1c0.1" , 1, 0.1 , 0.1));
    expected.push_back(makeIngredient("r1p0.1c0"   , 1, 0.1 ,   0));
    expected.push_back(makeIngredient("r2p0.2c0"   , 2, 0.2 ,   0));
    expected.push_back(makeIngredient("r1p0.2c0.1" , 1, 0.2 , 0.1));
    expected.push_back(makeIngredient("r1p0.2c0"   , 1, 0.2 ,   0));
    expected.push_back(makeIngredient("r1p0c0.1"   , 1, 0.05, 0.1));
    expected.push_back(makeIngredient("r1p0c0"     , 1, 0.05,   0));
    expected.push_back(makeIngredient("r2p0c0"     , 2, 0.1 ,   0));

    //test the order and the new packing priority computed for priority of 0 ->radius based
    for(int i=0; i<dispatcher.activeIngr.size(); ++i) {
            std::cout << dispatcher.activeIngr[i]->name << std::endl;
            BOOST_CHECK_EQUAL(expected[i].name, dispatcher.activeIngr[i]->name);
            BOOST_CHECK_EQUAL(expected[i].packingPriority, dispatcher.activeIngr[i]->packingPriority);
    }
    //there is no more ingredient with priority = 0.

    //Now drop ingredient. Negative ingredient should be pack first
    for (int i = 0 ; i< 6;i++ ) dispatcher.dropIngredient(dispatcher.activeIngr[0]);

    std::vector<Ingredient> expected_after_drop;
    expected_after_drop.push_back(makeIngredient("r1p0c0.1"   , 1, 0.05, 0.1));
    expected_after_drop.push_back(makeIngredient("r1p0c0"     , 1, 0.05,   0));
    expected_after_drop.push_back(makeIngredient("r2p0.1c0"   , 2, 0.1 ,   0));
    expected_after_drop.push_back(makeIngredient("r2p0c0"     , 2, 0.1 ,   0));
    expected_after_drop.push_back(makeIngredient("r1p0.1c0.1" , 1, 0.1 , 0.1));
    expected_after_drop.push_back(makeIngredient("r1p0.1c0"   , 1, 0.1 ,   0));
    expected_after_drop.push_back(makeIngredient("r2p0.2c0"   , 2, 0.2 ,   0));
    expected_after_drop.push_back(makeIngredient("r1p0.2c0.1" , 1, 0.2 , 0.1));
    expected_after_drop.push_back(makeIngredient("r1p0.2c0"   , 1, 0.2 ,   0));
    //test the order and the new packing priority computed for priority of 0 ->radius based
    for(int i=0; i<dispatcher.activeIngr.size(); ++i) {
            //std::cout << grid.activeIngr[i]->name << " " << expected_after_drop[i].name << std::endl;
            //std::cout << grid.activeIngr[i]->packingPriority << " " << expected_after_drop[i].packingPriority << std::endl;
            BOOST_CHECK_EQUAL(expected_after_drop[i].name, dispatcher.activeIngr[i]->name);
            BOOST_CHECK_EQUAL(expected_after_drop[i].packingPriority, dispatcher.activeIngr[i]->packingPriority);
    }
    //Continue to Drop and see how it goes, also with change of completion
    //after negative droppin the inrdent picking is based on the weight.
    //r2p0.2c0 should be the mot probable
    dispatcher.dropIngredient(dispatcher.activeIngr[6]);//
    dispatcher.activeIngr[1]->completion = 0.6f;

    std::vector<Ingredient> expected_after_drop_and_completion;
    expected_after_drop_and_completion.push_back(makeIngredient("r1p0c0.1"   , 1, 0.05, 0.1));//shouldn't be r1p0c0 here ? verify with Graham if behaviour is correct here
    expected_after_drop_and_completion.push_back(makeIngredient("r1p0c0"     , 1, 0.05, 0.6));
    expected_after_drop_and_completion.push_back(makeIngredient("r2p0.1c0"   , 2, 0.1 ,   0));
    expected_after_drop_and_completion.push_back(makeIngredient("r2p0c0"     , 2, 0.1 ,   0));
    expected_after_drop_and_completion.push_back(makeIngredient("r1p0.1c0.1" , 1, 0.1 , 0.1));
    expected_after_drop_and_completion.push_back(makeIngredient("r1p0.1c0"   , 1, 0.1 ,   0));
    expected_after_drop_and_completion.push_back(makeIngredient("r1p0.2c0.1" , 1, 0.2 , 0.1));
    expected_after_drop_and_completion.push_back(makeIngredient("r1p0.2c0"   , 1, 0.2 ,   0));

    for(int i=0; i<dispatcher.activeIngr.size(); ++i) {
            std::cout << dispatcher.activeIngr[i]->name << " " << expected_after_drop_and_completion[i].name << std::endl;
            std::cout << dispatcher.activeIngr[i]->packingPriority << " " << expected_after_drop_and_completion[i].packingPriority << std::endl;
            BOOST_CHECK_EQUAL(expected_after_drop_and_completion[i].name, dispatcher.activeIngr[i]->name);
            BOOST_CHECK_EQUAL(expected_after_drop_and_completion[i].packingPriority, dispatcher.activeIngr[i]->packingPriority);
    }

    //for( Ingradient * ing : grid.activeIngr) {
        //std::cout << "expected.push_back(makeIngredient(\"" << ing->name << "\", " << ing->minRadius << "," << ing->packingPriority << "," << ing->completion << "));" << std::endl;
    //}
}
