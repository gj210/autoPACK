#define BOOST_TEST_MODULE PrioritiesTesting

#if !defined( WIN32 )
    #define BOOST_TEST_DYN_LINK
#endif
#include <boost/test/unit_test.hpp>

#include <string>
#include "../src/BigGrid.h"
#include "../src/Types.h"

float stepsize = float(15*1.1547); 

static Ingradient makeIngradient(std::string name, float minRadius, float packingPriority, float completition) {
    Ingradient ing1;    
    ing1.minRadius=minRadius;
    ing1.packingPriority=packingPriority;
    ing1.completion=completition;
    ing1.name = name;
    return ing1;
};

BOOST_AUTO_TEST_CASE( PrioritiestSortingTest ) {
    const openvdb::Vec3d ibotleft(-3,-3,-3);
    const openvdb::Vec3d iupright(3,3,3);
    big_grid grid(2,ibotleft,iupright,123);    
    
    std::vector<Ingradient> ingradients;    
    ingradients.push_back(makeIngradient("r1p0.1c0" , 1,  0.1, 0));
    ingradients.push_back(makeIngradient("r1p0.2c0" , 1,  0.2, 0));
    ingradients.push_back(makeIngradient("r1p0c0"   , 1,    0, 0));
    ingradients.push_back(makeIngradient("r1p-0.1c0", 1, -0.1, 0));
    ingradients.push_back(makeIngradient("r1p-0.2c0", 1, -0.2, 0));
    
    ingradients.push_back(makeIngradient("r1p0.1c0.1" , 1,  0.1, 0.1));
    ingradients.push_back(makeIngradient("r1p0.2c0.1" , 1,  0.2, 0.1));
    ingradients.push_back(makeIngradient("r1p0c0.1"   , 1,    0, 0.1));
    ingradients.push_back(makeIngradient("r1p-0.1c0.1", 1, -0.1, 0.1));
    ingradients.push_back(makeIngradient("r1p-0.2c0.1", 1, -0.2, 0.1));
    
    ingradients.push_back(makeIngradient("r2p0.1c0" , 2,  0.1, 0));
    ingradients.push_back(makeIngradient("r2p0.2c0" , 2,  0.2, 0));
    ingradients.push_back(makeIngradient("r2p0c0"   , 2,    0, 0));
    ingradients.push_back(makeIngradient("r2p-0.1c0", 2, -0.1, 0));
    ingradients.push_back(makeIngradient("r2p-0.2c0", 2, -0.2, 0));
    
    grid.setIngredients(ingradients);
    grid.prepareIngredient();

    std::vector<Ingradient> expected;
    expected.push_back(makeIngradient("r2p-0.1c0"  , 2,-0.1 ,   0));
    expected.push_back(makeIngradient("r1p-0.1c0.1", 1,-0.1 , 0.1));
    expected.push_back(makeIngradient("r1p-0.1c0"  , 1,-0.1 ,   0));
    expected.push_back(makeIngradient("r2p-0.2c0"  , 2,-0.2 ,   0));
    expected.push_back(makeIngradient("r1p-0.2c0.1", 1,-0.2 , 0.1));
    expected.push_back(makeIngradient("r1p-0.2c0"  , 1,-0.2 ,   0));
    expected.push_back(makeIngradient("r2p0.1c0"   , 2, 0.1 ,   0));
    expected.push_back(makeIngradient("r1p0.1c0.1" , 1, 0.1 , 0.1));
    expected.push_back(makeIngradient("r1p0.1c0"   , 1, 0.1 ,   0));
    expected.push_back(makeIngradient("r2p0.2c0"   , 2, 0.2 ,   0));
    expected.push_back(makeIngradient("r1p0.2c0.1" , 1, 0.2 , 0.1));
    expected.push_back(makeIngradient("r1p0.2c0"   , 1, 0.2 ,   0));
    expected.push_back(makeIngradient("r1p0c0.1"   , 1, 0.05, 0.1));
    expected.push_back(makeIngradient("r1p0c0"     , 1, 0.05,   0));
    expected.push_back(makeIngradient("r2p0c0"     , 2, 0.1 ,   0));

    for(int i=0; i<grid.activeIngr.size(); ++i) {
            BOOST_CHECK_EQUAL(expected[i].name, grid.activeIngr[i]->name);
            BOOST_CHECK_EQUAL(expected[i].packingPriority, grid.activeIngr[i]->packingPriority);
    }
    
    //for( Ingradient * ing : grid.activeIngr) {
        //std::cout << "expected.push_back(makeIngradient(\"" << ing->name << "\", " << ing->minRadius << "," << ing->packingPriority << "," << ing->completion << "));" << std::endl;
    //}
}
