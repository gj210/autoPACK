
#include "Types.h"
#include <numeric>

namespace Geometric {
    
    openvdb::Vec3d calculateGeometricCenter( std::vector<openvdb::Vec3d>::const_iterator beginIter, std::vector<openvdb::Vec3d>::const_iterator endIter )
    {
        auto geometricCenterOffset = std::accumulate(beginIter, endIter, openvdb::Vec3d::zero() , 
            [] (openvdb::Vec3d const& acc, openvdb::Vec3d const& item) { return item + acc; });

        geometricCenterOffset /= std::distance(beginIter,endIter);
        return geometricCenterOffset;
    }
}