
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

    openvdb::Vec3d transformPossition( openvdb::Vec3d const& position, openvdb::Vec3d const& offset, openvdb::math::Mat4d const& rotMatj )
    {	
        openvdb::math::Transform::Ptr targetXform = openvdb::math::Transform::createLinearTransform();
        targetXform->preMult(rotMatj);
        targetXform->postTranslate(offset);
        openvdb::math::Mat4d mat = targetXform->baseMap()->getAffineMap()->getMat4();
        const openvdb::Vec3d pos = mat.transform(position);
        return pos;
    }

    double countNormalizedDistanceToAllPositions( std::vector<openvdb::Vec3d> const& rpossitions, openvdb::Vec3d const& transformedPosition)
    {
        if (rpossitions.empty())
            return 0;

        double sum = std::accumulate(std::begin(rpossitions), std::end(rpossitions), 0.0, [ &transformedPosition ](double acc, openvdb::Vec3d const& position) 
        { return acc + std::abs((transformedPosition - position).lengthSqr()); } );

        return sum/rpossitions.size();
    }
}