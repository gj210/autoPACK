#pragma once

#ifdef _MSC_VER
//Disable warnings from openvdb in Visual Studio
#pragma warning(push, 0)   
#endif

#include <openvdb/openvdb.h>
#include <openvdb/Types.h>
#include <openvdb/tools/GridTransformer.h>

#ifdef _MSC_VER
//Disable warnings from openvdb in Visual Studio
#pragma warning(pop)
#endif


namespace Geometric {
    openvdb::Vec3d calculateGeometricCenter( std::vector<openvdb::Vec3d>::const_iterator beginIter, std::vector<openvdb::Vec3d>::const_iterator endIter );


    openvdb::Vec3d transformPossition( openvdb::Vec3d const& position, openvdb::Vec3d const& offset, openvdb::math::Mat4d const& rotMatj );

    double countNormalizedDistanceToAllPositions( std::vector<openvdb::Vec3d> const& rpossitions, openvdb::Vec3d const& transformedPosition);
}