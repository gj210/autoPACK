#pragma once


#ifdef _MSC_VER
#pragma warning(disable:4146)
#pragma warning(disable:4503) // OpenVDB "warning decorated name length exceeded, name was truncated"
#define _SCL_SECURE_NO_WARNINGS
#endif

#include <vector>
//Disable warnings from openvdb in Visual Studio
#ifdef _MSC_VER
    #pragma warning(push, 0)   
#endif

/*
openvdb includes
*/
#include <openvdb/openvdb.h>
#include <openvdb/Types.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tree/Tree.h>
#include <openvdb/tools/GridTransformer.h>

#ifdef _MSC_VER
    #pragma warning(pop)
#endif

