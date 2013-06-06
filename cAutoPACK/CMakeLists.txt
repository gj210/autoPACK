CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

cmake_policy(SET CMP0015 OLD)

PROJECT(cAutoPack)

SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

find_path(ILMBASE_INCLUDE_DIR NAMES Half.h OpenEXR/half.h OpenEXR/Half.h PATH_SUFFIXES OpenEXR/include/ REQUIRED)
message(STATUS "OpenEXR includes found in: ${ILMBASE_INCLUDE_DIR}")
find_library(ILMBASE_LIBRARY Half PATH_SUFFIXES OpenEXR/lib REQUIRED)
message(STATUS "OpenEXR library found in: ${ILMBASE_LIBRARY}")

find_path(TBB_INCLUDE_DIR tbb/tbb.h PATH_SUFFIXES tbb/include REQUIRED)
message(STATUS "OpenEXR includes found in: ${TBB_INCLUDE_DIR}")
find_library(TBB_LIBRARY NAMES tbb PATH_SUFFIXES tbb/lib REQUIRED)
message(STATUS "TBB library found in: ${TBB_LIBRARY}")

find_path(PUGIXML_INCLUDE_DIR pugixml.hpp PATH_SUFFIXES pugixml/include REQUIRED)
message(STATUS "pugiXML includes found in: ${PUGIXML_INCLUDE_DIR}")
find_library(PUGIXML_LIBRARY NAMES pugixml PATH_SUFFIXES pugixml/lib REQUIRED)
message(STATUS "pugiXML library found in: ${PUGIXML_LIBRARY}")

find_path(OPENVDB_INCLUDE_DIR openvdb/openvdb.h PATH_SUFFIXES OpenVDB/include REQUIRED)
message(STATUS "OpenVDB includes found in: ${OPENVDB_INCLUDE_DIR}")
find_library(OPENVDB_LIBRARY NAMES openvdb vdb PATH_SUFFIXES OpenVDB/lib REQUIRED)
message(STATUS "OpenVDB library found in: ${OPENVDB_LIBRARY}")

find_package(Boost REQUIRED)
FIND_PACKAGE(ZLIB REQUIRED)
	
add_executable(autoPack autopack.cpp)

include_directories(${Boost_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS}
	${ILMBASE_INCLUDE_DIR} ${TBB_INCLUDE_DIR}
	${PUGIXML_INCLUDE_DIR} ${OPENVDB_INCLUDE_DIR})

target_link_libraries(autoPack ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${ILMBASE_LIBRARY} ${TBB_LIBRARY} ${PUGIXML_LIBRARY} ${OPENVDB_LIBRARY})

if(MSVC)
	ADD_DEFINITIONS(-DNOMINMAX)	
	ADD_DEFINITIONS(-DZLIB_WINAPI)
	ADD_DEFINITIONS(-D__TBB_NO_IMPLICIT_LINKAGE)
endif()
if (CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11")
    if (APPLE)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif()
