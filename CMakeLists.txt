cmake_minimum_required(VERSION 3.10)
project(HierarchicalFaceClustering)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

find_package(Eigen3 3.3 REQUIRED)
if(EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    message(STATUS "Eigen is not found.")
endif()

set(LIBIGL_ROOT_DIR "$ENV{LIBIGL_ROOT}")
if(LIBIGL_ROOT_DIR)
    include_directories("${LIBIGL_ROOT_DIR}/include")
else()
    message(FATAL_ERROR "Cannot find libigl at LIBIGL_ROOT")
endif()

add_executable(HierarchicalFaceClustering main.cpp ClusteringAlgorithm.cpp ClusteringAlgorithm.h)
