cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)
 
project(2DHolstein CXX)

find_package(Eigen3 REQUIRED)

message(STATUS ${EIGEN3_INCLUDE_DIR})

add_executable(main holsteinDQMC.cpp)
target_include_directories(main PUBLIC ${EIGEN3_INCLUDE_DIR})

