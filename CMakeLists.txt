cmake_minimum_required(VERSION 3.12.0)
project(navier_stokes LANGUAGES CXX C)
find_package(MPI REQUIRED)

include(cmake-common.cmake)

add_executable(navier_stokes src/main.cpp src/NavierStokes.cpp)
deal_ii_setup_target(navier_stokes)
add_executable(navier_stokes3D src/main3D.cpp src/NavierStokes3D.cpp)
deal_ii_setup_target(navier_stokes3D)