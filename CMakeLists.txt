cmake_minimum_required(VERSION 3.10)
project(GravityInversion)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(Programm main.cpp utils.h utils.cpp MeasurementData.h MeasurementData.cpp GravimetryInversion.h GravimetryInversion.cpp )
target_link_libraries(Programm Eigen3::Eigen)