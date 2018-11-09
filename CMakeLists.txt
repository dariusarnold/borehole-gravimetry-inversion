cmake_minimum_required(VERSION 3.10)
project(GravityInversion)

set(CMAKE_CXX_STANDARD 17)

if (WIN32)
    set(EIGEN_DIR "C:\\Eigen3")
    include_directories(${EIGEN_DIR})
endif (WIN32)
 if (UNIX)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
endif (UNIX)

add_executable(Programm main.cpp utils.h utils.cpp MeasurementData.h MeasurementData.cpp GravimetryInversion.h GravimetryInversion.cpp FileWriter.cpp FileWriter.h)
if (UNIX)
    target_link_libraries(Programm Eigen3::Eigen)
endif (UNIX)