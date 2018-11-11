cmake_minimum_required(VERSION 3.10)
project(GravityInversion)

set(CMAKE_CXX_STANDARD 17)

# enable to show complite compilation commands
#set(CMAKE_VERBOSE_MAKEFILE ON)

# set compiler warnings
add_compile_options(-Wall -Wextra -pedantic -Weffc++ -Wshadow)
# enable optimization for release build
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message("Activating optimizations (-O3) for ${CMAKE_BUILD_TYPE} build type")
    add_compile_options(-O3)
endif()

if (WIN32)
    set(EIGEN_DIR "C:\\Eigen3")
    include_directories(${EIGEN_DIR})
endif (WIN32)
 if (UNIX)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
endif (UNIX)

add_executable(Programm main.cpp utils.h utils.cpp MeasurementData.h MeasurementData.cpp GravimetryInversion.h GravimetryInversion.cpp FileWriter.cpp FileWriter.h Representant.cpp Representant.h)
if (UNIX)
    target_link_libraries(Programm Eigen3::Eigen)
endif (UNIX)
#link library for filesystem support
target_link_libraries(Programm stdc++fs)