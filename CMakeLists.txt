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

# compile under windows
if (CYGWIN)
    set(Eigen3_DIR "C:\\Eigen3")
    # use system directive to supress warnings from Eigen library
    include_directories(SYSTEM "${Eigen3_DIR}")
    link_directories(${Eigen3_DIR})
endif (CYGWIN)
if (UNIX AND NOT CYGWIN)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
endif (UNIX AND NOT CYGWIN)

add_executable(Inversion main_inversion.cpp)
add_executable(Interpolation main_interpolator.cpp)

# set include directory so compiler finds my headers
include_directories(source)

# include the CMakeLists that specifiy the sources for the targets
include(source/CMakeLists.txt)

# Link liberaries
if (UNIX AND NOT CYGWIN)
    target_link_libraries(Inversion Eigen3::Eigen)
    target_link_libraries(Interpolation Eigen3::Eigen)
endif (UNIX AND NOT CYGWIN)
#link library for filesystem support
target_link_libraries(Inversion stdc++fs)
target_link_libraries(Interpolation stdc++fs)