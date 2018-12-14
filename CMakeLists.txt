cmake_minimum_required(VERSION 3.10)
project(GravityInversion)

set(CMAKE_CXX_STANDARD 17)

# enable to show complite compilation commands
#set(CMAKE_VERBOSE_MAKEFILE ON)

# set compiler warnings
add_compile_options(-Wall -Wextra -pedantic -Weffc++ -Wshadow)

# include Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)


# Create the inversion library that I want to wrap
add_library(GravInvLibrary "")
# Enable -fPIC
set_target_properties(GravInvLibrary PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

#add_subdirectory(libGravimetryInversion)
include(libGravimetryInversion/CMakeLists.txt)

# Link liberaries
target_link_libraries(GravInvLibrary PRIVATE Eigen3::Eigen)
#link library for filesystem support
target_link_libraries(GravInvLibrary PRIVATE stdc++fs)

add_subdirectory(python-bindings)