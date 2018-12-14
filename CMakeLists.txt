cmake_minimum_required(VERSION 3.10)
project(GravityInversion)

set(CMAKE_CXX_STANDARD 17)

# enable to show complite compilation commands
#set(CMAKE_VERBOSE_MAKEFILE ON)

# set compiler warnings
add_compile_options(-Wall -Wextra -pedantic -Weffc++ -Wshadow)

# include Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# include pybind11
set(PYBIND11_CPP_STANDARD -std=gnu++17)
# set the correct python version to use
set(PYBIND11_PYTHON_VERSION 3.6)
find_package(pybind11 REQUIRED)

# Create the inversion library that I want to wrap
add_library(Library source/MeasurementData.cpp source/MeasurementData.h source/GravimetryInversion.h source/FileIO.cpp source/FileIO.h source/Result.cpp source/Result.h source/Library.cpp source/Library.h source/Norms.cpp source/Norms.h source/ErrorNorms.cpp source/utils.h source/utils.cpp)
# Enable -fPIC
set_target_properties(Library PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

# Create the python extension module
pybind11_add_module(PyLibrary source/PyLibrary.cpp)

# Link liberaries
target_link_libraries(Library PRIVATE Eigen3::Eigen)
#link library for filesystem support
target_link_libraries(Library PRIVATE stdc++fs)
target_link_libraries(PyLibrary PRIVATE Eigen3::Eigen Library)
