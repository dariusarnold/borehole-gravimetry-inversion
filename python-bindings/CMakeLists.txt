#cmake_minimum_required(VERSION 3.10)
#project(GravityInversion)

# include pybind11
set(PYBIND11_CPP_STANDARD -std=gnu++17)
# set the correct python version to use
set(PYBIND11_PYTHON_VERSION 3.6)
find_package(pybind11 REQUIRED)

# Create the python extension module
pybind11_add_module(PyLibrary PyLibrary.cpp)

target_link_libraries(PyLibrary PRIVATE Eigen3::Eigen GravInvLibrary)