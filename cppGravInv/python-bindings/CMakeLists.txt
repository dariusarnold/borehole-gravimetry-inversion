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

# add install rule
install(TARGETS PyLibrary DESTINATION ${CMAKE_INSTALL_PREFIX}/pyGravInv)

# execute install after build
add_custom_command(
        TARGET PyLibrary
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target install
)