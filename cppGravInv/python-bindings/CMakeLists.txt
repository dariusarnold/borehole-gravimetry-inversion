# include pybind11
set(PYBIND11_CPP_STANDARD -std=gnu++17)
# set the correct python version to use
set(PYBIND11_PYTHON_VERSION 3.6)

if (CYGWIN)
    # on windows, use local version
    add_subdirectory(pybind11)
endif (CYGWIN)
if (UNIX AND NOT CYGWIN)
    find_package(pybind11 REQUIRED)
endif (UNIX AND NOT CYGWIN)

# Create the python extension module
pybind11_add_module(PyLibrary PyLibrary.cpp)

if (UNIX AND NOT CYGWIN)
    target_link_libraries(PyLibrary PRIVATE Eigen3::Eigen GravInvLibrary)
else()
    target_link_libraries(PyLibrary PRIVATE GravInvLibrary ${PYTHON_LIBRARIES})
endif (UNIX AND NOT CYGWIN)

# add install rule
install(TARGETS PyLibrary DESTINATION ${CMAKE_INSTALL_PREFIX}/pyGravInv)

# execute install after build
add_custom_command(
        TARGET PyLibrary
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target install
)