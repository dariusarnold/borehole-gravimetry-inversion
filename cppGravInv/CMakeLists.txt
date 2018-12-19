cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 17)

# set compiler warnings
add_compile_options(-Wall -Wextra -pedantic -Weffc++ -Wshadow)

# include Eigen
# compile under windows
if (CYGWIN)
    message(STATUS "Using manually installed Eigen library on windows")
    set(Eigen3_DIR "C:\\Eigen3\\")
    # use system directive to supress warnings from Eigen library
    include_directories(SYSTEM "${Eigen3_DIR}")
    link_directories(${Eigen3_DIR})
endif (CYGWIN)
if (UNIX AND NOT CYGWIN)
    message(STATUS "Using Eigen installed from package manager")
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
else()
    message(STATUS "On Unix but haven't found Eigen3")
endif (UNIX AND NOT CYGWIN)

# Create the inversion library that I want to wrap
add_library(GravInvLibrary libGravimetryInversion/GravimetryInversion.h)
# Enable -fPIC
set_target_properties(GravInvLibrary PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
# get the sources from the CMakeLists.txt in the subdirectory
add_subdirectory(libGravimetryInversion)

# add CMakeLists.txt that creates that python bindings for the GravInvLibrary
add_subdirectory(python-bindings)

# add these exectables that run the program and save the output to file
add_executable(Interpolation main_interpolator.cpp)
add_executable(Inversion main_inversion.cpp)
add_executable(Inversion_Error main_inversion_with_errors.cpp)

# Link libraries
if (UNIX AND NOT CYGWIN)
    target_link_libraries(GravInvLibrary PUBLIC Eigen3::Eigen)
endif (UNIX AND NOT CYGWIN)
#link library for filesystem support
target_link_libraries(GravInvLibrary PUBLIC stdc++fs)
target_link_libraries(Interpolation PRIVATE GravInvLibrary)
target_link_libraries(Inversion PRIVATE GravInvLibrary)
target_link_libraries(Inversion_Error PRIVATE GravInvLibrary)

# add include directories as interface to the library so targets linking against it will find the headers
target_include_directories(GravInvLibrary INTERFACE libGravimetryInversion)