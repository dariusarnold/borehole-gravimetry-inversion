set(CMAKE_CXX_STANDARD 17)

# enable to show complite compilation commands
set(CMAKE_VERBOSE_MAKEFILE ON)

# set compiler warnings
add_compile_options(-Wall -Wextra -pedantic -Weffc++ -Wshadow)

# include Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)


# Create the inversion library that I want to wrap
add_library(GravInvLibrary libGravimetryInversion/GravimetryInversion.h)
# Enable -fPIC
set_target_properties(GravInvLibrary PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
# get the sources from the CMakeLists.txt in the subdirectory
add_subdirectory(libGravimetryInversion)

# add CMakeLists.txt that createes that python bindings for the GravInvLibrary
add_subdirectory(python-bindings)

# add these exectables that run the program and save the output to file
add_executable(Interpolation main_interpolator.cpp)
add_executable(Inversion main_inversion.cpp)
add_executable(Inversion_Error main_inversion_with_errors.cpp)

# Link liberaries
target_link_libraries(GravInvLibrary PUBLIC Eigen3::Eigen)
#link library for filesystem support
target_link_libraries(GravInvLibrary PUBLIC stdc++fs)
target_link_libraries(Interpolation PRIVATE GravInvLibrary)
target_link_libraries(Inversion PRIVATE GravInvLibrary)
target_link_libraries(Inversion_Error PRIVATE GravInvLibrary)

# add include directories as interface to the library so targets linking against it will find the headers
target_include_directories(GravInvLibrary INTERFACE libGravimetryInversion)