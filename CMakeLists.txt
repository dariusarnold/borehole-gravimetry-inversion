cmake_minimum_required(VERSION 3.13)
project(GravityInversion)

# enable to show complite compilation commands
#set(CMAKE_VERBOSE_MAKEFILE ON)

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
endif (UNIX AND NOT CYGWIN)

add_subdirectory(cppGravInv)