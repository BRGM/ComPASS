cmake_minimum_required(VERSION 3.2)

project(ComPASS Fortran C CXX)

# don't set CMAKE_BUILD_TYPE here, it is already in the cache
# This statically specifies what build type (configuration) will be built in this build tree. Possible values are empty, Debug, Release, RelWithDebInfo and MinSizeRel. This variable is only meaningful to single-configuration generators (such as make and Ninja)
# set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Release)

# make sure our own modules are found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules")

include(ComPASSCMakeUtils)

if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(FATAL_ERROR
"-------------------------------
In-source compilation forbidden.
Build a directory build: mkdir build.
Then: cd build; cmake ..
--------------------------------")
endif(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)

## FIXME: to be removed ! not portable
# set(CMAKE_Fortran_COMPILER mpifort)
# set(CMAKE_C_COMPILER mpicc)
# set(CMAKE_CXX_COMPILER mpicxx)

########### Configuration options

option(DEF_THERMIQUE "Thermique" ON)
if(DEF_THERMIQUE)
  add_definitions(-D_THERMIQUE_)
endif(DEF_THERMIQUE)

option(DEF_VISU "Visualization with VTK" ON)
if(DEF_VISU)
  add_definitions(-D_VISU_)
endif(DEF_VISU)

option(DEF_HDF5 "Checkpointing with HDF5" OFF)

option(DEF_DISPMODULE "Third Party debug using displaymodule" OFF)
if(DEF_DISPMODULE)
  add_definitions(-D_DISPMODULE_)
endif(DEF_DISPMODULE)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-D_DEBUG_)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

get_directory_property( DirDefs COMPILE_DEFINITIONS )
message(STATUS "COMPILE_DEFINITIONS ${DirDefs}" )

##################

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/bin/fortran CACHE PATH "Single Directory for all fortran modules.")

#CHECKME ???? -> portability + test constants are compiled as double
set(CMAKE_Fortran_FLAGS "-I${LIB_PREFIX}/include/shared -fdefault-real-8 -fdefault-double-8" CACHE STRING "Fortran FLAGS") # make sure that the constants (ex. 1.0) are compiled as double (ex. 1.d0) -fdefault-real-8 -fdefault-double-8
# add_definitions(-fdefault-real-8 -fdefault-double-8)

#CHECKME: is this the place ? portability -> Windows ?
add_definitions(-O0 -g)

## CHECKME: CMake -> cmake ?
set(CMAKE_FIND_DIR ${CMAKE_SOURCE_DIR}/CMake)

find_package(MPI REQUIRED)
find_package(METIS REQUIRED)
find_package(PETSc REQUIRED)
# cache PETSC version
set(PETSC_VERSION ${PETSC_VERSION} CACHE INTERNAL "PETSC version") 
if(${PETSC_VERSION} VERSION_GREATER 3.6)
    message(FATAL_ERROR "PETSC version should not be greater than 3.5.x")
endif()
find_package(LAPACK REQUIRED)
# vtk
# set(VTK_DIR ${LIB_PREFIX}/include/vtk-7.1)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# hdf5
if(DEF_HDF5)
  add_definitions(-D_HDF5_)
set(HDF5_USE_STATIC_LIBRARIES OFF)
option(HDF5_PREFER_PARALLEL "HDF5 prefers parallel" ON)
  # FIXME: If option is changed find_package will not reload the correct hdf5 version
  # FIXME: inforce parallel ? h5pset_fapl_mpio_f_ is needed ?
  find_package(HDF5 REQUIRED COMPONENTS Fortran Fortran_HL)
  if(HDF5_FOUND)
      message(STATUS "HDF5 version: " ${HDF5_VERSION})
  else()
      message(FATAL_ERROR "HDF5 was required but was not found on the system, CMake will exit." )
  endif()
endif(DEF_HDF5)

enable_testing()

# CHECKME: Should languages be enabled here or at the beginning as project properties
# enable_language(Fortran)
# enable_language(C)
# enable_language(CXX)

add_subdirectory(src)