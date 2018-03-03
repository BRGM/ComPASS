cmake_minimum_required(VERSION 3.2)

project(MeshTools CXX)

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
set( CMAKE_CXX_STANDARD 14)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

# Check if meshtools is being used directly or via add_subdirectory
set(MESHTOOLS_MASTER_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(MESHTOOLS_MASTER_PROJECT ON)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(MESHTOOLS_CONFIG_DIRECTORY ${PROJECT_SOURCE_DIR}/config CACHE PATH "path to configuration scripts")
set(MESHTOOLS_PYTHON_DIRECTORY ${PROJECT_SOURCE_DIR} CACHE PATH "python package directory")

include_directories(mapbox/include) 

set(PYBIND11_PYTHON_VERSION 3 CACHE STRING "Python version to use for compiling modules")
# Relying on PYBIND11_VERSION_MAJOR variable which is cached to detect pybind11 installation
if(NOT (DEFINED PYBIND11_VERSION_MAJOR))
  find_package(pybind11 REQUIRED)
  if(${pybind11_VERSION} VERSION_LESS 2.2)
    message(FATAL_ERROR "pybind11 version should 2.2.x at least")
  endif(${pybind11_VERSION} VERSION_LESS 2.2)
endif(NOT (DEFINED PYBIND11_VERSION_MAJOR))
# The following test is not mandatory as python interpreter comes as a by-product of requiring pybind11
if(NOT ${PYTHONINTERP_FOUND})
  message(FATAL_ERROR "python interpreter is mandatory, CMake will stop.")
endif(NOT ${PYTHONINTERP_FOUND})
# The following test is not mandatory as python interpreter comes as a by-product of requiring pybind11
if(NOT ${PYTHONLIBS_FOUND})
  message(FATAL_ERROR "python libraries are mandatory, CMake will stop.")
endif(NOT ${PYTHONLIBS_FOUND})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if(CGAL_FOUND)
    # include helper file
    include( ${CGAL_USE_FILE} )
    # Boost and its components
    # FIXME should fail silently without including CGAL
    find_package( Boost REQUIRED )
    if ( NOT Boost_FOUND )
      message(STATUS "This project requires the Boost library, and will not be compiled.")
      return()  
    endif()
    # Find Eigen3 (requires 3.1.0 or greater)
    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
        include( ${EIGEN3_USE_FILE} )
    endif()
    set(MESHTOOLS_CGAL_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/cgal)
    add_subdirectory(${MESHTOOLS_CGAL_DIRECTORY})
    set(MESHTOOLS_USES_CGAL TRUE)
    add_definitions(-DMESHTOOLS_USES_CGAL)
else(CGAL_FOUND)
    set(MESHTOOLS_USES_CGAL FALSE)
endif(CGAL_FOUND)

pybind11_add_module(_MeshTools meshtools-module.cpp meshtools-wrapper.cpp)
set_target_properties(_MeshTools
    PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    LIBRARY_OUTPUT_DIRECTORY ${MESHTOOLS_PYTHON_DIRECTORY}
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${MESHTOOLS_PYTHON_DIRECTORY}
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${MESHTOOLS_PYTHON_DIRECTORY}
)

if(${MESHTOOLS_MASTER_PROJECT})
    add_executable(test-meshtools test-meshtools.cpp)
    pybind11_add_module(consumer consumer-module.cpp)
    set_target_properties(consumer
        PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        LIBRARY_OUTPUT_DIRECTORY ${MESHTOOLS_PYTHON_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY_DEBUG ${MESHTOOLS_PYTHON_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY_RELEASE ${MESHTOOLS_PYTHON_DIRECTORY}
    )
endif(${MESHTOOLS_MASTER_PROJECT})


# make module available through user site packages
install(CODE "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" add_to_site.py MeshTools \"${MESHTOOLS_PYTHON_DIRECTORY}\" WORKING_DIRECTORY \"${MESHTOOLS_CONFIG_DIRECTORY}\")")
