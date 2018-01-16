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

include_directories(mapbox/include) 

set(PYBIND11_PYTHON_VERSION 3 CACHE STRING "Python version to use for compiling modules")
if((NOT (DEFINED pybind11_FOUND)) OR (NOT ${pybind11_FOUND}))
  find_package(pybind11 REQUIRED)
  if(${pybind11_VERSION} VERSION_LESS 2.2)
    message(FATAL_ERROR "pybind11 version should 2.2.x at least")
  endif(${pybind11_VERSION} VERSION_LESS 2.2)
endif((NOT (DEFINED pybind11_FOUND)) OR (NOT ${pybind11_FOUND}))
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

pybind11_add_module(_MeshTools meshtools-module.cpp meshtools-wrapper.cpp)
set_target_properties(_MeshTools
	PROPERTIES
	CXX_STANDARD 14
	CXX_STANDARD_REQUIRED YES
	CXX_EXTENSIONS NO
	LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}
	LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PROJECT_SOURCE_DIR}
	LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PROJECT_SOURCE_DIR}
)

if(${MESHTOOLS_MASTER_PROJECT})
    add_executable(test-meshtools test-meshtools.cpp)
    pybind11_add_module(consumer consumer-module.cpp)
    set_target_properties(consumer
        PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}
        LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PROJECT_SOURCE_DIR}
        LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PROJECT_SOURCE_DIR}
    )
endif(${MESHTOOLS_MASTER_PROJECT})