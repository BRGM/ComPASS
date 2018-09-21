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

pybind11_add_module(CGALWrappers CGALWrappers.cpp C3t3Wrapper-module.cpp DTMWrapper.cpp
                    implicit_functions.cpp mesh_implicit_domains.cpp)
target_compile_definitions(CGALWrappers PRIVATE BOOST_PARAMETER_MAX_ARITY=12)
set_target_properties(CGALWrappers
  PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  LIBRARY_OUTPUT_DIRECTORY ${MESHTOOLS_PYTHON_MODULES_OUTPUT_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY_DEBUG ${MESHTOOLS_PYTHON_MODULES_OUTPUT_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY_RELEASE ${MESHTOOLS_PYTHON_MODULES_OUTPUT_DIRECTORY}
)

pybind11_add_module(Corefinement corefinement.cpp)
#target_compile_definitions(Corefinement PRIVATE BOOST_PARAMETER_MAX_ARITY=12)
set_target_properties(Corefinement
  PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  LIBRARY_OUTPUT_DIRECTORY ${MESHTOOLS_PYTHON_MODULES_OUTPUT_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY_DEBUG ${MESHTOOLS_PYTHON_MODULES_OUTPUT_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY_RELEASE ${MESHTOOLS_PYTHON_MODULES_OUTPUT_DIRECTORY}
)

option(ACTIVATE_CONCURRENCY "Enable concurrency" ON)

if( ACTIVATE_CONCURRENCY )
    find_package( TBB )
    if( NOT TBB_FOUND )
        message( STATUS "NOTICE: Intel TBB NOT found! The example is faster if TBB is linked." )
    endif()
endif()

if(ACTIVATE_CONCURRENCY AND TBB_FOUND)
    CGAL_target_use_TBB(CGALWrappers)
endif()
