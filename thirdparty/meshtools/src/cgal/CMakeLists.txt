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

pybind11_add_module(C3t3Wrapper C3t3Wrapper-module.cpp)
set_target_properties(C3t3Wrapper
  PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  LIBRARY_OUTPUT_DIRECTORY ${MESHTOOLS_PYTHON_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY_DEBUG ${MESHTOOLS_PYTHON_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY_RELEASE ${MESHTOOLS_PYTHON_DIRECTORY}
)

pybind11_add_module(DTMWrapper DTMWrapper.cpp)
set_target_properties(DTMWrapper
  PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  LIBRARY_OUTPUT_DIRECTORY ${MESHTOOLS_PYTHON_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY_DEBUG ${MESHTOOLS_PYTHON_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY_RELEASE ${MESHTOOLS_PYTHON_DIRECTORY}
)
