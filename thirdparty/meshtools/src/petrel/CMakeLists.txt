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

#option(ACTIVATE_CONCURRENCY "Enable concurrency" ON)

pybind11_add_module(PetrelMesh petrel-mesh.cpp)
#target_compile_definitions(PetrelMesh PRIVATE BOOST_PARAMETER_MAX_ARITY=12)
set_target_properties(PetrelMesh
  PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  LIBRARY_OUTPUT_DIRECTORY ${MESHTOOLS_PYTHON_MODULES_OUTPUT_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY_DEBUG ${MESHTOOLS_PYTHON_MODULES_OUTPUT_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY_RELEASE ${MESHTOOLS_PYTHON_MODULES_OUTPUT_DIRECTORY}
)

#if(ACTIVATE_CONCURRENCY)
#  find_package( TBB REQUIRED )
#  if( TBB_FOUND )
#    add_definitions(-DCGAL_CONCURRENT_MESH_3)
#    include(${TBB_USE_FILE})
#    target_link_libraries(SharpEdges ${TBB_LIBRARIES})
#  endif(TBB_FOUND)
#endif(ACTIVATE_CONCURRENCY)
