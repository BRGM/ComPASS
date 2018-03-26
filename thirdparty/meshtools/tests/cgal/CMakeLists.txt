include( ${CGAL_USE_FILE} )
include_directories(${MESHTOOLS_CGAL_DIRECTORY})

add_executable(test-DTM test-DTM.cpp)
