# CGAL and its components
if(MESHTOOLS_USES_CGAL)
    add_subdirectory(${MESHTOOLS_CGAL_DIRECTORY})
    add_definitions(-DMESHTOOLS_USES_CGAL)
    add_subdirectory(petrel)
endif(MESHTOOLS_USES_CGAL)

add_subdirectory(core)
