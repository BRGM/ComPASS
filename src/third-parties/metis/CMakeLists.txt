file(GLOB metis_sources *.c)

add_library(dep-metis STATIC ${metis_sources})
target_compile_definitions(dep-metis PUBLIC IDXTYPEWIDTH=32 REALTYPEWIDTH=64)
set_target_properties(dep-metis PROPERTIES POSITION_INDEPENDENT_CODE True)
target_include_directories(dep-metis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(dep-metis PRIVATE GKlib)
