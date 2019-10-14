pybind11_add_module(bind_petsc bind_petsc.cpp VariousDumps.F90)
set_target_properties(bind_petsc
    PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    LIBRARY_OUTPUT_DIRECTORY ${COMPASS_PYTHON_DIRECTORY}
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${COMPASS_PYTHON_DIRECTORY}
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${COMPASS_PYTHON_DIRECTORY}
)
target_include_directories(bind_petsc PUBLIC  ${PROJECT_SOURCE_DIR}/../wrappers)
target_include_directories(bind_petsc PUBLIC ${pybind11_INCLUDE_DIR})
target_include_directories(bind_petsc PUBLIC ${MPI_CXX_INCLUDE_DIRS})
target_include_directories(bind_petsc PUBLIC ${PETSc_INCLUDE_DIRS})
target_link_libraries(bind_petsc PRIVATE ${PETSc_LIBRARIES} ${LAPACK_LIBRARIES})
