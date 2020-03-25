add_executable(scorep-python scorep-python.cpp)
target_link_libraries(scorep-python PRIVATE pybind11::embed)
set_target_properties(scorep-python
    PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

install(TARGETS scorep-python DESTINATION ${PROJECT_SOURCE_DIR}/bin EXPORT ComPASS)
