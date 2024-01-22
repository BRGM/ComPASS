pybind11_add_module(metis metis.cpp)
target_link_libraries(metis PRIVATE dep-metis)
install(
  TARGETS metis
  DESTINATION ${COMPASS_PYTHON_PACKAGE_NAME}
  EXPORT ComPASS
)
