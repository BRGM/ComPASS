add_custom_target(
  sitesconf
  COMMAND python3 parse_conf.py
  DEPENDS conf.yaml
  BYPRODUCTS generated-sites.h generated-sites.F90
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(sites-fortran test-sites.F90 sites.F90)
add_dependencies(sites-fortran sitesconf)

add_executable(sites-cpp test-sites.cpp sites.h)
add_dependencies(sites-cpp sitesconf)
