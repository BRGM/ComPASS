function(append_sources)
  file(RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}"
       "${CMAKE_CURRENT_SOURCE_DIR}"
  )
  # message(STATUS "arguments " ${ARGV}) message(STATUS "arguments " ${ARGN})
  # message(STATUS "arguments " ${ARGC}) set(TOTO A B C) message(STATUS
  # "variable " ${TOTO}) foreach(_elem ${ARGV}) message(STATUS "iter: "
  # ${_elem}) endforeach()
  list(GET ARGV 0 COLLECTION_NAME)
  # message(STATUS "extracted varname: " ${COLLECTION_NAME})
  file(RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}"
       "${CMAKE_CURRENT_SOURCE_DIR}"
  )
  list(REMOVE_AT ARGV 0)
  foreach(_src ${ARGV})
    if(_relPath)
      list(APPEND COLLECTION "${_relPath}/${_src}")
    else()
      list(APPEND COLLECTION "${_src}")
    endif()
  endforeach()
  if(_relPath)
    # propagate collection to parent directory
    set(${COLLECTION_NAME}
        ${${COLLECTION_NAME}} ${COLLECTION}
        PARENT_SCOPE
    )
  endif()
endfunction(append_sources)
