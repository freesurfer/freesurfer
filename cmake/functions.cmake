
function(install_symlinks)
  set(flags DESTINATION TYPE)
  cmake_parse_arguments(INSTALL "" "${flags}" "" ${ARGN})
  string(TOUPPER ${INSTALL_TYPE} INSTALL_TYPE)
  foreach(arg ${INSTALL_UNPARSED_ARGUMENTS})
    get_filename_component(BASENAME ${arg} NAME)
    get_filename_component(ABS_PATH ${arg} REALPATH)
    install(${INSTALL_TYPE} ${ABS_PATH} DESTINATION ${INSTALL_DESTINATION} RENAME ${BASENAME})
  endforeach()
endfunction()


function(install_tarball TARBALL DESTINATION_DIR)
  get_filename_component(ABS_TARBALL ${TARBALL} ABSOLUTE)
  install(DIRECTORY DESTINATION ${DESTINATION_DIR})
  install(CODE "execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar -xzf ${ABS_TARBALL}
    WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${DESTINATION_DIR})"
  )
endfunction()