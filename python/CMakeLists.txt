project(pysurfer)

install(CODE "
  message(STATUS \"Installing: ${CMAKE_INSTALL_PREFIX}/python/freesurfer\")
  execute_process(COMMAND bash -c \"pip install ${CMAKE_CURRENT_SOURCE_DIR} --upgrade --target=${CMAKE_INSTALL_PREFIX}/python\" RESULT_VARIABLE retcode)
  if(NOT \${retcode} STREQUAL 0)
    message(FATAL_ERROR \"Could not install freesurfer python module\")
  endif()"
)
