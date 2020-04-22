project(pysurfer)

add_subdirectories(bindings)

if(NOT INSTALL_PYTHON_DEPENDENCIES)
  set(PIP_FLAGS "--no-dependencies")
endif()

# run 'pip install' and direct package installs to freesufer/python/packages
install(CODE "
  message(STATUS \"Installing freesurfer python package\")
  execute_process(COMMAND bash -c \"${PYTHON_EXECUTABLE} -m pip install ${CMAKE_CURRENT_SOURCE_DIR} ${PIP_FLAGS} --upgrade --target=${CMAKE_INSTALL_PREFIX}/python/packages\" RESULT_VARIABLE retcode)
  if(NOT \${retcode} STREQUAL 0)
    # message(FATAL_ERROR \"Could not install freesurfer python module\")
    message(SEND_ERROR \"Could not install freesurfer python module\")
  endif()"
)

install(PROGRAMS fspython DESTINATION bin)
