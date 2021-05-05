project(qatools)

if(NOT INSTALL_PYTHON_DEPENDENCIES)
  set(PIP_FLAGS "--no-dependencies")
endif()

# run 'pip install' and direct package installs to freesufer/python/packages
# requirements are currently defined by python/requirements-extra.txt

install(CODE "
  message(STATUS \"Installing QA tools\")
  execute_process(COMMAND bash -c \"${PYTHON_EXECUTABLE} -m pip install ${CMAKE_CURRENT_SOURCE_DIR} --no-dependencies --upgrade --target=${CMAKE_INSTALL_PREFIX}/python/packages\" RESULT_VARIABLE retcode)
  if(NOT \${retcode} STREQUAL 0)
    # message(FATAL_ERROR \"Could not install qatools python module\")
    message(SEND_ERROR \"Could not install qatools python module\")
  endif()"
)

install_pyscript(qatools.py)
