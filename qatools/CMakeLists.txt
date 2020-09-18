project(qatools)

if(NOT INSTALL_PYTHON_DEPENDENCIES)
  set(PIP_FLAGS "--no-dependencies")
endif()

# run 'pip install' and direct package installs to freesufer/python/packages
install(CODE "
  message(STATUS \"Installing QA tools\")
  execute_process(COMMAND bash -c \"${PYTHON_EXECUTABLE} -m pip install ${CMAKE_CURRENT_SOURCE_DIR} --no-dependencies --upgrade --target=${CMAKE_INSTALL_PREFIX}/python/packages\" RESULT_VARIABLE retcode)
  if(NOT \${retcode} STREQUAL 0)
    # message(FATAL_ERROR \"Could not install qatools python module\")
    message(SEND_ERROR \"Could not install qatools python module\")
  endif()"
)

install_pyscript(qatools.py)

# install any package requirements
if(INSTALL_PYTHON_DEPENDENCIES)
  if(NOT DISTRIBUTE_FSPYTHON)
    set(PKG_TARGET "--target=${CMAKE_INSTALL_PREFIX}/python/packages")
  endif()
  install(CODE "
    message(STATUS \"Installing QA tools package dependencies\")
    execute_process(COMMAND bash -c \"${CMAKE_INSTALL_PREFIX}/python/bin/python3 -m pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt ${PKG_TARGET}\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not install QA tools dependencies\")
    endif()
    execute_process(COMMAND bash -c \"find ${CMAKE_INSTALL_PREFIX}/python/bin ! -name 'python3*' -type f -exec rm -f {} +\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not remove pip-installed scripts\")
    endif()"
  )
endif()
