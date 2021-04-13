project(fspython)

add_subdirectories(bindings)

# run 'pip install' and direct package installs to freesurfer/python/packages
install(CODE "
  message(STATUS \"Installing freesurfer python package\")
  execute_process(COMMAND bash -c \"${PYTHON_EXECUTABLE} -m pip install ${CMAKE_CURRENT_SOURCE_DIR} --no-dependencies --upgrade --target=${CMAKE_INSTALL_PREFIX}/python/packages\" RESULT_VARIABLE retcode)
  if(NOT \${retcode} STREQUAL 0)
    message(FATAL_ERROR \"Could not install freesurfer python module\")
  endif()"
)

# install the fspython script
install(PROGRAMS fspython DESTINATION bin)

# install any package requirements
if(INSTALL_PYTHON_DEPENDENCIES)

  # packages not in requirements.txt
  set(EXTRA_PACKAGES " \
    tensorflow==2.4.1 \
    git+git://github.com/adalca/pystrum.git@3befbd1166595476e8315fef99bea373bfdc302c \
    git+git://github.com/adalca/neurite.git@ac74518f367113ddadca90b70acf12f6a8526839 \
    git+git://github.com/voxelmorph/voxelmorph.git@80d0c489febfb4fa32b4a247629e79720fbb4c14"
  )

  if(NOT DISTRIBUTE_FSPYTHON)
    set(PKG_TARGET "--target=${CMAKE_INSTALL_PREFIX}/python/packages")
  endif()
  install(CODE "
    message(STATUS \"Installing freesurfer python package dependencies\")
    execute_process(COMMAND bash -c \"${CMAKE_INSTALL_PREFIX}/python/bin/python3 -m pip install --disable-pip-version-check ${EXTRA_PACKAGES} -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt ${PKG_TARGET}\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not install freesurfer python dependencies\")
    endif()
    execute_process(COMMAND bash -c \"find ${CMAKE_INSTALL_PREFIX}/python/bin ! -name 'python3*' -type f -exec rm -f {} +\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not remove pip-installed scripts\")
    endif()"
  )
endif()
