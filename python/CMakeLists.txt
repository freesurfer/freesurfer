project(fspython)

# build the cxx/python bindings for the below modules
add_subdirectories(fsbindings)
add_subdirectories(gems/cxx)

# run 'pip install' and direct package installs to freesurfer/python/packages
install(CODE "
  message(STATUS \"Installing python packages\")
  execute_process(COMMAND bash -c \"${PYTHON_EXECUTABLE} -m pip install ${CMAKE_CURRENT_SOURCE_DIR} --no-dependencies --upgrade --target=${CMAKE_INSTALL_PREFIX}/python/packages\" RESULT_VARIABLE retcode)
  if(NOT \${retcode} STREQUAL 0)
    message(FATAL_ERROR \"Could not install python modules\")
  endif()"
  COMPONENT fspy
)

# install the fspython script
if(APPLE_ARM64)
   # Allow script to use local packages with the system or Apple recommended Anaconda python
   install(PROGRAMS fspython.local DESTINATION bin RENAME fspython)
else()
   install(PROGRAMS fspython DESTINATION bin)
endif()

# install any package requirements
if(INSTALL_PYTHON_DEPENDENCIES)

  if(NOT DISTRIBUTE_FSPYTHON)
    set(PKG_TARGET "--target=${CMAKE_INSTALL_PREFIX}/python/packages")
  endif()

  install(CODE "
    message(STATUS \"Installing python package dependencies\")
    execute_process(COMMAND bash -c \"${CMAKE_INSTALL_PREFIX}/python/bin/python3 -m pip install --upgrade --disable-pip-version-check ${EXTRA_PACKAGES} -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-extra.txt ${PKG_TARGET}\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not install python dependencies\")
    endif()
    execute_process(COMMAND bash -c \"find ${CMAKE_INSTALL_PREFIX}/python/bin ! -name 'python3*' -type f -exec rm -f {} +\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not remove pip-installed scripts\")
    endif()"
  )

  if(PATCH_FSPYTHON)
     # Use to workaround a python module not installing, e.g., problem with setup.py.
     # Enable with -DPATCH_FSPYTHON=ON as arg to cmake command to copy install from previous build.
     # The path to the previous build must be available/mounted in the build container.
     set(FSPYTHON_PREVIOUS_BUILD_BASE "/space/freesurfer/build/garnet/sandboxes/utility/freesurfer/dev.20221207")
     install(CODE "
       message(STATUS \"patching fspython distribution\")
       execute_process(COMMAND bash -c \" ( cd ${FSPYTHON_PREVIOUS_BUILD_BASE}/python/lib/python3.8/site-packages && tar cpf - surfa* xxhash* ) | ( cd ${CMAKE_INSTALL_PREFIX}/python/lib/python3.8/site-packages && tar xpf - ) \" RESULT_VARIABLE patch)
       if(NOT \${patch} STREQUAL 0)
         message(FATAL_ERROR \"Could not patch fspython\")
       endif()"
     )
  endif()

endif()
