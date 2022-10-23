project(fspython)

# build the cxx/python bindings for the below modules
add_subdirectories(fsbindings)
add_subdirectories(gems/cxx)

if(DISTRIBUTE_FSPYTHON)
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
   install(PROGRAMS fspython DESTINATION bin)

   # install any package requirements
   if(INSTALL_PYTHON_DEPENDENCIES)

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
   endif()
else()
   if(INSTALL_PYTHON_DEPENDENCIES)
      set(PKG_TARGET "--target=${CMAKE_INSTALL_PREFIX}/python/packages")
   endif()
endif()
