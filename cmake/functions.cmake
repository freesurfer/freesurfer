# This file defines a few custom cmake utility functions to
# simplify the freesurfer build process


# add_subdirectories(<subdirs>)
# Simple utility to add multiple subdirectories at once
function(add_subdirectories)
  foreach(SUBDIR ${ARGN})
    add_subdirectory(${SUBDIR})
  endforeach()
endfunction()


# install_configured(<files>)
# Installs a file configured with cmake variables
function(install_configured)
  cmake_parse_arguments(INSTALL "" "DESTINATION" "" ${ARGN})
  foreach(FILE ${INSTALL_UNPARSED_ARGUMENTS})
    install(CODE "
      message(STATUS \"Configuring: ${CMAKE_INSTALL_PREFIX}/${INSTALL_DESTINATION}/${FILE}\")
      file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${INSTALL_DESTINATION})
      set(FS_VERSION ${FS_VERSION})
      set(BUILD_STAMP ${BUILD_STAMP})
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${FILE} ${CMAKE_INSTALL_PREFIX}/${INSTALL_DESTINATION} @ONLY)"
    )
  endforeach()
endfunction()


# install_symlinks(<files> DESTINATION <dir> TYPE <type>)
# Unfortunately, cmake does not follow symlinks when installing files,
# so this is a wrapper for the install() function that will actually follow
# links. The TYPE argument should specify the install type that would normally
# be provided to the install function (i.e. PROGRAMS, FILES, etc.)
function(install_symlinks)
  cmake_parse_arguments(INSTALL "NMR_ONLY" "DESTINATION;TYPE" "" ${ARGN})
  string(TOUPPER ${INSTALL_TYPE} INSTALL_TYPE)
  if(INSTALL_NMR_ONLY)
    set(COMPONENT_ARGS COMPONENT nmr EXCLUDE_FROM_ALL)
  endif()
  foreach(arg ${INSTALL_UNPARSED_ARGUMENTS})
    get_filename_component(BASENAME ${arg} NAME)
    get_filename_component(ABS_PATH ${arg} REALPATH)
    install(${INSTALL_TYPE} ${ABS_PATH} DESTINATION ${INSTALL_DESTINATION} RENAME ${BASENAME} ${COMPONENT_ARGS})
  endforeach()
endfunction()


# install_tarball(<tarball> DESTINATION <dir>)
# Use this to extract a tarball to the provided install destination
function(install_tarball)
  cmake_parse_arguments(INSTALL "" "DESTINATION" "" ${ARGN})
  get_filename_component(TARBALL ${INSTALL_UNPARSED_ARGUMENTS} ABSOLUTE)
  install(CODE "
    message(STATUS \"Extracting: ${TARBALL} to ${CMAKE_INSTALL_PREFIX}/${INSTALL_DESTINATION}\")
    file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${INSTALL_DESTINATION})
    execute_process(
      COMMAND bash -c \"tar -xzf ${TARBALL}\"
      WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${INSTALL_DESTINATION}
      RESULT_VARIABLE retcode
    )
    if(NOT \${retcode} STREQUAL 0)
      # message(FATAL_ERROR \"Could not extract ${INSTALL_UNPARSED_ARGUMENTS} - perhaps it hasn't been downloaded from the annex?\")
      message(SEND_ERROR \"Could not extract ${INSTALL_UNPARSED_ARGUMENTS} - perhaps it hasn't been downloaded from the annex?\")
    endif()"
  )
endfunction()


# mac_deploy_qt(TARGET <target> BUNDLE <bundle> PLIST <plist> ICONS <icons>)
# Creates a mac app bundle from a given target. The ICONS argument is optional.
function(mac_deploy_qt)
  cmake_parse_arguments(APP "" "TARGET;BUNDLE;PLIST;ICONS" "" ${ARGN})
  # install binary
  install(TARGETS ${APP_TARGET} DESTINATION ${APP_BUNDLE}/Contents/MacOS)
  # install the plist
  install(FILES ${APP_PLIST} DESTINATION ${APP_BUNDLE}/Contents)
  # install the resources
  if(APP_ICONS)
    install_symlinks(${APP_ICONS} TYPE files DESTINATION ${APP_BUNDLE}/Contents/Resources)
  endif()
  # run the qt deployment script
  install(CODE "
    message(STATUS \"Deploying ${APP_BUNDLE}\")
    execute_process(COMMAND bash -c \"${CMAKE_SOURCE_DIR}/qt/mac_deploy ${Qt5_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/${APP_BUNDLE}\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      # message(FATAL_ERROR \"Could not deploy ${APP_TARGET}\")
      message(SEND_ERROR \"Could not deploy ${APP_TARGET}\")
    endif()"
  )
endfunction()


# add_help(<binary> <xml>)
# Link an xml helptext to a target binary. This will create a target dependency on
# the help file and will run xxd to create the xml header during the build
function(add_help BINARY HELPTEXT)
  add_custom_command(COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} &&
                             xxd -i ${HELPTEXT} ${CMAKE_CURRENT_BINARY_DIR}/${HELPTEXT}.h
                     OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${HELPTEXT}.h
                     DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${HELPTEXT})
  include_directories(${CMAKE_CURRENT_BINARY_DIR})
  target_sources(${BINARY} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/${HELPTEXT}.h)
  install(FILES ${HELPTEXT} DESTINATION docs/xml)
  # make sure to validate the xml as well
  add_test(${BINARY}_help_test bash -c "xmllint --noout ${CMAKE_CURRENT_SOURCE_DIR}/${HELPTEXT}")
endfunction(add_help)


# install_append_help(<script> <xml> <destination>)
# Some shell scripts also utilize a help.xml file by appending the helptext to the end of it.
# This function will setup this build command and dependency
function(install_append_help SCRIPT HELPTEXT DESTINATION)
  install_configured(${SCRIPT} DESTINATION ${DESTINATION})
  install(CODE "
    message(STATUS \"Appending helptext: ${CMAKE_INSTALL_PREFIX}/${DESTINATION}/${SCRIPT}\")
    execute_process(
      COMMAND bash -c \"
        ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}/utils --target fsPrintHelp &&
        ${CMAKE_BINARY_DIR}/utils/fsPrintHelp ${CMAKE_CURRENT_SOURCE_DIR}/${HELPTEXT} >> ${CMAKE_INSTALL_PREFIX}/${DESTINATION}/${SCRIPT}\"
      OUTPUT_QUIET
      RESULT_VARIABLE retcode
    )
    if(NOT \${retcode} STREQUAL 0)
      # message(FATAL_ERROR \"Could not append help text to ${SCRIPT}\")
      message(SEND_ERROR \"Could not append help text to ${SCRIPT}\")
    endif()"
  )
  install(FILES ${HELPTEXT} DESTINATION docs/xml)
  # make sure to validate the xml as well
  add_test(${SCRIPT}_help_test bash -c "xmllint --noout ${CMAKE_CURRENT_SOURCE_DIR}/${HELPTEXT}")
endfunction()

# install_osx_app(<app>)
# This copies a pre-built os x application to bin (for Eugenio's subfield applications)
function(install_osx_app APP_PATH)
  get_filename_component(APP_NAME ${APP_PATH} NAME)
  install(CODE "
    message(STATUS \"Copying OS X Application: ${APP_NAME} to ${CMAKE_INSTALL_PREFIX}/bin/${APP_NAME}\")
    execute_process(
      COMMAND bash -c \"
        cp -RL ${CMAKE_CURRENT_SOURCE_DIR}/${APP_PATH} ${CMAKE_INSTALL_PREFIX}/bin &&
        chmod -R 755 ${CMAKE_INSTALL_PREFIX}/bin/${APP_NAME}\"
      RESULT_VARIABLE retcode
    )
    if(NOT \${retcode} STREQUAL 0)
      # message(FATAL_ERROR \"Could not install ${APP_NAME}\")
      message(SEND_ERROR \"Could not install ${APP_NAME}\")
    endif()"
  )
endfunction()


# symlink(<target> <linkname>)
# Creates a symlink at install time
function(symlink TARGET LINKNAME)
  install(CODE "
    message(STATUS \"Symlinking: ${LINKNAME} to ${TARGET}\")
    execute_process(COMMAND bash -c \"mkdir -p $(dirname ${LINKNAME}) && rm -f ${LINKNAME} && ln -s ${TARGET} ${LINKNAME}\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      # message(FATAL_ERROR \"Could not create symlink to ${TARGET}\")
      message(SEND_ERROR \"Could not create symlink to ${TARGET}\")
    endif()"
  )
endfunction()


# add_test_script(NAME <testname> SCRIPT <testscript> DEPENDS <targets>)
# Adds a script to the test framework. If the script calls a freesurfer binary, it
# should be specified with DEPENDS to guarantee it gets built beforehand
function(add_test_script)
  cmake_parse_arguments(TEST "" "NAME;SCRIPT" "DEPENDS" ${ARGN})
  # foreach(TARGET ${TEST_DEPENDS})
  #  set(TEST_CMD "${TEST_CMD} ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR} --target ${TARGET} &&")
  # endforeach()
  add_test(${TEST_NAME} bash -c "${TEST_CMD} ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_SCRIPT}")
endfunction()


# add_test_executable(<target> <sources>)
# Adds an executable to the test framework that only gets built by the test target.
# This function takes the same input as add_executable()
function(add_test_executable)
  set(TARGET ${ARGV0})
  LIST(REMOVE_AT ARGV 0)
  add_executable(${TARGET} EXCLUDE_FROM_ALL ${ARGV})
  # set(TEST_CMD "${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR} --target ${TARGET} &&")
  # add_test(${TARGET} bash -c "${TEST_CMD} ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}")
endfunction()


# library_paths(NAME <var> LIBDIR <library-dir> LIBRARIES <libs>)
# This function is meant to find valid paths of libraries in a given library directory.
# It's mostly a utility for locating VTK and PETSC libs that may or may not exist on a
# given platform.
function(library_paths)
  cmake_parse_arguments(ARG "" "NAME;LIBDIR" "LIBRARIES" ${ARGN})
  unset(LIBPATH)
  foreach(LIBRARY ${ARG_LIBRARIES})
    find_library(LIBPATH PATHS ${ARG_LIBDIR} NAMES ${LIBRARY} NO_DEFAULT_PATH)
    if(LIBPATH)
      set(LIB_LIST ${LIB_LIST} ${LIBPATH})
    endif()
    unset(LIBPATH CACHE)
  endforeach()
  set(${ARG_NAME} ${LIB_LIST} PARENT_SCOPE)
endfunction()


# pyscript(<scripts>)
# This function installs and wraps python scripts that are meant to be run
# with fspython. It will install the actual scripts to python/scripts and create
# a wrapper with the same names in the bin directory. For example:
#     pyscript(samseg)
# creates the bash script $FREESURFER_HOME/bin/samseg that calls the real samseg
# with the correct fspython distribution and packages.
function(install_pyscript)
  foreach(SCRIPT ${ARGN})
    install(FILES ${SCRIPT} DESTINATION python/scripts)
    install(CODE "
      message(STATUS \"Configuring python wrapper: ${CMAKE_INSTALL_PREFIX}/bin/${SCRIPT}\")
      set(SCRIPTNAME ${SCRIPT})
      configure_file(${CMAKE_SOURCE_DIR}/python/wrapper ${CMAKE_INSTALL_PREFIX}/bin/${SCRIPT} @ONLY)"
    )
  endforeach()
endfunction()
